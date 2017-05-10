/*
   Copyright 2017 University of Washington                          3-clause BSD license

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "Common.h"
#include "mango_Preprocess.h"
#include "mango_DataInternal.h"

//std::vector<Query*>           g_pvQuery;
//std::vector<InputFileInfo *>  g_pvInputFiles;
//StaticParams                  g_staticParams;
//MassRange                     g_massRange;


bool mango_preprocess::_bDoneProcessingAllSpectra;
bool mango_preprocess::_bFirstScan;
bool *mango_preprocess::pbMemoryPool;
double **mango_preprocess::ppdTmpRawDataArr;
double **mango_preprocess::ppdTmpFastXcorrDataArr;
double **mango_preprocess::ppdTmpCorrelationDataArr;

mango_preprocess::mango_preprocess()
{
}


mango_preprocess::~mango_preprocess()
{
}

void mango_preprocess::Reset()
{
    _bFirstScan = true;
    _bDoneProcessingAllSpectra = false;
}

void mango_preprocess::LoadAndPreprocessSpectra(Spectrum *mstSpectrum)
{
   g_massRange.iMaxFragmentCharge = 0;

   if (mstSpectrum->getScanNumber() != 0)   // should not be needed by quick sanity check to make sure scan is read
   {
      int iNumClearedPeaks = 0;

      // Clear out m/z range if clear_mz_range parameter is specified
      // Accomplish this by setting corresponding intensity to 0
      if (g_staticParams.options.clearMzRange.dEnd > 0.0
            && g_staticParams.options.clearMzRange.dStart <= g_staticParams.options.clearMzRange.dEnd)
      {
         int i=0;

         while (true)
         {
            if (i >= mstSpectrum->size() || mstSpectrum->at(i).mz > g_staticParams.options.clearMzRange.dEnd)
               break;

            if (mstSpectrum->at(i).mz >= g_staticParams.options.clearMzRange.dStart
                  && mstSpectrum->at(i).mz <= g_staticParams.options.clearMzRange.dEnd)
            {
               mstSpectrum->at(i).intensity = 0.0;
               iNumClearedPeaks++;
            }

            i++;
         }
      }

      if (mstSpectrum->size()-iNumClearedPeaks >= g_staticParams.options.iMinPeaks)
      {
         if (CheckActivationMethodFilter(mstSpectrum->getActivationMethod()))
         {
            int i=0;

            PreprocessSpectrum(*mstSpectrum,
                  ppdTmpRawDataArr[i],
                  ppdTmpFastXcorrDataArr[i],
                  ppdTmpCorrelationDataArr[i]);
         }
      }
   }
}


bool mango_preprocess::DoneProcessingAllSpectra()
{
   return _bDoneProcessingAllSpectra;
}


bool mango_preprocess::Preprocess(struct Query *pScoring,
                                  Spectrum mstSpectrum,
                                  double *pdTmpRawData,
                                  double *pdTmpFastXcorrData,
                                  double *pdTmpCorrelationData)
{
   int i;
   int x;
   int y;
   struct PreprocessStruct pPre;

   pPre.iHighestIon = 0;
   pPre.dHighestIntensity = 0;

   //MH: Find appropriately sized array cushion based on user parameters. Fixes error found by Patrick Pedrioli for
   // very wide mass tolerance searches (i.e. 500 Da).
   double dCushion = 0.0;
   if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
   {
      dCushion = g_staticParams.tolerances.dTolerancePeptide;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
        dCushion *= 8; //MH: hope +8 is large enough charge because g_staticParams.options.iEndCharge can be overridden.
      }
   }
   else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
   {
      dCushion = g_staticParams.tolerances.dTolerancePeptide * 0.001;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         dCushion *= 8; //MH: hope +8 is large enough charge because g_staticParams.options.iEndCharge can be overridden.
      }
   }
   else // ppm
   {
      dCushion = g_staticParams.tolerances.dTolerancePeptide * g_staticParams.options.dPeptideMassHigh / 1000000.0;
   }

   // initialize these temporary arrays before re-using
   size_t iTmp = (size_t)((g_staticParams.options.dPeptideMassHigh + dCushion + 2.0) * g_staticParams.dInverseBinWidth)*sizeof(double);

   memset(pdTmpRawData, 0, iTmp);
   memset(pdTmpFastXcorrData, 0, iTmp);
   memset(pdTmpCorrelationData, 0, iTmp);

   // pdTmpRawData is a binned array holding raw data
   if (!LoadIons(pScoring, pdTmpRawData, mstSpectrum, &pPre))
   {
      return false;
   }

   float pfFastXcorrData[pScoring->_spectrumInfoInternal.iArraySize];

   // Create data for correlation analysis.
   // pdTmpRawData intensities are normalized to 100; pdTmpCorrelationData is windowed
   MakeCorrData(pdTmpRawData, pdTmpCorrelationData, pScoring, &pPre);

   // Make fast xcorr spectrum.
   double dSum=0.0;
   int iTmpRange = 2*g_staticParams.iXcorrProcessingOffset + 1;
   double dTmp = 1.0 / (double)(iTmpRange - 1);

   dSum=0.0;
   for (i=0; i<g_staticParams.iXcorrProcessingOffset; i++)
      dSum += pdTmpCorrelationData[i];
   for (i=g_staticParams.iXcorrProcessingOffset; i < pScoring->_spectrumInfoInternal.iArraySize + g_staticParams.iXcorrProcessingOffset; i++)
   {
      if (i<pScoring->_spectrumInfoInternal.iArraySize)
         dSum += pdTmpCorrelationData[i];
      if (i>=iTmpRange && i-iTmpRange < pScoring->_spectrumInfoInternal.iArraySize)
         dSum -= pdTmpCorrelationData[i-iTmpRange];
      if (i-g_staticParams.iXcorrProcessingOffset >=0 && i-g_staticParams.iXcorrProcessingOffset < pScoring->_spectrumInfoInternal.iArraySize)
         pdTmpFastXcorrData[i-g_staticParams.iXcorrProcessingOffset] = (dSum - pdTmpCorrelationData[i-g_staticParams.iXcorrProcessingOffset])* dTmp;
   }

   pfFastXcorrData[0] = 0.0;
   for (i=1; i<pScoring->_spectrumInfoInternal.iArraySize; i++)
   {
      double dTmp = pdTmpCorrelationData[i] - pdTmpFastXcorrData[i];

      pfFastXcorrData[i] = (float)dTmp;

      // Add flanking peaks if used
      if (g_staticParams.ionInformation.iTheoreticalFragmentIons == 0)
      {
         int iTmp;

         iTmp = i-1;
         pfFastXcorrData[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp])*0.5);

         iTmp = i+1;
         if (iTmp < pScoring->_spectrumInfoInternal.iArraySize)
            pfFastXcorrData[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp])*0.5);
      }
   }

   pScoring->iFastXcorrData=pScoring->_spectrumInfoInternal.iArraySize/SPARSE_MATRIX_SIZE+1;

   //MH: Fill sparse matrix
   try
   {
      pScoring->ppfSparseFastXcorrData = new float*[pScoring->iFastXcorrData]();
   }
   catch (std::bad_alloc& ba)
   {
      fprintf(stderr, " Error - new(pScoring->ppfSparseFastXcorrData[%d]). bad_alloc: %s.\n", pScoring->iFastXcorrData, ba.what());
      fprintf(stderr, " mango ran out of memory. Look into \"spectrum_batch_size\"\n");
      fprintf(stderr, " parameters to address mitigate memory use.\n");
      return false;
   }

   for (i=1; i<pScoring->_spectrumInfoInternal.iArraySize; i++)
   {
      if (pfFastXcorrData[i]>FLOAT_ZERO || pfFastXcorrData[i]<-FLOAT_ZERO)
      {
         x=i/SPARSE_MATRIX_SIZE;
         if (pScoring->ppfSparseFastXcorrData[x]==NULL)
         {
            try
            {
               pScoring->ppfSparseFastXcorrData[x] = new float[SPARSE_MATRIX_SIZE]();
            }
            catch (std::bad_alloc& ba)
            {
               fprintf(stderr, " Error - new(pScoring->ppfSparseFastXcorrData[%d][%d]). bad_alloc: %s.\n", x, SPARSE_MATRIX_SIZE, ba.what());
               fprintf(stderr, " mango ran out of memory. Look into \"spectrum_batch_size\"\n");
               fprintf(stderr, " parameters to address mitigate memory use.\n");
               return false;
            }
            for (y=0; y<SPARSE_MATRIX_SIZE; y++)
               pScoring->ppfSparseFastXcorrData[x][y]=0;
         }
         y=i-(x*SPARSE_MATRIX_SIZE);
         pScoring->ppfSparseFastXcorrData[x][y] = pfFastXcorrData[i];
      }
   }

   return true;
}


bool mango_preprocess::CheckActivationMethodFilter(MSActivation act)
{
   bool bSearchSpectrum = true;

   // Check possible activation method filter.
   if (strcmp(g_staticParams.options.szActivationMethod, "ALL")!=0 && (act != mstNA))
   {
      if (!strcmp(g_staticParams.options.szActivationMethod, "CID") && (act != mstCID))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "HCD") && (act != mstHCD))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "ETD") && (act != mstETD))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "ECD") && (act != mstECD))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "PQD") && (act != mstPQD))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "IRMPD") && (act != mstIRMPD))
      {
         bSearchSpectrum = 0;
      }
   }

   return bSearchSpectrum;
}

bool mango_preprocess::CheckExit(int iAnalysisType,
                                int iScanNum,
                                int iTotalScans,
                                int iLastScan,
                                int iReaderLastScan,
                                int iNumSpectraLoaded)
{
   if (iAnalysisType == AnalysisType_SpecificScan)
   {
      _bDoneProcessingAllSpectra = true;
      return true;
   }

   if (iAnalysisType == AnalysisType_SpecificScanRange)
   {
      if (iLastScan > 0)
      {
         if (iScanNum >= iLastScan)
         {
            _bDoneProcessingAllSpectra = true;
            return true;
         }
      }
   }

   if (iAnalysisType == AnalysisType_EntireFile
         && IsValidInputType(g_staticParams.inputFile.iInputType)
         && iScanNum == 0)
   {
      _bDoneProcessingAllSpectra = true;
      return true;
   }

   // Horrible way to exit as this typically requires a quick cycle through
   // while loop but not sure what else to do when getScanNumber() returns 0
   // for non MS/MS scans.
   if (IsValidInputType(g_staticParams.inputFile.iInputType)
         && iTotalScans > iReaderLastScan)
   {
      _bDoneProcessingAllSpectra = true;
      return true;
   }

   if ((g_staticParams.options.iSpectrumBatchSize != 0)
         && (iNumSpectraLoaded >= g_staticParams.options.iSpectrumBatchSize))
   {
      return true;
   }

   return false;
}


bool mango_preprocess::PreprocessSpectrum(Spectrum &spec,
                                         double *pdTmpRawData,
                                         double *pdTmpFastXcorrData,
                                         double *pdTmpCorrelationData)
{
   int z;
   int zStop;

   int iScanNumber = spec.getScanNumber();

   // Set our boundaries for multiple z lines.
   zStop = spec.sizeZ();

   for (z=0; z<zStop; z++)
   {
      int iPrecursorCharge = spec.atZ(z).z;  // I need this before iChargeState gets assigned.
      double dMass = spec.atZ(z).mh;
      Query *pScoring = new Query();

      if (dMass >= g_staticParams.options.dPeptideMassHigh)
         continue;

      pScoring->_pepMassInfo.dExpPepMass = dMass;
      pScoring->_spectrumInfoInternal.iChargeState = iPrecursorCharge;
      pScoring->_spectrumInfoInternal.dTotalIntensity = 0.0;
      pScoring->_spectrumInfoInternal.dRTime = 60.0*spec.getRTime();;
      pScoring->_spectrumInfoInternal.iScanNumber = iScanNumber;

      if (iPrecursorCharge == 1)
         pScoring->_spectrumInfoInternal.iMaxFragCharge = 1;
      else
      {
         pScoring->_spectrumInfoInternal.iMaxFragCharge = iPrecursorCharge - 1;

         if (pScoring->_spectrumInfoInternal.iMaxFragCharge > g_staticParams.options.iMaxFragmentCharge)
            pScoring->_spectrumInfoInternal.iMaxFragCharge = g_staticParams.options.iMaxFragmentCharge;
      }

      //MH: Find appropriately sized array cushion based on user parameters. Fixes error found by Patrick Pedrioli for
      // very wide mass tolerance searches (i.e. 500 Da).
      double dCushion = 0.0;
      if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
      {
         dCushion = g_staticParams.tolerances.dTolerancePeptide;

         if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
         {
            dCushion *= pScoring->_spectrumInfoInternal.iChargeState;
         }
      }
      else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
      {
         dCushion = g_staticParams.tolerances.dTolerancePeptide * 0.001;

         if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
         {
            dCushion *= pScoring->_spectrumInfoInternal.iChargeState;
         }
      }
      else // ppm
      {
         dCushion = g_staticParams.tolerances.dTolerancePeptide * dMass / 1000000.0;
      }
      pScoring->_spectrumInfoInternal.iArraySize = (int)((dMass + dCushion + 2.0) * g_staticParams.dInverseBinWidth);

      // g_massRange.iMaxFragmentCharge is global maximum fragment ion charge across all spectra.
      if (pScoring->_spectrumInfoInternal.iMaxFragCharge > g_massRange.iMaxFragmentCharge)
      {
         g_massRange.iMaxFragmentCharge = pScoring->_spectrumInfoInternal.iMaxFragCharge;
      }

      if (!AdjustMassTol(pScoring))
      {
         return false;
      }

      // Populate pdCorrelation data.
      // NOTE: there must be a good way of doing this just once per spectrum instead
      //       of repeating for each charge state.
      if (!Preprocess(pScoring, spec, pdTmpRawData, pdTmpFastXcorrData, pdTmpCorrelationData))
      {
         return false;
      }

      g_pvQuery.push_back(pScoring);

      //delete pScoring;  //Why does this this causes a segfault??
   }

   return true;
}


// Skip repeating a search if output exists only works for .out files
bool mango_preprocess::CheckExistOutFile(int iCharge,
                                        int iScanNum)
{
   bool bSearchSpectrum = 1;

   if (g_staticParams.options.bOutputOutFiles
         && g_staticParams.options.bSkipAlreadyDone
         && !g_staticParams.options.bOutputPepXMLFile
         && !g_staticParams.options.bOutputPercolatorFile)
   {
      char szOutputFileName[SIZE_FILE];
      char *pStr;
      FILE *fpcheck;

      if ( (pStr = strrchr(g_staticParams.inputFile.szBaseName, '\\')) == NULL
            && (pStr = strrchr(g_staticParams.inputFile.szBaseName, '/')) == NULL)
      {
         pStr = g_staticParams.inputFile.szBaseName;
      }
      else
         (*pStr)++;

      sprintf(szOutputFileName, "%s/%s.%.5d.%.5d.%d.out",
            g_staticParams.inputFile.szBaseName,
            pStr,
            iScanNum,
            iScanNum,
            iCharge);

      // Check existence of .out file.
      if ((fpcheck = fopen(szOutputFileName, "r")) != NULL)
      {
         bSearchSpectrum = 0;
         fclose(fpcheck);
      }
   }

   return bSearchSpectrum;
}


bool mango_preprocess::AdjustMassTol(struct Query *pScoring)
{
   if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
   {
      pScoring->_pepMassInfo.dPeptideMassTolerance = g_staticParams.tolerances.dTolerancePeptide;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         pScoring->_pepMassInfo.dPeptideMassTolerance *= pScoring->_spectrumInfoInternal.iChargeState;
      }
   }
   else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
   {
      pScoring->_pepMassInfo.dPeptideMassTolerance = g_staticParams.tolerances.dTolerancePeptide * 0.001;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         pScoring->_pepMassInfo.dPeptideMassTolerance *= pScoring->_spectrumInfoInternal.iChargeState;
      }
   }
   else // ppm
   {
      pScoring->_pepMassInfo.dPeptideMassTolerance = g_staticParams.tolerances.dTolerancePeptide
         * pScoring->_pepMassInfo.dExpPepMass / 1000000.0;
   }

   if (g_staticParams.tolerances.iIsotopeError == 0)
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 1) // search 0, +1 isotope windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - C13_DIFF * PROTON_MASS;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 2) // search 0, +1, +2 isotope windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - 2.0 * C13_DIFF * PROTON_MASS;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 3) // search 0, +1, +2, +3 isotope windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - 3.0 * C13_DIFF * PROTON_MASS;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 4) // search -1 to +1, +2, +3 isotope windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - 3.0 * C13_DIFF * PROTON_MASS;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance + 1.0 * C13_DIFF * PROTON_MASS;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 5) // search -8, -4, 0, 4, 8 windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - 8.1;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance + 8.1;
   }
   else  // Should not get here.
   {
      fprintf(stderr,  " Error - iIsotopeError=%d\n",  g_staticParams.tolerances.iIsotopeError);
      return false;
   }

   if (g_staticParams.vectorMassOffsets.size() > 0)
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus -= g_staticParams.vectorMassOffsets[g_staticParams.vectorMassOffsets.size()-1];

   if (pScoring->_pepMassInfo.dPeptideMassTolerancePlus > g_staticParams.options.dPeptideMassHigh)
      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = g_staticParams.options.dPeptideMassHigh;

   if (pScoring->_pepMassInfo.dPeptideMassToleranceMinus < g_staticParams.options.dPeptideMassLow)
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = g_staticParams.options.dPeptideMassLow;

   if (pScoring->_pepMassInfo.dPeptideMassToleranceMinus < 100.0)
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = 100.0;

   return true;
}


//  Reads MSMS data file as ASCII mass/intensity pairs.
bool mango_preprocess::LoadIons(struct Query *pScoring,
                               double *pdTmpRawData,
                               Spectrum mstSpectrum,
                               struct PreprocessStruct *pPre)
{
   int  i;
   double dIon,
          dIntensity;

   i = 0;
   while(true)
   {
      if (i >= mstSpectrum.size())
         break;

      dIon = mstSpectrum.at(i).mz;
      dIntensity = mstSpectrum.at(i).intensity;
      i++;


      pScoring->_spectrumInfoInternal.dTotalIntensity += dIntensity;

      if ((dIntensity >= g_staticParams.options.dMinIntensity) && (dIntensity > 0.0))
      {
         if (dIon < (pScoring->_pepMassInfo.dExpPepMass + 50.0))
         {
            int iBinIon = BIN(dIon);

            dIntensity = sqrt(dIntensity);

            if (iBinIon > pPre->iHighestIon)
               pPre->iHighestIon = iBinIon;

            if ((iBinIon < pScoring->_spectrumInfoInternal.iArraySize) && (dIntensity > pdTmpRawData[iBinIon]))
            {
               if (     !(dIon > 310.16-0.1 && dIon < 310.16+0.1)        // clear out contaminant peaks
                     && !(dIon > 311.16-0.1 && dIon < 311.16+0.1)
                     && !(dIon > 430.21-0.1 && dIon < 430.21+0.1)
                     && !(dIon > 431.22-0.1 && dIon < 431.22+0.1)
                     && !(dIon > 655.36-0.1 && dIon < 655.36+0.1)
                     && !(dIon > 677.38-0.1 && dIon < 677.38+0.1)
                     && !(dIon > g_staticParams.options.dReporterMass+PROTON_MASS-0.1    // also remove 752 reporter ion and isotopes
                       && dIon < g_staticParams.options.dReporterMass+PROTON_MASS+0.1)
                     && !(dIon > g_staticParams.options.dReporterMass+2*PROTON_MASS-0.1
                       && dIon < g_staticParams.options.dReporterMass+2*PROTON_MASS+0.1)
                     && !(dIon > g_staticParams.options.dReporterMass+3*PROTON_MASS-0.1
                       && dIon < g_staticParams.options.dReporterMass+3*PROTON_MASS+0.1))
               {
                  if (dIntensity > pdTmpRawData[iBinIon])
                     pdTmpRawData[iBinIon] = dIntensity;

                  if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity)
                     pPre->dHighestIntensity = pdTmpRawData[iBinIon];
               }

            }
         }
      }
   }

   return true;
}


// pdTmpRawData now holds raw data, pdTmpCorrelationData is windowed data after this function
void mango_preprocess::MakeCorrData(double *pdTmpRawData,
                                   double *pdTmpCorrelationData,
                                   struct Query *pScoring,
                                   struct PreprocessStruct *pPre)
{
   int  i,
        ii,
        iBin,
        iWindowSize,
        iNumWindows=10;
   double dMaxWindowInten,
          dTmp1,
          dTmp2;

   iWindowSize = (int)((pPre->iHighestIon)/iNumWindows) + 1;

   for (i=0; i<iNumWindows; i++)
   {
      dMaxWindowInten = 0.0;

      for (ii=0; ii<iWindowSize; ii++)    // Find max inten. in window.
      {
         iBin = i*iWindowSize+ii;
         if (iBin < pScoring->_spectrumInfoInternal.iArraySize)
         {
            if (pdTmpRawData[iBin] > dMaxWindowInten)
               dMaxWindowInten = pdTmpRawData[iBin];
         }
      }

      if (dMaxWindowInten > 0.0)
      {
         dTmp1 = 50.0 / dMaxWindowInten;
         dTmp2 = 0.05 * pPre->dHighestIntensity;

         for (ii=0; ii<iWindowSize; ii++)    // Normalize to max inten. in window.
         {
            iBin = i*iWindowSize+ii;
            if (iBin < pScoring->_spectrumInfoInternal.iArraySize)
            {
               if (pdTmpRawData[iBin] > dTmp2)
                  pdTmpCorrelationData[iBin] = pdTmpRawData[iBin]*dTmp1;
            }
         }
      }
   }
}


//MH: This function allocates memory to be shared by threads for spectral processing
bool mango_preprocess::AllocateMemory(int maxNumThreads)
{
   int i;

   //MH: Find appropriately sized array cushion based on user parameters. Fixes error found by Patrick Pedrioli for
   // very wide mass tolerance searches (i.e. 500 Da).
   double dCushion = 0.0;
   if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
   {
      dCushion = g_staticParams.tolerances.dTolerancePeptide;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
        dCushion *= 8; //MH: hope +8 is large enough charge because g_staticParams.options.iEndCharge can be overridden.
      }
   }
   else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
   {
      dCushion = g_staticParams.tolerances.dTolerancePeptide * 0.001;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         dCushion *= 8; //MH: hope +8 is large enough charge because g_staticParams.options.iEndCharge can be overridden.
      }
   }
   else // ppm
   {
      dCushion = g_staticParams.tolerances.dTolerancePeptide * g_staticParams.options.dPeptideMassHigh / 1000000.0;
   }

   //MH: Must be equal to largest possible array
   int iArraySize = (int)((g_staticParams.options.dPeptideMassHigh + dCushion + 2.0) * g_staticParams.dInverseBinWidth);

   //MH: Initally mark all arrays as available (i.e. false=not inuse).
   pbMemoryPool = new bool[maxNumThreads];
   for (i=0; i<maxNumThreads; i++)
   {
      pbMemoryPool[i] = false;
   }

   //MH: Allocate arrays
   ppdTmpRawDataArr = new double*[maxNumThreads]();
   for (i=0; i<maxNumThreads; i++)
   {
      try
      {
         ppdTmpRawDataArr[i] = new double[iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pdTmpRawData[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         fprintf(stderr, "Mango ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   //MH: Allocate arrays
   ppdTmpFastXcorrDataArr = new double*[maxNumThreads]();
   for (i=0; i<maxNumThreads; i++)
   {
      try
      {
         ppdTmpFastXcorrDataArr[i] = new double[iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pdTmpFastXcorrData[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         fprintf(stderr, "Mango ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   //MH: Allocate arrays
   ppdTmpCorrelationDataArr = new double*[maxNumThreads]();
   for (i=0; i<maxNumThreads; i++)
   {
      try
      {
         ppdTmpCorrelationDataArr[i] = new double[iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pdTmpCorrelationData[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         fprintf(stderr, "Mango ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   return true;
}


//MH: Deallocates memory shared by threads during spectral processing.
bool mango_preprocess::DeallocateMemory(int maxNumThreads)
{
   int i;

   delete[] pbMemoryPool;

   for (i=0; i<maxNumThreads; i++)
   {
      delete[] ppdTmpRawDataArr[i];
      delete[] ppdTmpFastXcorrDataArr[i];
      delete[] ppdTmpCorrelationDataArr[i];
   }

   delete[] ppdTmpRawDataArr;
   delete[] ppdTmpFastXcorrDataArr;
   delete[] ppdTmpCorrelationDataArr;

   return true;
}

bool mango_preprocess::IsValidInputType(int inputType)
{
   return (inputType == InputType_MZXML || inputType == InputType_RAW);
}
