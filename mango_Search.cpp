/*
   Copyright 2016 University of Washington

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "Common.h"
#include "mango.h"
#include "mango_Search.h"
#include "mango_DataInternal.h"
#include "mango_Preprocess.h"
#include "CometDecoys.h"

// Generate data for both sp scoring (pfSpScoreData) and xcorr analysis (FastXcorr).
mango_Search::mango_Search()
{
}
#define NUMPEPTIDES 10

mango_Search::~mango_Search()
{
}

void insert_pep_pq(char *pepArray[], char *proArray[], float xcorrArray[], char *ins_pep, char *ins_pro, float ins_xcorr)
{
   int i; 

   // check for duplicates
   for (i = 0; i< NUMPEPTIDES - 1; i++) {
      if (pepArray[i] != NULL) {
         if (!strcmp(pepArray[i], ins_pep)) {
             //cout << "" << pepArray[i] << " " << ins_pep << " Found duplicate" << endl;
             return;
         }
      }
   }

   // Insert first into the array
   if (xcorrArray[NUMPEPTIDES - 1] < ins_xcorr) {
      xcorrArray[NUMPEPTIDES - 1] = ins_xcorr;

      if (pepArray[NUMPEPTIDES - 1]) delete [] pepArray[NUMPEPTIDES - 1];
      pepArray[NUMPEPTIDES - 1] = ins_pep;

      if (proArray[NUMPEPTIDES - 1]) delete [] proArray[NUMPEPTIDES - 1];
      proArray[NUMPEPTIDES - 1] = ins_pro;
   } else return;

   // shiffle
   for (i = NUMPEPTIDES - 1; i > 0; i--) {
      if (xcorrArray[i] > xcorrArray[i-1]) {
         // Swap
         float temp = xcorrArray[i];
         char *temp_pep = pepArray[i];
         char *temp_pro = proArray[i];
         xcorrArray[i] = xcorrArray[i-1];
         pepArray[i] = pepArray[i-1];
         proArray[i] = proArray[i-1];
         xcorrArray[i-1] = temp;
         pepArray[i-1] = temp_pep;
         proArray[i-1] = temp_pro;
      } else break;
   }
}

#define HISTOGRAM_BIN_SIZE 0.1
#define MAX_XCORR_VALUE 20
#define NUM_BINS (int)(MAX_XCORR_VALUE/HISTOGRAM_BIN_SIZE + 1)

inline int mango_get_histogram_bin_num(float value)
{
    if (value > MAX_XCORR_VALUE)
       value = MAX_XCORR_VALUE;
    else if (value < 0)
       value = 0;
    return value/HISTOGRAM_BIN_SIZE;
}

void mango_print_histogram(int hist_pep[])
{
   for (int i = 0; i <NUM_BINS; i++)
      printf("%d ", hist_pep[i]);
   printf("\n");
}

#define DECOY_SIZE 3000
#define MAX_DECOY_PEP_LEN 40

bool mango_Search::CalculateEValue(int *hist_pep,
                                    int iMatchPepCount,
                                    double *dSlope,
                                    double *dIntercept,
                                    double dNeutralPepMass,
                                    int iScanNumber)
{
   int iMaxCorr;
   int iStartCorr;
   int iNextCorr;

   if (iMatchPepCount < DECOY_SIZE)
   {
      if (!GenerateXcorrDecoys(dNeutralPepMass, iMatchPepCount, hist_pep, iScanNumber))
      {
         return false;
      }
   }

   LinearRegression(hist_pep, dSlope, dIntercept, &iMaxCorr, &iStartCorr, &iNextCorr);
   *dSlope *= 10.0; // Used in pow() function so do multiply outside of for loop.

   return true;
}


void mango_Search::LinearRegression(int *piHistogram,
                                     double *slope,
                                     double *intercept,
                                     int *iMaxXcorr,
                                     int *iStartXcorr,
                                     int *iNextXcorr)
{
   double Sx, Sxy;      // Sum of square distances.
   double Mx, My;       // means
   double b, a;
   double SumX, SumY;   // Sum of X and Y values to calculate mean.

   double dCummulative[HISTO_SIZE];  // Cummulative frequency at each xcorr value.

   int i;
   int iNextCorr;    // 2nd best xcorr index
   int iMaxCorr=0;   // max xcorr index
   int iStartCorr;
   int iNumPoints;

   // Find maximum correlation score index.
   for (i=HISTO_SIZE-2; i>=0; i--)
   {
      if (piHistogram[i] > 0)
         break;
   }
   iMaxCorr = i;

   iNextCorr = 0;
   for (i=0; i<iMaxCorr; i++)
   {
      if (piHistogram[i]==0)
      {
         // register iNextCorr if there's a histo value of 0 consecutively
         if (piHistogram[i+1]==0 || i+1 == iMaxCorr)
         {
            if (i>0)
               iNextCorr = i-1;
            break;
         }
      }
   }

   if (i==iMaxCorr)
   {
      iNextCorr = iMaxCorr;
      if (iMaxCorr>12)
         iNextCorr = iMaxCorr-2;
   }


   // Create cummulative distribution function from iNextCorr down, skipping the outliers.
   dCummulative[iNextCorr] = piHistogram[iNextCorr];
   for (i=iNextCorr-1; i>=0; i--)
   {
      dCummulative[i] = dCummulative[i+1] + piHistogram[i];
      if (piHistogram[i+1] == 0)
         dCummulative[i+1] = 0.0;
   }

   // log10
   for (i=iNextCorr; i>=0; i--)
   {
      piHistogram[i] = (int)dCummulative[i];  // First store cummulative in histogram.
      dCummulative[i] = log10(dCummulative[i]);
   }

   iStartCorr = 1;
   if (iNextCorr >= 30)
      iStartCorr = (int)(iNextCorr - iNextCorr*0.25);
   else if (iNextCorr >= 15)
      iStartCorr = (int)(iNextCorr - iNextCorr*0.5);

   Mx=My=a=b=0.0;

   while (iStartCorr >= 0)
   {
      Sx=Sxy=SumX=SumY=0.0;
      iNumPoints=0;

      // Calculate means.
      for (i=iStartCorr; i<=iNextCorr; i++)
      {
         if (piHistogram[i] > 0)
         {
            SumY += (float)dCummulative[i];
            SumX += i;
            iNumPoints++;
         }
      }

      if (iNumPoints > 0)
      {
         Mx = SumX / iNumPoints;
         My = SumY / iNumPoints;
      }
      else
         Mx = My = 0.0;

      // Calculate sum of squares.
      for (i=iStartCorr; i<=iNextCorr; i++)
      {
         if (dCummulative[i] > 0)
         {
            double dX;
            double dY;

            dX = i - Mx;
            dY = dCummulative[i] - My;

            Sx  += dX*dX;
            Sxy += dX*dY;
         }
      }

      if (Sx > 0)
         b = Sxy / Sx;   // slope
      else
         b = 0;

      if (b < 0.0)
         break;
      else
         iStartCorr--;
   }

   a = My - b*Mx;  // y-intercept

   *slope = b;
   *intercept = a;
   *iMaxXcorr = iMaxCorr;
   *iStartXcorr = iStartCorr;
   *iNextXcorr = iNextCorr;
}


// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
bool mango_Search::GenerateXcorrDecoys(double dNeutralPepMass,
                                        int iMatchPepCount,
                                        int *hist_pep,
                                        int iScanNumber)
{
   int i;
   int ii;
   int j;
   int bin_num;
   int iMaxFragCharge;
   int ctCharge;
   double dBion;
   double dYion;
   double dFastXcorr;
   double dFragmentIonMass;

   int *piHistogram;

   int iFragmentIonMass;
   int iWhichQuery;

   for (iWhichQuery=0; iWhichQuery<(int)g_pvQuery.size(); iWhichQuery++)
   {
      if (g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iScanNumber == iScanNumber)
         break;
   }
   if (iWhichQuery < (int)g_pvQuery.size() && g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iScanNumber == iScanNumber)
   {
      Query* pQuery = g_pvQuery.at(iWhichQuery);

      piHistogram = hist_pep;

      //iMaxFragCharge = pQuery->_spectrumInfoInternal.iMaxFragCharge;
      iMaxFragCharge = 1;  //FIX only considering 1+ charges now

      // DECOY_SIZE is the minimum # of decoys required or else this function is
      // called.  So need generate iLoopMax more xcorr scores for the histogram.
      int iLoopMax = DECOY_SIZE - iMatchPepCount;
      int iLastEntry;

      iLastEntry = iMatchPepCount;

      if (iLastEntry > g_staticParams.options.iNumStored)
         iLastEntry = g_staticParams.options.iNumStored;

      j=0;
      for (i=0; i<iLoopMax; i++)  // iterate through required # decoys
      {
         dFastXcorr = 0.0;

         for (j=0; j<MAX_DECOY_PEP_LEN; j++)  // iterate through decoy fragment ions
         {
            dBion = decoyIons[i].pdIonsN[j];
            dYion = decoyIons[i].pdIonsC[j];

            for (ii=0; ii<2; ii++)
            {
               dFragmentIonMass =  0.0;
               switch (ii)
               {
                  case 0:
                     dFragmentIonMass = dBion;
                     break;
                  case 1:
                     dFragmentIonMass = dYion;
                     break;
               }

               for (ctCharge=1; ctCharge<=iMaxFragCharge; ctCharge++)
               {
                  dFragmentIonMass = (dFragmentIonMass + (ctCharge-1)*PROTON_MASS)/ctCharge;

                  if (dFragmentIonMass < dNeutralPepMass)
                  {
                     iFragmentIonMass = BIN(dFragmentIonMass);

                     if (iFragmentIonMass < pQuery->_spectrumInfoInternal.iArraySize && iFragmentIonMass >= 0)
                     {
                        int x = iFragmentIonMass / SPARSE_MATRIX_SIZE;
                        if (pQuery->ppfSparseFastXcorrData[x]!=NULL)
                        {
                           int y = iFragmentIonMass - (x*SPARSE_MATRIX_SIZE);
                           dFastXcorr += pQuery->ppfSparseFastXcorrData[x][y];
                        }
                     }
                     else
                     {
                        char szErrorMsg[256];
                        sprintf(szErrorMsg,  " Error - XCORR DECOY: dFragMass %f, iFragMass %d, ArraySize %d, InputMass %f, scan %d, z %d",
                              dFragmentIonMass,
                              iFragmentIonMass,
                              pQuery->_spectrumInfoInternal.iArraySize,
                              pQuery->_pepMassInfo.dExpPepMass,
                              pQuery->_spectrumInfoInternal.iScanNumber,
                              ctCharge);

                        string strErrorMsg(szErrorMsg);
                        logerr(szErrorMsg);
                        return false;
                     }
                  }

               }
            }
         }

         dFastXcorr *= 0.005;
         bin_num = mango_get_histogram_bin_num(dFastXcorr);
         piHistogram[bin_num] += 1;
      }

      return true;
   }

   return false;

}


void mango_Search::SearchForPeptides(char *szMZXML,
                                      const char *protein_file,
                                      enzyme_cut_params params,
                                      const char *pep_hash_file)
{
   int i;
   int ii;
   double dTolerance;
   double dPPM = 20.0;  // use 20ppm tolerance for now
   int hist_pep1[NUM_BINS],
       hist_pep2[NUM_BINS],
       hist_combined[NUM_BINS],
       num_pep1,
       num_pep2,
       num_pep_combined;

   char szOutputTxt[SIZE_FILE];

   char *toppep1[NUMPEPTIDES], *toppep2[NUMPEPTIDES], *toppepcombined[NUMPEPTIDES];
   char *toppro1[NUMPEPTIDES], *toppro2[NUMPEPTIDES], *topprocombined[NUMPEPTIDES];
   float xcorrPep1[NUMPEPTIDES], xcorrPep2[NUMPEPTIDES], xcorrCombined[NUMPEPTIDES];

   for (ii = 0; ii < NUMPEPTIDES; ii++)
   {
      toppep1[ii] = toppep2[ii] = toppepcombined[ii] = NULL;
      toppro1[ii] = toppro2[ii] = topprocombined[ii] = NULL;
   }

   strcpy(szOutputTxt, szMZXML);
   szOutputTxt[strlen(szOutputTxt)-5]='\0';
   strcat(szOutputTxt, "txt");
   FILE *fptxt;
   if ((fptxt=fopen(szOutputTxt, "w")) == NULL)
   {
      printf(" Error - cannot write txt output %s\n", szOutputTxt);
      return;
   }

   MSReader mstReader;
   Spectrum mstSpectrum;
   // We want to read only MS2 scans.
   vector<MSSpectrumType> msLevel;
   msLevel.push_back(MS2);

   mstReader.setFilter(msLevel);
   mstReader.readFile(szMZXML, mstSpectrum, 1);
 
   mango_preprocess::AllocateMemory(1);

   // If PeptideHash not present, generate it now; otherwise open the hash file.
   protein_hash_db_t phdp = phd_retrieve_hash_db(protein_file, params, pep_hash_file);

   fprintf(fptxt, "scan\texp_mass1\texp_mass2\tpeptide1\txcorr1\tevalue1\tcalcmass1\tpeptide2\txcorr2\tevalue2\tcalcmass2\tcombinedxcorr\tcombinedevalue\n"); 

   FILE *fpxml;
   char szOutput[1024];
   char szBaseName[1024];
   char *combinedPep;
   char *szPeptide;
   char *szProtein;

   strcpy(szBaseName, szMZXML);
   if (!strcmp(szBaseName+strlen(szBaseName)-6, ".mzXML"))
      szBaseName[strlen(szBaseName)-6]='\0';
   else if (!strcmp(szBaseName+strlen(szBaseName)-6, ".mzML"))
      szBaseName[strlen(szBaseName)-5]='\0';
   sprintf(szOutput, "%s.pep.xml", szBaseName);

   if ((fpxml=fopen(szOutput, "w")) == NULL)
   {
      printf(" Error - cannot write pepXML output %s\n", szOutput);
      return;
   }
   int iIndex=0;

   WritePepXMLHeader(fpxml, szBaseName, protein_file);

   if (!g_staticParams.options.bVerboseOutput)
   {
      printf("percent complete: ");
      fflush(stdout);
   }

/*
   int iHistoSize[20];
   int iNonEmptyScan = 0;
   int iEmptyScan = 0;
   for (int x=0; x<=11; x++)
      iHistoSize[x]=0;
   for (i=0; i<(int)pvSpectrumList.size(); i++)
   {
      int x = (int)pvSpectrumList.at(i).pvdPrecursors.size();

      if (x>=20)
         x=11;
      else if (x>=10)
         x=10;

      iHistoSize[x] += 1;
      
      if (x>0)
         iNonEmptyScan++;
      else
         iEmptyScan++;
   }
   printf("\n");
   for (int x=0; x<=11; x++)
      printf("%d  %d\n", x, iHistoSize[x]);
   printf("iNonEmptyScans = %d\n", iNonEmptyScan);
   printf("iEmptyScans = %d\n", iEmptyScan);
   exit(1);
*/



//g_staticParams.options.bVerboseOutput = true;

   bool bSilac = false;
   int y=0;

/*
printf("MS1_scan\tMS2_scan\tintact_precursor_neutral_mass\tintact_precursor_mz\tintact_precursor_charge\treleased_neutral1\trelease_neutral2\treleased_mz1\treleased_mz2\trelease_charge1\treleased_charge2\n");
*/

   for (i=0; i<(int)pvSpectrumList.size(); i++)
   {

/*
for (ii=0; ii<(int)pvSpectrumList.at(i).pvdPrecursors.size(); ii++)
{
   double dMZ1 =  (pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1
                  + pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1 * PROTON_MASS)/ pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1;
   double dMZ2 =  (pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2
                  + pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2 * PROTON_MASS)/ pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2;

   printf("%d\t%d\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%d\t%d\n",
         pvSpectrumList.at(i).iPrecursorScanNumber,
         pvSpectrumList.at(i).iScanNumber,
         pvSpectrumList.at(i).dPrecursorMZ * pvSpectrumList.at(i).iPrecursorCharge - (pvSpectrumList.at(i).iPrecursorCharge * PROTON_MASS),
         pvSpectrumList.at(i).dPrecursorMZ,
         pvSpectrumList.at(i).iPrecursorCharge,
         pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1,
         pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2,
         dMZ1,
         dMZ2,
         pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1,
         pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2);

}
*/

      if (1) //pvSpectrumList.at(i).iScanNumber >=35475 && pvSpectrumList.at(i).iScanNumber<=35475) // limit analysis range during dev/testing
      {
         mango_preprocess::LoadAndPreprocessSpectra(mstReader, pvSpectrumList.at(i).iScanNumber, 0, 0);

         for (ii=0; ii<(int)pvSpectrumList.at(i).pvdPrecursors.size(); ii++)
         {
            for (int j = 0; j < NUM_BINS; j++)
               hist_pep1[j] = hist_pep2[j] = hist_combined[j] = 0;

            num_pep1 = num_pep2 = num_pep_combined = 0;

            double dMZ1 =  (pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1
                  + pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1 * PROTON_MASS)/ pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1;
            double dMZ2 =  (pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2
                  + pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2 * PROTON_MASS)/ pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2;

            for (int li = 0; li < NUMPEPTIDES; li++)
            {
               xcorrPep1[li] = xcorrPep2[li] = xcorrCombined[li] = -99999;
               toppep1[li] = toppep2[li] = toppepcombined[li] = NULL;
            }

            if (g_staticParams.options.bVerboseOutput)
            {
               printf("Scan %d (i=%d), retrieving peptides of mass %0.4f (%d+ %0.4f) and %0.4f (%d+ %0.4f)\n",
                     pvSpectrumList.at(i).iScanNumber,
                     i,
                     pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1,
                     pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1,
                     dMZ1,
                     pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2,
                     pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2,
                     dMZ2);
            }

            double pep_mass1 = pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 - LYSINE_MOD - g_staticParams.precalcMasses.dOH2;
            double pep_mass2 = pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 - LYSINE_MOD - g_staticParams.precalcMasses.dOH2;

            if (pep_mass1 <= 0)
            {
               cout << "Peptide mass1 is coming out to be zero after removing Lysine residue" << endl;
               exit(1);
            }

            if (pep_mass2 <= 0)
            {
               cout << "Peptide mass2 is coming out to be zero after removing Lysine residue" << endl;
               exit(1);
            }

            fprintf(fptxt, "%d\t%f\t%f", pvSpectrumList.at(i).iScanNumber,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2);

            double dXcorr = 0.0;
            vector<double> vdXcorr_pep1;  // store xcorr scores to be used in combined histogram
            vector<double> vdXcorr_pep2;

            if (g_staticParams.options.bVerboseOutput)
            {
               cout << "After Lysine residue reduction the peptide of mass " << pep_mass1 << " are being extracted";
               cout << " (" << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 << ")" << endl;
            }
            dTolerance = (dPPM * pep_mass1) / 1e6;

            double dSilacMass = 0.0;

  
            // bSilac
            for (y=0; y<2; y++)
            {
               if (y==0)
                  dSilacMass = 8.014199 + 8.014199;  //...K...K
               else
                  dSilacMass = 6.020129 + 8.014199;  //...K...R

               if (!bSilac)
               {
                  y=2;   // break out of y for-loop
                  dSilacMass = 0.0;
               }

               for (int x=0; x<2; x++)
               {
                  vector<peptide_hash_database::phd_peptide> *peptides1 = phdp->phd_get_peptides_ofmass_tolerance(pep_mass1 - dSilacMass - x*1.003355, dTolerance);
      
                  for (peptide_hash_database::phd_peptide peptide : *peptides1)
                  {
                     szPeptide = new char[peptide.phdpep_sequence().length() + 1];
                     strcpy(szPeptide, (peptide.phdpep_sequence()).c_str() );
      
                     szProtein = new char[peptide.phdpep_protein_list(0).phdpro_name().length() + 1];
                     strcpy(szProtein, (peptide.phdpep_protein_list(0).phdpro_name()).c_str() );
      
                     // sanity check to ignore peptides w/unknown AA residues
                     // should not be needed now that this is addressed in the hash building
                     if (strchr(szPeptide, 'B') || strchr(szPeptide, 'X') || strchr(szPeptide, 'J') || strchr(szPeptide, 'Z'))
                        dXcorr = 0.0;
                     else
                     {
      
                        if (bSilac)
                        {
                           if ((y==0 && szPeptide[strlen(szPeptide)-1]=='K') || (y==1 && szPeptide[strlen(szPeptide)-1]=='R')) // SILAC
                              dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber, true);
                        }
                        else
                           dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber, false);
      
                     }
      
                     vdXcorr_pep1.push_back(dXcorr);
      
                     hist_pep1[mango_get_histogram_bin_num(dXcorr)]++;
                     insert_pep_pq(toppep1, toppro1, xcorrPep1, szPeptide, szProtein, dXcorr);
                     num_pep1++;
                     if (g_staticParams.options.bVerboseOutput)
                        cout << "pep1: " << szPeptide << "  xcorr " << dXcorr << "  protein " << szProtein << endl;
                  }
      
                  (*peptides1).clear();
               }
            } // SILAC

            if (g_staticParams.options.bVerboseOutput)
            {
               cout << "After Lysine residue reduction the peptide of mass " << pep_mass2 << " are being extracted";
               cout << " (" << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 << ")" << endl;
            }
            dTolerance = (dPPM * pep_mass2) / 1e6;

            // bSilac
            for (y=0; y<2; y++)
            {
               if (y==0)
                  dSilacMass = 8.014199 + 8.014199;  //...K...K
               else
                  dSilacMass = 6.020129 + 8.014199;  //...K...R

               if (!bSilac)
               {
                  y=2;   // break out of y for-loop
                  dSilacMass = 0.0;
               }

               for (int x=0; x<2; x++)
               {
                  vector<peptide_hash_database::phd_peptide> *peptides2 = phdp->phd_get_peptides_ofmass_tolerance(pep_mass2 - dSilacMass - x*1.003355, dTolerance);
      
                  for (peptide_hash_database::phd_peptide peptide : *peptides2)
                  {
                     szPeptide = new char[peptide.phdpep_sequence().length() + 1];
                     strcpy(szPeptide, (peptide.phdpep_sequence()).c_str() );
      
                     szProtein = new char[peptide.phdpep_protein_list(0).phdpro_name().length() + 1];
                     strcpy(szProtein, (peptide.phdpep_protein_list(0).phdpro_name()).c_str() );
      
                     // sanity check to ignore peptides w/unknown AA residues
                     // should not be needed now that this is addressed in the hash building
                     if (strchr(szPeptide, 'B') || strchr(szPeptide, 'X') || strchr(szPeptide, 'J') || strchr(szPeptide, 'Z'))
                        dXcorr = 0.0;
                     else
                     {
                        if (bSilac)
                        {
                           if ((y==0 && szPeptide[strlen(szPeptide)-1]=='K') || (y==1 && szPeptide[strlen(szPeptide)-1]=='R')) //SILAC
                             dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber, true);
                        }
                        else
                           dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber, false);
                     }
      
                     vdXcorr_pep2.push_back(dXcorr);
      
                     hist_pep2[mango_get_histogram_bin_num(dXcorr)]++;
                     insert_pep_pq(toppep2, toppro2, xcorrPep2, szPeptide, szProtein, dXcorr);
                     num_pep2++;
                     if (g_staticParams.options.bVerboseOutput)
                        cout << "pep2: " << szPeptide << "  xcorr " << dXcorr << "  protein " << szProtein << endl;
                  }
                  (*peptides2).clear();
               }
            } //SILAC

            if (toppep1[0] == NULL || toppep2[0] == NULL)
               continue;

            double dSlope;
            double dIntercept;
            double dExpect;

            // return dSlope and dIntercept for histogram
            CalculateEValue(hist_pep1, num_pep1, &dSlope, &dIntercept,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1, pvSpectrumList.at(i).iScanNumber);

            if (g_staticParams.options.bVerboseOutput)
            {
               mango_print_histogram(hist_pep1);
               cout << "Top "<< NUMPEPTIDES << " pep1 peptides for this scan are " << endl;
            }

            double dExpect1 = 999;;
            for (int li = 0 ; li < NUMPEPTIDES; li++)
            {
               if (toppep1[li] != NULL)
               {
                  if (dSlope > 0)
                     dExpect = 999;
                  else
                     dExpect = pow(10.0, dSlope * xcorrPep1[li] + dIntercept);

                  if (li == 0)
                     dExpect1 = dExpect;

                  if (g_staticParams.options.bVerboseOutput)
                     cout << "pep1_top: " << toppep1[li] << " xcorr " << xcorrPep1[li] << " expect " << dExpect << endl;
               }
            }

            if (toppep1[0] != NULL)
               fprintf(fptxt, "\t%s\t%f\t%0.3E\t%f", toppep1[0], xcorrPep1[0], dExpect1, phdp->phd_calculate_mass_peptide(string(toppep1[0])));
            else
               fprintf(fptxt, "\t-\t0\t999\t0");

            CalculateEValue(hist_pep2, num_pep2, &dSlope, &dIntercept,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2, pvSpectrumList.at(i).iScanNumber);

            if (g_staticParams.options.bVerboseOutput)
            {
               mango_print_histogram(hist_pep2);
               cout << "Top "<< NUMPEPTIDES << " pep2 peptides for this scan are " << endl;
            }
            double dExpect2 = 999;
            for (int li = 0; li < NUMPEPTIDES; li++)
            {
               if (toppep2[li] != NULL)
               {
                  if (dSlope > 0)
                     dExpect = 999;
                  else
                     dExpect = pow(10.0, dSlope * xcorrPep2[li] + dIntercept);

                  if (li == 0)
                     dExpect2 = dExpect;

                  if (g_staticParams.options.bVerboseOutput)
                     cout << "pep2_top: " << toppep2[li] << " xcorr " << xcorrPep2[li] << " expect " << dExpect << endl;
               }
            }

            if (toppep2[0] != NULL)
               fprintf(fptxt, "\t%s\t%f\t%0.3E\t%f", toppep2[0], xcorrPep2[0], dExpect2, phdp->phd_calculate_mass_peptide(string(toppep2[0])));
            else
               fprintf(fptxt, "\t-\t0\t999\t0");

            if (g_staticParams.options.bVerboseOutput)
               cout << "Size of peptide1 list is " << num_pep1 << " and the size of peptide2 list is " << num_pep2 << endl;

            // Compute histogram of combined scores;
            for (int x=0; x<NUM_BINS; x++)
               hist_combined[x] = 0;
    
            for (vector<double>::iterator x = vdXcorr_pep1.begin(); x != vdXcorr_pep1.end(); ++x)
            {
               for (vector<double>::iterator y = vdXcorr_pep2.begin(); y != vdXcorr_pep2.end(); ++y)
               {
                  hist_combined[mango_get_histogram_bin_num(*x + *y)]++;
               }
            }

            CalculateEValue(hist_combined, num_pep_combined, &dSlope, &dIntercept,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2, pvSpectrumList.at(i).iScanNumber);

            if (g_staticParams.options.bVerboseOutput)
               mango_print_histogram(hist_combined);

            // take all combinations of top pep1 and pep2 and store best
            for (int x = 0; x< NUMPEPTIDES - 1; x++)
            {
               if (toppep1[x] != NULL)
               {
                  for (int y = 0; y< NUMPEPTIDES - 1; y++)
                  {
                     if (toppep2[y] != NULL)
                     {
                        combinedPep = new char[strlen(toppep1[x]) + strlen(toppep2[y]) + 4];
                        sprintf(combinedPep, "%s + %s", toppep1[x], toppep2[y]);

                        double dCombinedXcorr = xcorrPep1[x] + xcorrPep2[y];

                        insert_pep_pq(toppepcombined, topprocombined, xcorrCombined, combinedPep, NULL, dCombinedXcorr);
                     }
                  }
               }
            }

            double dExpectCombined = 999;
            for (int li = 0; li < NUMPEPTIDES; li++)
            {
               if (toppepcombined[li] != NULL)
               {
                  if (dSlope > 0)
                     dExpect = 999;
                  else
                     dExpect = pow(10.0, dSlope * xcorrCombined[li] + dIntercept);

                  if (g_staticParams.options.bVerboseOutput)
                     cout << "combined: " << toppepcombined[li] << " xcorr " << xcorrCombined[li] << " expect " << dExpect << endl;

                  if (li == 0)
                  {
                     fprintf(fptxt, "\t%f\t%0.3E\n",  xcorrCombined[li], dExpect);
                     dExpectCombined = dExpect;
                  }
               }
            }

            int iCharge = (pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1>pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2
                  ? pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1
                  : pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2);

            WriteSpectrumQuery(fpxml, szBaseName,
                  pep_mass1, pep_mass2,
                  xcorrPep1[0], xcorrPep2[0],
                  dExpect1, dExpect2,
                  phdp->phd_calculate_mass_peptide(string(toppep1[0])), phdp->phd_calculate_mass_peptide(string(toppep2[0])),
                  xcorrCombined[0], dExpectCombined,
                  toppep1[0], toppep2[0],
                  toppro1[0], toppro2[0],
   //             pvSpectrumList.at(i).iPrecursorCharge,
                  iCharge,                                     // report largest charge of the two released peptides
                  iIndex, pvSpectrumList.at(i).iScanNumber);
         }

         for (int y=0; y<(int)g_pvQuery.size(); y++)
         {
            // need to free processed spectrum data here
            for (int x=0;x<g_pvQuery.at(y)->iFastXcorrData;x++)
            {
               if (g_pvQuery.at(y)->ppfSparseFastXcorrData[x] != NULL)
                  delete[] g_pvQuery.at(y)->ppfSparseFastXcorrData[x];
            }
         }

         g_pvQuery.clear();

         if (!g_staticParams.options.bVerboseOutput)
         {
            printf("%5.1f%%", (float)(100.0*i/pvSpectrumList.size()));
            fflush(stdout);
            printf("\b\b\b\b\b\b");
         }
      } //scan range restriction
   }

   mango_preprocess::DeallocateMemory(1);

   fprintf(fpxml, "  </msms_run_summary>\n");
   fprintf(fpxml, "</msms_pipeline_analysis>\n");

   if (!g_staticParams.options.bVerboseOutput)
     printf("\n\nDone.\n");

   fclose(fptxt);
   fclose(fpxml);
}


double mango_Search::XcorrScore(const char *szPeptide,
                                int iScanNumber,
                                bool bSilacHeavy)
{

   int iWhichQuery;
   int iLenPeptide = strlen(szPeptide);
   double dXcorr = 0.0;

   // find which
   for (iWhichQuery=0; iWhichQuery<(int)g_pvQuery.size(); iWhichQuery++)
   {
      if (g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iScanNumber == iScanNumber)
         break;
   }

   if (iWhichQuery < (int)g_pvQuery.size() && g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iScanNumber == iScanNumber)
   {
      int bin, x, y;
      int iMax = g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iArraySize/SPARSE_MATRIX_SIZE + 1;

      double dBion = g_staticParams.precalcMasses.dNtermProton;
      double dYion = g_staticParams.precalcMasses.dCtermOH2Proton;

      bool bBionLysine = false; // set to true after first b-ion lysine is modified
      bool bYionLysine = false; // set to true after first y-ion lysine is modified

      for (int i=0; i<iLenPeptide-1; i++) // will ignore multiple fragment ion charge states for now
      {
         dBion += g_staticParams.massUtility.pdAAMassFragment[(int)szPeptide[i]];
         if (szPeptide[i] == 'K' && !bBionLysine)
         {
            dBion += LYSINE_MOD;
            bBionLysine = true;
         }
         if (bSilacHeavy)
         {
            if (szPeptide[i] == 'K')
               dBion += 8.014199;
            if (szPeptide[i] == 'R')
               dBion += 6.020129;
         }

         bin = BIN(dBion);
         x =  bin / SPARSE_MATRIX_SIZE;
         if (!(g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x]==NULL || x>iMax)) // x should never be > iMax so this is just a safety check
         {
            y = bin - (x*SPARSE_MATRIX_SIZE);
            dXcorr += g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x][y];
         }

         dYion += g_staticParams.massUtility.pdAAMassFragment[(int)szPeptide[iLenPeptide -1 - i]];
         if (szPeptide[iLenPeptide -1 - i] == 'K' && !bYionLysine && i>0)
         {
            dYion += LYSINE_MOD;
            bYionLysine = true;
         }
         if (bSilacHeavy)
         {
            if (szPeptide[iLenPeptide -1 - i] == 'K')
               dYion += 8.014199;
            if (szPeptide[iLenPeptide -1 - i] == 'R')
               dYion += 6.020129;
         }

         bin = BIN(dYion);
         x =  bin / SPARSE_MATRIX_SIZE;
         if (!(g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x]==NULL || x>iMax)) // x should never be > iMax so this is just a safety check
         {
            y = bin - (x*SPARSE_MATRIX_SIZE);
            dXcorr += g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x][y];
         }
      }

      if (!bBionLysine && !bYionLysine) // sanity check
      {
         cout << " Error, no internal lysine: " << szPeptide << endl;
         exit(1);
      }

      dXcorr *= 0.005;
   }

   if (dXcorr < 0.0)
      dXcorr = 0.0;

   return dXcorr;
}


void mango_Search::WritePepXMLHeader(FILE *fpxml,
                                     char *szBaseName,
                                     const char *szFastaFile)
{
   char szCWD[1024];
   if (getcwd(szCWD, sizeof(szCWD)) == NULL)
   {
      perror("getcwd() error");
      strcpy(szCWD, "/pwd");
   }

   fprintf(fpxml, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
   fprintf(fpxml, "<msms_pipeline_analysis date= \"2016-11-16T15:59:37\" summary_xml=\"%s/%s.pep.xml\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://regis-web.systemsbiology.net/pepXML http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v120.xsd\">\n", szCWD, szBaseName);
   fprintf(fpxml, " <msms_run_summary base_name=\"%s/%s\" raw_data_type=\"raw\" raw_data=\".mzXML\">\n", szCWD, szBaseName);
   fprintf(fpxml, "  <sample_enzyme name=\"trypsin\">\n");
   fprintf(fpxml, "   <specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>\n");
   fprintf(fpxml, "  </sample_enzyme>\n");
   fprintf(fpxml, "  <search_summary base_name=\"%s/%s\" search_engine=\"Kojak\" search_engine_version=\"1.0.0\" precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" search_id=\"1\">\n", szCWD, szBaseName);
   fprintf(fpxml, "   <search_database local_path=\"%s\" type=\"AA\"/>\n", szFastaFile);
   fprintf(fpxml, "   <enzymatic_search_contstraint enzyme=\"Trypsin\" max_num_internal_cleavages=\"1\" min_number_termini=\"2\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_A\" value=\"0\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_B\" value=\"1\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_C\" value=\"0\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_X\" value=\"0\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_Y\" value=\"1\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_Z\" value=\"0\"/>\n");
   fprintf(fpxml, "  </search_summary>\n");
}


void mango_Search::WriteSpectrumQuery(FILE *fpxml,
                                       char *szBaseName,
                                       double dExpMass1,
                                       double dExpMass2,
                                       double dXcorr1,
                                       double dXcorr2,
                                       double dExpect1,
                                       double dExpect2,
                                       double dCalcMass1,
                                       double dCalcMass2,
                                       double dXcorrCombined,
                                       double dExpectCombined,
                                       char *szPep1,
                                       char *szPep2,
                                       char *szProt1,
                                       char *szProt2,
                                       int iCharge,
                                       int iIndex,
                                       int iScan) 
{                         
   int i;

   // add terminal masses OH + H
   dCalcMass1 += LYSINE_MOD + g_staticParams.massUtility.pdAAMassFragment['o'] + 2*g_staticParams.massUtility.pdAAMassFragment['h'];
   dCalcMass2 += LYSINE_MOD + g_staticParams.massUtility.pdAAMassFragment['o'] + 2*g_staticParams.massUtility.pdAAMassFragment['h'];

   fprintf(fpxml, "  <spectrum_query spectrum=\"%s.%d.%d.%d\" start_scan=\"%d\" end_scan=\"%d\" precursor_neutral_mass=\"%0.6f\" assumed_charge=\"%d\" index=\"%d\">\n",
         szBaseName, iScan, iScan, iCharge, iScan, iScan, dExpMass1+dExpMass2+g_staticParams.options.dReporterMass, iCharge, ++iIndex);
   fprintf(fpxml, "   <search_result>\n");
   fprintf(fpxml, "    <search_hit hit_rank=\"1\" peptide=\"-\" peptide_prev_aa=\"-\" peptide_next_aa=\"-\" protein=\"-\" num_tot_proteins=\"1\" calc_neutral_pep_mass=\"%0.6f\" massdiff=\"%0.6f\" xlink_type=\"xl\">\n",
         dCalcMass1+dCalcMass2+g_staticParams.options.dReporterMass, (dCalcMass1+dCalcMass2)-(dExpMass1+dExpMass2));
   fprintf(fpxml, "     <xlink identifier=\"BDP-NHP\" mass=\"200.00\">\n");
   fprintf(fpxml, "      <linked_peptide peptide=\"%s\" peptide_prev_aa=\"-\" peptide_next_aa=\"-\" protein=\"%s\" num_tot_proteins=\"1\" calc_neutral_pep_mass=\"%0.6f\" complement_mass=\"%0.6f\" designation=\"alpha\">\n",
         szPep1, szProt1, dCalcMass1, dCalcMass2+g_staticParams.options.dReporterMass);
   fprintf(fpxml, "       <modification_info>\n");
   for (i=0; i<(int)strlen(szPep1); i++)
      if (szPep1[i]=='K')
         break;
   fprintf(fpxml, "        <mod_aminoacid_mass position=\"%d\" mass=\"325.127385\"/>\n", i+1);
   for (i=0; i<(int)strlen(szPep1); i++)
      if (szPep1[i]=='C')
         fprintf(fpxml, "        <mod_aminoacid_mass position=\"%d\" mass=\"160.03064805\"/>\n", i+1);
   fprintf(fpxml, "       </modification_info>\n");
   fprintf(fpxml, "       <xlink_score name=\"score\" value=\"%0.3E\"/>\n", dExpect1);
// fprintf(fpxml, "       <xlink_score name=\"xcorr\" value=\"%0.3f\"/>\n", dXcorr1);
   fprintf(fpxml, "      </linked_peptide>\n");
   fprintf(fpxml, "      <linked_peptide peptide=\"%s\" peptide_prev_aa=\"-\" peptide_next_aa=\"-\" protein=\"%s\" num_tot_proteins=\"1\" calc_neutral_pep_mass=\"%0.6f\" complement_mass=\"%0.6f\" designation=\"beta\">\n",
         szPep2, szProt2, dCalcMass2, dCalcMass1+g_staticParams.options.dReporterMass);
   fprintf(fpxml, "       <modification_info>\n");
   for (i=0; i<(int)strlen(szPep2); i++)
      if (szPep2[i]=='K')
         break;
   fprintf(fpxml, "        <mod_aminoacid_mass position=\"%d\" mass=\"325.127385\"/>\n", i+1);
   for (i=0; i<(int)strlen(szPep2); i++)
      if (szPep2[i]=='C')
         fprintf(fpxml, "        <mod_aminoacid_mass position=\"%d\" mass=\"160.03046805\"/>\n", i+1);
   fprintf(fpxml, "       </modification_info>\n");
   fprintf(fpxml, "       <xlink_score name=\"score\" value=\"%0.3E\"/>\n", dExpect2);
// fprintf(fpxml, "       <xlink_score name=\"delta_score\" value=\"%0.3f\"/>\n", dXcorr2);
   fprintf(fpxml, "      </linked_peptide>\n");
   fprintf(fpxml, "     </xlink>\n");
double dScore =  dExpect1 > dExpect2 ? dExpect1 : dExpect2;
//double dScore = dExpectCombined;
   fprintf(fpxml, "     <search_score name=\"kojak_score\" value=\"%0.3E\"/>\n", dScore);
   fprintf(fpxml, "     <search_score name=\"delta_score\" value=\"%0.3f\"/>\n", dExpectCombined); //dXcorrCombined);
   fprintf(fpxml, "     <search_score name=\"ppm_error\" value=\"0.0\"/>\n");
   fprintf(fpxml, "    </search_hit>\n");
   fprintf(fpxml, "   </search_result>\n");
   fprintf(fpxml, "  </spectrum_query>\n");

}

