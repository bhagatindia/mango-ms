/*
   Copyright 2017 University of Washington                          3-clause BSD license

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _MANGODATAINTERNAL_H_
#define _MANGODATAINTERNAL_H_

#include "mango_Data.h"

class MangoSearchManager;

#define PROTON_MASS                 1.00727646688
#define C13_DIFF                    1.00335483

#define FLOAT_ZERO                  1e-6     // 0.000001
#define MAX_PEPTIDE_LEN             64       // max # of AA for a peptide
#define MAX_PEPTIDE_LEN_P2          66       // max # of AA for a peptide plus 2 for N/C-term
#define SIZE_MASS                   128      // ascii value size
#define NUM_SP_IONS                 200      // num ions for preliminary scoring
#define NUM_ION_SERIES              9

#define WIDTH_REFERENCE             512      // size of the protein accession field to store

#define HISTO_SIZE                  152      // some number greater than 150; chose 152 for byte alignment?

#define VMODS                       9
#define VMOD_1_INDEX                0
#define VMOD_2_INDEX                1
#define VMOD_3_INDEX                2
#define VMOD_4_INDEX                3
#define VMOD_5_INDEX                4
#define VMOD_6_INDEX                5
#define VMOD_7_INDEX                6
#define VMOD_8_INDEX                7
#define VMOD_9_INDEX                8

#define ENZYME_SINGLE_TERMINI       1
#define ENZYME_DOUBLE_TERMINI       2
#define ENZYME_N_TERMINI            8
#define ENZYME_C_TERMINI            9

#define ION_SERIES_A                0
#define ION_SERIES_B                1
#define ION_SERIES_C                2
#define ION_SERIES_X                3
#define ION_SERIES_Y                4
#define ION_SERIES_Z                5

#define XCORR_CUTOFF                1E-8   // some near-zero cutoff

struct msdata                    // used in the preprocessing
{
   double dIon;
   double dIntensity;
};

struct Options             // output parameters
{
   int iNumPeptideOutputLines;
   int iWhichReadingFrame;
   int iEnzymeTermini;
   int iNumStored;               // # of search results to store for xcorr analysis
   int iSpectrumBatchSize;       // # of spectra to search at a time within the scan range
   int iStartCharge;
   int iEndCharge;
   int iMaxFragmentCharge;
   int iMaxPrecursorCharge;
   int iMSLevel;                 // mzXML only
   int iMinPeaks;
   int iRemovePrecursor;         // 0=no, 1=yes, 2=ETD precursors
   int iDecoySearch;             // 0=no, 1=concatenated search, 2=separate decoy search
   int iNumThreads;              // 0=poll CPU else set # threads to spawn
   int bOutputTxtFile;
   int bOutputPepXMLFile;
   int bOutputPercolatorFile;
   int bOutputOutFiles;
   int bClipNtermMet;            // 0=leave sequences alone; 1=also consider w/o N-term methionine
   int bSkipAlreadyDone;         // 0=search everything; 1=don't re-search if .out exists
   int bVerboseOutput;
   int bNoEnzymeSelected;
   int bShowFragmentIons;
   int bPrintExpectScore;
   int bOverrideCharge;
   int iMimicCometPepXML;
   int iReportedScore;
   int iSilacHeavy;
   int iDumpRelationshipData;
   double dMinIntensity;
   double dRemovePrecursorTol;
   double dPeptideMassLow;       // MH+ mass
   double dPeptideMassHigh;      // MH+ mass
   double dReporterMass;
   double dLysineStumpMass;

   IntRange scanRange;
   DoubleRange clearMzRange;
   char szActivationMethod[24];  // mzXML only

   Options& operator=(Options& a)
   {
      iNumPeptideOutputLines = a.iNumPeptideOutputLines;
      iWhichReadingFrame = a.iWhichReadingFrame;
      iEnzymeTermini = a.iEnzymeTermini;
      iNumStored = a.iNumStored;
      scanRange = a.scanRange;
      iSpectrumBatchSize = a.iSpectrumBatchSize;
      iStartCharge = a.iStartCharge;
      iEndCharge = a.iEndCharge;
      iMaxFragmentCharge = a.iMaxFragmentCharge;
      iMaxPrecursorCharge = a.iMaxPrecursorCharge;
      iMSLevel = a.iMSLevel;
      iMinPeaks = a.iMinPeaks;
      dMinIntensity = a.dMinIntensity;
      iRemovePrecursor = a.iRemovePrecursor;
      iDecoySearch = a.iDecoySearch;
      iNumThreads = a.iNumThreads;
      bOutputTxtFile = a.bOutputTxtFile;
      bOutputPepXMLFile = a.bOutputPepXMLFile;
      bOutputPercolatorFile = a.bOutputPercolatorFile;
      bOutputOutFiles = a.bOutputOutFiles;
      bClipNtermMet = a.bClipNtermMet;
      bSkipAlreadyDone = a.bSkipAlreadyDone;
      bVerboseOutput = a.bVerboseOutput;
      bNoEnzymeSelected = a.bNoEnzymeSelected;
      bShowFragmentIons = a.bShowFragmentIons;
      bPrintExpectScore = a.bPrintExpectScore;
      dRemovePrecursorTol = a.dRemovePrecursorTol;
      clearMzRange = a.clearMzRange;
      dPeptideMassLow = a.dPeptideMassLow;
      dPeptideMassHigh = a.dPeptideMassHigh;
      dReporterMass = a.dReporterMass;
      dLysineStumpMass = a.dLysineStumpMass;
      iMimicCometPepXML = a.iMimicCometPepXML;
      iReportedScore = a.iReportedScore;
      iSilacHeavy = a.iSilacHeavy;
      iDumpRelationshipData = a.iDumpRelationshipData;
      strcpy(szActivationMethod, a.szActivationMethod);

      return *this;
   }
};

struct Results
{
   double dPepMass;
   double dExpect;
   float  fScoreSp;
   float  fXcorr;
   int    iDuplicateCount;
   int    iLenPeptide;
   int    iRankSp;
   int    iMatchedIons;
   int    iTotalIons;
   long   lProteinFilePosition;
   int    piVarModSites[MAX_PEPTIDE_LEN_P2];   // store variable mods encoding, +2 to accomodate N/C-term
   double pdVarModSites[MAX_PEPTIDE_LEN_P2];   // store variable mods mass diffs, +2 to accomodate N/C-term
// char szProtein[WIDTH_REFERENCE];
   char   szPeptide[MAX_PEPTIDE_LEN];
   char   szPrevNextAA[2];                    // [0] stores prev AA, [1] stores next AA
};

struct PepMassInfo
{
   double dCalcPepMass;
   double dExpPepMass;
   double dPeptideMassTolerance;
   double dPeptideMassToleranceMinus;
   double dPeptideMassTolerancePlus;
};

struct SpectrumInfoInternal
{
   int iArraySize;         // m/z versus intensity array
   int iHighestIon;
   int iScanNumber;
   int iChargeState;
   int iMaxFragCharge;
   double dTotalIntensity;
   double dRTime;
   char szNativeID[128];   // nativeID string from mzML
};

// The minimum and maximum mass range of all peptides to consider
// i.e. lowestPepMass - tolerance to highestPepMass + tolerance
struct MassRange
{
   double dMinMass;
   double dMaxMass;
   int    iMaxFragmentCharge;  // global maximum fragment charge
};

extern MassRange g_massRange;

// PreprocessStruct stores information used in preprocessing
// each spectrum.  Information not kept around otherwise
struct PreprocessStruct
{
   int    iHighestIon;
   double dHighestIntensity;
   struct msdata pTmpSpData[NUM_SP_IONS];
};


//-->MH
typedef struct sDBEntry
{
   string strName;           // might be able to delete this here
   string strSeq;
   long lProteinFilePosition;
} sDBEntry;

struct DBInfo
{
   char szDatabase[SIZE_FILE];
   char szHash[SIZE_FILE];
   char szFileName[SIZE_FILE];
   int  iTotalNumProteins;
   unsigned long int uliTotAACount;

   DBInfo& operator=(DBInfo& a)
   {
      strcpy(szDatabase, a.szDatabase);
      strcpy(szHash, a.szHash);
      strcpy(szFileName, a.szFileName);
      iTotalNumProteins = a.iTotalNumProteins;
      uliTotAACount = a.uliTotAACount;

      return *this;
   }
};

struct StaticMod
{
   double dAddCterminusPeptide;
   double dAddNterminusPeptide;
   double dAddCterminusProtein;
   double dAddNterminusProtein;
   double pdStaticMods[SIZE_MASS];

   StaticMod& operator=(StaticMod& a)
   {
      dAddCterminusPeptide = a.dAddCterminusPeptide;
      dAddNterminusPeptide = a.dAddNterminusPeptide;
      dAddCterminusProtein = a.dAddCterminusProtein;
      dAddNterminusProtein = a.dAddNterminusProtein;

      for (int i = 0; i < SIZE_MASS; i++)
      {
         pdStaticMods[i] = a.pdStaticMods[i];
      }

      return *this;
   }
};

struct PrecalcMasses
{
   double dNtermProton;          // dAddNterminusPeptide + PROTON_MASS
   double dCtermOH2Proton;       // dAddCterminusPeptide + dOH2fragment + PROTON_MASS
   double dOH2ProtonCtermNterm;  // dOH2parent + PROTON_MASS + dAddCterminusPeptide + dAddNterminusPeptide
   double dOH2;                  // mass of OH2
   int iMinus17;                 // BIN'd value of mass(NH3)
   int iMinus18;                 // BIN'd value of mass(H2O)

   PrecalcMasses& operator=(PrecalcMasses& a)
   {
      dNtermProton = a.dNtermProton;
      dCtermOH2Proton = a.dCtermOH2Proton;
      dOH2ProtonCtermNterm = a.dOH2ProtonCtermNterm;
      dOH2 = a.dOH2 ;
      iMinus17 = a.iMinus17;
      iMinus18 = a.iMinus18;

      return *this;
   }
};

struct VarModParams
{
   bool    bVarModSearch;            // set to true if variable mods are specified
   bool    bBinaryModSearch;         // set to true if any of the variable mods are of binary mod variety
   int     bRequireVarMod;           // also set to true if any individual bRequireThisMod is true
   int     iMaxVarModPerPeptide;
   VarMods varModList[VMODS];

   VarModParams& operator=(VarModParams& a)
   {
      bVarModSearch = a.bVarModSearch;
      iMaxVarModPerPeptide = a.iMaxVarModPerPeptide;

      for (int i = 0; i < VMODS; i++)
      {
         varModList[i] = a.varModList[i];
      }

      return *this;
   }
};

struct MassUtil
{
   int    bMonoMassesParent;
   int    bMonoMassesFragment;
   double dCO;
   double dNH3;
   double dNH2;
   double dH2O;
   double dCOminusH2;
   double dOH2fragment;
   double dOH2parent;
   double pdAAMassParent[SIZE_MASS];
   double pdAAMassFragment[SIZE_MASS];

   MassUtil& operator=(MassUtil& a)
   {
      bMonoMassesParent = a.bMonoMassesParent;
      bMonoMassesFragment = a.bMonoMassesFragment;
      dCO = a.dCO;
      dNH3 = a.dNH3;
      dNH2 = a.dNH2;
      dH2O = a.dH2O;
      dCOminusH2 = a.dCOminusH2;
      dOH2fragment = a.dOH2fragment;
      dOH2parent = a.dOH2parent;

      for (int i = 0; i < SIZE_MASS; i++)
      {
         pdAAMassParent[i] = a.pdAAMassParent[i];
         pdAAMassFragment[i] = a.pdAAMassFragment[i];
      }

      return *this;
   }
};

struct ToleranceParams
{
   int    iMassToleranceUnits;    // 0=ppm, 1=da (default)
   int    iMassToleranceType;     // 0=MH+ (default), 1=precursor m/z; only valid if iMassToleranceUnits > 0
   int    iIsotopeError;
   double dTolerancePeptide;        // tolerance from param file
   double dToleranceRelationship;   // tolerance from param file
   double dFragmentBinSize;
   double dFragmentBinStartOffset;

   ToleranceParams& operator=(ToleranceParams& a)
   {
      iMassToleranceUnits = a.iMassToleranceUnits;
      iMassToleranceType = a.iMassToleranceType;
      iIsotopeError = a.iIsotopeError;
      dTolerancePeptide = a.dTolerancePeptide;
      dToleranceRelationship = a.dToleranceRelationship;
      dFragmentBinSize = a.dFragmentBinSize;
      dFragmentBinStartOffset = a.dFragmentBinStartOffset;

      return *this;
   }
};


struct IonInfo
{
   int iNumIonSeriesUsed;
   int piSelectedIonSeries[NUM_ION_SERIES];
   int bUseNeutralLoss;
   int iTheoreticalFragmentIons;
   int iIonVal[NUM_ION_SERIES];

   IonInfo& operator=(IonInfo& a)
   {
      iNumIonSeriesUsed = a.iNumIonSeriesUsed;
      bUseNeutralLoss = a.bUseNeutralLoss;
      iTheoreticalFragmentIons = a.iTheoreticalFragmentIons;

      for (int i = 0; i < NUM_ION_SERIES; i++)
      {
         piSelectedIonSeries[i] = a.piSelectedIonSeries[i];
         iIonVal[i] = a.iIonVal[i];
      }

      return *this;
   }
};

// static user params, won't change per thread - can make global!
struct StaticParams
{
   char            szHostName[SIZE_FILE];
   char            szOutFileTimeString[256];
   char            szIonSeries[256];   // used for .out files
   char            szDisplayLine[256]; // used for .out files
   char            szMod[512];         // used for .out files
   char            szDecoyPrefix[256]; // used for prefix to indicate decoys
   char            szOutputSuffix[256]; // used for suffix to append to output file base names
   int             iElapseTime;
   char            szDate[32];
   Options         options;
   DBInfo          databaseInfo;
   InputFileInfo   inputFile;
   int             bPrintDuplReferences;
   VarModParams    variableModParameters;
   ToleranceParams tolerances;
   StaticMod       staticModifications;
   PrecalcMasses   precalcMasses;
   EnzymeInfo      enzymeInformation;
   MassUtil        massUtility;
   double          dInverseBinWidth;    // this is used in BIN() many times so use inverse binWidth to do multiply vs. divide
   double          dOneMinusBinOffset;  // this is used in BIN() many times so calculate once
   IonInfo         ionInformation;
   int             iXcorrProcessingOffset;
   vector<double>  vectorMassOffsets;

   StaticParams()
   {
      int i;

      inputFile.iInputType = InputType_MS2;

      szMod[0] = '\0';

      iXcorrProcessingOffset = 75;

      databaseInfo.szDatabase[0] = '\0';
      databaseInfo.szHash[0] = '\0';

      strcpy(szDecoyPrefix, "DECOY_");
      szOutputSuffix[0] = '\0';

      for (i=0; i<SIZE_MASS; i++)
      {
         massUtility.pdAAMassParent[i] = 999999.;
         massUtility.pdAAMassFragment[i] = 999999.;
      }

      massUtility.bMonoMassesFragment = 1;
      massUtility.bMonoMassesParent = 1;

      for (int i=0; i<SIZE_MASS; i++)
      {
         staticModifications.pdStaticMods[i] = 0.0;
      }

//    staticModifications.pdStaticMods[(int)'C'] = 57.021464;

      enzymeInformation.iAllowedMissedCleavage = 2;

      for (i=0; i<VMODS; i++)
      {
         variableModParameters.varModList[i].iMaxNumVarModAAPerMod = 3;
         variableModParameters.varModList[i].iBinaryMod = 0;
         variableModParameters.varModList[i].bRequireThisMod = 0;
         if (i==0)
         {
            variableModParameters.varModList[i].dVarModMass = 15.9949;
            strcpy(variableModParameters.varModList[i].szVarModChar, "M");
         }
         else
         {
            variableModParameters.varModList[i].dVarModMass = 0.0;
            strcpy(variableModParameters.varModList[i].szVarModChar, "X");
         }
         variableModParameters.varModList[i].iVarModTermDistance = -1;   // distance from N or C-term distance
         variableModParameters.varModList[i].iWhichTerm = 0;             // specify N (0) or C-term (1)

      }

      variableModParameters.iMaxVarModPerPeptide = 5;

      ionInformation.bUseNeutralLoss = 0;
      ionInformation.iTheoreticalFragmentIons = 0;      // 0 = flanking peaks; 1 = no flanking peaks
      ionInformation.iIonVal[ION_SERIES_A] = 0;
      ionInformation.iIonVal[ION_SERIES_B] = 1;
      ionInformation.iIonVal[ION_SERIES_C] = 0;
      ionInformation.iIonVal[ION_SERIES_X] = 0;
      ionInformation.iIonVal[ION_SERIES_Y] = 1;
      ionInformation.iIonVal[ION_SERIES_Z] = 0;

      options.iNumPeptideOutputLines = 5;
      options.iWhichReadingFrame = 0;
      options.iEnzymeTermini = 2;
      options.iNumStored = 100;                         // default # of search results to store for xcorr analysis.

      options.bNoEnzymeSelected = 1;
      options.bShowFragmentIons = 0;
      options.bPrintExpectScore = 1;
      options.bOverrideCharge = 0;
      options.iRemovePrecursor = 0;
      options.dRemovePrecursorTol = 1.0;

      options.bOutputTxtFile = 0;
      options.bOutputPepXMLFile = 1;
      options.bOutputPercolatorFile = 0;
      options.bOutputOutFiles = 0;

      options.bSkipAlreadyDone = 1;
      options.bVerboseOutput = 0;
      options.iDecoySearch = 0;
      options.iNumThreads = 0;
      options.bClipNtermMet = 0;

      // These parameters affect mzXML/RAMP spectra only.
      options.scanRange.iStart = 0;
      options.scanRange.iEnd = 0;
      options.iSpectrumBatchSize = 0;
      options.iMinPeaks = 10;
      options.iStartCharge = 0;
      options.iEndCharge = 0;
      options.iMaxFragmentCharge = 3;
      options.iMaxPrecursorCharge = 6;
      options.iMSLevel = 2;
      options.dMinIntensity = 0.0;
      options.dPeptideMassLow = 2000.0;
      options.dPeptideMassHigh = 8000.0;
      strcpy(options.szActivationMethod, "ALL");
      // End of mzXML specific parameters.

      options.dReporterMass = 751.40508;
      options.dLysineStumpMass = 197.032422;
      options.iMimicCometPepXML = 0;
      options.iReportedScore = 0;
      options.iSilacHeavy = 0;
      options.iDumpRelationshipData= 0;

      options.clearMzRange.dStart = 0.0;
      options.clearMzRange.dEnd = 0.0;

      staticModifications.dAddCterminusPeptide = 0.0;
      staticModifications.dAddNterminusPeptide = 0.0;
      staticModifications.dAddCterminusProtein = 0.0;
      staticModifications.dAddNterminusProtein = 0.0;

      tolerances.iMassToleranceUnits = 2;
      tolerances.iMassToleranceType = 1;
      tolerances.iIsotopeError = 0;
      tolerances.dTolerancePeptide = 20.0;                    // peptide_mass_tolerance
      tolerances.dToleranceRelationship = 50.0;                    // peptide_mass_tolerance
      tolerances.dFragmentBinSize = 0.02;
      tolerances.dFragmentBinStartOffset = 0.0;
   }

   StaticParams& operator=(StaticParams& a)
   {
      strcpy(szHostName, a.szHostName);
      strcpy(szOutFileTimeString, a.szOutFileTimeString);
      strcpy(szIonSeries, a.szIonSeries);
      strcpy(szDisplayLine, a.szDisplayLine);
      strcpy(szMod, a.szMod);
      strcpy(szDecoyPrefix, a.szDecoyPrefix);
      strcpy(szOutputSuffix, a.szOutputSuffix);
      vectorMassOffsets = a.vectorMassOffsets;
      iElapseTime = a.iElapseTime;
      strcpy(szDate, a.szDate);
      options = a.options;
      databaseInfo = a.databaseInfo;
      inputFile = a.inputFile;
      bPrintDuplReferences = a.bPrintDuplReferences;
      variableModParameters = a.variableModParameters;
      tolerances = a.tolerances;
      staticModifications = a.staticModifications;
      precalcMasses = a.precalcMasses;
      enzymeInformation = a.enzymeInformation;
      massUtility = a.massUtility;
      dInverseBinWidth = a.dInverseBinWidth;
      dOneMinusBinOffset = a.dOneMinusBinOffset;
      iXcorrProcessingOffset = a.iXcorrProcessingOffset;
      ionInformation = a.ionInformation;
      return *this;
   }

};

extern StaticParams g_staticParams;

struct SparseMatrix
{
   int bin;
   float fIntensity;
};

// Query stores information for peptide scoring and results
// This struct is allocated for each spectrum/charge combination
struct Query
{
   int   iXcorrHistogram[HISTO_SIZE];
   int   iHistogramCount;   // # of entries in histogram
   float fPar[4];           // parameters of LMA regression

   int iMatchPeptideCount;        // # of peptide matches
   int iDecoyMatchPeptideCount;   // # of decoy peptide matches

   short siMaxXcorr;        // index of maximum correlation score in iXcorrHistogram

   short siLowestSpScoreIndex;
   short siLowestDecoySpScoreIndex;

   float fLowestSpScore;
   float fLowestDecoySpScore;

   float fLowestCorrScore;
   float fLowestDecoyCorrScore;

   unsigned long int  _uliNumMatchedPeptides;
   unsigned long int  _uliNumMatchedDecoyPeptides;

   // Sparse matrix representation of data
   int iFastXcorrData;  //MH: I believe these are all the same size now.
   float **ppfSparseFastXcorrData;

   PepMassInfo          _pepMassInfo;
   SpectrumInfoInternal _spectrumInfoInternal;
   Results              *_pResults;
   Results              *_pDecoys;

   Query()
   {
      for (int i=0; i < HISTO_SIZE; i++)
         iXcorrHistogram[i] = 0;

      iMatchPeptideCount= 0;
      iDecoyMatchPeptideCount= 0;
      iHistogramCount = 0;

      fPar[0]=0.0;
      fPar[1]=0.0;
      fPar[2]=0.0;
      fPar[3]=0.0;

      siMaxXcorr = 0;                        // index of maximum correlation score in iXcorrHistogram
      siLowestSpScoreIndex = 0;
      siLowestDecoySpScoreIndex = 0;

      fLowestSpScore = 0.0;
      fLowestDecoySpScore = 0.0;

      fLowestCorrScore = XCORR_CUTOFF;
      fLowestDecoyCorrScore = XCORR_CUTOFF;

      _uliNumMatchedPeptides = 0;
      _uliNumMatchedDecoyPeptides = 0;

      ppfSparseFastXcorrData = NULL;

      _pepMassInfo.dCalcPepMass = 0.0;
      _pepMassInfo.dExpPepMass = 0.0;
      _pepMassInfo.dPeptideMassTolerance = 0.0;
      _pepMassInfo.dPeptideMassToleranceMinus = 0.0;
      _pepMassInfo.dPeptideMassTolerancePlus = 0.0;

      _spectrumInfoInternal.dTotalIntensity = 0.0;
      _spectrumInfoInternal.iArraySize = 0;
      _spectrumInfoInternal.iHighestIon = 0;
      _spectrumInfoInternal.iScanNumber = 0;
      _spectrumInfoInternal.dTotalIntensity = 0.0;

      _pResults = NULL;
      _pDecoys = NULL;
   }

   ~Query()
   {
      int i;

      for (i=0;i<iFastXcorrData;i++)
      {
         if (ppfSparseFastXcorrData[i] != NULL)
            delete[] ppfSparseFastXcorrData[i];
      }
      delete[] ppfSparseFastXcorrData;
      ppfSparseFastXcorrData = NULL;

      delete[] _pResults;
      _pResults = NULL;

      if (g_staticParams.options.iDecoySearch==2)
      {
         delete[] _pDecoys;
         _pDecoys = NULL;
      }
   }
};

extern vector<Query*>          g_pvQuery;
extern vector<InputFileInfo*>  g_pvInputFiles;

struct IonSeriesStruct         // defines which fragment ion series are considered
{
   int bPreviousMatch[8];
};

#endif // _XLINKSDATAINTERNAL_H_
