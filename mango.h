#ifndef _XLINKX_H_
#define _XLINKX_H_

#define SIZE_BUF    8192
#define SIZE_FILE   512
#define MAX_PEAKS   1000

using namespace std;

#include "MSReader.h"
#include "Spectrum.h"
#include "MSObject.h"
#include "math.h"
#include <vector>
#include <cfloat>
#include "mango_Interfaces.h"
#include "mango_SearchManager.h"

using namespace MSToolkit;

struct PrecursorsStruct
{
   int    iCharge1;
   int    iCharge2;
   int    iIntensity1;
   int    iIntensity2;
   double dNeutralMass1;
   double dNeutralMass2;
};

// MS/MS scan information
struct ScanDataStruct
{
   int    iScanNumber;
   int    iPrecursorScanNumber;
   int    iPrecursorCharge;   // charge state of intact precursor
   double dPrecursorMZ;       // intact precursor m/z
   double dPrecursorNeutralMass;       // intact precursor Neutralmass
   double dHardklorPrecursorNeutralMass;  // precursor m/z found from Hardklor

   vector<PrecursorsStruct> pvdPrecursors;
};

struct ParamsStruct
{
   int iResolution1;
   int iResolution2;
   int iChargeMin1;
   int iChargeMin2;
   int iChargeMax1;
   int iChargeMax2;

   double dToleranceMS1;
   double dToleranceMS2;
   double dMassReporter;
   double dMassMod;
};

extern vector<ScanDataStruct> pvSpectrumList;

void Usage(char *pszCmd);
void ProcessCmdLine(int argc,
                    char *argv[],
                    char *szParamsFile,
                    vector<InputFileInfo*> &pvInputFiles,
                    IMangoSearchManager *pSearchMgr);
void SetOptions(char *arg,
                char *szParamsFile,
                bool *bPrintParams,
                IMangoSearchManager *pSearchMgr);
void LoadParameters(char *pszParamsFile,
                    IMangoSearchManager *pSearchMgr);
void PrintParams(void);




#endif // _XLINKX_H_
