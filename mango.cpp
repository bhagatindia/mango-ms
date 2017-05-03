#include <stdio.h>
#include <stdlib.h>

#include "Common.h"
#include "mango.h"
#include "mango_DataInternal.h"
#include "mango_Preprocess.h"
#include "mango_Search.h"
#include "mango_MassSpecUtils.h"
#include "mango_SearchManager.h"


// This program takes in an mzXML as input.
// If corresponding .hk file does not exist, it runs Hardklor to create it.
// Then reading each MS/MS scan in mzXML, find 2 peaks in .hk file
// that add up to ReACT relationship; need mzXML for precursor mass.

vector<ScanDataStruct> pvSpectrumList;
struct ParamsStruct pParams;

int main(int argc, char **argv)
{
   printf("\n Mango version \"%s\"\n\n", mango_version);

   if (argc < 2)
      Usage(argv[0]);

   vector<InputFileInfo*> pvInputFiles;
   IMangoSearchManager* pMangoSearchMgr = GetMangoSearchManager();
   char szParamsFile[SIZE_FILE];

   ProcessCmdLine(argc, argv, szParamsFile, pvInputFiles, pMangoSearchMgr);
   pMangoSearchMgr->AddInputFiles(pvInputFiles);

   bool bSearchSucceeded = pMangoSearchMgr->DoSearch();

   ReleaseMangoSearchManager();

   if (!bSearchSucceeded)
      exit(1);

   return 0;

} //main


void Usage(char *pszCmd)
{
   printf(" Mango usage:  %s [options] <input_files>\n", pszCmd);
   printf("\n");
   printf("       options:  -p         to print out a mango.params file (named mango.params.new)\n");
   printf("\n");
   printf(" Supported input formats include mzXML, mzML\n");
   printf("\n");
   printf(" Example:  %s file1.mzXML file2.mzXML\n", pszCmd);
   printf("\n");

   exit(1);
}


void ProcessCmdLine(int argc,
                    char *argv[],
                    char *szParamsFile,
                    vector<InputFileInfo*> &pvInputFiles,
                    IMangoSearchManager *pSearchMgr)
{
   bool bPrintParams = false;
   int iStartInputFile = 1;
   char *arg;

   if (iStartInputFile == argc)
   {  
      char szErrorMsg[256];
      sprintf(szErrorMsg+strlen(szErrorMsg), " Error - no input files specified so nothing to do.\n");
      logerr(szErrorMsg);
      exit(1);
   }

   strcpy(szParamsFile, "mango.params");

   arg = argv[iStartInputFile];

   // First process the command line options; do this only to see if an alternate
   // params file is specified before loading params file first.
   while ((iStartInputFile < argc) && (NULL != arg))
   {
      if (arg[0] == '-')
         SetOptions(arg, &bPrintParams);

      arg = argv[++iStartInputFile];
   }

   if (bPrintParams)
   {
      PrintParams();
      exit(0);
   }

   // Loads search parameters from mango.params file. This has to happen
   // after parsing command line arguments in case alt file is specified.
   LoadParameters(szParamsFile, pSearchMgr);

   // Now go through input arguments again.  Command line options will
   // override options specified in params file.
   iStartInputFile = 1;
   arg = argv[iStartInputFile];
   while ((iStartInputFile < argc) && (NULL != arg))
   {
      if (arg[0] == '-')
      {
         SetOptions(arg, &bPrintParams);
      }
      else if (arg != NULL)
      {
         // Get the file name. Because Windows can have ":" in the file path,
         // we can't just use "strtok" to grab the filename.
         int i;
         int iArgLen = strlen(arg);
         for (i=0; i < iArgLen; i++)
         {
            if (arg[i] == ':')
            {
               if ((i + 1) < iArgLen)
               {
                  if (arg[i+1] != '\\' && arg[i+1] != '/')
                  {
                     break;
                  }
               }
            }
         }

         InputFileInfo *pInputFileInfo = new InputFileInfo();
         strncpy(pInputFileInfo->szFileName, arg, i);
         pInputFileInfo->szFileName[i] = '\0';
         pvInputFiles.push_back(pInputFileInfo);
      }
      else
      {
         break;
      }

      arg = argv[++iStartInputFile];
   }
} // ProcessCmdLine


void SetOptions(char *arg,
                bool *bPrintParams)
{
   switch (arg[1])
   {
      case 'p':
         *bPrintParams = true;
         break;
      default:
         break;
   }
}


void PrintParams(void)
{
   FILE *fp;

   if ( (fp=fopen("mango.params.new", "w"))==NULL)
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg+strlen(szErrorMsg), " Error - cannot write file mango.params.new\n");
      logerr(szErrorMsg);
      exit(1);
   }

   fprintf(fp, "# mango_version %s\n", mango_version);
   fprintf(fp, "fasta_file = /some/path/db.fasta\n");
   fprintf(fp, "fasta_hash = /some/path/db.fasta.hash            # if exits, use hash; if not create hash from fasta\n");
   fprintf(fp, "mass_tolerance_relationship = %0.2f              # PPM units; intact crosslink vs. 751 + mass1 + mass2\n", g_staticParams.tolerances.dToleranceRelationship);
   fprintf(fp, "mass_tolerance_peptide = %0.2f                   # PPM units; mass1 vs. retrieved peptides from database\n", g_staticParams.tolerances.dTolerancePeptide);
   fprintf(fp, "mass_tolerance_fragment = %0.2f                   # Da bin size\n", g_staticParams.tolerances.dFragmentBinSize);
   fprintf(fp, "reporter_neutral_mass = %f\n", g_staticParams.options.dReporterMass);
   fprintf(fp, "lysine_stump_mass = %f\n", g_staticParams.options.dLysineStumpMass);
   fprintf(fp, "mimic_comet_pepxml = %d                          # if 1, will write out IDs as separate spectrum_query entries\n", g_staticParams.options.iMimicCometPepXML);
   fprintf(fp, "reported_score = %d                              # # 0=worst E-value; 1=combined E-value\n", g_staticParams.options.iReportedScore);
   fprintf(fp, "silac_heavy = %d                                 # 0=normal/light search; 1=SILAC heavy search\n", g_staticParams.options.iSilacHeavy);
   fprintf(fp, "dump_relationship_data = %d                      # 0=no, 1=yes, 2=yes but do not do search\n", g_staticParams.options.iDumpRelationshipData);
      
   logout("\n Created:  mango.params.new\n\n");
   fclose(fp);
}


// Reads comet.params parameter file.
void LoadParameters(char *pszParamsFile,
                    IMangoSearchManager *pSearchMgr)
{
   double dDoubleParam;
   int   iIntParam;
   char  szParamBuf[SIZE_BUF],
         szParamName[128],
         szParamVal[512],
         szParamStringVal[512],
         szVersion[128],
         szErrorMsg[512];
   FILE  *fp;
   bool  bValidParamsFile;
   char *pStr;
   VarMods varModsParam;
   IntRange intRangeParam;
   DoubleRange doubleRangeParam;

   if ((fp=fopen(pszParamsFile, "r")) == NULL)
   {
      sprintf(szErrorMsg+strlen(szErrorMsg), " Error - cannot open parameter file \"%s\".\n", pszParamsFile);
      logerr(szErrorMsg);
      exit(1);
   }

   // validate not using incompatible params file by checking "# mango_version" in first line of file
   strcpy(szVersion, "unknown");
   bValidParamsFile = false;
   while (!feof(fp))
   {
      if (fgets(szParamBuf, SIZE_BUF, fp) != NULL)
      {
         if (!strncmp(szParamBuf, "# mango_version ", 16))
         {
            char szRev1[12],
                 szRev2[12];

            sscanf(szParamBuf, "%*s %*s %128s %12s %12s", szVersion, szRev1, szRev2);

            if (pSearchMgr->IsValidMangoVersion(string(szVersion)))
            {
               bValidParamsFile = true;
               sprintf(szVersion, "%s %s %s", szVersion, szRev1, szRev2);
               pSearchMgr->SetParam("# mango_version ", szVersion, szVersion);
               break;
            }
         }
      }
   }

   if (!bValidParamsFile)
   {
      sprintf(szErrorMsg+strlen(szErrorMsg), " The mango.params file is from version %s\n", szVersion);
      sprintf(szErrorMsg+strlen(szErrorMsg), " Please update your mango.params file.  You can generate\n");
      sprintf(szErrorMsg+strlen(szErrorMsg), " a new parameters file using \"mango -p\"\n\n");
      logerr(szErrorMsg);
      exit(1);
   }

   rewind(fp);

   // now parse the parameter entries
   while (!feof(fp))
   {
      if (fgets(szParamBuf, SIZE_BUF, fp) != NULL)
      {
         if ( (pStr = strchr(szParamBuf, '#')) != NULL)  // take care of comments
            *pStr = 0;

         if ( (pStr = strchr(szParamBuf, '=')) != NULL)
         {
            strcpy(szParamVal, pStr + 1);  // Copy over value.
            *pStr = 0;                     // Null rest of szParamName at equal char.

            sscanf(szParamBuf, "%128s", szParamName);

            if (!strcmp(szParamName, "fasta_file"))
            {
               char szFile[512];

               // Support parsing a database string from params file that
               // includes spaces in the path.

               // Remove white spaces at beginning/end of szParamVal
               int iLen = strlen(szParamVal);
               char *szTrimmed = szParamVal;

               while (isspace(szTrimmed[iLen -1]))  // trim end
                  szTrimmed[--iLen] = 0;
               while (*szTrimmed && isspace(*szTrimmed))  // trim beginning
               {
                  ++szTrimmed;
                  --iLen;
               }

               memmove(szParamVal, szTrimmed, iLen+1);

               strcpy(szFile, szParamVal);
               pSearchMgr->SetParam("fasta_file", szFile, szFile);
            }
            else if (!strcmp(szParamName, "fasta_hash"))
            {
               char szFile[512];
               // Remove white spaces at beginning/end of szParamVal
               int iLen = strlen(szParamVal);
               char *szTrimmed = szParamVal;

               while (isspace(szTrimmed[iLen -1]))  // trim end
                  szTrimmed[--iLen] = 0;
               while (*szTrimmed && isspace(*szTrimmed))  // trim beginning
               {
                  ++szTrimmed;
                  --iLen;
               }

               memmove(szParamVal, szTrimmed, iLen+1);

               strcpy(szFile, szParamVal);
               pSearchMgr->SetParam("fasta_hash", szFile, szFile);
            }
            else if (!strcmp(szParamName, "mass_tolerance_relationship"))
            {  
               sscanf(szParamVal, "%lf", &dDoubleParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%lf", dDoubleParam);
               pSearchMgr->SetParam("mass_tolerance_relationship", szParamStringVal, dDoubleParam);
            }
            else if (!strcmp(szParamName, "mass_tolerance_peptide"))
            {  
               sscanf(szParamVal, "%lf", &dDoubleParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%lf", dDoubleParam);
               pSearchMgr->SetParam("mass_tolerance_peptide", szParamStringVal, dDoubleParam);
            }
            else if (!strcmp(szParamName, "mass_tolerance_fragment"))
            {  
               sscanf(szParamVal, "%lf", &dDoubleParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%lf", dDoubleParam);
               pSearchMgr->SetParam("mass_tolerance_fragment", szParamStringVal, dDoubleParam);
            }
            else if (!strcmp(szParamName, "reporter_neutral_mass"))
            {  
               sscanf(szParamVal, "%lf", &dDoubleParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%lf", dDoubleParam);
               pSearchMgr->SetParam("reporter_neutral_mass", szParamStringVal, dDoubleParam);
            }
            else if (!strcmp(szParamName, "lysine_stump_mass"))
            {  
               sscanf(szParamVal, "%lf", &dDoubleParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%lf", dDoubleParam);
               pSearchMgr->SetParam("lysine_stump_mass", szParamStringVal, dDoubleParam);
            }
            else if (!strcmp(szParamName, "mimic_comet_pepxml"))
            {  
               sscanf(szParamVal, "%d", &iIntParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%d", iIntParam); 
               pSearchMgr->SetParam("mimic_comet_pepxml", szParamStringVal, iIntParam);
            }
            else if (!strcmp(szParamName, "reported_score"))
            {  
               sscanf(szParamVal, "%d", &iIntParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%d", iIntParam); 
               pSearchMgr->SetParam("reported_score", szParamStringVal, iIntParam);
            }
            else if (!strcmp(szParamName, "silac_heavy"))
            {  
               sscanf(szParamVal, "%d", &iIntParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%d", iIntParam); 
               pSearchMgr->SetParam("silac_heavy", szParamStringVal, iIntParam);
            }
            else if (!strcmp(szParamName, "dump_relationship_data"))
            {  
               sscanf(szParamVal, "%d", &iIntParam);
               szParamStringVal[0] = '\0';
               sprintf(szParamStringVal, "%d", iIntParam); 
               pSearchMgr->SetParam("dump_relationship_data", szParamStringVal, iIntParam);
            }
            else
            {
               sprintf(szErrorMsg, " Warning - invalid parameter found: %s.  Parameter will be ignored.\n", szParamName);
               logout(szErrorMsg);
            }
         }
      }
   }

   fclose(fp);
}
