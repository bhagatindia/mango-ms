/*
   Copyright 2017 University of Washington                          3-clause BSD license

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _MANGOREADFASTA_H_
#define _MANGOREADFASTA_H_

#include "Common.h"
#include "mango_DataInternal.h"

#include "hash/mango-hash.h"

#define NUMPEPTIDES 10

class mango_Search
{
public:
   mango_Search();
   ~mango_Search();

   static void SearchForPeptides(char *szMZXML,
                                 const char *,
                                 enzyme_cut_params,
                                 const char *);

private:

   static void ScorePeptides(protein_hash_db_t phdp,
                             double pep_mass,
                             char *toppep[NUMPEPTIDES],
                             char *toppro[NUMPEPTIDES],
                             float *xcorrPep,
                             vector<double> &vdXcorr_pep,
                             int *hist_pep,
                             int *num_pep,
                             int iScanNumber);

   static double XcorrScore(const char *szPeptide,
                            int iScanNumber);

   static bool CalculateEValue(int *hist_pep,
                               int iMatchPepCount,
                               double *dSlope,
                               double *dIntercept,
                               double dNeutralPepMass,
                               int iScanNumber);

   static void LinearRegression(int *piHistogram,
                                double *slope,
                                double *intercept,
                                int *iMaxXcorr,
                                int *iStartXcorr,
                                int *iNextXcorr);

   static bool GenerateXcorrDecoys(double dNeutralPepMass,
                                   int iMatchPepCount,
                                   int *hist_pep,
                                   int iScanNumber);

   static void WritePepXMLHeader(FILE *fpxml,
                                 char *szBaseName,
                                 const char *szFastaFile,
                                 bool bMimicComet);

   static void WriteSpectrumQuery(FILE *fpxml,
                                  char *szBaseName,
                                  double dExpMass1,
                                  double dExpMass2,
                                  double dXcorr1,
                                  double dXcorr2,
                                  double dDeltaCn1,
                                  double dDeltaCn2,
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
                                  int iScan);

   static void WriteSplitSpectrumQuery(FILE *fpxml,
                                       char *szBaseName,
                                       double dExpMass1,
                                       double dExpMass2,
                                       double dXcorr1,
                                       double dXcorr2,
                                       double dDeltaCn1,
                                       double dDeltaCn2,
                                       double dExpect1,
                                       double dExpect2,
                                       double dCalcMass1,
                                       double dCalcMass2,
                                       char *szPep1,
                                       char *szPep2,
                                       char *szProt1,
                                       char *szProt2,
                                       int iCharge1,
                                       int iCharge2,
                                       int iIndex,
                                       int iScan,
                                       int iWhichDuplicatePrecursor);

   // Private static methods

   // Private member variables
};

#endif // _MANGOREADFASTA_H_
