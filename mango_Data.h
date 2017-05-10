/*
   Copyright 2017 University of Washington                          3-clause BSD license

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <string>

#ifndef _MANGODATA_H_
#define _MANGODATA_H_

#define SIZE_BUF                    8192
#define SIZE_FILE                   512

#define MAX_THREADS                 64

#define MAX_ENZYME_AA               20       // max # of AA for enzyme break point
#define MAX_VARMOD_AA               20       // max # of modified AAs in a peptide per variable modification

#define ENZYME_NAME_LEN             48

#define MAX_FRAGMENT_CHARGE         5
#define MAX_PRECURSOR_CHARGE        9

#define MAX_PERMUTATIONS            100000

#define SPARSE_MATRIX_SIZE          100

struct DoubleRange
{
   double dStart;
   double dEnd;

   DoubleRange()
   {
      dStart = 0.0;
      dEnd = 0.0;
   }

   DoubleRange(const DoubleRange& a)
   {
      dStart = a.dStart;
      dEnd = a.dEnd;
   }

   DoubleRange(double dStart_in, double dEnd_in)
   {
      dStart = dStart_in;
      dEnd = dEnd_in;
   }

   DoubleRange& operator=(DoubleRange& a)
   {
      dStart = a.dStart;
      dEnd = a.dEnd;
      return *this;
   }
};

struct IntRange
{
   int iStart;
   int iEnd;

   IntRange()
   {
      iStart = 0;
      iEnd = 0;
   }

   IntRange(const IntRange& a)
   {
      iStart = a.iStart;
      iEnd = a.iEnd;
   }

   IntRange(int iStart_in, int iEnd_in)
   {
      iStart = iStart_in;
      iEnd = iEnd_in;
   }

   IntRange& operator=(IntRange& a)
   {
      iStart = a.iStart;
      iEnd = a.iEnd;
      return *this;
   }
};

struct VarMods
{
   double dVarModMass;
   int    iBinaryMod;
   int    bRequireThisMod;
   int    iMaxNumVarModAAPerMod;
   int    iVarModTermDistance;
   int    iWhichTerm;
   char   szVarModChar[MAX_VARMOD_AA];

   VarMods()
   {
      iBinaryMod = 0;
      bRequireThisMod = 0;
      iMaxNumVarModAAPerMod = 0;
      iVarModTermDistance = -1;
      iWhichTerm = 0;
      dVarModMass = 0.0;
      szVarModChar[0] = '\0';
   }

   VarMods(const VarMods& a)
   {
      iBinaryMod = a.iBinaryMod;
      bRequireThisMod = a.bRequireThisMod;
      iMaxNumVarModAAPerMod = a.iMaxNumVarModAAPerMod;
      iVarModTermDistance = a.iVarModTermDistance;
      iWhichTerm = a.iWhichTerm;
      dVarModMass = a.dVarModMass;
      strcpy(szVarModChar, a.szVarModChar);
   }

   VarMods& operator=(VarMods& a)
   {
      iBinaryMod = a.iBinaryMod;
      bRequireThisMod = a.bRequireThisMod;
      iMaxNumVarModAAPerMod = a.iMaxNumVarModAAPerMod;
      iVarModTermDistance = a.iVarModTermDistance;
      iWhichTerm = a.iWhichTerm;
      dVarModMass = a.dVarModMass;
      strcpy(szVarModChar, a.szVarModChar);

      return *this;
   }
};

struct EnzymeInfo
{
   int  iAllowedMissedCleavage;

   int  iSearchEnzymeOffSet;
   char szSearchEnzymeName[ENZYME_NAME_LEN];
   char szSearchEnzymeBreakAA[MAX_ENZYME_AA];
   char szSearchEnzymeNoBreakAA[MAX_ENZYME_AA];

   int  iSampleEnzymeOffSet;
   char szSampleEnzymeName[ENZYME_NAME_LEN];
   char szSampleEnzymeBreakAA[MAX_ENZYME_AA];
   char szSampleEnzymeNoBreakAA[MAX_ENZYME_AA];

   EnzymeInfo()
   {
      iAllowedMissedCleavage = 0;
      iSearchEnzymeOffSet = 0;
      iSampleEnzymeOffSet = 0;

      szSearchEnzymeName[0] = '\0';
      szSearchEnzymeBreakAA[0] = '\0';
      szSearchEnzymeNoBreakAA[0] = '\0';

      szSampleEnzymeName[0] = '\0';
      szSampleEnzymeBreakAA[0] = '\0';
      szSampleEnzymeNoBreakAA[0] = '\0';
   }

   EnzymeInfo(const EnzymeInfo& a)
   {
      iAllowedMissedCleavage = a.iAllowedMissedCleavage;
      iSearchEnzymeOffSet = a.iSearchEnzymeOffSet;
      iSampleEnzymeOffSet = a.iSampleEnzymeOffSet;

      int i;

      for (i = 0; i < ENZYME_NAME_LEN; i++)
      {
         szSearchEnzymeName[i] = a.szSearchEnzymeName[i];
         szSampleEnzymeName[i] = a.szSampleEnzymeName[i];
      }

      for (i = 0; i < MAX_ENZYME_AA; i++)
      {
         szSearchEnzymeBreakAA[i] = a.szSearchEnzymeBreakAA[i];
         szSearchEnzymeNoBreakAA[i] = a.szSearchEnzymeNoBreakAA[i];
         szSampleEnzymeBreakAA[i] = a.szSampleEnzymeBreakAA[i];
         szSampleEnzymeNoBreakAA[i] = a.szSampleEnzymeNoBreakAA[i];
      }
   }

   EnzymeInfo& operator=(EnzymeInfo& a)
   {
      iAllowedMissedCleavage = a.iAllowedMissedCleavage;
      iSearchEnzymeOffSet = a.iSearchEnzymeOffSet;
      iSampleEnzymeOffSet = a.iSampleEnzymeOffSet;

      int i;

      for (i = 0; i < ENZYME_NAME_LEN; i++)
      {
         szSearchEnzymeName[i] = a.szSearchEnzymeName[i];
         szSampleEnzymeName[i] = a.szSampleEnzymeName[i];
      }

      for (i = 0; i < MAX_ENZYME_AA; i++)
      {
         szSearchEnzymeBreakAA[i] = a.szSearchEnzymeBreakAA[i];
         szSearchEnzymeNoBreakAA[i] = a.szSearchEnzymeNoBreakAA[i];
         szSampleEnzymeBreakAA[i] = a.szSampleEnzymeBreakAA[i];
         szSampleEnzymeNoBreakAA[i] = a.szSampleEnzymeNoBreakAA[i];
      }

      return *this;
   }
};

// *IMPORTANT* If you change this enum, please also change the corresponding
// enum in MangoDataWrapper.h in the MangoWrapper namespace.
enum AnalysisType
{
   AnalysisType_Unknown = 0,
   AnalysisType_DTA,
   AnalysisType_SpecificScan,
   AnalysisType_SpecificScanRange,
   AnalysisType_EntireFile
};

// *IMPORTANT* If you change this enum, please also change the corresponding
// enum in MangoDataWrapper.h in the MangoWrapper namespace.
enum InputType
{
   InputType_UNKNOWN = -1,
   InputType_MS2 = 0,           // ms2, cms2, bms2, etc.
   InputType_MZXML,
   InputType_MZML,
   InputType_RAW,
   InputType_MGF
};

struct InputFileInfo
{
   int  iInputType;
   int  iAnalysisType;
   int  iFirstScan;
   int  iLastScan;
   char szFileName[SIZE_FILE];
   char szBaseName[SIZE_FILE];

   InputFileInfo()
   {
      iInputType = 0;
      iAnalysisType = AnalysisType_Unknown;
      iFirstScan = 0;
      iLastScan = 0;

      szFileName[0] = '\0';
      szBaseName[0] = '\0';
   }

   InputFileInfo(const InputFileInfo& inputObj)
   {
      iInputType = inputObj.iInputType;
      iAnalysisType = inputObj.iAnalysisType;
      iFirstScan = inputObj.iFirstScan;
      iLastScan = inputObj.iLastScan;

      szBaseName[0] = '\0';
      strcpy(szBaseName, inputObj.szBaseName);

      szFileName[0] = '\0';
      strcpy(szFileName, inputObj.szFileName);
   }

   InputFileInfo(char *pszFileName)
   {
      iInputType = 0;
      iAnalysisType = AnalysisType_Unknown;
      iFirstScan = 0;
      iLastScan = 0;

      szBaseName[0] = '\0';

      pszFileName[0] = '\0';
      strcpy(szFileName, pszFileName);
   }

   InputFileInfo& operator = (InputFileInfo &inputObj)
   {
      iInputType = inputObj.iInputType;
      iAnalysisType = inputObj.iAnalysisType;
      iFirstScan = inputObj.iFirstScan;
      iLastScan = inputObj.iLastScan;

      szBaseName[0] = '\0';
      strcpy(szBaseName, inputObj.szBaseName);

      szFileName[0] = '\0';
      strcpy(szFileName, inputObj.szFileName);
      return *this;
   }
};

enum MangoParamType
{
   MangoParamType_Unknown = 0,
   MangoParamType_Int,
   MangoParamType_Double,
   MangoParamType_String,
   MangoParamType_VarMods,
   MangoParamType_DoubleRange,
   MangoParamType_IntRange,
   MangoParamType_EnzymeInfo,
   MangoParamType_DoubleVector
};


// A virtual class that provides a generic data structure to store any Mango
// parameter so that we can store all parameters in one data container
// (e.g. std::map). The specific type of parameter will use the TypedMangoParam
// class which inherits from this class and specifies _paramType and
// _strValue, a string representation of the value of the param

class MangoParam
{
public:
   MangoParam(MangoParamType paramType, const string& strValue)
      : _paramType(paramType), _strValue(strValue) {}
   virtual ~MangoParam() {}
   string& GetStringValue() { return _strValue; }
private:
   MangoParamType _paramType;
   string _strValue;
};


// A templated class to store Mango parameters of any type, specifying the type
// T upon creation. It inherits from MangoParam so after creation, an object of
// this class type can be stored as a MangoParam and cast back to
// TypedMangoParam to access the GetValue() method when needed.

template< typename T >
class TypedMangoParam : public MangoParam
{
public:
   TypedMangoParam (MangoParamType paramType, const string& strValue, const T& value)
      : MangoParam(paramType, strValue), _value(value) {}

   T& GetValue() { return _value; }

private:
   T _value;
};

#endif // _MANGODATA_H_
