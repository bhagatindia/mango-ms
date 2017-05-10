/*
   Copyright 2017 University of Washington                          3-clause BSD license

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _MANGOINTERFACES_H_
#define _MANGOINTERFACES_H_

#include "Common.h"
#include "mango_Data.h"

using namespace std;

namespace MangoInterfaces
{
   class IMangoSearchManager
   {
public:
      virtual ~IMangoSearchManager() {}
      virtual bool DoSearch() = 0;
      virtual void AddInputFiles(vector<InputFileInfo*> &pvInputFiles) = 0;
      virtual void SetOutputFileBaseName(const char *pszBaseName) = 0;
      virtual void SetParam(const string &name, const string &strValue, const string &value) = 0;
      virtual bool GetParamValue(const string &name, string &value) = 0;
      virtual void SetParam(const string &name, const string &strValue, const int &value) = 0;
      virtual bool GetParamValue(const string &name, int &value) = 0;
      virtual void SetParam(const string &name, const string &strValue, const double &value) = 0;
      virtual bool GetParamValue(const string &name, double &value) = 0;
      virtual void SetParam(const string &name, const string &strValue, const VarMods &value) = 0;
      virtual bool GetParamValue(const string &name, VarMods &value) = 0;
      virtual void SetParam(const string &name, const string &strValue, const DoubleRange &value) = 0;
      virtual bool GetParamValue(const string &name, DoubleRange &value) = 0;
      virtual void SetParam(const string &name, const string &strValue, const IntRange &value) = 0;
      virtual bool GetParamValue(const string &name, IntRange &value) = 0;
      virtual void SetParam(const string &name, const string &strValue, const EnzymeInfo &value) = 0;
      virtual bool GetParamValue(const string &name, EnzymeInfo &value) = 0;
      virtual void SetParam(const string &name, const string &strValue, const vector<double> &value) = 0;
      virtual bool GetParamValue(const string &name, vector<double> &value) = 0;
      virtual bool IsValidMangoVersion(const string &version) = 0;
   };

   IMangoSearchManager *GetMangoSearchManager();
   void ReleaseMangoSearchManager();
}

#endif // _MANGOINTERFACES_H_
