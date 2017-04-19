/*
   Copyright 2013 University of Washington

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

#ifndef _MANGOSEARCHMANAGER_H_
#define _MANGOSEARCHMANAGER_H_

#include "mango_Data.h"
#include "mango_Interfaces.h"

using namespace MangoInterfaces;

class MangoSearchManager : public IMangoSearchManager
{
public:
   MangoSearchManager();
   ~MangoSearchManager();

   std::map<std::string, MangoParam*>& GetParamsMap();

   // Methods inherited from IMangoSearchManager
   virtual bool DoSearch();
   virtual void AddInputFiles(vector<InputFileInfo*> &pvInputFiles);
   virtual void SetOutputFileBaseName(const char *pszBaseName);
   virtual void SetParam(const string &name, const string &strValue, const string &value);
   virtual bool GetParamValue(const string &name, string &value);
   virtual void SetParam(const string &name, const string &strValue, const int &value);
   virtual bool GetParamValue(const string &name, int &value);
   virtual void SetParam(const string &name, const string &strValue, const double &value);
   virtual bool GetParamValue(const string &name, double &value);
   virtual void SetParam(const string &name, const string &strValue, const VarMods &value);
   virtual bool GetParamValue(const string &name, VarMods &value);
   virtual void SetParam(const string &name, const string &strValue, const DoubleRange &value);
   virtual bool GetParamValue(const string &name, DoubleRange &value);
   virtual void SetParam(const string &name, const string &strValue, const IntRange &value);
   virtual bool GetParamValue(const string &name, IntRange &value);
   virtual void SetParam(const string &name, const string &strValue, const EnzymeInfo &value);
   virtual bool GetParamValue(const string &name, EnzymeInfo &value);
   virtual void SetParam(const string &name, const string &strValue, const vector<double> &value);
   virtual bool GetParamValue(const string &name, vector<double> &value);
   virtual bool IsValidMangoVersion(const string &version);


private:
   bool InitializeStaticParams();

   std::map<std::string, MangoParam*> _mapStaticParams;

   void READ_MZXMLSCANS(char *szMZXML);
   void READ_HK1(char *szHK);
   void READ_HK2(char *szHK);
   void GENERATE_HK(char *szHK);
   int WithinTolerance(double dMass1,
                       double dMass2,
                       double dPPM);
};

#endif
