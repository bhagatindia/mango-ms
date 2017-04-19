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

#include "mango_Interfaces.h"
#include "mango_SearchManager.h"

using namespace MangoInterfaces;

MangoSearchManager *g_pMangoSearchManager = NULL;

IMangoSearchManager *MangoInterfaces::GetMangoSearchManager()
{
   if (NULL == g_pMangoSearchManager)
   {
      g_pMangoSearchManager = new MangoSearchManager();
   }

   IMangoSearchManager *pMangoSearchMgr = static_cast<IMangoSearchManager*>(g_pMangoSearchManager);
   return pMangoSearchMgr;
}

void MangoInterfaces::ReleaseMangoSearchManager()
{
   if (NULL != g_pMangoSearchManager)
   {
      delete g_pMangoSearchManager;
      g_pMangoSearchManager = NULL;
   }
}


