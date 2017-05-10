/*
   Copyright 2017 University of Washington                          3-clause BSD license

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>
#include <string>
#include <ctime>

#ifdef _WIN32
#include <direct.h>
#include <errno.h>
#define STRCMP_IGNORE_CASE(a,b) _strcmpi(a,b)
#else
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <pthread.h>
#define STRCMP_IGNORE_CASE(a,b) strcasecmp(a,b)
#endif

using namespace std;

#include "MSReader.h"
#include "Spectrum.h"
#include "MSObject.h"
#include <vector>
#include <cfloat>

#define mango_version "2017.01 rev. 0alpha"
#define copyright "(c) University of Washington"

// Redefined how the bin offset is interpreted and applied.  The valid range for the offset is
// now between 0.0 and 1.0 and scales to the binWidth.
#define BIN(dMass) (int)(dMass*g_staticParams.dInverseBinWidth + g_staticParams.dOneMinusBinOffset)

#define isEqual(x, y) (std::abs(x-y) <= ( (std::abs(x) > std::abs(y) ? std::abs(y) : std::abs(x)) * FLT_EPSILON))

using namespace MSToolkit;
#define logout(szString) fputs(szString, stdout)
#define logerr(szString) fputs(szString, stderr)

#endif // _COMMON_H_
