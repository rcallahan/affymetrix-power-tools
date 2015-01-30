////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

//
// cvs:affy/sdk/stats/affy_random_sample_test.cpp ---
//
// $Id: affy_random_sample_test.cpp,v 1.7 2009-07-31 15:21:02 awilli Exp $
//
// g++ -o affy_random_sample_test{,.cpp}

#include <vector>
#ifdef WIN32
#include <process.h>
#define GETPID _getpid
#else
#include <unistd.h>
#define GETPID getpid
#endif

// switch between the two.
#define USE_AFFY_RANDOM_SAMPLE 1
#ifdef USE_AFFY_RANDOM_SAMPLE
  #include "affy_random_sample.h"
#else
  #include <ext/algorithm>
#endif

using namespace std;

int main() {
  vector<int> vi;
  vector<int> vo;

  // comment this out if you want the same results.
  srand((unsigned int)GETPID());

  // fill
  for (int i=0;i<100;i++) {
    vi.push_back(i);
  }
  vo.resize(10);

  // do
  random_sample(vi.begin(),vi.end(),vo.begin(),vo.end());

  // out
  printf("VO=[");
  for (int i=0;i<vo.size();i++) {
    printf("%d,",vo[i]);
  }
  printf("]\n");

}
