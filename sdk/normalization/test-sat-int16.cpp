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

// $Id: test-sat-int16.cpp,v 1.9 2009-09-18 03:37:30 mspald Exp $

/// @file   test-sat-int16.cpp
/// @brief  test functions of the sat_int16 datatype.

#include "normalization/sat_uint16.h"
//
#include <cassert>
#include <cstdio>
#include <ctime>
#include <vector>
//

time_t bench_ts;
time_t bench_te;

#define BENCH_START() time(&bench_ts)
#define BENCH_STOP()  time(&bench_te); printf("TIME=%6d\n",(int)(bench_te-bench_ts))

#define PRINTF_ISI() printf("i=%d si=%d\n",int(i),int(si))

int main(int argc,char *argv[]) {

  sat_uint16 si=0;
  int i,j,k;

  printf("sizeof(sat_uint16)=%d\n",(int)sizeof(sat_uint16));
  assert(sizeof(sat_uint16)==2); // crash otherwise

  si=i=-1;
  PRINTF_ISI();
  si=i=SAT_UINT16_MAX+1;
  PRINTF_ISI();

#define TEST_STEPS 10
#define TEST_INC   (int(SAT_UINT16_MAX*1.2/TEST_STEPS))

  i=0;
  for (j=0;j<TEST_STEPS;j++) {
    i+=TEST_INC;
    si=i;
    PRINTF_ISI();
  }

  // check that ++/-- work
  si=0;
  si++;
  assert(si==1);
  si--;
  assert(si==0);

#define BENCH_REPS (100000)
#define BENCH_MAX  (int(SAT_UINT16_MAX*1.5/5))

  printf("Native ints: ");
  BENCH_START();
  for (k=0;k<BENCH_REPS;k++) {
    j=0;
    j++;
    j--;
    for (i=0;i<BENCH_MAX;i++) {
      j=j+1;
    }
  }
  BENCH_STOP();
  int secs_int16=bench_te-bench_ts;

  printf("sat ints: ");
  BENCH_START();
  for (k=0;k<BENCH_REPS;k++) {
    si=0;
    si++;
    si--;
    for (i=0;i<BENCH_MAX;i++) {
      si=j+1;
    }
  }
  BENCH_STOP();
  int secs_sat_uint16=bench_te-bench_ts;

  printf("*** Slowdown = %f\n",double(secs_sat_uint16)/secs_int16);
}
