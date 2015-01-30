////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////
#include "util/Err.h"
#include "util/Util.h"
//
#include "portability/affy-base-types.h"
//
#include <iostream>
//
using namespace std;

int main(int argc, char *argv[]) {
  int megaByte = 1048576;
  uint64_t free=0, total=0, swapAvail=0, memAvail=0;
  bool retVal = Util::memInfo(free, total, swapAvail, memAvail);
  if(retVal != true) {
    Err::errAbort("Can't determine memory.");
  }
  else { 
    cout << "free: " << free / megaByte << endl;
    cout << "total: " << total / megaByte << endl;
    cout << "swap avail: " << swapAvail / megaByte << endl;
    cout << "mem avail: " << memAvail / megaByte << endl;
  }
  return 0;
}
