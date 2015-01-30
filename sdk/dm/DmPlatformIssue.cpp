////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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


#include "dm/DM.h"
//
#include "portability/affy-base-types.h"
//
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/* Template for printing out a generic vector or array. */
template<class T> inline void PrintDmExample(const std::string &descript, const std::string &name, T start, T end) {
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(20);
  std::cout << name << "\t" << descript;
  for(; start != end; ++start) {

    std::cout << "\t" << *start;
  }
  std::cout << std::endl;
}

/** Little template function to make string conversion easy. 
    this isn't the fastest way to do things, but it is easy. */
template <class T> 
std::string ToDmStr(const T &t) {
  std::ostringstream s;
  s.setf(std::ios_base::boolalpha);
  s << t; 
  return s.str();
};

int main(int argc, char *argv[]) {
  std::vector<CQuartet> qVec(6);
  for(uint32_t qIx = 0; qIx < 6; qIx++) {
    for(uint32_t pIx = 0; pIx < 4; pIx++) {
      qVec[qIx].pixels[pIx] = 9;
    }
  }
  int pos = 0;
  int qCount = 0;
  qVec[qCount].intensity[pos++] = 585.00000000000000000000; 
  qVec[qCount].intensity[pos++] = 230.00000000000000000000;
  qVec[qCount].intensity[pos++] = 616.00000000000000000000;
  qVec[qCount].intensity[pos++] = 274.00000000000000000000;
  pos = 0;
  qVec[qCount].variance[pos++] = 3893.76025390625000000000;
  qVec[qCount].variance[pos++] = 219.04000854492187500000;
  qVec[qCount].variance[pos++] = 6021.75976562500000000000;
  qVec[qCount].variance[pos++] = 1892.25000000000000000000;
  pos = 0;
  qCount++;
  qVec[qCount].intensity[pos++] = 377.00000000000000000000;
  qVec[qCount].intensity[pos++] = 266.00000000000000000000;
  qVec[qCount].intensity[pos++] = 741.00000000000000000000;
  qVec[qCount].intensity[pos++] = 284.00000000000000000000;
  pos = 0;
  qVec[qCount].variance[pos++] = 967.21002197265625000000;
  qVec[qCount].variance[pos++] = 146.41000366210937500000;
  qVec[qCount].variance[pos++] = 17582.76171875000000000000;
  qVec[qCount].variance[pos++] = 691.68994140625000000000;
  pos = 0;
  qCount++;
  qVec[qCount].intensity[pos++] = 623.00000000000000000000;
  qVec[qCount].intensity[pos++] = 384.00000000000000000000;
  qVec[qCount].intensity[pos++] = 1062.00000000000000000000;
  qVec[qCount].intensity[pos++] = 449.00000000000000000000;
  pos = 0;
  qVec[qCount].variance[pos++] = 10342.88964843750000000000;
  qVec[qCount].variance[pos++] = 1989.15991210937500000000;
  qVec[qCount].variance[pos++] = 9623.60937500000000000000;
  qVec[qCount].variance[pos++] = 2959.36010742187500000000;
  pos = 0;
  qCount++;
  qVec[qCount].intensity[pos++] = 382.00000000000000000000;
  qVec[qCount].intensity[pos++] = 310.00000000000000000000;
  qVec[qCount].intensity[pos++] = 1048.00000000000000000000;
  qVec[qCount].intensity[pos++] = 358.00000000000000000000;
  pos = 0;
  qVec[qCount].variance[pos++] = 2735.29003906250000000000;
  qVec[qCount].variance[pos++] = 2672.89013671875000000000;
  qVec[qCount].variance[pos++] = 20534.89062500000000000000;
  qVec[qCount].variance[pos++] = 2470.09008789062500000000;
  pos = 0;
  qCount++;
  qVec[qCount].intensity[pos++] = 374.00000000000000000000;
  qVec[qCount].intensity[pos++] = 177.00000000000000000000;
  qVec[qCount].intensity[pos++] = 1592.00000000000000000000;
  qVec[qCount].intensity[pos++] = 379.00000000000000000000;
  pos = 0;
  qVec[qCount].variance[pos++] = 2777.29003906250000000000;
  qVec[qCount].variance[pos++] = 595.35998535156250000000;
  qVec[qCount].variance[pos++] = 63705.75781250000000000000;
  qVec[qCount].variance[pos++] = 1980.25000000000000000000;
  pos = 0;
  qCount++;
  qVec[qCount].intensity[pos++] = 385.00000000000000000000;
  qVec[qCount].intensity[pos++] = 333.00000000000000000000;
  qVec[qCount].intensity[pos++] = 1290.00000000000000000000;
  qVec[qCount].intensity[pos++] = 676.00000000000000000000;
  pos = 0;
  qVec[qCount].variance[pos++] = 729.00000000000000000000;
  qVec[qCount].variance[pos++] = 784.00000000000000000000;
  qVec[qCount].variance[pos++] = 17662.40820312500000000000;
  qVec[qCount].variance[pos++] = 5700.25000000000000000000;

  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(20);

  std::cout << "Data In:" << std::endl;
  for(uint32_t qIx = 0; qIx < qVec.size(); qIx++) {
    PrintDmExample(ToDmStr("pixels-") + ToDmStr(qIx), "SNP_A-2245218", &qVec[qIx].pixels[0], &qVec[qIx].pixels[NUMCELLS]);
    PrintDmExample(ToDmStr("intensity-") + ToDmStr(qIx), "SNP_A-2245218", &qVec[qIx].intensity[0], &qVec[qIx].intensity[NUMCELLS]);
    PrintDmExample(ToDmStr("variance-") + ToDmStr(qIx), "SNP_A-2245218", &qVec[qIx].variance[0], &qVec[qIx].variance[NUMCELLS]);
  }
  std::pair<float,int> call;
  CallDM(qVec, call, 0.0);
  std::cout << "SNP_A-2245218" << "\t" << call.first << "\t" << call.second << std::endl;
  return 0;
}
