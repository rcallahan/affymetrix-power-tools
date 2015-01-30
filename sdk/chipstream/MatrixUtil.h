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

#ifndef _MATRIXUTIL_H_
#define _MATRIXUTIL_H_

//
#include "newmat.h"
//
#include <cassert>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
//


////
////
////

/** Conversion utility for transforming vectors to column vectors. */
template <typename T1>
ColumnVector columnVectorFromVector(const std::vector<T1>& v) {
  ColumnVector x;
  x.ReSize(v.size());
  for(unsigned int i = 0; i < v.size(); i++) {
    x.element(i) = v[i];
  }
  return x;
}

/** Conversion utility for transforming arrays to column vectors. */
template <typename T1>
ColumnVector columnVectorFromArray(const T1 * const array, const int count) {
  ColumnVector x;
  x.ReSize(count);
  for(unsigned int i = 0; i < count; i++) {
    x.element(i) = array[i];
  }
  return x;
}

/** utility: print out newmat column vector */
void printColumnVector(ColumnVector &v, std::ostream *out, const std::string& delim);

/** utility: */
std::string stringColumnVector(const ColumnVector &v, const std::string& delim);

////
////
////

void printVec(const std::string& lbl,const std::vector<char  > vec);
void printVec(const std::string& lbl,const std::vector<unsigned char> vec);
void printVec(const std::string& lbl,const std::vector<int   > vec);
void printVec(const std::string& lbl,const std::vector<float > vec);
void printVec(const std::string& lbl,const std::vector<double> vec);

////
////
////

/** utility: */

/** Conversion utility for transforming arrays to matrices. */
template <typename T1>
Matrix matrixFromArray(const T1 * const array, int count, int nRow, int nCol) {
  Matrix m(nRow, nCol);
  assert(nRow * nCol == count);
  for(int i = 0; i < count; i++) {
    int col = i / nRow;
    int row = i - col * nRow;
    m.element(row, col) = array[i];
  }
  return m;
}

/** utility: */
Matrix elementMultiplication(Matrix &x, Matrix &y);

/** utility: print out newmat Matrix */
void printMatrix(Matrix &m, std::ostream *out, const std::string& delim);

/** utility: print out newmat Matrix */
std::string stringMatrix(Matrix &m, const std::string& delim);

#endif /* _MATRIXUTIL_H_ */
