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

#include "chipstream/MatrixUtil.h"
//
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdio.h>
//


////
//// ColumnVectors
////

/** Utility printing function. */
void printColumnVector(ColumnVector &v, std::ostream *out, const std::string& delim) {
  int nRow = v.Nrows();
  int i = 0;
  if(out == NULL) 
    out = &cout;
  for(i = 0; i < nRow - 1; i++)
    (*out) << v.element(i) << delim;
  (*out) << v.element(i);
}

std::string stringColumnVector(ColumnVector &v, const std::string& delim)
{
  std::ostringstream ostm;
  printColumnVector(v,&ostm,delim);
  return ostm.str();
}

////
//// Matrix
////

/** Mutilply each element individually and return matrix. */
Matrix elementMultiplication(Matrix &x, Matrix &y) {
  assert(x.Nrows() == y.Nrows());
  assert(x.Ncols() == y.Ncols());
  Matrix m(x.Nrows(), x.Ncols());
  for(unsigned int i = 0; i < x.Nrows(); i++) {
    for(unsigned int j = 0; j < x.Ncols(); j++) {
      m.element(i,j) = x.element(i,j) * y.element(i,j);
    }
  }
  return m;
}

/** Utility printing function. */
void printMatrix(Matrix &m, std::ostream *out, const std::string& delim) {
  int nRow = m.Nrows();
  int nCol = m.Ncols();
  int i = 0;
  int j = 0;
  if(out == NULL) 
    out = &cout;

  for(i = 0; i < nRow-1; i++) {
    for(j = 0; j < nCol - 1; j++) {
      (*out) << m.element(i,j) << delim;
    }
    (*out) << m.element(i,j) << delim;
  }
  for(j = 0; j < nCol - 1; j++) {
    (*out) << m.element(i,j) << delim;
  }
  (*out) << m.element(i,j);
}

std::string stringMatrix(Matrix &m, const std::string& delim)
{
  std::ostringstream ostm;
  printMatrix(m,&ostm,delim);
  return ostm.str();
}


////
////
////


template <typename T1>
void printVec_tmpl(const std::string& lbl,const std::string& fmt,const std::vector<T1> vec)
{
  int s=vec.size();
  
  printf("%-15s(%4d):",lbl.c_str(),s);
  for (int i=0;i<s;i++) {
    //cout << vec[i];
    if (i!=0) {
      printf(",");
    }
    printf(fmt.c_str(),vec[i]);
  }
  printf("\n");
}

//
void printVec(const std::string& lbl,const std::vector<char> vec)
{
  printVec_tmpl(lbl,"%d",vec);
}
void printVec(const std::string& lbl,const std::vector<unsigned char> vec)
{
  printVec_tmpl(lbl,"%d",vec);
}
void printVec(const std::string& lbl,const std::vector<int> vec)
{
  printVec_tmpl(lbl,"%d",vec);
}
void printVec(const std::string& lbl,const std::vector<float> vec)
{
  printVec_tmpl(lbl,"%.10f",vec);
}
void printVec(const std::string& lbl,const std::vector<double> vec)
{
  printVec_tmpl(lbl,"%.12f",vec);
}
