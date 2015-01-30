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

/**
 * @file   AnalysisStreamTest.cpp
 * @author Chuck Sugnet
 * @date   Fri Nov 10 15:19:44 2006
 * 
 * @brief  
 * 
 * 
 */

#ifndef ANALYSISSTREAMTEST_H
#define ANALYSISSTREAMTEST_H
#include <iostream>
#include <string>
#include <vector>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "algorithm/spectclust/SpectClust.h"
#include "util/Convert.h"
#include "util/Util.h"
#include "newmatio.h"

using namespace std;

class AnalysisStreamTest : public CppUnit::TestFixture {

public:
  CPPUNIT_TEST_SUITE( AnalysisStreamTest );
  CPPUNIT_TEST( testCor );
  CPPUNIT_TEST( testCov );
  CPPUNIT_TEST( testMaxEigen );
  CPPUNIT_TEST_SUITE_END();

  void testCov();
  void testCor();
  void testMaxEigen();
  Matrix getMatrix();
};

CPPUNIT_TEST_SUITE_REGISTRATION( AnalysisStreamTest );

Matrix AnalysisStreamTest::getMatrix() {
  Real col1[] = {4.2833088, -1.6252153, -2.8697821};
  Real col2[] = {-0.2283135, -0.4979154, 0.9406485};
  Real col3[] = {-1.6465658, 1.4138764, 0.2557616};
  Matrix M(3,3);
  M.Column(1) << col1;
  M.Column(2) << col2;
  M.Column(3) << col3;
  return M;
}

void AnalysisStreamTest::testCov() {
  Real col1[] = {14.604386, -1.42652183, -5.04147784};
  Real col2[] = {-1.426522, 0.58477058, -0.04456243};
  Real col3[] = {-5.041478, -0.04456243, 2.38773105};
  Matrix M = getMatrix();
  Matrix Gold(M.Ncols(), M.Ncols());
  Gold.Column(1) << col1;
  Gold.Column(2) << col2;
  Gold.Column(3) << col3;
  cout << "Gold: " << endl << Gold << endl;
  Matrix Test(M.Ncols(), M.Ncols());
  SpectClust::MatrixCov(M, Test);
  cout << "Gold: " << endl << Gold << endl;
  cout << "Test: " << endl << Test << endl;
  Matrix Diff = Gold - Test;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
}

void AnalysisStreamTest::testCor() {
  Real col1[] = { 1.0000000, -0.48813958, -0.85373583 };
  Real col2[] = {-0.4881396,  1.00000000, -0.03771236};
  Real col3[] = {-0.8537358, -0.03771236,  1.00000000};
  Matrix M = getMatrix();
  Matrix Gold(M.Ncols(), M.Ncols());
  Gold.Column(1) << col1;
  Gold.Column(2) << col2;
  Gold.Column(3) << col3;
  cout << "Gold: " << endl << Gold << endl;
  Matrix Test(M.Ncols(), M.Ncols());
  SpectClust::MatrixCor(M, Test);
  cout << "Gold: " << endl << Gold << endl;
  cout << "Test: " << endl << Test << endl;
  Matrix Diff = Gold - Test;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
}

void AnalysisStreamTest::testMaxEigen() {
  Matrix M = getMatrix();
  Matrix C(M.Ncols(), M.Ncols());
  SpectClust::MatrixCov(M, C);
  ColumnVector T(M.Ncols());
  double eVal = 0;
  SpectClust::MaxEigen(C, eVal, T);
  ColumnVector Gold(M.Ncols());
  double goldVal = 1.652681e+01;
  Real goldVec[] = {0.9387423, -0.0830654, -0.3344593};
  Gold.Column(1) << goldVec;
  cout << "Gold: " << goldVal << endl << Gold << endl;
  cout << "Test: " << eVal << endl << T << endl;
  ColumnVector Diff = T - Gold;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > 0.00001 || !Convert::doubleCloseEnough(eVal, goldVal, 4)) {
    CPPUNIT_ASSERT(false);
  }
}

#endif 
