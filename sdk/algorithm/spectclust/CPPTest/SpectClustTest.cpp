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
 * @file   SpectClustTest.cpp
 * @author Chuck Sugnet
 * @date   Mon Dec 18 10:04:50 2006
 * 
 * @brief  Routines to test the spectral cluster algorithms.
 */


#ifndef SPECTCLUSTTEST_H
#define SPECTCLUSTTEST_H
#include "algorithm/spectclust/SpectClust.h"
//
#include "util/Convert.h"
#include "util/Util.h"
#include "util/RowFile.h"
//
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <iostream>
#include <string>
#include <vector>
//

#include "newmatio.h" // needs to be after iostream.h in order to compile.

using namespace std;

class SpectClustTest : public CppUnit::TestFixture {

public:
  CPPUNIT_TEST_SUITE( SpectClustTest );
  CPPUNIT_TEST( testNCutDynamicProgram );
  CPPUNIT_TEST( testAngleDist );
  //  CPPUNIT_TEST( testDoCluster );
   CPPUNIT_TEST( testCov );
  CPPUNIT_TEST( testMaxEigen );
  // CPPUNIT_TEST( testFillInDistance );
  CPPUNIT_TEST( testFindNLargestEvals );
  CPPUNIT_TEST_SUITE_END();

  void testNCutDynamicProgram();
  void testCov();
  void testMaxEigen();
  void testFillInDistance();
  void testFindNLargestEvals();
  void testDoCluster();
  void testAngleDist();
  Matrix getMatrix();
  Matrix getClustMatrix();
  void RFileToMatrix(Matrix &M, const char *fileName);
};

void SpectClustTest::RFileToMatrix(Matrix &M, const char *fileName) {
  vector<vector<Real> > data;
  RowFile::matrixFromFile(fileName, data, 1, 1);
  assert(data.size() > 0 && data[0].size() > 0);
  M.ReSize(data.size(), data[0].size());
  for(int i = 0; i < data.size(); i++) {
    M.Row(i+1) << &(data[i][0]);
  }
}

void SpectClustTest::testNCutDynamicProgram() {
  Matrix Dist;
  RFileToMatrix(Dist, "data/spike-in.norm.angleDist.b12.txt");
  bool converged;
  vector<double> eVals;
  Matrix EVec;
  converged = SpectClust::findNLargestEvals(Dist, 2, eVals, EVec, 200);
  vector<double> cutVals;
  std::vector<std::pair<double,int> > indices;
  for(int i = 0; i < Dist.Ncols(); i++) {
    std::pair<double,int> p;
    p.first = EVec.element(i,1);// / E.element(i,0) ;
    p.second = i;
    indices.push_back(p);
  }
  std::sort(indices.begin(), indices.end(), SpectClust::PairLess());
  SpectClust::calcNormalizedCutValues(Dist, indices, cutVals);
  for(int i = 0; i < cutVals.size(); i++) {
    double nCut = SpectClust::calcNormalizedCut(Dist, indices, i);
    double diff = nCut - cutVals[i];
    //    cout << "\t" << nCut << "\t" << cutVals[i] << "\t" << diff << endl;
    if(fabs(diff) > .000001) {
      CPPUNIT_ASSERT(false);
    }
  }
}

void SpectClustTest::testAngleDist() {
  Matrix M, Gold;
  SymmetricMatrix Dist;
  RFileToMatrix(M, "data/spike-in.b12.txt");
  RFileToMatrix(Gold, "data/spike-in.norm.angleDist.b12.0-2.txt");
  SpectClust::rowMedianDivide(M);
  M = M.t();
  AngleMetric metric(M);
  SpectClust::fillInDistance(Dist, M, metric, false);
  Matrix Diff = Gold - Dist;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
  Matrix GoldVec;
  RFileToMatrix(GoldVec, "data/spike-in.norm.angleDist.eVec.b12.0-2.txt");
  SpectClust::normalizeSum(Dist);
  bool converged = false;
  vector<double> eVals;
  Matrix EVec;
  converged = SpectClust::findNLargestEvals(Dist, 2, eVals, EVec, 200);
  // just so happens that our eigen vector finder gives a different
  // sign than R's.
  GoldVec.Column(1) = GoldVec.Column(1) * -1;
  //  GoldVec.Column(2) = GoldVec.Column(2) * -1;
  Diff = EVec - GoldVec;
  maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
  vector<int> clusters;
  SpectClust::partitionClusters(Dist, EVec, eVals, 2, clusters, 4, 0, 1);
  vector<int> pos, neg;
  for(int i = 0; i < clusters.size(); i++) {
    if(clusters[i] == 0) {
      neg.push_back(i);
    }
    else if(clusters[i] == 1) {
      pos.push_back(i);
    }
    else {
      Err::errAbort("Only expecting 0 or 1");
    }
  }
  Matrix Pos(pos.size(), 1), Neg(neg.size(),1);
  for(int i = 0; i < pos.size(); i++) {
    Pos.element(i,0) = pos[i] + 1;
  }
  for(int i = 0; i < neg.size(); i++) {
    Neg.element(i,0) = neg[i] + 1;
  }
  Matrix GoldPos, GoldNeg;
  RFileToMatrix(GoldPos, "data/spike-in.norm.angleDist.pos.b12.txt");
  RFileToMatrix(GoldNeg, "data/spike-in.norm.angleDist.neg.b12.txt");
  Diff = Pos - GoldPos;
  maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
  Diff = Neg - GoldNeg;
  maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
  vector<double> clusterVals;
  SpectClust::orderIntraClusterCorr(M, clusters, 2, clusterVals);
}

Matrix SpectClustTest::getClustMatrix() {
  // columns 1-3 are in the same cluster 4 and 5 are decoys.
  Real col1[] = {  5.274081,  4.794793,  1.1374331,  1.2308186,  5.707615,  3.206105,  0.8317420};
  Real col2[] = {  4.789248,  5.079418,  0.9303852,  0.3509492,  5.726203,  3.573452,  1.1412697};
  Real col3[] = {  4.937514,  4.889457,  1.6881732,  1.8803931,  6.095611,  2.868552,  1.0375580};
  Real col4[] = { -5.436639, -5.241824, -1.0230441, -1.9654199, -6.482444, -2.529698, -0.9388572};
  Real col5[] = { -5.174804, -4.906279, -2.1491716, -1.1618790, -5.094437, -4.392525, -0.5924768};
  Matrix M(7,5);
  // note I'm swapping col4 and col2 so application doesn't expect them in right order...
  M.Column(1) << col1;
  M.Column(4) << col2;
  M.Column(3) << col3;
  M.Column(2) << col4;
  M.Column(5) << col5;
  return M;
}

Matrix SpectClustTest::getMatrix() {
  Real col1[] = {4.2833088, -1.6252153, -2.8697821};
  Real col2[] = {-0.2283135, -0.4979154, 0.9406485};
  Real col3[] = {-1.6465658, 1.4138764, 0.2557616};
  Matrix M(3,3);
  M.Column(1) << col1;
  M.Column(2) << col2;
  M.Column(3) << col3;
  return M;
}

void SpectClustTest::testCov() {
  Real col1[] = {14.604386, -1.42652183, -5.04147784};
  Real col2[] = {-1.426522, 0.58477058, -0.04456243};
  Real col3[] = {-5.041478, -0.04456243, 2.38773105};
  Matrix M = getMatrix();
  Matrix Gold(M.Ncols(), M.Ncols());
  Gold.Column(1) << col1;
  Gold.Column(2) << col2;
  Gold.Column(3) << col3;
  //  cout << "Gold: " << endl << Gold << endl;
  Matrix Test(M.Ncols(), M.Ncols());
  SpectClust::MatrixCov(M, Test);
  //  cout << "Gold: " << endl << Gold << endl;
  //  cout << "Test: " << endl << Test << endl;
  Matrix Diff = Gold - Test;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
}

void SpectClustTest::testMaxEigen() {
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
  //  cout << "Gold: " << goldVal << endl << Gold << endl;
  //  cout << "Test: " << eVal << endl << T << endl;
  ColumnVector Diff = T - Gold;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > 0.00001 || !Convert::doubleCloseEnough(eVal, goldVal, 4)) {
    CPPUNIT_ASSERT(false);
  }

}

void SpectClustTest::testFindNLargestEvals() {
  Matrix M = getMatrix();
  Matrix C(M.Ncols(), M.Ncols());
  SpectClust::MatrixCov(M, C);
  Matrix Evec;
  vector<double> eVals;
  SpectClust::findNLargestEvals(C, 2, eVals, Evec);
  Matrix Gold(M.Ncols(), 2);
  double goldVal1 = 1.652681e+01, goldVal2 = 1.050075e+00;
  Real goldVec1[] = {0.9387423, -0.0830654, -0.3344593};
  Real goldVec2[] = {0.1960120, -0.6695446, 0.7164422};
  Gold.Column(1) << goldVec1;
  Gold.Column(2) << goldVec2;
//   cout << "Gold: " << goldVal1 << " " << goldVal2 << endl << Gold << endl;
//   cout << "Test: " << eVals[0] << " " << eVals[1] << endl << Evec << endl;
  Matrix Diff = Evec - Gold;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > 0.00001 || !Convert::doubleCloseEnough(eVals[0], goldVal1, 4) ||
     !Convert::doubleCloseEnough(eVals[1], goldVal2, 4)) {
    CPPUNIT_ASSERT(false);
  }

}

void SpectClustTest::testFillInDistance() {
  Real col1[] = {8.478264e-18, 0.026075522, 9.996101e-01};
  Real col2[] = {2.607552e-02, 0.016092326, 1.617114e-03};
  Real col3[] = {9.996101e-01, 0.001617114, 3.517769e-07};
//   double col1[] = {1.000000e-00, 1.243566e-08, 1.721017e-12};
//   double col2[] = {1.243566e-08, 7.706867e-01, 4.068833e-02};
//   double col3[] = {1.721017e-12, 4.068833e-02, 9.923325e-01};
  Matrix M = getMatrix();
  Matrix Gold(M.Ncols(), M.Ncols());
  Gold.Column(1) << col1;
  Gold.Column(2) << col2;
  Gold.Column(3) << col3;
  SymmetricMatrix A(M.Nrows());
  A = 0;
  DistanceMetric dist;
  SpectClust::fillInDistance(A, M, dist);
//   cout << "Orig: " <<  M << endl;
//   cout << "Distance: " << A << endl;
  Matrix Diff = Gold - A;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
}

void SpectClustTest::testDoCluster() {
  Real col1[] = {0.476427,0.401074,0.474626,0.474408,0.402273};
  Real col2[] = {-0.331007 ,0.584172,-0.327181,-0.325684,0.579706};
  Real evalGold[] = {1.000000000,0.726870496};
  Matrix M = getClustMatrix();
  //  SpectClust::printMatrix(M);
  SymmetricMatrix D;
  Matrix EVec;
  Matrix EGold(M.Ncols(), 2);
  EGold.Column(1) << col1;
  EGold.Column(2) << col2;
  std::vector<double> eVals;
  CorrelationMetric dist;
  std::vector<int> clusters;
  SpectClust::fillInDistance(D, M, dist);
  bool converged = false;
  converged = SpectClust::findNLargestEvals(D, 2, eVals, EVec);
  SpectClust::partitionClusters(D, EVec, eVals, 2, clusters, 2, 0, 1);
  Matrix Diff = EGold - EVec;
  double maxDiff = MaximumAbsoluteValue(Diff);
  if(maxDiff > .00001) {
    CPPUNIT_ASSERT(false);
  }
  if(!Convert::doubleCloseEnough(evalGold[0], eVals[0], 4) ||
     !Convert::doubleCloseEnough(evalGold[1], eVals[1], 4)) {
    CPPUNIT_ASSERT(false);
  }
  for(int i = 0; i < clusters.size(); i++) {
    if(i == 0 && clusters[i] != 0) {
      CPPUNIT_ASSERT(false);
    }
    if(i == 1 && clusters[i] != 1) {
      CPPUNIT_ASSERT(false);
    }
    if((i == 2 || i == 3) && clusters[i] != 0) {
      CPPUNIT_ASSERT(false);
    }
    if(i == 4 && clusters[i] != 1) {
      CPPUNIT_ASSERT(false);
    }
  }
}

CPPUNIT_TEST_SUITE_REGISTRATION( SpectClustTest );

/* 
Here is some R code used to get the correct solutions above.

dotDistance <- function(M) {
  D = matrix(nrow=ncol(M),ncol=ncol(M),0);
  for(i in 1:ncol(M)) {
    for(j in 1:ncol(M)) {
      D[i,j] = exp(-1 * (M[,i] %*% M[,j])); 
    }
  }
  return(D);
}

dat = matrix(nrow=3,c(4.2833088, -1.6252153, -2.8697821,-0.2283135, -0.4979154, 0.9406485,-1.6465658, 1.4138764, 0.2557616))

A = dotDistance(dat);
rsum = apply(A, 1, sum);
D = matrix(nrow=nrow(dat),ncol=ncol(dat),0);
diag(D) = rsum;
D = solve(D^{.5})
dist = D %*% A %*% D

val1 <- rnorm(7, mean=3, sd=2)
val2 <- val1 + rnorm(7, mean=0,sd=.5)
val3 <- val1 + rnorm(7, mean=0,sd=.5)
val4 <- (-1 * val1) + rnorm(7, mean=0,sd=.5)
val5 <- (-1 * val1) + rnorm(7, mean=0,sd=.5)
vals <- cbind(val1,val2,val3,val4,val5)

         val1      val2     val3       val4       val5
[1,] 5.274081 4.7892479 4.937514 -5.4366394 -5.1748044
[2,] 4.794793 5.0794183 4.889457 -5.2418239 -4.9062786
[3,] 1.137433 0.9303852 1.688173 -1.0230441 -2.1491716
[4,] 1.230819 0.3509492 1.880393 -1.9654199 -1.1618790
[5,] 5.707615 5.7262031 6.095611 -6.4824439 -5.0944373
[6,] 3.206105 3.5734521 2.868552 -2.5296975 -4.3925248
[7,] 0.831742 1.1412697 1.037558 -0.9388572 -0.5924768

A = exp(-1 * (1 - cor(vals)));
rsum = apply(A, 1, sum);
D = matrix(nrow=nrow(A),ncol=ncol(A),0);
diag(D) = rsum;
D = solve(D^{.5})
dist = D %*% A %*% D

valsX = vals[,c(1,4,2,5,3)]
A = exp(-1 * (1 - cor(valsX)));
rsum = apply(A, 1, sum);
D = matrix(nrow=nrow(A),ncol=ncol(A),0);
diag(D) = rsum;
D = solve(D^{.5})
dist = D %*% A %*% D
*/

#endif 
