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
 * @file   QuantBRLMMTest.cpp
 * @author Chuck Sugnet
 * @date   Wed Feb 22 14:43:23 2006
 * 
 * @brief  Testing out the BRLMM quantification method for classifying genotypes.
 * 
 */
#ifndef QUANTBRLMMTEST_H
#define QUANTBRLMMTEST_H


#include <iostream>
#include <string>
#include <vector>
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include "chipstream/QuantBRLMM.h"
#include "util/Util.h"

using namespace std;

/**
 * @class QuantBRLMMTest
 * @brief cppunit class for testing conversion functions.
 */
class QuantBRLMMTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( QuantBRLMMTest );
  CPPUNIT_TEST( brlmmCallTest );
  CPPUNIT_TEST( computeClusterStatsTest );
  CPPUNIT_TEST( shrinkVarianceTest );
  CPPUNIT_TEST( formatAsBigVarMatTest );
  CPPUNIT_TEST( bayesianMeanUpdateTest );
  CPPUNIT_TEST_SUITE_END();

public:
  // blank test.
  void computeClusterStatsTest();
  void shrinkVarianceTest();
  void formatAsBigVarMatTest();
  void bayesianMeanUpdateTest();
  void brlmmCallTest();
  ClusterPrior *makeClusterPrior();
  void readPriorsFromFile();
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( QuantBRLMMTest );

ClusterPrior *QuantBRLMMTest::makeClusterPrior() {
    double priorCenter[] = {10.13775762,-1.55991681,10.43226400,-0.01907997,10.16756143,1.49258271};
    double priorVars[] = {0.0277327446,0.0792960731,0.0114842020,0.0226200533,0.0256513244,-0.0002765137,0.0270509638,0.0731970281,-0.0100890109};
    double priorCenterVar[] =  {0.4212383466,-0.0547361986,0.4194632311,-0.0185995712,0.3955449182
                               ,0.0414846647,-0.0547361986,0.2451460962,-0.1030744211,0.0578217536
                               ,-0.0509736399,-0.1172055347,0.4194632311,-0.1030744211,0.4536137153
                               ,-0.0001502033,0.4259356205,0.0928455723,-0.0185995712,0.0578217536
                               ,-0.0001502033,0.0742498556,0.0212272364,0.0434582042,0.3955449182
                               ,-0.0509736399,0.4259356205,0.0212272364,0.4308954107,0.0554855210
                               ,0.0414846647,-0.1172055347,0.0928455723,0.0434582042,0.0554855210
                               ,0.1972760895};
    ClusterPrior *prior = new ClusterPrior();
    prior->covars = QuantBRLMM::matrixFromArray(priorCenterVar, ArraySize(priorCenterVar), 6, 6);
    prior->centers = QuantBRLMM::columnVectorFromArray(priorCenter, ArraySize(priorCenter));  
    prior->centerVars = QuantBRLMM::columnVectorFromArray(priorVars, ArraySize(priorVars));  
    return prior;
}
    
void QuantBRLMMTest::readPriorsFromFile() {
}    

void QuantBRLMMTest::computeClusterStatsTest() {
  double aValues[] = {
    11.361066,11.495905,11.499846,9.471878,12.089583,11.518948,12.151334,
    11.398744,9.450799,11.659728,12.066863,8.745842,11.531966,9.844392,
    9.404716,11.658256,11.609871,9.369379,11.315602,12.341380,9.305834,
    11.907416,9.115564,9.272863,11.485880,9.456970,9.152285,11.640832,
    11.321083,11.601214,11.532259,9.510566,9.002815,11.594278,9.691744,
    12.305406,11.732040,9.030667,11.483211,9.220378,8.848936,11.505117,
    11.585855,11.562815,9.084277,9.447083,9.348507,9.276822,11.469744,
    11.566530,11.513086,12.227135,12.172459,12.202828,11.582048,11.352650,
    12.191491,12.094770,11.459227,12.219199,11.574404,11.613651,12.149938,
    11.957283,9.564531,11.192293,11.494406,12.262800,11.575445,10.939653,
    12.034627,11.447445,9.197462,11.539740,9.411299,11.547039,8.896938,
    8.609179,11.750163,9.362601,9.537607,11.506308,11.371014,11.275252,
    11.531089,11.357662,9.343630,11.542983,11.328563,12.310868,11.119849,
    11.608486,11.829604,11.357002,11.330525,9.329796 };
  double bValues[] = {
    11.423694,11.404130,11.492454,11.998379,9.834945,11.667245,8.864805,
    11.644262,12.498001,11.810290,9.256916,11.911130,11.677543,12.053417,
    12.398824,11.642909,11.815944,12.200286,11.602374,10.112048,12.331561,
    11.879698,12.230201,12.095562,11.429825,12.253729,12.376397,11.245672,
    11.379920,11.763835,11.653248,12.333463,12.133431,11.642323,11.828057,
    9.627899,11.471726,12.184844,11.491402,12.442762,12.259773,11.671497,
    11.795796,11.660040,12.246741,12.089583,12.363259,12.365885,11.666978,
    11.854790,11.689649,9.908242,9.293012,9.318317,11.616457,11.212800,
    9.444808,9.366541,11.708006,9.654278,11.732083,11.627853,9.347178,
    9.885696,12.435983,11.258978,11.573174,9.470456,11.706884,11.175300,
    10.176423,11.503130,12.241953,11.975884,12.514640,11.713387,12.049304,
    9.683346,11.913450,12.250150,12.352485,11.631041,11.557751,11.408436,
    11.505812,11.565102,11.820658,11.506407,11.146632,9.480992,11.553389,
    9.530406,9.402159,11.461786,11.345128,12.081683 };
  char calls[] = { 
    1,1,1,2,0,1,0,1,2,1,0,2,1,2,2,1,1,2,1,0,2,1,2,2,1,2,2,1,1,1,1,2,2,1,2,0,1,2,
    1,2,2,1,1,1,2,2,2,2,1,1,1,0,0,0,1,1,0,0,1,0,1,1,0,0,2,1,1,0,1,1,0,1,2,1,2,1,
    2,2,1,2,2,1,1,1,1,1,2,1,1,0,1,0,0,1,1,2 };
  double goldCounts[] = {18, 48, 30};
  double goldMMeans[] = {10.83864375, 11.53349809, 10.70551412};
  double goldAMeans[] = {-2.56894072, 0.09450856, 2.85733757};
  double goldMVars[] = {0.03686030, 0.02599595, 0.11171006 };
  double goldAVars[] = {0.1353275, 0.0240022, 0.1991781 };
  double goldCovars[] = {0.037628698, 0.004908917, 0.086610288};
  double goldVars[] = {0.036860300,0.135327544,0.037628703,0.025995957,0.024002224,0.004908904,0.111710095,0.199178067,0.086610329};
  vector<double> aValuesVec(&aValues[0], &aValues[0]+(ArraySize(aValues)));
  vector<double> bValuesVec(&bValues[0], &bValues[0]+(ArraySize(bValues)));
  for(unsigned int i = 0; i < aValuesVec.size(); i++) {
    aValuesVec[i] = pow(2, aValuesVec[i]);
    bValuesVec[i] = pow(2, bValuesVec[i]);
  }
  vector<affx::GType> callsVec;
  for(unsigned int i = 0; i < ArraySize(calls); i++) {
    callsVec.push_back(affx::GType_from_int(calls[i]));
  }
  QuantBRLMM brlmm(1, 0, 0, 1, QuantBRLMM::MvA, 20, 2, true);
  QuantBRLMM::BrlmmParam param;
  param.m_Transform = QuantBRLMM::MvA;
  QuantBRLMM::transformData(aValuesVec, bValuesVec, param);
  ClusterStats stats = QuantBRLMM::computeClusterStats(aValuesVec, bValuesVec, callsVec, QuantBRLMM::MvA, 1.0);
  bool pass = true;
  for(unsigned int i = 0; i < 3; i++) {
    pass &= Convert::doubleCloseEnough(goldCounts[i], stats.counts[i], 6);
    pass &= Convert::doubleCloseEnough(goldMMeans[i], stats.mMeans[i], 6);
    pass &= Convert::doubleCloseEnough(goldMVars[i], stats.mVars[i], 6);
    pass &= Convert::doubleCloseEnough(goldAMeans[i], stats.aMeans[i], 6);
    pass &= Convert::doubleCloseEnough(goldAVars[i], stats.aVars[i], 6);
    pass &= Convert::doubleCloseEnough(goldCovars[i], stats.covars[i], 6);
  }
  for(unsigned int i = 0; i < stats.vars.size(); i++) {
    pass &= Convert::doubleCloseEnough(goldVars[i], stats.vars[i], 6);
  }
  CPPUNIT_ASSERT(pass);
}


void QuantBRLMMTest::shrinkVarianceTest() {
  int goldCounts[] = {18, 48, 30};
  double empiricalVars[] = {0.036860300,0.135327544,0.037628703,0.025995957,0.024002224,0.004908904,0.111710095,0.199178067,0.086610329};
  double priorVars[] = {0.0277327446,0.0792960731,0.0114842020,0.0226200533,0.0256513244,-0.0002765137,0.0270509638,0.0731970281,-0.0100890109};
  double goldShrinkVars[] = {0.032056323,0.105837296,0.023868439,0.025003044,0.024487254,0.003383781,0.077846442,0.148785652,0.047930593};
  vector<int> counts(&goldCounts[0], &goldCounts[0]+ArraySize(goldCounts));
  vector<double> eVars(&empiricalVars[0], &empiricalVars[0]+ArraySize(empiricalVars));
  vector<double> pVars(&priorVars[0], &priorVars[0]+ArraySize(priorVars));
  
  ClusterStats stats;
  ClusterPrior prior;
  stats.counts = counts;
  stats.vars = eVars;
  prior.centerVars.ReSize(pVars.size());
  for(unsigned int i = 0; i < pVars.size(); i++)
    prior.centerVars.element(i) = pVars[i];

  ColumnVector result = QuantBRLMM::shrinkVariance(stats, prior, 20);
  bool pass = true;
  for(unsigned int i = 0; i < ArraySize(goldShrinkVars); i++) {
    pass &= Convert::doubleCloseEnough(result.element(i), goldShrinkVars[i], 5);
  }
  CPPUNIT_ASSERT(pass);
}

void QuantBRLMMTest::formatAsBigVarMatTest() {
  double shrinkVars[] = {0.032056323,0.105837296,0.023868439,0.025003044,0.024487254,0.003383781,0.077846442,0.148785652,0.047930593};
  ColumnVector x = QuantBRLMM::columnVectorFromArray(shrinkVars, ArraySize(shrinkVars));
  Matrix m = QuantBRLMM::formatAsBigVarMat(x);
  bool pass = true;
  pass &= Convert::doubleCloseEnough(m.element(0,0), 0.03205632);
  pass &= Convert::doubleCloseEnough(m.element(0,1), 0.02386844);
  pass &= Convert::doubleCloseEnough(m.element(1,0), 0.02386844);

  pass &= Convert::doubleCloseEnough(m.element(1,1), 0.10583730);

  pass &= Convert::doubleCloseEnough(m.element(2,2), 0.025003044);
  pass &= Convert::doubleCloseEnough(m.element(2,3), 0.003383781);
  pass &= Convert::doubleCloseEnough(m.element(3,2), 0.003383781);

  pass &= Convert::doubleCloseEnough(m.element(3,3), 0.024487254);
  pass &= Convert::doubleCloseEnough(m.element(4,4), 0.07784644);
  pass &= Convert::doubleCloseEnough(m.element(4,5), 0.04793059);
  pass &= Convert::doubleCloseEnough(m.element(5,4), 0.04793059);
  pass &= Convert::doubleCloseEnough(m.element(5,5), 0.14878565);

  /* Make sure other things are zero... */
  pass &= Convert::doubleCloseEnough(m.element(1,5), 0.0);
  pass &= Convert::doubleCloseEnough(m.element(3,1), 0.0);
  pass &= Convert::doubleCloseEnough(m.element(1,4), 0.0);
  pass &= Convert::doubleCloseEnough(m.element(1,2), 0.0);
  
  CPPUNIT_ASSERT(pass);
}

void QuantBRLMMTest::bayesianMeanUpdateTest() {
  double newVar[] = {0.032056323,0.105837296,0.023868439,0.025003044,0.024487254,0.003383781,0.077846442,0.148785652,0.047930593};
  int counts[] = {18,48,30};
  double means[] = {10.83864374,-2.56894060,11.53349809,0.09450861,10.70551410,2.85733747};
  double priorCenterVar[] =  {0.4212383466,-0.0547361986,0.4194632311,-0.0185995712,0.3955449182
                        ,0.0414846647,-0.0547361986,0.2451460962,-0.1030744211,0.0578217536
                        ,-0.0509736399,-0.1172055347,0.4194632311,-0.1030744211,0.4536137153
                        ,-0.0001502033,0.4259356205,0.0928455723,-0.0185995712,0.0578217536
                        ,-0.0001502033,0.0742498556,0.0212272364,0.0434582042,0.3955449182
                        ,-0.0509736399,0.4259356205,0.0212272364,0.4308954107,0.0554855210
                        ,0.0414846647,-0.1172055347,0.0928455723,0.0434582042,0.0554855210
                        ,0.1972760895};
  double priorCenter[] = {10.13775762,-1.55991681,10.43226400,-0.01907997,10.16756143,1.49258271};
  double goldMeans[] = {10.84589,-2.569764,11.52277,0.09218012,10.73853,2.866197};

  vector<int> countsV(&counts[0], &counts[0] + ArraySize(counts));
  ColumnVector newVarV = QuantBRLMM::columnVectorFromArray(newVar, ArraySize(newVar));
  ColumnVector meansV = QuantBRLMM::columnVectorFromArray(means, ArraySize(means));
  Matrix priorCenterVarM = QuantBRLMM::matrixFromArray(priorCenterVar, ArraySize(priorCenterVar), 6, 6);
  ColumnVector priorCenterV = QuantBRLMM::columnVectorFromArray(priorCenter, ArraySize(priorCenter));
  
  Matrix result = QuantBRLMM::bayesianMeanUpdate(newVarV, countsV, meansV, priorCenterVarM, priorCenterV);
  bool pass = true;
  for(int i = 0; i < result.Nrows(); i++) {
    pass &= Convert::doubleCloseEnough(result.element(i,0), goldMeans[i],4);
  }
  CPPUNIT_ASSERT(pass);
}

void QuantBRLMMTest::brlmmCallTest() {
  double aValues[] = {
    11.361066,11.495905,11.499846,9.471878,12.089583,11.518948,12.151334,
    11.398744,9.450799,11.659728,12.066863,8.745842,11.531966,9.844392,
    9.404716,11.658256,11.609871,9.369379,11.315602,12.341380,9.305834,
    11.907416,9.115564,9.272863,11.485880,9.456970,9.152285,11.640832,
    11.321083,11.601214,11.532259,9.510566,9.002815,11.594278,9.691744,
    12.305406,11.732040,9.030667,11.483211,9.220378,8.848936,11.505117,
    11.585855,11.562815,9.084277,9.447083,9.348507,9.276822,11.469744,
    11.566530,11.513086,12.227135,12.172459,12.202828,11.582048,11.352650,
    12.191491,12.094770,11.459227,12.219199,11.574404,11.613651,12.149938,
    11.957283,9.564531,11.192293,11.494406,12.262800,11.575445,10.939653,
    12.034627,11.447445,9.197462,11.539740,9.411299,11.547039,8.896938,
    8.609179,11.750163,9.362601,9.537607,11.506308,11.371014,11.275252,
    11.531089,11.357662,9.343630,11.542983,11.328563,12.310868,11.119849,
    11.608486,11.829604,11.357002,11.330525,9.329796 };
  double bValues[] = {
    11.423694,11.404130,11.492454,11.998379,9.834945,11.667245,8.864805,
    11.644262,12.498001,11.810290,9.256916,11.911130,11.677543,12.053417,
    12.398824,11.642909,11.815944,12.200286,11.602374,10.112048,12.331561,
    11.879698,12.230201,12.095562,11.429825,12.253729,12.376397,11.245672,
    11.379920,11.763835,11.653248,12.333463,12.133431,11.642323,11.828057,
    9.627899,11.471726,12.184844,11.491402,12.442762,12.259773,11.671497,
    11.795796,11.660040,12.246741,12.089583,12.363259,12.365885,11.666978,
    11.854790,11.689649,9.908242,9.293012,9.318317,11.616457,11.212800,
    9.444808,9.366541,11.708006,9.654278,11.732083,11.627853,9.347178,
    9.885696,12.435983,11.258978,11.573174,9.470456,11.706884,11.175300,
    10.176423,11.503130,12.241953,11.975884,12.514640,11.713387,12.049304,
    9.683346,11.913450,12.250150,12.352485,11.631041,11.557751,11.408436,
    11.505812,11.565102,11.820658,11.506407,11.146632,9.480992,11.553389,
    9.530406,9.402159,11.461786,11.345128,12.081683 };
  char calls[] = { 
    1,1,1,2,0,1,0,1,2,1,0,2,1,2,2,1,1,2,1,0,2,1,2,2,1,2,2,1,1,1,1,2,2,1,2,0,1,2,
    1,2,2,1,1,1,2,2,2,2,1,1,1,0,0,0,1,1,0,0,1,0,1,1,0,0,2,1,1,0,1,1,0,1,2,1,2,1,
    2,2,1,2,2,1,1,1,1,1,2,1,1,0,1,0,0,1,1,2 };
  double priorCenter[] = {10.13775762,-1.55991681,10.43226400,-0.01907997,10.16756143,1.49258271};
  double priorVars[] = {0.0277327446,0.0792960731,0.0114842020,0.0226200533,0.0256513244,-0.0002765137,0.0270509638,0.0731970281,-0.0100890109};
  double priorCenterVar[] =  {0.4212383466,-0.0547361986,0.4194632311,-0.0185995712,0.3955449182
                        ,0.0414846647,-0.0547361986,0.2451460962,-0.1030744211,0.0578217536
                        ,-0.0509736399,-0.1172055347,0.4194632311,-0.1030744211,0.4536137153
                        ,-0.0001502033,0.4259356205,0.0928455723,-0.0185995712,0.0578217536
                        ,-0.0001502033,0.0742498556,0.0212272364,0.0434582042,0.3955449182
                        ,-0.0509736399,0.4259356205,0.0212272364,0.4308954107,0.0554855210
                        ,0.0414846647,-0.1172055347,0.0928455723,0.0434582042,0.0554855210
                        ,0.1972760895};
  char goldCalls [] = { 
    1,1,1,2,0,1,0,1,2,1,0,2,1,2,2,1,1,2,1,0,2,1,2,2,1,2,2,1,1,1,1,2,2,1,2,0,1,2,
    1,2,2,1,1,1,2,2,2,2,1,1,1,0,0,0,1,1,0,0,1,0,1,1,0,0,2,1,1,0,1,1,0,1,2,1,2,1,
    2,2,1,2,2,1,1,1,1,1,2,1,1,0,1,0,0,1,1,2};
  double goldConfs [] = {
    0.0104627971,0.0253592842,0.0066143650,0.0032563085,0.0044687386,
    0.0040885965,0.0198213663,0.0130856275,0.0020924352,0.0252604580,
    0.0044293554,0.0108611933,0.0048713903,0.0269553792,0.0010063871,
    0.0204642714,0.0235970665,0.0001830164,0.0243919767,0.0207208298,
    0.0005146705,0.0964470691,0.0019625785,0.0001120247,0.0163574147,
    0.0010511073,0.0025405105,0.2151426645,0.0181641688,0.0154759806,
    0.0030578748,0.0019697654,0.0038046979,0.0075840987,0.0215393117,
    0.0033092216,0.1094023516,0.0034193120,0.0049465483,0.0022498087,
    0.0092662183,0.0048985680,0.0201268184,0.0046377524,0.0026663130,
    0.0015927753,0.0006390624,0.0009174543,0.0067404200,0.0332148568,
    0.0065785811,0.0067987384,0.0035323600,0.0033875779,0.0065589449,
    0.0704277272,0.0011206503,0.0018101697,0.0141033883,0.0011668982,
    0.0108280477,0.0113290184,0.0021779722,0.0122898107,0.0033562668,
    0.0525212042,0.0002117166,0.0023778657,0.0083012119,0.1527278470,
    0.0329797787,0.0019608318,0.0008282237,0.0716562181,0.0020555060,
    0.0084928031,0.0062930866,0.1284246114,0.0509312245,0.0001831626,
    0.0025419431,0.0016374792,0.0079845593,0.0214020640,0.0092492419,
    0.0107873157,0.0034673627,0.0113463379,0.1016985406,0.0038423193,
    0.0991572186,0.0401892269,0.0121096908,0.0080807419,0.0233845699,
    0.0002629919};
  
  vector<double> aValuesVec(&aValues[0], &aValues[0]+(ArraySize(aValues)));
  vector<double> bValuesVec(&bValues[0], &bValues[0]+(ArraySize(bValues)));
  vector<affx::GType> callsVec;
  for(unsigned int i = 0; i < ArraySize(calls); i++) {
    callsVec.push_back(affx::GType_from_int(calls[i]));
  }
  vector<vector<double> > distances;
  for(unsigned int i = 0; i < aValuesVec.size(); i++) {
    aValuesVec[i] = pow(2, aValuesVec[i]);
    bValuesVec[i] = pow(2, bValuesVec[i]);
  }
  ClusterPrior prior;
  prior.covars = QuantBRLMM::matrixFromArray(priorCenterVar, ArraySize(priorCenterVar), 6, 6);
  prior.centers = QuantBRLMM::columnVectorFromArray(priorCenter, ArraySize(priorCenter));  
  prior.centerVars = QuantBRLMM::columnVectorFromArray(priorVars, ArraySize(priorVars));  
  QuantBRLMM::BrlmmParam param;
  param.m_Transform = QuantBRLMM::MvA;
  QuantBRLMM::transformData(aValuesVec, bValuesVec, param);
  ClusterModel model = QuantBRLMM::fitModel(aValuesVec, bValuesVec, distances, callsVec, prior, "", param);
                                            
  vector<affx::GType> genotypes;
  vector<double> confidences;
  QuantBRLMM::brlmmCall(aValuesVec, bValuesVec, model, genotypes, confidences, distances, param);
  bool callPass = true;

  CPPUNIT_ASSERT(genotypes.size() == ArraySize(goldCalls));
  for(unsigned int i = 0; i < genotypes.size(); i++) {
    callPass &= genotypes[i] == goldCalls[i];
  }
  CPPUNIT_ASSERT(callPass);

  bool confPass = true;
  CPPUNIT_ASSERT(confidences.size() == ArraySize(goldConfs));
  for(unsigned int i = 0; i < confidences.size(); i++) 
    confPass &= Convert::doubleCloseEnough(confidences[i], goldConfs[i], 5);
  CPPUNIT_ASSERT(confPass);
}


#endif
