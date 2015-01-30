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
 * @file   QuantExprMethodTest.cpp
 * @author csugne
 * @date   Mon Nov  7 14:07:42 PST 2005
 * 
 * @brief  Testing the QuantExprMethod functions.
 * 
 */
#ifndef QUANTEXPRMETHODTEST_H
#define QUANTEXPRMETHODTEST_H

#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantRma.h"
#include "chipstream/QuantPlier.h"
#include "chipstream/QuantIterPlier.h"
//
#include "util/Convert.h"
#include "util/TableFile.h"
#include "util/Err.h"
//
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <string>
#include <vector>

using namespace std;
/**
 * @class QuantExprMethodTest
 * @brief cppunit class for testing conversion functions.
 */
class QuantExprMethodTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( QuantExprMethodTest );
  CPPUNIT_TEST( testQuantRma );
  CPPUNIT_TEST( testQuantPlier );
  CPPUNIT_TEST( testQuantIterPlier );
  CPPUNIT_TEST_SUITE_END();

public:
  // blank test.
  void testQuantRma();
  void testQuantPlier();
  void testQuantIterPlier();
  bool testQuantExprMethod(QuantExprMethod &qMethod, const char *fileIn, 
                       const char *expectedOut, const char *expectedFeatOut,
                       const char *expectedResOut, bool transpose=false);
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( QuantExprMethodTest );

bool QuantExprMethodTest::testQuantExprMethod(QuantExprMethod &qMethod, const char *fileIn, 
                                      const char *expectedOut, const char *expectedFeatOut,
                                      const char *expectedResOut, bool transpose) {
  TableFile input('\t','#',false,false);
  TableFile expected('\t','#',false,false);
  TableFile expectedFeat('\t','#',false,false);
  TableFile expectedRes('\t','#',false,false);

  input.open(fileIn);
  expected.open(expectedOut);

  int rowIx = 0, colIx = 0;
  /* Load data into quantifier. */
  if(transpose)
    qMethod.setBounds(input.numRows(), input.numCols());
  else
    qMethod.setBounds(input.numCols(), input.numRows());
  for(rowIx = 0; rowIx < input.numRows(); rowIx++) {
    for(colIx = 0; colIx < input.numCols(); colIx++) {
      if(transpose) {
        qMethod.setPMDataAt(rowIx, colIx, Convert::toFloat(input.getData(rowIx, colIx).c_str()));
        qMethod.setMMDataAt(rowIx, colIx, 0);
      }
      else {
        qMethod.setPMDataAt(colIx, rowIx, Convert::toFloat(input.getData(rowIx, colIx).c_str()));
        qMethod.setMMDataAt(colIx, rowIx, 0);
      }
    }
  }
  /* Heavy lifting. */
  qMethod.computeEstimate();

  /* Check residuals */
  if(qMethod.haveResiduals() && expectedResOut != NULL) {
    expectedRes.open(expectedResOut);
    int count = 0;
    for(int probeIx = 0; probeIx < qMethod.getNumFeatures(); probeIx++) {
      if(qMethod.featureUsed(probeIx)) {
        for(int chipIx = 0; chipIx < qMethod.getNumTargets(); chipIx++) {
          float estimated = qMethod.getResidual(probeIx, chipIx);
          float gold = Convert::toFloat(expectedRes.getData(chipIx, count).c_str());
          CPPUNIT_ASSERT(Convert::doubleCloseEnough(estimated, gold,4));
        }
        count++;
      }
    }
  }
  /* Check feature effects */
  if(qMethod.haveFeatureEffects() && expectedFeatOut != NULL) {
    int count = 0;
    expectedFeat.open(expectedFeatOut);   
    for(int probeIx = 0; probeIx < qMethod.getNumFeatures(); probeIx++) {
      if(qMethod.featureUsed(probeIx)) {
        float estimated = qMethod.getFeatureEffect(probeIx);
        float gold = Convert::toFloat(expectedFeat.getData(0, count++).c_str());
        CPPUNIT_ASSERT(Convert::doubleCloseEnough(estimated, gold, 4));
      }
    }
  }
  /* Check target effects. */
  for(int chipIx = 0; chipIx < qMethod.getNumTargets(); chipIx++) {
    float estimated = qMethod.getSignalEstimate(chipIx);
    float gold = Convert::toFloat(expected.getData(0, chipIx).c_str());
    //cout << "Target: " << gold << "\t" << estimated << "\n";
    CPPUNIT_ASSERT(Convert::doubleCloseEnough(estimated, gold, 3));
  }
  // whew, all ok!
  return true;
}

void QuantExprMethodTest::testQuantRma() {
  QuantRma rma;
  CPPUNIT_ASSERT( testQuantExprMethod(rma, "input/quant-rma.txt", "expected/quant-rma-expected.txt", 
                                  "expected/quant-rma-features-expected.txt",
                                  "expected/quant-rma-residuals.txt", true) );
}

void QuantExprMethodTest::testQuantPlier() {
  std::map<std::string,std::string> param;
  QuantPlierParams plierParams;
  SelfDoc doc = QuantPlier::explainSelf();
  QuantPlierBase::fillPlierParams(param, plierParams, doc);
  plierParams.setFixPrecomputed(false);
  plierParams.setNumericalTolerance(0.0);
  plierParams.setSafetyZero(0.0);
  QuantPlier plier(plierParams);
  CPPUNIT_ASSERT( testQuantExprMethod(plier, "input/quant-rma.txt", "expected/quant-plier-expected.txt", 
                                  "expected/quant-plier-features-expected.txt",
                                  "expected/quant-plier-residuals.txt", true) );
}

void QuantExprMethodTest::testQuantIterPlier() {
  std::map<std::string,std::string> param;
  QuantPlierParams plierParams;
  SelfDoc doc = QuantIterPlier::explainSelf();
  QuantPlierBase::fillPlierParams(param, plierParams, doc);
  plierParams.setFixPrecomputed(false);
  plierParams.setNumericalTolerance(0.0);
  plierParams.setSafetyZero(0.0);
  QuantIterPlier iPlier(plierParams);
  vector<int> iterations;
  iterations.push_back(22);
  iterations.push_back(11);
  iPlier.setIterations(iterations);
  CPPUNIT_ASSERT( testQuantExprMethod(iPlier, "input/quant-iterplier.txt",
                                     "expected/quant-iterplier-expected.txt",
                                     "expected/quant-iterplier-features-expected.txt",
                                     "expected/quant-iterplier-residuals.txt", false) );
  
  CPPUNIT_ASSERT( testQuantExprMethod(iPlier, "input/quant-iterplier-rev.txt",
                                     "expected/quant-iterplier-expected.txt",
                                     "expected/quant-iterplier-features-expected.txt",
                                     "expected/quant-iterplier-residuals.txt", false) );

}


#endif 
