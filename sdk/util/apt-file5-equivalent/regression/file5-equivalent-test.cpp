////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

/**
 * @file   file5-equivalent-test.cpp
 * @brief  Program for doing regression tests on file5-equivalent.
 */

/// @todo we should really have a test for negative result

#include "util/RegressionTest.h"
#include "util/FsTestDir.h"

using namespace std;

class File5EquivalentTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  File5EquivalentTest() {
    numPassed = 0;
    numFailed = 0;
  }
  void testDifferencesPositive();
  void testEpsilonNegative();
  void testIntegerNegative();
  void testStringNegative();
  void testNonFiniteNegative1();
  void testNonFiniteNegative2();
  void testSignAllowNegationPositive1();
  void testSignAllowNegationPositive2();
  void testSignAllowNegationNegative1();
  void testSignAllowNegationNegative2();
  void testReportNanNumDiffPositive();
  void testReportNanNumDiffNegative();
  void testIgnoreNanNumDifferences();
  void testParametersNegative();
  void testCorrelationPositive();
  void testCorrelationNegative();
  void testIgnoreListPositive();
};

// allowable differences only in checked datasets and header parameters
void File5EquivalentTest::testDifferencesPositive()
{
  string command = "./apt-file5-equivalent --epsilon 0.005 --correlation 1.1 --datasets-file ./data/datasets_similar.txt ./data/demean_false_diff1.ref.a5 ./data/demean_false_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDifferencesPositive", command, checks);
  Verbose::out (1, "\nDoing testDifferencesPositive()");
  if (differencesTest.run()) {
    ++numPassed;
  }
  else {
    Verbose::out (1, "Error in File5EquivalentTest::testDifferencesPositive(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// datasets SketchCN, AntigenomicProbes, MedianSignals and header parameters have non-allowed differences
void File5EquivalentTest::testEpsilonNegative()
{
  string command = "./apt-file5-equivalent -e 0.000001 -d ./data/datasets_similar.txt ./data/demean_false_diff1.ref.a5 ./data/demean_false_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testEpsilonNegative", command, checks);
  Verbose::out (1, "\nDoing testEpsilonNegative()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testEpsilonNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// dataset AntigenomicProbes.ProbeID int32 values have non-allowed single digit differences: 1037413 vs 1037412, etc.
void File5EquivalentTest::testIntegerNegative()
{
  string command = "./apt-file5-equivalent -e 0.005 -c 0.99 -d ./data/datasets_integer_diff.txt ./data/demean_false_diff1.ref.a5 ./data/demean_false_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testIntegerNegative", command, checks);
  Verbose::out (1, "\nDoing testIntegerNegative()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testIntegerNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// dataset MedianSignals.probeset_id has string value has non-allowed difference: 'C-10XYZ' vs 'C-10T0Q'
void File5EquivalentTest::testStringNegative()
{
  string command = "./apt-file5-equivalent -e 0.005 -c 0.999 -d ./data/datasets_string_diff.txt ./data/demean_false_diff1.ref.a5 ./data/demean_false_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testStringNegative", command, checks);
  Verbose::out (1, "\nDoing testStringNegative()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testStringNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// dataset SketchCN.Sketch has non-allowed difference: NaN vs 1782.38155470
void File5EquivalentTest::testNonFiniteNegative1()
{
  string command = "./apt-file5-equivalent -e 0.005 -c 1.0 -d ./data/datasets_non-finite_diff.txt ./data/demean_false_diff1.ref.a5 ./data/demean_false_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testNonFiniteNegative1", command, checks);
  Verbose::out (1, "\nDoing testNonFiniteNegative1()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testNonFiniteNegative1(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// dataset SketchCN.Sketch has non-allowed difference: Inf vs 1782.38155470
void File5EquivalentTest::testNonFiniteNegative2()
{
  string command = "./apt-file5-equivalent -e 0.005 -c 1.1 -d ./data/datasets_non-finite_diff.txt ./data/demean_false_diff1.ref.a5 ./data/demean_true_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testNonFiniteNegative2", command, checks);
  Verbose::out (1, "\nDoing testNonFiniteNegative2()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testNonFiniteNegative2(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// some parameters are different with all values treated as text strings
void File5EquivalentTest::testParametersNegative()
{
  string command = "./apt-file5-equivalent -e 0 -d ./data/datasets_params_diff.txt ./data/demean_false_diff1.ref.a5 ./data/demean_false_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testCorrelationNegative", command, checks);
  Verbose::out (1, "\nDoing testCorrelationNegative()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testCorrelationNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// datasets WaveCorrection X1,X2 have all same values with changed sign, X3,X4,X5,X6 have only some values with changed sign
void File5EquivalentTest::testSignAllowNegationPositive1()
{
  string command = "./apt-file5-equivalent -s -e 0.0001 -c 0.9999 -d ./data/datasets_correlation.txt ./data/demean_true_diff2.ref.a5 ./data/demean_true_diff3.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testSignAllowNegationPositive1", command, checks);
  Verbose::out (1, "\nDoing testSignAllowNegationPositive1()");
  if (differencesTest.run()) {
    ++numPassed;
  }
  else {
    Verbose::out (1, "Error in File5EquivalentTest::testSignAllowNegationPositive1(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// datasets WaveCorrection X1,X2 have all same values with changed sign, equivalent to epsilon=0 without -s based on correlation with X3,X4,X5,X6 excluded
void File5EquivalentTest::testSignAllowNegationPositive2()
{
  string command = "./apt-file5-equivalent -e 0 -c 0.999999 -d ./data/datasets_correlation.txt -i ./data/ignore_columns_wc-except-x1-x2.txt ./data/demean_true_diff2.ref.a5 ./data/demean_true_diff3.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testSignAllowNegationPositive2", command, checks);
  Verbose::out (1, "\nDoing testSignAllowNegationPositive2()");
  if (differencesTest.run()) {
    ++numPassed;
  }
  else {
    Verbose::out (1, "Error in File5EquivalentTest::testSignAllowNegationPositive2(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// datasets WaveCorrection X3,X4 have 10+ values with changed sign plus small differences larger than epsilon=0.00001
void File5EquivalentTest::testSignAllowNegationNegative1()
{
  string command = "./apt-file5-equivalent --sign-allow-negation -e 0.00001 -c 0.9999 --datasets-file ./data/datasets_correlation.txt ./data/demean_true_diff2.ref.a5 ./data/demean_true_diff3.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testSignAllowNegationNegative1", command, checks);
  Verbose::out (1, "\nDoing testSignAllowNegationNegative1()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testSignAllowNegationNegative1(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// datasets WaveCorrection X1,X2,X3,X4,X5,X6 have values with changed sign; differences fail larger tolerance without -s option
void File5EquivalentTest::testSignAllowNegationNegative2()
{
  string command = "./apt-file5-equivalent -e 0.001 -c 0.999 --datasets-file ./data/datasets_correlation.txt ./data/demean_true_diff2.ref.a5 ./data/demean_true_diff3.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testSignAllowNegationNegative2", command, checks);
  Verbose::out (1, "\nDoing testSignAllowNegationNegative2()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testSignAllowNegationNegative2(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// demean_true_diff3.ref.a5 WaveCorrection datasets X7,X8 have NaN values, verify passes if NaN status is same by comparing file to itself
void File5EquivalentTest::testReportNanNumDiffPositive()
{
  string command = "./apt-file5-equivalent -n true -e 0.00001 -c 0.999999 -d ./data/datasets_correlation.txt -i ./data/ignore_columns_wc-except-x7-x8.txt ./data/demean_true_diff3.ref.a5 ./data/demean_true_diff3.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testReportNanNumDiffPositive", command, checks);
  Verbose::out (1, "\nDoing testReportNanNumDiffPositive()");
  if (differencesTest.run()) {
    ++numPassed;
  }
  else {
    Verbose::out (1, "Error in File5EquivalentTest::testReportNanNumDiffPositive(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// demean_true_diff3.ref.a5 WaveCorrection datasets X7,X8 have NaN values, verify fails if NaN status is different by comparing file without NaN values
void File5EquivalentTest::testReportNanNumDiffNegative()
{
  string command = "./apt-file5-equivalent --report-nan-num-diff true -e 0.00001 -c 0.9999 --datasets-file ./data/datasets_correlation.txt -i ./data/ignore_columns_wc-except-x7-x8.txt ./data/demean_true_diff2.ref.a5 ./data/demean_true_diff3.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testReportNanNumDiffNegative", command, checks);
  Verbose::out (1, "\nDoing testReportNanNumDiffNegative()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testReportNanNumDiffNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// ddemean_true_diff3.ref.a5 WaveCorrection datasets X7,X8 have NaN values, verify that it ignores differences if the -n option is false
void File5EquivalentTest::testIgnoreNanNumDifferences()
{
  string command = "./apt-file5-equivalent -n false -e 0.00001 -c 0.999999 -d ./data/datasets_correlation.txt -i ./data/ignore_columns_wc-except-x7-x8.txt ./data/demean_true_diff2.ref.a5 ./data/demean_true_diff3.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testIgnoreNanNumDifferences", command, checks);
  Verbose::out (1, "\nDoing testIgnoreNanNumDifferences()");
  if (differencesTest.run()) {
    ++numPassed;
  }
  else {
    Verbose::out (1, "Error in File5EquivalentTest::testIgnoreNanNumDifferences(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// datasets WaveCorrection has corr=0.99999999820667 above threshold
void File5EquivalentTest::testCorrelationPositive()
{
  string command = "./apt-file5-equivalent -e 0.00000001 -c 0.999 --datasets-file ./data/datasets_correlation.txt ./data/demean_false_diff1.ref.a5 ./data/demean_false_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testCorrelationPositive", command, checks);
  Verbose::out (1, "\nDoing testCorrelationPositive()");
  if (differencesTest.run()) {
    ++numPassed;
  }
  else {
    Verbose::out (1, "Error in File5EquivalentTest::testCorrelationPositive(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// datasets WaveCorrection has corr=0.1946 below threshold
void File5EquivalentTest::testCorrelationNegative()
{
  string command = "./apt-file5-equivalent -e 0.00000001 -c 0.90 --datasets-file ./data/datasets_correlation.txt ./data/demean_false_diff2.ref.a5 ./data/demean_true_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testCorrelationNegative", command, checks);
  Verbose::out (1, "\nDoing testCorrelationNegative()");
  if (!differencesTest.run()) {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in File5EquivalentTest::testCorrelationNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// ignore column with integer differences Cyto2.AntigenomicProbes.ProbeID
void File5EquivalentTest::testIgnoreListPositive()
{
  string command = "./apt-file5-equivalent -e 0.005 -c 1.1 -d ./data/datasets_integer_diff.txt -i ./data/ignore_columns.txt ./data/demean_false_diff1.ref.a5 ./data/demean_false_diff2.ref.a5";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testIgnoreListPositive", command, checks);
  Verbose::out (1, "\nDoing testIgnoreListPositive()");
  if (differencesTest.run()) {
    ++numPassed;
  }
  else {
    Verbose::out (1, "Error in File5EquivalentTest::testIgnoreListPositive(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

int main (int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("util/file5-equivalent", true);
    
    File5EquivalentTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);
    test.testDifferencesPositive();
    test.testEpsilonNegative();
  //test.testIntegerNegative();       //rsatin TODO fix APT-312
    test.testStringNegative();
  //test.testNonFiniteNegative1();    //rsatin TODO fix APT-312
  //test.testNonFiniteNegative2();    //rsatin TODO fix APT-312
    test.testParametersNegative();
    test.testSignAllowNegationPositive1();
    test.testSignAllowNegationPositive2();
    test.testSignAllowNegationNegative1();
    test.testSignAllowNegationNegative2();
    test.testReportNanNumDiffPositive();
    test.testReportNanNumDiffNegative();
	test.testIgnoreNanNumDifferences();
    test.testCorrelationPositive();
    test.testCorrelationNegative();
  //test.testIgnoreListPositive();    //rsatin TODO fix APT-312
   Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
