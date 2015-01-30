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
 * @file   matrix-diff-test.cpp
 * @brief  Program for doing regression tests on matrix-diff.
 */

/// @todo we should really have a test for negative result

#include "util/RegressionTest.h"

using namespace std;

class MatrixDiffTest
{

public:
  int numPassed, numFailed;
  MatrixDiffTest()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void testTextDifferences1();
  void testTextDifferences2();
  void testNumericDifferences();
  void testNumericDifferencesFractional();
  void testMixedDifferences();
  void testMixedDifferencesFractional();
  void testMixedDiffFracNegative();
};

// TextFileCheck POSITIVE test case: skip header and first (and only) row with differences
void MatrixDiffTest::testTextDifferences1()
{
  string command = "./apt-matrix-diff -t -l 2 data/testMatrixDiff2.txt data/testMatrixDiff3.txt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testTextDifferences1", command.c_str(), checks);
  Verbose::out (1, "\nDoing testTextDifferences1()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MatrixDiffTest::testTextDifferences1(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// TextFileCheck NEGATIVE test case: skip header only and not first row with differences
void MatrixDiffTest::testTextDifferences2()
{
  string command = "./apt-matrix-diff --text --line-skip 1 --print-mismatch false data/testMatrixDiff2.txt data/testMatrixDiff3.txt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testTextDifferences2", command.c_str(), checks);
  //differencesTest.m_NegTest = true;  //works only with RegressionTest::pass() method
  Verbose::out (1, "\nDoing testTextDifferences2()");
  if (!differencesTest.run())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in MatrixDiffTest::testTextDifferences2(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// MatrixDiff POSITIVE test case: skip header only and ignore difference < epsilon
void MatrixDiffTest::testNumericDifferences()
{
  string command = "./apt-matrix-diff --version -l 1 -e 0.1 data/test-numeric-gold.txt data/test-numeric-generated.txt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testNumericDifferences", command.c_str(), checks);
  Verbose::out (1, "\nDoing testNumericDifferences()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MatrixDiffTest::testNumericDifferences(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// MatrixDiff POSITIVE test case: skip header only and ignore difference > epsilon based on fraction difference threshold
void MatrixDiffTest::testNumericDifferencesFractional()
{
  string command = "./apt-matrix-diff --line-skip 1 --allowed-mismatch 2 --epsilon 0.001 --fraction 0.00005 data/testMatrixDiff1.txt data/testMatrixDiff2.txt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testNumericDifferencesFractional", command.c_str(), checks);
  Verbose::out (1, "\nDoing testNumericDifferencesFractional()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MatrixDiffTest::testNumericDifferencesFractional(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// MixedFileCheck POSITIVE test case: skip first two data rows with differences in numeric table data
void MatrixDiffTest::testMixedDifferences()
{
  string command = "./apt-matrix-diff --version --mixed -l 2 -p false data/test-mixed-gold.txt data/test-mixed-generated.txt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testMixedDifferences", command.c_str(), checks);
  Verbose::out (1, "Doing testMixedDifferences()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MatrixDiffTest::testMixedDifferences(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// MixedFileCheck POSITIVE test case: allow 2 of 2 large numeric differences per epsilon and fraction criteria
void MatrixDiffTest::testMixedDifferencesFractional()
{
  string command = "./apt-matrix-diff --mixed -l 0 -a 2 -e 0.001 -f 0.00005 --print-mismatch false --print-mismatch-max 1 data/testMixedDiff1.txt data/testMixedDiff2.txt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testMixedDifferencesFractional", command.c_str(), checks);
  Verbose::out (1, "Doing testMixedDifferencesFractional()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MatrixDiffTest::testMixedDifferencesFractional(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// MixedFileCheck NEGATIVE test case: allow 1 of 2 large numeric differences per epsilon and fraction criteria
void MatrixDiffTest::testMixedDiffFracNegative()
{
  string command = "./apt-matrix-diff --mixed -l 0 -a 1 -e 0.001 -f 0.00005 --print-mismatch false --print-mismatch-max 1 data/testMixedDiff1.txt data/testMixedDiff2.txt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testMixedDifferencesFractional", command.c_str(), checks);
  Verbose::out (1, "Doing testMixedDiffFracNegative()");
  if (!differencesTest.run())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in MatrixDiffTest::testMixedDiffFracNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

int main (int argc, char* argv[])
{
  try {
    MatrixDiffTest test;
    Verbose::setLevel(2);
    test.testTextDifferences1();
    test.testTextDifferences2();
    test.testNumericDifferences();
    test.testNumericDifferencesFractional();
    test.testMixedDifferences();
    test.testMixedDifferencesFractional();
    test.testMixedDiffFracNegative();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
