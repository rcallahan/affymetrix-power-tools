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
 * @file   calvin-equivalent-test.cpp
 * @brief  Program for doing regression tests on calvin-equivalent.
 */

/// @todo we should really have a test for negative result

#include "util/RegressionTest.h"

using namespace std;

class CalvinEquivalent
{

public:
  int numPassed, numFailed;
  CalvinEquivalent()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void testDifferences();
  void testDifferencesNegative();
  void testDifferencesCorrelation();
  void testDifferencesFractional();
  void testDifferencesFractionalNegative();
  void testDiffHeaderIgnore();
  void testDiffHeaderIgnoreNegative();
  void testDiffColumnsIgnore();
  void testDiffColumnsIgnoreNegative();
  void testDiffDatasetIgnore();
  void testDiffMapEpsilon();
  void testDiffMapEpsilonNegative();
  void testNaNNumDiffPositive();
  void testNaNNumDiffNegative();
  void testIgnoreNaNNumDifferences();
};

void CalvinEquivalent::testDifferences()
{
  string command = "./apt-calvin-equivalent -e 0.005 -c 1.0 -k 0 ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDifferences", command, checks);
  Verbose::out (1, "\nDoing testDifferences()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testDifferences(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void CalvinEquivalent::testDifferencesNegative()
{
  string command = "./apt-calvin-equivalent -e 0.0001 -c 1.0 -k 0 ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDifferencesNegative", command, checks);
  Verbose::out (1, "\nDoing testDifferencesNegative()");
  if (!differencesTest.run())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in CalvinEquivalent::testDifferencesNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

void CalvinEquivalent::testDifferencesCorrelation()
{
  string command = "./apt-calvin-equivalent --epsilon 0.0001 --correlation 0.98 --check-header false ./data/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDifferencesCorrelation", command, checks);
  Verbose::out (1, "\nDoing testDifferencesCorrelation()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testDifferencesCorrelation(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void CalvinEquivalent::testDifferencesFractional()
{
  string command = "./apt-calvin-equivalent --epsilon 0.001 --fraction 0.001 -k 0 ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDifferencesFractional", command, checks);
  Verbose::out (1, "\nDoing testDifferencesFractional()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testDifferencesFractional(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void CalvinEquivalent::testDifferencesFractionalNegative()
{
  string command = "./apt-calvin-equivalent --epsilon 0.001 --fraction 0.00005 --check-header false ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDifferencesFractionalNegative", command, checks);
  Verbose::out (1, "\nDoing testDifferencesFractionalNegative()");
  if (!differencesTest.run())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in CalvinEquivalent::testDifferencesFractionalNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

void CalvinEquivalent::testDiffHeaderIgnore()
{
  string command = "./apt-calvin-equivalent -k 1 -e 0.005 -i ./data/setIgnore1.txt ./data/NA06985_AxiomGWASHuSNP1_HeaderDiff1.chp ./data/NA06985_AxiomGWASHuSNP1_HeaderDiff2.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDiffHeaderIgnore", command, checks);
  Verbose::out (1, "\nDoing testDiffHeaderIgnore()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testDiffHeaderIgnore(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void CalvinEquivalent::testDiffHeaderIgnoreNegative()
{
  string command = "./apt-calvin-equivalent -k true -e 0.005 --ignore-params-file ./data/setIgnore2.txt ./data/NA06985_AxiomGWASHuSNP1_HeaderDiff1.chp ./data/NA06985_AxiomGWASHuSNP1_HeaderDiff2.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDiffHeaderIgnoreNegative", command, checks);
  Verbose::out (1, "\nDoing testDiffHeaderIgnoreNegative()");
  if (!differencesTest.run())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in CalvinEquivalent::testDiffHeaderIgnoreNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

void CalvinEquivalent::testDiffColumnsIgnore()
{
  string command = "./apt-calvin-equivalent -k 1 -e 0.0002 -c 1.0 -f 0.0 -i ./data/setIgnore3.txt ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDiffColumnsIgnore", command, checks);
  Verbose::out (1, "\nDoing testDiffColumnsIgnore()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testDiffColumnsIgnore(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void CalvinEquivalent::testDiffColumnsIgnoreNegative()
{
  string command = "./apt-calvin-equivalent -k 1 -e 0.0001 -c 1.0 -f 0.0 -i ./data/setIgnore3.txt ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDiffColumnsIgnoreNegative", command, checks);
  Verbose::out (1, "\nDoing testDiffColumnsIgnoreNegative()");
  if (!differencesTest.run())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in CalvinEquivalent::testDiffColumnsIgnoreNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

void CalvinEquivalent::testDiffDatasetIgnore()
{
  string command = "./apt-calvin-equivalent -k 1 -e 0.0 -c 1.0 -f 0.0 -i ./data/setIgnore4.txt ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDiffDatasetIgnore", command, checks);
  Verbose::out (1, "\nDoing testDiffDatasetIgnore()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testDiffDatasetIgnore(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void CalvinEquivalent::testDiffMapEpsilon()
{
  string command = "./apt-calvin-equivalent -k 1 -e 0 -m ./data/mapEpsilon1.txt ./data/NA06985_GW6_C.birdseed-v2_Gold_small.chp ./data/NA06985_GW6_C.birdseed-v2_bad8c_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDiffMapEpsilon", command, checks);
  Verbose::out (1, "\nDoing testDiffMapEpsilon()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testDiffMapEpsilon(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void CalvinEquivalent::testDiffMapEpsilonNegative()
{
  string command = "./apt-calvin-equivalent --check-header true --epsilon 0 --epsilon-map-file ./data/mapEpsilon2.txt ./data/NA06985_GW6_C.birdseed-v2_Gold_small.chp ./data/NA06985_GW6_C.birdseed-v2_bad8c_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testDiffMapEpsilonNegative", command, checks);
  Verbose::out (1, "\nDoing testDiffMapEpsilonNegative()");
  if (!differencesTest.run())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in CalvinEquivalent::testDiffMapEpsilonNegative(): Negative Test SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

// file NA06985_GW6_C.birdseed-v2_ValueDiff1_NaN_small.chp has NaN values, verify passes if NaN status is same by comparing file to itself
void CalvinEquivalent::testNaNNumDiffPositive()
{
  string command = "./apt-calvin-equivalent -n true -e 0.005 -c 1.0 -k 0 ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_NaN_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_NaN_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testNaNNumDiffPositive", command, checks);
  Verbose::out (1, "\nDoing testNaNNumDiffPositive()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testNaNNumDiffPositive(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

// file NA06985_GW6_C.birdseed-v2_ValueDiff1_NaN_small.chp has NaN values, verify fails if NaN status is different by comparing file without NaN values
void CalvinEquivalent::testNaNNumDiffNegative()
{
  string command = "./apt-calvin-equivalent --report-nan-num-diff true -e 0.005 -c 1.0 -k 0 ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_NaN_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testNaNNumDiffNegative", command, checks);
  Verbose::out (1, "\nDoing testNaNNumDiffNegative()");
  if (!differencesTest.run())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in CalvinEquivalent::testNaNNumDiffNegative(): Negative Test SHOULD HAVE FAILED!!");
    ++numFailed;
  }
}

// file NA06985_GW6_C.birdseed-v2_ValueDiff1_NaN_small.chp has NaN values, verify that NaN differences are ignored if -n option is false
void CalvinEquivalent::testIgnoreNaNNumDifferences()
{
  string command = "./apt-calvin-equivalent -n false -e 0.005 -c 1.0 -k 0 ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_NaN_small.chp ./data/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testIgnoreNaNNumDifferences", command, checks);
  Verbose::out (1, "\nDoing testIgnoreNaNNumDifferences()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in CalvinEquivalent::testIgnoreNaNNumDifferences(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

int main (int argc, char* argv[])
{
  try {
    CalvinEquivalent test;
    Verbose::setLevel(2);
    test.testDifferences();
    test.testDifferencesNegative();
    test.testDifferencesCorrelation();
    test.testDifferencesFractional();
    test.testDifferencesFractionalNegative();
    test.testDiffHeaderIgnore();
    test.testDiffHeaderIgnoreNegative();
    test.testDiffColumnsIgnore();
    test.testDiffColumnsIgnoreNegative();
    test.testDiffDatasetIgnore();
    test.testDiffMapEpsilon();
    test.testDiffMapEpsilonNegative();
	test.testNaNNumDiffPositive();
	test.testNaNNumDiffNegative();
	test.testIgnoreNaNNumDifferences();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
