////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
 * @file   cel-extract-test.cpp
 * @brief  Program for doing regression tests on cel-extract.
 */
#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/MixedFileCheck.h"
#include "util/RegressionTest.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>
//

using namespace std;

///@todo it would be good to add a lot more tests. w/ analysis string specifically.

string tissueCels = "  ../../../regression-data/data/idata/cel/HG-U133_Plus_2/pancrease-rep1.cel ../../../regression-data/data/idata/cel/HG-U133_Plus_2/pancrease-rep2.cel ../../../regression-data/data/idata/cel/HG-U133_Plus_2/pancrease-rep3.cel";

class CelExtractTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  CelExtractTest()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void doExtractIntensities();
  void doExtractIntensitiesSubset();
  void doExtractIntensitiesRaw();

private:
  string errorMsg;
};

void CelExtractTest::doExtractIntensitiesRaw()
{
  string command = "./apt-cel-extract -o " + testDir + "/qt.extract.intensities.raw.txt " + tissueCels;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(testDir + "/qt.extract.intensities.raw.txt",
                                      "../../../regression-data/data/idata/cel-extract/HG-U133_Plus_2/qt.extract.intensities.raw.txt",
                                      0.1,
                                      0,0));
  RegressionTest extractTest ("qt-doExtractIntensitiesRaw",command.c_str(), checks);
  Verbose::out (1, "Doing qt-doExtractIntensitiesRaw()");
  if (!extractTest.pass()) {
    Verbose::out(1, "Error in qt-doExtractIntensitiesRaw():" + extractTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}
void CelExtractTest::doExtractIntensities()
{
  string command = "./apt-cel-extract -o " + testDir + "/qt.extract.intensities.txt -s ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf " + tissueCels;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(testDir + "/qt.extract.intensities.txt",
                                      "../../../regression-data/data/idata/cel-extract/HG-U133_Plus_2/qt.extract.intensities.txt",
                                      0.1,
                                      0,0));
  RegressionTest extractTest ("qt-doExtractIntensities", command.c_str(), checks);
  Verbose::out (1, "Doing qt-doExtractIntensities()");
  if (!extractTest.pass()) {
    Verbose::out(1, "Error in qt-doExtractIntensities():" + extractTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}
void CelExtractTest::doExtractIntensitiesSubset()
{
  string command = "./apt-cel-extract --probeset-ids=../../../regression-data/data/idata/cel-extract/HG-U133_Plus_2/ps.txt -o " + testDir + "/qt.extract.intensities.sub.txt -s ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf " + tissueCels;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(testDir + "/qt.extract.intensities.sub.txt",
                                      "../../../regression-data/data/idata/cel-extract/HG-U133_Plus_2/qt.extract.intensities.sub.txt",
                                      0.1,
                                      0,0));
  RegressionTest extractTest ("qt-doExtractIntensitiesSubset", command.c_str(), checks);
  Verbose::out (1, "Doing qt-doExtractIntensitiesSubset()");
  if (!extractTest.pass()) {
    Verbose::out(1, "Error in qt-doExtractIntensitiesSubset():" + extractTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

int main (int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("chipstream/cel-extract-qt", true);
    
    CelExtractTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);
    test.doExtractIntensities();
    test.doExtractIntensitiesRaw();
    test.doExtractIntensitiesSubset();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 0;
}
