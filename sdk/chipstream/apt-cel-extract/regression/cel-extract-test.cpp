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

string tissueCels = "  ../../../regression-data/data/cel/HG-U133_Plus_2/pancrease-rep1.cel ../../../regression-data/data/cel/HG-U133_Plus_2/pancrease-rep2.cel ../../../regression-data/data/cel/HG-U133_Plus_2/pancrease-rep3.cel";

string axiomCels = " ../../../regression-data/data/cel/Axiom_GW_Hu_SNP/NA18526.CEL ../../../regression-data/data/cel/Axiom_GW_Hu_SNP/NA18856.CEL ../../../regression-data/data/cel/Axiom_GW_Hu_SNP/NA18862.CEL";

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
  void doExtractMultiChannelIntensitiesSubset();
  void doExtractMultiChannelIntensitiesRaw();

private:
  string errorMsg;
};

void CelExtractTest::doExtractIntensitiesRaw()
{
  string command = "./apt-cel-extract -o " + testDir + "/extract.intensities.raw.txt " + tissueCels;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(testDir + "/extract.intensities.raw.txt",
                                      "../../../regression-data/data/chipstream/cel-extract/HG-U133_Plus_2/extract.intensities.raw.txt",
                                      0.1,
                                      0,0));
  RegressionTest extractTest ("doExtractIntensitiesRaw",command.c_str(), checks);
  Verbose::out (1, "Doing doExtractIntensitiesRaw()");
  if (!extractTest.pass()) {
    Verbose::out(1, "Error in doExtractIntensitiesRaw():" + extractTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void CelExtractTest::doExtractMultiChannelIntensitiesRaw()
{
  string command = "./apt-cel-extract -o " + testDir + "/extract.mc.intensities.raw.txt " + axiomCels;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(testDir + "/extract.mc.intensities.raw.txt.channel0",
                                      "../../../regression-data/data/chipstream/cel-extract/Axiom_GW_Hu_SNP/extract.mc.intensities.raw.txt.channel0",
                                      0.1,
                                      0,0));
  checks.push_back(new MixedFileCheck(testDir + "/extract.mc.intensities.raw.txt.channel1",
                                      "../../../regression-data/data/chipstream/cel-extract/Axiom_GW_Hu_SNP/extract.mc.intensities.raw.txt.channel1",
                                      0.1,
                                      0,0));
  RegressionTest extractTest ("doExtractMultiChannelIntensitiesRaw",command.c_str(), checks);
  Verbose::out (1, "Doing doExtractMultiChannelIntensitiesRaw()");
  if (!extractTest.pass()) {
    Verbose::out(1, "Error in doExtractMultiChannelIntensitiesRaw():" + extractTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void CelExtractTest::doExtractIntensities()
{
  string command = "./apt-cel-extract -o " + testDir + "/extract.intensities.txt -d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf " + tissueCels;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(testDir + "/extract.intensities.txt",
                                      "../../../regression-data/data/chipstream/cel-extract/HG-U133_Plus_2/extract.intensities.txt",
                                      0.1,
                                      0,0));
  RegressionTest extractTest ("doExtractIntensities", command.c_str(), checks);
  Verbose::out (1, "Doing doExtractIntensities()");
  if (!extractTest.pass()) {
    Verbose::out(1, "Error in doExtractIntensities():" + extractTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}
void CelExtractTest::doExtractIntensitiesSubset()
{
  string command = "./apt-cel-extract --probeset-ids=../../../regression-data/data/chipstream/cel-extract/HG-U133_Plus_2/ps.txt -o " + testDir + "/extract.intensities.sub.txt -d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf " + tissueCels;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(testDir + "/extract.intensities.sub.txt",
                                      "../../../regression-data/data/chipstream/cel-extract/HG-U133_Plus_2/extract.intensities.sub.txt",
                                      0.1,
                                      0,0));
  RegressionTest extractTest ("doExtractIntensitiesSubset", command.c_str(), checks);
  Verbose::out (1, "Doing doExtractIntensitiesSubset()");
  if (!extractTest.pass()) {
    Verbose::out(1, "Error in doExtractIntensitiesSubset():" + extractTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void CelExtractTest::doExtractMultiChannelIntensitiesSubset()
{
  string command = "./apt-cel-extract --probeset-ids=../../../regression-data/data/chipstream/cel-extract/Axiom_GW_Hu_SNP/ps.txt -o " + testDir + "/extract.mc.intensities.sub.txt -d ../../../regression-data/data/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.s3.cdf" + axiomCels;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(testDir + "/extract.mc.intensities.sub.txt.channel0",
                                      "../../../regression-data/data/chipstream/cel-extract/Axiom_GW_Hu_SNP/extract.mc.intensities.sub.txt.channel0",
                                      0.1,
                                      0,0));
  checks.push_back(new MixedFileCheck(testDir + "/extract.mc.intensities.sub.txt.channel1",
                                      "../../../regression-data/data/chipstream/cel-extract/Axiom_GW_Hu_SNP/extract.mc.intensities.sub.txt.channel1",
                                      0.1,
                                      0,0));
  RegressionTest extractTest ("doExtractMultiChannelIntensitiesSubset", command.c_str(), checks);
  Verbose::out (1, "Doing doExtractMultiChannelIntensitiesSubset()");
  if (!extractTest.pass()) {
    Verbose::out(1, "Error in doExtractMultiChannelIntensitiesSubset():" + extractTest.getErrorMsg());
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
    testDir.setTestDir("chipstream/cel-extract", true);
    
    CelExtractTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);
    test.doExtractIntensities();
    test.doExtractIntensitiesRaw();
    test.doExtractIntensitiesSubset();
    test.doExtractMultiChannelIntensitiesRaw();
    test.doExtractMultiChannelIntensitiesSubset();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
