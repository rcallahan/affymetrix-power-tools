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
 * @file   probeset-genotype-test.cpp
 * @author Chuck Sugnet
 * @date   Thu Mar  2 16:39:16 2006
 * 
 * @brief  Program for doing regression tests on probeset-genotype data.
 * 
 * 
 */

#include "util/Convert.h"
#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/RegressionTest.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>
//

using namespace std;

class SnpSummaryTest {

public:
  int numPassed, numFailed;
  std::string testDir;
  SnpSummaryTest() {
    numPassed = 0;
    numFailed = 0;
  }
  void doCalvinChp();
  void doXdaChp();
  void doText();

};

void SnpSummaryTest::doCalvinChp() {
  string name = "qt-doCalvinChp";
  string command = "./apt-snp-summary -s " + testDir + "/qt-doCalvinChp/table.txt --chp-files=../../../regression-data/data/idata/snp-summary/doBrlmmpSnp5Chp.chpfiles.txt";
  vector<RegressionCheck *> checks;

  checks.push_back(new MixedFileCheck(testDir + "/qt-doCalvinChp/table.txt",
                                   "../../../regression-data/data/idata/snp-summary/qt-doCalvinChp/table.txt", 
                                   0.0000001,
                                   0, 0));

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SnpSummaryTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void SnpSummaryTest::doXdaChp() {
// engine does not support xda chp files
/*
  string name = "qt-doXdaChp";
  string command = "./apt-snp-summary -s " + testDir + "/qt-doXdaChp/table.txt --chp-files=../../../regression-data/data/idata/snp-summary/doChpFiles.chpfiles.txt";
  vector<RegressionCheck *> checks;

  checks.push_back(new MixedFileCheck(testDir + "/qt-doXdaChp/table.txt",
                                   "../../../regression-data/data/idata/snp-summary/qt-doXdaChp/table.txt", 
                                   0.0000001,
                                   0, 0));

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SnpSummaryTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
*/
}

void SnpSummaryTest::doText() {
  string name = "qt-doText";
  string command = "./apt-snp-summary -s " + testDir + "/qt-doText/table.txt --call-file=../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.calls.txt";
  vector<RegressionCheck *> checks;

  checks.push_back(new MixedFileCheck(testDir + "/qt-doText/table.txt",
                                   "../../../regression-data/data/idata/snp-summary/qt-doText/table.txt", 
                                   0.0000001,
                                   0, 0));

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SnpSummaryTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

/** Everybody's favorite function. */
int main(int argc, char* argv[]) {
  try {
    FsTestDir testDir;
    testDir.setTestDir("chipstream/snp-summary-qt", true);
    
    SnpSummaryTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);

    test.doXdaChp();
    test.doCalvinChp();
    test.doText();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

