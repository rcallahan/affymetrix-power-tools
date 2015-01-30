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
 * @file   canary-test.cpp
 * @author Chuck Sugnet
 * @date   Thu Mar  2 16:39:16 2006
 * 
 * @brief  Program for doing regression tests on canary data.
 * 
 * 
 */

#include "util/CalvinChpCheck.h"
#include "util/ChpCheck.h"
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

class CanaryTest {

public:
  int numPassed, numFailed;
  std::string testDir;
  
  CanaryTest() {
    numPassed = 0;
    numFailed = 0;
  }
  void doChpFilesDefault();
};

void CanaryTest::doChpFilesDefault() {
  string command = "./apt-canary "
    "--out-dir " + testDir + "/qt-doChpFilesDefault "
    "--spf-file ../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.spf "
    "--cnv-region-file ../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.canary-v1.region "
    "--cnv-normalization-file ../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.canary-v1.normalization "
    "--cnv-prior-file ../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.canary-v1.prior "
    "--cnv-map-file ../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.canary-v1.bed "
    "--analysis-name canary-v1 "
    "--cel-files ../../regression-data/data/idata/canary/GenomeWideSNP_6/fas_cel_files.txt "
    "--cc-chp-output=true "
    "--table-output=true";

  string chpFilesSuffix = ".canary-v1";
  const char *chpFiles[] = {
    "NA06985_GW6_C", "NA06993_GW6_C", "NA07000_GW6_C", 
    "NA07019_GW6_C", "NA07029_GW6_C",
    NULL
  };

  string name = "qt-doChpFilesDefault";
  vector<RegressionCheck *> checks;

  // checks for chp files.
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "../../regression-data/data/idata/canary/GenomeWideSNP_6/qt-doChpFilesDefault/", chpFilesSuffix + ".cnvchp");
  gen =  Util::addPrefixSuffix(chpFiles, testDir + "/qt-doChpFilesDefault/", chpFilesSuffix + ".cnvchp");
  
  checks.push_back(new CalvinChpCheck(gen,gold));

  // text output checks
  checks.push_back(new MatrixCheck(testDir + "/qt-doChpFilesDefault/canary-v1.allele_average_summary.txt",
                                   "../../regression-data/data/idata/canary/GenomeWideSNP_6/qt-doChpFilesDefault/canary-v1.allele_average_summary.txt", 
                                   0.0001,
                                   1, 1, false, 0));

  checks.push_back(new MatrixCheck(testDir + "/qt-doChpFilesDefault/canary-v1.calls.txt",
                                   "../../regression-data/data/idata/canary/GenomeWideSNP_6/qt-doChpFilesDefault/canary-v1.calls.txt", 
                                   0.0001,
                                   1, 1, false, 0));

  checks.push_back(new MatrixCheck(testDir + "/qt-doChpFilesDefault/canary-v1.cnv_region_summary.txt",
                                   "../../regression-data/data/idata/canary/GenomeWideSNP_6/qt-doChpFilesDefault/canary-v1.cnv_region_summary.txt", 
                                   0.0001,
                                   1, 1, false, 0));

  checks.push_back(new MatrixCheck(testDir + "/qt-doChpFilesDefault/canary-v1.confidences.txt",
                                   "../../regression-data/data/idata/canary/GenomeWideSNP_6/qt-doChpFilesDefault/canary-v1.confidences.txt", 
                                   0.0001,
                                   1, 1, false, 0));

  checks.push_back(new MatrixCheck(testDir + "/qt-doChpFilesDefault/canary-v1.summary.txt",
                                   "../../regression-data/data/idata/canary/GenomeWideSNP_6/qt-doChpFilesDefault/canary-v1.summary.txt", 
                                   0.0001,
                                   1, 1, false, 0));

  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in CanaryTest::" + name + "(): " + brlmmTest.getErrorMsg());
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
    testDir.setTestDir("canary-qt", true);
    
    CanaryTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);

    test.doChpFilesDefault();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

