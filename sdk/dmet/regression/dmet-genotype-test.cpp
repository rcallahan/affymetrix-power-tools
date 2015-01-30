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
 * @file   dmet-genotype-test.cpp
 * @author Alan Williams
 * @date   Fri Sep 12 12:07:04 PDT 2008
 * 
 * @brief  Program for doing regression tests on dmet-genotype data.
 */

#include "util/CalvinChpCheck.h"
#include "util/ChpCheck.h"
#include "util/Convert.h"
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

string m_celsSuffix = ".CEL";
const char* m_cels[]={
  "ZQ_C01_NA19193",    "ZQ_C02_NA12004",    "ZQ_C03_NA19194",
  "ZQ_C04_NA18856",    "ZQ_C05_NA19137",    "ZQ_C06_NA19144",
  "ZQ_C07_NA12003",    "ZQ_C08_NA19143",    "ZQ_C09_NA10838",
  "ZQ_C10_NA19130",    "ZQ_C11_NA18855",    "ZQ_C12_NA19139",
  "ZQ_D01_NA18856",    "ZQ_D02_NA19194",    "ZQ_D03_NA19138",
  "ZQ_D04_NA19131",    "ZQ_D05_NA10838",    "ZQ_D06_NA19192",
  "ZQ_D07_NA19145",    "ZQ_D08_NA19132",    "ZQ_D09_NA18857",
  "ZQ_D10_gCtrl_1",    "ZQ_D11_gCtrl_2",    "ZQ_D12_gCtrl_3",
  NULL
};



class DmetGenotypeTest {

public:
  int numPassed, numFailed;
  DmetGenotypeTest() {
    numPassed = 0;
    numFailed = 0;
  }
  void doRefBuild( const FsTestDir & );
  void doRefRun( const FsTestDir & );
  void doPlasmidRun( const FsTestDir & );
  string createString(const char **items, int size, const char *prefix, const char *suffix);
  // vector<string> m_cels;
  
};

string DmetGenotypeTest::createString(const char **items, int size, const char *prefix, const char *suffix) {
    string val = "";
    for(int i = 0; i < size; i++) {
        val = val + prefix + items[i] + suffix;
    }
    return val;
}

void DmetGenotypeTest::doRefBuild( const FsTestDir & testDir ) {
  
  std::string doRefBuildPath = Fs::join(testDir.asString(), "doRefBuild");
  
  string command = "./apt-dmet-genotype "
    "--batch-info=true "
    "--null-context=false "
    "--reference-input=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.genomic.ref.a5 "
    "--out-dir="
    + doRefBuildPath +
    " "
    "--cel-files=../../regression-data/data/idata/dmet-genotype/DMET_Plus/cels.small.txt "
    "--cdf-file=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cdf "
    "--chrX-probes=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.chrXprobes "
    "--chrY-probes=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.chrYprobes "
    "--batch-name=dynamic-batch-1 "
    "--reference-output="
    + Fs::join(doRefBuildPath, "reference.a5 ") +
    "--region-model=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cn-region-models.txt "
    "--probeset-model=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cn-probeset-models.txt "
    "--probeset-ids=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.genomic.gt.ps "
    "--probeset-ids-reported=../../regression-data/data/idata/dmet-genotype/DMET_Plus/consent1.txt "
    "--cn-region-gt-probeset-file=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cn-gt.ps "
    "--cc-chp-output=true "
    "-- "
    "--summaries=true "
    "--feat-effects "
    "-- "
    "--text-output "
    "-- "
    "--table-output "
    "--feat-effects "
    "--summaries ";

  string name = "doRefBuild";
  vector<RegressionCheck *> checks;

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(m_cels, "../../regression-data/data/idata/dmet-genotype/DMET_Plus/doRefBuild/chp/", ".dmet.chp");
  gen = Util::addPrefixSuffix(m_cels, Fs::join(doRefBuildPath, "chp") + "/", ".dmet.chp");
  
  checks.push_back(new CalvinChpCheck(gen, gold));

  // Run the test
  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in DmetGenotypeTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void DmetGenotypeTest::doRefRun( const FsTestDir & testDir ) {
  std::string doRefRunPath = Fs::join(testDir.asString(), "doRefRun");

  string command = "./apt-dmet-genotype "
    "--batch-info=false "
    "--null-context=true "
    "--reference-input=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.genomic.ref.a5 "
    "--out-dir="
    + doRefRunPath +
    " "
    "--cel-files=../../regression-data/data/idata/dmet-genotype/DMET_Plus/cels.small.txt "
    "--spf-file=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.spf "
    "--chrX-probes=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.chrXprobes "
    "--chrY-probes=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.chrYprobes "
    "--region-model=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cn-region-models.txt "
    "--probeset-model=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cn-probeset-models.txt "
    "--probeset-ids=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.genomic.gt.ps "
    "--probeset-ids-reported=../../regression-data/data/idata/dmet-genotype/DMET_Plus/consent2.txt "
    "--cn-region-gt-probeset-file=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cn-gt.ps "
    "--cc-chp-output=true "
    "--geno-call-thresh=0.01 "
    "-- "
    "--summaries=true "
    "--feat-effects "
    "-- "
    "--text-output "
    "-- "
    "--table-output "
    "--feat-effects "
    "--summaries";

  string name = "doRefRun";
  vector<RegressionCheck *> checks;

  //Chp File Checks

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(m_cels, "../../regression-data/data/idata/dmet-genotype/DMET_Plus/doRefRun/chp/", ".dmet.chp");
  gen = Util::addPrefixSuffix(m_cels, doRefRunPath + "/chp/", ".dmet.chp");

  checks.push_back(new CalvinChpCheck(gen, gold));

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in DmetGenotypeTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void DmetGenotypeTest::doPlasmidRun( const FsTestDir & testDir ) {

  std::string doPlasmidRunPath = Fs::join(testDir.asString(), "doPlasmidRun");

  string command = "./apt-dmet-genotype "
    "--reference-input=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.plasmid.ref.a5 "
    "--out-dir="
    + doPlasmidRunPath +
    " "
    "--cel-files=../../regression-data/data/idata/dmet-genotype/DMET_Plus/cels.small.txt "
    "--cdf-file=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cdf "
    "--chrX-probes=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.chrXprobes "
    "--chrY-probes=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.chrYprobes "
    "--region-model=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.cn-region-models.txt "
    "--probeset-ids=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.plasmid.gt.ps "
    "--probeset-ids-reported=../../regression-data/data/idata/lib/DMET_Plus/DMET_Plus.r3.plasmid.gt.ps "
    "--cc-chp-output=true "
    "--geno-call-thresh=0.2 "
    "--run-cn-engine=false";

  string name = "doPlasmidRun";
  vector<RegressionCheck *> checks;

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(m_cels, "../../regression-data/data/idata/dmet-genotype/DMET_Plus/doPlasmidRun/chp/", ".dmet.chp");
  gen = Util::addPrefixSuffix(m_cels, doPlasmidRunPath + "/chp/", ".dmet.chp");

  checks.push_back(new CalvinChpCheck(gen, gold));

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in DmetGenotypeTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

/** Everybody's favorite function. */
int main(int argc, char* argv[]) {
  try {
    DmetGenotypeTest test;
    Verbose::setLevel(2);
    FsTestDir testDir;
    testDir.setTestDir("dmet", true);
    test.doRefBuild( testDir );
    test.doRefRun( testDir );
    test.doPlasmidRun( testDir );

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

