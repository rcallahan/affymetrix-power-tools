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

class SummaryGenotypeTest {

public:
  int numPassed, numFailed;
  std::string testDir;
  SummaryGenotypeTest() {
    numPassed = 0;
    numFailed = 0;
  }
  void doBrlmmpSnp5Two();
  void doBrlmmpSnp5TwoA5();
  void doBrlmmpSnp5();

};

void SummaryGenotypeTest::doBrlmmpSnp5Two() {
  string command = "./apt-summary-genotype "
    "-p CM=1.bins=100.K=2.SB=0.003.MS=0.05 "
    "-o " + testDir + "/qt-doBrlmmpSnp5Two "
    "--summaries-a=../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/alleleA.plier.summary.txt "
    "--summaries-b=../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/alleleB.plier.summary.txt "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-genders ../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/brlmm-p.genders.txt";

  string name = "qt-doBrlmmpSnp5Two";
  vector<RegressionCheck *> checks;

  checks.push_back(new MatrixCheck(testDir + "/qt-doBrlmmpSnp5Two/brlmm-p.calls.txt",
                                   "../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/qt-doBrlmmpSnp5Two/brlmm-p.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));

  checks.push_back(new MatrixCheck(testDir + "/qt-doBrlmmpSnp5Two/brlmm-p.confidences.txt",
                                   "../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/qt-doBrlmmpSnp5Two/brlmm-p.confidences.txt", 
                                   0.0001,
                                   1, 1, false, 0));

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SummaryGenotypeTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void SummaryGenotypeTest::doBrlmmpSnp5TwoA5() {
  string command = "./apt-summary-genotype "
    "-p CM=1.bins=100.K=2.SB=0.003.MS=0.05 "
    "-o " + testDir + "/qt-doBrlmmpSnp5TwoA5 "
    "--file5-summaries=../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/summaries.a5 "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-genders ../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/brlmm-p.genders.txt";

  string name = "qt-doBrlmmpSnp5TwoA5";
  vector<RegressionCheck *> checks;

  checks.push_back(new MatrixCheck(testDir + "/qt-doBrlmmpSnp5TwoA5/brlmm-p.calls.txt",
                                   "../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/qt-doBrlmmpSnp5Two/brlmm-p.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));

  checks.push_back(new MatrixCheck(testDir + "/qt-doBrlmmpSnp5Two/brlmm-p.confidences.txt",
                                   "../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/qt-doBrlmmpSnp5Two/brlmm-p.confidences.txt", 
                                   0.0001,
                                   1, 1, false, 0));

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SummaryGenotypeTest::" + name + "(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void SummaryGenotypeTest::doBrlmmpSnp5() {
  string command = "./apt-summary-genotype "
    "-p CM=1.bins=100.K=2.SB=0.003.MS=0.05 "
    "-o " + testDir + "/qt-doBrlmmpSnp5 "
    "-s ../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/brlmm-p.plier.summary.txt "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-genders ../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/brlmm-p.genders.txt";

  string name = "qt-doBrlmmpSnp5";
  vector<RegressionCheck *> checks;

  checks.push_back(new MatrixCheck(testDir + "/qt-doBrlmmpSnp5/brlmm-p.calls.txt",
                                   "../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/qt-doBrlmmpSnp5/brlmm-p.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));

  checks.push_back(new MatrixCheck(testDir + "/qt-doBrlmmpSnp5/brlmm-p.confidences.txt",
                                   "../../../regression-data/data/idata/sum-geno/GenomeWideSNP_5/qt-doBrlmmpSnp5/brlmm-p.confidences.txt", 
                                   0.0001,
                                   1, 1, false, 0));

  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass()) {
    Verbose::out(1, "Error in SummaryGenotypeTest::" + name + "(): " + test.getErrorMsg());
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
    testDir.setTestDir("chipstream/summary-genotype-qt", true);
    
    SummaryGenotypeTest test;
    test.testDir = testDir.asString();

    Verbose::setLevel(2);

    test.doBrlmmpSnp5Two();
    test.doBrlmmpSnp5TwoA5();
    test.doBrlmmpSnp5();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

