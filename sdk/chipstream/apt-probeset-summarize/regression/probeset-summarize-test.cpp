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
 * @file   probeset-summarize-test.cpp
 * @author Chuck Sugnet
 * @date   Mon Dec  5 14:12:14 2005
 *
 * @brief  Program for doing regression tests on probeset-summarize data.
 *
 */
#include "util/CalvinChpCheck.h"
#include "util/ChpCheck.h"
#include "util/Convert.h"
#include "util/FsTestDir.h"
#include "util/RegressionSuite.h"
#include "util/RegressionTest.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>
//

using namespace std;

string celFileList =" ../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/cel-files.txt";
string celFileList2=" ../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/cel-files.txt";

string tissueCelsPrefix = "../../../regression-data/data/cel/HG-U133_Plus_2/";
string tissueCelsSuffix = ".cel";
const char* tissueCels[]={
  "heart-rep1",     "heart-rep2",     "heart-rep3",
  "hela-rep1",      "hela-rep2",      "hela-rep3",
  "pancrease-rep1", "pancrease-rep2", "pancrease-rep3",
  NULL
};

string tissueCcCelsPrefix = "../../../regression-data/data/cel/HG-U133_Plus_2/";
string tissueCcCelsSuffix = ".agcc.cel";
const char* tissueCcCels[]={
  "heart-rep1",     "heart-rep2",     "heart-rep3",
  "hela-rep1",      "hela-rep2",      "hela-rep3",
  "pancrease-rep1", "pancrease-rep2", "pancrease-rep3",
  NULL
};

string nvissaCelsPrefix = " ../../../regression-data/data/cel/HG-U133_Plus_2/";
string nvissaCelsSuffix = ".CEL";
const char* nvissaCels[]={
  "4221B_a01",     "4221B_b01",     "4221B_c01",
  "4221H_a02",     "4221H_b02",     "4221H_c02",
  "4221HL_a04",    "4221HL_b04",    "4221HL_c04",
  "4221P_a03",     "4221P_b03",     "4221P_c03",
  NULL
};

string huexCelsPrefix = " ../../../regression-data/data/cel/HuEx-1_0-st-v2/";
string huexCelsSuffix = ".CEL";
const char* huexCels[]={
  "huex_wta_cb_A",        "huex_wta_cb_B",        "huex_wta_cb_C",
  "huex_wta_heart_A",     "huex_wta_heart_B",     "huex_wta_heart_C",
  "huex_wta_muscle_A",    "huex_wta_muscle_B",    "huex_wta_muscle_C",
  "huex_wta_testes_A",    "huex_wta_testes_B",    "huex_wta_testes_C",
  NULL
};

class ProbeSetSummarizeTest : public RegressionSuite {

public:
  int numPassed, numFailed;
  std::string testDir;
  ProbeSetSummarizeTest() {
    numPassed = 0;
    numFailed = 0;
  }
  void doRmaTissueSketchSuppliedTest();
  void doPlierPrecompTissueTest();
  void doRmaNvissaTest();
  void doRmaNvissaReadWriteFeatureEffectsTest();
  void doRmaNvissaSketchTest();
  void doRmaTissueTest();
  void doRmaTissueSpfTest();
  void doRmaTissueTestCC();
  void doRmaTissueSketchTest();
  void doPlierGcBgTissueTest();
  void doPlierMMTissueSketchTest();
  void doPlierWtaRefSeqTest();
  void doHumanGeneTest();
  void doHumanGeneSpfTest();
  void doHumanGeneKillListTest();
  void doSNP6CopyNumberWorkflow1a();
  void doSNP6CopyNumberWorkflow1b();
  void doSNP6CopyNumberWorkflow2();
  void doPlierMMTissueMedianNormTest();
  void doDabgSubsetTest();
  void doDabgU133Test();
};

void ProbeSetSummarizeTest::doHumanGeneKillListTest() {

  if ( !Fs::dirExists(testDir + "/doHumanGeneKillListTest") ){
    Fs::mkdirPath(testDir + "/doHumanGeneKillListTest", false);
  }

  // run using kill list, pgf, mps file
  string outdir1 = testDir + "/doHumanGeneKillListTest/mask1";
  string command1 = "./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "--kill-list ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list.txt "
    "-m ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list.mps "
    "-o " + outdir1 + " "
    "--cel-files " + celFileList2;

  // compare against pgf w/ probes removed
  string outdir2 = testDir + "/doHumanGeneKillListTest/mask2";
  string command2 = "./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list.pgf "
    "-m ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list2.mps "
    "-o " + outdir2 + " "
    "--cel-files " + celFileList2;

  // and against CDF masking
  string outdir3 = testDir + "/doHumanGeneKillListTest/mask3";
  string command3 = "./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.cdf "
    "--kill-list ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list.txt "
    "-m ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list.mps "
    "-o " + outdir3 + " "
    "--force "
    "--cel-files " + celFileList2;

  // and when using x/y rather than probe ID
  string outdir4 = testDir + "/doHumanGeneKillListTest/mask4";
  string command4 = "./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "--kill-list ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list.xy.txt "
    "-m ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list.mps "
    "-o " + outdir4 + " " 
    "--cel-files " + celFileList2;

  // no mps file
  string outdir5 = testDir + "/doHumanGeneKillListTest/mask5";
  string command5 = "./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "--kill-list ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list.txt "
    "-o " + outdir5 + " "
    "--cel-files " + celFileList2;

  // this one should fail -- check to see that if masking was not applied, the data would differ
  string outdir6 = testDir + "/doHumanGeneKillListTest/mask6";
  string command6 = "./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "-m ../../../regression-data/data/lib/HuGene-1_0-st-v1/test-kill-list2.mps "
    "-o " + outdir6 + " "
    "--cel-files " + celFileList2;

  RegressionTest test1(
          "doHumanGeneKillListTest-phase1",
          testDir + "/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command1.c_str(), 1, 1, true, 0);
  test1.setSuite(*this, outdir1, outdir1 + "/apt-probeset-summarize.log", outdir1 + "/valgrind.log");
  RegressionTest test2(
          "doHumanGeneKillListTest-phase2",
          testDir + "/doHumanGeneKillListTest/mask2/quant-norm.pm-only.med-polish.summary.txt",
          testDir + "/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command2.c_str(), 1, 1, true, 0);
  test2.setSuite(*this, outdir2, outdir2 + "/apt-probeset-summarize.log", outdir2 + "/valgrind.log");
  RegressionTest test3(
          "doHumanGeneKillListTest-phase3",
          testDir + "/doHumanGeneKillListTest/mask3/quant-norm.pm-only.med-polish.summary.txt",
          testDir + "/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command3.c_str(), 1, 1, true, 0);
  test3.setSuite(*this, outdir3, outdir3 + "/apt-probeset-summarize.log", outdir3 + "/valgrind.log");
  RegressionTest test4(
          "doHumanGeneKillListTest-phase4",
          testDir + "/doHumanGeneKillListTest/mask4/quant-norm.pm-only.med-polish.summary.txt",
          testDir + "/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command4.c_str(), 1, 1, true, 0);
  test4.setSuite(*this, outdir4, outdir4 + "/apt-probeset-summarize.log", outdir4 + "/valgrind.log");
  RegressionTest test5(
          "doHumanGeneKillListTest-phase5",
          testDir + "/doHumanGeneKillListTest/mask5/quant-norm.pm-only.med-polish.summary.txt",
          "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHumanGeneKillListTest/mask5/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command5.c_str(), 1, 1, true, 0);
  test5.setSuite(*this, outdir5, outdir5 + "/apt-probeset-summarize.log", outdir5 + "/valgrind.log");
  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
          testDir + "/doHumanGeneKillListTest/mask6/quant-norm.pm-only.med-polish.summary.txt",
          "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, 1, 1, true, 0));
  checks[0]->m_NegTest = true;
  RegressionTest test6( "doHumanGeneKillListTest-phase6", command6.c_str(), checks);
  test6.setSuite(*this, outdir6, outdir6 + "/apt-probeset-summarize.log", outdir6 + "/valgrind.log");
  Verbose::out(1, "Doing doHumanGeneKillListTest(): Phase 1");
  if(!test1.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneKillListTest(): Phase 1: " + test1.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

  Verbose::out(1, "Doing doHumanGeneKillListTest(): Phase 2");
  if(!test2.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneKillListTest(): Phase 2: " + test2.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

  Verbose::out(1, "Doing doHumanGeneKillListTest(): Phase 3");
  if(!test3.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneKillListTest(): Phase 3: " + test3.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

  Verbose::out(1, "Doing doHumanGeneKillListTest(): Phase 4");
  if(!test4.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneKillListTest(): Phase 4: " + test4.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

  Verbose::out(1, "Doing doHumanGeneKillListTest(): Phase 5");
  if(!test5.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneKillListTest(): Phase 5: " + test5.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

  Verbose::out(1, "Doing doHumanGeneKillListTest(): Phase 6 -- expect this to fail");
  if(!test6.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneKillListTest(): Phase 6: " + test6.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

}

void ProbeSetSummarizeTest::doRmaNvissaTest() {
  string outdir = testDir + "/doRmaNvissaTest";
  string command = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + outdir + " "
    "--use-disk=false "
    "--xda-chp-output "
    "--cc-chp-output ";
  command += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  string chpFilesSuffix = ".rma-bg.quant-norm.pm-only.med-polish.chp";
  // checks for chp files.
  const char *chpFiles[] = {"4221B_a01", "4221B_b01",
                "4221B_c01", "4221HL_a04",
                "4221HL_b04", "4221HL_c04",
                "4221H_a02", "4221H_b02",
                "4221H_c02", "4221P_a03",
                "4221P_b03", "4221P_c03",
                NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  vector<string> goldCC,genCC;
  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaTest/chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, testDir + "/doRmaNvissaTest/chp/", chpFilesSuffix);
  goldCC = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaTest/cc-chp/", chpFilesSuffix);
  genCC = Util::addPrefixSuffix(chpFiles, testDir + "/doRmaNvissaTest/cc-chp/", chpFilesSuffix);

  checks.push_back(new ChpCheck(gen, gold));
  checks.push_back(new CalvinChpCheck(genCC, goldCC));
  checks.push_back(new MatrixCheck(testDir + "/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));

  RegressionTest rmaTest("doRmaNvissaTest", command.c_str(), checks);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doRmaNvissaTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doRmaNvissaTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doRmaNvissaReadWriteFeatureEffectsTest() {

  //First Test
  string outdir = testDir + "/doRmaNvissaReadWriteFeatureEffectsTest";
  string commandMakeFeature = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + outdir + " "
    "--feat-effects ";
    commandMakeFeature += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.0001,
                                   1, 1, false, 0));

  ///@todo need to tighten up
  checks.push_back(new MatrixCheck(testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt",
                                   0.001,
                                   1, 1, false, 0));
  RegressionTest rmaMakeTest("doRmaNvissaReadWriteFeatureEffectsTest-part1", commandMakeFeature.c_str(), checks);
  rmaMakeTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doRmaNvissaReadWriteFeatureEffectsTest(), creating feature effects file.  ");
  bool passed = rmaMakeTest.pass();

  //Second Test
  string outdir2 = testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/useFeatures";
  string commandUseFeature = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + outdir2 + " "
    "--use-feat-eff " + testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt ";
    commandUseFeature += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  vector<RegressionCheck *> useChecks;
  ///@todo need to tighten up once we have binary/a5 input/output
  useChecks.push_back(new MatrixCheck(testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/useFeatures/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.001,
                                   1, 1, false, 0));
  RegressionTest rmaUseTest("doRmaNvissaReadWriteFeatureEffectsTest-part2", commandUseFeature.c_str(), useChecks);
  rmaUseTest.setSuite(*this, outdir2, outdir2 + "/apt-probeset-summarize.log", outdir2 + "/valgrind.log");

  Verbose::out(1, "Doing doRmaNvissaReadWriteFeatureEffectsTest(), using the just computed feature effects file. ");
  passed = passed && rmaUseTest.pass();

  //  Third Test - this case checks that code which was designed to use a new style of feature effects file, continues to function
  //  when the old style of feature effects file in text format is used as input.
  string outdir3 = testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/useText";
  string commandUseOldStyleTextFeature = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + outdir3 + " "
    "--use-feat-eff ../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt ";
    commandUseOldStyleTextFeature += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  vector<RegressionCheck *> useOldStyleTextFeatureChecks;
  useOldStyleTextFeatureChecks.push_back(new MatrixCheck(testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/useText/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.001,
                                   1, 1, false, 0));
  RegressionTest rmaUseOldStyleTextFeatureTest("doRmaNvissaReadWriteFeatureEffectsTest-part3", commandUseOldStyleTextFeature.c_str(), useOldStyleTextFeatureChecks);
  rmaUseOldStyleTextFeatureTest.setSuite(*this, outdir3, outdir3 + "/apt-probeset-summarize.log", outdir3 + "/valgrind.log");
  Verbose::out(1, "Doing doRmaNvissaReadWriteFeatureEffectsTest(), checking use of old style text format feature effects file.  ");
   passed = passed && rmaUseOldStyleTextFeatureTest.pass();

  //  Fourth Test - this case checks that code which was designed to use a new style of feature effects file, continues to function
  //  when the old style of feature effects file in HDF5 format is used as input.
  string outdir4 = testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/useHDF5";
  string commandUseOldStyleHDF5Feature = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + outdir4 + " "
    "--a5-feature-effects-input-file ../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.a5 "
    "--a5-feature-effects-input-name  rma-bg.quant-norm.pm-only.med-polish.feature-response ";
    commandUseOldStyleHDF5Feature += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

#ifndef _WIN32    
  vector<RegressionCheck *> useOldStyleHDF5FeatureChecks;
  useOldStyleHDF5FeatureChecks.push_back(new MatrixCheck(testDir + "/doRmaNvissaReadWriteFeatureEffectsTest/useHDF5/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.001,
                                   1, 1, false, 0));
  RegressionTest rmaUseOldStyleHDF5FeatureTest("doRmaNvissaReadWriteFeatureEffectsTest-part4", commandUseOldStyleHDF5Feature.c_str(), useOldStyleHDF5FeatureChecks);
  rmaUseOldStyleHDF5FeatureTest.setSuite(*this, outdir4, outdir4 + "/apt-probeset-summarize.log", outdir4 + "/valgrind.log");
  passed = passed && rmaUseOldStyleHDF5FeatureTest.pass();
#endif
  
  Verbose::out(1, "Doing doRmaNvissaReadWriteFeatureEffectsTest(), checking use of old style HDF5 format feature effects file.  ");

  if(!passed) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doRmaNvissaReadWriteFeatureEffectsTest(): " + rmaMakeTest.getErrorMsg()
                 + " " + rmaUseTest.getErrorMsg());
    numFailed++;

  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doRmaNvissaSketchTest() {
  string outdir = testDir + "/doRmaNvissaSketchTest";
  string command = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=100000.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + outdir + " ";
    command += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  RegressionTest rmaTest("doRmaNvissaSketchTest", testDir + "/doRmaNvissaSketchTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaNvissaSketchTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 225);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doRmaNvissaSketchTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doRmaNvissaSketchTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doRmaTissueTest() {
  string outdir = testDir + "/doRmaTissueTest";
  string command = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + testDir + "/doRmaTissueTest "
    "--cel-files " + celFileList;

  RegressionTest rmaTest("doRmaTissueTest", testDir + "/doRmaTissueTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaTissueTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doRmaTissueTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doRmaTissueTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doRmaTissueSpfTest() {
  string outdir = testDir + "/doRmaTissueSpfTest";
  string command = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o " + outdir + " "
    "--cel-files " + celFileList;

  RegressionTest rmaTest("doRmaTissueSpfTest", testDir + "/doRmaTissueSpfTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaTissueSpfTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doRmaTissueSpfTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doRmaTissueSpfTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doRmaTissueTestCC() {
  string outdir = testDir + "/doRmaTissueTestCC";
  string command = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + outdir + " ";
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCcCels, tissueCcCelsPrefix, tissueCcCelsSuffix), " ");

  RegressionTest rmaTest("doRmaTissueTestCC", testDir + "/doRmaTissueTestCC/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaTissueTestCC/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doRmaTissueTestCC()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doRmaTissueTestCC(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doPlierPrecompTissueTest() {
  string outdir = testDir + "/doPlierPrecompTissueTest";
  string commandMakeFeat = "./apt-probeset-summarize "
    "-a pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-p ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-b ../../../regression-data/data/lib/HG-U133_Plus_2/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp "
    "-o " + outdir + " "
    "--use-disk=false "
    "--feat-effects "
    "-x 5 ";
    commandMakeFeat += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest plierTest("doPlierPrecompTissueTest", testDir + "/doPlierPrecompTissueTest/pm-gcbg.plier.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doPlierPrecompTissueTest/pm-gcbg.plier.summary.txt",
                         0.0001,
                         commandMakeFeat.c_str(),
                         1, 1, true, 0);
  plierTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");

  Verbose::out(1, "Doing doPlierPrecompTissueTest() phase 1");
  string outdir2 = testDir + "/doPlierPrecompTissueTest/useFeat";
  string commandUseFeat = "./apt-probeset-summarize "
    "-a pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-p ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-b ../../../regression-data/data/lib/HG-U133_Plus_2/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp "
    "-o " + outdir2 + " "
    "--use-disk=false "
    "--use-feat-eff " + testDir + "/doPlierPrecompTissueTest/pm-gcbg.plier.feature-response.txt "
    "-x 5 ";
  commandUseFeat += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ") ;

  RegressionTest plierPrecompTest("doPlierPrecompTissueTest", testDir + "/doPlierPrecompTissueTest/useFeat/pm-gcbg.plier.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doPlierPrecompTissueTest/useFeat/pm-gcbg.plier.summary.txt",
                         0.0001,
                         commandUseFeat.c_str(),
                         1, 1, true, 0);
  plierPrecompTest.setSuite(*this, outdir2, outdir2 + "/apt-probeset-summarize.log", outdir2 + "/valgrind.log");

  bool passed = plierTest.pass();
  Verbose::out(1, "Doing doPlierPrecompTissueTest() phase 2");
  passed = passed && plierPrecompTest.pass();
  if(!passed) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doPlierPrecompTissueTest(): " + plierTest.getErrorMsg()
                 + " " + plierPrecompTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doRmaTissueSketchSuppliedTest() {
  // NOTE: The allowable error on this test has been raised to 0.021
  // to account for solaris sometimes inaccurate "log" function.
  // On ppc and amd64, we get the same answers, but on sparc,
  // the answers are different.  Normally it is in the low bits (8-10th place)
  // but sometimes not.  (Like in the 5th place)
  string outdir = testDir + "/doRmaTissueSketchSuppliedTest";
  string command = "./apt-probeset-summarize "
    "-a quant-norm.sketch=100000.lowprecision=true,pm-only,med-polish.expon=true "
    "--target-sketch ../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/quant-norm.100000.txt "
    "-p ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-o " + outdir + " "
    "-x 6 ";
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest("doRmaTissueSketchSuppliedTest", testDir + "/doRmaTissueSketchSuppliedTest/quant-norm.pm-only.med-polish.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaTissueSketchSuppliedTest/quant-norm.pm-only.med-polish.summary.txt",
                         0.021,
                         command.c_str(),
                         1, 1, true, 0);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");

  Verbose::out(1, "Doing doRmaTissueSketchSuppliedTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTets::doRmaTissueSketchSuppliedTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doRmaTissueSketchTest() {
  string outdir = testDir + "/doRmaTissueSketchTest";
  string command = "./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=100000.bioc=true.usepm=true,pm-only,med-polish "
    "-d ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf "
    "-x 5 "
    "-o " + outdir + " ";
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest("doRmaTissueSketchTest", testDir + "/doRmaTissueSketchTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",

                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doRmaTissueSketchTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.01,
                         command.c_str(),
                         1, 1, true, 1306);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doRmaTissueSketchTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doRmaTissueSketchTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doPlierGcBgTissueTest() {
  string outdir = testDir + "/doPlierGcBgTissueTest";
  string command = "./apt-probeset-summarize "
    "-a pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-b ../../../regression-data/data/lib/HG-U133_Plus_2/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp  "
    "-p ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-x 5 "
    "-o " + outdir + " "
    "--use-disk=false ";
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest("doPlierGcBgTissueTest", testDir + "/doPlierGcBgTissueTest/pm-gcbg.plier.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doPlierGcBgTissueTest/pm-gcbg.plier.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doPlierGcBgTissueTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doPlierGcBgTissueTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doPlierMMTissueSketchTest() {
  string outdir = testDir + "/doPlierMMTissueSketchTest";
  string command = "./apt-probeset-summarize "
    "-a quant-norm.sketch=100000.lowprecision=true,pm-mm,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-p ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-x 6 "
    "-o " + outdir + " ";
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest("doPlierMMTissueSketchTest", testDir + "/doPlierMMTissueSketchTest/quant-norm.pm-mm.plier.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doPlierMMTissueSketchTest/quant-norm.pm-mm.plier.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");

  Verbose::out(1, "Doing doPlierTissueTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doPlierTissueTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doPlierMMTissueMedianNormTest() {
  string outdir = testDir + "/doPlierMMTissueMedianNormTest";
  string command  = "./apt-probeset-summarize "
    "-a med-norm.lowprecision=true,pm-mm,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-p ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-x 5 "
    "-o " + outdir + " ";
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest("doPlierMMTissueMedianNormTest", testDir + "/doPlierMMTissueMedianNormTest/med-norm.pm-mm.plier.summary.txt",
                         "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doPlierMMTissueMedianNormTest/med-norm.pm-mm.plier.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");

  Verbose::out(1, "Doing doPlierMMTissueMedianTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doPlierMMTissueMedianTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doPlierWtaRefSeqTest() {
  string outdir = testDir + "/doPWRST";
  string command = "./apt-probeset-summarize "
    "-a pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-b ../../../regression-data/data/lib/HuEx-1_0-st-v2/antigenomic.probe_id-only.bgp "
    "-p ../../../regression-data/data/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.pgf "
    "-c ../../../regression-data/data/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.clf "
    "-x 5 "
    "-o " + outdir + " "
    "-m ../../../regression-data/data/lib/HuEx-1_0-st-v2/map.refseq.txt  "
    "--cc-chp-output "
    "-a gc-bg,pm-only,med-polish "
    "-a pm-gcbg,med-polish "
    "-a rma-bg,quant-norm.bioc=true,pm-only,med-polish "
    "-a rma-bg,quant-norm.bioc=true,pm-only,med-polish,pca-select.hard-min=5.min-percent=0.qnorm-only=false.info-criterion=bic "
    "-a rma-bg,quant-norm.bioc=true,pm-only,med-polish,spect-select.metric=angle.log2=false.cut-val=zero.margin=.9.min-percent=0.0.info-criterion=bic ";
  command += Util::joinVectorString(Util::addPrefixSuffix(huexCels, huexCelsPrefix, huexCelsSuffix), " ");

  string chpFilesSuffix = ".chp";
  const char *chpFiles[] = {"huex_wta_cb_A.gc-bg.pm-only.med-polish", "huex_wta_cb_A.pm-gcbg.med-polish", "huex_wta_cb_A.pm-gcbg.plier",
                "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish.pca-select",
                "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_cb_B.gc-bg.pm-only.med-polish",
                "huex_wta_cb_B.pm-gcbg.med-polish", "huex_wta_cb_B.pm-gcbg.plier", "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish",
                "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish.spect-select",
                "huex_wta_cb_C.gc-bg.pm-only.med-polish", "huex_wta_cb_C.pm-gcbg.med-polish", "huex_wta_cb_C.pm-gcbg.plier",
                "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish.pca-select",
                "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_heart_A.gc-bg.pm-only.med-polish",
                "huex_wta_heart_A.pm-gcbg.med-polish", "huex_wta_heart_A.pm-gcbg.plier", "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish",
                "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish.spect-select",
                "huex_wta_heart_B.gc-bg.pm-only.med-polish", "huex_wta_heart_B.pm-gcbg.med-polish", "huex_wta_heart_B.pm-gcbg.plier",
                "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish.pca-select",
                "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_heart_C.gc-bg.pm-only.med-polish",
                "huex_wta_heart_C.pm-gcbg.med-polish", "huex_wta_heart_C.pm-gcbg.plier", "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish",
                "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish.spect-select",
                "huex_wta_muscle_A.gc-bg.pm-only.med-polish", "huex_wta_muscle_A.pm-gcbg.med-polish", "huex_wta_muscle_A.pm-gcbg.plier",
                "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish.pca-select",
                "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_muscle_B.gc-bg.pm-only.med-polish",
                "huex_wta_muscle_B.pm-gcbg.med-polish", "huex_wta_muscle_B.pm-gcbg.plier", "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish",
                "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish.spect-select",
                "huex_wta_muscle_C.gc-bg.pm-only.med-polish", "huex_wta_muscle_C.pm-gcbg.med-polish", "huex_wta_muscle_C.pm-gcbg.plier",
                "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish.pca-select",
                "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_testes_A.gc-bg.pm-only.med-polish",
                "huex_wta_testes_A.pm-gcbg.med-polish", "huex_wta_testes_A.pm-gcbg.plier", "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish",
                "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish.spect-select",
                "huex_wta_testes_B.gc-bg.pm-only.med-polish", "huex_wta_testes_B.pm-gcbg.med-polish", "huex_wta_testes_B.pm-gcbg.plier",
                "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish.pca-select",
                "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_testes_C.gc-bg.pm-only.med-polish",
                "huex_wta_testes_C.pm-gcbg.med-polish", "huex_wta_testes_C.pm-gcbg.plier", "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish",
                "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish.spect-select",
                NULL
};


  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-summarize/HuEx-1_0-st-v2/doPlierWtaRefSeqTest/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles,testDir + "/doPWRST/cc-chp/" , chpFilesSuffix);

  checks.push_back(new CalvinChpCheck(gen, gold));

  checks.push_back(new MatrixCheck(testDir + "/doPWRST/pm-gcbg.plier.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HuEx-1_0-st-v2/doPlierWtaRefSeqTest/pm-gcbg.plier.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck(testDir + "/doPWRST/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HuEx-1_0-st-v2/doPlierWtaRefSeqTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck(testDir + "/doPWRST/rma-bg.quant-norm.pm-only.med-polish.pca-select.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HuEx-1_0-st-v2/doPlierWtaRefSeqTest/rma-bg.quant-norm.pm-only.med-polish.pca-select.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck(testDir + "/doPWRST/rma-bg.quant-norm.pm-only.med-polish.spect-select.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HuEx-1_0-st-v2/doPlierWtaRefSeqTest/rma-bg.quant-norm.pm-only.med-polish.spect-select.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck(testDir + "/doPWRST/pm-gcbg.med-polish.summary.txt",
                                   testDir + "/doPWRST/gc-bg.pm-only.med-polish.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  RegressionTest rmaTest("doPlierWtaRefSeqTest", command.c_str(), checks);
  Verbose::out(1, "Doing doPlierWtaRefSeqTest()");
  rmaTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doPlierWtaRefSeqTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doDabgSubsetTest() {
  string outdir = testDir + "/doDabgSubsetTest";
  string command = "./apt-probeset-summarize "
    "-a pm-only,dabg "
    "-c ../../../regression-data/data/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.clf "
    "-p ../../../regression-data/data/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.pgf "
    "-b ../../../regression-data/data/lib/HuEx-1_0-st-v2/antigenomic.probe_id-only.bgp  "
    "-o " + outdir + " "
    "--use-disk=false "
    "-cc-chp-output "
    "-s ../../../regression-data/data/lib/HuEx-1_0-st-v2/dabg.subset.txt "
    "-x 5 ";
  command += Util::joinVectorString(Util::addPrefixSuffix(huexCels, huexCelsPrefix, huexCelsSuffix), " ");

  string chpFilesSuffix = ".pm-only.dabg.chp";
  const char *chpFiles[] = {"huex_wta_cb_A", "huex_wta_cb_B", "huex_wta_cb_C",
                "huex_wta_heart_A", "huex_wta_heart_B", "huex_wta_heart_C",
                "huex_wta_muscle_A", "huex_wta_muscle_B", "huex_wta_muscle_C",
                "huex_wta_testes_A", "huex_wta_testes_B", "huex_wta_testes_C",
                NULL
  };


  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-summarize/HuEx-1_0-st-v2/doDabgSubsetTest/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, testDir + "/doDabgSubsetTest/cc-chp/", chpFilesSuffix);

  checks.push_back(new CalvinChpCheck(gen, gold));
  checks.push_back(new MatrixCheck(testDir + "/doDabgSubsetTest/pm-only.dabg.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HuEx-1_0-st-v2/doDabgSubsetTest/pm-only.dabg.summary.txt",
                                   0.0001, 1, 1, true, 0));
  RegressionTest dabgTest( "doDabgSubsetTest", command.c_str(), checks);
  dabgTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doDabgSubsetTest()");
  if(!dabgTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doDabgSubsetTest(): " + dabgTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doDabgU133Test() {
  string outdir = testDir + "/doDabgU133Test";
  string command = "./apt-probeset-summarize "
    "-a pm-only,dabg "
    "-b ../../../regression-data/data/lib/HG-U133_Plus_2/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp  "
    "-p ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-x 5 "
    "-o " + outdir + " ";
  command +=  Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest dabgTest("doDabgU133Test", testDir + "/doDabgU133Test/pm-only.dabg.summary.txt",
                          "../../../regression-data/data/chipstream/probeset-summarize/HG-U133_Plus_2/doDabgU133Test/pm-only.dabg.summary.txt",
                          0.0001,
                          command.c_str(),
                          1, 1, false, 0);
  dabgTest.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doDabgU133Test()");
  if(!dabgTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doDabgU133Test(): " + dabgTest.getErrorMsg());
   numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doHumanGeneSpfTest() {
  string outdir = testDir + "/doHumanGeneSpfTest";
  string command = "./apt-probeset-summarize "
    "-a rma-sketch "
    "--spf-file ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.spf "
    "--qc-probesets ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.qcc "
    "-m ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.mps "
    "-o " + outdir + " "
    "--cel-files " + celFileList2 ;
  vector<RegressionCheck *> checks;

  // checks for text matrix files
  checks.push_back(new MatrixCheck(testDir + "/doHumanGeneSpfTest/rma-sketch.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHumanGeneSpfTest/rma-sketch.summary.txt",
                                   0.0001, 1, 0, false, 0));

  RegressionTest test("doHumanGeneSpfTest", command.c_str(), checks);
  test.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doHumanGeneSpfTest()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneSpfTest(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doHumanGeneTest() {
  string outdir = testDir + "/doHGT";
  string command = "./apt-probeset-summarize "
    "-a plier-gcbg-sketch "
    "-a rma-sketch "
    "-a quant-norm,pm-only,med-polish,pca-select "
    "-a quant-norm,pm-only,med-polish,spect-select "
    "-b ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.bgp "
    "-c ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "--qc-probesets ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.qcc "
    "-m ../../../regression-data/data/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.mps "
    "-o " + outdir + " "
    "--use-disk=false "
    "--cel-files " + celFileList2 + " "
    "--feat-effects "
    "--feat-details "
    "--write-sketch "
    "--cc-chp-output";
  vector<RegressionCheck *> checks;

  string chpFilesSuffix = ".chp";
  // checks for chp files.
  const char *chpFiles[] = {"TisMap_Brain_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Brain_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Brain_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Brain_01_v1_WTGene1.rma-sketch",
                "TisMap_Breast_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Breast_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Breast_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Breast_01_v1_WTGene1.rma-sketch",
                "TisMap_Heart_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Heart_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Heart_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Heart_01_v1_WTGene1.rma-sketch",
                "TisMap_Kidney_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Kidney_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Kidney_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Kidney_01_v1_WTGene1.rma-sketch",
                "TisMap_Liver_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Liver_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Liver_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Liver_01_v1_WTGene1.rma-sketch",
                "TisMap_Panc_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Panc_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Panc_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Panc_01_v1_WTGene1.rma-sketch",
                "TisMap_Prost_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Prost_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Prost_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Prost_01_v1_WTGene1.rma-sketch",
                "TisMap_SkMus_01_v1_WTGene1.plier-gcbg-sketch","TisMap_SkMus_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_SkMus_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_SkMus_01_v1_WTGene1.rma-sketch",
                "TisMap_Spleen_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Spleen_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Spleen_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Spleen_01_v1_WTGene1.rma-sketch",
                "TisMap_Testis_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Testis_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Testis_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Testis_01_v1_WTGene1.rma-sketch",
                "TisMap_Thyroid_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Thyroid_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
                "TisMap_Thyroid_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Thyroid_01_v1_WTGene1.rma-sketch",
                NULL
};
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, testDir + "/doHGT/cc-chp/", chpFilesSuffix);

  checks.push_back(new CalvinChpCheck(gen, gold));

  // checks for text matrix files
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/plier-gcbg-sketch.feature-response.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/plier-gcbg-sketch.feature-response.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.pca-select.feature-response.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.pca-select.feature-response.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.spect-select.feature-response.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.spect-select.feature-response.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/rma-sketch.feature-response.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/rma-sketch.feature-response.txt",
                                    0.001, 1, 0, false, 0));


  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/plier-gcbg-sketch.report.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/plier-gcbg-sketch.report.txt",
                                    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/rma-sketch.report.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/rma-sketch.report.txt",
                                    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.pca-select.report.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.pca-select.report.txt",
                                    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.spect-select.report.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.spect-select.report.txt",
                                    0.0001, 1, 1, true, 0));

  checks.push_back(new MixedFileCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt",
                                    0.0001, 0, 0));
  checks.push_back(new MixedFileCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.spect-select.spect-select.report.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.spect-select.spect-select.report.txt",
                                    0.0001, 0, 0));

  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/plier-gcbg-sketch.residuals.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/plier-gcbg-sketch.residuals.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/rma-sketch.residuals.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/rma-sketch.residuals.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.pca-select.residuals.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.pca-select.residuals.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.spect-select.residuals.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.spect-select.residuals.txt",
                                    0.001, 1, 0, false, 0));

  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/plier-gcbg-sketch.summary.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/plier-gcbg-sketch.summary.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/rma-sketch.summary.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/rma-sketch.summary.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.pca-select.summary.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.pca-select.summary.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.spect-select.summary.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.spect-select.summary.txt",
                                    0.0001, 1, 0, false, 0));

  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.normalization-target.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.normalization-target.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/rma-bg.quant-norm.normalization-target.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/rma-bg.quant-norm.normalization-target.txt",
                                    0.0001, 1, 0, false, 0));

  checks.push_back(new MatrixCheck(
                testDir + "/doHGT/quant-norm.pm-only.med-polish.pca-select.self-qnormquant-norm.normalization-target.txt",
                "../../../regression-data/data/chipstream/probeset-summarize/HuGene-1_0-st-v1/doHGT/quant-norm.pm-only.med-polish.pca-select.self-qnormquant-norm.normalization-target.txt",
                                    0.0001, 1, 0, false, 0));

  RegressionTest test("doHumanGeneTest", command.c_str(), checks);
  test.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doHumanGeneTest()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneTest(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}



void ProbeSetSummarizeTest::doSNP6CopyNumberWorkflow1a() {
  string outdir = testDir + "/doSNP6CopyNumberWorkflow1a";
  string command = "./apt-probeset-summarize "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf "
    "--analysis quant-norm.sketch=-1,pm-sum,med-polish,expr.genotype=true.allele-a=true "
    "--use-disk=false "
    "--out-dir " + outdir + " ../../../regression-data/data/cel/GenomeWideSNP_6/NA04626_3X.rep3.CEL ../../../regression-data/data/cel/GenomeWideSNP_6/NA10851_1X_fosmid_ref.rep1.CEL ../../../regression-data/data/cel/GenomeWideSNP_6/NA15510_2X_fosmid.rep2.CEL";

vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
    testDir + "/doSNP6CopyNumberWorkflow1a/quant-norm.pm-sum.med-polish.expr.report.txt",
    "../../../regression-data/data/chipstream/probeset-summarize/GenomeWideSNP_6/doSNP6CopyNumberWorkflow1a/quant-norm.pm-sum.med-polish.expr.report.txt",
    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
    testDir + "/doSNP6CopyNumberWorkflow1a/quant-norm.pm-sum.med-polish.expr.summary.txt",
    "../../../regression-data/data/chipstream/probeset-summarize/GenomeWideSNP_6/doSNP6CopyNumberWorkflow1a/quant-norm.pm-sum.med-polish.expr.summary.txt",
    0.0001, 1, 1, true, 0));
  RegressionTest test("doSNP6CopyNumberWorkflow1a", command.c_str(), checks);
  test.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doSNP6CopyNumberWorkflow1a()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doSNP6CopyNumberWorkflow1a(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doSNP6CopyNumberWorkflow1b() {
  string outdir = testDir + "/doSNP6CopyNumberWorkflow1b";
  string command = "./apt-probeset-summarize "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf "
    "--analysis quant-norm.sketch=-1,pm-only,med-polish,expr.genotype=true "
    "--out-dir " + outdir + " ../../../regression-data/data/cel/GenomeWideSNP_6/NA04626_3X.rep3.CEL ../../../regression-data/data/cel/GenomeWideSNP_6/NA10851_1X_fosmid_ref.rep1.CEL ../../../regression-data/data/cel/GenomeWideSNP_6/NA15510_2X_fosmid.rep2.CEL";

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
    testDir + "/doSNP6CopyNumberWorkflow1b/quant-norm.pm-only.med-polish.expr.report.txt",
    "../../../regression-data/data/chipstream/probeset-summarize/GenomeWideSNP_6/doSNP6CopyNumberWorkflow1b/quant-norm.pm-only.med-polish.expr.report.txt",
    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
    testDir + "/doSNP6CopyNumberWorkflow1b/quant-norm.pm-only.med-polish.expr.summary.txt",
    "../../../regression-data/data/chipstream/probeset-summarize/GenomeWideSNP_6/doSNP6CopyNumberWorkflow1b/quant-norm.pm-only.med-polish.expr.summary.txt",
    0.0001, 1, 1, true, 0));
  RegressionTest test("doSNP6CopyNumberWorkflow1b", command.c_str(), checks);
  test.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doSNP6CopyNumberWorkflow1b()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doSNP6CopyNumberWorkflow1b(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}
void ProbeSetSummarizeTest::doSNP6CopyNumberWorkflow2() {
  string outdir = testDir + "/doSNP6CopyNumberWorkflow2";
  string command = "./apt-probeset-summarize "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf "
    "--analysis quant-norm.sketch=-1,pm-only,med-polish,pca-select "
    "--meta-probesets ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na22.dgv-cnvMay07.mps "
    "--out-dir " + outdir + " ../../../regression-data/data/cel/GenomeWideSNP_6/NA04626_3X.rep3.CEL ../../../regression-data/data/cel/GenomeWideSNP_6/NA10851_1X_fosmid_ref.rep1.CEL ../../../regression-data/data/cel/GenomeWideSNP_6/NA15510_2X_fosmid.rep2.CEL";

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
    testDir + "/doSNP6CopyNumberWorkflow2/quant-norm.pm-only.med-polish.pca-select.summary.txt",
    "../../../regression-data/data/chipstream/probeset-summarize/GenomeWideSNP_6/doSNP6CopyNumberWorkflow2/quant-norm.pm-only.med-polish.pca-select.summary.txt",
    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
    testDir + "/doSNP6CopyNumberWorkflow2/quant-norm.pm-only.med-polish.pca-select.report.txt",
    "../../../regression-data/data/chipstream/probeset-summarize/GenomeWideSNP_6/doSNP6CopyNumberWorkflow2/quant-norm.pm-only.med-polish.pca-select.report.txt",
    0.0001, 1, 1, true, 0));
  checks.push_back(new MixedFileCheck(
    testDir + "/doSNP6CopyNumberWorkflow2/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt",
    "../../../regression-data/data/chipstream/probeset-summarize/GenomeWideSNP_6/doSNP6CopyNumberWorkflow2/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt",
    0.0001, 0, 0));
  RegressionTest test("doSNP6CopyNumberWorkflow2", command.c_str(), checks);
  test.setSuite(*this, outdir, outdir + "/apt-probeset-summarize.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing doSNP6CopyNumberWorkflow2()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doSNP6CopyNumberWorkflow2(): " + test.getErrorMsg());
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
    testDir.setTestDir("chipstream/summarize", true);
    
    ProbeSetSummarizeTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);

    test.parseArgv(argv);
    bool doValgrind = test.doValgrind();
    string database = test.getDatabase();
    cout << "Valgrind: " + ToStr(doValgrind) + " database: " + database << endl;
    test.doRmaNvissaReadWriteFeatureEffectsTest();
    // HuGene-1_0-st-v1 based tests
    test.doHumanGeneTest();
    test.doHumanGeneSpfTest();
    test.doHumanGeneKillListTest();

    // HG-U133_Plus_2 based tests
    test.doRmaTissueTest();
    test.doRmaTissueSpfTest();
    test.doRmaTissueSketchTest();
    test.doRmaTissueTestCC();
    test.doRmaNvissaTest();
    test.doRmaNvissaSketchTest();
    test.doRmaNvissaReadWriteFeatureEffectsTest();


    test.doPlierGcBgTissueTest();
    test.doPlierMMTissueMedianNormTest();
    test.doPlierMMTissueSketchTest();
    test.doDabgU133Test();

    test.doPlierPrecompTissueTest();
    test.doRmaTissueSketchSuppliedTest();
    // HuEx-1_0-st-v2 based tests
    test.doPlierWtaRefSeqTest();
    
    test.doDabgSubsetTest();

    // GenomeWideSNP_6 based tests
    test.doSNP6CopyNumberWorkflow1a();
    test.doSNP6CopyNumberWorkflow1b();
    test.doSNP6CopyNumberWorkflow2();

    /// @todo add tests for --chip-type and --force options, including ones that should fail
    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
