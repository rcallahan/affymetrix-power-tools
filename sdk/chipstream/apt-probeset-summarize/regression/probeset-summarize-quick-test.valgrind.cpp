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

string valgrindString1 = "valgrind --leak-check=yes --log-file=";
string valgrindString2 = ".valgrind.txt ";


string celFileList =" ../../../regression-data/data/idata/p-sum/HG-U133_Plus_2/cel-files.txt";
string celFileList2=" ../../../regression-data/data/idata/p-sum/HuGene-1_0-st-v1/cel-files.txt";

string tissueCelsPrefix = "../../../regression-data/data/idata/cel/HG-U133_Plus_2/";
string tissueCelsSuffix = ".cel";
const char* tissueCels[]={
  "heart-rep1",     "heart-rep2",     "heart-rep3", 
  "hela-rep1",      "hela-rep2",      "hela-rep3",
  "pancrease-rep1", "pancrease-rep2", "pancrease-rep3",
  NULL
};

string tissueCcCelsPrefix = "../../../regression-data/data/idata/cel/HG-U133_Plus_2/";
string tissueCcCelsSuffix = ".agcc.cel";
const char* tissueCcCels[]={
  "heart-rep1",     "heart-rep2",     "heart-rep3",
  "hela-rep1",      "hela-rep2",      "hela-rep3",
  "pancrease-rep1", "pancrease-rep2", "pancrease-rep3",
  NULL
};

string nvissaCelsPrefix = " ../../../regression-data/data/idata/cel/HG-U133_Plus_2/";
string nvissaCelsSuffix = ".CEL";
const char* nvissaCels[]={
  "4221B_a01",	 "4221B_b01",	 "4221B_c01",			
  "4221H_a02",	 "4221H_b02",	 "4221H_c02",			
  "4221HL_a04",	 "4221HL_b04",	 "4221HL_c04",			
  "4221P_a03",	 "4221P_b03",	 "4221P_c03",			
  NULL
};

string huexCelsPrefix = " ../../../regression-data/data/idata/cel/HuEx-1_0-st-v2/";
string huexCelsSuffix = ".CEL";
const char* huexCels[]={
  "huex_wta_cb_A",	 "huex_wta_cb_B",	 "huex_wta_cb_C",	
  "huex_wta_heart_A",	 "huex_wta_heart_B",	 "huex_wta_heart_C",	
  "huex_wta_muscle_A",	 "huex_wta_muscle_B",	 "huex_wta_muscle_C",		
  "huex_wta_testes_A",	 "huex_wta_testes_B",	 "huex_wta_testes_C",	
  NULL
};

class ProbeSetSummarizeTest {

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
  void doSNP6CnWf1a();
  void doSNP6CnWf1b();
  void doSNP6CnWf2();
  void doPlierMMTissueMedianNormTest();
  void doDabgSubsetTest();
  void doDabgU133Test();
};

void ProbeSetSummarizeTest::doHumanGeneKillListTest() {

  if ( !Fs::dirExists("test-generated/doHumanGeneKillListTest") ) {
    Fs::mkdirPath("test-generated/doHumanGeneKillListTest", false);
  }

  string name1 = "qt-doHumanGeneKillListTest-part1"; 
  string name2 = "qt-doHumanGeneKillListTest-part2"; 
  string name3 = "qt-doHumanGeneKillListTest-part3"; 
  string name4 = "qt-doHumanGeneKillListTest-part4"; 
  string name5 = "qt-doHumanGeneKillListTest-part5"; 
  string name6 = "qt-doHumanGeneKillListTest-part6"; 


  // run using kill list, pgf, mps file
  string command1 = valgrindString1 + name1 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "--kill-list ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.txt "
    "-m ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.mps "
    "-o test-generated/doHumanGeneKillListTest/mask1 "
    "--cel-files " + celFileList2;
 
  // compare against pgf w/ probes removed
  string command2 = valgrindString1 + name2 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.pgf "
    "-m ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list2.mps "
    "-o test-generated/doHumanGeneKillListTest/mask2 "
    "--cel-files " + celFileList2;

  // and against CDF masking
  string command3 = valgrindString1 + name3 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.pgf "
    "--kill-list ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.txt "
    "-m ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.mps "
    "-o test-generated/doHumanGeneKillListTest/mask3  "
    "--force "
    "--cel-files " + celFileList2;

  // and when using x/y rather than probe ID
  string command4 = valgrindString1 + name4 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "--kill-list ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.xy.txt "
    "-m ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.mps "
    "-o test-generated/doHumanGeneKillListTest/mask4 "
    "--cel-files " + celFileList2;

  // no mps file
  string command5 = valgrindString1 + name5 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "--kill-list ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list.txt "
    "-o test-generated/doHumanGeneKillListTest/mask5 "
    "--cel-files " + celFileList2;

  // this one should fail -- check to see that if masking was not applied, the data would differ
  string command6 = valgrindString1 + name6 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a quant-norm.sketch=-1,pm-only,med-polish "
    "-c ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "-m ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/test-kill-list2.mps "
    "-o test-generated/doHumanGeneKillListTest/mask6 "
    "--cel-files " + celFileList2;



  RegressionTest test1(
          name1.c_str(),
          "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          "../../../regression-data/data/idata/p-sum/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command1.c_str(), 1, 1, true, 0);
  RegressionTest test2(
          name2.c_str(),
          "test-generated/doHumanGeneKillListTest/mask2/quant-norm.pm-only.med-polish.summary.txt",
          "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command2.c_str(), 1, 1, true, 0);
  RegressionTest test3(
          name3.c_str(),          
          "test-generated/doHumanGeneKillListTest/mask3/quant-norm.pm-only.med-polish.summary.txt",
          "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command3.c_str(), 1, 1, true, 0);
  RegressionTest test4(
          name4.c_str(),
          "test-generated/doHumanGeneKillListTest/mask4/quant-norm.pm-only.med-polish.summary.txt",
          "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command4.c_str(), 1, 1, true, 0);
  RegressionTest test5(
          name5.c_str(),
          "test-generated/doHumanGeneKillListTest/mask5/quant-norm.pm-only.med-polish.summary.txt",
          "../../../regression-data/data/idata/p-sum/doHumanGeneKillListTest/mask5/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command5.c_str(), 1, 1, true, 0);
  RegressionTest test6(
          name6.c_str(),
          "test-generated/doHumanGeneKillListTest/mask6/quant-norm.pm-only.med-polish.summary.txt",
          "../../../regression-data/data/idata/p-sum/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt",
          0.0001, command6.c_str(), 1, 1, true, 0);

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

  // @TODO - Fix this
//   Verbose::out(1, "Doing doHumanGeneKillListTest(): Phase 6 -- expect this to fail");
//   if(!test6.pass()) {
//     numPassed++;
//   }
//   else {
//     Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneKillListTest(): Phase 6: " + test6.getErrorMsg());
//     numFailed++;
//   }

}

void ProbeSetSummarizeTest::doRmaNvissaTest() {
  string name = "qt-doRmaNvissaTest";
  
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doRmaNvissaTest "
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
  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-sum/doRmaNvissaTest/chp/", chpFilesSuffix); 
  gen = Util::addPrefixSuffix(chpFiles, "test-generated/doRmaNvissaTest/chp/", chpFilesSuffix);
  goldCC = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-sum/doRmaNvissaTest/cc-chp/", chpFilesSuffix);
  genCC = Util::addPrefixSuffix(chpFiles, "test-generated/doRmaNvissaTest/cc-chp/", chpFilesSuffix);
  
  checks.push_back(new ChpCheck(gen, gold));
  checks.push_back(new CalvinChpCheck(genCC, goldCC));
  checks.push_back(new MatrixCheck("test-generated/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt", 
                                   "../../../regression-data/data/idata/p-sum/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));

  RegressionTest rmaTest(name.c_str(), command.c_str(), checks);

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
  string name1 = "qt-doRmaNvissaReadWriteFeatureEffectsTest-part1";
  string name2 = "qt-doRmaNvissaReadWriteFeatureEffectsTest-part2";
  string name3 = "qt-doRmaNvissaReadWriteFeatureEffectsTest-part3";
  //First Test 
  string commandMakeFeature = valgrindString1 + name1 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doFE "
    "--feat-effects "; 
    commandMakeFeature += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck("test-generated/doFE/rma-bg.quant-norm.pm-only.med-polish.summary.txt", 
                                   "../../../regression-data/data/idata/p-sum/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.0001,
                                   1, 1, false, 0));

  ///@todo need to tighten up
  checks.push_back(new MatrixCheck("test-generated/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt", 
                                   "../../../regression-data/data/idata/p-sum/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt",
                                   0.001,
                                   1, 1, false, 0));
  RegressionTest rmaMakeTest(name1.c_str(), commandMakeFeature.c_str(), checks);



  //Second Test
  string commandUseFeature = valgrindString1 + name2 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doFE/2 "
    "--use-feat-eff test-generated/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt "; 
    commandUseFeature += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  vector<RegressionCheck *> useChecks;
  ///@todo need to tighten up once we have binary/a5 input/output
  useChecks.push_back(new MatrixCheck("test-generated/doFE/2/rma-bg.quant-norm.pm-only.med-polish.summary.txt", 
                                   "../../../regression-data/data/idata/p-sum/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.001,
                                   1, 1, false, 0));
  RegressionTest rmaUseTest(name2.c_str(), commandUseFeature.c_str(), useChecks);


  //  Third Test - this case checks that code which was designed to use a new style of feature effects file, continues to function
  //  when the old style of feature effects file in text format is used as input.  
  string commandUseOldStyleTextFeature = valgrindString1 + name3 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doFE/3 "
    "--use-feat-eff ../../../regression-data/data/idata/p-sum/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt "; 
    commandUseOldStyleTextFeature += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  vector<RegressionCheck *> useOldStyleTextFeatureChecks;
  useOldStyleTextFeatureChecks.push_back(new MatrixCheck("test-generated/doFE/3/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   "../../../regression-data/data/idata/p-sum/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.001,
                                   1, 1, false, 0));
  RegressionTest rmaUseOldStyleTextFeatureTest(name3.c_str(), commandUseOldStyleTextFeature.c_str(), useOldStyleTextFeatureChecks);

  //  Fourth Test - this case checks that code which was designed to use a new style of feature effects file, continues to function
  //  when the old style of feature effects file in HDF5 format is used as input.  
//   string commandUseOldStyleHDF5Feature = "./apt-probeset-summarize "
//     "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
//     "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
//     "-x 5 "
//     "-o test-generated/doFE/4 "
//     "--a5-feature-effects-input-file ../../../regression-data/data/idata/p-sum/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.a5 "
//     "--a5-feature-effects-input-name  rma-bg.quant-norm.pm-only.med-polish.feature-response "; 
//     commandUseOldStyleHDF5Feature += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

//   vector<RegressionCheck *> useOldStyleHDF5FeatureChecks;
//   useOldStyleHDF5FeatureChecks.push_back(new MatrixCheck("test-generated/doFE/4/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
//                                    "../../../regression-data/data/idata/p-sum/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
//                                    0.001,
//                                    1, 1, false, 0));
//   RegressionTest rmaUseOldStyleHDF5FeatureTest("qt-doRmaNvissaReadWriteFeatureEffectsTest-part4", commandUseOldStyleHDF5Feature.c_str(), useOldStyleHDF5FeatureChecks);

  Verbose::out(1, "Doing doRmaNvissaReadWriteFeatureEffectsTest(), creating feature effects file.  ");
  bool passed = rmaMakeTest.pass();

  Verbose::out(1, "Doing doRmaNvissaReadWriteFeatureEffectsTest(), using the just computed feature effects file. ");
  passed = passed && rmaUseTest.pass();

  Verbose::out(1, "Doing doRmaNvissaReadWriteFeatureEffectsTest(), checking use of old style text format feature effects file.  ");
   passed = passed && rmaUseOldStyleTextFeatureTest.pass();
//   Verbose::out(1, "Doing doRmaNvissaReadWriteFeatureEffectsTest(), checking use of old style HDF5 format feature effects file.  ");
//   passed = passed && rmaUseOldStyleHDF5FeatureTest.pass();

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
  string name = "qt-doRmaNvissaSketchTest"; 
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=100000.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doRmaNvissaSketchTest "; 
    command += Util::joinVectorString(Util::addPrefixSuffix(nvissaCels, nvissaCelsPrefix, nvissaCelsSuffix), " ");

  RegressionTest rmaTest(name.c_str(), "test-generated/doRmaNvissaSketchTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doRmaNvissaSketchTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 225);
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
  string name = "qt-doRmaTissueTest";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doRmaTissueTest "
    "--cel-files " + celFileList;
  
  RegressionTest rmaTest(name.c_str(), "test-generated/doRmaTissueTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doRmaTissueTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
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
  string name = "qt-doRmaTissueSpfTest";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doRmaTissueSpfTest "
    "--cel-files " + celFileList;

  RegressionTest rmaTest(name.c_str(), "test-generated/doRmaTissueSpfTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doRmaTissueSpfTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
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
  string name = "qt-doRmaTissueTestCC";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doRmaTissueTestCC "; 
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCcCels, tissueCcCelsPrefix, tissueCcCelsSuffix), " ");

  RegressionTest rmaTest(name.c_str(), "test-generated/doRmaTissueTestCC/rma-bg.quant-norm.pm-only.med-polish.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doRmaTissueTestCC/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
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
  string name1 = "qt-doPlierPrecompTissueTest-part1";
  string commandMakeFeat = valgrindString1 + name1 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-p ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-b ../../../regression-data/data/idata/lib/HG-U133_Plus_2/pooled-mm-probes.rand-1000-per-bin.bgp "
    "-o test-generated/doPlierPrecompTissueTest "
    "--feat-effects "
    "-x 5 "; 
    commandMakeFeat += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest plierTest(name1.c_str(), "test-generated/doPlierPrecompTissueTest/pm-gcbg.plier.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doPlierPrecompTissueTest/pm-gcbg.plier.summary.txt", 
                         0.0001,
                         commandMakeFeat.c_str(),
                         1, 1, true, 0);
  Verbose::out(1, "Doing doPlierPrecompTissueTest() phase 1");


  string name2 = "qt-doPlierPrecompTissueTest-part2";
  string commandUseFeat = valgrindString1 + name2 + valgrindString2 +  " ./apt-probeset-summarize "
    "-a pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-p ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-b ../../../regression-data/data/idata/lib/HG-U133_Plus_2/pooled-mm-probes.rand-1000-per-bin.bgp "
    "-o test-generated/doPlierPrecompTissueTest/useFeat "
    "--use-feat-eff test-generated/doPlierPrecompTissueTest/pm-gcbg.plier.feature-response.txt "
    "-x 5 ";
  commandUseFeat += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ") ;

  RegressionTest plierPrecompTest(name2.c_str(), "test-generated/doPlierPrecompTissueTest/useFeat/pm-gcbg.plier.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doPlierPrecompTissueTest/useFeat/pm-gcbg.plier.summary.txt",
                         0.0001,
                         commandUseFeat.c_str(),
                         1, 1, true, 0);

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
  string name = "qt-doRmaTissueSketchSuppliedTest"; 
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a quant-norm.sketch=100000.lowprecision=true,pm-only,med-polish.expon=true "
    "--target-sketch ../../../regression-data/data/idata/p-sum/HG-U133_Plus_2/quant-norm.100000.txt "
    "-p ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-o test-generated/doRmaTissueSketchSuppliedTest "
    "-x 6 ";  
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest(name.c_str(), "test-generated/doRmaTissueSketchSuppliedTest/quant-norm.pm-only.med-polish.summary.txt",
                         "../../../regression-data/data/idata/p-sum/doRmaTissueSketchSuppliedTest/quant-norm.pm-only.med-polish.summary.txt",
                         0.021,
                         command.c_str(),
                         1, 1, true, 0);
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
  string name = "qt-doRmaTissueSketchTest";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-bg,quant-norm.sketch=100000.bioc=true.usepm=true,pm-only,med-polish "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf "
    "-x 5 "
    "-o test-generated/doRmaTissueSketchTest "; 
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest(name.c_str(), "test-generated/doRmaTissueSketchTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doRmaTissueSketchTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                         0.01,
                         command.c_str(),
                         1, 1, true, 1306);
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
  string name = "qt-doPlierGcBgTissueTest"; 
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-b ../../../regression-data/data/idata/lib/HG-U133_Plus_2/pooled-mm-probes.rand-1000-per-bin.bgp  "
    "-p ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-x 5 "
    "-o test-generated/doPlierGcBgTissueTest "; 
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest(name.c_str(), "test-generated/doPlierGcBgTissueTest/pm-gcbg.plier.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doPlierGcBgTissueTest/pm-gcbg.plier.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
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
  string name = "qt-doPlierMMTissueSketchTest"; 
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a quant-norm.sketch=100000.lowprecision=true,pm-mm,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-p ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-x 6 "
    "-o test-generated/doPlierMMTissueSketchTest "; 
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest(name.c_str(), "test-generated/doPlierMMTissueSketchTest/quant-norm.pm-mm.plier.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doPlierMMTissueSketchTest/quant-norm.pm-mm.plier.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
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
  string name = "qt-doPlierMMTissueMedianNormTest"; 
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a med-norm.lowprecision=true,pm-mm,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-p ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-x 5 "
    "-o test-generated/doPlierMMTissueMedianNormTest "; 
    command += Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest rmaTest(name.c_str(), "test-generated/doPlierMMTissueMedianNormTest/med-norm.pm-mm.plier.summary.txt", 
                         "../../../regression-data/data/idata/p-sum/doPlierMMTissueMedianNormTest/med-norm.pm-mm.plier.summary.txt",
                         0.0001,
                         command.c_str(),
                         1, 1, true, 0);
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
  string name = "qt-doPlierWtaRefSeqTest"; 
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0 "
    "-b ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/antigenomic.bgp "
    "-p ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.pgf "
    "-c ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.clf "
    "-x 5 "
    "-o test-generated/doRS "
    "-m ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/map.refseq.txt  "
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

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-sum/doRS/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles,"test-generated/doRS/cc-chp/" , chpFilesSuffix);

  checks.push_back(new CalvinChpCheck(gen, gold));

  checks.push_back(new MatrixCheck("test-generated/doRS/pm-gcbg.plier.summary.txt",
                                   "../../../regression-data/data/idata/p-sum/doRS/pm-gcbg.plier.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck("test-generated/doRS/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   "../../../regression-data/data/idata/p-sum/doRS/rma-bg.quant-norm.pm-only.med-polish.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck("test-generated/doRS/rma-bg.quant-norm.pm-only.med-polish.pca-select.summary.txt",
                                   "../../../regression-data/data/idata/p-sum/doRS/rma-bg.quant-norm.pm-only.med-polish.pca-select.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck("test-generated/doRS/rma-bg.quant-norm.pm-only.med-polish.spect-select.summary.txt",
                                   "../../../regression-data/data/idata/p-sum/doRS/rma-bg.quant-norm.pm-only.med-polish.spect-select.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck("test-generated/doRS/pm-gcbg.med-polish.summary.txt",
                                   "test-generated/doRS/gc-bg.pm-only.med-polish.summary.txt",
                                   0.0001,
                                   1, 1, true, 0));
  RegressionTest rmaTest(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing doPlierWtaRefSeqTest()");
  if(!rmaTest.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doPlierWtaRefSeqTest(): " + rmaTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doDabgSubsetTest() {
  string name = "qt-doDabgSubsetTest"; 
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a pm-only,dabg "
    "-c ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.clf "
    "-p ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.pgf "
    "-b ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/antigenomic.bgp  "
    "-o test-generated/doDabgSubsetTest "
    "-cc-chp-output "
    "-s ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/dabg.subset.txt "
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

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-sum/doDabgSubsetTest/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, "test-generated/doDabgSubsetTest/cc-chp/", chpFilesSuffix);

  checks.push_back(new CalvinChpCheck(gen, gold));
  checks.push_back(new MatrixCheck("test-generated/doDabgSubsetTest/pm-only.dabg.summary.txt",
                                   "../../../regression-data/data/idata/p-sum/doDabgSubsetTest/pm-only.dabg.summary.txt",
                                   0.0001, 1, 1, true, 0));
  RegressionTest dabgTest(name.c_str(), command.c_str(), checks);
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
  string name = "qt-doDabgU133Test";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a pm-only,dabg "
    "-b ../../../regression-data/data/idata/lib/HG-U133_Plus_2/pooled-mm-probes.rand-1000-per-bin.bgp  "
    "-p ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "-c ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "-x 5 "
    "-o test-generated/doDabgU133Test ";
  command +=  Util::joinVectorString(Util::addPrefixSuffix(tissueCels, tissueCelsPrefix, tissueCelsSuffix), " ");

  RegressionTest dabgTest(name.c_str(), "test-generated/doDabgU133Test/pm-only.dabg.summary.txt",
                          "../../../regression-data/data/idata/p-sum/doDabgU133Test/pm-only.dabg.summary.txt",
                          0.0001,
                          command.c_str(),
                          1, 1, false, 0);
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
  string name = "qt-doHumanGeneSpfTest";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a rma-sketch "
    "--spf-file ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.spf "
    "--qc-probesets ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.qcc "
    "-m ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.mps "
    "-o test-generated/doHumanGeneSpfTest "
    "--cel-files " + celFileList2 ;
  vector<RegressionCheck *> checks;

  // checks for text matrix files
  checks.push_back(new MatrixCheck("test-generated/doHumanGeneSpfTest/rma-sketch.summary.txt",
                                   "../../../regression-data/data/idata/p-sum/doHumanGeneSpfTest/rma-sketch.summary.txt",
                                   0.0001, 1, 0, false, 0));

  RegressionTest test(name.c_str(), command.c_str(), checks);
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
  string name = "qt-doHumanGeneTest";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "-a plier-gcbg-sketch "
    "-a rma-sketch "
    "-a quant-norm,pm-only,med-polish,pca-select "
    "-a quant-norm,pm-only,med-polish,spect-select "
    "-b ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.bgp "
    "-c ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf "
    "-p ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.pgf "
    "--qc-probesets ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.qcc "
    "-m ../../../regression-data/data/idata/lib/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.mps "
    "-o test-generated/doHGT "
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

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-sum/doHGT/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, "test-generated/doHGT/cc-chp/", chpFilesSuffix);

  checks.push_back(new CalvinChpCheck(gen, gold));

  // checks for text matrix files
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/plier-gcbg-sketch.feature-response.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/plier-gcbg-sketch.feature-response.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.pca-select.feature-response.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.pca-select.feature-response.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.spect-select.feature-response.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.spect-select.feature-response.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/rma-sketch.feature-response.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/rma-sketch.feature-response.txt",
                                    0.001, 1, 0, false, 0));


  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/plier-gcbg-sketch.report.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/plier-gcbg-sketch.report.txt",
                                    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/rma-sketch.report.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/rma-sketch.report.txt",
                                    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.pca-select.report.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.pca-select.report.txt",
                                    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.spect-select.report.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.spect-select.report.txt",
                                    0.0001, 1, 1, true, 0));

  checks.push_back(new MixedFileCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt",
                                    0.0001, 0, 0));
  checks.push_back(new MixedFileCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.spect-select.spect-select.report.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.spect-select.spect-select.report.txt",
                                    0.0001, 0, 0));

  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/plier-gcbg-sketch.residuals.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/plier-gcbg-sketch.residuals.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/rma-sketch.residuals.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/rma-sketch.residuals.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.pca-select.residuals.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.pca-select.residuals.txt",
                                    0.001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.spect-select.residuals.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.spect-select.residuals.txt",
                                    0.001, 1, 0, false, 0));

  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/plier-gcbg-sketch.summary.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/plier-gcbg-sketch.summary.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/rma-sketch.summary.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/rma-sketch.summary.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.pca-select.summary.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.pca-select.summary.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.spect-select.summary.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.spect-select.summary.txt",
                                    0.0001, 1, 0, false, 0));

  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.normalization-target.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.normalization-target.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/rma-bg.quant-norm.normalization-target.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/rma-bg.quant-norm.normalization-target.txt",
                                    0.0001, 1, 0, false, 0));
  checks.push_back(new MatrixCheck(
                "test-generated/doHGT/quant-norm.pm-only.med-polish.pca-select.self-qnormquant-norm.normalization-target.txt",
                "../../../regression-data/data/idata/p-sum/doHGT/quant-norm.pm-only.med-polish.pca-select.self-qnormquant-norm.normalization-target.txt",
                                    0.0001, 1, 0, false, 0));
  RegressionTest test("qt-doHumanGeneTest", command.c_str(), checks);
  Verbose::out(1, "Doing doHumanGeneTest()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doHumanGeneTest(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}



void ProbeSetSummarizeTest::doSNP6CnWf1a() {
  string name = "qt-doSNP6CnWf1a";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6_cn/GenomeWideSNP_6.spf "
    "--analysis quant-norm.sketch=-1,pm-sum,med-polish,expr.genotype=true.allele-a=true "
    "--out-dir test-generated/doSNP6CnWf1a ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA04626_3X.rep3.CEL ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA10851_1X_fosmid_ref.rep1.CEL ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA15510_2X_fosmid.rep2.CEL";
  
vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
    "test-generated/doSNP6CnWf1a/quant-norm.pm-sum.med-polish.expr.report.txt",
    "../../../regression-data/data/idata/p-sum/doSNP6CnWf1a/quant-norm.pm-sum.med-polish.expr.report.txt",
    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
    "test-generated/doSNP6CnWf1a/quant-norm.pm-sum.med-polish.expr.summary.txt",
    "../../../regression-data/data/idata/p-sum/doSNP6CnWf1a/quant-norm.pm-sum.med-polish.expr.summary.txt",
    0.0001, 1, 1, true, 0));
  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing doSNP6CnWf1a()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doSNP6CnWf1a(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetSummarizeTest::doSNP6CnWf1b() {
  string name = "qt-doSNP6CnWf1b";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6_cn/GenomeWideSNP_6.spf "
    "--analysis quant-norm.sketch=-1,pm-only,med-polish,expr.genotype=true "
    "--out-dir test-generated/doSNP6CnWf1b ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA04626_3X.rep3.CEL ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA10851_1X_fosmid_ref.rep1.CEL ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA15510_2X_fosmid.rep2.CEL";

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
    "test-generated/doSNP6CnWf1b/quant-norm.pm-only.med-polish.expr.report.txt",
    "../../../regression-data/data/idata/p-sum/doSNP6CnWf1b/quant-norm.pm-only.med-polish.expr.report.txt",
    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
    "test-generated/doSNP6CnWf1b/quant-norm.pm-only.med-polish.expr.summary.txt",
    "../../../regression-data/data/idata/p-sum/doSNP6CnWf1b/quant-norm.pm-only.med-polish.expr.summary.txt",
    0.0001, 1, 1, true, 0));
  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing doSNP6CnWf1b()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doSNP6CnWf1b(): " + test.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}
void ProbeSetSummarizeTest::doSNP6CnWf2() {
  string name = "qt-doSNP6CnWf2";
  string command = valgrindString1 + name + valgrindString2 +  " ./apt-probeset-summarize "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6_cn/GenomeWideSNP_6.spf "
    "--analysis quant-norm.sketch=-1,pm-only,med-polish,pca-select "
    "--meta-probesets ../../../regression-data/data/idata/lib/GenomeWideSNP_6_cn/GenomeWideSNP_6.na22.dgv-cnvMay07.mps "
    "--out-dir test-generated/doSNP6CnWf2 ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA04626_3X.rep3.CEL ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA10851_1X_fosmid_ref.rep1.CEL ../../../regression-data/data/idata/cel/GenomeWideSNP_6_cn/NA15510_2X_fosmid.rep2.CEL";

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
    "test-generated/doSNP6CnWf2/quant-norm.pm-only.med-polish.pca-select.summary.txt",
    "../../../regression-data/data/idata/p-sum/doSNP6CnWf2/quant-norm.pm-only.med-polish.pca-select.summary.txt",
    0.0001, 1, 1, true, 0));
  checks.push_back(new MatrixCheck(
    "test-generated/doSNP6CnWf2/quant-norm.pm-only.med-polish.pca-select.report.txt",
    "../../../regression-data/data/idata/p-sum/doSNP6CnWf2/quant-norm.pm-only.med-polish.pca-select.report.txt",
    0.0001, 1, 1, true, 0));
  checks.push_back(new MixedFileCheck(
    "test-generated/doSNP6CnWf2/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt",
    "../../../regression-data/data/idata/p-sum/doSNP6CnWf2/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt",
    0.0001, 0, 0));
  RegressionTest test(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing doSNP6CnWf2()");
  if(!test.pass()) {
    Verbose::out(1, "Error in ProbeSetSummarizeTest::doSNP6CnWf2(): " + test.getErrorMsg());
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
    testDir.setTestDir("chipstream/probeset-summarize-valgrind", true);
    
    ProbeSetSummarizeTest test;
    test.testDir = testDir.asString();
    valgrindString1 = valgrindString1 + testDir.asString() + "/";
    
    Verbose::setLevel(2);

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
    test.doSNP6CnWf1a();
    test.doSNP6CnWf1b();
    test.doSNP6CnWf2(); 

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
