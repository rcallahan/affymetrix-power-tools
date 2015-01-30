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

//
#include "chipstream/apt-probeset-genotype/regression/TrustedProbesCheck.h"
//

#include "util/CalvinChpCheck.h"
#include "util/ChpCheck.h"
#include "util/Convert.h"
#include "util/FsTestDir.h"
#include "util/RegressionSuite.h"
#include "util/RegressionTest.h"
#include "util/TextFileCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string.h>
#include <string>
#include <vector>
//
#include "util/Fs.h"

using namespace std;

const char *sty20cels[] = {
  "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
  "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
  "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
  "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
  "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
  "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
  "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
  NULL
};

class ProbeSetGenotypeTest : public RegressionSuite {

public:
  int numPassed, numFailed;
  std::string testDir;
  ProbeSetGenotypeTest() {
    numPassed = 0;
    numFailed = 0;
  }

  void doBirdseedSnp6Spf();
  void doBirdseedSnp6Full();
  void doBirdseedSnp6Chp();
  void doBirdseed2Snp6Chp();
  void doBrlmmpSnp6Chp();
  void doBrlmmpSnp5Chp();
  void doBrlmmpSnp5Spf();
  void doBrlmmpSnp5KillList();
  void doMultiChannelAxiom();
  void doReadGenotypesIn();
  void doChpFiles();
  void doSpf();
  void doWritePriorFileStyCCS_2_20_1_0();
  void doStyCCS_2_40_0_9_0_033();
  void doStyCCS_2_20_1_0();
  void doStyCES_2_20_1_0();
  void doStyRVT_2_20_1_0();
  void doStyMva_2_20_1_2();
  void doStyMva_2_20_1_0();
  void doStyMva_2_20_0_8_0();
  void doFileStyCCS_2_20_1_0();
  void doSet1();
  void doLabelZStyCCS();
  void doLabelZStyTestPriors();
  void doLabelZStyCCSseveral();
  void doLabelZStyHints();
  void doTrustedProbes();
  string createString(const char **items, int size, const char *prefix, const char *suffix);
  void doBrlmmpSnp5CdfReadWriteFeatureEffectsA5();


};

string ProbeSetGenotypeTest::createString(const char **items, int size, const char *prefix, const char *suffix) {
    string val = "";
    for(int i = 0; i < size; i++) {
        val = val + prefix + items[i] + suffix;
    }
    return val;
}

void ProbeSetGenotypeTest::doBirdseedSnp6Spf() {
  string outdir = testDir + "/doBirdseedSnp6Spf";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "--analysis birdseed "
    "--use-disk=false "
    "--special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--spf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.spf "
    "--read-models-birdseed ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o " + outdir + " "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/fas_cel_files.txt";

  string name = "doBirdseedSnp6Spf";
  vector<RegressionCheck *> checks;

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck(testDir + "/doBirdseedSnp6Spf/birdseed.report.txt",
  //                                 "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(testDir + "/doBirdseedSnp6Spf/birdseed.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.calls.txt",
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck(testDir + "/doBirdseedSnp6Spf/birdseed.confidences.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.confidences.txt",
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
  birdseedTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!birdseedTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + birdseedTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doBirdseedSnp6Full() {
  string outdir = testDir + "/doBirdseedSnp6Full";
  string command = "./apt-probeset-genotype "
    "--cc-chp-output "
    "--table-output "
    "--analysis birdseed "
    "--special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.Full.specialSNPs "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.Full.cdf "
    "--read-models-birdseed ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o " + outdir + " " 
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/fas_cel_files.txt";

  const char *chpFiles[] = {"NA06985_GW6_C.birdseed", "NA06993_GW6_C.birdseed", "NA07000_GW6_C.birdseed",
                "NA07019_GW6_C.birdseed", "NA07029_GW6_C.birdseed", "NA07056_GW6_C.birdseed",
                "NA07345_GW6_C.birdseed", "NA10830_GW6_C.birdseed", "NA10831_GW6_C.birdseed",
                "NA10839_GW6_C.birdseed", "NA10855_GW6_C.birdseed", "NA10857_GW6_C.birdseed",
                "NA10860_GW6_C.birdseed", "NA10863_GW6_C.birdseed", "NA11881_GW6_C.birdseed",
                "NA11882_GW6_C.birdseed", "NA11993_GW6_C.birdseed", "NA12003_GW6_C.birdseed",
                "NA12004_GW6_C.birdseed",
                NULL
};

  string name = "doBirdseedSnp6Full";
  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Full/cc-chp/", ".chp");
  gen = Util::addPrefixSuffix(chpFiles, testDir + "/doBirdseedSnp6Full/cc-chp/", ".chp");

  // checks for chp files.
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001));

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck(testDir + "/doBirdseedSnp6Full/birdseed.report.txt",
  //                                 "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Full/birdseed.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(testDir + "/doBirdseedSnp6Full/birdseed.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Full/birdseed.calls.txt",
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck(testDir + "/doBirdseedSnp6Full/birdseed.confidences.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Full/birdseed.confidences.txt",
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
  birdseedTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!birdseedTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + birdseedTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doBirdseedSnp6Chp() {
  string outdir = testDir + "/doBirdseedSnp6Chp";
  string command = "./apt-probeset-genotype "
    "--cc-chp-output "
    "--table-output "
    "--analysis birdseed "
    "--use-disk=false "
    "--special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf "
    "--read-models-birdseed ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o " + outdir + " "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/fas_cel_files.txt";

  const char *chpFiles[] = {"NA06985_GW6_C.birdseed", "NA06993_GW6_C.birdseed", "NA07000_GW6_C.birdseed",
                "NA07019_GW6_C.birdseed", "NA07029_GW6_C.birdseed", "NA07056_GW6_C.birdseed",
                "NA07345_GW6_C.birdseed", "NA10830_GW6_C.birdseed", "NA10831_GW6_C.birdseed",
                "NA10839_GW6_C.birdseed", "NA10855_GW6_C.birdseed", "NA10857_GW6_C.birdseed",
                "NA10860_GW6_C.birdseed", "NA10863_GW6_C.birdseed", "NA11881_GW6_C.birdseed",
                "NA11882_GW6_C.birdseed", "NA11993_GW6_C.birdseed", "NA12003_GW6_C.birdseed",
                "NA12004_GW6_C.birdseed",
                NULL
  };

  string name = "doBirdseedSnp6Chp";
  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Chp/cc-chp/", ".chp");
  gen = Util::addPrefixSuffix(chpFiles, testDir + "/doBirdseedSnp6Chp/cc-chp/", ".chp");

  // checks for chp files.
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001));
  // checks for report file
  checks.push_back(new MixedFileCheck(testDir + "/doBirdseedSnp6Chp/birdseed.report.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Chp/birdseed.report.txt",
                                   0.0001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(testDir + "/doBirdseedSnp6Chp/birdseed.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Chp/birdseed.calls.txt",
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck(testDir + "/doBirdseedSnp6Chp/birdseed.confidences.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseedSnp6Chp/birdseed.confidences.txt",
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
  birdseedTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!birdseedTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + birdseedTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doBirdseed2Snp6Chp() {
  string outdir = testDir + "/doBirdseed2Snp6Chp";
  string command = "./apt-probeset-genotype "
    "--cc-chp-output "
    "--table-output "
    "--analysis birdseed-v2 "
    "--special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf "
    "--read-models-birdseed ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o " + outdir + " "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/fas_cel_files.txt "
    "--set-gender-method cn-probe-chrXY-ratio "
    "--chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes "
    "--chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes ";

  const char *chpFiles[] = {"NA06985_GW6_C", "NA06993_GW6_C", "NA07000_GW6_C",
                "NA07019_GW6_C", "NA07029_GW6_C", "NA07056_GW6_C",
                "NA07345_GW6_C", "NA10830_GW6_C", "NA10831_GW6_C",
                "NA10839_GW6_C", "NA10855_GW6_C", "NA10857_GW6_C",
                "NA10860_GW6_C", "NA10863_GW6_C", "NA11881_GW6_C",
                "NA11882_GW6_C", "NA11993_GW6_C", "NA12003_GW6_C",
                "NA12004_GW6_C",
                NULL
  };

  string name = "doBirdseed2Snp6Chp";
  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseed2Snp6Chp/cc-chp/", ".birdseed-v2.chp");
  gen = Util::addPrefixSuffix(chpFiles, testDir + "/doBirdseed2Snp6Chp/cc-chp/", ".birdseed-v2.chp");


  // checks for chp files.
  ///@todo we need to relax Eps due to insignificant difference in extra (non-call, non-confidence) values for one SNP between win32 and other platforms.
  ///      the interface does not (yet) allow us to set different eps for each field -- which would be ideal. With this change we are not
  ///      doing as tight a check on confidences as we probably should. Of course the test on the matrix file mitigates some risk here.
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001));
  // checks for report file
  checks.push_back(new MixedFileCheck(testDir + "/doBirdseed2Snp6Chp/birdseed-v2.report.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseed2Snp6Chp/birdseed-v2.report.txt",
                                   0.0001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(testDir + "/doBirdseed2Snp6Chp/birdseed-v2.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseed2Snp6Chp/birdseed-v2.calls.txt",
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck(testDir + "/doBirdseed2Snp6Chp/birdseed-v2.confidences.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBirdseed2Snp6Chp/birdseed-v2.confidences.txt",
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
  birdseedTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");

  Verbose::out(1, "Doing " + name + "()");
  if(!birdseedTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + birdseedTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doBrlmmpSnp6Chp() {
  string outdir = testDir + "/doBrlmmpSnp6Chp";
  string command = "./apt-probeset-genotype "
    "--cc-chp-output "
    "--table-output "
    "--analysis brlmm-p "
    "--use-disk=false "
    "--special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf "
    "--read-models-brlmmp  ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.brlmm-p.models  "
    "-o " + outdir + " "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/fas_cel_files.txt";

  const char *chpFiles[] = {"NA06985_GW6_C.brlmm-p", "NA06993_GW6_C.brlmm-p", "NA07000_GW6_C.brlmm-p",
                "NA07019_GW6_C.brlmm-p", "NA07029_GW6_C.brlmm-p", "NA07056_GW6_C.brlmm-p",
                "NA07345_GW6_C.brlmm-p", "NA10830_GW6_C.brlmm-p", "NA10831_GW6_C.brlmm-p",
                "NA10839_GW6_C.brlmm-p", "NA10855_GW6_C.brlmm-p", "NA10857_GW6_C.brlmm-p",
                "NA10860_GW6_C.brlmm-p", "NA10863_GW6_C.brlmm-p", "NA11881_GW6_C.brlmm-p",
                "NA11882_GW6_C.brlmm-p", "NA11993_GW6_C.brlmm-p", "NA12003_GW6_C.brlmm-p",
                "NA12004_GW6_C.brlmm-p",
                NULL
  };

  string name = "doBrlmmpSnp6Chp";
  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBrlmmpSnp6Chp/cc-chp/", ".chp");
  gen = Util::addPrefixSuffix(chpFiles, testDir + "/doBrlmmpSnp6Chp/cc-chp/", ".chp");

  // checks for chp files.
  checks.push_back(new CalvinChpCheck(gen,gold, 0, L"apt-", 0.001,false));
  // checks for report file
  ///@todo relaxed due to difference in chr-x-het-rate metric on win32
  checks.push_back(new MixedFileCheck(testDir + "/doBrlmmpSnp6Chp/brlmm-p.report.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBrlmmpSnp6Chp/brlmm-p.report.txt",
                                   0.0001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp6Chp/brlmm-p.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBrlmmpSnp6Chp/brlmm-p.calls.txt",
                                   0.0000001,
                                   1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks.push_back(new MatrixCheck("../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_6/doBrlmmpSnp6Chp/brlmm-p.confidences.txt",
                                   testDir + "/doBrlmmpSnp6Chp/brlmm-p.confidences.txt",
                                   0.0001,
                                   1, 1, false, 0));
  RegressionTest brlmmpTest(name.c_str(), command.c_str(), checks);
  brlmmpTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");

  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmpTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmpTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doBrlmmpSnp5Spf() {
  string outdir = testDir + "/doBrlmmpSnp5Spf";
  std::string command = "./apt-probeset-genotype "
    "--spf-file ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.spf "
    "--chrX-snps ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir " + outdir;

  const char *chpFiles[] = {
    "NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C","NA06994_GW5_C","NA07000_GW5_C",
    "NA07019_GW5_C","NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C","NA07048_GW5_C",
    "NA07055_GW5_C","NA07056_GW5_C","NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
    "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C","NA10838_GW5_C","NA10839_GW5_C",
    NULL
  };

  string name = "doBrlmmpSnp5Spf";
  vector<RegressionCheck *> checks;


  command += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/cel/GenomeWideSNP_5/", ".CEL"), " ");

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5Spf/brlmm-p.report.txt",
  //                                 "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Spf/brlmm-p.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5Spf/brlmm-p.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Spf/brlmm-p.calls.txt",
                                   0.0000001,
                                   1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5Spf/brlmm-p.confidences.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Spf/brlmm-p.confidences.txt",
                                   0.0001,
                                   1, 1, false, 0));
  RegressionTest brlmmpTest(name.c_str(), command.c_str(), checks);
  brlmmpTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");

  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmpTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmpTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}


/*  This is the only test within the genotyping regression suite which tests feature effects.  It uses a cdf file as input and does a write and read of feature effects in two stages.  At both steps it validates that the feature effects are correct by checking that the summary values created using the feature effects match the golden values.  It would be nice to check the feature effects themselves but for this testcase they are in A5 format. The A5 format was chosen since we do not even have an testcase in the summarization regression suite which tests the A5 input/output of feature effects */

void ProbeSetGenotypeTest::doBrlmmpSnp5CdfReadWriteFeatureEffectsA5() {
  string outdir = testDir + "/doBrlmmpSnp5CdfWriteFeatureEffectsA5";
  string commandMakeFeatureEffects = "./apt-probeset-genotype "
    "--qmethod-spec med-polish "
    "--analysis brlmm-p "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.cdf "
    "--chrX-snps ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--out-dir " + outdir + "  "
    "--a5-feature-effects "
    "--summaries";

  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
                "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
                "NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C",
                "NA07048_GW5_C","NA07055_GW5_C","NA07056_GW5_C",
                "NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
                "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C",
                "NA10838_GW5_C","NA10839_GW5_C",
                NULL
  };

  commandMakeFeatureEffects += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/cel/GenomeWideSNP_5/", ".CEL"), " ");

  //First Test  Feature effects are created and written into A5 format.
  vector<RegressionCheck *> makeChecks;

  makeChecks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5CdfWriteFeatureEffectsA5/brlmm-p.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.summary.txt",
                                   0.00001,
                                   1, 1, false, 0));

  RegressionTest makeFeatureEffectsTest("doBrlmmpSnp5CdfReadWriteFeatureEffectsA5-part1", commandMakeFeatureEffects.c_str(), makeChecks);
  makeFeatureEffectsTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");

  //Second Test  The feature effects just created are read in and used to create summary values.
  //TODO  This test presently does not read in the just created feature effects file. This is because the present version of the FE file
  //      needs to have suffixes from the probeset names removed.  The final solution will be to remove the addition of the suffixes
  //      from the QuantMethods.  An interim solution would be to have an HDF5 comparison of the present FE files.  Or as is done now,
  //      just drop the check that the FE file was correctly written.  This second check still checks that translated FE files are
  //      correctly read.
  string outdir2 = testDir + "/doBrlmmpSnp5CdfReadFeatureEffectsA5";
  string commandUseFeatureEffects = "./apt-probeset-genotype "
    "--qmethod-spec med-polish "
    "--analysis brlmm-p "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.cdf "
    "--chrX-snps ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
//    "--a5-feature-effects-input-file " + testDir + "/doBrlmmpSnp5CdfWriteFeatureEffectsA5/brlmm-p.feature-response.a5 "
    "--a5-feature-effects-input-file ../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.feature-response.a5.translated "
    "--a5-feature-effects-input-name  brlmm-p.feature-response "
    "--out-dir " + outdir2 + " "
    "--summaries";

  commandUseFeatureEffects += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/cel/GenomeWideSNP_5/", ".CEL"), " ");



  vector<RegressionCheck *> useChecks;

  useChecks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5CdfReadFeatureEffectsA5/brlmm-p.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.summary.txt",
                                   0.00002,
                                   1, 1, false, 0));


  RegressionTest useFeatureEffectsTest("doBrlmmpSnp5CdfReadWriteFeatureEffectsA5-part2", commandUseFeatureEffects.c_str(), useChecks);
  useFeatureEffectsTest.setSuite(*this, outdir2, outdir2 + "/apt-probeset-genotype.log", outdir2 + "/valgrind.log");

  Verbose::out(1, "Doing doBrlmmpSnp5CdfReadWriteFeatureEffectsA5Test(), creating feature effects file.  ");
  bool passed = makeFeatureEffectsTest.pass();

  Verbose::out(1, "Doing doBrlmmpSnp5CdfReadWriteFeatureEffectsA5Test(), using the just computed feature effects file. ");
  passed = passed && useFeatureEffectsTest.pass();

  if(!passed) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::doBrlmmpSnp5CdfReadWriteFeatureEffectsA5Test(): " + makeFeatureEffectsTest.getErrorMsg()
                 + " " + useFeatureEffectsTest.getErrorMsg());
    numFailed++;

  }
  else {
    numPassed++;
  }

}




void ProbeSetGenotypeTest::doBrlmmpSnp5Chp() {
  string outdir = testDir + "/doBrlmmpSnp5Chp";
  string command = "./apt-probeset-genotype "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.cdf "
    "--chrX-snps ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir " + outdir + " "
    "--use-disk=false "
    "--summaries "
    "--cc-chp-output";

  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
                "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
                "NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C",
                "NA07048_GW5_C","NA07055_GW5_C","NA07056_GW5_C",
                "NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
                "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C",
                "NA10838_GW5_C","NA10839_GW5_C",
                NULL
  };

  string name = "doBrlmmpSnp5Chp";
  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Chp/cc-chp/", ".brlmm-p.chp");
  gen = Util::addPrefixSuffix(chpFiles,testDir + "/doBrlmmpSnp5Chp/cc-chp/", ".brlmm-p.chp");
  command += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/cel/GenomeWideSNP_5/", ".CEL"), " ");

  // checks for chp files. [ignore header because of command string]
  checks.push_back(new CalvinChpCheck(gen,gold, 0, L"apt-", 0.001,false));
  // checks for report file
  checks.push_back(new MixedFileCheck(testDir + "/doBrlmmpSnp5Chp/brlmm-p.report.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.report.txt",
                                   0.00001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5Chp/brlmm-p.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.calls.txt",
                                   0.0000001,
                                   1, 1, false, 0));
  // checks for the table files.
  ///@todo Is this too loose?
  checks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5Chp/brlmm-p.normalized-summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.normalized-summary.txt",
                                   0.001,
                                   1, 2, false, 0));
  ///@todo Is this too loose?
  checks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5Chp/brlmm-p.summary.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.summary.txt",
                                   0.001,
                                   1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks.push_back(new MatrixCheck(testDir + "/doBrlmmpSnp5Chp/brlmm-p.confidences.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.confidences.txt",
                                   0.0001,
                                   1, 1, false, 0));
  RegressionTest brlmmpTest(name.c_str(), command.c_str(), checks);
  brlmmpTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmpTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmpTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doBrlmmpSnp5KillList() {
  string outdir1 = testDir + "/doBrlmmpSnp5KillList";
  string command1 ="./apt-probeset-genotype "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.cdf "
    "--kill-list ../../../regression-data/data/lib/GenomeWideSNP_5/test-kill-list.txt "
    "-s ../../../regression-data/data/lib/GenomeWideSNP_5/test-kill-list.snp-list.txt "
    "--chrX-snps ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir " + outdir1 + " "
    "--summaries ";
  string outdir2 = testDir + "/doBrlmmpSnp5KillList/mask2";
  string command2 = "./apt-probeset-genotype "
    "--force "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_5/test-kill-list.cdf "
    "-s ../../../regression-data/data/lib/GenomeWideSNP_5/test-kill-list.snp-list.txt "
    "--chrX-snps ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir " + outdir2 + " "
    "--summaries ";
  string outdir3 = testDir + "/doBrlmmpSnp5KillList/mask3";
  string command3 = "./apt-probeset-genotype "
    "--cdf-file ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.cdf "
    "-s ../../../regression-data/data/lib/GenomeWideSNP_5/test-kill-list.snp-list.txt "
    "--chrX-snps ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir " + outdir3 + " "
    "--summaries ";

  //const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C","NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C","NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C","NA07048_GW5_C","NA07055_GW5_C","NA07056_GW5_C","NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C","NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C","NA10838_GW5_C","NA10839_GW5_C"};
  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
                "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
                NULL
  };

  string name = "doBrlmmpSnp5KillList";
  vector<RegressionCheck *> checks1;
  vector<RegressionCheck *> checks2;
  vector<RegressionCheck *> checks3;

  command1 += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/cel/GenomeWideSNP_5/", ".CEL"), " ");
  command2 += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/cel/GenomeWideSNP_5/", ".CEL"), " ");
  command3 += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/cel/GenomeWideSNP_5/", ".CEL"), " ");

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  //checks2.push_back(new MatrixCheck(
  //                            testDir + "/doBrlmmpSnp5KillList/mask2/brlmm-p.report.txt",
  //                            testDir + "/doBrlmmpSnp5KillList/brlmm-p.report.txt",
  //                            0.00001,
  //                            1,2, false, 0));
  checks2.push_back(new MatrixCheck(
                              testDir + "/doBrlmmpSnp5KillList/mask2/brlmm-p.calls.txt",
                              testDir + "/doBrlmmpSnp5KillList/brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, false, 0));
  ///@todo Is this too loose?
  checks2.push_back(new MatrixCheck(
                              testDir + "/doBrlmmpSnp5KillList/mask2/brlmm-p.normalized-summary.txt",
                              testDir + "/doBrlmmpSnp5KillList/brlmm-p.normalized-summary.txt",
                              0.001,
                              1, 2, false, 0));
  ///@todo Is this too loose?
  checks2.push_back(new MatrixCheck(
                              testDir + "/doBrlmmpSnp5KillList/mask2/brlmm-p.summary.txt",
                              testDir + "/doBrlmmpSnp5KillList/brlmm-p.summary.txt",
                              0.001,
                              1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks2.push_back(new MatrixCheck(
                              testDir + "/doBrlmmpSnp5KillList/mask2/brlmm-p.confidences.txt",
                              testDir + "/doBrlmmpSnp5KillList/brlmm-p.confidences.txt",
                              0.0001,
                              1, 1, false, 0));

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  //checks3.push_back(new MatrixCheck(
  //                            testDir + "/doBrlmmpSnp5KillList/mask3/brlmm-p.report.txt",
  //                            testDir + "/doBrlmmpSnp5KillList/brlmm-p.report.txt",
  //                            0.00001,
  //                            1,2, false, 0));
  checks3.push_back(new MatrixCheck(
                              testDir + "/doBrlmmpSnp5KillList/mask3/brlmm-p.calls.txt",
                              testDir + "/doBrlmmpSnp5KillList/brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, false, 0));
  checks3[checks3.size()-1]->m_NegTest = true;
  ///@todo Is this too loose?
  checks3.push_back(new MatrixCheck(
                              testDir + "/doBrlmmpSnp5KillList/mask3/brlmm-p.normalized-summary.txt",
                              testDir + "/doBrlmmpSnp5KillList/brlmm-p.normalized-summary.txt",
                              0.001,
                              1, 2, false, 0));
  checks3[checks3.size()-1]->m_NegTest = true;
  ///@todo Is this too loose?
  checks3.push_back(new MatrixCheck(
                              testDir + "/doBrlmmpSnp5KillList/mask3/brlmm-p.summary.txt",
                              testDir + "/doBrlmmpSnp5KillList/brlmm-p.summary.txt",
                              0.001,
                              1, 1, false, 0));
  checks3[checks3.size()-1]->m_NegTest = true;
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks3.push_back(new MatrixCheck(
                              testDir + "/doBrlmmpSnp5KillList/mask3/brlmm-p.confidences.txt",
                              testDir + "/doBrlmmpSnp5KillList/brlmm-p.confidences.txt",
                              0.0001,
                              1, 1, false, 0));
  checks3[checks3.size()-1]->m_NegTest = true;

  string subName;
  subName = name + "-part1";
  RegressionTest test1(subName.c_str(), command1.c_str(), checks1);
  test1.setSuite(*this, outdir1,  outdir1 + "/apt-probeset-genotype.log", outdir1 + "/valgrind.log");
  subName = name + "-part2";
  RegressionTest test2(subName.c_str(), command2.c_str(), checks2);
  test2.setSuite(*this, outdir2, outdir2 + "/apt-probeset-genotype.log", outdir2 + "/valgrind.log");
  subName = name + "-part3";
  RegressionTest test3(subName.c_str(), command3.c_str(), checks3);
  test3.setSuite(*this, outdir3, outdir3 + "/apt-probeset-genotype.log", outdir3 + "/valgrind.log");

  Verbose::out(1, "Doing " + name + "(): Phase 1");
  if(!test1.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): Phase 1: " + test1.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

  Verbose::out(1, "Doing " + name + "(): Phase 2");
  if(!test2.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): Phase 2: " + test2.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

  Verbose::out(1, "Doing " + name + "(): Phase 3 (expect differences -- error if same)");
  if(!test3.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): Phase 2: " + test2.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }

}

void ProbeSetGenotypeTest::doMultiChannelAxiom() {
  const char *file_prefixes[] = {
    "NA07029_AxiomGWASHuSNP1_20090812_MCKI_T01_B11_v1",
    "NA07034_AxiomGWASHuSNP1_20090812_MCKI_T01_A05_v1",
    "NA10831_AxiomGWASHuSNP1_20090812_MCKI_T01_D06_v1",
    "NA11831_AxiomGWASHuSNP1_20090812_MCKI_T01_G03_v1",
    "NA11839_AxiomGWASHuSNP1_20090812_MCKI_T01_A03_v1",
    "NA11882_AxiomGWASHuSNP1_20090812_MCKI_T01_H04_v1",
    "NA11993_AxiomGWASHuSNP1_20090812_MCKI_T01_G06_v1",
    "NA12716_AxiomGWASHuSNP1_20090812_MCKI_T01_E01_v1",
    "NA12801_AxiomGWASHuSNP1_20090812_MCKI_T01_D09_v1",
    "NA12813_AxiomGWASHuSNP1_20090812_MCKI_T01_E08_v1",
    "NA12874_AxiomGWASHuSNP1_20090812_MCKI_T01_E04_v1",
    "NA12892_AxiomGWASHuSNP1_20090812_MCKI_T01_B12_v1",
  NULL
  };
  string outdir = testDir + "/doMultiChannelAxiom";
  string name = "doMultiChannelAxiom";
  string command = "./apt-probeset-genotype "
    "-a artifact-reduction.ResType=2.Clip=0.4.Close=2.Open=2.Fringe=4.CC=2,quant-norm.target=1000.sketch=50000,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.lambda=1.0.HARD=3.SB=0.75.transform=MVA.copyqc=0.00000.wobble=0.05.MS=0.15.copytype=-1.clustertype=2.ocean=0.00001.CSepPen=0.1.CSepThr=4 "
    "--qmethod-spec med-polish.expon=true "
    "--read-models-brlmmp ../../../regression-data/data/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.models "
    "--cdf-file ../../../regression-data/data/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.cdf "
    "--special-snps ../../../regression-data/data/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.specialSNPs "
    "--chrX-probes ../../../regression-data/data/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.chrXprobes "
    "--chrY-probes ../../../regression-data/data/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.chrYprobes "
    "--target-sketch ../../../regression-data/data/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.sketch "
    "--use-feat-eff ../../../regression-data/data/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.feature-effects "
    "--set-gender-method cn-probe-chrXY-ratio "
    "--em-gender false "
    "--female-thresh 0.54 "
    "--male-thresh 1.0 "
    "--cc-chp-output "
    "--chip-type Axiom_GW_Hu_SNP "
    "--chip-type Axiom_GW_Hu_SNP.r2 "
    "--set-analysis-name AxiomGT1 "
    "--out-dir " + outdir + " ";
  command += Util::joinVectorString(Util::addPrefixSuffix(file_prefixes," ../../../regression-data/data/cel/Axiom_GW_Hu_SNP/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(file_prefixes, "../../../regression-data/data/chipstream/probeset-genotype/Axiom_GW_Hu_SNP/doMultiChannelAxiom/cc-chp/", ".AxiomGT1.chp");
  gen = Util::addPrefixSuffix(file_prefixes, testDir + "/doMultiChannelAxiom/cc-chp/", ".AxiomGT1.chp");
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001));
  checks.push_back(new MatrixCheck(testDir + "/doMultiChannelAxiom/AxiomGT1.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/Axiom_GW_Hu_SNP/doMultiChannelAxiom/AxiomGT1.calls.txt",
                                   0.0000001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck(testDir + "/doMultiChannelAxiom/AxiomGT1.confidences.txt", 
                                   "../../../regression-data/data/chipstream/probeset-genotype/Axiom_GW_Hu_SNP/doMultiChannelAxiom/AxiomGT1.confidences.txt",
                                   0.0001,
                                   1, 1, true, 0));

  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doLabelZStyHints() {
  string name = "doStyCCS_BRLMM-P_Hints";
  string outdir =  testDir + "/doLabelZStyHints";
  string command = "./apt-probeset-genotype "
    " --genotypes ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/call.hints "
    "--table-output -s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--summaries "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm-p.lowprecision=true.transform=ccs.K=4.MS=2.hints=1.CP=8 "
    "--write-models "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(testDir + "/doLabelZStyHints/quant-norm.pm-only.brlmm-p.calls.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doLabelZStyHints/quant-norm.pm-only.brlmm-p.calls.txt",
                                   0.0000001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck(testDir + "/doLabelZStyHints/quant-norm.pm-only.brlmm-p.confidences.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doLabelZStyHints/quant-norm.pm-only.brlmm-p.confidences.txt",
                                   0.0001,
                                   1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doLabelZStyTestPriors() {
  string name = "doStyCCS_BRLMM-P-priors";
  string outdir = testDir + "/doLabelZStyTestPriors";
  string command = "./apt-probeset-genotype "
    "--read-models-brlmmp ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/test.snp-prior.ref "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--summaries "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
    command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doLabelZStyTestPriors/quant-norm.pm-only.brlmm-p.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doLabelZStyTestPriors/quant-norm.pm-only.brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doLabelZStyTestPriors/quant-norm.pm-only.brlmm-p.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doLabelZStyTestPriors/quant-norm.pm-only.brlmm-p.confidences.txt",
                              0.001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doLabelZStyCCS() {
  string outdir = testDir + "/doLabelZStyCCS";
  string name = "doStyCCS_BRLMM-P-default";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--summaries "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
    command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doLabelZStyCCS/quant-norm.pm-only.brlmm-p.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doLabelZStyCCS/quant-norm.pm-only.brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doLabelZStyCCS/quant-norm.pm-only.brlmm-p.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doLabelZStyCCS/quant-norm.pm-only.brlmm-p.confidences.txt",
                              0.001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}
void ProbeSetGenotypeTest::doLabelZStyCCSseveral() {
  string name = "doStyCCS_BRLMM-P-severalparameters";
  string outdir = testDir + "/doLabelZStyCCSseveral";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--summaries "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm-p.lowprecision=true.transform=ccs.K=4.MS=2.bins=10.KX=1 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doLabelZStyCCSseveral/quant-norm.pm-only.brlmm-p.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doLabelZStyCCSseveral/quant-norm.pm-only.brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doLabelZStyCCSseveral/quant-norm.pm-only.brlmm-p.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doLabelZStyCCSseveral/quant-norm.pm-only.brlmm-p.confidences.txt",
                              0.0001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doStyCCS_2_20_1_0() {
  string name = "doStyCCS_2_20_1_0";
  string outdir = testDir + "/doStyCCS_2_20_1_0";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
    command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doStyCCS_2_40_0_9_0_033() {
  string name = "doStyCCS_2_40_0_9_0_033";
  string outdir = testDir + "/doStyCCS_2_40_0_9_0_033";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=40.lowprecision=true.het-mult=0.9.transform=ccs.K=4.MS=2 --dm-thresh 0.33 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyCCS_2_40_0_9_0_033/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyCCS_2_40_0_9_0_033/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyCCS_2_40_0_9_0_033/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyCCS_2_40_0_9_0_033/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doSpf() {
  string outdir = testDir + "/doSpf";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "--chrX-snps ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.chrx "
    "--spf-file ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o " + outdir + " "
    "--use-disk=false "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/cel-files.txt";

  string name = "doSpf";
  vector<RegressionCheck *> checks;

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck(
  //                                 testDir + "/doSpf/brlmm.report.txt",
  //                                 "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doSpf/brlmm.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(
                              testDir + "/doSpf/brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doSpf/brlmm.calls.txt",
                              0.0000001,
                              1, 1, false, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doSpf/brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doSpf/brlmm.confidences.txt",
                              0.00001,
                              1, 1, false, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doChpFiles() {
  string command = "./apt-probeset-genotype "
    "--xda-chp-output "
    "--cc-chp-output "
    "--table-output "
    "--use-disk=false "
    "--chrX-snps ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.chrx  "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + testDir + "/doChpFiles "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/cel-files.txt";

  const char *chpFiles[] = { "NA06985_B01_Sty_Plate1.brlmm","NA06991_B03_Sty_Plate1.brlmm","NA06993_B02_Sty_Plate1.brlmm",
                 "NA06994_A11_Sty_Plate1.brlmm","NA07000_A10_Sty_Plate1.brlmm","NA07019_A09_Sty_Plate1.brlmm",
                 "NA07022_A08_Sty_Plate1.brlmm","NA07029_A12_Sty_Plate1.brlmm","NA07034_B05_Sty_Plate1.brlmm",
                 "NA07048_B06_Sty_Plate1.brlmm","NA07055_B04_Sty_Plate1.brlmm","NA07056_A07_Sty_Plate1.brlmm",
                 "NA07345_B10_Sty_Plate1.brlmm","NA07348_B12_Sty_Plate1.brlmm","NA07357_B11_Sty_Plate1.brlmm",
                 "NA10846_A06_Sty_Plate1.brlmm","NA10847_A03_Sty_Plate1.brlmm","NA10851_B09_Sty_Plate1.brlmm",
                 "NA10854_C09_Sty_Plate1.brlmm","NA10855_C12_Sty_Plate1.brlmm",
                 NULL
  };

  string name = "doChpFiles";
  vector<RegressionCheck *> checks;

  vector<string> gold,gen,goldX,genX;
  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doChpFiles/chp/",".chp");
  gen = Util::addPrefixSuffix(chpFiles, testDir + "/doChpFiles/chp/", ".chp");
  goldX = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doChpFiles/cc-chp/", ".chp");
  genX = Util::addPrefixSuffix(chpFiles, testDir + "/doChpFiles/cc-chp/",".chp");


  // checks for chp files.
  checks.push_back(new ChpCheck(gen,gold));
  checks.push_back(new CalvinChpCheck(genX,goldX));

  // checks for report file
  checks.push_back(new MixedFileCheck(
                                   testDir + "/doChpFiles/brlmm.report.txt",
                                   "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doChpFiles/brlmm.report.txt",
                                   0.00001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(
                              testDir + "/doChpFiles/brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doChpFiles/brlmm.calls.txt",
                              0.0000001,
                              1, 1, false, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doChpFiles/brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doChpFiles/brlmm.confidences.txt",
                              0.00001,
                              1, 1, false, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}


void ProbeSetGenotypeTest::doFileStyCCS_2_20_1_0() {
  string outdir = testDir + "/doFileStyCCS_2_20_1_0";
  string name = "doFileStyCCS_2_20_1_0";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/cel-files.txt";

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}


void ProbeSetGenotypeTest::doWritePriorFileStyCCS_2_20_1_0() {
  string outdir = testDir + "/doWritePriorFileStyCCS_2_20_1_0";
  string name = "doWritePriorFileStyCCS_2_20_1_0";
  // shuld get the same answer if using prior generated on same run.
  string command ="./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/cel-files.txt "
    "--write-prior";

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "() Phase 1" );
  if(brlmmTest.pass()) {
    string outdir2 = testDir + "/doWritePriorFileStyCCS_2_20_1_0/phase2";
    string command2 = "./apt-probeset-genotype "
      "--table-output "
      "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
      "--list-sample "
      "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
      "--no-gender-force "
      "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
      "-o " + outdir2 + " "
      "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/cel-files.txt "
      "--read-priors-brlmm " + testDir + "/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.prior.txt";

    vector<RegressionCheck *> checks2;
    checks2.push_back(new MatrixCheck(testDir + "/doWritePriorFileStyCCS_2_20_1_0/phase2/quant-norm.pm-only.brlmm.calls.txt",
                                testDir + "/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                                0.0000001,
                                1, 1, true, 0));
    ///@todo relaxed test from 10-5 to 10-4 due to failures. probably due to changes in rounding of reports
    checks2.push_back(new MatrixCheck(testDir + "/doWritePriorFileStyCCS_2_20_1_0/phase2/quant-norm.pm-only.brlmm.confidences.txt",
                                testDir + "/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                                0.0001,
                                1, 1, true, 0));
    string name2 = "doWritePriorFileStyCCS_2_20_1_0-phase2";
    RegressionTest brlmmTest2(name2.c_str(), command2.c_str(), checks2);
    brlmmTest2.setSuite(*this, outdir2, outdir2 + "/apt-probeset-genotype.log", outdir2 + "/valgrind.log");
    Verbose::out(1, "Doing " + name2 + "() Phase 2" );
    if(!brlmmTest2.pass()) {
      Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name2 + "(): " + brlmmTest.getErrorMsg());
      numFailed++;
    }
    else {
      numPassed++;
    }
  }
  else {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
}

void ProbeSetGenotypeTest::doReadGenotypesIn() {
  string outdir = testDir + "/doReadGenotypesIn";
  string name = "doReadGenotypesIn";
  // should get the same answer if using prior generated on same run.
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/cel-files.txt "
    "--dm-out "
    "--dm-thresh .33";

  vector<RegressionCheck *> checks;

  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "() Phase 1" );
  if(brlmmTest.pass()) {
    string outdir2 = testDir + "/doReadGenotypesIn/phase2";
    string name2 = "doReadGenotypesIn-phase2";
    string command2 = "./apt-probeset-genotype "
      "--table-output "
      "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
      "--list-sample "
      "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
      "--no-gender-force "
      "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
      "-o " + outdir2 + " "
      "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/cel-files.txt "
      "--genotypes " + testDir + "/doReadGenotypesIn/dm.calls.txt";

    vector<RegressionCheck *> checks2;
    checks2.push_back(new MatrixCheck(
                                testDir + "/doReadGenotypesIn/phase2/quant-norm.pm-only.brlmm.calls.txt",
                                testDir + "/doReadGenotypesIn/quant-norm.pm-only.brlmm.calls.txt",
                                0.0000001,
                                1, 1, true, 0));
    checks2.push_back(new MatrixCheck(
                                testDir + "/doReadGenotypesIn/phase2/quant-norm.pm-only.brlmm.confidences.txt",
                                testDir + "/doReadGenotypesIn/quant-norm.pm-only.brlmm.confidences.txt",
                                0.00001,
                                1, 1, true, 0));
    RegressionTest brlmmTest2(name2.c_str(), command2.c_str(), checks2);
    brlmmTest2.setSuite(*this, outdir2, outdir2 + "/apt-probeset-genotype.log", outdir2 + "/valgrind.log");
    Verbose::out(1, "Doing " + name + "() Phase 2" );
    if(!brlmmTest2.pass()) {
      Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
      numFailed++;
    }
    else {
      numPassed++;
    }
  }
  else {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
}

void ProbeSetGenotypeTest::doStyRVT_2_20_1_0() {
  string name = "doStyRvT_2_20_1_0";
  string outdir = testDir + "/doStyRVT_2_20_1_0";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample -a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=rvt.MS=2 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " " ;
    command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyRVT_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyRVT_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyRVT_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyRVT_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doStyCES_2_20_1_0() {
  string outdir = testDir + "/doStyCES_2_20_1_0";
  string name = "doStyCES_2_20_1_0";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ces.K=1.MS=2 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyCES_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyCES_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyCES_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyCES_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeTest::doStyMva_2_20_1_0() {
  string name = "doStyMva_2_20_1_0";
  string outdir = testDir + "/doStyMva_2_20_1_0";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=1.MS=2 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyMva_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyMva_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyMva_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyMva_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}


void ProbeSetGenotypeTest::doStyMva_2_20_1_2() {
  string outdir = testDir + "/doStyMva_2_20_1_2";
  string name = "doStyMva_2_20_1_2";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=1.iterations=1.iter-thresh=.6.MS=2 --no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyMva_2_20_1_2/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyMva_2_20_1_2/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyMva_2_20_1_2/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyMva_2_20_1_2/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}


void ProbeSetGenotypeTest::doStyMva_2_20_0_8_0() {
  string name = "doStyMva_2_20_0_8_0";
  string outdir = testDir + "/doStyMva_2_20_0_8_0";
  string command = "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/lib/Mapping250K_Sty/set1snps.txt "
    "--list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=.8.MS=2 "
    "--no-gender-force "
    "-c ../../../regression-data/data/lib/Mapping250K_Sty/Mapping250K_Sty.cdf "
    "-o " + outdir + " ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyMva_2_20_0_8_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyMva_2_20_0_8_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              testDir + "/doStyMva_2_20_0_8_0/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/chipstream/probeset-genotype/Mapping250K_Sty/doStyMva_2_20_0_8_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  brlmmTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + brlmmTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}


void ProbeSetGenotypeTest::doTrustedProbes() {
  string name = "doTrustedProbes";
  string outdir = testDir + "/doTrustedProbes";
  string command = "./apt-probeset-genotype "
    "--analysis-files-path=../../../regression-data/data/lib/AxiomLib_Release "
    "--xml-file=../../../regression-data/data/chipstream/probeset-genotype/Axiom_Release/doTrustedProbes/trusted.probes.regression.xml "
    "--read-models-brlmmp ../../../regression-data/data/lib/AxiomLib_Release/AxiomGeneric.models "
    "--out-dir " + testDir + "/doTrustedProbes "
    "--summaries --write-models "
    "--cel-files ../../../regression-data/data/chipstream/probeset-genotype/Axiom_Release/doTrustedProbes/CELregression.list.txt "
    //    "&& perl -i -p -e 's{ ,  }{\t}gxms' " + testDir + "/doTrustedProbes/AxiomGT1.snp-posteriors.txt "
    ;

  std::string golddir = "../../../regression-data/data/chipstream/probeset-genotype/Axiom_Release/doTrustedProbes";

  const TrustedProbesCheck trusted_reports[] = {
    {"AxiomGT1.calls.txt", false, false, 74, 1, 0.0001 },
    {"AxiomGT1.confidences.txt", false, true, 71, 1, 0.0001 },
    {"AxiomGT1.normalized-summary.txt", false, true, 3, 2, 0.0001 },
    {"AxiomGT1.report.txt", false, false, 73, 1, 0.0001 },
    {"AxiomGT1.snp-posteriors.txt", true, true, 73, 1, 0.0003 },
    {"AxiomGT1.summary.txt", false, true,  72, 2, 0.06 },
  };
  
  int trusted_reports_size = sizeof( trusted_reports ) / (sizeof(trusted_reports[0]));

  vector<RegressionCheck *> checks;
  for ( int i =0; i < trusted_reports_size; i++ ) {
    if ( trusted_reports[i].posteriorMatrixCheck  ) {
      checks.push_back(new PosteriorMatrixCheck( outdir + "/" + trusted_reports[i].name, 
                                        golddir + "/" + trusted_reports[i].name,
                                        trusted_reports[i].eps,
                                        trusted_reports[i].skiplines,
                                        trusted_reports[i].skipcols,
                                        false, 0 ) );

    }
    else if ( trusted_reports[i].matrixCheck  ) {
      checks.push_back(new MatrixCheck( outdir + "/" + trusted_reports[i].name, 
                                        golddir + "/" + trusted_reports[i].name,
                                        trusted_reports[i].eps,
                                        trusted_reports[i].skiplines,
                                        trusted_reports[i].skipcols,
                                        false, 0 ) );

    }
    else {
      checks.push_back( new TextFileCheck( outdir + "/" + trusted_reports[i].name, 
                                           golddir + "/" + trusted_reports[i].name,
                                           std::string("#")));

    }
  }
  
  RegressionTest trustedProbesTest(name.c_str(), command.c_str(), checks);
  trustedProbesTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");
  Verbose::out(1, "Doing " + name + "()");
  if(!trustedProbesTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + trustedProbesTest.getErrorMsg());
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
    testDir.setTestDir("chipstream/genotype", true);
    
    ProbeSetGenotypeTest test;
    test.testDir = testDir.asString();
    
    test.parseArgv(argv);
    bool doValgrind = test.doValgrind();
    string database = test.getDatabase();
    cout << "Valgrind: " + ToStr(doValgrind) + " database: " + database << endl;
    int testbirdseedflag, testbrlmmflag,testbrlmmpflag;

    testbirdseedflag=testbrlmmflag=testbrlmmpflag = 0;

    Verbose::setLevel(2);

    if (test.leftoverArgs() < 1)
    {
        testbrlmmflag = 1;
        testbrlmmpflag = 1;
        testbirdseedflag=1;
    }
    else
    {
      if (!strcmp(argv[1],"brlmm"))
          testbrlmmflag=1;
      if (!strcmp(argv[1],"brlmmp"))
          testbrlmmpflag=1;
      if (!strcmp(argv[1],"birdseed"))
          testbirdseedflag=1;
    }
    if (testbrlmmflag)
    {
        test.doFileStyCCS_2_20_1_0();
        test.doStyCCS_2_20_1_0();
        test.doChpFiles();
        test.doStyCCS_2_40_0_9_0_033();
        test.doReadGenotypesIn();
        test.doWritePriorFileStyCCS_2_20_1_0();
        test.doStyMva_2_20_1_2();
        test.doStyCES_2_20_1_0();
        test.doStyRVT_2_20_1_0();
        test.doStyMva_2_20_1_0();
        test.doStyMva_2_20_0_8_0();

        test.doSpf();
    }
    // labelz = brlmm-p
    if (testbrlmmpflag)
    {
        #ifdef __linux__
        #ifdef __LP64__
        test.doBrlmmpSnp5CdfReadWriteFeatureEffectsA5();
        #endif
        #endif
        test.doTrustedProbes();
        test.doLabelZStyCCS();
        test.doLabelZStyTestPriors();
        test.doLabelZStyCCSseveral();
        test.doLabelZStyHints();

        test.doBrlmmpSnp5Chp();
        test.doBrlmmpSnp5Spf();
        test.doBrlmmpSnp5KillList();
        ///@todo broken test Martin is workin on fix
        ///@todo He only did this hack for now.
        test.doMultiChannelAxiom();
    }

    if (testbirdseedflag)
    {
        test.doBirdseedSnp6Chp();
        test.doBrlmmpSnp6Chp();
        test.doBirdseedSnp6Spf();
        test.doBirdseedSnp6Full();
        test.doBirdseed2Snp6Chp();
    }

        
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

