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

#include "util/CalvinChpCheck.h"
#include "util/ChpCheck.h"
#include "util/Convert.h"
#include "util/Fs.h"
#include "util/RegressionTest.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>
//

using namespace std;

string valgrindString1 = "valgrind --leak-check=yes --log-file=./test-generated/valgrind/";
string valgrindString2 = ".valgrind.txt ";

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

class ProbeSetGenotypeTest {

public:
  int numPassed, numFailed;
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
  string name = "qt-doBirdseedSnp6Spf";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "--analysis birdseed "
    "--special-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.SNP_A-4232288.spf "
    "--read-models-birdseed ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o test-generated/doBirdseedSnp6Spf "
    "--cel-files ../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/fas_cel_files.txt";

  vector<RegressionCheck *> checks;

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck("test-generated/doBirdseedSnp6Spf/birdseed.report.txt",
  //                                 "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck("test-generated/doBirdseedSnp6Spf/birdseed.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck("test-generated/doBirdseedSnp6Spf/birdseed.confidences.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.confidences.txt", 
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
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
  string name = "qt-doBirdseedSnp6Full";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--cc-chp-output "
    "--table-output "
    "--analysis birdseed "
    "--special-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.Full.specialSNPs "
    "--cdf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.Full.cdf "
    "--read-models-birdseed ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o test-generated/doBirdseedSnp6Full "
    "--cel-files ../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/fas_cel_files.txt";

  const char *chpFiles[] = {"NA06985_GW6_C.birdseed", "NA06993_GW6_C.birdseed", "NA07000_GW6_C.birdseed", 
			    "NA07019_GW6_C.birdseed", "NA07029_GW6_C.birdseed", "NA07056_GW6_C.birdseed", 
			    "NA07345_GW6_C.birdseed", "NA10830_GW6_C.birdseed", "NA10831_GW6_C.birdseed", 
			    "NA10839_GW6_C.birdseed", "NA10855_GW6_C.birdseed", "NA10857_GW6_C.birdseed", 
			    "NA10860_GW6_C.birdseed", "NA10863_GW6_C.birdseed", "NA11881_GW6_C.birdseed", 
			    "NA11882_GW6_C.birdseed", "NA11993_GW6_C.birdseed", "NA12003_GW6_C.birdseed", 
			    "NA12004_GW6_C.birdseed",
			    NULL
};

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Full/cc-chp/", ".chp");
  gen = Util::addPrefixSuffix(chpFiles, "test-generated/doBirdseedSnp6Full/cc-chp/", ".chp");

  // checks for chp files.
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001));

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck("test-generated/doBirdseedSnp6Full/birdseed.report.txt",
  //                                 "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Full/birdseed.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck("test-generated/doBirdseedSnp6Full/birdseed.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Full/birdseed.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck("test-generated/doBirdseedSnp6Full/birdseed.confidences.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Full/birdseed.confidences.txt", 
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
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
  string name = "qt-doBirdseedSnp6Chp";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--cc-chp-output "
    "--table-output "
    "--analysis birdseed "
    "--special-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.SNP_A-4232288.spf "
    "--read-models-birdseed ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o test-generated/doBirdseedSnp6Chp "
    "--cel-files ../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/fas_cel_files.txt";

  const char *chpFiles[] = {"NA06985_GW6_C.birdseed", "NA06993_GW6_C.birdseed", "NA07000_GW6_C.birdseed", 
			    "NA07019_GW6_C.birdseed", "NA07029_GW6_C.birdseed", "NA07056_GW6_C.birdseed", 
			    "NA07345_GW6_C.birdseed", "NA10830_GW6_C.birdseed", "NA10831_GW6_C.birdseed", 
			    "NA10839_GW6_C.birdseed", "NA10855_GW6_C.birdseed", "NA10857_GW6_C.birdseed", 
			    "NA10860_GW6_C.birdseed", "NA10863_GW6_C.birdseed", "NA11881_GW6_C.birdseed", 
			    "NA11882_GW6_C.birdseed", "NA11993_GW6_C.birdseed", "NA12003_GW6_C.birdseed", 
			    "NA12004_GW6_C.birdseed",
			    NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Chp/cc-chp/", ".chp");
  gen = Util::addPrefixSuffix(chpFiles, "test-generated/doBirdseedSnp6Chp/cc-chp/", ".chp");

  // checks for chp files.
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001));
  // checks for report file
  checks.push_back(new MixedFileCheck("test-generated/doBirdseedSnp6Chp/birdseed.report.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Chp/birdseed.report.txt",
                                   0.0001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck("test-generated/doBirdseedSnp6Chp/birdseed.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Chp/birdseed.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck("test-generated/doBirdseedSnp6Chp/birdseed.confidences.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Chp/birdseed.confidences.txt", 
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
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
  string name = "qt-doBirdseed2Snp6Chp";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--cc-chp-output "
    "--table-output "
    "--analysis birdseed-v2 "
    "--special-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.spf "
    "--read-models-birdseed ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o test-generated/doBirdseed2Snp6Chp "
    "--cel-files ../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/fas_cel_files.txt "
    "--set-gender-method cn-probe-chrXY-ratio "
    "--chrX-probes ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes "
    "--chrY-probes ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes ";

  const char *chpFiles[] = {"NA06985_GW6_C", "NA06993_GW6_C", "NA07000_GW6_C", 
			    "NA07019_GW6_C", "NA07029_GW6_C", "NA07056_GW6_C", 
			    "NA07345_GW6_C", "NA10830_GW6_C", "NA10831_GW6_C", 
			    "NA10839_GW6_C", "NA10855_GW6_C", "NA10857_GW6_C", 
			    "NA10860_GW6_C", "NA10863_GW6_C", "NA11881_GW6_C", 
			    "NA11882_GW6_C", "NA11993_GW6_C", "NA12003_GW6_C", 
			    "NA12004_GW6_C",
			    NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseed2Snp6Chp/cc-chp/", ".birdseed-v2.chp");
  gen = Util::addPrefixSuffix(chpFiles, "test-generated/doBirdseed2Snp6Chp/cc-chp/", ".birdseed-v2.chp");

  
  // checks for chp files.
  ///@todo we need to relax Eps due to insignificant difference in extra (non-call, non-confidence) values for one SNP between win32 and other platforms.
  ///      the interface does not (yet) allow us to set different eps for each field -- which would be ideal. With this change we are not
  ///      doing as tight a check on confidences as we probably should. Of course the test on the matrix file mitigates some risk here.
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001));
  // checks for report file
  checks.push_back(new MixedFileCheck("test-generated/doBirdseed2Snp6Chp/birdseed-v2.report.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseed2Snp6Chp/birdseed-v2.report.txt",
                                   0.0001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck("test-generated/doBirdseed2Snp6Chp/birdseed-v2.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseed2Snp6Chp/birdseed-v2.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck("test-generated/doBirdseed2Snp6Chp/birdseed-v2.confidences.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseed2Snp6Chp/birdseed-v2.confidences.txt", 
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
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
  string name = "qt-doBrlmmpSnp6Chp";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--cc-chp-output "
    "--table-output "
    "--analysis brlmm-p "
    "--special-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.spf "
    "--read-models-brlmmp  ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.brlmm-p.models  "
    "-o test-generated/doBrlmmpSnp6Chp "
    "--cel-files ../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/fas_cel_files.txt";
  
  const char *chpFiles[] = {"NA06985_GW6_C.brlmm-p", "NA06993_GW6_C.brlmm-p", "NA07000_GW6_C.brlmm-p", 
			    "NA07019_GW6_C.brlmm-p", "NA07029_GW6_C.brlmm-p", "NA07056_GW6_C.brlmm-p", 
			    "NA07345_GW6_C.brlmm-p", "NA10830_GW6_C.brlmm-p", "NA10831_GW6_C.brlmm-p", 
			    "NA10839_GW6_C.brlmm-p", "NA10855_GW6_C.brlmm-p", "NA10857_GW6_C.brlmm-p", 
			    "NA10860_GW6_C.brlmm-p", "NA10863_GW6_C.brlmm-p", "NA11881_GW6_C.brlmm-p", 
			    "NA11882_GW6_C.brlmm-p", "NA11993_GW6_C.brlmm-p", "NA12003_GW6_C.brlmm-p", 
			    "NA12004_GW6_C.brlmm-p", 
			    NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBrlmmpSnp6Chp/cc-chp/", ".chp");
  gen = Util::addPrefixSuffix(chpFiles, "test-generated/doBrlmmpSnp6Chp/cc-chp/", ".chp");

  // checks for chp files.
  checks.push_back(new CalvinChpCheck(gen,gold, 0, L"apt-", 0.001,false));
  // checks for report file
  ///@todo relaxed due to difference in chr-x-het-rate metric on win32
  checks.push_back(new MixedFileCheck("test-generated/doBrlmmpSnp6Chp/brlmm-p.report.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBrlmmpSnp6Chp/brlmm-p.report.txt",
                                   0.0001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp6Chp/brlmm-p.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBrlmmpSnp6Chp/brlmm-p.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks.push_back(new MatrixCheck("../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBrlmmpSnp6Chp/brlmm-p.confidences.txt", 
                                   "test-generated/doBrlmmpSnp6Chp/brlmm-p.confidences.txt",
                                   0.0001,
                                   1, 1, false, 0));
  RegressionTest brlmmpTest(name.c_str(), command.c_str(), checks);
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
  string name = "qt-doBrlmmpSnp5Spf";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.spf "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir test-generated/doBrlmmpSnp5Spf";
          
  const char *chpFiles[] = {
    "NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C","NA06994_GW5_C","NA07000_GW5_C",
    "NA07019_GW5_C","NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C","NA07048_GW5_C",
    "NA07055_GW5_C","NA07056_GW5_C","NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
    "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C","NA10838_GW5_C","NA10839_GW5_C",
    NULL
  };

  vector<RegressionCheck *> checks;


  command += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/idata/cel/GenomeWideSNP_5/", ".CEL"), " ");

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5Spf/brlmm-p.report.txt",
  //                                 "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Spf/brlmm-p.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5Spf/brlmm-p.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Spf/brlmm-p.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5Spf/brlmm-p.confidences.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Spf/brlmm-p.confidences.txt", 
                                   0.0001,
                                   1, 1, false, 0));
  RegressionTest brlmmpTest(name.c_str(), command.c_str(), checks);
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
  string name1 = "qt-doBrlmmpSnp5CdfReadWriteFeatureEffectsA5-part1";  

  string commandMakeFeatureEffects = valgrindString1 + name1 + valgrindString2 + "./apt-probeset-genotype "
    "--qmethod-spec med-polish "
    "--analysis brlmm-p "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.spf "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--out-dir test-generated/doBrlmmpSnp5CdfWriteFeatureEffectsA5 "
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
  
  commandMakeFeatureEffects += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/idata/cel/GenomeWideSNP_5/", ".CEL"), " ");

  //First Test  Feature effects are created and written into A5 format. 
  vector<RegressionCheck *> makeChecks;

  makeChecks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5CdfWriteFeatureEffectsA5/brlmm-p.summary.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.summary.txt",
                                   0.00001,
                                   1, 1, false, 0));

  RegressionTest makeFeatureEffectsTest(name1.c_str(), commandMakeFeatureEffects.c_str(), makeChecks);


  //Second Test  The feature effects just created are read in and used to create summary values.  
  string name2 = "qt-doBrlmmpSnp5CdfReadWriteFeatureEffectsA5-part2"; 

  string commandUseFeatureEffects = valgrindString1 + name2 + valgrindString2 + " ./apt-probeset-genotype "  
    "--qmethod-spec med-polish "
    "--analysis brlmm-p "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.spf "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--a5-feature-effects-input-file test-generated/doBrlmmpSnp5CdfWriteFeatureEffectsA5/brlmm-p.feature-response.a5 "
    "--a5-feature-effects-input-name  brlmm-p.feature-response "
    "--out-dir test-generated/doBrlmmpSnp5CdfReadFeatureEffectsA5 "
    "--summaries";
  
  commandUseFeatureEffects += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/idata/cel/GenomeWideSNP_5/", ".CEL"), " ");
  


  vector<RegressionCheck *> useChecks;

  useChecks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5CdfReadFeatureEffectsA5/brlmm-p.summary.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.summary.txt",
                                   0.00002,
                                   1, 1, false, 0));
 

  RegressionTest useFeatureEffectsTest(name2.c_str(), commandUseFeatureEffects.c_str(), useChecks);

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
  string name = "qt-doBrlmmpSnp5Chp";
    
  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.spf "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir test-generated/doBrlmmpSnp5Chp "
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

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Chp/cc-chp/", ".brlmm-p.chp");
  gen = Util::addPrefixSuffix(chpFiles,"test-generated/doBrlmmpSnp5Chp/cc-chp/", ".brlmm-p.chp");
  command += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/idata/cel/GenomeWideSNP_5/", ".CEL"), " ");

  // checks for chp files. [ignore header because of command string]
  checks.push_back(new CalvinChpCheck(gen,gold, 0, L"apt-", 0.001,false));
  // checks for report file
  checks.push_back(new MixedFileCheck("test-generated/doBrlmmpSnp5Chp/brlmm-p.report.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.report.txt",
                                   0.00001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5Chp/brlmm-p.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.calls.txt", 
                                   0.0000001,
                                   1, 1, false, 0));
  // checks for the table files.
  ///@todo Is this too loose?
  checks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5Chp/brlmm-p.normalized-summary.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.normalized-summary.txt", 
                                   0.001,
                                   1, 2, false, 0));
  ///@todo Is this too loose?
  checks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5Chp/brlmm-p.summary.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.summary.txt", 
                                   0.001,
                                   1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks.push_back(new MatrixCheck("test-generated/doBrlmmpSnp5Chp/brlmm-p.confidences.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_5/doBrlmmpSnp5Chp/brlmm-p.confidences.txt", 
                                   0.0001,
                                   1, 1, false, 0));
  RegressionTest brlmmpTest(name.c_str(), command.c_str(), checks);
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
  string name = "qt-doBrlmmpSnp5KillList";
 
  string command1 = valgrindString1 + name + "-part1"+ valgrindString2 + "./apt-probeset-genotype "   
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.spf "
    "--kill-list ../../../regression-data/data/idata/lib/GenomeWideSNP_5/test-kill-list.txt "
    "-s ../../../regression-data/data/idata/lib/GenomeWideSNP_5/test-kill-list.snp-list.txt "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir test-generated/doBrlmmpSnp5KillList "
    "--summaries ";

  string command2 = valgrindString1 + name + "-part2" + valgrindString2 + "./apt-probeset-genotype "   
    "--force "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/test-kill-list.spf "
    "-s ../../../regression-data/data/idata/lib/GenomeWideSNP_5/test-kill-list.snp-list.txt "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir test-generated/doBrlmmpSnp5KillList/mask2 "
    "--summaries ";

  string command3 = valgrindString1 + name + "-part3" + valgrindString2 + "./apt-probeset-genotype "   
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.spf "
    "-s ../../../regression-data/data/idata/lib/GenomeWideSNP_5/test-kill-list.snp-list.txt "
    "--chrX-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.chrx "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.models "
    "--analysis brlmm-p "
    "--out-dir test-generated/doBrlmmpSnp5KillList/mask3 "
    "--summaries ";

  //const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C","NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C","NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C","NA07048_GW5_C","NA07055_GW5_C","NA07056_GW5_C","NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C","NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C","NA10838_GW5_C","NA10839_GW5_C"};
  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
			    "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
			    NULL
  };

  vector<RegressionCheck *> checks1;
  vector<RegressionCheck *> checks2;
  vector<RegressionCheck *> checks3;

  command1 += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/idata/cel/GenomeWideSNP_5/", ".CEL"), " ");
  command2 += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/idata/cel/GenomeWideSNP_5/", ".CEL"), " ");
  command3 += Util::joinVectorString(Util::addPrefixSuffix(chpFiles, " ../../../regression-data/data/idata/cel/GenomeWideSNP_5/", ".CEL"), " ");

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  //checks2.push_back(new MatrixCheck(
  //                            "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.report.txt",
  //                            "test-generated/doBrlmmpSnp5KillList/brlmm-p.report.txt",
  //                            0.00001,
  //                            1,2, false, 0));
  checks2.push_back(new MatrixCheck(
                              "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.calls.txt",
                              "test-generated/doBrlmmpSnp5KillList/brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, false, 0));
  ///@todo Is this too loose?
  checks2.push_back(new MatrixCheck(
                              "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.normalized-summary.txt",
                              "test-generated/doBrlmmpSnp5KillList/brlmm-p.normalized-summary.txt",
                              0.001,
                              1, 2, false, 0));
  ///@todo Is this too loose?
  checks2.push_back(new MatrixCheck(
                              "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.summary.txt",
                              "test-generated/doBrlmmpSnp5KillList/brlmm-p.summary.txt",
                              0.001,
                              1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks2.push_back(new MatrixCheck(
                              "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.confidences.txt",
                              "test-generated/doBrlmmpSnp5KillList/brlmm-p.confidences.txt",
                              0.0001,
                              1, 1, false, 0));

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  //checks3.push_back(new MatrixCheck(
  //                            "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.report.txt",
  //                            "test-generated/doBrlmmpSnp5KillList/brlmm-p.report.txt",
  //                            0.00001,
  //                            1,2, false, 0));
  checks3.push_back(new MatrixCheck(
                              "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.calls.txt",
                              "test-generated/doBrlmmpSnp5KillList/brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, false, 0));
  ///@todo Is this too loose?
  checks3.push_back(new MatrixCheck(
                              "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.normalized-summary.txt",
                              "test-generated/doBrlmmpSnp5KillList/brlmm-p.normalized-summary.txt",
                              0.001,
                              1, 2, false, 0));
  ///@todo Is this too loose?
  checks3.push_back(new MatrixCheck(
                              "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.summary.txt",
                              "test-generated/doBrlmmpSnp5KillList/brlmm-p.summary.txt",
                              0.001,
                              1, 1, false, 0));
  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  checks3.push_back(new MatrixCheck(
                              "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.confidences.txt",
                              "test-generated/doBrlmmpSnp5KillList/brlmm-p.confidences.txt",
                              0.0001,
                              1, 1, false, 0));

  string subName;
  subName = name + "-part1";
  RegressionTest test1(subName.c_str(), command1.c_str(), checks1);
  subName = name + "-part2";
  RegressionTest test2(subName.c_str(), command2.c_str(), checks2);
  subName = name + "-part3";
  RegressionTest test3(subName.c_str(), command3.c_str(), checks3);

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
    numPassed++;
  }
  else {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): Phase 2: " + test2.getErrorMsg());
    numFailed++;
  }

}

void ProbeSetGenotypeTest::doLabelZStyHints() {
  string name = "qt-doStyCCS_BRLMM-P_Hints";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    " --genotypes  ../../../regression-data/data/idata/lib/Mapping250K_Sty/call.hints "
    "--table-output -s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt " 
    "--summaries " 
    "--prior-size 1000 --list-sample " 
    "-a quant-norm.sketch=50000,pm-only,brlmm-p.lowprecision=true.transform=ccs.K=4.MS=2.hints=1.CP=8 " 
    "--write-models " 
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doLabelZStyHints ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck("test-generated/doLabelZStyHints/quant-norm.pm-only.brlmm-p.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doLabelZStyHints/quant-norm.pm-only.brlmm-p.calls.txt",
                                   0.0000001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck("test-generated/doLabelZStyHints/quant-norm.pm-only.brlmm-p.confidences.txt", 
                                   "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doLabelZStyHints/quant-norm.pm-only.brlmm-p.confidences.txt",
                                   0.0001,
                                   1, 1, true, 0));
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

void ProbeSetGenotypeTest::doLabelZStyTestPriors() {
  string name = "qt-doStyCCS_BRLMM-P-priors";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--read-models-brlmmp ../../../regression-data/data/idata/lib/Mapping250K_Sty/test.snp-prior.ref "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--summaries "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doLabelZStyTestPriors ";
    command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doLabelZStyTestPriors/quant-norm.pm-only.brlmm-p.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doLabelZStyTestPriors/quant-norm.pm-only.brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doLabelZStyTestPriors/quant-norm.pm-only.brlmm-p.confidences.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doLabelZStyTestPriors/quant-norm.pm-only.brlmm-p.confidences.txt",
                              0.001,
                              1, 1, true, 0));
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

void ProbeSetGenotypeTest::doLabelZStyCCS() {
  string name = "qt-doStyCCS_BRLMM-P-default";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--summaries "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doLabelZStyCCS ";
    command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doLabelZStyCCS/quant-norm.pm-only.brlmm-p.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doLabelZStyCCS/quant-norm.pm-only.brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doLabelZStyCCS/quant-norm.pm-only.brlmm-p.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doLabelZStyCCS/quant-norm.pm-only.brlmm-p.confidences.txt",
                              0.001,
                              1, 1, true, 0));
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
void ProbeSetGenotypeTest::doLabelZStyCCSseveral() {
  string name = "qt-doStyCCS_BRLMM-P-severalparameters";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--summaries "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm-p.lowprecision=true.transform=ccs.K=4.MS=2.bins=10.KX=1 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doLabelZStyCCSserveral ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doLabelZStyCCSserveral/quant-norm.pm-only.brlmm-p.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doLabelZStyCCSserveral/quant-norm.pm-only.brlmm-p.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doLabelZStyCCSserveral/quant-norm.pm-only.brlmm-p.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doLabelZStyCCSserveral/quant-norm.pm-only.brlmm-p.confidences.txt",
                              0.0001,
                              1, 1, true, 0));
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

void ProbeSetGenotypeTest::doStyCCS_2_20_1_0() {
  string name = "qt-doStyCCS_2_20_1_0";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doStyCCS_2_20_1_0 ";
    command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
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

void ProbeSetGenotypeTest::doStyCCS_2_40_0_9_0_033() {
  string name = "qt-doStyCCS_2_40_0_9_0_033";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=40.lowprecision=true.het-mult=0.9.transform=ccs.K=4.MS=2 --dm-thresh 0.33 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doStyCCS_2_40_0_9_0_033 "; 
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels," ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyCCS_2_40_0_9_0_033/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyCCS_2_40_0_9_0_033/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyCCS_2_40_0_9_0_033/quant-norm.pm-only.brlmm.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyCCS_2_40_0_9_0_033/quant-norm.pm-only.brlmm.confidences.txt", 
                              0.00001,
                              1, 1, true, 0));
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

void ProbeSetGenotypeTest::doSpf() {
  string name = "qt-doSpf";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "--chrX-snps ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.chrx "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doSpf "
    "--cel-files ../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt";

  vector<RegressionCheck *> checks;

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck(
  //                                 "test-generated/doSpf/brlmm.report.txt",
  //                                 "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doSpf/brlmm.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(
                              "test-generated/doSpf/brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doSpf/brlmm.calls.txt",
                              0.0000001,
                              1, 1, false, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doSpf/brlmm.confidences.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doSpf/brlmm.confidences.txt",
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

void ProbeSetGenotypeTest::doChpFiles() {
  string name = "qt-doChpFiles";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--xda-chp-output "
    "--cc-chp-output "
    "--table-output "
    "--chrX-snps ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.chrx  "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doChpFiles "
    "--cel-files ../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt";

  const char *chpFiles[] = { "NA06985_B01_Sty_Plate1.brlmm","NA06991_B03_Sty_Plate1.brlmm","NA06993_B02_Sty_Plate1.brlmm",
			     "NA06994_A11_Sty_Plate1.brlmm","NA07000_A10_Sty_Plate1.brlmm","NA07019_A09_Sty_Plate1.brlmm",
			     "NA07022_A08_Sty_Plate1.brlmm","NA07029_A12_Sty_Plate1.brlmm","NA07034_B05_Sty_Plate1.brlmm",
			     "NA07048_B06_Sty_Plate1.brlmm","NA07055_B04_Sty_Plate1.brlmm","NA07056_A07_Sty_Plate1.brlmm",
			     "NA07345_B10_Sty_Plate1.brlmm","NA07348_B12_Sty_Plate1.brlmm","NA07357_B11_Sty_Plate1.brlmm",
			     "NA10846_A06_Sty_Plate1.brlmm","NA10847_A03_Sty_Plate1.brlmm","NA10851_B09_Sty_Plate1.brlmm",
			     "NA10854_C09_Sty_Plate1.brlmm","NA10855_C12_Sty_Plate1.brlmm",
			     NULL
  };

  vector<RegressionCheck *> checks;
  
  vector<string> gold,gen,goldX,genX;
  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doChpFiles/chp/",".chp");
  gen = Util::addPrefixSuffix(chpFiles, "test-generated/doChpFiles/chp/", ".chp");
  goldX = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doChpFiles/cc-chp/", ".chp");
  genX = Util::addPrefixSuffix(chpFiles, "test-generated/doChpFiles/cc-chp/",".chp");

  
  // checks for chp files.
  checks.push_back(new ChpCheck(gen,gold));
  checks.push_back(new CalvinChpCheck(genX,goldX));

  // checks for report file
  checks.push_back(new MixedFileCheck(
                                   "test-generated/doChpFiles/brlmm.report.txt",
                                   "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doChpFiles/brlmm.report.txt",
                                   0.00001,
                                   0,0));
  // checks for the table files.
  checks.push_back(new MatrixCheck( 
                              "test-generated/doChpFiles/brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doChpFiles/brlmm.calls.txt",
                              0.0000001,
                              1, 1, false, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doChpFiles/brlmm.confidences.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doChpFiles/brlmm.confidences.txt",
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
  string name = "qt-doFileStyCCS_2_20_1_0";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doFileStyCCS_2_20_1_0 "
    "--cel-files ../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt";

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              0.00001,
                              1, 1, true, 0));
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


void ProbeSetGenotypeTest::doWritePriorFileStyCCS_2_20_1_0() {
  string name = "qt-doWritePriorFileStyCCS_2_20_1_0";
  // shuld get the same answer if using prior generated on same run.

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doWritePriorFileStyCCS_2_20_1_0 "
    "--cel-files ../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt "
    "--write-prior";

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                              0.00001,
                              1, 1, true, 0));
  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "() Phase 1" );
  if(brlmmTest.pass()) {
    string command2 = "./apt-probeset-genotype "
      "--table-output "
      "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
      "--prior-size 1000 --list-sample "
      "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
      "--no-gender-force "
      "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
      "-o test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2 "
      "--cel-files ../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt "
      "--read-priors-brlmm test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.prior.txt";
    
    vector<RegressionCheck *> checks2;
    checks2.push_back(new MatrixCheck("test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2/quant-norm.pm-only.brlmm.calls.txt",
                                "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                                0.0000001,
                                1, 1, true, 0));
    ///@todo relaxed test from 10-5 to 10-4 due to failures. probably due to changes in rounding of reports
    checks2.push_back(new MatrixCheck("test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2/quant-norm.pm-only.brlmm.confidences.txt",
                                "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt",
                                0.0001,
                                1, 1, true, 0));
    RegressionTest brlmmTest2(name.c_str(), command2.c_str(), checks2);
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

void ProbeSetGenotypeTest::doReadGenotypesIn() {
  string name = "qt-doReadGenotypesIn";
    // should get the same answer if using prior generated on same run.

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doReadGenotypesIn "
    "--cel-files ../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt "
    "--dm-out "
    "--dm-thresh .33";

  vector<RegressionCheck *> checks;

  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "() Phase 1" );
  if(brlmmTest.pass()) {
    string command2 = "./apt-probeset-genotype "
      "--table-output "
      "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
      "--prior-size 1000 --list-sample "
      "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2 "
      "--no-gender-force "
      "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
      "-o test-generated/doReadGenotypesIn/phase2 "
      "--cel-files ../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt "
      "--genotypes test-generated/doReadGenotypesIn/dm.calls.txt";

    vector<RegressionCheck *> checks2;
    checks2.push_back(new MatrixCheck(
                                "test-generated/doReadGenotypesIn/phase2/quant-norm.pm-only.brlmm.calls.txt",
                                "test-generated/doReadGenotypesIn/quant-norm.pm-only.brlmm.calls.txt", 
                                0.0000001,
                                1, 1, true, 0));
    checks2.push_back(new MatrixCheck(
                                "test-generated/doReadGenotypesIn/phase2/quant-norm.pm-only.brlmm.confidences.txt", 
                                "test-generated/doReadGenotypesIn/quant-norm.pm-only.brlmm.confidences.txt", 
                                0.00001,
                                1, 1, true, 0));
    RegressionTest brlmmTest2(name.c_str(), command2.c_str(), checks2);
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
  string name = "qt-doStyRvT_2_20_1_0";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample -a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=rvt.MS=2 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doStyRVT_2_20_1_0 "; 
    command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyRVT_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyRVT_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyRVT_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyRVT_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              0.00001,
                              1, 1, true, 0));
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

void ProbeSetGenotypeTest::doStyCES_2_20_1_0() {
  string name = "qt-doStyCES_2_20_1_0";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ces.K=1.MS=2 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doStyCES_2_20_1_0 ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyCES_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyCES_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyCES_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyCES_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              0.00001,
                              1, 1, true, 0));
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

void ProbeSetGenotypeTest::doStyMva_2_20_1_0() {
  string name = "qt-doStyMva_2_20_1_0";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=1.MS=2 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doStyMva_2_20_1_0 ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyMva_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyMva_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyMva_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyMva_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              0.00001,
                              1, 1, true, 0));
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


void ProbeSetGenotypeTest::doStyMva_2_20_1_2() {
  string name = "qt-doStyMva_2_20_1_2";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=1.iterations=1.iter-thresh=.6.MS=2 --no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doStyMva_2_20_1_2 ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyMva_2_20_1_2/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyMva_2_20_1_2/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyMva_2_20_1_2/quant-norm.pm-only.brlmm.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyMva_2_20_1_2/quant-norm.pm-only.brlmm.confidences.txt", 
                              0.00001,
                              1, 1, true, 0));
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


void ProbeSetGenotypeTest::doStyMva_2_20_0_8_0() {
  string name = "qt-doStyMva_2_20_0_8_0";

  string command = valgrindString1 + name + valgrindString2 + "./apt-probeset-genotype "
    "--table-output "
    "-s ../../../regression-data/data/idata/lib/Mapping250K_Sty/set1snps.txt "
    "--prior-size 1000 --list-sample "
    "-a quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=.8.MS=2 "
    "--no-gender-force "
    "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
    "-o test-generated/doStyMva_2_20_0_8_0 ";
  command += Util::joinVectorString(Util::addPrefixSuffix(sty20cels, " ../../../regression-data/data/idata/cel/Mapping250K_Sty/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyMva_2_20_0_8_0/quant-norm.pm-only.brlmm.calls.txt",
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyMva_2_20_0_8_0/quant-norm.pm-only.brlmm.calls.txt",
                              0.0000001,
                              1, 1, true, 0));
  checks.push_back(new MatrixCheck(
                              "test-generated/doStyMva_2_20_0_8_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/doStyMva_2_20_0_8_0/quant-norm.pm-only.brlmm.confidences.txt", 
                              0.00001,
                              1, 1, true, 0));
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


/** Everybody's favorite function. */
int main(int argc, char* argv[]) {
  try {
    ProbeSetGenotypeTest test;
    int testbirdseedflag, testbrlmmflag,testbrlmmpflag;

    testbirdseedflag=testbrlmmflag=testbrlmmpflag = 0;
  
    Verbose::setLevel(2);

    if ( !Fs::dirExists("test-generated/valgrind") ) {
      Fs::mkdirPath("test-generated/valgrind", false);
    }

    if (argc<2)
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
        test.doChpFiles();
        test.doStyCCS_2_20_1_0();
        test.doStyCCS_2_40_0_9_0_033();
        test.doReadGenotypesIn();
        test.doWritePriorFileStyCCS_2_20_1_0();
        test.doFileStyCCS_2_20_1_0();
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
        test.doLabelZStyCCS();
        test.doLabelZStyTestPriors();
        test.doLabelZStyCCSseveral();
        test.doLabelZStyHints();
        
        test.doBrlmmpSnp5Chp();
        test.doBrlmmpSnp5Spf();
        //test.doBrlmmpSnp5KillList();
        ///@todo broken test Martin is workin on fix
        //test.doBrlmmpSnp5CdfReadWriteFeatureEffectsA5();    
    }
    
    if (testbirdseedflag)
    {
        test.doBirdseedSnp6Chp();
        test.doBrlmmpSnp6Chp();
        test.doBirdseedSnp6Spf();
        //    test.doBirdseedSnp6Full();
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

