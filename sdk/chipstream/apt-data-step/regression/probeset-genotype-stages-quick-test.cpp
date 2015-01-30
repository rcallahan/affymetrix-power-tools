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
 * @file   probeset-genotype-stages-test.cpp
 * @author Chuck Sugnet
 * @date   Fri Jan 29 07:01:42 PST 2010
 *
 * @brief  Program for doing regression tests on probeset-genotype data.
 *
 *
 */

//
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

class ProbeSetGenotypeStagesTest {

public:
  int numPassed, numFailed;
  std::string testDir;
  ProbeSetGenotypeStagesTest() {
    numPassed = 0;
    numFailed = 0;
  }

  
  void doMultiChannelAxiom();
    void doBirdseedSnp6Spf();
  string createString(const char **items, int size, const char *prefix, const char *suffix);

};

string ProbeSetGenotypeStagesTest::createString(const char **items, int size, const char *prefix, const char *suffix) {
    string val = "";
    for(int i = 0; i < size; i++) {
        val = val + prefix + items[i] + suffix;
    }
    return val;
}

void ProbeSetGenotypeStagesTest::doBirdseedSnp6Spf() {
  string outdir = testDir + "/qt-doBirdseedSnp6Spf";
  string command = "./apt-probeset-genotype-stages "
    "--table-output "
    "--analysis birdseed "
    "--use-disk=false "
    "--special-snps ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs "
    "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.SNP_A-4232288.spf "
    "--read-models-birdseed ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.birdseed.models "
    "-o " + outdir + " "
    "--cel-files ../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/fas_cel_files.txt";

  string name = "qt-doBirdseedSnp6Spf";
  vector<RegressionCheck *> checks;

  // AW: enough flux in what is going in the report file that we are not going
  //     to check it for every test. just some of the tests.
  // checks for report file
  //checks.push_back(new MatrixCheck(testDir + "/qt-doBirdseedSnp6Spf/birdseed.report.txt",
  //                                 "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.report.txt",
  //                                 0.00001,
  //                                 1,2, false, 0));
  // checks for the table files.
  checks.push_back(new MatrixCheck(testDir + "/qt-doBirdseedSnp6Spf/birdseed.calls.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.calls.txt",
                                   0.0000001,
                                   1, 1, false, 0));
  checks.push_back(new MatrixCheck(testDir + "/qt-doBirdseedSnp6Spf/birdseed.confidences.txt",
                                   "../../../regression-data/data/idata/p-geno/GenomeWideSNP_6/doBirdseedSnp6Spf/birdseed.confidences.txt",
                                   0.0003,
                                   1, 1, false, 0)); // currently 173 nan's in file

  RegressionTest birdseedTest(name.c_str(), command.c_str(), checks);
  //  birdseedTest.setSuite(*this, outdir, outdir + "/apt-probeset-genotype.log", outdir + "/valgrind.log");

  Verbose::out(1, "Doing " + name + "()");
  if(!birdseedTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeTest::" + name + "(): " + birdseedTest.getErrorMsg());
    numFailed++;
  }
  else {
    numPassed++;
  }
}

void ProbeSetGenotypeStagesTest::doMultiChannelAxiom() {

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

  std::string idata="../../../regression-data/data/idata";

  string name = "doMultiChannelAxiom";
  string command = "./apt-probeset-genotype-stages"
    " -a artifact-reduction.ResType=2.Clip=0.4.Close=2.Open=2.Fringe=4.CC=2,quant-norm.target=1000.sketch=5000,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.lambda=1.0.HARD=3.SB=0.75.transform=MVA.copyqc=0.00000.wobble=0.05.MS=0.15.copytype=-1.clustertype=2.ocean=0.00001.CSepPen=0.1.CSepThr=4"
    " --qmethod-spec med-polish.expon=true"
    " --set-gender-method cn-probe-chrXY-ratio"
    " --em-gender false"
    " --female-thresh 0.54"
    " --male-thresh 1.0"
    " --cc-chp-output"
    " --chip-type Axiom_GW_Hu_SNP"
    " --chip-type Axiom_GW_Hu_SNP.r2"
    " --set-analysis-name AxiomGT1"
    " --out-dir " + testDir + "/doMultiChannelAxiom"
    " --force"
    " --read-models-brlmmp "+Fs::Unc(idata+"/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.models")+
    " --spf-file "+Fs::Unc(idata+"/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.spf")+
    " --special-snps "+Fs::Unc(idata+"/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.specialSNPs ")+
    " --chrX-probes "+Fs::Unc(idata+"/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.chrXprobes")+
    " --chrY-probes "+Fs::Unc(idata+"/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.chrYprobes")+
    " --target-sketch "+Fs::Unc(idata+"/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.sketch")+
    " --use-feat-eff "+Fs::Unc(idata+"/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.feature-effects");

  command += " ";
  command +=Util::joinVectorString(Util::addPrefixSuffix(file_prefixes,idata+"/cel/Axiom_GW_Hu_SNP/",".CEL"), " ");

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(file_prefixes,idata+"/p-geno-stages/Axiom_GW_Hu_SNP/doMultiChannelAxiom/cc-chp/", ".AxiomGT1.chp");
  gen = Util::addPrefixSuffix(file_prefixes, testDir + "/doMultiChannelAxiom/cc-chp/", ".AxiomGT1.chp");
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001));
  checks.push_back(new MatrixCheck(testDir + "/doMultiChannelAxiom/AxiomGT1.calls.txt",
                                   "../../../regression-data/data/idata/p-geno-stages/Axiom_GW_Hu_SNP/doMultiChannelAxiom/AxiomGT1.calls.txt",
                                   0.0000001,
                                   1, 1, true, 0));
  checks.push_back(new MatrixCheck(testDir + "/doMultiChannelAxiom/AxiomGT1.confidences.txt",
                                   "../../../regression-data/data/idata/p-geno-stages/Axiom_GW_Hu_SNP/doMultiChannelAxiom/AxiomGT1.confidences.txt",
                                   0.0001,
                                   1, 1, true, 0));

  RegressionTest brlmmTest(name.c_str(), command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!brlmmTest.pass()) {
    Verbose::out(1, "Error in ProbeSetGenotypeStagesTest::" + name + "(): " + brlmmTest.getErrorMsg());
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
    testDir.setTestDir("chipstream/data-step-qt", true);
    
    ProbeSetGenotypeStagesTest test;
    test.testDir = testDir.asString();
    
    int testbirdseedflag, testbrlmmflag,testbrlmmpflag;

    testbirdseedflag=testbrlmmflag=testbrlmmpflag = 0;

    Verbose::setLevel(2);

    // APT-617: fails on all platforms.
    //test.doBirdseedSnp6Spf(); 
    if(Fs::fileExists("JUnitTestResults.qt-doBirdseedSnp6Spf.xml"))
      Fs::rm("JUnitTestResults.qt-doBirdseedSnp6Spf.xml");
    test.doMultiChannelAxiom();
    
    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
