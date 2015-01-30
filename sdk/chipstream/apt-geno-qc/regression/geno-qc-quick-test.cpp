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
* @file   geno-qc-test.cpp
* @author Xiaojun Di
* @date   Nov 2, 2006
*
* @brief  Program for doing regression tests on apt-geno-qc on genotyping data.
*/

//
#include "util/ChpCheck.h"
#include "util/Convert.h"
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

class GenoQcTest {
public:
    int numPassed;
    int numFailed;
    std::string testDir;
  
    GenoQcTest() {
        numPassed = 0;
        numFailed = 0;
    };

    void doGenoQc500k();
    void doGenoQcSnp5();
    void doGenoQcSnp6();
  //
};

void GenoQcTest::doGenoQc500k()
{
    string name = "qt-doGenoQcTest500k";
    string command = "./apt-geno-qc "
      "--dm-out " + testDir + "/qt-doGenoQcTest500k "
      "--spf-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf "
      "--qca-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.qca "
      "--qcc-file ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.qcc "
      "--cel-files ../../../regression-data/data/idata/gqc/Mapping250K_Sty/cel-files.txt "
      "--out-file " + testDir + "/qt-doGenoQcTest500k/geno-qc-sty-test.gqc";

    string dmCallsPrefix = "qt-doGenoQcTest500k/";
    string dmCallsSuffix = ".dm.txt";
    const char *dmCalls[] = {
        "NA06985_B01_Sty_Plate1", "NA06991_B03_Sty_Plate1", "NA06993_B02_Sty_Plate1",
        "NA06994_A11_Sty_Plate1", "NA07000_A10_Sty_Plate1",
    NULL
    };

    vector<RegressionCheck *> checks;
    checks.push_back(new MixedFileCheck(
                testDir + "/qt-doGenoQcTest500k/geno-qc-sty-test.gqc",
                "../../../regression-data/data/idata/gqc/Mapping250K_Sty/qt-doGenoQcTest500k/geno-qc-sty-test.gqc",
                0.0001, 1, 0));
    vector<string> dmTest = Util::addPrefixSuffix(dmCalls, testDir + "/" + dmCallsPrefix, dmCallsSuffix);
    vector<string> dmGold = Util::addPrefixSuffix(dmCalls, "../../../regression-data/data/idata/gqc/Mapping250K_Sty/" + dmCallsPrefix, dmCallsSuffix);
    for(int i=0; i<dmTest.size(); i++)
        checks.push_back(new MixedFileCheck(dmTest[i], dmGold[i], 0.0001, 1, 0));

    RegressionTest qcTest(name,command, checks);
    Verbose::out(1, "Doing " + name + "()");
    if(!qcTest.pass())
    {
        Verbose::out(1, "Error in GenoQcTest::" + name + "(): " + qcTest.getErrorMsg());
        numFailed++;
    }
    else
    {
        numPassed++;
    }
}

void GenoQcTest::doGenoQcSnp5()
{
    string name = "qt-doGenoQcTestSnp5";
    string command = "./apt-geno-qc "
      "--dm-out " + testDir + "/qt-doGenoQcTestSnp5 "
      "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.spf "
      "--qca-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.qca "
      "--qcc-file ../../../regression-data/data/idata/lib/GenomeWideSNP_5/GenomeWideSNP_5.qcc "
      "--cel-files ../../../regression-data/data/idata/gqc/GenomeWideSNP_5/cel-files.txt "
      "--out-file " + testDir + "/qt-doGenoQcTestSnp5/geno-qc-snp5-test.gqc";

    string dmCallsPrefix = "qt-doGenoQcTestSnp5/";
    string dmCallsSuffix = ".dm.txt";
    const char *dmCalls[] = {
        "NA06985_GW5_C", "NA06991_GW5_C", "NA06993_GW5_C",
        "NA06994_GW5_C", "NA07000_GW5_C",
    NULL
    };

    vector<RegressionCheck *> checks;
    checks.push_back(new MixedFileCheck(
                testDir + "/qt-doGenoQcTestSnp5/geno-qc-snp5-test.gqc",
                "../../../regression-data/data/idata/gqc/GenomeWideSNP_5/qt-doGenoQcTestSnp5/geno-qc-snp5-test.gqc",
                0.0001, 1, 0));

    vector<string> dmTest = Util::addPrefixSuffix(dmCalls, testDir + "/" + dmCallsPrefix, dmCallsSuffix);
    vector<string> dmGold = Util::addPrefixSuffix(dmCalls, "../../../regression-data/data/idata/gqc/GenomeWideSNP_5/" + dmCallsPrefix, dmCallsSuffix);
    for(int i=0; i<dmTest.size(); i++)
        checks.push_back(new MixedFileCheck(dmTest[i], dmGold[i], 0.0001, 1, 0));

    RegressionTest qcTest(name, command, checks);
    Verbose::out(1, "Doing " + name + "()");
    if(!qcTest.pass())
    {
        Verbose::out(1, "Error in GenoQcTest::" + name + "(): " + qcTest.getErrorMsg());
        numFailed++;
    }
    else
    {
        numPassed++;
    }
}

void GenoQcTest::doGenoQcSnp6()
{
    string name = "qt-doGenoQcTestSnp6";
    string command = "./apt-geno-qc "
      "--dm-out " + testDir + "/qt-doGenoQcTestSnp6 "
      "--spf-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.spf "
      "--qca-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.r2.qca "
      "--qcc-file ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.r2.qcc "
      "--chrX-probes ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes "
      "--chrY-probes ../../../regression-data/data/idata/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes "
      "--cel-files  ../../../regression-data/data/idata/gqc/GenomeWideSNP_6/fas_cel_files.txt "
      "--out-file " + testDir + "/qt-doGenoQcTestSnp6/geno-qc-snp6-test.gqc";

    string dmCallsPrefix = "qt-doGenoQcTestSnp6/";
    string dmCallsSuffix = ".dm.txt";
    const char *dmCalls[] = {
        "NA06985_GW6_C", "NA06993_GW6_C", "NA07000_GW6_C",
        "NA07019_GW6_C", "NA07029_GW6_C",
    NULL
    };

    vector<RegressionCheck *> checks;
    checks.push_back(new MixedFileCheck(
                testDir + "/qt-doGenoQcTestSnp6/geno-qc-snp6-test.gqc",
                "../../../regression-data/data/idata/gqc/GenomeWideSNP_6/qt-doGenoQcTestSnp6/geno-qc-snp6-test.gqc",
                0.0001, 1, 0));

    vector<string> dmTest = Util::addPrefixSuffix(dmCalls, testDir + "/" + dmCallsPrefix, dmCallsSuffix);
    vector<string> dmGold = Util::addPrefixSuffix(dmCalls, "../../../regression-data/data/idata/gqc/GenomeWideSNP_6/" + dmCallsPrefix, dmCallsSuffix);
    for(int i=0; i<dmTest.size(); i++)
        checks.push_back(new MixedFileCheck(dmTest[i], dmGold[i], 0.0001, 1, 0));

    RegressionTest qcTest(name, command, checks);
    Verbose::out(1, "Doing " + name + "()");
    if(!qcTest.pass())
    {
        Verbose::out(1, "Error in GenoQcTest::" + name + "(): " + qcTest.getErrorMsg());
        numFailed++;
    }
    else
    {
        numPassed++;
    }
}

//////////

/** the simple and short main function. */
int main(int argc,const char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("chipstream/geno-qc-qt",true);
    
    GenoQcTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);

    // this one is the slowest -- as we are running over more probesets w/ mm probes
    test.doGenoQc500k();
    test.doGenoQcSnp5();
    test.doGenoQcSnp6();

    ///@todo add integration tests for axiom

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
