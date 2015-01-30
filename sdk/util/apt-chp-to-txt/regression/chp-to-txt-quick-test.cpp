////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

/**
 * @file   chp-to-txt-test.cpp
 * @brief  Program for doing tests on chp-to-txt.
 */
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

class ChpToTxtTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  ChpToTxtTest()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void doTest(string name, string file, string cdf = "");
  void doBadChipTypeTest(string name, string file, string cdf);

private:
  string errorMsg;
};

void ChpToTxtTest::doTest(string name, string file, string cdf)
{
  string command = "./apt-chp-to-txt -out-dir " + testDir + "/" + name + " ../../../regression-data/data/idata/chp-to-txt/chp/" + file;
  if(cdf != "")
      command += " --cdf-file " + cdf;
  vector<RegressionCheck *> checks;
  checks.push_back(new MixedFileCheck(Fs::join(testDir,name,file + ".txt"),
                                      Fs::join("../../../regression-data/data/idata/chp-to-txt",name,file + ".txt"),
                                      0.00001, 0,0));
  RegressionTest test (name.c_str(), command.c_str(), checks);
  Verbose::out (1, "\nDoing doTest(" + name + ", " + file + ")");
  if (!test.pass())
  {
    Verbose::out (1, "Error in ChpToTxtTest::doTest(" + name + ", " + file + "): " + test.getErrorMsg());
    ++numFailed;
  } else
      ++numPassed;
}

void ChpToTxtTest::doBadChipTypeTest(string name, string file, string cdf)
{
  string command = "./apt-chp-to-txt -v -1 -out-dir " + testDir + "/" + name + " ../../../regression-data/data/idata/chp-to-txt/chp/" + file + " --cdf-file " + cdf;
  vector<RegressionCheck *> checks;
  RegressionTest test (name.c_str(), command.c_str(), checks);
  test.m_NegTest = true;
  Verbose::out (1, "\nDoing doTest(" + name + ", " + file + ")");
  if (test.pass())
  {
    Verbose::setLevel(2);
    ++numPassed;
  } else {
    Verbose::setLevel(2);
    Verbose::out (1, "Error in ChpToTxtTest::doTest(" + name + ", " + file + "): " + test.getErrorMsg());
    Verbose::out (1, "SHOULD HAVE FAILED!!!");
    ++numFailed;
  }
}

int main (int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("util/chp-to-txt-qt", true);
    
    ChpToTxtTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);


    // aggcc checks
    test.doTest("qt-doAgccMdGeno","GenomeWideSNP_5.agcc.multidata--genotyping.brlmm-p.chp");
    test.doTest("qt-doAgccMdCnv","GenomeWideSNP_6.agcc.multidata--copynumber-region.canary-v1.chp");
    test.doTest("qt-doAgccMdCn5","GenomeWideSNP_6.agcc.multidata--copynumber.cn5.chp");
    test.doTest("qt-doAgccMdCn4Cn","Mapping250K.agcc.multidata--copynumber.cn4-cn.chp");
    test.doTest("qt-doAgccMdCn4Loh","Mapping250K.agcc.multidata--copynumber.cn4-loh.chp");
    test.doTest("qt-doAgccExpMas5","HG-U133_Plus_2.agcc.expression.mas5.chp");
    test.doTest("qt-doAgccQuant","HuEx-1_0-st-v2.agcc.quantification.rma.chp");
    test.doTest("qt-doAgccQuantDet","HuEx-1_0-st-v2.agcc.quantification-detection.rma-dabg.chp");

    // xda checks
    test.doTest("qt-doXdaGenoMpam","Mapping10K_Xba131.xda.genotyping.mpam.chp");
    test.doTest("qt-doXdaGenoBrlmm","Mapping250K_Sty.xda.genotyping.brlmm.chp");
    test.doTest("qt-doXdaGenoDm","Mapping50K_Hind240.xda.genotyping.dm.chp");
    test.doTest("qt-doXdaExpMas5","HG-U133_Plus_2.xda.expression.mas5.chp");

    ///@todo regression test for "no-header" option
    ///@todo regression test for "no-body" option
    ///@todo regression test for "chp-files" option

    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}


