////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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
 * @file   cel-transformer-test.cpp
 * @author Pete Klosterman
 * @date   Fri Mar 17 11:34:05 2006
 *
 * @brief  Program for doing regression test on cel-transformer data.
 *
 */

//
#include "util/CelCheck.h"
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

class CelTransformerTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  CelTransformerTest()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void doSpfTest();
  void doAgccTest();
  void doPgfTest();
  void doA5SketchTest();
  void doMultiChanTest();
};

void CelTransformerTest::doSpfTest()
{
  string command = "./apt-cel-transformer -c rma-bg,quant-norm -o " + testDir + "/qt-doSpfTest --spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf";

  string celFilesSuffix = ".cel";
  const char *celFiles[] = {"heart-rep1","heart-rep2","heart-rep3",
                NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(celFiles, "../../../regression-data/data/idata/cel-tran/qt-doSpfTest/", celFilesSuffix);
  gen = Util::addPrefixSuffix(celFiles, testDir + "/qt-doSpfTest/", celFilesSuffix);

  command += Util::joinVectorString(Util::addPrefixSuffix(celFiles, " ../../../regression-data/data/idata/cel/HG-U133_Plus_2/", celFilesSuffix), " ");

  checks.push_back(new CelCheck(gen, gold, .2));

  RegressionTest test("qt-doSpfTest", command.c_str(), checks);
  Verbose::out(1, "Doing doSpfTest()");
  if(!test.pass())
  {
    Verbose::out(1, "Error in CelTransformerTest::doSpfTest(): " + test.getErrorMsg());
    numFailed++;
  }
  else
  {
    numPassed++;
  }
}

void CelTransformerTest::doAgccTest()
{
  string command = "./apt-cel-transformer "
    "-c rma-bg,quant-norm "
    "-o " + testDir + "/qt-doAgccTest "
    "--spf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.spf ";

  string celFilesGoldSuffix = ".agcc.cel";
  string celFilesSuffix = ".agcc.cel";
  const char *celFilesGold[] = {"heart-rep1","heart-rep2","heart-rep3",
                NULL
  };
  const char *celFiles[] = {"heart-rep1","heart-rep2","heart-rep3",
                NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(celFilesGold, "../../../regression-data/data/idata/cel-tran/qt-doAgccTest/", celFilesGoldSuffix);
  gen = Util::addPrefixSuffix(celFiles, testDir + "/qt-doAgccTest/", celFilesSuffix);
  command += Util::joinVectorString(Util::addPrefixSuffix(celFiles, " ../../../regression-data/data/idata/cel/HG-U133_Plus_2/", celFilesSuffix), " ");

  checks.push_back(new CelCheck(gen, gold, .2));

  RegressionTest test("qt-doAgccTest", command.c_str(), checks);
  Verbose::out(1, "Doing doAgccTest()");
  if(!test.pass())
  {
    Verbose::out(1, "Error in CelTransformerTest::doAgccTest(): " + test.getErrorMsg());
    numFailed++;
  }
  else
  {
    numPassed++;
  }
}

void CelTransformerTest::doPgfTest()
{
  string command = "./apt-cel-transformer "
    "-c rma-bg,quant-norm "
    "-o " + testDir + "/qt-doPgfTest "
    "--pgf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.pgf "
    "--clf-file ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.clf "
    "--cel-files ../../../regression-data/data/idata/cel-tran/cel-files.txt";

  string celFilesSuffix = ".cel";
  const char *celFiles[] = {"heart-rep1","heart-rep2","heart-rep3",
                NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(celFiles, "../../../regression-data/data/idata/cel-tran/qt-doPgfTest/", celFilesSuffix);
  gen = Util::addPrefixSuffix(celFiles, testDir + "/qt-doPgfTest/", celFilesSuffix);

  checks.push_back(new CelCheck(gen, gold, .2));

  RegressionTest test("qt-doPgfTest", command.c_str(), checks);
  Verbose::out(1, "Doing doPgfTest()");
  if(!test.pass())
  {
    Verbose::out(1, "Error in CelTransformerTest::doPgfTest(): " + test.getErrorMsg());
    numFailed++;
  }
  else
  {
    numPassed++;
  }
}

void CelTransformerTest::doA5SketchTest()
{
  string command = "./apt-cel-transformer "
    "-c quant-norm.sketch=1000 "
    "-o " + testDir + "/qt-doA5SketchTest "
    "--write-sketch "
    "--a5-write-sketch ";

  string celFilesSuffix = ".cel";
  const char *celFiles[] = {"heart-rep1","heart-rep2","heart-rep3",
                NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  command += Util::joinVectorString(Util::addPrefixSuffix(celFiles, " ../../../regression-data/data/idata/cel/HG-U133_Plus_2/", celFilesSuffix), " ");

  gold = Util::addPrefixSuffix(celFiles, "../../../regression-data/data/idata/cel-tran/qt-doA5SketchTest/", celFilesSuffix);
  gen = Util::addPrefixSuffix(celFiles, testDir + "/qt-doA5SketchTest/", celFilesSuffix);
  checks.push_back(new CelCheck(gen, gold, .2));

  ///@todo add check for txt sketch output file
  ///@todo add check for a5 sketch output file

  RegressionTest test("qt-doA5SketchTest", command.c_str(), checks);
  Verbose::out(1, "Doing doA5SketchTest()");
  if(!test.pass())
  {
    Verbose::out(1, "Error in CelTransformerTest::doA5SketchTest(): " + test.getErrorMsg());
    numFailed++;
  }
  else
  {
    numPassed++;
  }
}

void CelTransformerTest::doMultiChanTest() {
  string command = "./apt-cel-transformer "
    "-c rma-bg,quant-norm "
    "-o " + testDir + "/qt-doMultiChanTest "
    "--spf-file ../../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.spf ";

  const char* celPrefixes[] = {
    "NA07029_AxiomGWASHuSNP1_20090812_MCKI_T01_B11_v1",
    "NA07034_AxiomGWASHuSNP1_20090812_MCKI_T01_A05_v1",
    "NA10831_AxiomGWASHuSNP1_20090812_MCKI_T01_D06_v1",
    NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;

  gold = Util::addPrefixSuffix(celPrefixes, "../../../regression-data/data/idata/cel-tran/qt-doMultiChanTest/", ".CEL");
  gen = Util::addPrefixSuffix(celPrefixes, testDir + "/qt-doMultiChanTest/", ".CEL");
  command += Util::joinVectorString(Util::addPrefixSuffix(celPrefixes, " ../../../regression-data/data/idata/cel/Axiom_GW_Hu_SNP/", ".CEL"), " ");
  // TODO: allow CelCheck to work on multi-channel CELs.  in it's current state
  // the test serendipitously passes
  checks.push_back(new CelCheck(gen, gold, .2));

  RegressionTest test("qt-doMultiChanTest", command.c_str(), checks);
  Verbose::out(1, "Doing doMultiChanTest()");
  if(!test.pass())
  {
    Verbose::out(1, "Error in CelTransformerTest::doMultiChanTest(): " + test.getErrorMsg());
    numFailed++;
  }
  else
  {
    numPassed++;
  }
}

int main(int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("chipstream/cel-transformer-qt", true);
    
    CelTransformerTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);

    test.doAgccTest();
    test.doSpfTest();
    test.doPgfTest();
    test.doA5SketchTest();
    test.doMultiChanTest();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
