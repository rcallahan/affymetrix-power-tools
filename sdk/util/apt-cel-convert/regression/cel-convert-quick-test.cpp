////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
 * @file   cel-convert-test.cpp
 * @brief  Program for doing regression tests on cel-extract.
 */

//
#include "util/CelCheck.h"
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

const char *xdaCels[] = {
    "heart-rep1.cel",
    "heart-rep2.cel",
};

const char *textCels[] = {
    "heart-rep1.text.cel",
    "heart-rep2.text.cel",
};

const char *agccCels[] = {
    "heart-rep1.agcc.cel",
    "heart-rep2.agcc.cel",
};

const char *goldPath = "../../../regression-data/data/idata/cel/HG-U133_Plus_2/";

class CelConvertTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  CelConvertTest()
  {
    numPassed = 0;
    numFailed = 0;
  }


  void doConvertRoundRobin();

  void doConvertXDAtoAGCC();
  void doConvertXDAtoTEXT();
  void doConvertAGCCtoXDA();
  void doConvertAGCCtoTEXT();
  void doConvertTEXTtoXDA();
  void doConvertTEXTtoAGCC();

private:
  string errorMsg;

  void runTestHarness(const char *cels[], int size,  const char *inFormat, const char *outFormat);
};

void CelConvertTest::doConvertXDAtoAGCC()  { runTestHarness(xdaCels,ArraySize(xdaCels),"xda", "agcc"); }
void CelConvertTest::doConvertXDAtoTEXT()  { runTestHarness(xdaCels,ArraySize(xdaCels),"xda", "text"); }
void CelConvertTest::doConvertAGCCtoXDA()  { runTestHarness(agccCels,ArraySize(xdaCels),"agcc","xda" ); }
void CelConvertTest::doConvertAGCCtoTEXT() { runTestHarness(agccCels,ArraySize(xdaCels),"agcc","text"); }
void CelConvertTest::doConvertTEXTtoXDA()  { runTestHarness(textCels,ArraySize(xdaCels),"text","xda" ); }
void CelConvertTest::doConvertTEXTtoAGCC() { runTestHarness(textCels,ArraySize(xdaCels),"text","agcc"); }

void CelConvertTest::runTestHarness(const char *cels[], int size, const char *inFormat, const char *outFormat) {

  string testName = "qt-doConvert_" + ToStr(inFormat) + "_to_" + ToStr(outFormat);

  string celString = "";
  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  for(int i = 0; i < size; i++) {
    gen.push_back(Fs::join(testDir , cels[i]));
    gold.push_back(Fs::join(goldPath,  cels[i]));
    celString += " " + Fs::join(goldPath, cels[i]);
  }
  checks.push_back(new CelCheck(gen, gold, .00001));

  string command = "./apt-cel-convert -out-dir " + testDir + " --format " + ToStr(outFormat) + " " + ToStr(celString);

  RegressionTest test(testName.c_str(), command.c_str(), checks);

  Verbose::out (1, "Doing " + testName );
  if (!test.pass()) {
    Verbose::out (1, "Error in CelConvertTest::" + testName + ": " + test.getErrorMsg());
    numFailed++;
  } else {
    numPassed++;
  }

}

void CelConvertTest::doConvertRoundRobin () {
  string testName = "doConvertRoundRobin";

  // Round 1: xda to text
  string command = "./apt-cel-convert -out-dir " + testDir + " --format text " + ToStr(goldPath) + "heart-rep1.cel";
  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  gen.push_back(Fs::join(testDir, "heart-rep1.cel"));
  gold.push_back(Fs::join(goldPath, "heart-rep1.cel"));
  checks.push_back(new CelCheck(gen, gold, .00001));
  RegressionTest test1(testName.c_str(), command.c_str(), checks);

  Verbose::out (1, "Doing " + testName + " (xda to text)");
  if (!test1.pass()) {
    Verbose::out (1, "Error in CelConvertTest::" + testName + ": " + test1.getErrorMsg());
    numFailed++;
  } else {
    numPassed++;
  }

  // Round 2: text to agcc
  command = "./apt-cel-convert --in-place --format agcc " + testDir + "/heart-rep1.cel";
  checks.clear();
  gold.clear();
  gen.clear();
  gen.push_back(Fs::join(testDir , "heart-rep1.cel"));
  gold.push_back(Fs::join(goldPath, "heart-rep1.cel"));
  checks.push_back(new CelCheck(gen, gold, .00001));
  RegressionTest test2(testName.c_str(), command.c_str(), checks);

  Verbose::out (1, "Doing " + testName + " (text to agcc)");
  if (!test2.pass()) {
    Verbose::out (1, "Error in CelConvertTest::" + testName + ": " + test2.getErrorMsg());
    numFailed++;
  } else {
    numPassed++;
  }

  // Round 3: agcc to xda
  command = "./apt-cel-convert --in-place --format xda " + testDir + "/heart-rep1.cel";
  checks.clear();
  gold.clear();
  gen.clear();
  gen.push_back(Fs::join(testDir, "heart-rep1.cel"));
  gold.push_back(Fs::join(goldPath, "heart-rep1.cel"));
  checks.push_back(new CelCheck(gen, gold, .00001));
  RegressionTest test3(testName.c_str(), command.c_str(), checks);

  Verbose::out (1, "Doing " + testName + " (agcc to xda)");
  if (!test3.pass()) {
    Verbose::out (1, "Error in CelConvertTest::" + testName + ": " + test3.getErrorMsg());
    numFailed++;
  } else {
    numPassed++;
  }

}

int main (int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("util/cel-convert-qt", true);
    
    CelConvertTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);

    test.doConvertRoundRobin();

    test.doConvertXDAtoAGCC();
    test.doConvertXDAtoTEXT();
    test.doConvertAGCCtoXDA();
    test.doConvertAGCCtoTEXT();
    test.doConvertTEXTtoXDA();
    test.doConvertTEXTtoAGCC();
    
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
