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
 * @file   dump-pgf-test.cpp
 * @brief  Program for doing regression tests on dump-pgf.
 */
#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/RegressionTest.h"
#include "util/TextFileCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>
//

using namespace std;

class DumpPgfTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  DumpPgfTest()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void doDumpPgfTest();

private:
  string errorMsg;
};

void DumpPgfTest::doDumpPgfTest()
{

  string command = "./apt-dump-pgf "
    "-out-file " + testDir + "/dump.probeset_ids_HG-U133_Plus_2.txt "
    "-p ../../../regression-data/data/idata/lib/HG-U133_Plus_2/HG-U133_Plus_2.for-tsv-dump.pgf "
    "-probeset-ids ../../../regression-data/data/idata/lib/HG-U133_Plus_2/probeset_ids_HG-U133_Plus_2.txt";

  vector<RegressionCheck *> checks;
  RegressionTest dumpTest ("doDumpPgfTest", command.c_str(), checks);
  Verbose::out (1, "Doing doDumpPgfTest()");
  if (!dumpTest.run())
  {
    Verbose::out (1, "Error in DumpPgfTest::doDumpPgfTest(): " + dumpTest.getErrorMsg());
    numFailed++;
  }
  else
  {
    TextFileCheck textCheck (testDir + "/dump.probeset_ids_HG-U133_Plus_2.txt",
                             "../../../regression-data/data/idata/dump-pgf/dump.probeset_ids_HG-U133_Plus_2.txt", 5);
    if (textCheck.check (errorMsg))
    {
      numPassed++;
    }
    else
    {
      Verbose::out (1, "Error in DumpPgfTest::doDumpPgfTest(): " + errorMsg);
      numFailed++;
    }
  }
}

int main (int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("file/TsvFile/dump-pgf", true);
    
    DumpPgfTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);
    test.doDumpPgfTest();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
