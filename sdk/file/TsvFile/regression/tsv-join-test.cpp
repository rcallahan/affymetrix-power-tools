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
 * @file   tsv-join-test.cpp
 * @brief  Program for doing regression tests on tsv-join.
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

class TsvJoinTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  TsvJoinTest()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void doTranscriptClusterTest();
  void doProbesetIdTest();

private:
  string errorMsg;
};

void TsvJoinTest::doTranscriptClusterTest()
{
  // apt-tsv-join does not make an output directory

  string command = "./apt-tsv-join -k transcript_cluster_id -o " + testDir + "/join-transcript-cluster.txt ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/HuEx-1_0-st-probeset-annot.85000.csv  ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/HuEx-1_0-st-transcript-annot.16000.csv";
  vector<RegressionCheck *> checks;
  RegressionTest joinTest ("doTranscriptClusterTest", command.c_str(), checks);
  Verbose::out (1, "Doing doTranscriptClusterTest()");
  if (!joinTest.run())
  {
    Verbose::out (1, "Error in TsvJoinTest::doTranscriptClusterTest(): " + joinTest.getErrorMsg());
    numFailed++;
  }
  else
  {
    TextFileCheck textCheck (testDir + "/join-transcript-cluster.txt",
                             "../../../regression-data/data/idata/tsv-join/join-transcript-cluster.txt", 34);
    if (textCheck.check (errorMsg))
    {
      numPassed++;
    }
    else
    {
      Verbose::out (1, "Error in TsvJoinTest::doTranscriptClusterTest(): " + errorMsg);
      numFailed++;
    }
  }
}

void TsvJoinTest::doProbesetIdTest()
{

  string command = "./apt-tsv-join -k probeset_id  -o " + testDir + "/join-probeset-id.txt ../../../regression-data/data/idata/lib/HuEx-1_0-st-v2/HuEx-1_0-st-transcript-annot.csv ../../../regression-data/data/idata/tsv-join/quant-norm.pm-gcbg.iter-plier.summary.txt";
  vector<RegressionCheck *> checks;
  RegressionTest joinTest ("doProbesetIdTest", command.c_str(), checks);
  Verbose::out (1, "Doing doProbesetIdTest()");
  if (!joinTest.run())
  {
    Verbose::out (1, "Error in TsvJoinTest::doProbesetIdTest(): " + joinTest.getErrorMsg());
    ++numFailed;
  }
  else
  {
    TextFileCheck textCheck (testDir + "/join-probeset-id.txt",
                             "../../../regression-data/data/idata/tsv-join/join-probeset-id.txt", 34);
    if (textCheck.check (errorMsg))
    {
      numPassed++;
    }
    else
    {
      Verbose::out (1, "Error in TsvJoinTest::doTranscriptClusterTest(): " + errorMsg);
      numFailed++;
    }
  }
}

int main (int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("file/TsvFile/tsv-join",true);
    
    TsvJoinTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);
    test.doProbesetIdTest();
    test.doTranscriptClusterTest();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
