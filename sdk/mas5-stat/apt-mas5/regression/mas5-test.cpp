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
 * @file   mas5-test.cpp
 * @author Pete Klosterman
 * @date   Fri Mar 17 11:34:05 2006
 *
 * @brief  Program for doing regression test on mas5 data.
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

class Mas5Test
{

public:
  int numPassed, numFailed;
  std::string testDir;
  Mas5Test()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void doU133Calvin();
  void doU133Xda();
};

void Mas5Test::doU133Calvin()
{
  string name = "doU133Calvin";
  string command = "./apt-mas5 -c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf --cel-files=../../../regression-data/data/mas5/cel-files.txt  -o " + testDir + "/" + name;

  string chpFilesSuffix = ".chp";
  const char *chpFiles[] = {
                "heart-rep1",
                "heart-rep2",
                "heart-rep3",
                "hela-rep1",
                "hela-rep2",
                "hela-rep3",
                "pancrease-rep1",
                "pancrease-rep2",
                "pancrease-rep3",
                NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/mas5/"+name+"/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, Fs::join( testDir ,  name) + "/", chpFilesSuffix);
  ///@todo APT-378: fix win32 differences in mean signal for P/A/M/All probesets
  // Check including headers at lower eps to deal with win32 differences
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.001, true, 0.0));
  // Check not including headers at higher stringency
  checks.push_back(new CalvinChpCheck(gen, gold, 0, L"apt-", 0.0001, false, 0.0));

  RegressionTest test(name, command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass())
  {
    Verbose::out(1, "Error in Mas5Test::"+name+"(): " + test.getErrorMsg());
    numFailed++;
  }
  else
  {
    numPassed++;
  }
}

void Mas5Test::doU133Xda()
{
  string name = "doU133Xda";
  string command = "./apt-mas5 -c ../../../regression-data/data/lib/HG-U133_Plus_2/HG-U133_Plus_2.cdf --cc-chp-output=false --xda-chp-output=true --cel-files=../../../regression-data/data/mas5/cel-files.txt  -o " + testDir + "/" + name;

  string chpFilesSuffix = ".chp";
  const char *chpFiles[] = {
                "heart-rep1",
                "heart-rep2",
                "heart-rep3",
                "hela-rep1",
                "hela-rep2",
                "hela-rep3",
                "pancrease-rep1",
                "pancrease-rep2",
                "pancrease-rep3",
                NULL
  };

  vector<RegressionCheck *> checks;
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "../../../regression-data/data/mas5/"+name+"/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, Fs::join(testDir, name) + "/", chpFilesSuffix);
  ///@todo APT-378: fix win32 differences in mean signal for P/A/M/All probesets
  // Check including headers at lower eps to deal with win32 differences
  checks.push_back(new ChpCheck(gen, gold, 0, "apt-",1,true,0.0));

  RegressionTest test(name, command.c_str(), checks);
  Verbose::out(1, "Doing " + name + "()");
  if(!test.pass())
  {
    Verbose::out(1, "Error in Mas5Test::"+name+"(): " + test.getErrorMsg());
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
    testDir.setTestDir("mas5-stat/mas5", true);
    Mas5Test test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);

    test.doU133Calvin();
    test.doU133Xda();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
