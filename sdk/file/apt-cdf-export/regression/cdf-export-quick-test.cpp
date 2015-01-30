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
 * @file   cdf-export-test.cpp
 * @brief  Program for doing regression tests on cdf-export.
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

class CdfExportTest
{

public:
  int numPassed, numFailed;
  std::string testDir;
  CdfExportTest()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void exportCdf();

private:
  string errorMsg;
};

void CdfExportTest::exportCdf()
{

  string command = "./apt-cdf-export -c ../../../rawq/test/data/Test3.CDF -ps Pae_16SrRNA_s_at -ps AFFX-Athal-Actin_5_r_at > "
    + testDir + "/export.Test3.txt";
  
  vector<RegressionCheck *> checks;
  RegressionTest exportTest ("exportCdf", command.c_str(), checks);
  Verbose::out (1, "\nDoing exportCdf()");
  if (!exportTest.run())
  {
    Verbose::out (1, "Error in CdfExportTest::exportCdf(): " + exportTest.getErrorMsg());
    numFailed++;
  }
  else
  {
    TextFileCheck exportCheck (testDir + "/export.Test3.txt",
                               "data/export.Test3.txt", 0);
    if (exportCheck.check (errorMsg))
    {
      numPassed++;
    }
    else
    {
      Verbose::out (1, "Error in CdfExportTest::exportCdf(): " + errorMsg);
      numFailed++;
    }
  }
}

int main (int argc, char* argv[])
{
  try {
    FsTestDir testDir;
    testDir.setTestDir("file/cdf-export-qt", true);
    
    CdfExportTest test;
    test.testDir = testDir.asString();
    
    Verbose::setLevel(2);
    test.exportCdf();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
