////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
 * @file   summary-vis-test.cpp
 * @brief  Program for doing regression tests on summary-vis.
 */
#include "util/Fs.h"
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

class SummaryVisTest
{

public:
  int numPassed, numFailed;
  SummaryVisTest()
  {
    numPassed = 0;
    numFailed = 0;
  }
  void doSummaryVisTest();

private:
  string errorMsg;
};

void SummaryVisTest::doSummaryVisTest()
{
  // apt-summary-vis does not make an output directory
  if ( !Fs::dirExists("summary-vis-test-generated") ) {
    Fs::mkdir("summary-vis-test-generated", false);
  }

  string command = "./apt-summary-vis -o summary-vis-test-generated/summaryvis.txt -g ../../../regression-data/data/lib/HuEx-1_0-st-v2/HuEx-1_0-st-v2.annot.hg16.csv ../../../regression-data/data/chipstream/summary-vis/HuEx-1_0-st-v2/exon-summary.txt  ../../../regression-data/data/chipstream/summary-vis/HuEx-1_0-st-v2/exon-summary0.txt  ../../../regression-data/data/chipstream/summary-vis/HuEx-1_0-st-v2/exon-summary1.txt";
  vector<RegressionCheck *> checks;
  RegressionTest summaryTest ("doSummaryVisTest", command.c_str(), checks);
  Verbose::out (1, "Doing doSummaryVisTest()");
  if (!summaryTest.run())
  {
    Verbose::out (1, "Error in SummaryVisTest::doSummaryVisTest(): " + summaryTest.getErrorMsg());
    numFailed++;
  }
  else
  {
    TextFileCheck textCheck ("summary-vis-test-generated/summaryvis.txt",
                             "../../../regression-data/data/chipstream/summary-vis/HuEx-1_0-st-v2/summaryvis.txt", 41);
    if (textCheck.check (errorMsg))
    {
      numPassed++;
    }
    else
    {
      Verbose::out (1, "Error in SummaryVisTest::doSummaryVisTest(): " + errorMsg);
      numFailed++;
    }
  }
}

int main (int argc, char* argv[])
{
  try {
    SummaryVisTest test;
    Verbose::setLevel(2);
    test.doSummaryVisTest();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
