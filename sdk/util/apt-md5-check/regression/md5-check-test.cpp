////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
 * @file   md5-check-test.cpp
 * @brief  Program for doing regression tests on md5-check.
 */

#include "util/RegressionTest.h"

using namespace std;

class MD5CheckTest
{

public:
  int numPassed, numFailed;
  string md5app;
  MD5CheckTest()
  {
    numPassed = 0;
    numFailed = 0;
    #ifdef _WIN32
    md5app = "./apt-md5-check.exe";
    #else
    md5app = "./apt-md5-check";
    #endif
  }
  void testWithAppParam();
  void testWithSaltParam();
  void testFailure();
  void testFileNotFound();



};

void MD5CheckTest::testWithAppParam()
{
  string command = md5app + " -v 2 -a \"DMET Console\" -i data/translations_md5.rpt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testWithAppParam", command.c_str(), checks);
  Verbose::out (1, "\nDoing testWithAppParam()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MD5CheckTest::testWithAppParam(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void MD5CheckTest::testWithSaltParam()
{
  string command = md5app + " -v 2 -s salt -i data/translations_md5.rpt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testWithSaltParam", command.c_str(), checks);
  Verbose::out (1, "\nDoing testWithSaltParam()");
  if (differencesTest.run())
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MD5CheckTest::testWithSaltParam(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void MD5CheckTest::testFailure()
{
  string command = md5app + " -v 2 -i data/translations_md5.rpt";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testFailure", command.c_str(), checks);
  Verbose::out (1, "\nDoing testFailure()");
  if (differencesTest.run() == false)
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MD5CheckTest::testFailure(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

void MD5CheckTest::testFileNotFound()
{
  string command =  md5app + " -v 2 -i file_does_not_exist";
  vector<RegressionCheck *> checks;
  RegressionTest differencesTest ("testFileNotFound", command.c_str(), checks);
  Verbose::out (1, "\nDoing testFileNotFound()");
  if (differencesTest.run() == false)
    ++numPassed;
  else
  {
    Verbose::out (1, "Error in MD5CheckTest::testFileNotFound(): " + differencesTest.getErrorMsg());
    ++numFailed;
  }
}

int main (int argc, char* argv[])
{
  try {
    MD5CheckTest test;
    Verbose::setLevel(2);
    test.testWithAppParam();
	test.testWithSaltParam();
	test.testFailure();
	//test.testFileNotFound();
    Verbose::out (1, "NumPassed: " + ToStr (test.numPassed) + " NumFailed: " + ToStr (test.numFailed));
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}
