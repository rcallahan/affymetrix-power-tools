////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
 * @file   test-regression-gender.cpp
 * @brief  Program for doing regression tests on apt-copynumber-gender.
 */

// This has to be last. Otherwise windows compile fails
//#include "calvin_files/utils/src/Calvin.h"
//#include "copynumber/apt-copynumber-cyto/regression/File5EquivalentCheck.h"
//#include "file5/File5.h"
#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/LogStream.h"
#include "util/RegressionCheck.h"
#include "util/RegressionSuite.h"
#include "util/RegressionTest.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
//
using namespace std;

class GenderTest : public RegressionSuite 
{
public:
    int numPassed;
    int numFailed;
    std::string testDir;
  
    GenderTest()
    {
        numPassed = 0;
        numFailed = 0;
    }
    void genderTest1();
};

void GenderTest::genderTest1()
{ 
        string outdir = Fs::join(testDir, "genderTest1");
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;
      
        string name = "genderTest1";
		
        string command = "./apt-copynumber-gender"
                         " --out-dir " + testDir + "/genderTest1"
                         " --xml-file-append-only  ./../../../regression-data/data/lib/gender/genderTest1/genderxml.txt " ;
        
	checks.push_back(new MatrixCheck(       testDir + "/genderTest1/genderCalls.txt",
                                                "../../../regression-data/data/copynumber/gender/genderTest1/genderCalls.txt",
                                                 0.0000001,
                                                 1, 2, false, 0));
 

        RegressionTest genderTest1(name.c_str(), command.c_str(), checks);
        genderTest1.setSuite(*this,  outdir, outdir + "/GenderTest.log", outdir + "/valgrind.log");
        Verbose::out(1, "Doing " + name + "()");
        if(!genderTest1.pass()) 
        {
            Verbose::out(1, "Error in genderTest1::" + name + "(): " + genderTest1.getErrorMsg());
            numFailed++;
        }
        else 
        {
            numPassed++;
        }
} 


int main(int argc, char* argv[]) {
    try {
        /* Set up the logging and message handlers. */
        FsTestDir testDir;
        testDir.setTestDir("copynumber/copynumber-gender", true);
      

        ofstream logOut;
        string logName;
        logName = Fs::join(testDir.asString(), "test-regression-copynumber-gender-run.log");
        Verbose::out(1, "Log file: " + logName);
        Fs::mustOpenToWrite(logOut, logName.c_str());
        LogStream log(3, &logOut);

        Verbose::pushMsgHandler(&log);
        Verbose::pushProgressHandler(&log);
        Verbose::pushWarnHandler(&log);

        Verbose::setLevel(3);
        GenderTest test;
        test.testDir = testDir.asString();
        test.parseArgv(argv);
        test.doValgrind();
        string database = test.getDatabase();

        Verbose::out(1, "Execute regression test: GenderTest");

		test.genderTest1();

        Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed) + " for test-regression-copynumber-gender");
        logOut.close();
        return test.numFailed != 0;
    }
    catch(...) {
        Verbose::out(1,"Unexpected Error: uncaught exception.");
        return 1;
    }
    return 1;
}

