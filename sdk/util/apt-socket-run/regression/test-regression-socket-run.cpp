////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
 * @file   test-regression-socket-run.cpp
 * @brief  Program for doing regression tests on sockets.
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

class SocketTest : public RegressionSuite 
{
public:
    int numPassed;
    int numFailed;
    std::string testDir;
  
    SocketTest()
    {
        numPassed = 0;
        numFailed = 0;
    }
    void socketGenotypeTest();
    void socketGenderTest();
};

void SocketTest::socketGenotypeTest()
{ 
        string outdir = testDir + "/socketGenotypeTest";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;
    

        string name = "socketGenotypeTest";
        const char * chpFiles [] = { "F1cross02.brlmm-p",  "F1cross03.brlmm-p",  "F1cross04.brlmm-p", "F1cross05.brlmm-p", 
                                     "Inbred01.brlmm-p", "Inbred02.brlmm-p", "Inbred03.brlmm-p",
                                     "Inbred04.brlmm-p", "Inbred05.brlmm-p", NULL
        };
        std::vector<std::string> result_files;
        result_files = Util::addPrefixSuffix(chpFiles, testDir , ".chp");
        
#ifdef _WIN32
        string command = "apt-socket-run apt-probeset-genotype --xml-file-append-only ..\\..\\..\\regression-data\\data\\lib\\sockets\\socketGenotypeTest\\jobInbred --out-dir " + outdir ;
#else
        string command = "./apt-socket-run ./apt-probeset-genotype --xml-file-append-only ../../../regression-data/data/lib/sockets/socketGenotypeTest/jobInbred --out-dir " + outdir;
#endif
        for ( int i =0; i < result_files.size(); i++ ) {
          command = command + " --result-files " + result_files[i];
        }
        checks.push_back(new MatrixCheck(Fs::join(outdir, "brlmm-p.calls.txt"),
                                                "../../../regression-data/data/util/sockets/socketGenotypeTest/brlmm-p.calls.txt",
                                                 0.0000001,
                                                 1, 1, false, 0));

          

        RegressionTest socketGenotypeTest(name.c_str(), command.c_str(), checks);
        socketGenotypeTest.setSuite(*this,  outdir, outdir + "/apt-socket-run.log", outdir + "/valgrind.log");
        Verbose::out(1, "Doing " + name + "()");
        if(!socketGenotypeTest.pass()) 
        {
            Verbose::out(1, "Error in socketGenotypeTest::" + name + "(): " + socketGenotypeTest.getErrorMsg());
            numFailed++;
        }
        else 
        {
            numPassed++;
        }
} 

void SocketTest::socketGenderTest()
{
        string outdir = testDir + "/socketGenderTest";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        string name = "socketGenderTest";

#ifdef _WIN32
        string command = "apt-socket-run  apt-copynumber-gender --xml-file-append-only ..\\..\\..\\regression-data\\data\\lib\\sockets\\socketGenderTest\\genderxml --out-dir " + outdir;
#else
        string command = "./apt-socket-run ./apt-copynumber-gender --xml-file-append-only ../../../regression-data/data/lib/sockets/socketGenderTest/genderxml.txt --out-dir " + outdir;
#endif

        checks.push_back(new MatrixCheck(Fs::join(outdir, "genderCalls.txt"),
                                         "../../../regression-data/data/util/sockets/socketGenderTest/genderCalls.txt",
                                         0.0000001,
                                         1, 2, false, 0));


        RegressionTest socketGenderTest(name.c_str(), command.c_str(), checks);
        socketGenderTest.setSuite(*this,  outdir, outdir + "/apt-socket-run.log", outdir + "/valgrind.log");
        Verbose::out(1, "Doing " + name + "()");
        if(!socketGenderTest.pass())
        {
            Verbose::out(1, "Error in socketGenderTest::" + name + "(): " + socketGenderTest.getErrorMsg());
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
        testDir.setTestDir("util/socket-run", true);

        ofstream logOut;
        string logName;
        logName = testDir.asString() + "/test-regression-socket-run.log";
        Verbose::out(1, "Log file: " + logName);
        Fs::mustOpenToWrite(logOut, logName.c_str());
        LogStream log(3, &logOut);

        Verbose::pushMsgHandler(&log);
        Verbose::pushProgressHandler(&log);
        Verbose::pushWarnHandler(&log);

        Verbose::setLevel(3);
        SocketTest test;
        test.testDir = testDir.asString();
        
        test.parseArgv(argv);
        test.doValgrind();
        string database = test.getDatabase();

        Verbose::out(1, "Execute regression test: test-regression-socket_run");

	test.socketGenotypeTest();
	test.socketGenderTest();

        Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed) + " for test-regression-socket_run");
        logOut.close();
        return test.numFailed != 0;
    }
    catch(...) {
        Verbose::out(1,"Unexpected Error: uncaught exception.");
        return 1;
    }
    return 1;
}

