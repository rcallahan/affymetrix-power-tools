////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
 * @file   test-regression-copynumber-wave.cpp
 * @brief  Program for doing regression tests on cyto data.
 */

#include "file5/File5.h"
#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/LogStream.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
//

using namespace std;

class test_regression_copynumber_wave
{
public:
    int numPassed;
    int numFailed;
    std::string testDir;
    std::set<std::string> setIgnore;
    std::map<std::string, float> mapEpsilon;

    test_regression_copynumber_wave()
    {
        numPassed = 0;
        numFailed = 0;

        // Header paramters to ignore
        setIgnore.insert("FileCreationTime");
        setIgnore.insert("FileIdentifier");
        setIgnore.insert("program-version");
        setIgnore.insert("create_date");
        setIgnore.insert("create-date");
        setIgnore.insert("affymetrix-algorithm-param-option-verbose");
        setIgnore.insert("affymetrix-algorithm-param-option-exec-guid");
        setIgnore.insert("affymetrix-algorithm-param-option-program-cvs-id");
        setIgnore.insert("affymetrix-algorithm-param-option-version-to-report");
        setIgnore.insert("affymetrix-algorithm-param-option-command-line");
        setIgnore.insert("affymetrix-algorithm-param-option-mem-usage");
        setIgnore.insert("affymetrix-algorithm-param-option-run-probeset-genotype");
        setIgnore.insert("affymetrix-algorithm-param-option-cels");
        setIgnore.insert("affymetrix-algorithm-param-option-out-dir");
        setIgnore.insert("affymetrix-algorithm-param-state-time-start");
        setIgnore.insert("affymetrix-algorithm-param-state-free-mem-at-start");
        setIgnore.insert("affymetrix-algorithm-param-option-temp-dir");

        // Data Set columns to ignore (Specify as Group.Set.Column)
//        setIgnore.insert("AlgorithmData.MarkerABSignal.SCAR");

        // Data Set column epsilon override (Specify as Group.Set.Column)
//        mapEpsilon.insert(std::pair<std::string, float>("AlgorithmData.MarkerABSignal.SCAR", (float)0.0005));
    }

    void execute(const std::string& str)  {
      std::string cmd = Fs::convertCommandToUnc(str);
      if (system(cmd.c_str()) != 0) {Err::errAbort("Execution failed for command: " + cmd);}
    }

    void run()
    {
        execute("./apt-copynumber-wave \
                    -v 3 \
                    --cn-reference-input ../../../regression-data/data/copynumber/copynumber-cyto/Cytogenetics_Array/wave/Cytogenetics_Array_input.REF_MODEL \
                    --cn-reference-output " + testDir + "/Cytogenetics_Array/wave/Cytogenetics_Array_output.REF_MODEL \
                    --analysis additional-waves-reference-method.additional-wave-count=3.trim=2.0.percentile=0.75.demean=false.cn-qc-cutoff=0.27.snp-qc-cutoff=1.1.waviness-seg-count-cutoff=10.force=true.keep-temp-data=false \
                    --cychp-files ../../../regression-data/data/copynumber/copynumber-cyto/Cytogenetics_Array/wave/cychp_files.txt");

        bool bPassed = true;

        std::set<std::string> setIgnore;
        AffxString strFileName1 = "../../../regression-data/data/copynumber/copynumber-cyto/Cytogenetics_Array/wave/Cytogenetics_Array_output.REF_MODEL";
        AffxString strFileName2 = "" + testDir + "/Cytogenetics_Array/wave/Cytogenetics_Array_output.REF_MODEL";
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "Parameters", setIgnore, 0, 1)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "AntigenomicProbes", setIgnore, 0, 1)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "MedianSignals", setIgnore, 0, 1)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "ProbeEffects", setIgnore, 0, 1)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "SNPPosteriors", setIgnore, 0, 1)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "SketchCN", setIgnore, 0, 1)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "SketchSNP", setIgnore, 0, 1)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "WaveCorrection", setIgnore, 0, 1)) {bPassed = false;}

        if (bPassed) {numPassed++;} else {numFailed++;}
    }

};

int main(int argc, char* argv[]) {
  try {

    FsTestDir testDir;
    testDir.setTestDir("copynumber/copynumber-wave", true);

    ofstream logOut;
    string logName;
    logName = Fs::join(testDir.asString(), "/test-regression-copynumber-wave.log");
    Verbose::out(1, "Log file: " + logName);
    Fs::mustOpenToWrite(logOut, logName.c_str());
    LogStream log(3, &logOut);

    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);

    Verbose::setLevel(3);
    test_regression_copynumber_wave test;
    test.testDir = testDir.asString();

    if ( !Fs::dirExists(test.testDir + "/Cytogenetics_Array") ) {
      Fs::mkdirPath(test.testDir + "/Cytogenetics_Array", false);
    }

    Verbose::out(1, "Execute regression test: test-regression-copynumber-wave");

    test.run();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed) + " for test-regression-copynumber-wave");
    logOut.close();
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

