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
 * @file   test-regression-copynumber-workflow-quick.cpp
 * @brief  Program for doing regression tests on probeset-genotype data.
 */

#include "calvin_files/utils/src/Calvin.h"

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

string valgrindString1 = "valgrind --leak-check=yes --log-file=";
string valgrindString2 = ".valgrind.txt ";

class test_regression_copynumber_workflow_quick
{
public:
    int numPassed;
    int numFailed;
    std::string testDir;
    std::map<std::string, float> mapEpsilon;
    std::set<std::string> setIgnore;
    std::vector<std::string> vSixFiles;
    std::vector<std::string> vOneFile;

    test_regression_copynumber_workflow_quick()
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

        // Data Set columns to ignore (Specify as Group.Set.Column)
        setIgnore.insert("MultiData.CopyNumber.SmoothSignal");
        setIgnore.insert("Segments.CN.SegmentID");
        setIgnore.insert("Segments.LOH.SegmentID");
        setIgnore.insert("Segments.CNNeutralLOH.SegmentID");
        setIgnore.insert("Segments.NormalDiploid.SegmentID");
        setIgnore.insert("Segments.Mosaicism.SegmentID");
        setIgnore.insert("Segments.NoCall.SegmentID");

        // File tags to check for equivalency
        vSixFiles.push_back("NA06985_GW6_C");
        vSixFiles.push_back("NA06993_GW6_C");
        vSixFiles.push_back("NA07000_GW6_C");
        vSixFiles.push_back("NA07019_GW6_C");
        vSixFiles.push_back("NA07029_GW6_C");
        vSixFiles.push_back("NA07056_GW6_C");

        vOneFile.push_back("NA07029_GW6_C");
    }

    void execute(const std::string& str) {
      std::string cmd =Fs::convertCommandToUnc(str);
      if (system(cmd.c_str()) != 0) {Err::errAbort("Execution failed for command: " + cmd);}
    }

    void equivalency(const std::string& name, std::vector<string>& vFileTags)
    {
        bool bPassed = true;
        AffxString strFileName1;
        AffxString strFileName2;
        Verbose::out(1, "*");
        for(unsigned int i = 0; (i < vFileTags.size()); i++)
        {
            strFileName1 = "../../../regression-data/data/idata/copy-number/GenomeWideSNP_6/" + name + "/" + vFileTags[i] + ".CN5.CNCHP";
            strFileName2 = testDir + "/GenomeWideSNP_6/" + name + "/" + vFileTags[i] + ".CN5.CNCHP";
            if (!Calvin::equivalent(strFileName1, strFileName2, setIgnore, mapEpsilon, 0.0001, 1.0, false)) {bPassed = false;}

        }
        Verbose::out(1, "*");
        if (bPassed) {numPassed++;} else {numFailed++;}
    }

    void run()
    {

        copynumber_HapMap270_1_cels();
        copynumber_HapMap270_6_cels();
        reference_6_cels();
        copynumber_1_cels();

    }



    void copynumber_HapMap270_1_cels()
    {
                string name = "copynumber_HapMap270_1_cels";
                execute( valgrindString1 + name + valgrindString2 + " ./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/GenomeWideSNP_6/copynumber_HapMap270_1_cels \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ");

        equivalency("copynumber_HapMap270_1_cels", vOneFile);
    }


    void copynumber_HapMap270_6_cels()
    {
                string name = "copynumber_HapMap270_6_cels";
                execute( valgrindString1 + name + valgrindString2 + " ./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/GenomeWideSNP_6/copynumber_HapMap270_6_cels \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        equivalency("copynumber_HapMap270_6_cels", vSixFiles);
    }

    void reference_6_cels()
    {
                string name = "reference_6_cels";
        Fs::rm(testDir + "/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref", false);
                execute( valgrindString1 + name + valgrindString2 + " ./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-output " + testDir + "/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/GenomeWideSNP_6/reference_6_cels \
             --temp-dir " + testDir + "/GenomeWideSNP_6/reference_6_cels \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        equivalency("reference_6_cels", vSixFiles); // Reference building
    }


    void copynumber_1_cels()
    {
                string name = "copynumber_1_cels";
                execute( valgrindString1 + name + valgrindString2 + " ./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input " + testDir + "/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/GenomeWideSNP_6/copynumber_1_cels \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ");

        equivalency("copynumber_1_cels", vOneFile);
    }
};

int main(int argc, char* argv[])
{
  try {

    FsTestDir testDir;
    testDir.setTestDir("copynumber/copynumber-workflow-valgrind", true);

    ofstream logOut;
    string logName;
    logName = testDir.asString() + "/test-regression-copynumber-workflow.log";
    valgrindString1 = valgrindString1 + testDir.asString() + "/";
    
    Verbose::out(1, "Log file: " + logName);
    Fs::mustOpenToWrite(logOut, logName.c_str());
    LogStream log(3, &logOut);
    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);

    Verbose::setLevel(3);
    test_regression_copynumber_workflow_quick test;
    test.testDir = testDir.asString();

    if ( !Fs::dirExists(test.testDir + "/GenomeWideSNP_6") ) {
      Fs::mkdirPath(test.testDir + "/GenomeWideSNP_6", false);
    }

    Verbose::out(1, "Execute regression test: test-regression-copynumber-workflow");

    test.run();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed) + " for test-regression-copynumber-workflow");
    logOut.close();
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

