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

class test_regression_copynumber_workflow_quick  : public RegressionSuite, public RegressionCheck
{
public:
    int numPassed;
    int numFailed;
    std::string testDir;
    std::map<std::string, float> mapEpsilon;
    std::set<std::string> setIgnore;
    std::vector<std::string> vSixFiles;
    std::vector<std::string> vOneFile;

  std::string m_OutPrefix;
  std::vector<std::string> m_ToCompare;

  bool check(std::string &msg) {
    return equivalency(m_OutPrefix, m_ToCompare);
  }

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
      std::string cmd = Fs::convertCommandToUnc(str);
      if (system(cmd.c_str()) != 0) {Err::errAbort("Execution failed for command: " + cmd);}
    }

    bool equivalency(const std::string& name, std::vector<string>& vFileTags)
    {
        bool bPassed = true;
        AffxString strFileName1;
        AffxString strFileName2;
        Verbose::out(1, "*");
        for(unsigned int i = 0; (i < vFileTags.size()); i++)
        {
            strFileName1 = "../../../regression-data/data/idata/copy-number/GenomeWideSNP_6/" + name + "/" + vFileTags[i] + ".CN5.CNCHP";
            strFileName2 = Fs::join(testDir, "qt-GenomeWideSNP_6", name ) + "/" + vFileTags[i] + ".CN5.CNCHP";
            if (!Calvin::equivalent(strFileName1, strFileName2, setIgnore, mapEpsilon, 0.0001, 1.0, false)) {bPassed = false;}

//            strFileName1 = "../../../regression-data/data/idata/copy-number/GenomeWideSNP_6/" + name + "/" + vFileTags[i] + ".CN5.cyto2.cychp";
//            strFileName2 = "test-generated/qt-GenomeWideSNP_6/" + name + "/" + vFileTags[i] + ".CN5.cyto2.cychp";
//            if (!Calvin::equivalent(strFileName1, strFileName2, setIgnore, mapEpsilon, 0.0001, false)) {bPassed = false;}
        }
        Verbose::out(1, "*");
        if (bPassed) {numPassed++;} else {numFailed++;}
        return bPassed;
    }

    void run()
    {
//        AffxString strFileName1 = "test-generated/qt-GenomeWideSNP_6/reference_6_cels/NA06985_GW6_C.CN5.CNCHP";
//        AffxString strFileName2 = "../../../regression-data/data/idata/copy-number/GenomeWideSNP_6/reference_6_cels/NA06985_GW6_C.CN5.CNCHP";
//        AffxString strFileName1 = "\\\\10.80.199.50\\share\\TestData\\SNP6withARR\\3.0.2_Batch_NA28_Default_GC\\NA06985_GW6_C.CN5.CNCHP";
//        AffxString strFileName2 = "\\\\10.80.199.50\\share\\TestData\\SNP6withARR\\4.0_build296_x64_Batch_NA28_Default_GC\\NA06985_GW6_C.CN5.CNCHP";
//        AffxString strFileName1 = "\\\\10.80.199.50\\share\\TestData\\SNP6withARR\\3.0.2_Single_NA28_Default_GC\\NA06985_GW6_C.CN5.CNCHP";
//        AffxString strFileName2 = "\\\\10.80.199.50\\share\\TestData\\SNP6withARR\\4.0_build296_x64_Single_NA28_Default_GC\\NA06985_GW6_C.CN5.CNCHP";
//        Calvin::equivalent(strFileName1, strFileName2, setIgnore, mapEpsilon, 0.0001, 1.0, false);
//        return;
//        copynumber_HapMap270_1_cels_probeset_ids();
//        copynumber_HapMap270_1_cels_cancer();
//        copynumber_HapMap270_for_Carl();

//        reference_6_cels_med_norm();
//        copynumber_1_cels_med_norm();
//        return;

//APT-619
#ifndef _WIN32      
        copynumber_HapMap270_1_cels();
        copynumber_HapMap270_6_cels();
#endif        
        reference_6_cels();
        copynumber_1_cels();

//        copynumber_6_cels();
//        copynumber_6_cels_data_paging();
//        copynumber_6_cels_gc_correction_false();

//        component_6_cels();
    }

    void copynumber_HapMap270_1_cels_cancer()
    {
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output true \
             --cnchp-output true \
             --text-output true \
             --reference-input ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/qt-GenomeWideSNP_6/copynumber_HapMap270_1_cels_cancer \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ");

        equivalency("copynumber_HapMap270_1_cels_cancer", vOneFile);
    }

    void copynumber_HapMap270_for_Carl()
    {
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --run-probeset-genotype false \
             --run-geno-qc true \
             --cychp-output true \
             --cnchp-output true \
             --text-output false \
             --qca-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.qca \
             --qcc-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.qcc \
             --reference-input ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/qt-GenomeWideSNP_6/copynumber_HapMap270_for_Carl \
             --cel-files ../../../../../../cyto2/cyto2_cel-files/upd_trio/CELFileList.txt \
             ");
    }

    void copynumber_HapMap270_1_cels()
    {
      string outdir = testDir + "/qt-GenomeWideSNP_6/copynumber_HapMap270_1_cels";
      //Fs::rm_rf(outdir, false);
      vector<RegressionCheck *> checks;
      m_OutPrefix = "copynumber_HapMap270_1_cels";
      m_ToCompare = vOneFile;
      checks.push_back(this);
      string name = "qt-copynumber_HapMap270_1_cels";
      string command = "./apt-copynumber-workflow \
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
             --o " + outdir + " \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ";
      RegressionTest test(name.c_str(), command.c_str(), checks, false);
      test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
      if(!test.pass()) {
        Verbose::out(1, "Error in::" + name + "(): " + test.getErrorMsg());
      }
    }

    void copynumber_HapMap270_1_cels_probeset_ids()
    {
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --probeset-ids ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/ProbeSetRestrictList.txt \
             --run-geno-qc false \
             --cychp-output true \
             --cnchp-output true \
             --text-output true \
             --reference-input ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/qt-GenomeWideSNP_6/copynumber_HapMap270_1_cels_probeset_ids \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ");

        equivalency("copynumber_HapMap270_1_cels_probeset_ids", vOneFile);
    }

    void copynumber_HapMap270_6_cels()
    {
      string outdir = testDir + "/qt-GenomeWideSNP_6/copynumber_HapMap270_6_cels";
      //Fs::rm_rf(outdir, false);
      vector<RegressionCheck *> checks;
      m_OutPrefix =  "copynumber_HapMap270_6_cels";
      m_ToCompare = vSixFiles;
      checks.push_back(this);
      string name = "qt-copynumber_HapMap270_6_cels";
      string command = "./apt-copynumber-workflow \
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
             --o " + outdir + " \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ";
      RegressionTest test(name.c_str(), command.c_str(), checks, false);
      test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
      if(!test.pass()) {
        Verbose::out(1, "Error in::" + name + "(): " + test.getErrorMsg());
      }

    }

    void reference_6_cels()
    {
      string outdir = testDir + "/qt-GenomeWideSNP_6/reference_6_cels";
      //Fs::rm_rf(outdir, false);
      vector<RegressionCheck *> checks;
      m_OutPrefix = "reference_6_cels";
      m_ToCompare = vSixFiles;
      checks.push_back(this);
      string name = "qt-reference_6_cels";
        Fs::rm(testDir + "/qt-GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref", false);
       string command = "./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-output " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + outdir + " \
             --temp-dir " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ";

      RegressionTest test(name.c_str(), command.c_str(), checks, false);
      test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
      if(!test.pass()) {
        Verbose::out(1, "Error in::" + name + "(): " + test.getErrorMsg());
      }
    }

    void reference_6_cels_med_norm()
    {
        Fs::rm(testDir + "/qt-GenomeWideSNP_6/reference_6_cels_med_norm/CNReference.a5.ref", false);
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --normalization-type 2 \
             --run-geno-qc false \
             --cychp-output true \
             --cnchp-output true \
             --text-output true \
             --reference-output " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels_med_norm/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels_med_norm \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

//        equivalency("reference_6_cels_med_norm", vSixFiles); // Reference building
    }

    void copynumber_6_cels()
    {
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/qt-GenomeWideSNP_6/copynumber_6_cels \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        equivalency("copynumber_6_cels", vSixFiles);
    }

    void copynumber_6_cels_data_paging()
    {
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --run-probeset-genotype false \
             --run-geno-qc false \
             --cychp-output true \
             --cnchp-output true \
             --text-output true \
             --reference-input " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/qt-GenomeWideSNP_6/copynumber_6_cels \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        equivalency("copynumber_6_cels", vSixFiles);
    }

    void copynumber_1_cels()
    {
      string outdir = testDir + "/qt-GenomeWideSNP_6/copynumber_1_cels";
      //Fs::rm_rf(outdir, false);
      vector<RegressionCheck *> checks;
      m_OutPrefix = "copynumber_1_cels";
      m_ToCompare = vOneFile;
      checks.push_back(this);
      string name = "qt-copynumber_1_cels";
      string command = "./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + outdir + " \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ";
      RegressionTest test(name.c_str(), command.c_str(), checks, false);
      test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
      if(!test.pass()) {
        Verbose::out(1, "Error in::" + name + "(): " + test.getErrorMsg());
      }
    }

    void copynumber_1_cels_med_norm()
    {
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --normalization-type 2 \
             --run-geno-qc false \
             --cychp-output true \
             --cnchp-output true \
             --text-output true \
             --reference-input " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels_med_norm/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/qt-GenomeWideSNP_6/copynumber_1_cels_med_norm \
             ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ");

//        equivalency("copynumber_1_cels_med_norm", vOneFile);
    }

    void copynumber_6_cels_gc_correction_false()
    {
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --analysis log2-ratio.gc-correction=false \
             --analysis gaussian-smooth \
             --analysis allelic-difference \
             --analysis cn-state.hmmCN_mu=\"'\"-2,-0.552,0,0.339,0.543\"'\" \
             --analysis loh \
             --analysis cn-neutral-loh \
             --analysis normal-diploid \
             --cychp-output true --cnchp-output true --text-output true \
             --reference-input " + testDir + "/qt-GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + testDir + "/qt-GenomeWideSNP_6/copynumber_6_cels_gc_correction_false \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        equivalency("copynumber_6_cels_gc_correction_false", vSixFiles);
    }

    void component_6_cels()
    {
        execute("./apt-probeset-genotype \
             -v 1 \
             --all-types true \
             --set-gender-method cn-probe-chrXY-ratio \
             --a5-global-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CNReference.a5 \
             --analysis adapter-type-norm.StyOnlySnpRatio=0.8287.NspOnlySnpRatio=0.7960.NspOnlyCnRatio=1.4218.NspBothSnpRatio=0.9954.NspBothCnRatio=1.6392,quant-norm,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.wobble=0.05.MS=1 \
             --qmethod-spec plier.optmethod=1.FixFeatureEffect=true \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --out-dir " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             --table-output false \
             --summaries false \
             --xda-chp-output false \
             --cc-chp-output false \
             --residuals false \
             --set-analysis-name CN5 \
             --include-quant-in-report-file-name true \
             --a5-calls true \
             --a5-calls-use-global false \
             --a5-feature-effects true \
             --a5-feature-effects-use-global true \
             --a5-residuals false \
             --a5-residuals-use-global false \
             --a5-sketch true \
             --a5-sketch-use-global true \
             --a5-summaries true \
             --a5-summaries-use-global false \
             --a5-write-models true \
             --a5-write-models-use-global true \
             --feat-effects false \
             --write-sketch false \
             --write-models false \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        execute("./apt-geno-qc \
             -v 3 \
             --spf-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.spf \
             --chrX-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --qca-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.qca \
             --qcc-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.qcc \
             --out-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/idata/cel/copynumber-workflow/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        execute("./apt-copynumber-reference \
             -v 3 \
             --genotype-report-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.report.txt \
             --reference-text-output false \
             --reference-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CNReference.a5.ref \
             --log2-input false \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --expr-summary-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.plier.summary.a5 \
             --genotype-calls-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.calls.a5 \
             --genotype-confidences-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.confidences.a5 \
             --out-dir " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-log2ratio \
             -v 3 \
             --log2-input false \
             --log2ratio-hdf5-output true \
             --log2ratio-text-output true \
             --text-output false \
             --cnchp-output false \
             --cychp-output false \
             --call-copynumber-engine true \
             --annotation-file ../../../regression-data/data/idata/lib/copynumber-workflow/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --reference-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CNReference.a5.ref \
             --genotype-report-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.report.txt \
             --expr-summary-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.plier.summary.a5 \
             --genotype-calls-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.calls.a5 \
             --genotype-confidences-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.confidences.a5 \
             --out-dir " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             --geno-qc-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.NA06985_GW6_C.CEL.a5 \
             -o " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.NA06993_GW6_C.CEL.a5 \
             -o " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             ");

         execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.NA07000_GW6_C.CEL.a5 \
             -o " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.NA07019_GW6_C.CEL.a5 \
             -o " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.NA07029_GW6_C.CEL.a5 \
             -o " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file " + testDir + "/qt-GenomeWideSNP_6/component_6_cels/CN5.NA07056_GW6_C.CEL.a5 \
             -o " + testDir + "/qt-GenomeWideSNP_6/component_6_cels \
             ");

        equivalency("component_6_cels", vSixFiles);
    }
};

int main(int argc, char* argv[])
{
  try {

    FsTestDir testDir;
    testDir.setTestDir("copynumber/copynumber-workflow-qt", true);

    ofstream logOut;
    string logName;
    logName = Fs::join(testDir.asString(), "qt-test-regression-copynumber-workflow.log");
    Verbose::out(1, "Log file: " + logName);
    Fs::mustOpenToWrite(logOut, logName.c_str());
    LogStream log(3, &logOut);
    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);

    Verbose::setLevel(3);
    test_regression_copynumber_workflow_quick test;
    test.testDir = testDir.asString();
    test.parseArgv(argv);
    //bool doValgrind = 
    test.doValgrind();
    string database = test.getDatabase();

    if ( !Fs::dirExists(test.testDir + "/qt-GenomeWideSNP_6") ) {
      Fs::mkdirPath(test.testDir + "/qt-GenomeWideSNP_6", false);
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

