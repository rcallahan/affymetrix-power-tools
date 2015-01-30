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
 * @file   test-regression-copynumber-workflow.cpp
 * @brief  Program for doing regression tests on probeset-genotype data.
 */

// This has to be last. Otherwise windows compile fails
#include "calvin_files/utils/src/Calvin.h"

#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/CnchpCheck.cpp"
#include "copynumber/apt-copynumber-cyto/regression/File5EquivalentCheck.h"
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

class test_regression_copynumber_workflow :  public RegressionSuite
{
public:
    int m_numPassed;
    int m_numFailed;
    std::string m_testDir;
    std::set<std::string> m_setIgnore;
    std::vector<std::string> m_vSixFiles;
    std::vector<std::string> m_vOneFile;
    std::vector<std::string> m_vSNP6CytoAberrations;
    std::vector<std::string> m_vSNP6CytoAberrationsInternal;


test_regression_copynumber_workflow()
{
        m_numPassed = 0;
        m_numFailed = 0;

        // Header paramters to ignore
        m_setIgnore.insert("FileCreationTime");
        m_setIgnore.insert("FileIdentifier");
        m_setIgnore.insert("program-version");
        m_setIgnore.insert("create_date");
        m_setIgnore.insert("create-date");
        m_setIgnore.insert("affymetrix-algorithm-param-option-verbose");
        m_setIgnore.insert("affymetrix-algorithm-param-option-exec-guid");
        m_setIgnore.insert("affymetrix-algorithm-param-option-program-cvs-id");
        m_setIgnore.insert("affymetrix-algorithm-param-option-version-to-report");
        m_setIgnore.insert("affymetrix-algorithm-param-option-command-line");
        m_setIgnore.insert("affymetrix-algorithm-param-option-mem-usage");
        m_setIgnore.insert("affymetrix-algorithm-param-option-run-probeset-genotype");
        m_setIgnore.insert("affymetrix-algorithm-param-option-cels");
        m_setIgnore.insert("affymetrix-algorithm-param-option-out-dir");
        m_setIgnore.insert("affymetrix-algorithm-param-state-time-start");
        m_setIgnore.insert("affymetrix-algorithm-param-state-free-mem-at-start");

        // Data Set columns to ignore (Specify as Group.Set.Column)
        m_setIgnore.insert("MultiData.CopyNumber.SmoothSignal");
        m_setIgnore.insert("Segments.CN.SegmentID");
        m_setIgnore.insert("Segments.LOH.SegmentID");
        m_setIgnore.insert("Segments.CNNeutralLOH.SegmentID");
        m_setIgnore.insert("Segments.NormalDiploid.SegmentID");
        m_setIgnore.insert("Segments.Mosaicism.SegmentID");
        m_setIgnore.insert("Segments.NoCall.SegmentID");

        // File tags to check for equivalency
        m_vSixFiles.push_back("NA06985_GW6_C");
        m_vSixFiles.push_back("NA06993_GW6_C");
        m_vSixFiles.push_back("NA07000_GW6_C");
        m_vSixFiles.push_back("NA07019_GW6_C");
        m_vSixFiles.push_back("NA07029_GW6_C");
        m_vSixFiles.push_back("NA07056_GW6_C");

        m_vOneFile.push_back("NA07029_GW6_C");

        m_vSNP6CytoAberrations.push_back("080207_LC_F2_U141_BETA6.cn5");
        m_vSNP6CytoAberrations.push_back("20100218_CC4_CC4_GenomeWideSNP_6.cn5");
        m_vSNP6CytoAberrations.push_back("Beta1.cn5");
        m_vSNP6CytoAberrations.push_back("Beta2.cn5");
        m_vSNP6CytoAberrations.push_back("Beta3.cn5");
        m_vSNP6CytoAberrations.push_back("Beta4.cn5");
        m_vSNP6CytoAberrations.push_back("Beta7.cn5");
        m_vSNP6CytoAberrations.push_back("CO_LL060080-6.0-E2.cn5");
        m_vSNP6CytoAberrations.push_back("CO_LL070093-6.0-D2.cn5");
        m_vSNP6CytoAberrations.push_back("CO_LL070202-6.0-G1.cn5");
        m_vSNP6CytoAberrations.push_back("NA02269_GW6-R_CYTOKIN_A7.cn5");
        m_vSNP6CytoAberrations.push_back("NA02521_GW6-R_CYTOKIN_A8.cn5");
        m_vSNP6CytoAberrations.push_back("NA03091_GW6-R_CYTOKIN_A12.cn5");
        m_vSNP6CytoAberrations.push_back("NA03535_GW6-R_CYTOKIN_B3.cn5");
        m_vSNP6CytoAberrations.push_back("NA04126_GW6-R_CYTOKIN_B6.cn5");
        m_vSNP6CytoAberrations.push_back("NA07691_GW6-R_CYTOKIN_B10.cn5");
        m_vSNP6CytoAberrations.push_back("NA09024_GW6-R_CYTOKIN_B12.cn5");
        m_vSNP6CytoAberrations.push_back("NA09209_GW6-R_CYTOKIN_C1.cn5");
        m_vSNP6CytoAberrations.push_back("NA09888_GW6-R_CYTOKIN_C2.cn5");
        m_vSNP6CytoAberrations.push_back("NA11033_GW6-R_CYTOKIN_C6.cn5");
        m_vSNP6CytoAberrations.push_back("NA11167_GW6-R_CYTOKIN_C7.cn5");
        m_vSNP6CytoAberrations.push_back("NA11382_GW6-R_CYTOKIN_C8.cn5");
        m_vSNP6CytoAberrations.push_back("NA11385_GW6-R_CYTOKIN_C9.cn5");
        m_vSNP6CytoAberrations.push_back("NA11404_GW6-R_CYTOKIN_C10.cn5");
        m_vSNP6CytoAberrations.push_back("NA11419_GW6-R_CYTOKIN_C11.cn5");
        m_vSNP6CytoAberrations.push_back("NA11515_GW6-R_CYTOKIN_C12.cn5");
        m_vSNP6CytoAberrations.push_back("NA12016_GW6-R_CYTOKIN_D1.cn5");
        m_vSNP6CytoAberrations.push_back("NA13164_GW6-R_CYTOKIN_D2.cn5");
        m_vSNP6CytoAberrations.push_back("NA13415_GW6-R_CYTOKIN_D4.cn5");
        m_vSNP6CytoAberrations.push_back("NA13419_GW6-R_CYTOKIN_D5.cn5");
        m_vSNP6CytoAberrations.push_back("NA13476_GW6-R_CYTOKIN_D7.cn5");
        m_vSNP6CytoAberrations.push_back("NA13553_GW6-R_CYTOKIN_D8.cn5");
        m_vSNP6CytoAberrations.push_back("NA14124_GW6-R_CYTOKIN_D9.cn5");
        m_vSNP6CytoAberrations.push_back("NA14129_GW6-R_CYTOKIN_D10.cn5");
        m_vSNP6CytoAberrations.push_back("NA16447_GW6-R_CYTOKIN_D11.cn5");
        m_vSNP6CytoAberrations.push_back("NA16453_GW6-R_CYTOKIN_D12.cn5");
        m_vSNP6CytoAberrations.push_back("NA16455_GW6-R_CYTOKIN_E1.cn5");
        m_vSNP6CytoAberrations.push_back("NA16718_GW6-R_CYTOKIN_E2.cn5");
        m_vSNP6CytoAberrations.push_back("NA18319_GW6-R_CYTOKIN_E3.cn5");
        m_vSNP6CytoAberrations.push_back("UC_2_R2271.cn5");

        m_vSNP6CytoAberrationsInternal.push_back("09-0161EH.cn5");
        m_vSNP6CytoAberrationsInternal.push_back("3_10-0622MW_05-05-10.cn5");
        m_vSNP6CytoAberrationsInternal.push_back("5_09-0161EH_02-27-09.cn5");
}



void run()
{

#if ( (!defined(_WIN32)) || (defined(_WIN64)))

        reference_6_cels_na32();
        copynumber_1_cels_na32();

        reference_6_cels_na31();
        copynumber_1_cels_na31();

        reference_6_cels_na30();
        copynumber_1_cels_na30();

        reference_6_cels();
        copynumber_1_cels();

        copynumber_HapMap270_1_cels();
        copynumber_HapMap270_6_cels();

        snp6CytoAberrations();

        reference_6_cels_wave_correction();
        copynumber_1_cels_wave_correction();

        
#else
        Verbose::out(1,"APT-994: test-regression-copynumber-workflow.cpp commented out for Windows 32 bit.");
#endif        
        
}

void snp6CytoAberrations()
{
        Verbose::out(1, "Working on snp6CytoAberrations.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/snp6CytoAberrations";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/snp6CytoAberrations/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/snp6CytoAberrations";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        string command = "./apt-copynumber-workflow "
          " --out-dir " + newDirectory +
          "  --xml-file-append-only ../../../regression-data/data/lib/GenomeWideSNP_6/jobFile7 ";


        checks.push_back(new CnchpCheck(        m_vSNP6CytoAberrations,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));

        goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/snp6CytoAberrations/INTERNAL";


        checks.push_back(new CnchpCheck(        m_vSNP6CytoAberrationsInternal,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));





        RegressionTest test("snp6CytoAberrations",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: snp6CytoAberrations " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: snp6CytoAberrations " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }

}


void copynumber_HapMap270_1_cels()
{
        Verbose::out(1, "Working on copynumber_HapMap_1_cels.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/copynumber_HapMap270_1_cels";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/copynumber_HapMap270_1_cels/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/copynumber_HapMap270_1_cels";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + outdir + " \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ";


        checks.push_back(new CnchpCheck(        m_vOneFile,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));


        RegressionTest test("copynumber_HapMap270_1_cels",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: copynumber_HapMap270_1_cels " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: copynumber_HapMap270_1_cels " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }

}


void copynumber_HapMap270_6_cels()
{

        Verbose::out(1, "Working on copynumber_HapMap_6_cels.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/copynumber_HapMap270_6_cels";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/copynumber_HapMap270_6_cels/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/copynumber_HapMap270_6_cels";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;


        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + outdir + " \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ";
        checks.push_back(new CnchpCheck(        m_vSixFiles,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));


        RegressionTest test("copynumber_HapMap270_6_cels",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: copynumber_HapMap270_6_cels " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: copynumber_HapMap270_6_cels " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}


void reference_6_cels_na32()
{

        Verbose::out(1, "Working on reference_6_cels_na32.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/reference_6_cels_na32";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/reference_6_cels_na32/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/reference_6_cels_na32";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        Fs::rm( m_testDir + "/GenomeWideSNP_6/reference_6_cels_na32/CNReference.a5.ref", false);
        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-output " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na32/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/release/GenomeWideSNP_6.na32.annot.db \
             --o " + outdir + " \
             --temp-dir " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na32 \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ";
        checks.push_back(new CnchpCheck(        m_vSixFiles,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));
        std::string strFileName1 = "../../../regression-data/data/copynumber/GenomeWideSNP_6/reference_6_cels_na32/CNReference.a5.ref";
        std::string strFileName2 =  m_testDir + "/GenomeWideSNP_6/reference_6_cels_na32/CNReference.a5.ref";

        std::set<std::string> setLocalIgnore;
        float epsilon = 0.0000001;
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CopyNumber", "Reference", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "CN5.plier.feature-response", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "CN5.snp-posteriors", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "target-sketch", setLocalIgnore, epsilon, 1.0));


        RegressionTest test("reference_6_cels_na32",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: reference_6_cels_na32." + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: reference_6_cels_na32 " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}


void reference_6_cels_na31()
{

        Verbose::out(1, "Working on reference_6_cels_na31.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/reference_6_cels_na31";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/reference_6_cels_na31/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/reference_6_cels_na31";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        Fs::rm( m_testDir + "/GenomeWideSNP_6/reference_6_cels_na31/CNReference.a5.ref", false);
        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-output " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na31/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/release/GenomeWideSNP_6.na31.annot.db \
             --o " + outdir + " \
             --temp-dir " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na31 \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ";
        checks.push_back(new CnchpCheck(        m_vSixFiles,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));
        std::string strFileName1 = "../../../regression-data/data/copynumber/GenomeWideSNP_6/reference_6_cels_na31/CNReference.a5.ref";
        std::string strFileName2 =  m_testDir + "/GenomeWideSNP_6/reference_6_cels_na31/CNReference.a5.ref";

        std::set<std::string> setLocalIgnore;
        float epsilon = 0.0000001;
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CopyNumber", "Reference", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "CN5.plier.feature-response", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "CN5.snp-posteriors", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "target-sketch", setLocalIgnore, epsilon, 1.0));


        RegressionTest test("reference_6_cels_na31",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: reference_6_cels_na31." + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: reference_6_cels_na31 " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}


void reference_6_cels_na30()
{

        Verbose::out(1, "Working on reference_6_cels_na30.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/reference_6_cels_na30";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/reference_6_cels_na30/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/reference_6_cels_na30";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        Fs::rm( m_testDir + "/GenomeWideSNP_6/reference_6_cels_na30/CNReference.a5.ref", false);
        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-output " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na30/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na30.annot.db \
             --o " + outdir + " \
             --temp-dir " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na30 \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ";
        checks.push_back(new CnchpCheck(        m_vSixFiles,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));
        std::string strFileName1 = "../../../regression-data/data/copynumber/GenomeWideSNP_6/reference_6_cels_na30/CNReference.a5.ref";
        std::string strFileName2 =  m_testDir + "/GenomeWideSNP_6/reference_6_cels_na30/CNReference.a5.ref";

        std::set<std::string> setLocalIgnore;
        float epsilon = 0.0000001;
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CopyNumber", "Reference", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "CN5.plier.feature-response", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "CN5.snp-posteriors", setLocalIgnore, epsilon, 1.0));
        checks.push_back(new File5EquivalentCheck(strFileName1, strFileName2, "CN5", "target-sketch", setLocalIgnore, epsilon, 1.0));


        RegressionTest test("reference_6_cels_na30",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: reference_6_cels_na30 " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: reference_6_cels_na30 " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}

void reference_6_cels()
{

        Verbose::out(1, "Working on reference_6_cels.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/reference_6_cels";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/reference_6_cels/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/reference_6_cels";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        Fs::rm( m_testDir + "/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref", false);
        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-output " + m_testDir + "/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + outdir + " \
             --temp-dir " + m_testDir + "/GenomeWideSNP_6/reference_6_cels \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ";
        checks.push_back(new CnchpCheck(        m_vSixFiles,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));


        RegressionTest test("reference_6_cels",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: reference_6_cels " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: reference_6_cels " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}

 
void reference_6_cels_wave_correction()
{

        Verbose::out(1, "Working on reference_6_cels_wave_correction.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/r6wc";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/r6wc/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/r6wc";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        Fs::rm( m_testDir + "/GenomeWideSNP_6/r6wc/CNReference.a5.ref", false);

        string command = "./apt-copynumber-workflow \
                         --verbose 1                  \
                         --run-geno-qc false          \
                         --cychp-output false                           \
                         --cnchp-output true                            \
                         --text-output true                             \
                         --reference-output " + m_testDir + "/GenomeWideSNP_6/r6wc/CNReference.a5.ref \
                         --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
                         --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
                         --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
                         --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
                         --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
                         --o " + outdir + "                             \
                         --temp-dir " + m_testDir + "/GenomeWideSNP_6/r6wc \
                         --wave-correction-reference-method wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false.wave-count=6 \
                         ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
                         ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
                         ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
                         ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
                         ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
                         ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
                         ";
        checks.push_back(new CnchpCheck(        m_vSixFiles,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));




        RegressionTest test("reference_6_cels_wave_correction",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: reference_6_cels_wave_correction " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: reference_6_cels_wave_correction " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }

}


void copynumber_6_cels()
{

        Verbose::out(1, "Working on copynumber_6_cels.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/copynumber_6_cels";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/copynumber_6_cels/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/copynumber_6_cels";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;


        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input " + m_testDir + "/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + m_testDir + "/GenomeWideSNP_6/copynumber_6_cels \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ";

        checks.push_back(new CnchpCheck(        m_vSixFiles,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));

        RegressionTest test("copynumber_6_cels",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: copynumber_6_cels " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: copynumber_6_cels " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }

}


void copynumber_1_cels_wave_correction()
{
        Verbose::out(1, "Working on copynumber_1_cels_wave_correction.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/c1wc";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/c1wc/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/c1wc";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;


        string command = "./apt-copynumber-workflow \
                         --verbose 1 \
                         --run-geno-qc false \
                         --cychp-output false \
                         --cnchp-output true \
                         --text-output true \
                         --reference-input " + m_testDir + "/GenomeWideSNP_6/r6wc/CNReference.a5.ref \
                         --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
                         --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
                         --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
                         --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
                         --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
                         --o " + outdir + "                             \
                         --analysis log2-ratio \
                         --wave-correction-log2ratio-adjustment-method wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=4 \
                         ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
                         ";

        checks.push_back(new CnchpCheck(        m_vOneFile,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));



        RegressionTest test("copynumber_1_cels_wave_correction",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: copynumber_1_cels_wave_correction " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: copynumber_1_cels_wave_correction " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }

}


void copynumber_1_cels_na32()
{

        Verbose::out(1, "Working on copynumber_1_cels_na32.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/copynumber_1_cels_na32";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/copynumber_1_cels_na32/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/copynumber_1_cels_na32";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na32/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/release/GenomeWideSNP_6.na32.annot.db \
             --o " + outdir + " \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ";

        checks.push_back(new CnchpCheck(        m_vOneFile,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));



        RegressionTest test("copynumber_1_cels_na32",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: copynumber_1_cels_na32 " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: copynumber_1_cels_na32." + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}



void copynumber_1_cels_na31()
{

        Verbose::out(1, "Working on copynumber_1_cels_na31.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/copynumber_1_cels_na31";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/copynumber_1_cels_na31/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/copynumber_1_cels_na31";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na31/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/release/GenomeWideSNP_6.na31.annot.db \
             --o " + outdir + " \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ";

        checks.push_back(new CnchpCheck(        m_vOneFile,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));



        RegressionTest test("copynumber_1_cels_na31",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: copynumber_1_cels_na31 " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: copynumber_1_cels_na31." + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}




void copynumber_1_cels_na30()
{

        Verbose::out(1, "Working on copynumber_1_cels_na30.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/copynumber_1_cels_na30";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/copynumber_1_cels_na30/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/copynumber_1_cels_na30";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input " + m_testDir + "/GenomeWideSNP_6/reference_6_cels_na30/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na30.annot.db \
             --o " + outdir + " \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ";

        checks.push_back(new CnchpCheck(        m_vOneFile,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));



        RegressionTest test("copynumber_1_cels_na30",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: copynumber_1_cels_na30 " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: copynumber_1_cels_na30 " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}


void copynumber_1_cels()
{

        Verbose::out(1, "Working on copynumber_1_cels.");

        std::string  goldDirectory="../../../regression-data/data/copynumber/GenomeWideSNP_6/copynumber_1_cels";
        std::string  newDirectory= m_testDir + "/GenomeWideSNP_6/copynumber_1_cels/";
        std::string  infix="CN5";

        std::set<std::string> setSetIgnore;
        setSetIgnore.clear();
        std::map<std::string, float> mapEpsilon;

        string outdir =  m_testDir + "/GenomeWideSNP_6/copynumber_1_cels";
        Fs::rmdirPath(outdir, false);
        Fs::mkdir(outdir, false);
        vector<RegressionCheck *> checks;

        string command = "./apt-copynumber-workflow \
             --verbose 1 \
             --run-geno-qc false \
             --cychp-output false \
             --cnchp-output true \
             --text-output true \
             --reference-input " + m_testDir + "/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o " + outdir + " \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ";

        checks.push_back(new CnchpCheck(        m_vOneFile,
                                                goldDirectory,
                                                newDirectory,
                                                infix,
                                                m_setIgnore,
                                                setSetIgnore,
                                                mapEpsilon,
                                                0.0000001,
                                                1.0,
                                                false));



        RegressionTest test("copynumber_1_cels",command.c_str(), checks, false);
        test.setSuite(*this,  outdir, outdir + "/apt-copynumber-workflow.log", outdir + "/valgrind.log");
        if(!test.pass()) {
            Verbose::out(1, "Error in: copynumber_1_cels " + test.getErrorMsg());
        }
        else
        {
            Verbose::out(1, "Passed: copynumber_1_cels " + test.getErrorMsg());
        }

        {
            m_numPassed += test.getNumPassed();
            m_numFailed += test.getNumFailed();
        }
}


int getNumPassed()
{
    return m_numPassed;
}

int getNumFailed()
{
    return m_numFailed;
}

};


int main(int argc, char* argv[])
{
  try {
    double byteSize = (10.0 * 1073741824);
    FsTestDir testDir;
    testDir.setTestDir("copynumber/workflow", true, byteSize);

    ofstream logOut;
    string logName;
    logName =  Fs::join(testDir.asString() ,"test-regression-copynumber-workflow.log");
    Verbose::out(1, "Log file: " + logName);
    Fs::mustOpenToWrite(logOut, logName.c_str());
    LogStream log(3, &logOut);
    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);

    Verbose::setLevel(3);
    test_regression_copynumber_workflow test;
    test.m_testDir = testDir.asString();
    test.parseArgv(argv);
    //bool doValgrind = 
    test.doValgrind();
    string database = test.getDatabase();

    if ( !Fs::dirExists(test.m_testDir + "/GenomeWideSNP_6") ) {
      Fs::mkdirPath(test.m_testDir + "/GenomeWideSNP_6", false);
    }

    Verbose::out(1, "Execute regression test: test-regression-copynumber-workflow");

    test.run();

    Verbose::out(1, "NumPassed: " + ToStr(test.getNumPassed()) + " NumFailed: " + ToStr(test.getNumFailed()) + " for test-regression-copynumber-workflow");
    logOut.close();
    return test.getNumFailed() != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

/*   The following are regression testcases created by Walt that have not been run since at least Jan 10.  They have not been translated
 *   to the standard style used in all chipstream/copynumber regression testcases.


    void component_6_cels()
    {
        execute("./apt-probeset-genotype \
             -v 1 \
             --all-types true \
             --set-gender-method cn-probe-chrXY-ratio \
             --a5-global-file test-generated/GenomeWideSNP_6/component_6_cels/CNReference.a5 \
             --analysis adapter-type-norm.StyOnlySnpRatio=0.8287.NspOnlySnpRatio=0.7960.NspOnlyCnRatio=1.4218.NspBothSnpRatio=0.9954.NspBothCnRatio=1.6392,quant-norm,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.wobble=0.05.MS=1 \
             --qmethod-spec plier.optmethod=1.FixFeatureEffect=true \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --out-dir test-generated/GenomeWideSNP_6/component_6_cels \
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
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        execute("./apt-geno-qc \
             -v 3 \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --qca-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.qca \
             --qcc-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.qcc \
             --out-file test-generated/GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        execute("./apt-copynumber-reference \
             -v 3 \
             --genotype-report-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.report.txt \
             --reference-text-output false \
             --reference-file test-generated/GenomeWideSNP_6/component_6_cels/CNReference.a5.ref \
             --log2-input false \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --expr-summary-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.plier.summary.a5 \
             --genotype-calls-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.calls.a5 \
             --genotype-confidences-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.confidences.a5 \
             --out-dir test-generated/GenomeWideSNP_6/component_6_cels \
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
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --reference-file test-generated/GenomeWideSNP_6/component_6_cels/CNReference.a5.ref \
             --genotype-report-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.report.txt \
             --expr-summary-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.plier.summary.a5 \
             --genotype-calls-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.calls.a5 \
             --genotype-confidences-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.confidences.a5 \
             --out-dir test-generated/GenomeWideSNP_6/component_6_cels \
             --geno-qc-file test-generated/GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file test-generated/GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.NA06985_GW6_C.CEL.a5 \
             -o test-generated/GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file test-generated/GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.NA06993_GW6_C.CEL.a5 \
             -o test-generated/GenomeWideSNP_6/component_6_cels \
             ");

         execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file test-generated/GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.NA07000_GW6_C.CEL.a5 \
             -o test-generated/GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file test-generated/GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.NA07019_GW6_C.CEL.a5 \
             -o test-generated/GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file test-generated/GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.NA07029_GW6_C.CEL.a5 \
             -o test-generated/GenomeWideSNP_6/component_6_cels \
             ");

        execute("./apt-copynumber-analysis \
             -v 3 \
             --geno-qc-file test-generated/GenomeWideSNP_6/component_6_cels/apt-geno-qc.txt \
             --cnchp-output true \
             --cychp-output true \
             --text-output true \
             --log2ratio-file test-generated/GenomeWideSNP_6/component_6_cels/CN5.NA07056_GW6_C.CEL.a5 \
             -o test-generated/GenomeWideSNP_6/component_6_cels \
             ");

        equivalency("component_6_cels", vSixFiles);
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
             --reference-input test-generated/GenomeWideSNP_6/reference_6_cels_med_norm/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o test-generated/GenomeWideSNP_6/copynumber_1_cels_med_norm \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
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
             --reference-input test-generated/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o test-generated/GenomeWideSNP_6/copynumber_6_cels_gc_correction_false \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        equivalency("copynumber_6_cels_gc_correction_false", vSixFiles);
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
             --reference-input test-generated/GenomeWideSNP_6/reference_6_cels/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o test-generated/GenomeWideSNP_6/copynumber_6_cels \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

        equivalency("copynumber_6_cels", vSixFiles);
    }

    void reference_6_cels_med_norm()
    {
        Fs::rm("test-generated/GenomeWideSNP_6/reference_6_cels_med_norm/CNReference.a5.ref");
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --normalization-type 2 \
             --run-geno-qc false \
             --cychp-output true \
             --cnchp-output true \
             --text-output true \
             --reference-output test-generated/GenomeWideSNP_6/reference_6_cels_med_norm/CNReference.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o test-generated/GenomeWideSNP_6/reference_6_cels_med_norm \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06985_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA06993_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07000_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07019_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07056_GW6_C.CEL \
             ");

//        equivalency("reference_6_cels_med_norm", vSixFiles); // Reference building
    }


    void copynumber_HapMap270_1_cels_probeset_ids()
    {
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --probeset-ids ../../../regression-data/data/lib/GenomeWideSNP_6/ProbeSetRestrictList.txt \
             --run-geno-qc false \
             --cychp-output true \
             --cnchp-output true \
             --text-output true \
             --reference-input ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o test-generated/GenomeWideSNP_6/copynumber_HapMap270_1_cels_probeset_ids \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ");

        equivalency("copynumber_HapMap270_1_cels_probeset_ids", vOneFile);
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
             --qca-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.qca \
             --qcc-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.qcc \
             --reference-input ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o test-generated/GenomeWideSNP_6/copynumber_HapMap270_for_Carl \
             --cel-files ../../../../../../cyto2/cyto2_cel-files/upd_trio/CELFileList.txt \
             ");
    }

void copynumber_HapMap270_1_cels_cancer()
{
        execute("./apt-copynumber-workflow \
             --verbose 3 \
             --run-geno-qc false \
             --cychp-output true \
             --cnchp-output true \
             --text-output true \
             --reference-input ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.hapmap270.na26.1.r2.a5.ref \
             --cdf-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.cdf \
             --chrX-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrXprobes \
             --chrY-probes ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.chrYprobes \
             --special-snps ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.specialSNPs \
             --annotation-file ../../../regression-data/data/lib/GenomeWideSNP_6/GenomeWideSNP_6.na27.6.annot.db \
             --o test-generated/GenomeWideSNP_6/copynumber_HapMap270_1_cels_cancer \
             ../../../regression-data/data/cel/GenomeWideSNP_6/NA07029_GW6_C.CEL \
             ");

        equivalency("copynumber_HapMap270_1_cels_cancer", vOneFile);
}

*/












