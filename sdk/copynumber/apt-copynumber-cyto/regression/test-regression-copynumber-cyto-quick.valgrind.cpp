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
 * @file   test-regression-copynumber-cyto-quick.cpp
 * @brief  Program for doing regression tests on cyto data.
 */

#include "calvin_files/utils/src/Calvin.h"
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

string valgrindString1 = "valgrind --leak-check=yes --log-file=";
string valgrindString2 = ".valgrind.txt ";




class test_regression_copynumber_cyto_quick
{
public:
    int numPassed;
    int numFailed;
    std::string testDir;
    std::set<std::string> setIgnore;
    std::map<std::string, float> mapEpsilon;
    std::vector<std::string> vCytogenetics_Array_SixFiles;

    test_regression_copynumber_cyto_quick()
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
        mapEpsilon.insert(std::pair<string, float>("Segments.CN.MeanMarkerDistance", 376.0f));
        // File tags to check for equivalency
        vCytogenetics_Array_SixFiles.push_back("HapMap-As_NA18547_A08_01_NN_20081218");
        vCytogenetics_Array_SixFiles.push_back("HapMap-As_NA18603_B03_01_NN_20081218");
        vCytogenetics_Array_SixFiles.push_back("HapMap-As_NA18971_C02_01_NN_20081218");
        vCytogenetics_Array_SixFiles.push_back("HapMap-Cs_NA10857_A10_01_NN_20081218");
        vCytogenetics_Array_SixFiles.push_back("HapMap-Cs_NA12003_D04_01_NN_20081218");
    }

    void execute(const std::string& str) {
      std::string cmd = Fs::convertCommandToUnc(str);
      if (system(cmd.c_str()) != 0) {Err::errAbort("Execution failed for command: " + cmd);}
    }

    void equivalency(const std::string& name, const std::string& modifier, std::vector<string>& vFileTags)
    {
        bool bPassed = true;
        AffxString strFileName1;
        AffxString strFileName2;
        Verbose::out(1, "*");
        for(unsigned int i = 0; (i < vFileTags.size()); i++)
        {

            strFileName1 = "../../../regression-data/data/idata/copynumber-cyto/" + name + "/" + vFileTags[i] + "." + modifier + ".cychp";
            strFileName2 = testDir + "/" + name + "/" + vFileTags[i] + "." + modifier + ".cychp";
            if (!Calvin::equivalent(strFileName1, strFileName2, setIgnore, mapEpsilon, 0.0001, 0.9999989, false)) {bPassed = false;}
        }
        Verbose::out(1, "*");
        if (bPassed) {
                  Verbose::out(1, "***Passed*** " + name + " vs " + modifier);
                  numPassed++;
                }
                else {
                  Verbose::out(1, "***Failed*** " + name + " vs " + modifier);
                  numFailed++;
                }
    }

    void run()
    {
        single();
        batch();
    }

    void single()
    {
          string name = "qt-single";
          Fs::rmdirPath(testDir + "/single/", false);

                 execute( valgrindString1 + name + valgrindString2 + " ./apt-copynumber-cyto \
                    -v 3 --cyto2 true \
                    --force true \
                    --text-output true \
                    --cnchp-output false \
                    --cychp-output true \
                    --set-analysis-name ca \
                    --reference-input ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.na28.REF_MODEL \
                    --snp-reference-input-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.na28.snpref.a5 \
                    --spf-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.spf \
                    --chrX-probes ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.chrXprobes \
                    --chrY-probes ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.chrYprobes \
                    --annotation-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.na28.annot.db \
                    --probe-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.probe_tab \
                    --snp-qc-snp-list ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.snplist.txt \
                    --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-As_NA18547_A08_01_NN_20081218.CEL \
                    --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-As_NA18603_B03_01_NN_20081218.CEL \
                    --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-As_NA18971_C02_01_NN_20081218.CEL \
                    --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-Cs_NA10857_A10_01_NN_20081218.CEL \
                    --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-Cs_NA12003_D04_01_NN_20081218.CEL \
                    --out-dir " + testDir + "/single \
                    --male-gender-ratio-cutoff=1.3 \
                    --female-gender-ratio-cutoff=1.0 \
                    --xx-cutoff=0.8 \
                    --xx-cutoff-high=1.07 \
                    --y-cutoff=0.65 \
                    --analysis log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5 \
                    --analysis kernel-smooth.sigma_span=50 \
                    --analysis cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6 \
                    --analysis cn-cyto2-gender.cutoff=0.5 \
                    --analysis cn-segment \
                    --analysis loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0 \
                    --analysis mosaicism.marker-bandwidth=6000.confidence-window=251 \
                    --analysis allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05 \
                    --local-gc-background-correction-reference-method pdnn-reference-method \
                    --wave-correction-reference-method wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false \
                    --local-gc-background-intensity-adjustment-method pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0 \
                    --image-correction-intensity-adjustment-method high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001 \
                    --wave-correction-log2ratio-adjustment-method wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3 --reference-chromosome=22 \
            ");
        equivalency("single", "ca", vCytogenetics_Array_SixFiles);
    }

    void batch()
    {
          string name1 = "qt-batch-part1";
          string name2 = "qt-batch-part2";
          Fs::rmdirPath(testDir + "/batch/", false);

                 execute( valgrindString1 + name1 + valgrindString2 + " ./apt-copynumber-cyto \
                    -v 3 --cyto2 true \
                    --force true \
                    --text-output true \
                    --cnchp-output false \
                    --cychp-output true \
                    --set-analysis-name ca \
                    --reference-output " + testDir + "/batch/Cytogenetics_Array.10_cels.REF_MODEL \
                    --snp-reference-input-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.na28.snpref.a5 \
                    --snp-reference-output-file " + testDir + "/batch/Cytogenetics_Array.10_cels.snpref.a5 \
                    --spf-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.spf \
                    --chrX-probes ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.chrXprobes \
                    --chrY-probes ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.chrYprobes \
                    --annotation-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.na28.annot.db \
                    --probe-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.probe_tab \
                    --snp-qc-snp-list ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.snplist.txt \
                    --gender-override-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/HapMap-LOHReference-Full-Genders.txt \
                    --genotype-call-override-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/HapMapCalls_LOHReferenceSet-FullArray.txt \
                    --cel-files ../../../regression-data/data/idata/cel/Cytogenetics_Array/CELFileList_5F5M.txt \
                    --out-dir " + testDir + "/batch \
                    --male-gender-ratio-cutoff=1.3 \
                    --female-gender-ratio-cutoff=1.0 \
                    --xx-cutoff=0.8 \
                    --xx-cutoff-high=1.07 \
                    --y-cutoff=0.65 \
                    --analysis log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5 \
                    --analysis kernel-smooth.sigma_span=50 \
                    --analysis cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6 \
                    --analysis cn-cyto2-gender.cutoff=0.5 \
                    --analysis cn-segment \
                    --analysis loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0 \
                    --analysis mosaicism.marker-bandwidth=6000.confidence-window=251 \
                    --analysis allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05 \
                    --local-gc-background-correction-reference-method pdnn-reference-method \
                    --wave-correction-reference-method wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false \
                    --local-gc-background-intensity-adjustment-method pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0 \
                    --image-correction-intensity-adjustment-method high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001 \
                    --wave-correction-log2ratio-adjustment-method wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3 --reference-chromosome=22 \
            ");

        bool bPassed = true;
        std::set<std::string> setIgnore;
        std::string strFileName1 = "../../../regression-data/data/idata/copynumber-cyto/batch/Cytogenetics_Array.10_cels.snpref.a5";
        std::string strFileName2 = testDir + "/batch/Cytogenetics_Array.10_cels.snpref.a5";
//        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "SNPReference", "Parameters", setIgnore, 0.0001, 1.0)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "SNPReference", "EMSNParameters", setIgnore, 0.0001, 1.0)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "SNPReference", "EMSNArray", setIgnore, 0.0001, 1.0)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "SNPReference", "SNPReference", setIgnore, 0.0003, 1.0)) {bPassed = false;}

        strFileName1 = "../../../regression-data/data/idata/copynumber-cyto/batch/Cytogenetics_Array.10_cels.REF_MODEL";
        strFileName2 = testDir + "/batch/Cytogenetics_Array.10_cels.REF_MODEL";
//        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "Parameters", setIgnore, 0.0001, 1.0)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "AntigenomicProbes", setIgnore, 0.0001, 1.0)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "MedianSignals", setIgnore, 0.5, 0.999999)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "ProbeEffects", setIgnore, 0.0001, 0.999999)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "SNPPosteriors", setIgnore, 0.0001, 1.0)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "SketchCN", setIgnore, 0.0001, 1.0)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "SketchSNP", setIgnore, 0.0001, 1.0)) {bPassed = false;}
        if (!affx::File5_File::equivalent(strFileName1, strFileName2, "Cyto2", "WaveCorrection", setIgnore, 0.0001, 0.999999)) {bPassed = false;}

        if (bPassed)
        {
                 execute( valgrindString1 + name2 + valgrindString2 + " ./apt-copynumber-cyto \
                        -v 3 --cyto2 true \
                        --force true \
                        --text-output true \
                        --cnchp-output false \
                        --cychp-output true \
                        --set-analysis-name ca \
                        --reference-input " + testDir + "/batch/Cytogenetics_Array.10_cels.REF_MODEL \
                        --snp-reference-input-file " + testDir + "/batch/Cytogenetics_Array.10_cels.snpref.a5 \
                        --spf-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.spf \
                        --chrX-probes ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.chrXprobes \
                        --chrY-probes ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.chrYprobes \
                        --annotation-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.na28.annot.db \
                        --probe-file ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.probe_tab \
                        --snp-qc-snp-list ../../../regression-data/data/idata/lib/Cytogenetics_Array/Cytogenetics_Array.snplist.txt \
                        --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-As_NA18547_A08_01_NN_20081218.CEL \
                        --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-As_NA18603_B03_01_NN_20081218.CEL \
                        --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-As_NA18971_C02_01_NN_20081218.CEL \
                        --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-Cs_NA10857_A10_01_NN_20081218.CEL \
                        --cels ../../../regression-data/data/idata/cel/Cytogenetics_Array/HapMap-Cs_NA12003_D04_01_NN_20081218.CEL \
                        --out-dir " + testDir + "/batch \
                        --male-gender-ratio-cutoff=1.3 \
                        --female-gender-ratio-cutoff=1.0 \
                        --xx-cutoff=0.8 \
                        --xx-cutoff-high=1.07 \
                        --y-cutoff=0.65 \
                        --analysis log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5 \
                        --analysis kernel-smooth.sigma_span=50 \
                        --analysis cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6 \
                        --analysis cn-cyto2-gender.cutoff=0.5 \
                        --analysis cn-segment \
                        --analysis loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0 \
                        --analysis mosaicism.marker-bandwidth=6000.confidence-window=251 \
                        --analysis allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05 \
                        --local-gc-background-correction-reference-method pdnn-reference-method \
                        --wave-correction-reference-method wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false \
                        --local-gc-background-intensity-adjustment-method pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0 \
                        --image-correction-intensity-adjustment-method high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001 \
                        --wave-correction-log2ratio-adjustment-method wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3 --reference-chromosome=22 \
                ");

            equivalency("batch", "ca", vCytogenetics_Array_SixFiles);
        }
        else {numFailed++;}
    }
};

int main(int argc, char* argv[]) {
  try {
    FsTestDir testDir;
    testDir.setTestDir("copynumber/copynumber-cyto-valgrind",true);
    
    /* Set up the logging and message handlers. */

    ofstream logOut;
    string logName;
    logName = testDir.asString() + "/test-regression-copynumber-cyto-quick.log";
    valgrindString1 = valgrindString1 + testDir.asString() + "/";
    
    Verbose::out(1, "Log file: " + logName);
    Fs::mustOpenToWrite(logOut, logName.c_str());
    LogStream log(3, &logOut);

    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);

    Verbose::setLevel(3);
    test_regression_copynumber_cyto_quick test;
    test.testDir = testDir.asString();

    Verbose::out(1, "Execute regression test: test-regression-copynumber-cyto-valgrind");

    test.run();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed) + " for test-regression-copynumber-cyto-valgrind");
    logOut.close();
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  return 1;
}

