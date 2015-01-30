////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
// This has to be last. Otherwise windows compile fails
#include "calvin_files/utils/src/Calvin.h"
#include "file5/File5.h"
#include "rtest/RT_Args.h"
#include "rtest/RT_Check.h"
#include "rtest/RT_Cmd.h"
#include "rtest/RT_Main.h"
#include "rtest/RT_Test.h"
#include "rtest/RT_Util.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>
//
using namespace std;

RT_Test* qt_single();
RT_Test* qt_batch();
RT_Test* single();
RT_Test* batch();

const char* cytogenetics_array_sixfiles[] = {"HapMap-As_NA18547_A08_01_NN_20081218", "HapMap-As_NA18603_B03_01_NN_20081218",
                         "HapMap-As_NA18943_C11_01_NN_20081218", "HapMap-As_NA18971_C02_01_NN_20081218",
                         "HapMap-Cs_NA10857_A10_01_NN_20081218", "HapMap-Cs_NA12003_D04_01_NN_20081218",
                         NULL
};

const char* qt_cytogenetics_array_sixfiles[] = {"HapMap-As_NA18547_A08_01_NN_20081218", "HapMap-As_NA18603_B03_01_NN_20081218",
                        "HapMap-As_NA18971_C02_01_NN_20081218",
                        "HapMap-Cs_NA10857_A10_01_NN_20081218", "HapMap-Cs_NA12003_D04_01_NN_20081218",
                         NULL
};

// Data Set columns to ignore (Specify as Group.Set.Column)
//ignoreHeaders.insert("AlgorithmData.MarkerABSignal.SCAR");

//Data Set column epsilon override (Specify as Group.Set.Column)
//mapEpsilon.insert(std::pair<std::string, float>("AlgorithmData.MarkerABSignal.SCAR", (float)0.0005));


///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if ( !Fs::dirExists("test-generated") ) {
    Fs::mkdir("test-generated");
  }

  RT_Main* rMain;
  rMain = new RT_Main("apt-rt-copynumber-cyto");

  //Store names of tests for Help Message
  //Store it in our RT_Main Instance
  std::vector<std::string> testNameList;
  testNameList.push_back("single");
  testNameList.push_back("batch");
  testNameList.push_back("qt_single");
  testNameList.push_back("qt_batch");
  rMain->setTestNameList(testNameList);

  //Set Variables
  rMain->setVar("sdktop", "../../..");
  rMain->setVar("gold_idata_copynumbercyto_testname", "${gold_dir}/idata/copynumber-cyto/${testname}");
  rMain->setVar("gold_idata_copynumbercyto_cytogeneticsarray_testname", "${gold_dir}/idata/copynumber-cyto/Cytogenetics_Array/${testname}");
  rMain->setVar("gold_idata_copynumbercyto_batch", "${gold_dir}/idata/copynumber-cyto/batch");
  rMain->setVar("gold_idata_copynumbercyto_cytogeneticsarray", "${gold_dir}/idata/copynumber-cyto/Cytogenetics_Array");
  rMain->setVar("gold_idata_lib_cytogeneticsarray", "${gold_dir}/idata/lib/Cytogenetics_Array");
  rMain->setVar("gold_idata_cel_cytogeneticsarray", "${gold_dir}/idata/cel/Cytogenetics_Array");
  rMain->setVar("out_cytogeneticsarray", "${out_dir}/Cytogenetics_Array");
  rMain->setVar("out_testname", "${out_dir}/${testname}");
  rMain->setVar("out_batch", "${out_dir}/batch");
  rMain->setVar("out_cytogeneticsarray_testname", "${out_dir}/Cytogenetics_Array/${testname}");

  rMain->setVar("gold_testname", "${gold_dir}/${testname}");
  rMain->setVar("gold_lib_cytogeneticsarray", "${gold_dir}/lib/Cytogenetics_Array");
  rMain->setVar("gold_cel_cytogeneticsarray", "${gold_dir}/cel/Cytogenetics_Array");
  rMain->setVar("gold_copynumbercyto", "${gold_dir}/copynumber-cyto");
  rMain->setVar("gold_copynumbercyto_testname", "${gold_copynumbercyto}/${testname}");
  rMain->setVar("gold_copynumbercyto_cytogeneticsarray_testname", "${gold_copynumbercyto}/Cytogenetics_Array/${testname}");
  rMain->setVar("gold_copynumbercyto_cytogeneticsarray", "${gold_copynumbercyto}/Cytogenetics_Array");
  rMain->setVar("gold_copynumbercyto_batch", "${gold_copynumbercyto}/batch");
  rMain->parseArgv(argc, argv);
  if(argc < 1)
    {
      cout<<"Error: Incorrect number of arguments given"<<endl;
      exit(1);
    }
  else if(rMain->getVarVal("allIntegration") != "")
    {
      rMain->addTest(qt_single());
      rMain->addTest(qt_batch());
    }
  else if(rMain->getVarVal("allNormal") != "")
    {
      rMain->addTest(single());
      rMain->addTest(batch());
    }
  //They want us to run all default tests
  else if(argc == 1 || rMain->getVarVal("all") != "")
    {
      //RUN ALL TESTS OMG
      cout<<"Running all default Regression Tests!"<<endl;
      rMain->addTest(qt_single());
      rMain->addTest(qt_batch());
      rMain->addTest(single());
      rMain->addTest(batch());
    }
  else
    {
      cout<<"Running Built-In Regression Tests"<<endl;

      if(rMain->getVarVal("qt_single") != "")
    {
      rMain->addTest(qt_single());
    }
      if(rMain->getVarVal("qt_batch") != "")
    {
      rMain->addTest(qt_batch());
    }
      if(rMain->getVarVal("single") != "")
    {
      rMain->addTest(single());
    }
      if(rMain->getVarVal("batch") != "")
    {
      rMain->addTest(batch());
    }
    }

  rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_single()
{
  RT_Test* rt;
  rt = new RT_Test("single");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_copynumber_cyto}");
  cmd->addArg("-v", "3");
  cmd->addArg("--cyto2", "true");
  cmd->addArg("--force", "true");
  cmd->addArg("--text-output", "true");
  cmd->addArg("--cnchp-output", "false");
  cmd->addArg("--cychp-output", "true");
  cmd->addArg("--set-analysis-name", "ca");
  cmd->addArg("--reference-input", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.na28.REF_MODEL");
  cmd->addArg("--snp-reference-input-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.na28.snpref.a5");
  cmd->addArg("--spf-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.spf");
  cmd->addArg("--chrX-probes", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.chrXprobes");
  cmd->addArg("--chrY-probes", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.chrYprobes");
  cmd->addArg("--annotation-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.na28.annot.db");
  cmd->addArg("--probe-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.probe_tab");
  cmd->addArg("--snp-qc-snp-list", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.snplist.txt");
  cmd->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-As_NA18547_A08_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-As_NA18603_B03_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-As_NA18971_C02_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-Cs_NA10857_A10_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-Cs_NA12003_D04_01_NN_20081218.CEL");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("--male-gender-ratio-cutoff=1.3");
  cmd->addArg("--female-gender-ratio-cutoff=1.0");
  cmd->addArg("--xx-cutoff=0.8");
  cmd->addArg("--xx-cutoff-high=1.07");
  cmd->addArg("--y-cutoff=0.65");
  cmd->addArg("--analysis", "log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
  cmd->addArg("--analysis", "kernel-smooth.sigma_span=50");
  cmd->addArg("--analysis", "cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6");
  cmd->addArg("--analysis", "cn-cyto2-gender.cutoff=0.5");
  cmd->addArg("--analysis", "cn-segment");
  cmd->addArg("--analysis", "loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
  cmd->addArg("--analysis", "mosaicism.marker-bandwidth=6000.confidence-window=251");
  cmd->addArg("--analysis", "allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05");
  cmd->addArg("--local-gc-background-correction-reference-method", "pdnn-reference-method");
  cmd->addArg("--wave-correction-reference-method", "wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false");
  cmd->addArg("--local-gc-background-intensity-adjustment-method", "pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
  cmd->addArg("--image-correction-intensity-adjustment-method", "high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001");
  cmd->addArg("--wave-correction-log2ratio-adjustment-method", "wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3");
  cmd->addArg("--reference-chromosome=22");

  std::vector<std::string> gen, gold;
  gen = Util::addPrefixSuffix(qt_cytogenetics_array_sixfiles, "${out_testname}/", ".ca.cychp");
  gold = Util::addPrefixSuffix(qt_cytogenetics_array_sixfiles, "${gold_idata_copynumbercyto_testname}/", ".ca.cychp");

RT_Check* check = rt->newCheck();
  check->setExe("${apt_calvin_equivalent}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  check->setGenGoldFileIsArg(true);
  check->addArg("--ignore-params-file", "./ignoreParamFile.tsv");
  check->addArg("--epsilon-map-file", "./qt_mapEpsilonFile.tsv");
  check->addArg("--epsilon", "0.0001");
  check->addArg("--correlation", "0.9999989");
  check->addArg("--check-header", "false");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* single()
{
  RT_Test* rt;
  rt = new RT_Test("single");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_copynumber_cyto}");
  cmd->addArg("-v", "3");
  cmd->addArg("--cyto2", "true");
  cmd->addArg("--force", "true");
  cmd->addArg("--text-output", "true");
  cmd->addArg("--cnchp-output", "false");
  cmd->addArg("--cychp-output", "true");
  cmd->addArg("--set-analysis-name", "ca");
  cmd->addArg("--reference-input", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.na28.REF_MODEL");
  cmd->addArg("--snp-reference-input-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.na28.snpref.a5");
  cmd->addArg("--cdf-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.CDF");
  cmd->addArg("--chrX-probes", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.chrXprobes");
  cmd->addArg("--chrY-probes", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.chrYprobes");
  cmd->addArg("--annotation-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.na28.annot.db");
  cmd->addArg("--probe-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.probe_tab");
  cmd->addArg("--snp-qc-snp-list", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.snplist.txt");
  cmd->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-As_NA18547_A08_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-As_NA18603_B03_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-As_NA18943_C11_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-As_NA18971_C02_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-Cs_NA10857_A10_01_NN_20081218.CEL");
  cmd->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-Cs_NA12003_D04_01_NN_20081218.CEL");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("--male-gender-ratio-cutoff=1.3");
  cmd->addArg("--female-gender-ratio-cutoff=1.0");
  cmd->addArg("--xx-cutoff=0.8");
  cmd->addArg("--xx-cutoff-high=1.07");
  cmd->addArg("--y-cutoff=0.65");
  cmd->addArg("--analysis", "log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
  cmd->addArg("--analysis", "kernel-smooth.sigma_span=50");
  cmd->addArg("--analysis", "cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6");
  cmd->addArg("--analysis", "cn-cyto2-gender.cutoff=0.5");
  cmd->addArg("--analysis", "cn-segment");
  cmd->addArg("--analysis", "loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
  cmd->addArg("--analysis", "mosaicism.marker-bandwidth=6000.confidence-window=251");
  cmd->addArg("--analysis", "allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05");
  cmd->addArg("--local-gc-background-correction-reference-method", "pdnn-reference-method");
  cmd->addArg("--wave-correction-reference-method", "wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false");
  cmd->addArg("--local-gc-background-intensity-adjustment-method", "pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
  cmd->addArg("--image-correction-intensity-adjustment-method", "high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001");
  cmd->addArg("--wave-correction-log2ratio-adjustment-method", "wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3");

  std::vector<std::string> gen, gold;
  gen = Util::addPrefixSuffix(cytogenetics_array_sixfiles, "${out_testname}/", ".ca.cychp");
  gold = Util::addPrefixSuffix(cytogenetics_array_sixfiles, "${gold_copynumbercyto_testname}/", ".ca.cychp");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_calvin_equivalent}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  check->setGenGoldFileIsArg(true);
  check->addArg("--ignore-params-file", "./ignoreParamFile.tsv");
  check->addArg("--epsilon", "0.0001");
  check->addArg("--correlation", "0.9999989");
  check->addArg("--check-header", "false");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_batch()
{
  RT_Test* rt;
  rt = new RT_Test("batch");
  RT_Test* rt1;
  rt1 = new RT_Test("batch-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("batch-part2");


  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_copynumber_cyto}");
  cmd1->addArg("-v", "3");
  cmd1->addArg("--cyto2", "true");
  cmd1->addArg("--force", "true");
  cmd1->addArg("--text-output", "true");
  cmd1->addArg("--cnchp-output", "false");
  cmd1->addArg("--cychp-output", "true");
  cmd1->addArg("--set-analysis-name", "ca");
  cmd1->addArg("--reference-output", "${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  cmd1->addArg("--snp-reference-input-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.na28.snpref.a5");
  cmd1->addArg("--snp-reference-output-file", "${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  cmd1->addArg("--spf-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.spf");
  cmd1->addArg("--chrX-probes", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.chrXprobes");
  cmd1->addArg("--chrY-probes", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.chrYprobes");
  cmd1->addArg("--annotation-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.na28.annot.db");
  cmd1->addArg("--probe-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.probe_tab");
  cmd1->addArg("--snp-qc-snp-list", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.snplist.txt");
  cmd1->addArg("--gender-override-file", "${gold_idata_lib_cytogeneticsarray}/HapMap-LOHReference-Full-Genders.txt");
  cmd1->addArg("--genotype-call-override-file", "${gold_idata_lib_cytogeneticsarray}/HapMapCalls_LOHReferenceSet-FullArray.txt");
  cmd1->addArg("--cel-files", "${gold_idata_cel_cytogeneticsarray}/CELFileList_5F5M.txt");
  cmd1->addArg("--out-dir", "${out_batch}");
  cmd1->addArg("--male-gender-ratio-cutoff=1.3");
  cmd1->addArg("--female-gender-ratio-cutoff=1.0");
  cmd1->addArg("--xx-cutoff=0.8");
  cmd1->addArg("--xx-cutoff-high=1.07");
  cmd1->addArg("--y-cutoff=0.65");
  cmd1->addArg("--analysis", "log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
  cmd1->addArg("--analysis", "kernel-smooth.sigma_span=50");
  cmd1->addArg("--analysis", "cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6");
  cmd1->addArg("--analysis", "cn-cyto2-gender.cutoff=0.5");
  cmd1->addArg("--analysis", "cn-segment");
  cmd1->addArg("--analysis", "loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
  cmd1->addArg("--analysis", "mosaicism.marker-bandwidth=6000.confidence-window=251");
  cmd1->addArg("--analysis", "allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05");
  cmd1->addArg("--local-gc-background-correction-reference-method", "pdnn-reference-method");
  cmd1->addArg("--wave-correction-reference-method", "wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false");
  cmd1->addArg("--local-gc-background-intensity-adjustment-method", "pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
  cmd1->addArg("--image-correction-intensity-adjustment-method", "high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001");
  cmd1->addArg("--wave-correction-log2ratio-adjustment-method", "wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3");
  cmd1->addArg("--reference-chromosome=22");

  std::vector<std::string> gen;
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");

  std::vector<std::string> gold;
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_idata_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");

  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_file5_equivalent}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->setGenGoldFileIsArg(true);
  check1->addArg("--ignore-columns-file", "./ignoreColumnFile.tsv");
  check1->addArg("--datasets-file", "./datasetFile.tsv");
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--correlation", "1.0");
  /*
  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_copynumber_cyto}");
  cmd2->addArg("-v", "3");
  cmd2->addArg("--cyto2", "true");
  cmd2->addArg("--force", "true");
  cmd2->addArg("--text-output", "true");
  cmd2->addArg("--cnchp-output", "false");
  cmd2->addArg("--cychp-output", "true");
  cmd2->addArg("--set-analysis-name", "ca");
  cmd2->addArg("--reference-input", "${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  cmd2->addArg("--snp-reference-input-file", "${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  cmd2->addArg("--spf-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.spf");
  cmd2->addArg("--chrX-probes", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.chrXprobes");
  cmd2->addArg("--chrY-probes", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.chrYprobes");
  cmd2->addArg("--annotation-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.na28.annot.db");
  cmd2->addArg("--probe-file", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.probe_tab");
  cmd2->addArg("--snp-qc-snp-list", "${gold_idata_lib_cytogeneticsarray}/Cytogenetics_Array.snplist.txt");
  cmd2->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-As_NA18547_A08_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-As_NA18603_B03_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-As_NA18971_C02_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-Cs_NA10857_A10_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_idata_cel_cytogeneticsarray}/HapMap-Cs_NA12003_D04_01_NN_20081218.CEL");
  cmd2->addArg("--out-dir", "${out_batch}");
  cmd2->addArg("--male-gender-ratio-cutoff=1.3");
  cmd2->addArg("--female-gender-ratio-cutoff=1.0");
  cmd2->addArg("--xx-cutoff=0.8");
  cmd2->addArg("--xx-cutoff-high=1.07");
  cmd2->addArg("--y-cutoff=0.65");
  cmd2->addArg("--analysis", "log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
  cmd2->addArg("--analysis", "kernel-smooth.sigma_span=50");
  cmd2->addArg("--analysis", "cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6");
  cmd2->addArg("--analysis", "cn-cyto2-gender.cutoff=0.5");
  cmd2->addArg("--analysis", "cn-segment");
  cmd2->addArg("--analysis", "loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
  cmd2->addArg("--analysis", "mosaicism.marker-bandwidth=6000.confidence-window=251");
  cmd2->addArg("--analysis", "allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05");
  cmd2->addArg("--local-gc-background-correction-reference-method", "pdnn-reference-method");
  cmd2->addArg("--wave-correction-reference-method", "wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false");
  cmd2->addArg("--local-gc-background-intensity-adjustment-method", "pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
  cmd2->addArg("--image-correction-intensity-adjustment-method", "high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001");
  cmd2->addArg("--wave-correction-log2ratio-adjustment-method", "wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3");
  cmd2->addArg("--reference-chromosome=22");

  std::vector<std::string> gen2,gold2;
  gen2 = Util::addPrefixSuffix(qt_cytogenetics_array_sixfiles, "${out_batch}/", ".ca.cychp");
  gold2 = Util::addPrefixSuffix(qt_cytogenetics_array_sixfiles, "${gold_idata_copynumbercyto_batch}/", ".ca.cychp");

  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_calvin_equivalent}");
  check2->setGenFiles(gen2);
  check2->setGoldFiles(gold2);
  check2->setGenGoldFileIsArg(true);
  check2->addArg("--epsilon-map-file", "./qt_mapEpsilonFile.tsv");
  check2->addArg("--ignore-params-file", "./ignoreParamFile.tsv");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--correlation", "0.9999989");
  check2->addArg("--check-header", "false");
  */
  rt->addSubtest(rt1);
  //rt->addSubtest(rt2);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* batch()
{
  RT_Test* rt;
  rt = new RT_Test("batch");
  RT_Test* rt1;
  rt1 = new RT_Test("batch-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("batch-part2");

  RT_Cmd* cmd = rt1->newCmd();
  cmd->setExe("${apt_copynumber_cyto}");
  cmd->addArg("-v", "3");
  cmd->addArg("--cyto2", "true");
  cmd->addArg("--force", "true");
  cmd->addArg("--text-output", "true");
  cmd->addArg("--cnchp-output", "false");
  cmd->addArg("--cychp-output", "true");
  cmd->addArg("--set-analysis-name", "ca");
  cmd->addArg("--reference-output", "${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  cmd->addArg("--snp-reference-input-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.na28.snpref.a5");
  cmd->addArg("--snp-reference-output-file", "${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  cmd->addArg("--cdf-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.CDF");
  cmd->addArg("--chrX-probes", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.chrXprobes");
  cmd->addArg("--chrY-probes", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.chrYprobes");
  cmd->addArg("--annotation-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.na28.annot.db");
  cmd->addArg("--probe-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.probe_tab");
  cmd->addArg("--snp-qc-snp-list", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.snplist.txt");
  cmd->addArg("--gender-override-file", "${gold_lib_cytogeneticsarray}/HapMap-LOHReference-Full-Genders.txt");
  cmd->addArg("--genotype-call-override-file", "${gold_lib_cytogeneticsarray}/HapMapCalls_LOHReferenceSet-FullArray.txt");
  cmd->addArg("--cel-files", "${gold_cel_cytogeneticsarray}/CELFileList_5F5M.txt");
  cmd->addArg("--out-dir", "${out_batch}");
  cmd->addArg("--male-gender-ratio-cutoff=1.3");
  cmd->addArg("--female-gender-ratio-cutoff=1.0");
  cmd->addArg("--xx-cutoff=0.8");
  cmd->addArg("--xx-cutoff-high=1.07");
  cmd->addArg("--y-cutoff=0.65");
  cmd->addArg("--analysis", "log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
  cmd->addArg("--analysis", "kernel-smooth.sigma_span=50");
  cmd->addArg("--analysis", "cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6");
  cmd->addArg("--analysis", "cn-cyto2-gender.cutoff=0.5");
  cmd->addArg("--analysis", "cn-segment");
  cmd->addArg("--analysis", "loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
  cmd->addArg("--analysis", "mosaicism.marker-bandwidth=6000.confidence-window=251");
  cmd->addArg("--analysis", "allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05");
  cmd->addArg("--local-gc-background-correction-reference-method", "pdnn-reference-method");
  cmd->addArg("--wave-correction-reference-method", "wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false");
  cmd->addArg("--local-gc-background-intensity-adjustment-method", "pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
  cmd->addArg("--image-correction-intensity-adjustment-method", "high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001");
  cmd->addArg("--wave-correction-log2ratio-adjustment-method", "wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3");

  std::vector<std::string> gen;
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gen.push_back("${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");

  std::vector<std::string> gold;
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  gold.push_back("${gold_copynumbercyto_batch}/Cytogenetics_Array.10_cels.REF_MODEL");

  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_file5_equivalent}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->setGenGoldFileIsArg(true);
  check1->addArg("--ignore-columns-file", "./ignoreColumnFile.tsv");
  check1->addArg("--datasets-file", "./datasetFile.tsv");
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--correlation", "1.0");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_copynumber_cyto}");
  cmd2->addArg("-v", "3");
  cmd2->addArg("--cyto2", "true");
  cmd2->addArg("--force", "true");
  cmd2->addArg("--text-output", "true");
  cmd2->addArg("--cnchp-output", "false");
  cmd2->addArg("--cychp-output", "true");
  cmd2->addArg("--set-analysis-name", "ca");
  cmd2->addArg("--reference-input", "${out_batch}/Cytogenetics_Array.10_cels.REF_MODEL");
  cmd2->addArg("--snp-reference-input-file", "${out_batch}/Cytogenetics_Array.10_cels.snpref.a5");
  cmd2->addArg("--cdf-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.CDF");
  cmd2->addArg("--chrX-probes", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.chrXprobes");
  cmd2->addArg("--chrY-probes", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.chrYprobes");
  cmd2->addArg("--annotation-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.na28.annot.db");
  cmd2->addArg("--probe-file", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.probe_tab");
  cmd2->addArg("--snp-qc-snp-list", "${gold_lib_cytogeneticsarray}/Cytogenetics_Array.snplist.txt");
  cmd2->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-As_NA18547_A08_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-As_NA18603_B03_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-As_NA18943_C11_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-As_NA18971_C02_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-Cs_NA10857_A10_01_NN_20081218.CEL");
  cmd2->addArg("--cels", "${gold_cel_cytogeneticsarray}/HapMap-Cs_NA12003_D04_01_NN_20081218.CEL");
  cmd2->addArg("--out-dir", "${out_batch}");
  cmd2->addArg("--male-gender-ratio-cutoff=1.3");
  cmd2->addArg("--female-gender-ratio-cutoff=1.0");
  cmd2->addArg("--xx-cutoff=0.8");
  cmd2->addArg("--xx-cutoff-high=1.07");
  cmd2->addArg("--y-cutoff=0.65");
  cmd2->addArg("--analysis", "log2-ratio-cyto2.gc-correction=true.median-autosome-median-normalization=true.median-smooth-marker-count=5.trim-high=2.0.trim-low=-2.5");
  cmd2->addArg("--analysis", "kernel-smooth.sigma_span=50");
  cmd2->addArg("--analysis", "cn-cyto2.diagonal-weight=0.995.mapd-weight=0.22.min-segment-size=5.hmm-confidence-weight=0.6");
  cmd2->addArg("--analysis", "cn-cyto2-gender.cutoff=0.5");
  cmd2->addArg("--analysis", "cn-segment");
  cmd2->addArg("--analysis", "loh-cyto2.lohCNSegSeparation=1000000.minInformation=100.lambdaCritical=8.0");
  cmd2->addArg("--analysis", "mosaicism.marker-bandwidth=6000.confidence-window=251");
  cmd2->addArg("--analysis", "allele-peaks.step=100.window=250.point-count=128.bandwidth=0.45.cutoff=0.05");
  cmd2->addArg("--local-gc-background-correction-reference-method", "pdnn-reference-method");
  cmd2->addArg("--wave-correction-reference-method", "wave-correction-reference-method.trim=2.0.percentile=0.75.demean=false");
  cmd2->addArg("--local-gc-background-intensity-adjustment-method", "pdnn-intensity-adjustment-method.predicted-intensity-bin-count=20.gc-bin-count=20.residual-trim=2.0");
  cmd2->addArg("--image-correction-intensity-adjustment-method", "high-pass-filter-intensity-adjustment-method.data-block-rows=320.data-block-cols=2015.mini-block-rows=8.mini-block-cols=8.global-smooth-weight=256.local-smooth-weight=64.converged=0.0001");
  cmd2->addArg("--wave-correction-log2ratio-adjustment-method", "wave-correction-log2ratio-adjustment-method.bandwidth=101.bin-count=25.wave-count=3");

  std::vector<std::string> gen2,gold2;
  gen2 = Util::addPrefixSuffix(cytogenetics_array_sixfiles, "${out_batch}/", ".ca.cychp");
  gold2 = Util::addPrefixSuffix(cytogenetics_array_sixfiles, "${gold_copynumbercyto_batch}/", ".ca.cychp");

  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_calvin_equivalent}");
  check2->setGenFiles(gen2);
  check2->setGoldFiles(gold2);
  check2->setGenGoldFileIsArg(true);
  check2->addArg("--ignore-params-file", "./ignoreParamFile.tsv");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--correlation", "0.9999989");
  check2->addArg("--check-header", "false");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;
}
