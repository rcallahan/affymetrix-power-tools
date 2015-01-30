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
#include "rtest/RT_Args.h"
#include "rtest/RT_Check.h"
#include "rtest/RT_Cmd.h"
#include "rtest/RT_Main.h"
#include "rtest/RT_Test.h"
#include "rtest/RT_Util.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
//
using namespace std;

RT_Test* doRefBuild();
RT_Test* doRefRun();
RT_Test* doPlasmidRun();

string m_celsSuffix = ".CEL";
const char* m_cels[]={
  "ZQ_C01_NA19193",    "ZQ_C02_NA12004",    "ZQ_C03_NA19194",
  "ZQ_C04_NA18856",    "ZQ_C05_NA19137",    "ZQ_C06_NA19144",
  "ZQ_C07_NA12003",    "ZQ_C08_NA19143",    "ZQ_C09_NA10838",
  "ZQ_C10_NA19130",    "ZQ_C11_NA18855",    "ZQ_C12_NA19139",
  "ZQ_D01_NA18856",    "ZQ_D02_NA19194",    "ZQ_D03_NA19138",
  "ZQ_D04_NA19131",    "ZQ_D05_NA10838",    "ZQ_D06_NA19192",
  "ZQ_D07_NA19145",    "ZQ_D08_NA19132",    "ZQ_D09_NA18857",
  "ZQ_D10_gCtrl_1",    "ZQ_D11_gCtrl_2",    "ZQ_D12_gCtrl_3",
  NULL
};

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if( !Fs::dirExists("test-generated" ) ){
    Fs::mkdir("test-generated", false);
  }

  RT_Main* rMain = new RT_Main("apt-rt-dmet-genotype");

 //Store names of tests for Help Message
 //Store it in our RT_Main Instance
 std::vector<std::string> testNameList;
 testNameList.push_back("doRefBuild");
 testNameList.push_back("doRefRun");
 testNameList.push_back("doPlasmidRun");
 rMain->setTestNameList(testNameList);

 //Set Variables
 rMain->setVar("sdktop", "../..");
 rMain->setVar("celsuffix", ".CEL");
 rMain->setVar("gold_lib_dmetplus", "${gold_dir}/lib/DMET_Plus");
 rMain->setVar("gold_dmetgenotype_dmetplus", "${gold_dir}/dmet-genotype/DMET_Plus");
 rMain->setVar("gold_dmetgenotype_dmetplus_testname_chp", "${gold_dmetgenotype_dmetplus}/${testname}/chp");
 rMain->setVar("out_testname_chp", "${out_testname}/chp");
 rMain->setVar("dmetchpsuffix", ".dmet.chp");

 if(argc < 1)
   {
     cout<<"Error: Incorrect number of arguments given"<<endl;
     exit(1);
   }
 //They want us to run all default tests
 else if(argc == 1)
   {
     rMain->addTest(doRefBuild());
     rMain->addTest(doRefRun());
     rMain->addTest(doPlasmidRun());
   }
 else
   {
     rMain->parseArgv(argc, argv);

     cout<<"Running Built-In Regression Tests"<<endl;
     if(rMain->getVarVal("doRefBuild") != "")
       {
	 rMain->addTest(doRefBuild());
       }
     if(rMain->getVarVal("doRefRun") != "")
       {
	 rMain->addTest(doRefRun());
       }
     if(rMain->getVarVal("doPlasmidRun") != "")
       {
	 rMain->addTest(doPlasmidRun());
       }

   }
 rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRefBuild()
{
  RT_Test* rt;
  rt = new RT_Test("doRefBuild");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_dmet_genotype}");
  cmd->addArg("--batch-info=true");
  cmd->addArg("--null-context=false");
  cmd->addArg("--reference-input=${gold_lib_dmetplus}/DMET_Plus.r3.genomic.ref.a5");
  cmd->addArg("--out-dir=${out_testname}");
  cmd->addArg("--cel-files=${gold_dmetgenotype_dmetplus}/cels.small.txt");
  cmd->addArg("--cdf-file=${gold_lib_dmetplus}/DMET_Plus.r3.cdf");
  cmd->addArg("--chrX-probes=${gold_lib_dmetplus}/DMET_Plus.r3.chrXprobes");
  cmd->addArg("--chrY-probes=${gold_lib_dmetplus}/DMET_Plus.r3.chrYprobes");
  cmd->addArg("--batch-name=dynamic-batch-1");
  cmd->addArg("--reference-output=${out_testname}/reference.a5");
  cmd->addArg("--region-model=${gold_lib_dmetplus}/DMET_Plus.r3.cn-region-models.txt");
  cmd->addArg("--probeset-model=${gold_lib_dmetplus}/DMET_Plus.r3.cn-probeset-models.txt");
  cmd->addArg("--probeset-ids=${gold_lib_dmetplus}/DMET_Plus.r3.genomic.gt.ps");
  cmd->addArg("--probeset-ids-reported=${gold_dmetgenotype_dmetplus}/consent1.txt");
  cmd->addArg("--cn-region-gt-probeset-file=${gold_lib_dmetplus}/DMET_Plus.r3.cn-gt.ps");
  cmd->addArg("--cc-chp-output=true");
  cmd->addArg("--");
  cmd->addArg("--summaries=true");
  cmd->addArg("--feat-effects");
  cmd->addArg("--");
  cmd->addArg("--text-output");
  cmd->addArg("--");
  cmd->addArg("--table-output");
  cmd->addArg("--feat-effects");
  cmd->addArg("--summaries");

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(m_cels, "${gold_dmetgenotype_dmetplus_testname_chp}/", "${dmetchpsuffix}");
  gen = Util::addPrefixSuffix(m_cels, "${out_testname_chp}/", "${dmetchpsuffix}");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  /*
  //Calvin Chp Check
  rt->addCheckMult("${apt_check_calvinchp}", gen, gold);
  */
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRefRun()
{
  RT_Test* rt;
  rt = new RT_Test("doRefRun");
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_dmet_genotype}");
  cmd->addArg("--batch-info=false");
  cmd->addArg("--null-context=true");
  cmd->addArg("--reference-input=${gold_lib_dmetplus}/DMET_Plus.r3.genomic.ref.a5");
  cmd->addArg("--out-dir=${out_testname}");
  cmd->addArg("--cel-files=${gold_dmetgenotype_dmetplus}/cels.small.txt");
  cmd->addArg("--spf-file=${gold_lib_dmetplus}/DMET_Plus.r3.spf");
  cmd->addArg("--chrX-probes=${gold_lib_dmetplus}/DMET_Plus.r3.chrXprobes");
  cmd->addArg("--chrY-probes=${gold_lib_dmetplus}/DMET_Plus.r3.chrYprobes");
  cmd->addArg("--region-model=${gold_lib_dmetplus}/DMET_Plus.r3.cn-region-models.txt");
  cmd->addArg("--probeset-model=${gold_lib_dmetplus}/DMET_Plus.r3.cn-probeset-models.txt");
  cmd->addArg("--probeset-ids=${gold_lib_dmetplus}/DMET_Plus.r3.genomic.gt.ps");
  cmd->addArg("--probeset-ids-reported=${gold_dmetgenotype_dmetplus}/consent2.txt");
  cmd->addArg("--cn-region-gt-probeset-file=${gold_lib_dmetplus}/DMET_Plus.r3.cn-gt.ps");
  cmd->addArg("--cc-chp-output=true");
  cmd->addArg("--geno-call-thresh=0.01");
  cmd->addArg("--");
  cmd->addArg("--summaries=true");
  cmd->addArg("--feat-effects");
  cmd->addArg("--");
  cmd->addArg("--text-output");
  cmd->addArg("--");
  cmd->addArg("--table-output");
  cmd->addArg("--feat-effects");
  cmd->addArg("--summaries");

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(m_cels, "${gold_dmetgenotype_dmetplus_testname_chp}/", "${dmetchpsuffix}");
  gen = Util::addPrefixSuffix(m_cels, "${out_testname_chp}/", "${dmetchpsuffix}");

  //Calvin Chp Check
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);

  // rt->addCheckMult("${apt_check_calvinchp}", gen, gold);

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doPlasmidRun()
{
  RT_Test* rt;
  rt = new RT_Test("doPlasmidRun");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_dmet_genotype}");
  cmd->addArg("--reference-input=${gold_lib_dmetplus}/DMET_Plus.r3.plasmid.ref.a5");
  cmd->addArg("--out-dir=${out_testname}");
  cmd->addArg("--cel-files=${gold_dmetgenotype_dmetplus}/cels.small.txt");
  cmd->addArg("--cdf-file=${gold_lib_dmetplus}/DMET_Plus.r3.cdf");
  cmd->addArg("--chrX-probes=${gold_lib_dmetplus}/DMET_Plus.r3.chrXprobes");
  cmd->addArg("--chrY-probes=${gold_lib_dmetplus}/DMET_Plus.r3.chrYprobes");
  cmd->addArg("--region-model=${gold_lib_dmetplus}/DMET_Plus.r3.cn-region-models.txt");
  cmd->addArg("--probeset-ids=${gold_lib_dmetplus}/DMET_Plus.r3.plasmid.gt.ps");
  cmd->addArg("--probeset-ids-reported=${gold_lib_dmetplus}/DMET_Plus.r3.plasmid.gt.ps");
  cmd->addArg("--cc-chp-output=true");
  cmd->addArg("--geno-call-thresh=0.2");
  cmd->addArg("--run-cn-engine=false");
  //cmd->addArg("--run-cn-engine=false");

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(m_cels, "${gold_dmetgenotype_dmetplus_testname_chp}/", "${dmetchpsuffix}");
  gen = Util::addPrefixSuffix(m_cels, "${out_testname_chp}/", "${dmetchpsuffix}");
  
  //Calvin Chp Check
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);

  //  rt->addCheckMult("${apt_check_calvinchp}", gen, gold);

  return rt;
}
