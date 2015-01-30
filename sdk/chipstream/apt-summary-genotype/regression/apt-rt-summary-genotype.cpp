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

RT_Test* qt_doBrlmmpSnp5Two();
RT_Test* qt_doBrlmmpSnp5();
RT_Test* doBrlmmpSnp5Two();
RT_Test* doBrlmmpSnp5();

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  RT_Main* rMain;
  rMain = new RT_Main("apt-rt-summary-genotype");

  if ( !Fs::dirExists("test-generated") ) {
    Fs::mkdir("test-generated", false);
  }

  std::vector<std::string> testNameList;
  testNameList.push_back("qt_doBrlmmpSnp5Two");
  testNameList.push_back("qt_doBrlmmpSnp5");
  testNameList.push_back("doBrlmmpSnp5Two");
  testNameList.push_back("doBrlmmpSnp5");
  rMain->setTestNameList(testNameList);

  rMain->setVar("gold_lib_genomewidesnp5", "${gold_dir}/lib/GenomeWideSNP_5");
  rMain->setVar("gold_summarygenotype_genomewidesnp5", "${gold_dir}/chipstream/summary-genotype/GenomeWideSNP_5");
  rMain->setVar("gold_summarygenotype_genomewidesnp5_testname", "${gold_dir}/chipstream/summary-genotype/GenomeWideSNP_5/${testname}");
  rMain->setVar("gold_idata_sumgeno_genomewidesnp5", "${gold_dir}/idata/sum-geno/GenomeWideSNP_5");
  rMain->setVar("gold_idata_sumgeno_genomewidesnp5_testname", "${gold_dir}/idata/sum-geno/GenomeWideSNP_5/${testname}");
  rMain->setVar("gold_idata_lib_genomewidesnp5", "${gold_dir}/idata/lib/GenomeWideSNP_5");

  rMain->parseArgv(argc, argv);

  if(argc < 1)
    {
      cout<<"Error: Incorrect number of arguments given"<<endl;
      exit(1);
    }
  else if(rMain->getVarVal("allNormal") != "")
    {
      rMain->addTest(doBrlmmpSnp5Two());
      rMain->addTest(doBrlmmpSnp5());
    }
  else if(rMain->getVarVal("allIntegration") != "")
    {
      rMain->addTest(qt_doBrlmmpSnp5Two());
      rMain->addTest(qt_doBrlmmpSnp5());
    }
  //They want us to run all default tests
  else if(argc == 1 || rMain->getVarVal("all") != "")
    {
      rMain->addTest(qt_doBrlmmpSnp5Two());
      rMain->addTest(qt_doBrlmmpSnp5());
      rMain->addTest(doBrlmmpSnp5Two());
      rMain->addTest(doBrlmmpSnp5());
    }
  else
    {
      cout<<"Running Built-In Regression Tests"<<endl;

      if(rMain->getVarVal("qt_doBrlmmpSnp5Two") != "")
	{
	  rMain->addTest(qt_doBrlmmpSnp5Two());
	}
      if(rMain->getVarVal("qt_doBrlmmpSnp5") != "")
	{
	  rMain->addTest(qt_doBrlmmpSnp5());
	}

      if(rMain->getVarVal("doBrlmmpSnp5Two") != "")
	{
	  rMain->addTest(doBrlmmpSnp5Two());
	}
      if(rMain->getVarVal("doBrlmmpSnp5") != "")
	{
	  rMain->addTest(doBrlmmpSnp5());
	}
    }
  rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBrlmmpSnp5Two()
{
  RT_Test* rt;
  rt = new RT_Test("qt-doBrlmmpSnp5Two");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_summary_genotype}");
  cmd->addArg("-p", "CM=1.bins=100.K=2.SB=0.003.MS=0.05");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--summaries-a=${gold_idata_sumgeno_genomewidesnp5}/alleleA.plier.summary.txt");
  cmd->addArg("--summaries-b=${gold_idata_sumgeno_genomewidesnp5}/alleleB.plier.summary.txt");
  cmd->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd->addArg("--read-genders", "${gold_idata_sumgeno_genomewidesnp5}/brlmm-p.genders.txt");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_idata_sumgeno_genomewidesnp5_testname}/brlmm-p.calls.txt"); 
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");
  check1->dump();

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_idata_sumgeno_genomewidesnp5_testname}/brlmm-p.confidences.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");
  check2->dump();
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBrlmmpSnp5Two()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5Two");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_summary_genotype}");
  cmd->addArg("-p", "CM=1.bins=100.K=2.SB=0.003.MS=0.05");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--summaries-a=${gold_summarygenotype_genomewidesnp5}/alleleA.plier.summary.txt");
  cmd->addArg("--summaries-b=${gold_summarygenotype_genomewidesnp5}/alleleB.plier.summary.txt");
  cmd->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd->addArg("--read-genders", "${gold_summarygenotype_genomewidesnp5}/brlmm-p.genders.txt");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_summarygenotype_genomewidesnp5_testname}/brlmm-p.calls.txt"); 
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_summarygenotype_genomewidesnp5_testname}/brlmm-p.confidences.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBrlmmpSnp5()
{
  RT_Test* rt;
  rt = new RT_Test("qt-doBrlmmpSnp5");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_summary_genotype}");
  cmd->addArg("-p", "CM=1.bins=100.K=2.SB=0.003.MS=0.05");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("-s", "${gold_idata_sumgeno_genomewidesnp5}/brlmm-p.plier.summary.txt");
  cmd->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd->addArg("--read-genders", "${gold_idata_sumgeno_genomewidesnp5}/brlmm-p.genders.txt");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_idata_sumgeno_genomewidesnp5_testname}/brlmm-p.calls.txt"); 
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_idata_sumgeno_genomewidesnp5_testname}/brlmm-p.confidences.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBrlmmpSnp5()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_summary_genotype}");
  cmd->addArg("-p", "CM=1.bins=100.K=2.SB=0.003.MS=0.05");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("-s", "${gold_summarygenotype_genomewidesnp5}/brlmm-p.plier.summary.txt");
  cmd->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd->addArg("--read-genders", "${gold_summarygenotype_genomewidesnp5}/brlmm-p.genders.txt");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_summarygenotype_genomewidesnp5_testname}/brlmm-p.calls.txt"); 
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_summarygenotype_genomewidesnp5_testname}/brlmm-p.confidences.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}
