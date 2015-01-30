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
#include "rtest/RT_Check.h"
#include "rtest/RT_Cmd.h"
#include "rtest/RT_Main.h"
#include "rtest/RT_Test.h"
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

//Array containing the names example CHP files
const char *chpFiles[] = {
  "NA06985_GW6_C", "NA06993_GW6_C", "NA07000_GW6_C", 
  "NA07019_GW6_C", "NA07029_GW6_C", "NA07056_GW6_C", 
  "NA07345_GW6_C", "NA10830_GW6_C", "NA10831_GW6_C", 
  "NA10839_GW6_C", "NA10855_GW6_C", "NA10857_GW6_C", 
  "NA10860_GW6_C", "NA10863_GW6_C", "NA11881_GW6_C", 
  "NA11882_GW6_C", "NA11993_GW6_C", "NA12003_GW6_C", 
  "NA12004_GW6_C",
  NULL
};

//Function Forward Declaration
RT_Test* defineDoChpFiles();
RT_Test* defineDoChpFilesDefault();

///////////////////////////////////////////////////////////////////////////////////////////////////
//Main Function
int main(int argc, char* argv[])
{

  ///Create Output Directories
  if ( !Fs::dirExists("test-generated") ) {
    Fs::mkdir("test-generated", false);
  }

  //Create new instance or RT_Main that will contain our suite of tests
  RT_Main* rMain = new RT_Main("apt-rt-canary");

  //Store names of tests for Help Message
  //Store it in our RT_Main Instance
  std::vector<std::string> testNameList;
  testNameList.push_back("doChpFiles");
  testNameList.push_back("doChpFilesDefault");
  rMain->setTestNameList(testNameList);

  ///Set Variables for this test environment
  rMain->setVar("sdktop", "../..");
  rMain->setVar("chpfilessuffix", ".canary-v1.cnvchp"); 
  rMain->setVar("gold_lib_genomewidesnp6", "${gold_dir}/lib/GenomeWideSNP_6");
  rMain->setVar("gold_canary_genomewidesnp6", "${gold_dir}/canary/GenomeWideSNP_6");
  rMain->setVar("gold_canary_genomewidesnp6_testname", "${gold_canary_genomewidesnp6}/${testname}");

  //If no Arguments are Given
  ///They want us to run all default tests
  if(argc == 1)
    {
      cout<<"Running All Default Regression Tests"<<endl;
      rMain->addTest(defineDoChpFiles());
      rMain->addTest(defineDoChpFilesDefault());
    }
  //Parse the argv and run specified tests/options
  else
    {
      rMain->parseArgv(argc, argv);

      cout<<"Running Built-In Regression Tests"<<endl;
     
      if(rMain->getVarVal("doChpFiles") != "")
	{
	  rMain->addTest(defineDoChpFiles());
	}
      if(rMain->getVarVal("doChpFilesDefault") != "")
	{
	  rMain->addTest(defineDoChpFilesDefault());
	}
    }

  //Run All tests that have been built up into the framework
  rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//Creates new test called doChpFiles
RT_Test* defineDoChpFiles()
{
  //Create new RT_Test instance
  RT_Test* rt =  new RT_Test("doChpFiles");

  //Create new command that will be stored in our Test
  RT_Cmd* cmd = rt->newCmd();
  //Set the Executable
  cmd->setExe("${apt_canary}");
  //Add command line arguments for this command
  cmd->addArg("--apt-summarize-analysis", "quant-norm.target=1000,pm-only,plier.optmethod=1,expr.genotype=true");
  cmd->addArg("--apt-canary-analysis", "af-weight=1.0,TOL=1e-3,hwe_tol=1e-4,hwe_tol2=1e-11,fraction-giveaway-0=0.20,fraction-giveaway-1=0.10,fraction-giveaway-2=0.15,fraction-giveaway-3=0.20,fraction-giveaway-4=0.30,min-fill-prop=0.10,conf-interval-half-width=1.959963984540054,inflation=1.3,min-cluster-variance=0.0001,pseudopoint-factor=100");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--cnv-region-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.canary-v1.region");
  cmd->addArg("--cnv-normalization-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.canary-v1.normalization");
  cmd->addArg("--cnv-prior-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.canary-v1.prior");
  cmd->addArg("--cnv-map-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.canary-v1.bed");
  cmd->addArg("--analysis-name", "canary-v1");
  cmd->addArg("--cel-files", "${gold_canary_genomewidesnp6}/fas_cel_files.txt");
  cmd->addArg("--cc-chp-output", "true");
  cmd->addArg("--table-output", "true");
  
  //Produce a vector of generated and gold files by combining our array of files with a given prefix and suffix
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_canary_genomewidesnp6_testname}/", "${chpfilessuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/", "${chpfilessuffix}");

  //Create new Calvin Chp Checks
  //Creates multiple checks from the given gen and gold vectors
  //Breaks them down into individual pairs of files
  //length of gen must equal length of gold
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  //rt->addCheckMult("${apt_check_calvinchp}", gen, gold);

  //Create new Matrix Check
  RT_Check* check1 = rt->newCheck();
  //Set the exe
  check1->setExe("${apt_check_matrix}");
  //Add command line Arguments
  check1->addArg("--gen", "${out_testname}/canary-v1.allele_average_summary.txt");
  check1->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.allele_average_summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  //Create new Matrix Check
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/canary-v1.calls.txt");
  check2->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.calls.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  //Create new Matrix Check
  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/canary-v1.cnv_region_summary.txt");
  check3->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.cnv_region_summary.txt"); 
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  //Create new Matrix Check
  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/canary-v1.confidences.txt");
  check4->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.confidences.txt"); 
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //Create new Matrix Check
  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_testname}/canary-v1.summary.txt");
  check5->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.summary.txt"); 
  check5->addArg("--epsilon", "0.0001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");

  //Return a pointer to the test just created
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* defineDoChpFilesDefault()
{
  //Create new RT_Test instance
  RT_Test* rt= new RT_Test("doChpFilesDefault");

  //Create new command that will be stored in our Test
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_canary}");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--cnv-region-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.canary-v1.region");
  cmd->addArg("--cnv-normalization-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.canary-v1.normalization");
  cmd->addArg("--cnv-prior-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.canary-v1.prior");
  cmd->addArg("--cnv-map-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.canary-v1.bed");
  cmd->addArg("--analysis-name", "canary-v1");
  cmd->addArg("--cel-files", "${gold_canary_genomewidesnp6}/fas_cel_files.txt");
  cmd->addArg("--cc-chp-output", "true");
  cmd->addArg("--table-output", "true");

  //Produce a vector of generated and gold files by combining our array of files with a given prefix and suffix
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_canary_genomewidesnp6_testname}/", "${chpfilessuffix}");
  gen =  Util::addPrefixSuffix(chpFiles, "${out_testname}/", "${chpfilessuffix}");

  //Create new Calvin Chp Checks
  //Creates multiple checks from the given gen and gold vectors
  //Breaks them down into individual pairs of files
  //length of gen must equal length of gold
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  //rt->addCheckMult("${apt_check_calvinchp}", gen, gold);

  //Create new Matrix Check
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/canary-v1.allele_average_summary.txt");
  check1->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.allele_average_summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  //Create new Matrix Check
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/canary-v1.calls.txt");
  check2->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.calls.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  //Create new Matrix Check
  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/canary-v1.cnv_region_summary.txt");
  check3->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.cnv_region_summary.txt"); 
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  //Create new Matrix Check
  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/canary-v1.confidences.txt");
  check4->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.confidences.txt"); 
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //Create new Matrix Check
  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_testname}/canary-v1.summary.txt");
  check5->addArg("--gold", "${gold_canary_genomewidesnp6_testname}/canary-v1.summary.txt"); 
  check5->addArg("--epsilon", "0.0001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");

  //Return a pointer to the test just created
  return rt;
}
