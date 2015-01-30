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

const int negTest = 1;

RT_Test* qt_doBirdseedSnp6Spf();
//RT_Test* qt_doBirdseedSnp6Full();
RT_Test* qt_doBirdseedSnp6Chp();
RT_Test* qt_doBirdseed2Snp6Chp();
RT_Test* qt_doBrlmmpSnp6Chp();
RT_Test* qt_doBrlmmpSnp5Spf();
//RT_Test* qt_doBrlmmpSnp5CdfReadWriteFeatureEffectsA5();
RT_Test* qt_doBrlmmpSnp5Chp();
//RT_Test* qt_doBrlmmpSnp5KillList();
RT_Test* qt_doLabelZStyHints();
RT_Test* qt_doLabelZStyTestPriors();
RT_Test* qt_doLabelZStyCCS();
RT_Test* qt_doLabelZStyCCSseveral();
RT_Test* qt_doStyCCS_2_20_1_0();
RT_Test* qt_doStyCCS_2_40_0_9_0_033();
RT_Test* qt_doSpf();
RT_Test* qt_doChpFiles();
RT_Test* qt_doFileStyCCS_2_20_1_0();
RT_Test* qt_doWritePriorFileStyCCS_2_20_1_0();
RT_Test* qt_doReadGenotypesIn();
RT_Test* qt_doStyRVT_2_20_1_0();
RT_Test* qt_doStyCES_2_20_1_0();
RT_Test* qt_doStyMva_2_20_1_0();
RT_Test* qt_doStyMva_2_20_1_2();
RT_Test* qt_doStyMva_2_20_0_8_0();
RT_Test* doBirdseedSnp6Spf();
RT_Test* doBirdseedSnp6Full();
RT_Test* doBirdseedSnp6Chp();
RT_Test* doBirdseed2Snp6Chp();
RT_Test* doBrlmmpSnp6Chp();
RT_Test* doBrlmmpSnp5Spf();
RT_Test* doBrlmmpSnp5CdfReadWriteFeatureEffectsA5();
RT_Test* doBrlmmpSnp5Chp();
RT_Test* doBrlmmpSnp5KillList();
RT_Test* doLabelZStyHints();
RT_Test* doLabelZStyTestPriors();
RT_Test* doLabelZStyCCS();
RT_Test* doLabelZStyCCSseveral();
RT_Test* doStyCCS_2_20_1_0();
RT_Test* doStyCCS_2_40_0_9_0_033();
RT_Test* doSpf();
RT_Test* doChpFiles();
RT_Test* doFileStyCCS_2_20_1_0();
RT_Test* doWritePriorFileStyCCS_2_20_1_0();
RT_Test* doReadGenotypesIn();
RT_Test* doStyRVT_2_20_1_0();
RT_Test* doStyCES_2_20_1_0();
RT_Test* doStyMva_2_20_1_0();
RT_Test* doStyMva_2_20_1_2();
RT_Test* doStyMva_2_20_0_8_0();

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if ( !Fs::dirExists("test-generated") ) {
    Fs::mkdir("test-generated", false);
  }

  RT_Main* rMain;
  rMain = new RT_Main("apt-rt-probeset-genotype");
  
  //Store names of tests for Help Message
  //Store it in our RT_Main Instance
  std::vector<std::string> testNameList;
  testNameList.push_back("qt_doBirdseedSnp6Spf");
  //testNameList.push_back("qt_doBirdseedSnp6Full");
  testNameList.push_back("qt_doBirdseedSnp6Chp");
  testNameList.push_back("qt_doBirdseed2Snp6Chp");
  testNameList.push_back("qt_doBrlmmpSnp6Chp");
  testNameList.push_back("qt_doBrlmmpSnp5Spf");
  //testNameList.push_back("qt_doBrlmmpSnp5CdfReadWriteFeatureEffectsA5");
  testNameList.push_back("qt_doBrlmmpSnp5Chp");
  //testNameList.push_back("qt_doBrlmmpSnp5KillList");
  testNameList.push_back("qt_doLabelZStyHints");
  testNameList.push_back("qt_doLabelZStyTestPriors");
  testNameList.push_back("qt_doLabelZStyCCS");
  testNameList.push_back("qt_doLabelZStyCCSseveral");
  testNameList.push_back("qt_doStyCCS_2_20_1_0");
  testNameList.push_back("qt_doStyCCS_2_40_0_9_0_033");
  testNameList.push_back("qt_doSpf");
  testNameList.push_back("qt_doChpFiles");
  testNameList.push_back("qt_doFileStyCCS_2_20_1_0");
  testNameList.push_back("qt_doWritePriorFileStyCCS_2_20_1_0");
  testNameList.push_back("qt_doReadGenotypesIn");
  testNameList.push_back("qt_doStyRVT_2_20_1_0");
  testNameList.push_back("qt_doStyCES_2_20_1_0");
  testNameList.push_back("qt_doStyMva_2_20_1_0");
  testNameList.push_back("qt_doStyMva_2_20_1_2");
  testNameList.push_back("qt_doStyMva_2_20_0_8_0");
  testNameList.push_back("doBirdseedSnp6Spf");
  testNameList.push_back("doBirdseedSnp6Full");
  testNameList.push_back("doBirdseedSnp6Chp");
  testNameList.push_back("doBirdseed2Snp6Chp");
  testNameList.push_back("doBrlmmpSnp6Chp");
  testNameList.push_back("doBrlmmpSnp5Spf");
  testNameList.push_back("doBrlmmpSnp5CdfReadWriteFeatureEffectsA5");
  testNameList.push_back("doBrlmmpSnp5Chp");
  testNameList.push_back("doBrlmmpSnp5KillList");
  testNameList.push_back("doLabelZStyHints");
  testNameList.push_back("doLabelZStyTestPriors");
  testNameList.push_back("doLabelZStyCCS");
  testNameList.push_back("doLabelZStyCCSseveral");
  testNameList.push_back("doStyCCS_2_20_1_0");
  testNameList.push_back("doStyCCS_2_40_0_9_0_033");
  testNameList.push_back("doSpf");
  testNameList.push_back("doChpFiles");
  testNameList.push_back("doFileStyCCS_2_20_1_0");
  testNameList.push_back("doWritePriorFileStyCCS_2_20_1_0");
  testNameList.push_back("doReadGenotypesIn");
  testNameList.push_back("doStyRVT_2_20_1_0");
  testNameList.push_back("doStyCES_2_20_1_0");
  testNameList.push_back("doStyMva_2_20_1_0");
  testNameList.push_back("doStyMva_2_20_1_2");
  testNameList.push_back("doStyMva_2_20_0_8_0");
  rMain->setTestNameList(testNameList);

  //Set Variables
  rMain->setVar("sdktop", "../../..");
  rMain->setVar("gold_lib_genomewidesnp6", "${gold_dir}/lib/GenomeWideSNP_6");
  rMain->setVar("gold_probesetgenotype_genomewidesnp6", "${gold_dir}/probeset-genotype/GenomeWideSNP_6");
  rMain->setVar("gold_probesetgenotype_genomewidesnp6_testname", "${gold_dir}/probeset-genotype/GenomeWideSNP_6/${testname}");
  rMain->setVar("gold_lib_genomewidesnp5", "${gold_dir}/lib/GenomeWideSNP_5");
  rMain->setVar("gold_cel_genomewidesnp5", "${gold_dir}/cel/GenomeWideSNP_5");
  rMain->setVar("gold_probesetgenotype_genomewidesnp5", "${gold_dir}/probeset-genotype/GenomeWideSNP_5");
  rMain->setVar("gold_probesetgenotype_genomewidesnp5_testname", "${gold_probesetgenotype_genomewidesnp5}/${testname}");
  rMain->setVar("gold_cel_mapping250ksty", "${gold_dir}/cel/Mapping250K_Sty");
  rMain->setVar("gold_lib_mapping250ksty", "${gold_dir}/lib/Mapping250K_Sty");
  rMain->setVar("gold_probesetgenotype_mapping250ksty", "${gold_dir}/probeset-genotype/Mapping250K_Sty");
  rMain->setVar("gold_probesetgenotype_mapping250ksty_testname", "${gold_probesetgenotype_mapping250ksty}/${testname}");
  rMain->setVar("chpsuffix", ".chp");
  rMain->setVar("birdseedv2chpsuffix", ".birdseed-v2.chp");
  rMain->setVar("celsuffix", ".CEL");
  rMain->setVar("brlmmp_chp_suffix", ".brlmm-p.chp");

  rMain->setVar("gold_idata_lib_genomewidesnp6", "${gold_dir}/idata/lib/GenomeWideSNP_6");
  rMain->setVar("gold_idata_cel_genomewidesnp6", "${gold_dir}/idata/cel/GenomeWideSNP_6");
  rMain->setVar("gold_idata_pgeno_genomewidesnp6", "${gold_dir}/idata/p-geno/GenomeWideSNP_6");
  rMain->setVar("gold_idata_pgeno_genomewidesnp6_testname", "${gold_dir}/idata/p-geno/GenomeWideSNP_6/${testname}");
  rMain->setVar("gold_idata_cel_genomewidesnp5", "${gold_dir}/idata/cel/GenomeWideSNP_5");
  rMain->setVar("gold_idata_lib_genomewidesnp5", "${gold_dir}/idata/lib/GenomeWideSNP_5");
  rMain->setVar("gold_idata_pgeno_genomewidesnp5", "${gold_dir}/idata/p-geno/GenomeWideSNP_5");
  rMain->setVar("gold_idata_pgeno_genomewidesnp5_testname", "${gold_dir}/idata/p-geno/GenomeWideSNP_5/${testname}");
  rMain->setVar("gold_idata_cel_mapping250ksty", "${gold_dir}/idata/cel/Mapping250K_Sty");
  rMain->setVar("gold_idata_cel_mapping250ksty_testname", "${gold_dir}/idata/cel/Mapping250K_Sty/${testname}");
  rMain->setVar("gold_idata_pgeno_mapping250ksty", "${gold_dir}/idata/p-geno/Mapping250K_Sty");
  rMain->setVar("gold_idata_pgeno_mapping250ksty_testname", "${gold_dir}/idata/p-geno/Mapping250K_Sty/${testname}");
  rMain->setVar("gold_idata_lib_mapping250ksty", "${gold_dir}/idata/lib/Mapping250K_Sty");

  rMain->parseArgv(argc, argv);
  if(argc < 1)
    {
      cout<<"Error: Incorrect number of arguments given"<<endl;
      exit(1);
    }
  else if(rMain->getVarVal("allIntegration") != "")
    {
      rMain->addTest(qt_doBirdseedSnp6Spf());
      //rMain->addTest(qt_doBirdseedSnp6Full());
      rMain->addTest(qt_doBirdseedSnp6Chp());
      rMain->addTest(qt_doBirdseed2Snp6Chp());
      rMain->addTest(qt_doBrlmmpSnp6Chp());
      rMain->addTest(qt_doBrlmmpSnp5Spf());
      //rMain->addTest(qt_doBrlmmpSnp5CdfReadWriteFeatureEffectsA5());
      rMain->addTest(qt_doBrlmmpSnp5Chp());
      //rMain->addTest(qt_doBrlmmpSnp5KillList());
      rMain->addTest(qt_doLabelZStyHints());
      rMain->addTest(qt_doLabelZStyTestPriors());
      rMain->addTest(qt_doLabelZStyCCS());
      rMain->addTest(qt_doLabelZStyCCSseveral());
      rMain->addTest(qt_doStyCCS_2_20_1_0());
      rMain->addTest(qt_doStyCCS_2_40_0_9_0_033());
      rMain->addTest(qt_doSpf());
      rMain->addTest(qt_doChpFiles());
      rMain->addTest(qt_doFileStyCCS_2_20_1_0());
      rMain->addTest(qt_doWritePriorFileStyCCS_2_20_1_0());
      rMain->addTest(qt_doReadGenotypesIn());
      rMain->addTest(qt_doStyRVT_2_20_1_0());
      rMain->addTest(qt_doStyCES_2_20_1_0());
      rMain->addTest(qt_doStyMva_2_20_1_0());
      rMain->addTest(qt_doStyMva_2_20_1_2());
      rMain->addTest(qt_doStyMva_2_20_0_8_0());
    }
  else if(rMain->getVarVal("allNormal") != "")
    {
 rMain->addTest(doBirdseedSnp6Spf());
      rMain->addTest(doBirdseedSnp6Full());
      rMain->addTest(doBirdseedSnp6Chp());
      rMain->addTest(doBirdseed2Snp6Chp());
      rMain->addTest(doBrlmmpSnp6Chp());
      rMain->addTest(doBrlmmpSnp5Spf());
      rMain->addTest(doBrlmmpSnp5CdfReadWriteFeatureEffectsA5());
      rMain->addTest(doBrlmmpSnp5Chp());
      rMain->addTest(doBrlmmpSnp5KillList());
      rMain->addTest(doLabelZStyHints());
      rMain->addTest(doLabelZStyTestPriors());
      rMain->addTest(doLabelZStyCCS());
      rMain->addTest(doLabelZStyCCSseveral());
      rMain->addTest(doStyCCS_2_20_1_0());
      rMain->addTest(doStyCCS_2_40_0_9_0_033());
      rMain->addTest(doSpf());
      rMain->addTest(doChpFiles());
      rMain->addTest(doFileStyCCS_2_20_1_0());
      rMain->addTest(doWritePriorFileStyCCS_2_20_1_0());
      rMain->addTest(doReadGenotypesIn());
      rMain->addTest(doStyRVT_2_20_1_0());
      rMain->addTest(doStyCES_2_20_1_0());
      rMain->addTest(doStyMva_2_20_1_0());
      rMain->addTest(doStyMva_2_20_1_2());
      rMain->addTest(doStyMva_2_20_0_8_0());
    }
  //They want us to run all default tests
  else if(argc == 1 || rMain->getVarVal("all") != "")
    {
      //RUN ALL TESTS OMG
      cout<<"Running all default Regression Tests!"<<endl;
      rMain->addTest(qt_doBirdseedSnp6Spf());
      //rMain->addTest(qt_doBirdseedSnp6Full());
      rMain->addTest(qt_doBirdseedSnp6Chp());
      rMain->addTest(qt_doBirdseed2Snp6Chp());
      rMain->addTest(qt_doBrlmmpSnp6Chp());
      rMain->addTest(qt_doBrlmmpSnp5Spf());
      //rMain->addTest(qt_doBrlmmpSnp5CdfReadWriteFeatureEffectsA5());
      rMain->addTest(qt_doBrlmmpSnp5Chp());
      //rMain->addTest(qt_doBrlmmpSnp5KillList());
      rMain->addTest(qt_doLabelZStyHints());
      rMain->addTest(qt_doLabelZStyTestPriors());
      rMain->addTest(qt_doLabelZStyCCS());
      rMain->addTest(qt_doLabelZStyCCSseveral());
      rMain->addTest(qt_doStyCCS_2_20_1_0());
      rMain->addTest(qt_doStyCCS_2_40_0_9_0_033());
      rMain->addTest(qt_doSpf());
      rMain->addTest(qt_doChpFiles());
      rMain->addTest(qt_doFileStyCCS_2_20_1_0());
      rMain->addTest(qt_doWritePriorFileStyCCS_2_20_1_0());
      rMain->addTest(qt_doReadGenotypesIn());
      rMain->addTest(qt_doStyRVT_2_20_1_0());
      rMain->addTest(qt_doStyCES_2_20_1_0());
      rMain->addTest(qt_doStyMva_2_20_1_0());
      rMain->addTest(qt_doStyMva_2_20_1_2());
      rMain->addTest(qt_doStyMva_2_20_0_8_0());


      rMain->addTest(doBirdseedSnp6Spf());
      rMain->addTest(doBirdseedSnp6Full());
      rMain->addTest(doBirdseedSnp6Chp());
      rMain->addTest(doBirdseed2Snp6Chp());
      rMain->addTest(doBrlmmpSnp6Chp());
      rMain->addTest(doBrlmmpSnp5Spf());
      rMain->addTest(doBrlmmpSnp5CdfReadWriteFeatureEffectsA5());
      rMain->addTest(doBrlmmpSnp5Chp());
      rMain->addTest(doBrlmmpSnp5KillList());
      rMain->addTest(doLabelZStyHints());
      rMain->addTest(doLabelZStyTestPriors());
      rMain->addTest(doLabelZStyCCS());
      rMain->addTest(doLabelZStyCCSseveral());
      rMain->addTest(doStyCCS_2_20_1_0());
      rMain->addTest(doStyCCS_2_40_0_9_0_033());
      rMain->addTest(doSpf());
      rMain->addTest(doChpFiles());
      rMain->addTest(doFileStyCCS_2_20_1_0());
      rMain->addTest(doWritePriorFileStyCCS_2_20_1_0());
      rMain->addTest(doReadGenotypesIn());
      rMain->addTest(doStyRVT_2_20_1_0());
      rMain->addTest(doStyCES_2_20_1_0());
      rMain->addTest(doStyMva_2_20_1_0());
      rMain->addTest(doStyMva_2_20_1_2());
      rMain->addTest(doStyMva_2_20_0_8_0());
    }
  else
    {
      cout<<"Running Built-In Regression Tests"<<endl;

      if(rMain->getVarVal("qt_doBirdseedSnp6Spf") != "")
	{
	  rMain->addTest(qt_doBirdseedSnp6Spf());
	}

      /* 
	 if(rMain->getVarVal("qt_doBirdseedSnp6Full") != "")
	{
	  rMain->addTest(qt_doBirdseedSnp6Full());
	}
      */
      if(rMain->getVarVal("qt_doBirdseedSnp6Chp") != "")
	{
	  rMain->addTest(qt_doBirdseedSnp6Chp());
	}
	
      if(rMain->getVarVal("qt_doBirdseed2Snp6Chp") != "")
	{
	  rMain->addTest(qt_doBirdseed2Snp6Chp());
	}

      if(rMain->getVarVal("qt_doBrlmmpSnp6Chp") != "")
	{
	  rMain->addTest(qt_doBrlmmpSnp6Chp());
	}

      if(rMain->getVarVal("qt_doBrlmmpSnp5Spf") != "")
	{
	  rMain->addTest(qt_doBrlmmpSnp5Spf());
	}

      /*
      if(rMain->getVarVal("qt_doBrlmmpSnp5CdfReadWriteFeatureEffectsA5") != "")
	{
	  rMain->addTest(qt_doBrlmmpSnp5CdfReadWriteFeatureEffectsA5());
	}
      */

      if(rMain->getVarVal("qt_doBrlmmpSnp5Chp") != "")
	{
	  rMain->addTest(qt_doBrlmmpSnp5Chp());
	}

      /*
      if(rMain->getVarVal("qt_doBrlmmpSnp5KillList") != "")
	{
	  rMain->addTest(qt_doBrlmmpSnp5KillList());
	}
      */

      if(rMain->getVarVal("qt_doLabelZStyHints") != "")
	{
	  rMain->addTest(qt_doLabelZStyHints());
	}

      if(rMain->getVarVal("qt_doLabelZStyTestPriors") != "")
	{
	  rMain->addTest(qt_doLabelZStyTestPriors());
	}

      if(rMain->getVarVal("qt_doLabelZStyCCS") != "")
	{
	  rMain->addTest(qt_doLabelZStyCCS());
	}

      if(rMain->getVarVal("qt_doLabelZStyCCSseveral") != "")
	{
	  rMain->addTest(qt_doLabelZStyCCSseveral());
	}

      if(rMain->getVarVal("qt_doStyCCS_2_20_1_0") != "")
	{
	  rMain->addTest(qt_doStyCCS_2_20_1_0());
	}

      if(rMain->getVarVal("qt_doStyCCS_2_40_0_9_0_033") != "")
	{
	  rMain->addTest(qt_doStyCCS_2_40_0_9_0_033());
	}

      if(rMain->getVarVal("qt_doSpf") != "")
	{
	  rMain->addTest(qt_doSpf());
	}

      if(rMain->getVarVal("qt_doChpFiles") != "")
	{
	  rMain->addTest(qt_doChpFiles());
	}

      if(rMain->getVarVal("qt_doFileStyCCS_2_20_1_0") != "")
	{
	  rMain->addTest(qt_doFileStyCCS_2_20_1_0());
	}

      if(rMain->getVarVal("qt_doWritePriorFileStyCCS_2_20_1_0") != "")
	{
	  rMain->addTest(qt_doWritePriorFileStyCCS_2_20_1_0());
	}

      if(rMain->getVarVal("qt_doReadGenotypesIn") != "")
	{
	  rMain->addTest(qt_doReadGenotypesIn());
	}

      if(rMain->getVarVal("qt_doStyRVT_2_20_1_0") != "")
	{
	  rMain->addTest(qt_doStyRVT_2_20_1_0());
	}

      if(rMain->getVarVal("qt_doStyCES_2_20_1_0") != "")
	{
	  rMain->addTest(qt_doStyCES_2_20_1_0());
	}

      if(rMain->getVarVal("qt_doStyMva_2_20_1_0") != "")
	{
	  rMain->addTest(qt_doStyMva_2_20_1_0());
	}

      if(rMain->getVarVal("qt_doStyMva_2_20_1_2") != "")
	{
	  rMain->addTest(qt_doStyMva_2_20_1_2());
	}

      if(rMain->getVarVal("qt_doStyMva_2_20_0_8_0") != "")
	{
	  rMain->addTest(qt_doStyMva_2_20_0_8_0());
	}

      if(rMain->getVarVal("doBirdseedSnp6Spf") != "")
	{
	  rMain->addTest(doBirdseedSnp6Spf());
	}

      if(rMain->getVarVal("doBirdseedSnp6Full") != "")
	{
	  rMain->addTest(doBirdseedSnp6Full());
	}

      if(rMain->getVarVal("doBirdseedSnp6Chp") != "")
	{
	  rMain->addTest(doBirdseedSnp6Chp());
	}
	
      if(rMain->getVarVal("doBirdseed2Snp6Chp") != "")
	{
	  rMain->addTest(doBirdseed2Snp6Chp());
	}

      if(rMain->getVarVal("doBrlmmpSnp6Chp") != "")
	{
	  rMain->addTest(doBrlmmpSnp6Chp());
	}

      if(rMain->getVarVal("doBrlmmpSnp5Spf") != "")
	{
	  rMain->addTest(doBrlmmpSnp5Spf());
	}

      if(rMain->getVarVal("doBrlmmpSnp5CdfReadWriteFeatureEffectsA5") != "")
	{
	  rMain->addTest(doBrlmmpSnp5CdfReadWriteFeatureEffectsA5());
	}

      if(rMain->getVarVal("doBrlmmpSnp5Chp") != "")
	{
	  rMain->addTest(doBrlmmpSnp5Chp());
	}

      if(rMain->getVarVal("doBrlmmpSnp5KillList") != "")
	{
	  rMain->addTest(doBrlmmpSnp5KillList());
	}

      if(rMain->getVarVal("doLabelZStyHints") != "")
	{
	  rMain->addTest(doLabelZStyHints());
	}

      if(rMain->getVarVal("doLabelZStyTestPriors") != "")
	{
	  rMain->addTest(doLabelZStyTestPriors());
	}

      if(rMain->getVarVal("doLabelZStyCCS") != "")
	{
	  rMain->addTest(doLabelZStyCCS());
	}

      if(rMain->getVarVal("doLabelZStyCCSseveral") != "")
	{
	  rMain->addTest(doLabelZStyCCSseveral());
	}

      if(rMain->getVarVal("doStyCCS_2_20_1_0") != "")
	{
	  rMain->addTest(doStyCCS_2_20_1_0());
	}

      if(rMain->getVarVal("doStyCCS_2_40_0_9_0_033") != "")
	{
	  rMain->addTest(doStyCCS_2_40_0_9_0_033());
	}

      if(rMain->getVarVal("doSpf") != "")
	{
	  rMain->addTest(doSpf());
	}

      if(rMain->getVarVal("doChpFiles") != "")
	{
	  rMain->addTest(doChpFiles());
	}

      if(rMain->getVarVal("doFileStyCCS_2_20_1_0") != "")
	{
	  rMain->addTest(doFileStyCCS_2_20_1_0());
	}

      if(rMain->getVarVal("doWritePriorFileStyCCS_2_20_1_0") != "")
	{
	  rMain->addTest(doWritePriorFileStyCCS_2_20_1_0());
	}

      if(rMain->getVarVal("doReadGenotypesIn") != "")
	{
	  rMain->addTest(doReadGenotypesIn());
	}

      if(rMain->getVarVal("doStyRVT_2_20_1_0") != "")
	{
	  rMain->addTest(doStyRVT_2_20_1_0());
	}

      if(rMain->getVarVal("doStyCES_2_20_1_0") != "")
	{
	  rMain->addTest(doStyCES_2_20_1_0());
	}

      if(rMain->getVarVal("doStyMva_2_20_1_0") != "")
	{
	  rMain->addTest(doStyMva_2_20_1_0());
	}

      if(rMain->getVarVal("doStyMva_2_20_1_2") != "")
	{
	  rMain->addTest(doStyMva_2_20_1_2());
	}

      if(rMain->getVarVal("doStyMva_2_20_0_8_0") != "")
	{
	  rMain->addTest(doStyMva_2_20_0_8_0());
	}
    }

  rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBirdseedSnp6Spf()
{
  RT_Test* rt;
  rt = new RT_Test("doBirdseedSnp6Spf");
  //rt->setIntegrationTestBool(true);

  //APT PROBESET GENOTYPE
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis birdseed");
  cmd->addArg("--special-snps", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.specialSNPs");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.SNP_A-4232288.spf");
  cmd->addArg("--read-models-birdseed", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.birdseed.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_idata_pgeno_genomewidesnp6}/fas_cel_files.txt");
  
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/birdseed.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed.calls.txt"); 
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/birdseed.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed.confidences.txt"); 
  check2->addArg("--epsilon", "0.0003");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBirdseedSnp6Spf()
{
  RT_Test* rt;
  rt = new RT_Test("doBirdseedSnp6Spf");

  //APT PROBESET GENOTYPE
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis birdseed");
  cmd->addArg("--special-snps", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.specialSNPs");
  cmd->addArg("--spf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.spf");
  cmd->addArg("--read-models-birdseed", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.birdseed.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_probesetgenotype_genomewidesnp6}/fas_cel_files.txt");
  
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/birdseed.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed.calls.txt"); 
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/birdseed.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed.confidences.txt"); 
  check2->addArg("--epsilon", "0.0003");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

/* Broken/Not Used
///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBirdseedSnp6Full()
{
  RT_Test* rt;
  rt = new RT_Test("doBirdseedSnp6Full");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis", "birdseed");
  cmd->addArg("--special-snps", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.Full.specialSNPs");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.Full.cdf");
  cmd->addArg("--read-models-birdseed", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.birdseed.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_idata_pgeno_genomewidesnp6}/fas_cel_files.txt");

  const char *chpFiles[] = {"NA06985_GW6_C.birdseed", "NA06993_GW6_C.birdseed", "NA07000_GW6_C.birdseed", 
			    "NA07019_GW6_C.birdseed", "NA07029_GW6_C.birdseed", "NA07056_GW6_C.birdseed", 
			    "NA07345_GW6_C.birdseed", "NA10830_GW6_C.birdseed", "NA10831_GW6_C.birdseed", 
			    "NA10839_GW6_C.birdseed", "NA10855_GW6_C.birdseed", "NA10857_GW6_C.birdseed", 
			    "NA10860_GW6_C.birdseed", "NA10863_GW6_C.birdseed", "NA11881_GW6_C.birdseed", 
			    "NA11882_GW6_C.birdseed", "NA11993_GW6_C.birdseed", "NA12003_GW6_C.birdseed", 
			    "NA12004_GW6_C.birdseed",
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_pgeno_genomewidesnp6_testname}/cc-chp/", "${chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", "${chpsuffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/birdseed.calls.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed.calls.txt"); 
  check2->addArg("--epsilon", "0.0000001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/birdseed.confidences.txt");
  check3->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed.confidences.txt"); 
  check3->addArg("--epsilon", "0.0003");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  return rt;
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBirdseedSnp6Full()
{
  RT_Test* rt;
  rt = new RT_Test("doBirdseedSnp6Full");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis", "birdseed");
  cmd->addArg("--special-snps", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.Full.specialSNPs");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.Full.cdf");
  cmd->addArg("--read-models-birdseed", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.birdseed.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_probesetgenotype_genomewidesnp6}/fas_cel_files.txt");

  const char *chpFiles[] = {"NA06985_GW6_C.birdseed", "NA06993_GW6_C.birdseed", "NA07000_GW6_C.birdseed", 
			    "NA07019_GW6_C.birdseed", "NA07029_GW6_C.birdseed", "NA07056_GW6_C.birdseed", 
			    "NA07345_GW6_C.birdseed", "NA10830_GW6_C.birdseed", "NA10831_GW6_C.birdseed", 
			    "NA10839_GW6_C.birdseed", "NA10855_GW6_C.birdseed", "NA10857_GW6_C.birdseed", 
			    "NA10860_GW6_C.birdseed", "NA10863_GW6_C.birdseed", "NA11881_GW6_C.birdseed", 
			    "NA11882_GW6_C.birdseed", "NA11993_GW6_C.birdseed", "NA12003_GW6_C.birdseed", 
			    "NA12004_GW6_C.birdseed",
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetgenotype_genomewidesnp6_testname}/cc-chp/", "${chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", "${chpsuffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/birdseed.calls.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed.calls.txt"); 
  check2->addArg("--epsilon", "0.0000001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/birdseed.confidences.txt");
  check3->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed.confidences.txt"); 
  check3->addArg("--epsilon", "0.0003");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBirdseedSnp6Chp()
{
  RT_Test* rt;
  rt = new RT_Test("doBirdseedSnp6Chp");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis", "birdseed");
  cmd->addArg("--special-snps", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.specialSNPs");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.SNP_A-4232288.spf");
  cmd->addArg("--read-models-birdseed", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.birdseed.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_idata_pgeno_genomewidesnp6}/fas_cel_files.txt");

  const char *chpFiles[] = {"NA06985_GW6_C.birdseed", "NA06993_GW6_C.birdseed", "NA07000_GW6_C.birdseed", 
			    "NA07019_GW6_C.birdseed", "NA07029_GW6_C.birdseed", "NA07056_GW6_C.birdseed", 
			    "NA07345_GW6_C.birdseed", "NA10830_GW6_C.birdseed", "NA10831_GW6_C.birdseed", 
			    "NA10839_GW6_C.birdseed", "NA10855_GW6_C.birdseed", "NA10857_GW6_C.birdseed", 
			    "NA10860_GW6_C.birdseed", "NA10863_GW6_C.birdseed", "NA11881_GW6_C.birdseed", 
			    "NA11882_GW6_C.birdseed", "NA11993_GW6_C.birdseed", "NA12003_GW6_C.birdseed", 
			    "NA12004_GW6_C.birdseed",
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_pgeno_genomewidesnp6_testname}/cc-chp/", "${chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", "${chpsuffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");

  //MIXED FILE CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_testname}/birdseed.report.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed.report.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/birdseed.calls.txt");
  check3->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed.calls.txt"); 
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/birdseed.confidences.txt");
  check4->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed.confidences.txt"); 
  check4->addArg("--epsilon", "0.0003");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBirdseedSnp6Chp()
{
  RT_Test* rt;
  rt = new RT_Test("doBirdseedSnp6Chp");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis", "birdseed");
  cmd->addArg("--special-snps", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.specialSNPs");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--read-models-birdseed", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.birdseed.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_probesetgenotype_genomewidesnp6}/fas_cel_files.txt");

  const char *chpFiles[] = {"NA06985_GW6_C.birdseed", "NA06993_GW6_C.birdseed", "NA07000_GW6_C.birdseed", 
			    "NA07019_GW6_C.birdseed", "NA07029_GW6_C.birdseed", "NA07056_GW6_C.birdseed", 
			    "NA07345_GW6_C.birdseed", "NA10830_GW6_C.birdseed", "NA10831_GW6_C.birdseed", 
			    "NA10839_GW6_C.birdseed", "NA10855_GW6_C.birdseed", "NA10857_GW6_C.birdseed", 
			    "NA10860_GW6_C.birdseed", "NA10863_GW6_C.birdseed", "NA11881_GW6_C.birdseed", 
			    "NA11882_GW6_C.birdseed", "NA11993_GW6_C.birdseed", "NA12003_GW6_C.birdseed", 
			    "NA12004_GW6_C.birdseed",
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetgenotype_genomewidesnp6_testname}/cc-chp/", "${chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", "${chpsuffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");

  //MIXED FILE CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_testname}/birdseed.report.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed.report.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/birdseed.calls.txt");
  check3->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed.calls.txt"); 
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/birdseed.confidences.txt");
  check4->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed.confidences.txt"); 
  check4->addArg("--epsilon", "0.0003");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBirdseed2Snp6Chp()
{
  RT_Test* rt;
  rt = new RT_Test("doBirdseed2Snp6Chp");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis", "birdseed-v2");
  cmd->addArg("--special-snps", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.specialSNPs");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.spf");
  cmd->addArg("--read-models-birdseed", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.birdseed.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_idata_pgeno_genomewidesnp6}/fas_cel_files.txt");
  cmd->addArg("--set-gender-method", "cn-probe-chrXY-ratio");
  cmd->addArg("--chrX-probes", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.chrXprobes");
  cmd->addArg("--chrY-probes", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.chrYprobes");

  const char *chpFiles[] = {"NA06985_GW6_C", "NA06993_GW6_C", "NA07000_GW6_C", 
			    "NA07019_GW6_C", "NA07029_GW6_C", "NA07056_GW6_C", 
			    "NA07345_GW6_C", "NA10830_GW6_C", "NA10831_GW6_C", 
			    "NA10839_GW6_C", "NA10855_GW6_C", "NA10857_GW6_C", 
			    "NA10860_GW6_C", "NA10863_GW6_C", "NA11881_GW6_C", 
			    "NA11882_GW6_C", "NA11993_GW6_C", "NA12003_GW6_C", 
			    "NA12004_GW6_C",
			    NULL
  };
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_pgeno_genomewidesnp6_testname}/cc-chp/", "${birdseedv2chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", "${birdseedv2chpsuffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");

  
  //MIXED FILE CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_testname}/birdseed-v2.report.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed-v2.report.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/birdseed-v2.calls.txt");
  check3->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed-v2.calls.txt"); 
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/birdseed-v2.confidences.txt");
  check4->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/birdseed-v2.confidences.txt"); 
  check4->addArg("--epsilon", "0.0003");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBirdseed2Snp6Chp()
{
  RT_Test* rt;
  rt = new RT_Test("doBirdseed2Snp6Chp");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis", "birdseed-v2");
  cmd->addArg("--special-snps", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.specialSNPs");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--read-models-birdseed", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.birdseed.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_probesetgenotype_genomewidesnp6}/fas_cel_files.txt");
  cmd->addArg("--set-gender-method", "cn-probe-chrXY-ratio");
  cmd->addArg("--chrX-probes", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.chrXprobes");
  cmd->addArg("--chrY-probes", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.chrYprobes");

  const char *chpFiles[] = {"NA06985_GW6_C", "NA06993_GW6_C", "NA07000_GW6_C", 
			    "NA07019_GW6_C", "NA07029_GW6_C", "NA07056_GW6_C", 
			    "NA07345_GW6_C", "NA10830_GW6_C", "NA10831_GW6_C", 
			    "NA10839_GW6_C", "NA10855_GW6_C", "NA10857_GW6_C", 
			    "NA10860_GW6_C", "NA10863_GW6_C", "NA11881_GW6_C", 
			    "NA11882_GW6_C", "NA11993_GW6_C", "NA12003_GW6_C", 
			    "NA12004_GW6_C",
			    NULL
  };
  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetgenotype_genomewidesnp6_testname}/cc-chp/", "${birdseedv2chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", "${birdseedv2chpsuffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");

  
  //MIXED FILE CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_testname}/birdseed-v2.report.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed-v2.report.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/birdseed-v2.calls.txt");
  check3->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed-v2.calls.txt"); 
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/birdseed-v2.confidences.txt");
  check4->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/birdseed-v2.confidences.txt"); 
  check4->addArg("--epsilon", "0.0003");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBrlmmpSnp6Chp()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp6Chp");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis", "brlmm-p");
  cmd->addArg("--special-snps", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.specialSNPs");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.spf");
  cmd->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp6}/GenomeWideSNP_6.brlmm-p.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_idata_pgeno_genomewidesnp6}/fas_cel_files.txt");

  const char *chpFiles[] = {"NA06985_GW6_C.brlmm-p", "NA06993_GW6_C.brlmm-p", "NA07000_GW6_C.brlmm-p", 
			    "NA07019_GW6_C.brlmm-p", "NA07029_GW6_C.brlmm-p", "NA07056_GW6_C.brlmm-p", 
			    "NA07345_GW6_C.brlmm-p", "NA10830_GW6_C.brlmm-p", "NA10831_GW6_C.brlmm-p", 
			    "NA10839_GW6_C.brlmm-p", "NA10855_GW6_C.brlmm-p", "NA10857_GW6_C.brlmm-p", 
			    "NA10860_GW6_C.brlmm-p", "NA10863_GW6_C.brlmm-p", "NA11881_GW6_C.brlmm-p", 
			    "NA11882_GW6_C.brlmm-p", "NA11993_GW6_C.brlmm-p", "NA12003_GW6_C.brlmm-p", 
			    "NA12004_GW6_C.brlmm-p", 
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_pgeno_genomewidesnp6_testname}/cc-chp/", "${chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", "${chpsuffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");
  check1->addArg("--bCheckHeaders", "false");

  //MIXED FILE CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.report.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/brlmm-p.report.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check3->addArg("--gold", "${gold_idata_pgeno_genomewidesnp6_testname}/brlmm-p.calls.txt"); 
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  ///@todo Fairly loose stringency on brlmmp confidence values - chp output
  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${gold_idata_pgeno_genomewidesnp6_testname}/brlmm-p.confidences.txt");
  check4->addArg("--gold", "${out_testname}/brlmm-p.confidences.txt"); 
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBrlmmpSnp6Chp()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp6Chp");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--analysis", "brlmm-p");
  cmd->addArg("--special-snps", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.specialSNPs");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.brlmm-p.models");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_probesetgenotype_genomewidesnp6}/fas_cel_files.txt");

  const char *chpFiles[] = {"NA06985_GW6_C.brlmm-p", "NA06993_GW6_C.brlmm-p", "NA07000_GW6_C.brlmm-p", 
			    "NA07019_GW6_C.brlmm-p", "NA07029_GW6_C.brlmm-p", "NA07056_GW6_C.brlmm-p", 
			    "NA07345_GW6_C.brlmm-p", "NA10830_GW6_C.brlmm-p", "NA10831_GW6_C.brlmm-p", 
			    "NA10839_GW6_C.brlmm-p", "NA10855_GW6_C.brlmm-p", "NA10857_GW6_C.brlmm-p", 
			    "NA10860_GW6_C.brlmm-p", "NA10863_GW6_C.brlmm-p", "NA11881_GW6_C.brlmm-p", 
			    "NA11882_GW6_C.brlmm-p", "NA11993_GW6_C.brlmm-p", "NA12003_GW6_C.brlmm-p", 
			    "NA12004_GW6_C.brlmm-p", 
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetgenotype_genomewidesnp6_testname}/cc-chp/", "${chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", "${chpsuffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");
check1->addArg("--bCheckHeaders", "false");

  //MIXED FILE CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.report.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/brlmm-p.report.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check3->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/brlmm-p.calls.txt"); 
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check4->addArg("--gold", "${gold_probesetgenotype_genomewidesnp6_testname}/brlmm-p.confidences.txt"); 
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBrlmmpSnp5Spf()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5Spf");

  const char *chpFiles[] = {
    "NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C","NA06994_GW5_C","NA07000_GW5_C",
    "NA07019_GW5_C","NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C","NA07048_GW5_C",
    "NA07055_GW5_C","NA07056_GW5_C","NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
    "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C","NA10838_GW5_C","NA10839_GW5_C",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--spf-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.spf");
  cmd->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd->addArg("--analysis", "brlmm-p");
  cmd->addArg("--out-dir", "${out_testname}");
  std::vector<std::string> fullChpFiles = Util::addPrefixSuffix(chpFiles, " ${gold_cel_genomewidesnp5}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5_testname}/brlmm-p.calls.txt"); 
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5_testname}/brlmm-p.confidences.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");


  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBrlmmpSnp5Spf()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5Spf");
  //rt->setIntegrationTestBool(true);

  const char *chpFiles[] = {
    "NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C","NA06994_GW5_C","NA07000_GW5_C",
    "NA07019_GW5_C","NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C","NA07048_GW5_C",
    "NA07055_GW5_C","NA07056_GW5_C","NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
    "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C","NA10838_GW5_C","NA10839_GW5_C",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.spf");
  cmd->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd->addArg("--analysis", "brlmm-p");
  cmd->addArg("--out-dir", "${out_testname}");
  std::vector<std::string> fullChpFiles = Util::addPrefixSuffix(chpFiles, " ${gold_idata_cel_genomewidesnp5}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5_testname}/brlmm-p.calls.txt"); 
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5_testname}/brlmm-p.confidences.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");


  return rt;
}
/* Broken/Not Used
///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBrlmmpSnp5CdfReadWriteFeatureEffectsA5()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5CdfReadWriteFeatureEffectsA5");
  RT_Test* rt1;
  rt1 = new RT_Test("doBrlmmpSnp5CdfWriteFeatureEffectsA5");
  RT_Test* rt2;
  rt2 = new RT_Test("doBrlmmpSnp5CdfReadFeatureEffectsA5");

  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
			    "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
			    "NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C",
			    "NA07048_GW5_C","NA07055_GW5_C","NA07056_GW5_C",
			    "NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
			    "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C",
			    "NA10838_GW5_C","NA10839_GW5_C",
			    NULL
  };

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_genotype}");
  cmd1->addArg("--qmethod-spec", "med-polish");
  cmd1->addArg("--analysis", "brlmm-p");
  cmd1->addArg("--spf-file", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.spf");
  cmd1->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd1->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd1->addArg("--out-dir", "${out_testname}");
  cmd1->addArg("--a5-feature-effects");
  cmd1->addArg("--summaries");
  
  std::vector<std::string> fullChpFiles1 = Util::addPrefixSuffix(chpFiles, "${gold_idata_cel_genomewidesnp5}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullChpFiles1.begin();i < fullChpFiles1.end();i++)
    {
      std::string stringArg = *i;
      cmd1->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm-p.summary.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5}/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.summary.txt"); 
  check1->addArg("--epsilon", "0.00001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");
	       
  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_genotype}");
  cmd2->addArg("--qmethod-spec", "med-polish");
  cmd2->addArg("--analysis", "brlmm-p");
  cmd2->addArg("--spf-file", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.spf");
  cmd2->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd2->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  //cmd2->addArg("--a5-feature-effects-input-file", "${gold_probesetgenotype_genomewidesnp5}/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.feature-response.a5.translated");
  cmd2->addArg("--a5-feature-effects-input-file", "test-generated/doBrlmmpSnp5CdfWriteFeatureEffectsA5/brlmm-p.feature-response.a5");
  cmd2->addArg("--a5-feature-effects-input-name", "brlmm-p.feature-response");
  cmd2->addArg("--out-dir", "${out_testname}");
  cmd2->addArg("--summaries");

  std::vector<std::string> fullChpFiles2 = Util::addPrefixSuffix(chpFiles, " ${gold_idata_cel_genomewidesnp5}/", "${celsuffix}");

  for(std::vector<std::string>::iterator i = fullChpFiles2.begin();i < fullChpFiles2.end();i++)
    {
      std::string stringArg = *i;
      cmd2->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.summary.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5}/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.summary.txt"); 
  check2->addArg("--epsilon", "0.00002");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;								 
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBrlmmpSnp5CdfReadWriteFeatureEffectsA5()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5CdfReadWriteFeatureEffectsA5");
  RT_Test* rt1;
  rt1 = new RT_Test("doBrlmmpSnp5CdfWriteFeatureEffectsA5");
  RT_Test* rt2;
  rt2 = new RT_Test("doBrlmmpSnp5CdfReadFeatureEffectsA5");

  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
			    "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
			    "NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C",
			    "NA07048_GW5_C","NA07055_GW5_C","NA07056_GW5_C",
			    "NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
			    "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C",
			    "NA10838_GW5_C","NA10839_GW5_C",
			    NULL
  };
  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_genotype}");
  cmd1->addArg("--qmethod-spec", "med-polish");
  cmd1->addArg("--analysis", "brlmm-p");
  cmd1->addArg("--cdf-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.cdf");
  cmd1->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd1->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd1->addArg("--out-dir", "${out_testname}");
  cmd1->addArg("--a5-feature-effects");
  cmd1->addArg("--summaries");
  
  std::vector<std::string> fullChpFiles1 = Util::addPrefixSuffix(chpFiles, "${gold_cel_genomewidesnp5}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullChpFiles1.begin();i < fullChpFiles1.end();i++)
    {
      std::string stringArg = *i;
      cmd1->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm-p.summary.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5}/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.summary.txt"); 
  check1->addArg("--epsilon", "0.00001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");
	       
  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_genotype}");
  cmd2->addArg("--qmethod-spec", "med-polish");
  cmd2->addArg("--analysis", "brlmm-p");
  cmd2->addArg("--cdf-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.cdf");
  cmd2->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd2->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd2->addArg("--a5-feature-effects-input-file", "${gold_probesetgenotype_genomewidesnp5}/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.feature-response.a5.translated");
  //cmd2->addArg("--a5-feature-effects-input-file", "test-generated/doBrlmmpSnp5CdfWriteFeatureEffectsA5/brlmm-p.feature-response.a5");
  cmd2->addArg("--a5-feature-effects-input-name", "brlmm-p.feature-response");
  cmd2->addArg("--out-dir", "${out_testname}");
  cmd2->addArg("--summaries");
  std::vector<std::string> fullChpFiles2 = Util::addPrefixSuffix(chpFiles, " ${gold_cel_genomewidesnp5}/", "${celsuffix}");

  for(std::vector<std::string>::iterator i = fullChpFiles2.begin();i < fullChpFiles2.end();i++)
    {
      std::string stringArg = *i;
      cmd2->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.summary.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5}/doBrlmmpSnp5CdfReadWriteFeatureEffectsA5/brlmm-p.summary.txt"); 
  check2->addArg("--epsilon", "0.00002");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;								 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBrlmmpSnp5Chp()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5Chp");
  //rt->setIntegrationTestBool(true);

  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
			    "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
			    "NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C",
			    "NA07048_GW5_C","NA07055_GW5_C","NA07056_GW5_C",
			    "NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
			    "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C",
			    "NA10838_GW5_C","NA10839_GW5_C",
			    NULL
  };
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.spf");
  cmd->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd->addArg("--analysis", "brlmm-p");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("--summaries");
  cmd->addArg("--cc-chp-output");
 
  std::vector<std::string> fullChpFiles = Util::addPrefixSuffix(chpFiles, " ${gold_idata_cel_genomewidesnp5}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_pgeno_genomewidesnp5_testname}/cc-chp/", "${brlmmp_chp_suffix}");
  gen = Util::addPrefixSuffix(chpFiles,"${out_testname}/cc-chp/", "${brlmmp_chp_suffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");
  check1->addArg("--bCheckHeaders", "false");

  //MIXED FILE CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.report.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5_testname}/brlmm-p.report.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check3->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5_testname}/brlmm-p.calls.txt");
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1");
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/brlmm-p.normalized-summary.txt");
  check4->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5_testname}/brlmm-p.normalized-summary.txt");
  check4->addArg("--epsilon", "0.001");
  check4->addArg("--rowSkip", "1");
  check4->addArg("--columnSkip", "2");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_testname}/brlmm-p.summary.txt");
  check5->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5_testname}/brlmm-p.summary.txt");
  check5->addArg("--epsilon", "0.001");
  check5->addArg("--rowSkip", "1");
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check6 = rt->newCheck();
  check6->setExe("${apt_check_matrix}");
  check6->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check6->addArg("--gold", "${gold_idata_pgeno_genomewidesnp5_testname}/brlmm-p.confidences.txt");
  check6->addArg("--epsilon", "0.0001");
  check6->addArg("--rowSkip", "1");
  check6->addArg("--columnSkip", "1");
  check6->addArg("--matchNames", "false");
  check6->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBrlmmpSnp5Chp()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5Chp");

  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
			    "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
			    "NA07022_GW5_C","NA07029_GW5_C","NA07034_GW5_C",
			    "NA07048_GW5_C","NA07055_GW5_C","NA07056_GW5_C",
			    "NA07345_GW5_C","NA07348_GW5_C","NA07357_GW5_C",
			    "NA10830_GW5_C","NA10831_GW5_C","NA10835_GW5_C",
			    "NA10838_GW5_C","NA10839_GW5_C",
			    NULL
  };
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.cdf");
  cmd->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd->addArg("--analysis", "brlmm-p");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("--summaries");
  cmd->addArg("--cc-chp-output");

  std::vector<std::string> fullChpFiles = Util::addPrefixSuffix(chpFiles, " ${gold_cel_genomewidesnp5}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetgenotype_genomewidesnp5_testname}/cc-chp/", "${brlmmp_chp_suffix}");
  gen = Util::addPrefixSuffix(chpFiles,"${out_testname}/cc-chp/", "${brlmmp_chp_suffix}");

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_calvinchp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  check1->addArg("--diffAllowed", "0");
  check1->addArg("--prefix", "apt-");
  check1->addArg("--epsilon", "0.001");
  check1->addArg("--bCheckHeaders", "false");

  //MIXED FILE CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_testname}/brlmm-p.report.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5_testname}/brlmm-p.report.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/brlmm-p.calls.txt");
  check3->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5_testname}/brlmm-p.calls.txt");
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1");
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/brlmm-p.normalized-summary.txt");
  check4->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5_testname}/brlmm-p.normalized-summary.txt");
  check4->addArg("--epsilon", "0.001");
  check4->addArg("--rowSkip", "1");
  check4->addArg("--columnSkip", "2");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_testname}/brlmm-p.summary.txt");
  check5->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5_testname}/brlmm-p.summary.txt");
  check5->addArg("--epsilon", "0.001");
  check5->addArg("--rowSkip", "1");
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check6 = rt->newCheck();
  check6->setExe("${apt_check_matrix}");
  check6->addArg("--gen", "${out_testname}/brlmm-p.confidences.txt");
  check6->addArg("--gold", "${gold_probesetgenotype_genomewidesnp5_testname}/brlmm-p.confidences.txt");
  check6->addArg("--epsilon", "0.0001");
  check6->addArg("--rowSkip", "1");
  check6->addArg("--columnSkip", "1");
  check6->addArg("--matchNames", "false");
  check6->addArg("--allowedMismatch", "0");

  return rt;
}

/* Broken/Unused Test
///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doBrlmmpSnp5KillList()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5KillList");
  RT_Test* rt1;
  rt1 = new RT_Test("doBrlmmpSnp5KillList-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("doBrlmmpSnp5KillList-part2");
  RT_Test* rt3;
  rt3 = new RT_Test("doBrlmmpSnp5KillList-part3");

  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
			    "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
			    NULL
  };
  std::vector<std::string> fullChpFiles = Util::addPrefixSuffix(chpFiles, " ${gold_idata_cel_genomewidesnp5}/", "${celsuffix}");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_genotype}");
  cmd1->addArg("--spf-file", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.spf");
  cmd1->addArg("--kill-list", "${gold_idata_lib_genomewidesnp5}/test-kill-list.txt");
  cmd1->addArg("-s", "${gold_idata_lib_genomewidesnp5}/test-kill-list.snp-list.txt");
  cmd1->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd1->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd1->addArg("--analysis", "brlmm-p");
  cmd1->addArg("--out-dir", "test-generated/doBrlmmpSnp5KillList");
  cmd1->addArg("--summaries");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd1->addArg(stringArg);
    }
  
  //CMD
  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_genotype}");
  cmd2->addArg("--force");
  cmd2->addArg("--spf-file", "${gold_idata_lib_genomewidesnp5}/test-kill-list.spf");
  cmd2->addArg("-s", "${gold_idata_lib_genomewidesnp5}/test-kill-list.snp-list.txt");
  cmd2->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd2->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd2->addArg("--analysis", "brlmm-p");
  cmd2->addArg("--out-dir", "test-generated/doBrlmmpSnp5KillList/mask2");
  cmd2->addArg("--summaries");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd2->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt2->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.calls.txt");
  check1->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.normalized-summary.txt");
  check2->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.normalized-summary.txt");
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "2");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check3 = rt2->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.summary.txt");
  check3->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.summary.txt");
  check3->addArg("--epsilon", "0.001");
  check3->addArg("--rowSkip", "1");
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check4 = rt2->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.confidences.txt");
  check4->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.confidences.txt");
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1");
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //CMD
  RT_Cmd* cmd3 = rt3->newCmd();
  cmd3->setExe("${apt_probeset_genotype}");
  cmd3->addArg("--spf-file", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.spf");
  cmd3->addArg("-s", "${gold_idata_lib_genomewidesnp5}/test-kill-list.snp-list.txt");
  cmd3->addArg("--chrX-snps", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd3->addArg("--read-models-brlmmp", "${gold_idata_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd3->addArg("--analysis", "brlmm-p");
  cmd3->addArg("--out-dir", "test-generated/doBrlmmpSnp5KillList/mask3");
  cmd3->addArg("--summaries");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd3->addArg(stringArg);
    }

 
  //MATRIX CHECK
  //Negative Test
  RT_Check* check5 = rt3->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.calls.txt");
  check5->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.calls.txt");
  check5->addArg("--epsilon", "0.0000001");
  check5->addArg("--rowSkip", "1");
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");
  //Negative Test
  check5->setNegTest();
  //rt3->addCheck(MATRIX_CHECK, checkArgs5, 1);

  //MATRIX CHECK
  //Negative Test
  RT_Check* check6 = rt3->newCheck();
  check6->setExe("${apt_check_matrix}");
  check6->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.normalized-summary.txt");
  check6->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.normalized-summary.txt");
  check6->addArg("--epsilon", "0.001");
  check6->addArg("--rowSkip", "1");
  check6->addArg("--columnSkip", "2");
  check6->addArg("--matchNames", "false");
  check6->addArg("--allowedMismatch", "0");
  //Negative Test
  check6->setNegTest();
  //rt3->addCheck(MATRIX_CHECK, checkArgs6, 1);

  //MATRIX CHECK
  //Negative Test
  RT_Check* check7 = rt3->newCheck();
  check7->setExe("${apt_check_matrix}");
  check7->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.summary.txt");
  check7->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.summary.txt");
  check7->addArg("--epsilon", "0.001");
  check7->addArg("--rowSkip", "1");
  check7->addArg("--columnSkip", "1");
  check7->addArg("--matchNames", "false");
  check7->addArg("--allowedMismatch", "0");
  //Negative Test
  check7->setNegTest();
  //rt3->addCheck(MATRIX_CHECK, checkArgs7, 1);

  //MATRIX CHECK
  //Negative Test
  RT_Check* check8 = rt3->newCheck();
  check8->setExe("${apt_check_matrix}");
  check8->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.confidences.txt");
  check8->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.confidences.txt");
  check8->addArg("--epsilon", "0.0001");
  check8->addArg("--rowSkip", "1");
  check8->addArg("--columnSkip", "1");
  check8->addArg("--matchNames", "false");
  check8->addArg("--allowedMismatch", "0");
  //Negative Test
  check8->setNegTest();

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  rt->addSubtest(rt3);
  return rt;
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBrlmmpSnp5KillList()
{
  RT_Test* rt;
  rt = new RT_Test("doBrlmmpSnp5KillList");
  RT_Test* rt1;
  rt1 = new RT_Test("doBrlmmpSnp5KillList-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("doBrlmmpSnp5KillList-part2");
  RT_Test* rt3;
  rt3 = new RT_Test("doBrlmmpSnp5KillList-part3");

  const char *chpFiles[] = {"NA06985_GW5_C","NA06991_GW5_C","NA06993_GW5_C",
			    "NA06994_GW5_C","NA07000_GW5_C","NA07019_GW5_C",
			    NULL
  };
  std::vector<std::string> fullChpFiles = Util::addPrefixSuffix(chpFiles, " ${gold_cel_genomewidesnp5}/", "${celsuffix}");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_genotype}");
  cmd1->addArg("--cdf-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.cdf");
  cmd1->addArg("--kill-list", "${gold_lib_genomewidesnp5}/test-kill-list.txt");
  cmd1->addArg("-s", "${gold_lib_genomewidesnp5}/test-kill-list.snp-list.txt");
  cmd1->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd1->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd1->addArg("--analysis", "brlmm-p");
  cmd1->addArg("--out-dir", "test-generated/doBrlmmpSnp5KillList");
  cmd1->addArg("--summaries");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd1->addArg(stringArg);
    }
  
  //CMD
  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_genotype}");
  cmd2->addArg("--force");
  cmd2->addArg("--cdf-file", "${gold_lib_genomewidesnp5}/test-kill-list.cdf");
  cmd2->addArg("-s", "${gold_lib_genomewidesnp5}/test-kill-list.snp-list.txt");
  cmd2->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd2->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd2->addArg("--analysis", "brlmm-p");
  cmd2->addArg("--out-dir", "test-generated/doBrlmmpSnp5KillList/mask2");
  cmd2->addArg("--summaries");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd2->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt2->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.calls.txt");
  check1->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.normalized-summary.txt");
  check2->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.normalized-summary.txt");
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "2");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check3 = rt2->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.summary.txt");
  check3->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.summary.txt");
  check3->addArg("--epsilon", "0.001");
  check3->addArg("--rowSkip", "1");
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check4 = rt2->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask2/brlmm-p.confidences.txt");
  check4->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.confidences.txt");
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1");
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //CMD
  RT_Cmd* cmd3 = rt3->newCmd();
  cmd3->setExe("${apt_probeset_genotype}");
  cmd3->addArg("--cdf-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.cdf");
  cmd3->addArg("-s", "${gold_lib_genomewidesnp5}/test-kill-list.snp-list.txt");
  cmd3->addArg("--chrX-snps", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.chrx");
  cmd3->addArg("--read-models-brlmmp", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.models");
  cmd3->addArg("--analysis", "brlmm-p");
  cmd3->addArg("--out-dir", "test-generated/doBrlmmpSnp5KillList/mask3");
  cmd3->addArg("--summaries");
  for(std::vector<std::string>::iterator i = fullChpFiles.begin();i < fullChpFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd3->addArg(stringArg);
    }

 
  //MATRIX CHECK
  //Negative Test
  RT_Check* check5 = rt3->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.calls.txt");
  check5->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.calls.txt");
  check5->addArg("--epsilon", "0.0000001");
  check5->addArg("--rowSkip", "1");
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");
  //Negative Test
  check5->setNegTest();
  //rt3->addCheck(MATRIX_CHECK, checkArgs5, 1);

  //MATRIX CHECK
  //Negative Test
  RT_Check* check6 = rt3->newCheck();
  check6->setExe("${apt_check_matrix}");
  check6->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.normalized-summary.txt");
  check6->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.normalized-summary.txt");
  check6->addArg("--epsilon", "0.001");
  check6->addArg("--rowSkip", "1");
  check6->addArg("--columnSkip", "2");
  check6->addArg("--matchNames", "false");
  check6->addArg("--allowedMismatch", "0");
  //Negative Test
  check6->setNegTest();
  //rt3->addCheck(MATRIX_CHECK, checkArgs6, 1);

  //MATRIX CHECK
  //Negative Test
  RT_Check* check7 = rt3->newCheck();
  check7->setExe("${apt_check_matrix}");
  check7->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.summary.txt");
  check7->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.summary.txt");
  check7->addArg("--epsilon", "0.001");
  check7->addArg("--rowSkip", "1");
  check7->addArg("--columnSkip", "1");
  check7->addArg("--matchNames", "false");
  check7->addArg("--allowedMismatch", "0");
  //Negative Test
  check7->setNegTest();
  //rt3->addCheck(MATRIX_CHECK, checkArgs7, 1);

  //MATRIX CHECK
  //Negative Test
  RT_Check* check8 = rt3->newCheck();
  check8->setExe("${apt_check_matrix}");
  check8->addArg("--gen", "test-generated/doBrlmmpSnp5KillList/mask3/brlmm-p.confidences.txt");
  check8->addArg("--gold", "test-generated/doBrlmmpSnp5KillList/brlmm-p.confidences.txt");
  check8->addArg("--epsilon", "0.0001");
  check8->addArg("--rowSkip", "1");
  check8->addArg("--columnSkip", "1");
  check8->addArg("--matchNames", "false");
  check8->addArg("--allowedMismatch", "0");
  //Negative Test
  check8->setNegTest();

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  rt->addSubtest(rt3);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doLabelZStyHints()
{
  RT_Test* rt;
  rt = new RT_Test("doLabelZStyHints");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--genotypes", "${gold_idata_lib_mapping250ksty}/call.hints");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt"); 
  cmd->addArg("--summaries"); 
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample"); 
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm-p.lowprecision=true.transform=ccs.K=4.MS=2.hints=1.CP=8"); 
  cmd->addArg("--write-models");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  //cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels,"${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doLabelZStyHints()
{
  RT_Test* rt;
  rt = new RT_Test("doLabelZStyHints");

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--genotypes", "${gold_probesetgenotype_mapping250ksty}/call.hints");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt"); 
  cmd->addArg("--summaries"); 
  cmd->addArg("--list-sample"); 
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm-p.lowprecision=true.transform=ccs.K=4.MS=2.hints=1.CP=8"); 
  cmd->addArg("--write-models");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels,"${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doLabelZStyTestPriors()
{
  RT_Test* rt;
  rt = new RT_Test("doLabelZStyTestPriors");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--read-models-brlmmp", "${gold_idata_lib_mapping250ksty}/test.snp-prior.ref");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--summaries");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels," ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");


  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doLabelZStyTestPriors()
{
  RT_Test* rt;
  rt = new RT_Test("doLabelZStyTestPriors");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--read-models-brlmmp", "${gold_probesetgenotype_mapping250ksty}/test.snp-prior.ref");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--summaries");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels," ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");


  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doLabelZStyCCS()
{
 RT_Test* rt;
  rt = new RT_Test("doLabelZStyCCS");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--summaries");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullSty20celsFiles =Util::addPrefixSuffix(sty20cels, " ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doLabelZStyCCS()
{
 RT_Test* rt;
  rt = new RT_Test("doLabelZStyCCS");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--summaries");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm-p.CM=1.bins=100");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles =Util::addPrefixSuffix(sty20cels, " ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doLabelZStyCCSseveral()
{
  RT_Test* rt;
  rt = new RT_Test("doLabelZStyCCSseveral");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--summaries");
  cmd->addArg("--prior-size" ,"1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm-p.lowprecision=true.transform=ccs.K=4.MS=2.bins=10.KX=1");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels," ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doLabelZStyCCSseveral()
{
  RT_Test* rt;
  rt = new RT_Test("doLabelZStyCCSseveral");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--summaries");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm-p.lowprecision=true.transform=ccs.K=4.MS=2.bins=10.KX=1");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels," ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm-p.confidences.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doStyCCS_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyCCS_2_20_1_0");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels," ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doStyCCS_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyCCS_2_20_1_0");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels," ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doStyCCS_2_40_0_9_0_033()
{
   RT_Test* rt;
  rt = new RT_Test("doStyCCS_2_40_0_9_0_033");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=40.lowprecision=true.het-mult=0.9.transform=ccs.K=4.MS=2 --dm-thresh 0.33");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels," ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doStyCCS_2_40_0_9_0_033()
{
   RT_Test* rt;
  rt = new RT_Test("doStyCCS_2_40_0_9_0_033");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=40.lowprecision=true.het-mult=0.9.transform=ccs.K=4.MS=2 --dm-thresh 0.33");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels," ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doSpf()
{
  RT_Test* rt;
  rt = new RT_Test("doSpf");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("--chrX-snps", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.chrx");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_idata_pgeno_mapping250ksty}/cel-files.txt");

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doSpf()
{
  RT_Test* rt;
  rt = new RT_Test("doSpf");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("--chrX-snps", "${gold_lib_mapping250ksty}/Mapping250K_Sty.chrx");
  cmd->addArg("--spf-file", "${gold_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_probesetgenotype_mapping250ksty}/cel-files.txt");

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doChpFiles()
{
  RT_Test* rt;
  rt = new RT_Test("doChpFiles");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--xda-chp-output");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--chrX-snps", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.chrx");
  cmd->addArg("-spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_idata_pgeno_mapping250ksty}/cel-files.txt");

  const char *chpFiles[] = { "NA06985_B01_Sty_Plate1.brlmm","NA06991_B03_Sty_Plate1.brlmm","NA06993_B02_Sty_Plate1.brlmm",
			     "NA06994_A11_Sty_Plate1.brlmm","NA07000_A10_Sty_Plate1.brlmm","NA07019_A09_Sty_Plate1.brlmm",
			     "NA07022_A08_Sty_Plate1.brlmm","NA07029_A12_Sty_Plate1.brlmm","NA07034_B05_Sty_Plate1.brlmm",
			     "NA07048_B06_Sty_Plate1.brlmm","NA07055_B04_Sty_Plate1.brlmm","NA07056_A07_Sty_Plate1.brlmm",
			     "NA07345_B10_Sty_Plate1.brlmm","NA07348_B12_Sty_Plate1.brlmm","NA07357_B11_Sty_Plate1.brlmm",
			     "NA10846_A06_Sty_Plate1.brlmm","NA10847_A03_Sty_Plate1.brlmm","NA10851_B09_Sty_Plate1.brlmm",
			     "NA10854_C09_Sty_Plate1.brlmm","NA10855_C12_Sty_Plate1.brlmm",
			     NULL
  };

  //Append appropriate paths to CHP files
  std::vector<std::string> gold,gen,goldX,genX;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_pgeno_mapping250ksty_testname}/chp/","${chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/chp/", "${chpsuffix}");
  goldX = Util::addPrefixSuffix(chpFiles, "${gold_idata_pgeno_mapping250ksty_testname}/cc-chp/", "${chpsuffix}");
  genX = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/","${chpsuffix}");

  //CHP CHECK
  //Generate Check objects for ChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_chp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  
  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_calvinchp}");
  check2->setGenFiles(genX);
  check2->setGoldFiles(goldX);

  //MIXED FILE CHECK
  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_mixedfile}");
  check3->addArg("--gen", "${out_testname}/brlmm.report.txt");
  check3->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/brlmm.report.txt");
  check3->addArg("--epsilon", "0.00001");
  check3->addArg("--skipLines", "0");
  check3->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/brlmm.calls.txt");
  check4->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/brlmm.calls.txt");
  check4->addArg("--epsilon", "0.0000001");
  check4->addArg("--rowSkip", "1");
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_testname}/brlmm.confidences.txt");
  check5->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/brlmm.confidences.txt");
  check5->addArg("--epsilon", "0.00001");
  check5->addArg("--rowSkip", "1");
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doChpFiles()
{
  RT_Test* rt;
  rt = new RT_Test("doChpFiles");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--xda-chp-output");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("--table-output");
  cmd->addArg("--chrX-snps", "${gold_lib_mapping250ksty}/Mapping250K_Sty.chrx");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_probesetgenotype_mapping250ksty}/cel-files.txt");

  const char *chpFiles[] = { "NA06985_B01_Sty_Plate1.brlmm","NA06991_B03_Sty_Plate1.brlmm","NA06993_B02_Sty_Plate1.brlmm",
			     "NA06994_A11_Sty_Plate1.brlmm","NA07000_A10_Sty_Plate1.brlmm","NA07019_A09_Sty_Plate1.brlmm",
			     "NA07022_A08_Sty_Plate1.brlmm","NA07029_A12_Sty_Plate1.brlmm","NA07034_B05_Sty_Plate1.brlmm",
			     "NA07048_B06_Sty_Plate1.brlmm","NA07055_B04_Sty_Plate1.brlmm","NA07056_A07_Sty_Plate1.brlmm",
			     "NA07345_B10_Sty_Plate1.brlmm","NA07348_B12_Sty_Plate1.brlmm","NA07357_B11_Sty_Plate1.brlmm",
			     "NA10846_A06_Sty_Plate1.brlmm","NA10847_A03_Sty_Plate1.brlmm","NA10851_B09_Sty_Plate1.brlmm",
			     "NA10854_C09_Sty_Plate1.brlmm","NA10855_C12_Sty_Plate1.brlmm",
			     NULL
  };

  //Append appropriate paths to CHP files
  std::vector<std::string> gold,gen,goldX,genX;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetgenotype_mapping250ksty_testname}/chp/","${chpsuffix}");
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/chp/", "${chpsuffix}");
  goldX = Util::addPrefixSuffix(chpFiles, "${gold_probesetgenotype_mapping250ksty_testname}/cc-chp/", "${chpsuffix}");
  genX = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/","${chpsuffix}");

  //CHP CHECK
  //Generate Check objects for ChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_chp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  
  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_calvinchp}");
  check2->setGenFiles(genX);
  check2->setGoldFiles(goldX);

  //MIXED FILE CHECK
  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_mixedfile}");
  check3->addArg("--gen", "${out_testname}/brlmm.report.txt");
  check3->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/brlmm.report.txt");
  check3->addArg("--epsilon", "0.00001");
  check3->addArg("--skipLines", "0");
  check3->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/brlmm.calls.txt");
  check4->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/brlmm.calls.txt");
  check4->addArg("--epsilon", "0.0000001");
  check4->addArg("--rowSkip", "1");
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_testname}/brlmm.confidences.txt");
  check5->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/brlmm.confidences.txt");
  check5->addArg("--epsilon", "0.00001");
  check5->addArg("--rowSkip", "1");
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doFileStyCCS_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doFileStyCCS_2_20_1_0");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_idata_pgeno_mapping250ksty}/cel-files.txt");

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doFileStyCCS_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doFileStyCCS_2_20_1_0");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${gold_probesetgenotype_mapping250ksty}/cel-files.txt");

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doWritePriorFileStyCCS_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doWritePriorFileStyCCS_2_20_1_0");
  //rt->setIntegrationTestBool(true);  
  RT_Test* rt1;
  rt1 = new RT_Test("doWritePriorFileStyCCS_2_20_1_0-part1");
//rt1->setIntegrationTestBool(true);
  RT_Test* rt2;
  rt2 = new RT_Test("doWritePriorFileStyCCS_2_20_1_0-part2");
//rt2->setIntegrationTestBool(true);

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_genotype}");
  cmd1->addArg("--table-output");
  cmd1->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd1->addArg("--prior-size", "1000");
  cmd1->addArg("--list-sample");
  cmd1->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd1->addArg("--no-gender-force");
  cmd1->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd1->addArg("-o", "test-generated/doWritePriorFileStyCCS_2_20_1_0");
  cmd1->addArg("--cel-files", "${gold_idata_pgeno_mapping250ksty}/cel-files.txt");
  cmd1->addArg("--write-prior");

  //MATRIX CHECK
  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty}/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt1->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty}/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_genotype}");
  cmd2->addArg("--table-output");
  cmd2->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd2->addArg("--prior-size", "1000");
  cmd2->addArg("--list-sample");
  cmd2->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd2->addArg("--no-gender-force");
  cmd2->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd2->addArg("-o", "test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2");
  cmd2->addArg("--cel-files", "${gold_idata_pgeno_mapping250ksty}/cel-files.txt");
  cmd2->addArg("--read-priors-brlmm", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.prior.txt");

  //MATRIX CHECK
  RT_Check* check3 = rt2->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2/quant-norm.pm-only.brlmm.calls.txt");
  check3->addArg("--gold", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt");
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1");
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "true");
  check3->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check4 = rt2->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2/quant-norm.pm-only.brlmm.confidences.txt");
  check4->addArg("--gold", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt");
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1");
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "true");
  check4->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doWritePriorFileStyCCS_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doWritePriorFileStyCCS_2_20_1_0");
  RT_Test* rt1;
  rt1 = new RT_Test("doWritePriorFileStyCCS_2_20_1_0-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("doWritePriorFileStyCCS_2_20_1_0-part2");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_genotype}");
  cmd1->addArg("--table-output");
  cmd1->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd1->addArg("--list-sample");
  cmd1->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd1->addArg("--no-gender-force");
  cmd1->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd1->addArg("-o", "test-generated/doWritePriorFileStyCCS_2_20_1_0");
  cmd1->addArg("--cel-files", "${gold_probesetgenotype_mapping250ksty}/cel-files.txt");
  cmd1->addArg("--write-prior");

  //MATRIX CHECK
  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty}/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt1->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty}/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_genotype}");
  cmd2->addArg("--table-output");
  cmd2->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd2->addArg("--list-sample");
  cmd2->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd2->addArg("--no-gender-force");
  cmd2->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd2->addArg("-o", "test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2");
  cmd2->addArg("--cel-files", "${gold_probesetgenotype_mapping250ksty}/cel-files.txt");
  cmd2->addArg("--read-priors-brlmm", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.prior.txt");

  //MATRIX CHECK
  RT_Check* check3 = rt2->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2/quant-norm.pm-only.brlmm.calls.txt");
  check3->addArg("--gold", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.calls.txt");
  check3->addArg("--epsilon", "0.0000001");
  check3->addArg("--rowSkip", "1");
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "true");
  check3->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check4 = rt2->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "test-generated/doWritePriorFileStyCCS_2_20_1_0/phase2/quant-norm.pm-only.brlmm.confidences.txt");
  check4->addArg("--gold", "test-generated/doWritePriorFileStyCCS_2_20_1_0/quant-norm.pm-only.brlmm.confidences.txt");
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1");
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "true");
  check4->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doReadGenotypesIn()
{
  RT_Test* rt;
  rt = new RT_Test("doReadGenotypesIn");
  //rt->setIntegrationTestBool(true);
  RT_Test* rt1;
  rt1 = new RT_Test("doReadGenotypesIn-part1");
//rt1->setIntegrationTestBool(true);
  RT_Test* rt2;
  rt2 = new RT_Test("doReadGenotypesIn-part2");
//rt2->setIntegrationTestBool(true);

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_genotype}");
  cmd1->addArg("--table-output ");
  cmd1->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd1->addArg("--prior-size", "1000");
  cmd1->addArg("--list-sample");
  cmd1->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd1->addArg("--no-gender-force");
  cmd1->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd1->addArg("-o", "test-generated/doReadGenotypesIn");
  cmd1->addArg("--cel-files", "${gold_idata_pgeno_mapping250ksty}/cel-files.txt");
  cmd1->addArg("--dm-out");
  cmd1->addArg("--dm-thresh", ".33");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_genotype}");
  cmd2->addArg("--table-output");
  cmd2->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd2->addArg("--prior-size", "1000");
  cmd2->addArg("--list-sample");
  cmd2->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd2->addArg("--no-gender-force");
  cmd2->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd2->addArg("-o", "test-generated/doReadGenotypesIn/phase2");
  cmd2->addArg("--cel-files", "${gold_idata_pgeno_mapping250ksty}/cel-files.txt");
  cmd2->addArg("--genotypes", "test-generated/doReadGenotypesIn/dm.calls.txt");

  //MATRIX CHECK
  RT_Check* check1 = rt2->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "test-generated/doReadGenotypesIn/phase2/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "test-generated/doReadGenotypesIn/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "test-generated/doReadGenotypesIn/phase2/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "test-generated/doReadGenotypesIn/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doReadGenotypesIn()
{
  RT_Test* rt;
  rt = new RT_Test("doReadGenotypesIn");
  RT_Test* rt1;
  rt1 = new RT_Test("doReadGenotypesIn-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("doReadGenotypesIn-part2");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_genotype}");
  cmd1->addArg("--table-output ");
  cmd1->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd1->addArg("--list-sample");
  cmd1->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd1->addArg("--no-gender-force");
  cmd1->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd1->addArg("-o", "test-generated/doReadGenotypesIn");
  cmd1->addArg("--cel-files", "${gold_probesetgenotype_mapping250ksty}/cel-files.txt");
  cmd1->addArg("--dm-out ");
  cmd1->addArg("--dm-thresh", ".33");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_genotype}");
  cmd2->addArg("--table-output");
  cmd2->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd2->addArg("--list-sample");
  cmd2->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ccs.K=4.MS=2");
  cmd2->addArg("--no-gender-force");
  cmd2->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd2->addArg("-o", "test-generated/doReadGenotypesIn/phase2");
  cmd2->addArg("--cel-files", "${gold_probesetgenotype_mapping250ksty}/cel-files.txt");
  cmd2->addArg("--genotypes", "test-generated/doReadGenotypesIn/dm.calls.txt");

  //MATRIX CHECK
  RT_Check* check1 = rt2->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "test-generated/doReadGenotypesIn/phase2/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "test-generated/doReadGenotypesIn/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "test-generated/doReadGenotypesIn/phase2/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "test-generated/doReadGenotypesIn/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doStyRVT_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyRVT_2_20_1_0");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=rvt.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doStyRVT_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyRVT_2_20_1_0");

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=rvt.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doStyCES_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyCES_2_20_1_0");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ces.K=1.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doStyCES_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyCES_2_20_1_0");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.lowprecision=true.het-mult=1.transform=ces.K=1.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doStyMva_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyMva_2_20_1_0");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=1.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doStyMva_2_20_1_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyMva_2_20_1_0");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=1.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doStyMva_2_20_1_2()
{
  RT_Test* rt;
  rt = new RT_Test("doStyMva_2_20_1_2");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=1.iterations=1.iter-thresh=.6.MS=2 --no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doStyMva_2_20_1_2()
{
  RT_Test* rt;
  rt = new RT_Test("doStyMva_2_20_1_2");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=1.iterations=1.iter-thresh=.6.MS=2 --no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doStyMva_2_20_0_8_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyMva_2_20_0_8_0");
  //rt->setIntegrationTestBool(true);

  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_idata_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--prior-size", "1000");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=.8.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("--spf-file", "${gold_idata_lib_mapping250ksty}/Mapping250K_Sty.spf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_idata_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_idata_pgeno_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doStyMva_2_20_0_8_0()
{
  RT_Test* rt;
  rt = new RT_Test("doStyMva_2_20_0_8_0");
  const char *sty20cels[] = {
    "NA06985_B01_Sty_Plate1","NA06991_B03_Sty_Plate1","NA06993_B02_Sty_Plate1",
    "NA06994_A11_Sty_Plate1","NA07000_A10_Sty_Plate1","NA07019_A09_Sty_Plate1",
    "NA07022_A08_Sty_Plate1","NA07029_A12_Sty_Plate1","NA07034_B05_Sty_Plate1",
    "NA07048_B06_Sty_Plate1","NA07055_B04_Sty_Plate1","NA07056_A07_Sty_Plate1",
    "NA07345_B10_Sty_Plate1","NA07348_B12_Sty_Plate1","NA07357_B11_Sty_Plate1",
    "NA10846_A06_Sty_Plate1","NA10847_A03_Sty_Plate1","NA10851_B09_Sty_Plate1",
    "NA10854_C09_Sty_Plate1","NA10855_C12_Sty_Plate1",
    NULL
  };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_genotype}");
  cmd->addArg("--table-output");
  cmd->addArg("-s", "${gold_lib_mapping250ksty}/set1snps.txt");
  cmd->addArg("--list-sample");
  cmd->addArg("-a", "quant-norm.sketch=50000,pm-only,brlmm.prior-weight=20.transform=mva.lowprecision=true.het-mult=.8.MS=2");
  cmd->addArg("--no-gender-force");
  cmd->addArg("-c", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("-o", "${out_testname}");
  std::vector<std::string> fullSty20celsFiles = Util::addPrefixSuffix(sty20cels, " ${gold_cel_mapping250ksty}/","${celsuffix}");
  for(std::vector<std::string>::iterator i = fullSty20celsFiles.begin();i < fullSty20celsFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 //MATRIX CHECK
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.calls.txt");
  check1->addArg("--epsilon", "0.0000001");
  check1->addArg("--rowSkip", "1");
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  //MATRIX CHECK
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--gold", "${gold_probesetgenotype_mapping250ksty_testname}/quant-norm.pm-only.brlmm.confidences.txt");
  check2->addArg("--epsilon", "0.00001");
  check2->addArg("--rowSkip", "1");
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}
