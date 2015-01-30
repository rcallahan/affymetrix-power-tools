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

RT_Test* qt_doRmaTissueSketchSuppliedTest();
RT_Test* qt_doPlierPrecompTissueTest();
RT_Test* qt_doRmaNvissaTest();
RT_Test* qt_doRmaNvissaReadWriteFeatureEffectsTest(); 
RT_Test* qt_doRmaNvissaSketchTest();
RT_Test* qt_doRmaTissueTest();
RT_Test* qt_doRmaTissueSpfTest();
RT_Test* qt_doRmaTissueTestCC();
RT_Test* qt_doRmaTissueSketchTest();
RT_Test* qt_doPlierGcBgTissueTest();
RT_Test* qt_doPlierMMTissueSketchTest();
RT_Test* qt_doPlierWtaRefSeqTest();
RT_Test* qt_doHumanGeneTest();
RT_Test* qt_doHumanGeneSpfTest();
RT_Test* qt_doHumanGeneKillListTest();
RT_Test* qt_doSNP6CopyNumberWorkflow1a();
RT_Test* qt_doSNP6CopyNumberWorkflow1b();
RT_Test* qt_doSNP6CopyNumberWorkflow2();
RT_Test* qt_doPlierMMTissueMedianNormTest();
RT_Test* qt_doDabgSubsetTest();
RT_Test* qt_doDabgU133Test();
RT_Test* doRmaTissueSketchSuppliedTest();
RT_Test* doPlierPrecompTissueTest();
RT_Test* doRmaNvissaTest();
RT_Test* doRmaNvissaReadWriteFeatureEffectsTest(); 
RT_Test* doRmaNvissaSketchTest();
RT_Test* doRmaTissueTest();
RT_Test* doRmaTissueSpfTest();
RT_Test* doRmaTissueTestCC();
RT_Test* doRmaTissueSketchTest();
RT_Test* doPlierGcBgTissueTest();
RT_Test* doPlierMMTissueSketchTest();
RT_Test* doPlierWtaRefSeqTest();
RT_Test* doHumanGeneTest();
RT_Test* doHumanGeneSpfTest();
RT_Test* doHumanGeneKillListTest();
RT_Test* doSNP6CopyNumberWorkflow1a();
RT_Test* doSNP6CopyNumberWorkflow1b();
RT_Test* doSNP6CopyNumberWorkflow2();
RT_Test* doPlierMMTissueMedianNormTest();
RT_Test* doDabgSubsetTest();
RT_Test* doDabgU133Test();



const char* tissueCels[]={
  "heart-rep1",     "heart-rep2",     "heart-rep3", 
  "hela-rep1",      "hela-rep2",      "hela-rep3",
  "pancrease-rep1", "pancrease-rep2", "pancrease-rep3",
  NULL
};

const char* tissueCcCels[]={
  "heart-rep1",     "heart-rep2",     "heart-rep3",
  "hela-rep1",      "hela-rep2",      "hela-rep3",
  "pancrease-rep1", "pancrease-rep2", "pancrease-rep3",
  NULL
};

const char* nvissaCels[]={
  "4221B_a01",	 "4221B_b01",	 "4221B_c01",			
  "4221H_a02",	 "4221H_b02",	 "4221H_c02",			
  "4221HL_a04",	 "4221HL_b04",	 "4221HL_c04",			
  "4221P_a03",	 "4221P_b03",	 "4221P_c03",			
  NULL
};

const char* huexCels[]={
  "huex_wta_cb_A",	 "huex_wta_cb_B",	 "huex_wta_cb_C",	
  "huex_wta_heart_A",	 "huex_wta_heart_B",	 "huex_wta_heart_C",	
  "huex_wta_muscle_A",	 "huex_wta_muscle_B",	 "huex_wta_muscle_C",		
  "huex_wta_testes_A",	 "huex_wta_testes_B",	 "huex_wta_testes_C",	
  NULL
};


///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

  RT_Main* rMain;
  rMain = new RT_Main("apt-rt-probeset-summarize");

  if ( !Fs::dirExists("test-generated") ) {
    Fs::mkdir("test-generated", false);
  }

  //Store names of tests for Help Message
  //Store it in our RT_Main Instance
  std::vector<std::string> testNameList;
  testNameList.push_back("qt_doRmaTissueSketchSuppliedTest");
  testNameList.push_back("qt_doPlierPrecompTissueTest");
  testNameList.push_back("qt_doRmaNvissaTest");
  testNameList.push_back("qt_doRmaNvissaReadWriteFeatureEffectsTest");
  testNameList.push_back("qt_doRmaNvissaSketchTest");
  testNameList.push_back("qt_doRmaTissueTest");
  testNameList.push_back("qt_doRmaTissueSpfTest");
  testNameList.push_back("qt_doRmaTissueTestCC");
  testNameList.push_back("qt_doRmaTissueSketchTest");
  testNameList.push_back("qt_doPlierGcBgTissueTest");
  testNameList.push_back("qt_doPlierMMTissueSketchTest");
  testNameList.push_back("qt_doPlierWtaRefSeqTest");
  testNameList.push_back("qt_doHumanGeneTest");
  testNameList.push_back("qt_doHumanGeneSpfTest");
  testNameList.push_back("qt_doHumanGeneKillListTest");
  testNameList.push_back("qt_doSNP6CopyNumberWorkflow1a");
  testNameList.push_back("qt_doSNP6CopyNumberWorkflow1b");
  testNameList.push_back("qt_doSNP6CopyNumberWorkflow2");
  testNameList.push_back("qt_doPlierMMTissueMedianNormTest");
  testNameList.push_back("qt_doDabgSubsetTest");
  testNameList.push_back("qt_doDabgU133Test");
  testNameList.push_back("doRmaTissueSketchSuppliedTest");
  testNameList.push_back("doPlierPrecompTissueTest");
  testNameList.push_back("doRmaNvissaTest");
  testNameList.push_back("doRmaNvissaReadWriteFeatureEffectsTest");
  testNameList.push_back("doRmaNvissaSketchTest");
  testNameList.push_back("doRmaTissueTest");
  testNameList.push_back("doRmaTissueSpfTest");
  testNameList.push_back("doRmaTissueTestCC");
  testNameList.push_back("doRmaTissueSketchTest");
  testNameList.push_back("doPlierGcBgTissueTest");
  testNameList.push_back("doPlierMMTissueSketchTest");
  testNameList.push_back("doPlierWtaRefSeqTest");
  testNameList.push_back("doHumanGeneTest");
  testNameList.push_back("doHumanGeneSpfTest");
  testNameList.push_back("doHumanGeneKillListTest");
  testNameList.push_back("doSNP6CopyNumberWorkflow1a");
  testNameList.push_back("doSNP6CopyNumberWorkflow1b");
  testNameList.push_back("doSNP6CopyNumberWorkflow2");
  testNameList.push_back("doPlierMMTissueMedianNormTest");
  testNameList.push_back("doDabgSubsetTest");
  testNameList.push_back("doDabgU133Test");
  rMain->setTestNameList(testNameList);

  rMain->setVar("qt_celFileList", "${gold_dir}/idata/p-sum/HG-U133_Plus_2/cel-files.txt");
  rMain->setVar("qt_celFileList2", "${gold_dir}/idata/p-sum/HuGene-1_0-st-v1/cel-files.txt");
  rMain->setVar("celFileList", "${gold_dir}/chipstream/probeset-summarize/HG-U133_Plus_2/cel-files.txt");
  rMain->setVar("celFileList2", "${gold_dir}/chipstream/probeset-summarize/HuGene-1_0-st-v1/cel-files.txt");
  rMain->setVar("gold_idata_cel_hgu133plus2", "${gold_dir}/idata/cel/HG-U133_Plus_2");
  rMain->setVar("gold_cel_hgu133plus2", "${gold_dir}/cel/HG-U133_Plus_2");
  rMain->setVar("gold_cel_huex10stv2", "${gold_dir}/cel/HuEx-1_0-st-v2");
  rMain->setVar("sdk_top", "../../..");
  rMain->setVar("gold_lib_hgu133plus2", "${gold_dir}/lib/HG-U133_Plus_2");
  rMain->setVar("gold_probesetsummarize_hgu133plus2", "${gold_dir}/chipstream/probeset-summarize/HG-U133_Plus_2");
  rMain->setVar("gold_probesetsummarize_hgu133plus2_testname", "${gold_dir}/chipstream/probeset-summarize/HG-U133_Plus_2/${testname}");
  rMain->setVar("gold_lib_huex10stv2", "${gold_dir}/lib/HuEx-1_0-st-v2");
  rMain->setVar("gold_lib_hugene10stv1", "${gold_dir}/lib/HuGene-1_0-st-v1");
  rMain->setVar("gold_idata_lib_hugene10stv1", "${gold_dir}/idata/lib/HuGene-1_0-st-v1");
  rMain->setVar("gold_probesetsummarize_hugene10stv1", "${gold_dir}/chipstream/probeset-summarize/HuGene-1_0-st-v1");
  rMain->setVar("gold_probesetsummarize_hugene10stv1_testname", "${gold_dir}/chipstream/probeset-summarize/HuGene-1_0-st-v1/${testname}");
  rMain->setVar("gold_cel_genomewidesnp6", "${gold_dir}/cel/GenomeWideSNP_6");
  rMain->setVar("gold_lib_genomewidesnp6", "${gold_dir}/lib/GenomeWideSNP_6");
  rMain->setVar("gold_idata_cel_genomewidesnp6", "${gold_dir}/idata/cel/GenomeWideSNP_6");
  rMain->setVar("gold_idata_cel_genomewidesnp6cn", "${gold_dir}/idata/cel/GenomeWideSNP_6_cn");
  rMain->setVar("gold_idata_lib_genomewidesnp6cn", "${gold_dir}/idata/lib/GenomeWideSNP_6_cn");
  rMain->setVar("gold_probesetsummarize_genomewidesnp6_testname", "${gold_dir}/chipstream/probeset-summarize/GenomeWideSNP_6/${testname}");
  rMain->setVar("gold_probesetsummarize_huex10stv2", "${gold_dir}/chipstream/probeset-summarize/HuEx-1_0-st-v2");
  rMain->setVar("gold_probesetsummarize_huex10stv2_testname", "${gold_dir}/chipstream/probeset-summarize/HuEx-1_0-st-v2/${testname}");
  rMain->setVar("gold_lib_hugene10stv1", "${gold_dir}/lib/HuGene-1_0-st-v1");
  rMain->setVar("celsuffixLower", ".cel");
  rMain->setVar("celsuffixUp", ".CEL");
  rMain->setVar("agcccelsuffix", ".agcc.cel");
  rMain->setVar("gold_idata_psum_hgu133plus2", "${gold_dir}/idata/p-sum/HG-U133_Plus_2");
  rMain->setVar("gold_idata_lib_hgu133plus2", "${gold_dir}/idata/lib/HG-U133_Plus_2");
  rMain->setVar("gold_idata_psum_testname", "${gold_dir}/idata/p-sum/${testname}");
  rMain->setVar("gold_idata_lib", "${gold_dir}/idata/lib");
  rMain->setVar("gold_idata_psum", "${gold_dir}/idata/p-sum");
  rMain->setVar("out_doRS", "${out_dir}/doRS");
  rMain->setVar("gold_idata_lib_huex10stv2", "${gold_dir}/idata/lib/HuEx-1_0-st-v2");
  rMain->setVar("gold_idata_cel_huex10stv2", "${gold_dir}/idata/cel/HuEx-1_0-st-v2");
    


  rMain->parseArgv(argc, argv);

  if(argc < 1)
    {
      cout<<"Error: Incorrect number of arguments given"<<endl;
      exit(1);
    }
  //They want us to run all default tests
  else if(rMain->getVarVal("allintegration") != "")
    {
      rMain->addTest(qt_doRmaTissueSketchSuppliedTest());
      rMain->addTest(qt_doPlierPrecompTissueTest());
      rMain->addTest(qt_doRmaNvissaTest());
      rMain->addTest(qt_doRmaNvissaReadWriteFeatureEffectsTest()); 
      rMain->addTest(qt_doRmaNvissaSketchTest());
      rMain->addTest(qt_doRmaTissueTest());
      rMain->addTest(qt_doRmaTissueSpfTest());
      rMain->addTest(qt_doRmaTissueTestCC());
      rMain->addTest(qt_doRmaTissueSketchTest());
      rMain->addTest(qt_doPlierGcBgTissueTest());
      rMain->addTest(qt_doPlierMMTissueSketchTest());
      rMain->addTest(qt_doPlierWtaRefSeqTest());
      rMain->addTest(qt_doHumanGeneTest());
      rMain->addTest(qt_doHumanGeneSpfTest());
      rMain->addTest(qt_doHumanGeneKillListTest());
      rMain->addTest(qt_doSNP6CopyNumberWorkflow1a());
      rMain->addTest(qt_doSNP6CopyNumberWorkflow1b());
      rMain->addTest(qt_doSNP6CopyNumberWorkflow2());
      rMain->addTest(qt_doPlierMMTissueMedianNormTest());
      rMain->addTest(qt_doDabgSubsetTest());
      rMain->addTest(qt_doDabgU133Test());
    }
  //They want us to run all default tests
  else if(rMain->getVarVal("allnormal") != "")
    {
      rMain->addTest(doRmaTissueSketchSuppliedTest());
      rMain->addTest(doPlierPrecompTissueTest());
      rMain->addTest(doRmaNvissaTest());
      rMain->addTest(doRmaNvissaReadWriteFeatureEffectsTest()); 
      rMain->addTest(doRmaNvissaSketchTest());
      rMain->addTest(doRmaTissueTest());
      rMain->addTest(doRmaTissueSpfTest());
      rMain->addTest(doRmaTissueTestCC());
      rMain->addTest(doRmaTissueSketchTest());
      rMain->addTest(doPlierGcBgTissueTest());
      rMain->addTest(doPlierMMTissueSketchTest());
      rMain->addTest(doPlierWtaRefSeqTest());
      rMain->addTest(doHumanGeneTest());
      rMain->addTest(doHumanGeneSpfTest());
      rMain->addTest(doHumanGeneKillListTest());
      rMain->addTest(doSNP6CopyNumberWorkflow1a());
      rMain->addTest(doSNP6CopyNumberWorkflow1b());
      rMain->addTest(doSNP6CopyNumberWorkflow2());
      rMain->addTest(doPlierMMTissueMedianNormTest());
      rMain->addTest(doDabgSubsetTest());
      rMain->addTest(doDabgU133Test());
    }
  //They want us to run all default tests
  else if(argc == 1 || rMain->getVarVal("all") != "")
    {
      //RUN ALL TESTS OMG
      cout<<"Running All Default Regression Tests!"<<endl;
      rMain->addTest(qt_doRmaTissueSketchSuppliedTest());
      rMain->addTest(qt_doPlierPrecompTissueTest());
      rMain->addTest(qt_doRmaNvissaTest());
      rMain->addTest(qt_doRmaNvissaReadWriteFeatureEffectsTest()); 
      rMain->addTest(qt_doRmaNvissaSketchTest());
      rMain->addTest(qt_doRmaTissueTest());
      rMain->addTest(qt_doRmaTissueSpfTest());
      rMain->addTest(qt_doRmaTissueTestCC());
      rMain->addTest(qt_doRmaTissueSketchTest());
      rMain->addTest(qt_doPlierGcBgTissueTest());
      rMain->addTest(qt_doPlierMMTissueSketchTest());
      rMain->addTest(qt_doPlierWtaRefSeqTest());
      rMain->addTest(qt_doHumanGeneTest());
      rMain->addTest(qt_doHumanGeneSpfTest());
      rMain->addTest(qt_doHumanGeneKillListTest());
      rMain->addTest(qt_doSNP6CopyNumberWorkflow1a());
      rMain->addTest(qt_doSNP6CopyNumberWorkflow1b());
      rMain->addTest(qt_doSNP6CopyNumberWorkflow2());
      rMain->addTest(qt_doPlierMMTissueMedianNormTest());
      rMain->addTest(qt_doDabgSubsetTest());
      rMain->addTest(qt_doDabgU133Test());
      rMain->addTest(doRmaTissueSketchSuppliedTest());
      rMain->addTest(doPlierPrecompTissueTest());
      rMain->addTest(doRmaNvissaTest());
      rMain->addTest(doRmaNvissaReadWriteFeatureEffectsTest()); 
      rMain->addTest(doRmaNvissaSketchTest());
      rMain->addTest(doRmaTissueTest());
      rMain->addTest(doRmaTissueSpfTest());
      rMain->addTest(doRmaTissueTestCC());
      rMain->addTest(doRmaTissueSketchTest());
      rMain->addTest(doPlierGcBgTissueTest());
      rMain->addTest(doPlierMMTissueSketchTest());
      rMain->addTest(doPlierWtaRefSeqTest());
      rMain->addTest(doHumanGeneTest());
      rMain->addTest(doHumanGeneSpfTest());
      rMain->addTest(doHumanGeneKillListTest());
      rMain->addTest(doSNP6CopyNumberWorkflow1a());
      rMain->addTest(doSNP6CopyNumberWorkflow1b());
      rMain->addTest(doSNP6CopyNumberWorkflow2());
      rMain->addTest(doPlierMMTissueMedianNormTest());
      rMain->addTest(doDabgSubsetTest());
      rMain->addTest(doDabgU133Test());
    }
  else
    {
      cout<<"Running Built-In Regression Tests"<<endl;
	 
      if(rMain->getVarVal("qt_doRmaTissueSketchSuppliedTest") != "")
	{
	  rMain->addTest(qt_doRmaTissueSketchSuppliedTest());
	}
      if(rMain->getVarVal("qt_doPlierPrecompTissueTest") != "")
	{
	  rMain->addTest(qt_doPlierPrecompTissueTest());
	}
      if(rMain->getVarVal("qt_doRmaNvissaTest") != "")
	{
	  rMain->addTest(qt_doRmaNvissaTest());
	}
      if(rMain->getVarVal("qt_doRmaNvissaReadWriteFeatureEffectsTest") != "")
	{
	  rMain->addTest(qt_doRmaNvissaReadWriteFeatureEffectsTest());
	}
      if(rMain->getVarVal("qt_doRmaNvissaSketchTest") != "")
	{
	  rMain->addTest(qt_doRmaNvissaSketchTest());
	}
      if(rMain->getVarVal("qt_doRmaTissueTest") != "")
	{
	  rMain->addTest(qt_doRmaTissueTest());
	}
      if(rMain->getVarVal("qt_doRmaTissueSpfTest") != "")
	{
	  rMain->addTest(qt_doRmaTissueSpfTest());
	}
      if(rMain->getVarVal("qt_doRmaTissueTestCC") != "")
	{
	  rMain->addTest(qt_doRmaTissueTestCC());
	}
      if(rMain->getVarVal("qt_doRmaTissueSketchTest") != "")
	{
	  rMain->addTest(qt_doRmaTissueSketchTest());
	}
      if(rMain->getVarVal("qt_doPlierGcBgTissueTest") != "")
	{
	  rMain->addTest(qt_doPlierGcBgTissueTest());
	}
      if(rMain->getVarVal("qt_doPlierMMTissueSketchTest") != "")
	{
	  rMain->addTest(qt_doPlierMMTissueSketchTest());
	}
      if(rMain->getVarVal("qt_doPlierWtaRefSeqTest") != "")
	{
	  rMain->addTest(qt_doPlierWtaRefSeqTest());
	}
      if(rMain->getVarVal("qt_doHumanGeneTest") != "")
	{
	  rMain->addTest(qt_doHumanGeneTest());
	}
      if(rMain->getVarVal("qt_doHumanGeneSpfTest") != "")
	{
	  rMain->addTest(qt_doHumanGeneSpfTest());
	}
      if(rMain->getVarVal("qt_doHumanGeneKillListTest") != "")
	{
	  rMain->addTest(qt_doHumanGeneKillListTest());
	}
      if(rMain->getVarVal("qt_doSNP6CopyNumberWorkflow1a") != "")
	{
	  rMain->addTest(qt_doSNP6CopyNumberWorkflow1a());
	}
      if(rMain->getVarVal("qt_doSNP6CopyNumberWorkflow1b") != "")
	{
	  rMain->addTest(qt_doSNP6CopyNumberWorkflow1b());
	}
      if(rMain->getVarVal("qt_doSNP6CopyNumberWorkflow2") != "")
	{
	  rMain->addTest(qt_doSNP6CopyNumberWorkflow2());
	}
      if(rMain->getVarVal("qt_doPlierMMTissueMedianNormTest") != "")
	{
	  rMain->addTest(qt_doPlierMMTissueMedianNormTest());
	}
      if(rMain->getVarVal("qt_doDabgSubsetTest") != "")
	{
	  rMain->addTest(qt_doDabgSubsetTest());
	}
      if(rMain->getVarVal("qt_doDabgU133Test") != "")
	{
	  rMain->addTest(qt_doDabgU133Test());
	}
      if(rMain->getVarVal("doRmaTissueSketchSuppliedTest") != "")
	{
	  rMain->addTest(doRmaTissueSketchSuppliedTest());
	}
      if(rMain->getVarVal("doPlierPrecompTissueTest") != "")
	{
	  rMain->addTest(doPlierPrecompTissueTest());
	}
      if(rMain->getVarVal("doRmaNvissaTest") != "")
	{
	  rMain->addTest(doRmaNvissaTest());
	}
      if(rMain->getVarVal("doRmaNvissaReadWriteFeatureEffectsTest") != "")
	{
	  rMain->addTest(doRmaNvissaReadWriteFeatureEffectsTest());
	}
      if(rMain->getVarVal("doRmaNvissaSketchTest") != "")
	{
	  rMain->addTest(doRmaNvissaSketchTest());
	}
      if(rMain->getVarVal("doRmaTissueTest") != "")
	{
	  rMain->addTest(doRmaTissueTest());
	}
      if(rMain->getVarVal("doRmaTissueSpfTest") != "")
	{
	  rMain->addTest(doRmaTissueSpfTest());
	}
      if(rMain->getVarVal("doRmaTissueTestCC") != "")
	{
	  rMain->addTest(doRmaTissueTestCC());
	}
      if(rMain->getVarVal("doRmaTissueSketchTest") != "")
	{
	  rMain->addTest(doRmaTissueSketchTest());
	}
      if(rMain->getVarVal("doPlierGcBgTissueTest") != "")
	{
	  rMain->addTest(doPlierGcBgTissueTest());
	}
      if(rMain->getVarVal("doPlierMMTissueSketchTest") != "")
	{
	  rMain->addTest(doPlierMMTissueSketchTest());
	}
      if(rMain->getVarVal("doPlierWtaRefSeqTest") != "")
	{
	  rMain->addTest(doPlierWtaRefSeqTest());
	}
      if(rMain->getVarVal("doHumanGeneTest") != "")
	{
	  rMain->addTest(doHumanGeneTest());
	}
      if(rMain->getVarVal("doHumanGeneSpfTest") != "")
	{
	  rMain->addTest(doHumanGeneSpfTest());
	}
      if(rMain->getVarVal("doHumanGeneKillListTest") != "")
	{
	  rMain->addTest(doHumanGeneKillListTest());
	}
      if(rMain->getVarVal("doSNP6CopyNumberWorkflow1a") != "")
	{
	  rMain->addTest(doSNP6CopyNumberWorkflow1a());
	}
      if(rMain->getVarVal("doSNP6CopyNumberWorkflow1b") != "")
	{
	  rMain->addTest(doSNP6CopyNumberWorkflow1b());
	}
      if(rMain->getVarVal("doSNP6CopyNumberWorkflow2") != "")
	{
	  rMain->addTest(doSNP6CopyNumberWorkflow2());
	}
      if(rMain->getVarVal("doPlierMMTissueMedianNormTest") != "")
	{
	  rMain->addTest(doPlierMMTissueMedianNormTest());
	}
      if(rMain->getVarVal("doDabgSubsetTest") != "")
	{
	  rMain->addTest(doDabgSubsetTest());
	}
      if(rMain->getVarVal("doDabgU133Test") != "")
	{
	  rMain->addTest(doDabgU133Test());
	}
    }
  rMain->runTests();

  //DEBUG
  //rMain->dump();
}


///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doRmaTissueSketchSuppliedTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueSketchSuppliedTest");
  //Set testname variable for the purposes of pathname generation
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "quant-norm.sketch=100000.lowprecision=true,pm-only,med-polish.expon=true");
  cmd->addArg("--target-sketch", "${gold_idata_psum_hgu133plus2}/quant-norm.100000.txt");
  cmd->addArg("-p", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("-x", "6");
  
  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.021");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRmaTissueSketchSuppliedTest()
{
  // NOTE: The allowable error on this test has been raised to 0.021
  // to account for solaris sometimes inaccurate "log" function.
  // On ppc and amd64, we get the same answers, but on sparc,
  // the answers are different.  Normally it is in the low bits (8-10th place)
  // but sometimes not.  (Like in the 5th place)

  RT_Test* rt;
  rt = new RT_Test("doRmaTissueSketchSuppliedTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "quant-norm.sketch=100000.lowprecision=true,pm-only,med-polish.expon=true");
  cmd->addArg("--target-sketch", "${gold_probesetsummarize_hgu133plus2}/quant-norm.100000.txt");
  cmd->addArg("-p", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("-x", "6");
  
  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.021");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doPlierPrecompTissueTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierPrecompTissueTest");
  //rt->setIntegrationTestBool(true);
  RT_Test* rt1;
  rt1 = new RT_Test("doPlierPrecompTissueTest-part1");
//rt1->setIntegrationTestBool(true);
  RT_Test* rt2;
  rt2 = new RT_Test("doPlierPrecompTissueTest-part2");
//rt2->setIntegrationTestBool(true);

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_summarize}");
  cmd1->addArg("-a", "pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd1->addArg("-p", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd1->addArg("-c", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd1->addArg("-b", "${gold_idata_lib_hgu133plus2}/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp");
  cmd1->addArg("-o", "${out_dir}/doPlierPrecompTissueTest");
  cmd1->addArg("--feat-effects");
  cmd1->addArg("-x", "5");

  std::vector<std::string> fulltissueCels1 = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels1.begin();i < fulltissueCels1.end();i++)
    {
      std::string stringArg = *i;
      cmd1->addArg(stringArg);
    }


  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_dir}/doPlierPrecompTissueTest/pm-gcbg.plier.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum}/doPlierPrecompTissueTest/pm-gcbg.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_summarize}");
  cmd2->addArg("-a", "pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd2->addArg("-p", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd2->addArg("-c", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd2->addArg("-b", "${gold_idata_lib_hgu133plus2}/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp");
  cmd2->addArg("-o", "${out_dir}/doPlierPrecompTissueTest/useFeat");
  cmd2->addArg("--use-feat-eff", "${out_dir}/doPlierPrecompTissueTest/pm-gcbg.plier.feature-response.txt");
  cmd2->addArg("-x", "5");

  std::vector<std::string> fulltissueCels2 = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels2.begin();i < fulltissueCels2.end();i++)
    {
      std::string stringArg = *i;
      cmd2->addArg(stringArg);
    }

  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_dir}/doPlierPrecompTissueTest/useFeat/pm-gcbg.plier.summary.txt");
  check2->addArg("--gold", "${gold_idata_psum}/doPlierPrecompTissueTest/useFeat/pm-gcbg.plier.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doPlierPrecompTissueTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierPrecompTissueTest");
  RT_Test* rt1;
  rt1 = new RT_Test("doPlierPrecompTissueTest-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("doPlierPrecompTissueTest-part2");


  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_summarize}");
  cmd1->addArg("-a", "pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd1->addArg("-p", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd1->addArg("-c", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd1->addArg("-b", "${gold_lib_hgu133plus2}/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp");
  cmd1->addArg("-o", "${out_dir}/doPlierPrecompTissueTest");
  cmd1->addArg("--feat-effects");
  cmd1->addArg("-x", "5");

  std::vector<std::string> fulltissueCels1 = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels1.begin();i < fulltissueCels1.end();i++)
    {
      std::string stringArg = *i;
      cmd1->addArg(stringArg);
    }

  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_dir}/doPlierPrecompTissueTest/pm-gcbg.plier.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/doPlierPrecompTissueTest/pm-gcbg.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_summarize}");
  cmd2->addArg("-a", "pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd2->addArg("-p", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd2->addArg("-c", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd2->addArg("-b", "${gold_lib_hgu133plus2}/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp");
  cmd2->addArg("-o", "${out_dir}/doPlierPrecompTissueTest/useFeat");
  cmd2->addArg("--use-feat-eff", "${out_dir}/doPlierPrecompTissueTest/pm-gcbg.plier.feature-response.txt");
  cmd2->addArg("-x", "5");

  std::vector<std::string> fulltissueCels2 = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels2.begin();i < fulltissueCels2.end();i++)
    {
      std::string stringArg = *i;
      cmd2->addArg(stringArg);
    }

  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_dir}/doPlierPrecompTissueTest/useFeat/pm-gcbg.plier.summary.txt");
  check2->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/doPlierPrecompTissueTest/useFeat/pm-gcbg.plier.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doRmaNvissaTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaNvissaTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("--spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--xda-chp-output");
  cmd->addArg("--cc-chp-output");
   
  std::vector<std::string> fullNvissaCels = Util::addPrefixSuffix(nvissaCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixUp}");
  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  string chpFilesSuffix = ".rma-bg.quant-norm.pm-only.med-polish.chp";
  // checks for chp files.
  const char *chpFiles[] = {"4221B_a01", "4221B_b01", 
			    "4221B_c01", "4221HL_a04", 
			    "4221HL_b04", "4221HL_c04", 
			    "4221H_a02", "4221H_b02", 
			    "4221H_c02", "4221P_a03", 
			    "4221P_b03", "4221P_c03",
			    NULL
  };

  vector<string> gold,gen;
  vector<string> goldCC,genCC;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_psum_testname}/chp/", chpFilesSuffix); 
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/chp/", chpFilesSuffix);
  goldCC = Util::addPrefixSuffix(chpFiles, "${gold_idata_psum_testname}/cc-chp/", chpFilesSuffix);
  genCC = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", chpFilesSuffix);

  //CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_chp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  
  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_calvinchp}");
  check2->setGenFiles(genCC);
  check2->setGoldFiles(goldCC);

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check3->addArg("--gold", "${gold_idata_psum_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "true");
  check3->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRmaNvissaTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaNvissaTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--xda-chp-output");
  cmd->addArg("--cc-chp-output");
   
  std::vector<std::string> fullNvissaCels = Util::addPrefixSuffix(nvissaCels, "${gold_cel_hgu133plus2}/", "${celsuffixUp}");
  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  string chpFilesSuffix = ".rma-bg.quant-norm.pm-only.med-polish.chp";
  // checks for chp files.
  const char *chpFiles[] = {"4221B_a01", "4221B_b01", 
			    "4221B_c01", "4221HL_a04", 
			    "4221HL_b04", "4221HL_c04", 
			    "4221H_a02", "4221H_b02", 
			    "4221H_c02", "4221P_a03", 
			    "4221P_b03", "4221P_c03",
			    NULL
  };

  vector<string> gold,gen;
  vector<string> goldCC,genCC;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetsummarize_hgu133plus2_testname}/chp/", chpFilesSuffix); 
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/chp/", chpFilesSuffix);
  goldCC = Util::addPrefixSuffix(chpFiles, "${gold_probesetsummarize_hgu133plus2_testname}/cc-chp/", chpFilesSuffix);
  genCC = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", chpFilesSuffix);

  //CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_chp}");
  check1->setGenFiles(gen);
  check1->setGoldFiles(gold);
  
  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_calvinchp}");
  check2->setGenFiles(genCC);
  check2->setGoldFiles(goldCC);

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check3->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "true");
  check3->addArg("--allowedMismatch", "0");

  return rt;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doRmaNvissaReadWriteFeatureEffectsTest()
{
  std::vector<std::string> fullNvissaCels = Util::addPrefixSuffix(nvissaCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixUp}");  
  
  RT_Test* rt;
  rt = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest");
  //rt->setIntegrationTestBool(true);
  RT_Test* rt1;
  rt1 = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest-part1");
//rt1->setIntegrationTestBool(true);
  RT_Test* rt2;
  rt2 = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest-part2");
//rt2->setIntegrationTestBool(true);
  RT_Test* rt3;
  rt3 = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest-part3");
  //rt3->setIntegrationTestBool(true);

  // RT_Test* rt4;
  //rt4 = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest-part4");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_summarize}");
  cmd1->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd1->addArg("--spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd1->addArg("-x", "5");
  cmd1->addArg("-o", "${out_dir}/doFE");
  cmd1->addArg("--feat-effects");
  
  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd1->addArg(stringArg);
    }

  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_dir}/doFE/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum}/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt1->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_dir}/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt");
  check2->addArg("--gold", "${gold_idata_psum}/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt"); 
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_summarize}");
  cmd2->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd2->addArg("--spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd2->addArg("-x", "5");
  cmd2->addArg("-o", "${out_dir}/doFE/2");
  cmd2->addArg("--use-feat-eff", "${out_dir}/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt");

  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd2->addArg(stringArg);
    }

  ///@todo need to tighten up once we have binary/a5 input/output
  RT_Check* check3 = rt2->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_dir}/doFE/2/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check3->addArg("--gold", "${gold_idata_psum}/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check3->addArg("--epsilon", "0.001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  //  Third Test - this case checks that code which was designed to use a new style of feature effects file, continues to function
  //  when the old style of feature effects file in text format is used as input.  
  RT_Cmd* cmd3 = rt3->newCmd();
  cmd3->setExe("${apt_probeset_summarize}");
  cmd3->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd3->addArg("--spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd3->addArg("-x", "5");
  cmd3->addArg("-o", "${out_dir}/doFE/3");
  cmd3->addArg("--use-feat-eff", "${gold_idata_psum}/doFE/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt");

  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd3->addArg(stringArg);
    }

  RT_Check* check4 = rt3->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_dir}/doFE/3/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check4->addArg("--gold", "${gold_idata_psum}/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check4->addArg("--epsilon", "0.001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  /*
  //  Fourth Test - this case checks that code which was designed to use a new style of feature effects file, continues to function
  //  when the old style of feature effects file in HDF5 format is used as input.  
  RT_Cmd* cmd4 = rt4->newCmd();
  cmd4->setExe("${apt_probeset_summarize}");
  cmd4->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd4->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd4->addArg("-x", "5");
  cmd4->addArg("-o", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/useHDF5");
  cmd4->addArg("--a5-feature-effects-input-file", "${gold_probesetsummarize_hgu133plus2}/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.a5");
  cmd4->addArg("--a5-feature-effects-input-name", "rma-bg.quant-norm.pm-only.med-polish.feature-response");

  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd4->addArg(stringArg);
    }
  
  RT_Check* check5 = rt4->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/useHDF5/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check5->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check5->addArg("--epsilon", "0.001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");
  */
  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  rt->addSubtest(rt3);
  //rt->addSubtest(rt4);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRmaNvissaReadWriteFeatureEffectsTest()
{
  std::vector<std::string> fullNvissaCels = Util::addPrefixSuffix(nvissaCels, "${gold_cel_hgu133plus2}/", "${celsuffixUp}");  
  
  RT_Test* rt;
  rt = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest");
  RT_Test* rt1;
  rt1 = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest-part2");
  RT_Test* rt3;
  rt3 = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest-part3");
  RT_Test* rt4;
  rt4 = new RT_Test("doRmaNvissaReadWriteFeatureEffectsTest-part4");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_summarize}");
  cmd1->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd1->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd1->addArg("-x", "5");
  cmd1->addArg("-o", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest");
  cmd1->addArg("--feat-effects");
  
  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd1->addArg(stringArg);
    }

  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt1->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt");
  check2->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt"); 
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_summarize}");
  cmd2->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd2->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd2->addArg("-x", "5");
  cmd2->addArg("-o", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/useFeatures");
  cmd2->addArg("--use-feat-eff", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt");

  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd2->addArg(stringArg);
    }

  ///@todo need to tighten up once we have binary/a5 input/output
  RT_Check* check3 = rt2->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/useFeatures/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check3->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check3->addArg("--epsilon", "0.001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  //  Third Test - this case checks that code which was designed to use a new style of feature effects file, continues to function
  //  when the old style of feature effects file in text format is used as input.  
  RT_Cmd* cmd3 = rt3->newCmd();
  cmd3->setExe("${apt_probeset_summarize}");
  cmd3->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd3->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd3->addArg("-x", "5");
  cmd3->addArg("-o", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/useText");
  cmd3->addArg("--use-feat-eff", "${gold_probesetsummarize_hgu133plus2}/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.txt");

  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd3->addArg(stringArg);
    }

  RT_Check* check4 = rt3->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/useText/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check4->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check4->addArg("--epsilon", "0.001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  //  Fourth Test - this case checks that code which was designed to use a new style of feature effects file, continues to function
  //  when the old style of feature effects file in HDF5 format is used as input.  
  RT_Cmd* cmd4 = rt4->newCmd();
  cmd4->setExe("${apt_probeset_summarize}");
  cmd4->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd4->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd4->addArg("-x", "5");
  cmd4->addArg("-o", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/useHDF5");
  cmd4->addArg("--a5-feature-effects-input-file", "${gold_probesetsummarize_hgu133plus2}/doRmaNvissaReadWriteFeatureEffectsTest/rma-bg.quant-norm.pm-only.med-polish.feature-response.a5");
  cmd4->addArg("--a5-feature-effects-input-name", "rma-bg.quant-norm.pm-only.med-polish.feature-response");

  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd4->addArg(stringArg);
    }
  
  RT_Check* check5 = rt4->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_dir}/doRmaNvissaReadWriteFeatureEffectsTest/useHDF5/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check5->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/doRmaNvissaTest/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check5->addArg("--epsilon", "0.001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "false");
  check5->addArg("--allowedMismatch", "0");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  rt->addSubtest(rt3);
  rt->addSubtest(rt4);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doRmaNvissaSketchTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaNvissaSketchTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=100000.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("--spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullNvissaCels = Util::addPrefixSuffix(nvissaCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixUp}");
  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  std::vector<std::string> checkArgs1;
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "225");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRmaNvissaSketchTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaNvissaSketchTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=100000.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullNvissaCels = Util::addPrefixSuffix(nvissaCels, "${gold_cel_hgu133plus2}/", "${celsuffixUp}");
  for(std::vector<std::string>::iterator i = fullNvissaCels.begin();i < fullNvissaCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  std::vector<std::string> checkArgs1;
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "225");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doRmaTissueTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("-spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${qt_celFileList}");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRmaTissueTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueTest");


  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${celFileList}");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doRmaTissueSpfTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueSpfTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("--spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${qt_celFileList}");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRmaTissueSpfTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueSpfTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("--spf-file", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${celFileList}");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doRmaTissueTestCC()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueTestCC");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("--spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}"); 

  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${agcccelsuffix}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");
  
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRmaTissueTestCC()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueTestCC");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=0.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}"); 

  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${agcccelsuffix}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2}/doRmaTissueTestCC/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");
  
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doRmaTissueSketchTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueSketchTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=100000.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("--spf-file", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.spf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.01");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "1306");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doRmaTissueSketchTest()
{
  RT_Test* rt;
  rt = new RT_Test("doRmaTissueSketchTest");
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-bg,quant-norm.sketch=100000.bioc=true.usepm=true,pm-only,med-polish");
  cmd->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.01");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "1306");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doPlierGcBgTissueTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierGcBgTissueTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd->addArg("-b", "${gold_idata_lib_hgu133plus2}/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp");
  cmd->addArg("-p", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  
  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/pm-gcbg.plier.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/pm-gcbg.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doPlierGcBgTissueTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierGcBgTissueTest");
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd->addArg("-b", "${gold_lib_hgu133plus2}/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp");
  cmd->addArg("-p", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  
  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/pm-gcbg.plier.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/pm-gcbg.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doPlierMMTissueSketchTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierMMTissueSketchTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "quant-norm.sketch=100000.lowprecision=true,pm-mm,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd->addArg("-p", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-x", "6");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-mm.plier.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-mm.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doPlierMMTissueSketchTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierMMTissueSketchTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "quant-norm.sketch=100000.lowprecision=true,pm-mm,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd->addArg("-p", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-x", "6");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fulltissueCels = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fulltissueCels.begin();i < fulltissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-mm.plier.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/quant-norm.pm-mm.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doPlierWtaRefSeqTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierWtaRefSeqTest");
  //rt->setIntegrationTestBool(true);
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd->addArg("-b", "${gold_idata_lib_huex10stv2}/antigenomic.bgp");
  cmd->addArg("-p", "${gold_idata_lib_huex10stv2}/HuEx-1_0-st-v2.pgf");
  cmd->addArg("-c", "${gold_idata_lib_huex10stv2}/HuEx-1_0-st-v2.clf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_doRS}");
  cmd->addArg("-m", "${gold_idata_lib_huex10stv2}/map.refseq.txt");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("-a", "gc-bg,pm-only,med-polish");
  cmd->addArg("-a", "pm-gcbg,med-polish");
  cmd->addArg("-a", "rma-bg,quant-norm.bioc=true,pm-only,med-polish");
  cmd->addArg("-a", "rma-bg,quant-norm.bioc=true,pm-only,med-polish,pca-select.hard-min=5.min-percent=0.qnorm-only=false.info-criterion=bic");
  cmd->addArg("-a", "rma-bg,quant-norm.bioc=true,pm-only,med-polish,spect-select.metric=angle.log2=false.cut-val=zero.margin=.9.min-percent=0.0.info-criterion=bic");
  
  std::vector<std::string> fullhuexCels = Util::addPrefixSuffix(huexCels, "${gold_idata_cel_huex10stv2}/", "${celsuffixUp}");
  for(std::vector<std::string>::iterator i = fullhuexCels.begin();i < fullhuexCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  string chpFilesSuffix = ".chp";
  const char *chpFiles[] = {"huex_wta_cb_A.gc-bg.pm-only.med-polish", "huex_wta_cb_A.pm-gcbg.med-polish", "huex_wta_cb_A.pm-gcbg.plier", 
			    "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			    "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_cb_B.gc-bg.pm-only.med-polish", 
			    "huex_wta_cb_B.pm-gcbg.med-polish", "huex_wta_cb_B.pm-gcbg.plier", "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish", 
			    "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			    "huex_wta_cb_C.gc-bg.pm-only.med-polish", "huex_wta_cb_C.pm-gcbg.med-polish", "huex_wta_cb_C.pm-gcbg.plier", 
			    "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			    "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_heart_A.gc-bg.pm-only.med-polish", 
			    "huex_wta_heart_A.pm-gcbg.med-polish", "huex_wta_heart_A.pm-gcbg.plier", "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish", 
			    "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			    "huex_wta_heart_B.gc-bg.pm-only.med-polish", "huex_wta_heart_B.pm-gcbg.med-polish", "huex_wta_heart_B.pm-gcbg.plier", 
			    "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			    "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_heart_C.gc-bg.pm-only.med-polish", 
			    "huex_wta_heart_C.pm-gcbg.med-polish", "huex_wta_heart_C.pm-gcbg.plier", "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish", 
			    "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			    "huex_wta_muscle_A.gc-bg.pm-only.med-polish", "huex_wta_muscle_A.pm-gcbg.med-polish", "huex_wta_muscle_A.pm-gcbg.plier", 
			    "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			    "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_muscle_B.gc-bg.pm-only.med-polish", 
			    "huex_wta_muscle_B.pm-gcbg.med-polish", "huex_wta_muscle_B.pm-gcbg.plier", "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish", 
			    "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			    "huex_wta_muscle_C.gc-bg.pm-only.med-polish", "huex_wta_muscle_C.pm-gcbg.med-polish", "huex_wta_muscle_C.pm-gcbg.plier", 
			    "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			    "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_testes_A.gc-bg.pm-only.med-polish", 
			    "huex_wta_testes_A.pm-gcbg.med-polish", "huex_wta_testes_A.pm-gcbg.plier", "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish", 
			    "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			    "huex_wta_testes_B.gc-bg.pm-only.med-polish", "huex_wta_testes_B.pm-gcbg.med-polish", "huex_wta_testes_B.pm-gcbg.plier", 
			    "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			    "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_testes_C.gc-bg.pm-only.med-polish", 
			    "huex_wta_testes_C.pm-gcbg.med-polish", "huex_wta_testes_C.pm-gcbg.plier", "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish", 
			    "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish.spect-select",
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_psum}/doRS/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles,"${out_doRS}/cc-chp/" , chpFilesSuffix);

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_doRS}/pm-gcbg.plier.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum}/doRS/pm-gcbg.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_doRS}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
  check2->addArg("--gold", "${gold_idata_psum}/doRS/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_doRS}/rma-bg.quant-norm.pm-only.med-polish.pca-select.summary.txt");
  check3->addArg("--gold", "${gold_idata_psum}/doRS/rma-bg.quant-norm.pm-only.med-polish.pca-select.summary.txt"); 
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "true");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_doRS}/rma-bg.quant-norm.pm-only.med-polish.spect-select.summary.txt");
  check4->addArg("--gold", "${gold_idata_psum}/doRS/rma-bg.quant-norm.pm-only.med-polish.spect-select.summary.txt"); 
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "true");
  check4->addArg("--allowedMismatch", "0");

  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_doRS}/pm-gcbg.med-polish.summary.txt");
  check5->addArg("--gold", "${out_doRS}/gc-bg.pm-only.med-polish.summary.txt"); 
  check5->addArg("--epsilon", "0.0001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "true");
  check5->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doPlierWtaRefSeqTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierWtaRefSeqTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "pm-gcbg,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd->addArg("-b", "${gold_lib_huex10stv2}/antigenomic.probe_id-only.bgp");
  cmd->addArg("-p", "${gold_lib_huex10stv2}/HuEx-1_0-st-v2.pgf");
  cmd->addArg("-c", "${gold_lib_huex10stv2}/HuEx-1_0-st-v2.clf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("-m", "${gold_lib_huex10stv2}/map.refseq.txt");
  cmd->addArg("--cc-chp-output");
  cmd->addArg("-a", "gc-bg,pm-only,med-polish");
  cmd->addArg("-a", "pm-gcbg,med-polish");
  cmd->addArg("-a", "rma-bg,quant-norm.bioc=true,pm-only,med-polish");
  cmd->addArg("-a", "rma-bg,quant-norm.bioc=true,pm-only,med-polish,pca-select.hard-min=5.min-percent=0.qnorm-only=false.info-criterion=bic");
  cmd->addArg("-a", "rma-bg,quant-norm.bioc=true,pm-only,med-polish,spect-select.metric=angle.log2=false.cut-val=zero.margin=.9.min-percent=0.0.info-criterion=bic");
  
  std::vector<std::string> fullhuexCels = Util::addPrefixSuffix(huexCels, "${gold_cel_huex10stv2}/", "${celsuffixUp}");
  for(std::vector<std::string>::iterator i = fullhuexCels.begin();i < fullhuexCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

 string chpFilesSuffix = ".chp";
 const char *chpFiles[] = {"huex_wta_cb_A.gc-bg.pm-only.med-polish", "huex_wta_cb_A.pm-gcbg.med-polish", "huex_wta_cb_A.pm-gcbg.plier", 
			   "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			   "huex_wta_cb_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_cb_B.gc-bg.pm-only.med-polish", 
			   "huex_wta_cb_B.pm-gcbg.med-polish", "huex_wta_cb_B.pm-gcbg.plier", "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish", 
			   "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_cb_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			   "huex_wta_cb_C.gc-bg.pm-only.med-polish", "huex_wta_cb_C.pm-gcbg.med-polish", "huex_wta_cb_C.pm-gcbg.plier", 
			   "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			   "huex_wta_cb_C.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_heart_A.gc-bg.pm-only.med-polish", 
			   "huex_wta_heart_A.pm-gcbg.med-polish", "huex_wta_heart_A.pm-gcbg.plier", "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish", 
			   "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_heart_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			   "huex_wta_heart_B.gc-bg.pm-only.med-polish", "huex_wta_heart_B.pm-gcbg.med-polish", "huex_wta_heart_B.pm-gcbg.plier", 
			   "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			   "huex_wta_heart_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_heart_C.gc-bg.pm-only.med-polish", 
			   "huex_wta_heart_C.pm-gcbg.med-polish", "huex_wta_heart_C.pm-gcbg.plier", "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish", 
			   "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_heart_C.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			   "huex_wta_muscle_A.gc-bg.pm-only.med-polish", "huex_wta_muscle_A.pm-gcbg.med-polish", "huex_wta_muscle_A.pm-gcbg.plier", 
			   "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			   "huex_wta_muscle_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_muscle_B.gc-bg.pm-only.med-polish", 
			   "huex_wta_muscle_B.pm-gcbg.med-polish", "huex_wta_muscle_B.pm-gcbg.plier", "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish", 
			   "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_muscle_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			   "huex_wta_muscle_C.gc-bg.pm-only.med-polish", "huex_wta_muscle_C.pm-gcbg.med-polish", "huex_wta_muscle_C.pm-gcbg.plier", 
			   "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			   "huex_wta_muscle_C.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_testes_A.gc-bg.pm-only.med-polish", 
			   "huex_wta_testes_A.pm-gcbg.med-polish", "huex_wta_testes_A.pm-gcbg.plier", "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish", 
			   "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_testes_A.rma-bg.quant-norm.pm-only.med-polish.spect-select", 
			   "huex_wta_testes_B.gc-bg.pm-only.med-polish", "huex_wta_testes_B.pm-gcbg.med-polish", "huex_wta_testes_B.pm-gcbg.plier", 
			   "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish", "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish.pca-select", 
			   "huex_wta_testes_B.rma-bg.quant-norm.pm-only.med-polish.spect-select", "huex_wta_testes_C.gc-bg.pm-only.med-polish", 
			   "huex_wta_testes_C.pm-gcbg.med-polish", "huex_wta_testes_C.pm-gcbg.plier", "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish", 
			   "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish.pca-select", "huex_wta_testes_C.rma-bg.quant-norm.pm-only.med-polish.spect-select",
			   NULL
 };

 vector<string> gold,gen;


 gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetsummarize_huex10stv2_testname}/cc-chp/", chpFilesSuffix);
 gen = Util::addPrefixSuffix(chpFiles,"${out_testname}/cc-chp/" , chpFilesSuffix);

 //CALVIN CHP CHECK
 //Generate Check objects for CalvinChpCheck and store them in rt
 RT_Check* check = rt->newCheck();
 check->setExe("${apt_check_calvinchp}");
 check->setGenFiles(gen);
 check->setGoldFiles(gold);

 RT_Check* check1 = rt->newCheck();
 check1->setExe("${apt_check_matrix}");
 check1->addArg("--gen", "${out_testname}/pm-gcbg.plier.summary.txt");
 check1->addArg("--gold", "${gold_probesetsummarize_huex10stv2_testname}/pm-gcbg.plier.summary.txt"); 
 check1->addArg("--epsilon", "0.0001");
 check1->addArg("--rowSkip", "1"); 
 check1->addArg("--columnSkip", "1");
 check1->addArg("--matchNames", "true");
 check1->addArg("--allowedMismatch", "0");

 RT_Check* check2 = rt->newCheck();
 check2->setExe("${apt_check_matrix}");
 check2->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt");
 check2->addArg("--gold", "${gold_probesetsummarize_huex10stv2_testname}/rma-bg.quant-norm.pm-only.med-polish.summary.txt"); 
 check2->addArg("--epsilon", "0.0001");
 check2->addArg("--rowSkip", "1"); 
 check2->addArg("--columnSkip", "1");
 check2->addArg("--matchNames", "true");
 check2->addArg("--allowedMismatch", "0");

 RT_Check* check3 = rt->newCheck();
 check3->setExe("${apt_check_matrix}");
 check3->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.pca-select.summary.txt");
 check3->addArg("--gold", "${gold_probesetsummarize_huex10stv2_testname}/rma-bg.quant-norm.pm-only.med-polish.pca-select.summary.txt"); 
 check3->addArg("--epsilon", "0.0001");
 check3->addArg("--rowSkip", "1"); 
 check3->addArg("--columnSkip", "1");
 check3->addArg("--matchNames", "true");
 check3->addArg("--allowedMismatch", "0");

 RT_Check* check4 = rt->newCheck();
 check4->setExe("${apt_check_matrix}");
 check4->addArg("--gen", "${out_testname}/rma-bg.quant-norm.pm-only.med-polish.spect-select.summary.txt");
 check4->addArg("--gold", "${gold_probesetsummarize_huex10stv2_testname}/rma-bg.quant-norm.pm-only.med-polish.spect-select.summary.txt"); 
 check4->addArg("--epsilon", "0.0001");
 check4->addArg("--rowSkip", "1"); 
 check4->addArg("--columnSkip", "1");
 check4->addArg("--matchNames", "true");
 check4->addArg("--allowedMismatch", "0");

 RT_Check* check5 = rt->newCheck();
 check5->setExe("${apt_check_matrix}");
 check5->addArg("--gen", "${out_testname}/pm-gcbg.med-polish.summary.txt");
 check5->addArg("--gold", "${out_testname}/gc-bg.pm-only.med-polish.summary.txt"); 
 check5->addArg("--epsilon", "0.0001");
 check5->addArg("--rowSkip", "1"); 
 check5->addArg("--columnSkip", "1");
 check5->addArg("--matchNames", "true");
 check5->addArg("--allowedMismatch", "0");

 return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doHumanGeneTest()
{
  RT_Test* rt;
  //rt = new RT_Test("doHumanGeneTest");
  rt = new RT_Test("doHGT");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "plier-gcbg-sketch");
  cmd->addArg("-a", "rma-sketch");
  cmd->addArg("-a", "quant-norm,pm-only,med-polish,pca-select");
  cmd->addArg("-a", "quant-norm,pm-only,med-polish,spect-select");
  cmd->addArg("-b", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.bgp");
  cmd->addArg("-c", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd->addArg("-p", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd->addArg("--qc-probesets", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.qcc");
  cmd->addArg("-m", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.mps");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${qt_celFileList2}");
  cmd->addArg("--feat-effects");
  cmd->addArg("--feat-details");
  cmd->addArg("--write-sketch");
  cmd->addArg("--cc-chp-output");

  std::string chpFilesSuffix = ".chp";
  const char *chpFiles[] = {"TisMap_Brain_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Brain_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Brain_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Brain_01_v1_WTGene1.rma-sketch",
			    "TisMap_Breast_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Breast_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Breast_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Breast_01_v1_WTGene1.rma-sketch",
			    "TisMap_Heart_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Heart_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Heart_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Heart_01_v1_WTGene1.rma-sketch",
			    "TisMap_Kidney_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Kidney_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Kidney_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Kidney_01_v1_WTGene1.rma-sketch",
			    "TisMap_Liver_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Liver_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Liver_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Liver_01_v1_WTGene1.rma-sketch",
			    "TisMap_Panc_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Panc_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Panc_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Panc_01_v1_WTGene1.rma-sketch",
			    "TisMap_Prost_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Prost_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Prost_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Prost_01_v1_WTGene1.rma-sketch",
			    "TisMap_SkMus_01_v1_WTGene1.plier-gcbg-sketch","TisMap_SkMus_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_SkMus_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_SkMus_01_v1_WTGene1.rma-sketch",
			    "TisMap_Spleen_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Spleen_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Spleen_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Spleen_01_v1_WTGene1.rma-sketch",
			    "TisMap_Testis_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Testis_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Testis_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Testis_01_v1_WTGene1.rma-sketch",
			    "TisMap_Thyroid_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Thyroid_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Thyroid_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Thyroid_01_v1_WTGene1.rma-sketch",
			    NULL
  };

  std::vector<std::string> gold, gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_psum_testname}/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", chpFilesSuffix);

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/plier-gcbg-sketch.feature-response.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/plier-gcbg-sketch.feature-response.txt"); 
  check1->addArg("--epsilon", "0.001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "0");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.feature-response.txt");
  check2->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.feature-response.txt"); 
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "0");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.feature-response.txt");
  check3->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.feature-response.txt"); 
  check3->addArg("--epsilon", "0.001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "0");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/rma-sketch.feature-response.txt");
  check4->addArg("--gold", "${gold_idata_psum_testname}/rma-sketch.feature-response.txt"); 
  check4->addArg("--epsilon", "0.001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "0");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_testname}/plier-gcbg-sketch.report.txt");
  check5->addArg("--gold", "${gold_idata_psum_testname}/plier-gcbg-sketch.report.txt"); 
  check5->addArg("--epsilon", "0.0001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "true");
  check5->addArg("--allowedMismatch", "0");

  RT_Check* check6 = rt->newCheck();
  check6->setExe("${apt_check_matrix}");
  check6->addArg("--gen", "${out_testname}/rma-sketch.report.txt");
  check6->addArg("--gold", "${gold_idata_psum_testname}/rma-sketch.report.txt"); 
  check6->addArg("--epsilon", "0.0001");
  check6->addArg("--rowSkip", "1"); 
  check6->addArg("--columnSkip", "1");
  check6->addArg("--matchNames", "true");
  check6->addArg("--allowedMismatch", "0");

  RT_Check* check7 = rt->newCheck();
  check7->setExe("${apt_check_matrix}");
  check7->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.report.txt");
  check7->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.report.txt"); 
  check7->addArg("--epsilon", "0.0001");
  check7->addArg("--rowSkip", "1"); 
  check7->addArg("--columnSkip", "1");
  check7->addArg("--matchNames", "true");
  check7->addArg("--allowedMismatch", "0");

  RT_Check* check8 = rt->newCheck();
  check8->setExe("${apt_check_matrix}");
  check8->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.spect-select.report.txt");
  check8->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.spect-select.report.txt"); 
  check8->addArg("--epsilon", "0.0001");
  check8->addArg("--rowSkip", "1"); 
  check8->addArg("--columnSkip", "1");
  check8->addArg("--matchNames", "true");
  check8->addArg("--allowedMismatch", "0");

  //MIXED FILE CHECK
  RT_Check* check9 = rt->newCheck();
  check9->setExe("${apt_check_mixedfile}");
  check9->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt");
  check9->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt");
  check9->addArg("--epsilon", "0.0001");
  check9->addArg("--skipLines", "0");
  check9->addArg("--allowedMismatch", "0");

  //MIXED FILE CHECK
  RT_Check* check10 = rt->newCheck();
  check10->setExe("${apt_check_mixedfile}");
  check10->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.spect-select.spect-select.report.txt");
  check10->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.spect-select.spect-select.report.txt");
  check10->addArg("--epsilon", "0.0001");
  check10->addArg("--skipLines", "0");
  check10->addArg("--allowedMismatch", "0");

  std::vector<std::string> checkArgs11;
  RT_Check* check11 = rt->newCheck();
  check11->setExe("${apt_check_matrix}");
  check11->addArg("--gen", "${out_testname}/plier-gcbg-sketch.residuals.txt");
  check11->addArg("--gold", "${gold_idata_psum_testname}/plier-gcbg-sketch.residuals.txt"); 
  check11->addArg("--epsilon", "0.001");
  check11->addArg("--rowSkip", "1"); 
  check11->addArg("--columnSkip", "0");
  check11->addArg("--matchNames", "false");
  check11->addArg("--allowedMismatch", "0");

  RT_Check* check12 = rt->newCheck();
  check12->setExe("${apt_check_matrix}");
  check12->addArg("--gen", "${out_testname}/rma-sketch.residuals.txt");
  check12->addArg("--gold", "${gold_idata_psum_testname}/rma-sketch.residuals.txt"); 
  check12->addArg("--epsilon", "0.001");
  check12->addArg("--rowSkip", "1"); 
  check12->addArg("--columnSkip", "0");
  check12->addArg("--matchNames", "false");
  check12->addArg("--allowedMismatch", "0");

  RT_Check* check13 = rt->newCheck();
  check13->setExe("${apt_check_matrix}");
  check13->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.residuals.txt");
  check13->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.residuals.txt"); 
  check13->addArg("--epsilon", "0.001");
  check13->addArg("--rowSkip", "1"); 
  check13->addArg("--columnSkip", "0");
  check13->addArg("--matchNames", "false");
  check13->addArg("--allowedMismatch", "0");

  RT_Check* check14 = rt->newCheck();
  check14->setExe("${apt_check_matrix}");
  check14->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.spect-select.residuals.txt");
  check14->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.spect-select.residuals.txt"); 
  check14->addArg("--epsilon", "0.001");
  check14->addArg("--rowSkip", "1"); 
  check14->addArg("--columnSkip", "0");
  check14->addArg("--matchNames", "false");
  check14->addArg("--allowedMismatch", "0");

  RT_Check* check15 = rt->newCheck();
  check15->setExe("${apt_check_matrix}");
  check15->addArg("--gen", "${out_testname}/plier-gcbg-sketch.summary.txt");
  check15->addArg("--gold", "${gold_idata_psum_testname}/plier-gcbg-sketch.summary.txt"); 
  check15->addArg("--epsilon", "0.0001");
  check15->addArg("--rowSkip", "1"); 
  check15->addArg("--columnSkip", "0");
  check15->addArg("--matchNames", "false");
  check15->addArg("--allowedMismatch", "0");

  RT_Check* check16 = rt->newCheck();
  check16->setExe("${apt_check_matrix}");
  check16->addArg("--gen", "${out_testname}/rma-sketch.summary.txt");
  check16->addArg("--gold", "${gold_idata_psum_testname}/rma-sketch.summary.txt"); 
  check16->addArg("--epsilon", "0.0001");
  check16->addArg("--rowSkip", "1"); 
  check16->addArg("--columnSkip", "0");
  check16->addArg("--matchNames", "false");
  check16->addArg("--allowedMismatch", "0");

  RT_Check* check17 = rt->newCheck();
  check17->setExe("${apt_check_matrix}");
  check17->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.summary.txt");
  check17->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.summary.txt"); 
  check17->addArg("--epsilon", "0.0001");
  check17->addArg("--rowSkip", "1"); 
  check17->addArg("--columnSkip", "0");
  check17->addArg("--matchNames", "false");
  check17->addArg("--allowedMismatch", "0");

  RT_Check* check18 = rt->newCheck();
  check18->setExe("${apt_check_matrix}");
  check18->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.spect-select.summary.txt");
  check18->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.spect-select.summary.txt"); 
  check18->addArg("--epsilon", "0.0001");
  check18->addArg("--rowSkip", "1"); 
  check18->addArg("--columnSkip", "0");
  check18->addArg("--matchNames", "false");
  check18->addArg("--allowedMismatch", "0");

  RT_Check* check19 = rt->newCheck();
  check19->setExe("${apt_check_matrix}");
  check19->addArg("--gen", "${out_testname}/quant-norm.normalization-target.txt");
  check19->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.normalization-target.txt"); 
  check19->addArg("--epsilon", "0.0001");
  check19->addArg("--rowSkip", "1"); 
  check19->addArg("--columnSkip", "0");
  check19->addArg("--matchNames", "false");
  check19->addArg("--allowedMismatch", "0");

  RT_Check* check20 = rt->newCheck();
  check20->setExe("${apt_check_matrix}");
  check20->addArg("--gen", "${out_testname}/rma-bg.quant-norm.normalization-target.txt");
  check20->addArg("--gold", "${gold_idata_psum_testname}/rma-bg.quant-norm.normalization-target.txt"); 
  check20->addArg("--epsilon", "0.0001");
  check20->addArg("--rowSkip", "1"); 
  check20->addArg("--columnSkip", "0");
  check20->addArg("--matchNames", "false");
  check20->addArg("--allowedMismatch", "0");

  // windows chokes on file names > MAX_PATH (typically 260 characters)
  RT_Check* check21 = rt->newCheck();
  check21->setExe("${apt_check_matrix}");
  check21->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.self-qnormquant-norm.normalization-target.txt");
  check21->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.self-qnormquant-norm.normalization-target.txt"); 
  check21->addArg("--epsilon", "0.0001");
  check21->addArg("--rowSkip", "1"); 
  check21->addArg("--columnSkip", "0");
  check21->addArg("--matchNames", "false");
  check21->addArg("--allowedMismatch", "0");

  return rt;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doHumanGeneTest()
{
  RT_Test* rt;
  //rt = new RT_Test("doHumanGeneTest");
  rt = new RT_Test("doHGT");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "plier-gcbg-sketch");
  cmd->addArg("-a", "rma-sketch");
  cmd->addArg("-a", "quant-norm,pm-only,med-polish,pca-select");
  cmd->addArg("-a", "quant-norm,pm-only,med-polish,spect-select");
  cmd->addArg("-b", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.bgp");
  cmd->addArg("-c", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd->addArg("-p", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd->addArg("--qc-probesets", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.qcc");
  cmd->addArg("-m", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.mps");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${celFileList2}");
  cmd->addArg("--feat-effects");
  cmd->addArg("--feat-details");
  cmd->addArg("--write-sketch");
  cmd->addArg("--cc-chp-output");

  std::string chpFilesSuffix = ".chp";
  const char *chpFiles[] = {"TisMap_Brain_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Brain_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Brain_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Brain_01_v1_WTGene1.rma-sketch",
			    "TisMap_Breast_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Breast_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Breast_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Breast_01_v1_WTGene1.rma-sketch",
			    "TisMap_Heart_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Heart_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Heart_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Heart_01_v1_WTGene1.rma-sketch",
			    "TisMap_Kidney_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Kidney_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Kidney_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Kidney_01_v1_WTGene1.rma-sketch",
			    "TisMap_Liver_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Liver_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Liver_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Liver_01_v1_WTGene1.rma-sketch",
			    "TisMap_Panc_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Panc_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Panc_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Panc_01_v1_WTGene1.rma-sketch",
			    "TisMap_Prost_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Prost_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Prost_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Prost_01_v1_WTGene1.rma-sketch",
			    "TisMap_SkMus_01_v1_WTGene1.plier-gcbg-sketch","TisMap_SkMus_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_SkMus_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_SkMus_01_v1_WTGene1.rma-sketch",
			    "TisMap_Spleen_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Spleen_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Spleen_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Spleen_01_v1_WTGene1.rma-sketch",
			    "TisMap_Testis_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Testis_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Testis_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Testis_01_v1_WTGene1.rma-sketch",
			    "TisMap_Thyroid_01_v1_WTGene1.plier-gcbg-sketch","TisMap_Thyroid_01_v1_WTGene1.quant-norm.pm-only.med-polish.pca-select",
			    "TisMap_Thyroid_01_v1_WTGene1.quant-norm.pm-only.med-polish.spect-select","TisMap_Thyroid_01_v1_WTGene1.rma-sketch",
			    NULL
  };

  std::vector<std::string> gold, gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetsummarize_hugene10stv1_testname}/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", chpFilesSuffix);

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/plier-gcbg-sketch.feature-response.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/plier-gcbg-sketch.feature-response.txt"); 
  check1->addArg("--epsilon", "0.001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "0");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.feature-response.txt");
  check2->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.pca-select.feature-response.txt"); 
  check2->addArg("--epsilon", "0.001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "0");
  check2->addArg("--matchNames", "false");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.feature-response.txt");
  check3->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.pca-select.feature-response.txt"); 
  check3->addArg("--epsilon", "0.001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "0");
  check3->addArg("--matchNames", "false");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "${out_testname}/rma-sketch.feature-response.txt");
  check4->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/rma-sketch.feature-response.txt"); 
  check4->addArg("--epsilon", "0.001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "0");
  check4->addArg("--matchNames", "false");
  check4->addArg("--allowedMismatch", "0");

  RT_Check* check5 = rt->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "${out_testname}/plier-gcbg-sketch.report.txt");
  check5->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/plier-gcbg-sketch.report.txt"); 
  check5->addArg("--epsilon", "0.0001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "true");
  check5->addArg("--allowedMismatch", "0");

  RT_Check* check6 = rt->newCheck();
  check6->setExe("${apt_check_matrix}");
  check6->addArg("--gen", "${out_testname}/rma-sketch.report.txt");
  check6->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/rma-sketch.report.txt"); 
  check6->addArg("--epsilon", "0.0001");
  check6->addArg("--rowSkip", "1"); 
  check6->addArg("--columnSkip", "1");
  check6->addArg("--matchNames", "true");
  check6->addArg("--allowedMismatch", "0");

  RT_Check* check7 = rt->newCheck();
  check7->setExe("${apt_check_matrix}");
  check7->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.report.txt");
  check7->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.pca-select.report.txt"); 
  check7->addArg("--epsilon", "0.0001");
  check7->addArg("--rowSkip", "1"); 
  check7->addArg("--columnSkip", "1");
  check7->addArg("--matchNames", "true");
  check7->addArg("--allowedMismatch", "0");

  RT_Check* check8 = rt->newCheck();
  check8->setExe("${apt_check_matrix}");
  check8->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.spect-select.report.txt");
  check8->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.spect-select.report.txt"); 
  check8->addArg("--epsilon", "0.0001");
  check8->addArg("--rowSkip", "1"); 
  check8->addArg("--columnSkip", "1");
  check8->addArg("--matchNames", "true");
  check8->addArg("--allowedMismatch", "0");

  //MIXED FILE CHECK
  RT_Check* check9 = rt->newCheck();
  check9->setExe("${apt_check_mixedfile}");
  check9->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt");
  check9->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt");
  check9->addArg("--epsilon", "0.0001");
  check9->addArg("--skipLines", "0");
  check9->addArg("--allowedMismatch", "0");

  //MIXED FILE CHECK
  RT_Check* check10 = rt->newCheck();
  check10->setExe("${apt_check_mixedfile}");
  check10->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.spect-select.spect-select.report.txt");
  check10->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.spect-select.spect-select.report.txt");
  check10->addArg("--epsilon", "0.0001");
  check10->addArg("--skipLines", "0");
  check10->addArg("--allowedMismatch", "0");

  std::vector<std::string> checkArgs11;
  RT_Check* check11 = rt->newCheck();
  check11->setExe("${apt_check_matrix}");
  check11->addArg("--gen", "${out_testname}/plier-gcbg-sketch.residuals.txt");
  check11->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/plier-gcbg-sketch.residuals.txt"); 
  check11->addArg("--epsilon", "0.001");
  check11->addArg("--rowSkip", "1"); 
  check11->addArg("--columnSkip", "0");
  check11->addArg("--matchNames", "false");
  check11->addArg("--allowedMismatch", "0");

  RT_Check* check12 = rt->newCheck();
  check12->setExe("${apt_check_matrix}");
  check12->addArg("--gen", "${out_testname}/rma-sketch.residuals.txt");
  check12->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/rma-sketch.residuals.txt"); 
  check12->addArg("--epsilon", "0.001");
  check12->addArg("--rowSkip", "1"); 
  check12->addArg("--columnSkip", "0");
  check12->addArg("--matchNames", "false");
  check12->addArg("--allowedMismatch", "0");

  RT_Check* check13 = rt->newCheck();
  check13->setExe("${apt_check_matrix}");
  check13->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.residuals.txt");
  check13->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.pca-select.residuals.txt"); 
  check13->addArg("--epsilon", "0.001");
  check13->addArg("--rowSkip", "1"); 
  check13->addArg("--columnSkip", "0");
  check13->addArg("--matchNames", "false");
  check13->addArg("--allowedMismatch", "0");

  RT_Check* check14 = rt->newCheck();
  check14->setExe("${apt_check_matrix}");
  check14->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.spect-select.residuals.txt");
  check14->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.spect-select.residuals.txt"); 
  check14->addArg("--epsilon", "0.001");
  check14->addArg("--rowSkip", "1"); 
  check14->addArg("--columnSkip", "0");
  check14->addArg("--matchNames", "false");
  check14->addArg("--allowedMismatch", "0");

  RT_Check* check15 = rt->newCheck();
  check15->setExe("${apt_check_matrix}");
  check15->addArg("--gen", "${out_testname}/plier-gcbg-sketch.summary.txt");
  check15->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/plier-gcbg-sketch.summary.txt"); 
  check15->addArg("--epsilon", "0.0001");
  check15->addArg("--rowSkip", "1"); 
  check15->addArg("--columnSkip", "0");
  check15->addArg("--matchNames", "false");
  check15->addArg("--allowedMismatch", "0");

  RT_Check* check16 = rt->newCheck();
  check16->setExe("${apt_check_matrix}");
  check16->addArg("--gen", "${out_testname}/rma-sketch.summary.txt");
  check16->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/rma-sketch.summary.txt"); 
  check16->addArg("--epsilon", "0.0001");
  check16->addArg("--rowSkip", "1"); 
  check16->addArg("--columnSkip", "0");
  check16->addArg("--matchNames", "false");
  check16->addArg("--allowedMismatch", "0");

  RT_Check* check17 = rt->newCheck();
  check17->setExe("${apt_check_matrix}");
  check17->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.summary.txt");
  check17->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.pca-select.summary.txt"); 
  check17->addArg("--epsilon", "0.0001");
  check17->addArg("--rowSkip", "1"); 
  check17->addArg("--columnSkip", "0");
  check17->addArg("--matchNames", "false");
  check17->addArg("--allowedMismatch", "0");

  RT_Check* check18 = rt->newCheck();
  check18->setExe("${apt_check_matrix}");
  check18->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.spect-select.summary.txt");
  check18->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.pm-only.med-polish.spect-select.summary.txt"); 
  check18->addArg("--epsilon", "0.0001");
  check18->addArg("--rowSkip", "1"); 
  check18->addArg("--columnSkip", "0");
  check18->addArg("--matchNames", "false");
  check18->addArg("--allowedMismatch", "0");

  RT_Check* check19 = rt->newCheck();
  check19->setExe("${apt_check_matrix}");
  check19->addArg("--gen", "${out_testname}/quant-norm.normalization-target.txt");
  check19->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/quant-norm.normalization-target.txt"); 
  check19->addArg("--epsilon", "0.0001");
  check19->addArg("--rowSkip", "1"); 
  check19->addArg("--columnSkip", "0");
  check19->addArg("--matchNames", "false");
  check19->addArg("--allowedMismatch", "0");

  RT_Check* check20 = rt->newCheck();
  check20->setExe("${apt_check_matrix}");
  check20->addArg("--gen", "${out_testname}/rma-bg.quant-norm.normalization-target.txt");
  check20->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/rma-bg.quant-norm.normalization-target.txt"); 
  check20->addArg("--epsilon", "0.0001");
  check20->addArg("--rowSkip", "1"); 
  check20->addArg("--columnSkip", "0");
  check20->addArg("--matchNames", "false");
  check20->addArg("--allowedMismatch", "0");


  RT_Check* check21 = rt->newCheck();
  check21->setExe("${apt_check_matrix}");
  check21->addArg("--gen", "${out_testname}/rma-bg.quant-norm.normalization-target.txt");
  check21->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/rma-bg.quant-norm.normalization-target.txt"); 
  check21->addArg("--epsilon", "0.0001");
  check21->addArg("--rowSkip", "1"); 
  check21->addArg("--columnSkip", "0");
  check21->addArg("--matchNames", "false");
  check21->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doHumanGeneSpfTest()
{
  RT_Test* rt;
  rt = new RT_Test("doHumanGeneSpfTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-sketch");
  cmd->addArg("--spf-file", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.spf");
  cmd->addArg("--qc-probesets", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.qcc");
  cmd->addArg("-m", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.mps");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${qt_celFileList2}");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-sketch.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/rma-sketch.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "0");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doHumanGeneSpfTest()
{
  RT_Test* rt;
  rt = new RT_Test("doHumanGeneSpfTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "rma-sketch");
  cmd->addArg("--spf-file", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.spf");
  cmd->addArg("--qc-probesets", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.qcc");
  cmd->addArg("-m", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.mps");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("--cel-files", "${celFileList2}");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/rma-sketch.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hugene10stv1_testname}/rma-sketch.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "0");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doHumanGeneKillListTest()
{
  if ( !Fs::dirExists("test-generated/doHumanGeneKillListTest") ) {
    Fs::mkdirPath("test-generated/doHumanGeneKillListTest", false);
  }
  RT_Test* rt;
  rt = new RT_Test("doHumanGeneKillListTest");
  //rt->setIntegrationTestBool(true);
  RT_Test* rt1;
  rt1 = new RT_Test("doHumanGeneKillListTest-phase1");
//rt1->setIntegrationTestBool(true);
  RT_Test* rt2;
  rt2 = new RT_Test("doHumanGeneKillListTest-phase2");
  //rt2->setIntegrationTestBool(true);
  RT_Test* rt3;
  rt3 = new RT_Test("doHumanGeneKillListTest-phase3");
  //rt3->setIntegrationTestBool(true);
  RT_Test* rt4;
  rt4 = new RT_Test("doHumanGeneKillListTest-phase4");
  //rt4->setIntegrationTestBool(true);
  RT_Test* rt5;
  rt5 = new RT_Test("doHumanGeneKillListTest-phase5");
  //rt5->setIntegrationTestBool(true);
  RT_Test* rt6;
  rt6 = new RT_Test("doHumanGeneKillListTest-phase6");
  //rt6->setIntegrationTestBool(true);

  //Util::makeDir("test-generated/doHumanGeneKillListTest");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_summarize}");
  cmd1->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd1->addArg("-c", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd1->addArg("-p", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd1->addArg("--kill-list", "${gold_idata_lib_hugene10stv1}/test-kill-list.txt");
  cmd1->addArg("-m", "${gold_idata_lib_hugene10stv1}/test-kill-list.mps");
  cmd1->addArg("-o", "test-generated/doHumanGeneKillListTest/mask1");
  cmd1->addArg("--cel-files", "${qt_celFileList2}");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_summarize}");
  cmd2->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd2->addArg("-c", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd2->addArg("-p", "${gold_idata_lib_hugene10stv1}/test-kill-list.pgf");
  cmd2->addArg("-m", "${gold_idata_lib_hugene10stv1}/test-kill-list2.mps");
  cmd2->addArg("-o", "test-generated/doHumanGeneKillListTest/mask2");
  cmd2->addArg("--cel-files", "${qt_celFileList2}");

  RT_Cmd* cmd3 = rt3->newCmd();
  cmd3->setExe("${apt_probeset_summarize}");
  cmd3->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd3->addArg("-c", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd3->addArg("-p", "${gold_idata_lib_hugene10stv1}/test-kill-list.pgf");
  cmd3->addArg("--kill-list", "${gold_idata_lib_hugene10stv1}/test-kill-list.txt");
  cmd3->addArg("-m", "${gold_idata_lib_hugene10stv1}/test-kill-list.mps");
  cmd3->addArg("-o", "test-generated/doHumanGeneKillListTest/mask3");
  cmd3->addArg("--force");
  cmd3->addArg("--cel-files", "${qt_celFileList2}");

  RT_Cmd* cmd4 = rt4->newCmd();
  cmd4->setExe("${apt_probeset_summarize}");
  cmd4->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd4->addArg("-c", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd4->addArg("-p", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd4->addArg("--kill-list", "${gold_idata_lib_hugene10stv1}/test-kill-list.xy.txt");
  cmd4->addArg("-m", "${gold_idata_lib_hugene10stv1}/test-kill-list.mps");
  cmd4->addArg("-o", "test-generated/doHumanGeneKillListTest/mask4");
  cmd4->addArg("--cel-files", "${qt_celFileList2}");

  RT_Cmd* cmd5 = rt5->newCmd();
  cmd5->setExe("${apt_probeset_summarize}");
  cmd5->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd5->addArg("-c", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd5->addArg("-p", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd5->addArg("--kill-list", "${gold_idata_lib_hugene10stv1}/test-kill-list.txt");
  cmd5->addArg("-o", "test-generated/doHumanGeneKillListTest/mask5");
  cmd5->addArg("--cel-files", "${qt_celFileList2}");

  RT_Cmd* cmd6 = rt6->newCmd();
  cmd6->setExe("${apt_probeset_summarize}");
  cmd6->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd6->addArg("-c", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd6->addArg("-p", "${gold_idata_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd6->addArg("-m", "${gold_idata_lib_hugene10stv1}/test-kill-list2.mps");
  cmd6->addArg("-o", "test-generated/doHumanGeneKillListTest/mask6");
  cmd6->addArg("--cel-files", "${qt_celFileList2}");

  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum}/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask2/quant-norm.pm-only.med-polish.summary.txt");
  check2->addArg("--gold", "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt3->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask3/quant-norm.pm-only.med-polish.summary.txt");
  check3->addArg("--gold", "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "true");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt4->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask4/quant-norm.pm-only.med-polish.summary.txt");
  check4->addArg("--gold", "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "true");
  check4->addArg("--allowedMismatch", "0");

  RT_Check* check5 = rt5->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask5/quant-norm.pm-only.med-polish.summary.txt");
  check5->addArg("--gold", "${gold_idata_psum}/doHumanGeneKillListTest/mask5/quant-norm.pm-only.med-polish.summary.txt"); 
  check5->addArg("--epsilon", "0.0001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "true");
  check5->addArg("--allowedMismatch", "0");

  RT_Check* check6 = rt6->newCheck();
  check6->setExe("${apt_check_matrix}");
  check6->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask6/quant-norm.pm-only.med-polish.summary.txt");
  check6->addArg("--gold", "${gold_idata_psum}/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check6->addArg("--epsilon", "0.0001");
  check6->addArg("--rowSkip", "1"); 
  check6->addArg("--columnSkip", "1");
  check6->addArg("--matchNames", "true");
  check6->addArg("--allowedMismatch", "0");
  //check6->setNegTest();

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  rt->addSubtest(rt3);
  rt->addSubtest(rt4);
  rt->addSubtest(rt5);
  rt->addSubtest(rt6);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doHumanGeneKillListTest()
{
  if ( !Fs::dirExists("test-generated/doHumanGeneKillListTest") ) {
    Fs::mkdirPath("test-generated/doHumanGeneKillListTest", false);
  }
  RT_Test* rt;
  rt = new RT_Test("doHumanGeneKillListTest");
  RT_Test* rt1;
  rt1 = new RT_Test("doHumanGeneKillListTest-phase1");
  RT_Test* rt2;
  rt2 = new RT_Test("doHumanGeneKillListTest-phase2");
  RT_Test* rt3;
  rt3 = new RT_Test("doHumanGeneKillListTest-phase3");
  RT_Test* rt4;
  rt4 = new RT_Test("doHumanGeneKillListTest-phase4");
  RT_Test* rt5;
  rt5 = new RT_Test("doHumanGeneKillListTest-phase5");
  RT_Test* rt6;
  rt6 = new RT_Test("doHumanGeneKillListTest-phase6");

  //Util::makeDir("test-generated/doHumanGeneKillListTest");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_probeset_summarize}");
  cmd1->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd1->addArg("-c", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd1->addArg("-p", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd1->addArg("--kill-list", "${gold_lib_hugene10stv1}/test-kill-list.txt");
  cmd1->addArg("-m", "${gold_lib_hugene10stv1}/test-kill-list.mps");
  cmd1->addArg("-o", "test-generated/doHumanGeneKillListTest/mask1");
  cmd1->addArg("--cel-files", "${celFileList2}");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_probeset_summarize}");
  cmd2->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd2->addArg("-c", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd2->addArg("-p", "${gold_lib_hugene10stv1}/test-kill-list.pgf");
  cmd2->addArg("-m", "${gold_lib_hugene10stv1}/test-kill-list2.mps");
  cmd2->addArg("-o", "test-generated/doHumanGeneKillListTest/mask2");
  cmd2->addArg("--cel-files", "${celFileList2}");

  RT_Cmd* cmd3 = rt3->newCmd();
  cmd3->setExe("${apt_probeset_summarize}");
  cmd3->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd3->addArg("-d", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.cdf");
  cmd3->addArg("--kill-list", "${gold_lib_hugene10stv1}/test-kill-list.txt");
  cmd3->addArg("-m", "${gold_lib_hugene10stv1}/test-kill-list.mps");
  cmd3->addArg("-o", "test-generated/doHumanGeneKillListTest/mask3");
  cmd3->addArg("--force");
  cmd3->addArg("--cel-files", "${celFileList2}");

  RT_Cmd* cmd4 = rt4->newCmd();
  cmd4->setExe("${apt_probeset_summarize}");
  cmd4->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd4->addArg("-c", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd4->addArg("-p", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd4->addArg("--kill-list", "${gold_lib_hugene10stv1}/test-kill-list.xy.txt");
  cmd4->addArg("-m", "${gold_lib_hugene10stv1}/test-kill-list.mps");
  cmd4->addArg("-o", "test-generated/doHumanGeneKillListTest/mask4");
  cmd4->addArg("--cel-files", "${celFileList2}");

  RT_Cmd* cmd5 = rt5->newCmd();
  cmd5->setExe("${apt_probeset_summarize}");
  cmd5->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd5->addArg("-c", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd5->addArg("-p", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd5->addArg("--kill-list", "${gold_lib_hugene10stv1}/test-kill-list.txt");
  cmd5->addArg("-o", "test-generated/doHumanGeneKillListTest/mask5");
  cmd5->addArg("--cel-files", "${celFileList2}");

  RT_Cmd* cmd6 = rt6->newCmd();
  cmd6->setExe("${apt_probeset_summarize}");
  cmd6->addArg("-a", "quant-norm.sketch=-1,pm-only,med-polish");
  cmd6->addArg("-c", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.clf");
  cmd6->addArg("-p", "${gold_lib_hugene10stv1}/HuGene-1_0-st-v1.r3.pgf");
  cmd6->addArg("-m", "${gold_lib_hugene10stv1}/test-kill-list2.mps");
  cmd6->addArg("-o", "test-generated/doHumanGeneKillListTest/mask6");
  cmd6->addArg("--cel-files", "${celFileList2}");

  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hugene10stv1}/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask2/quant-norm.pm-only.med-polish.summary.txt");
  check2->addArg("--gold", "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  RT_Check* check3 = rt3->newCheck();
  check3->setExe("${apt_check_matrix}");
  check3->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask3/quant-norm.pm-only.med-polish.summary.txt");
  check3->addArg("--gold", "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--rowSkip", "1"); 
  check3->addArg("--columnSkip", "1");
  check3->addArg("--matchNames", "true");
  check3->addArg("--allowedMismatch", "0");

  RT_Check* check4 = rt4->newCheck();
  check4->setExe("${apt_check_matrix}");
  check4->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask4/quant-norm.pm-only.med-polish.summary.txt");
  check4->addArg("--gold", "test-generated/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check4->addArg("--epsilon", "0.0001");
  check4->addArg("--rowSkip", "1"); 
  check4->addArg("--columnSkip", "1");
  check4->addArg("--matchNames", "true");
  check4->addArg("--allowedMismatch", "0");

  RT_Check* check5 = rt5->newCheck();
  check5->setExe("${apt_check_matrix}");
  check5->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask5/quant-norm.pm-only.med-polish.summary.txt");
  check5->addArg("--gold", "${gold_probesetsummarize_hugene10stv1}/doHumanGeneKillListTest/mask5/quant-norm.pm-only.med-polish.summary.txt"); 
  check5->addArg("--epsilon", "0.0001");
  check5->addArg("--rowSkip", "1"); 
  check5->addArg("--columnSkip", "1");
  check5->addArg("--matchNames", "true");
  check5->addArg("--allowedMismatch", "0");

  RT_Check* check6 = rt6->newCheck();
  check6->setExe("${apt_check_matrix}");
  check6->addArg("--gen", "test-generated/doHumanGeneKillListTest/mask6/quant-norm.pm-only.med-polish.summary.txt");
  check6->addArg("--gold", "${gold_probesetsummarize_hugene10stv1}/doHumanGeneKillListTest/mask1/quant-norm.pm-only.med-polish.summary.txt"); 
  check6->addArg("--epsilon", "0.0001");
  check6->addArg("--rowSkip", "1"); 
  check6->addArg("--columnSkip", "1");
  check6->addArg("--matchNames", "true");
  check6->addArg("--allowedMismatch", "0");
  check6->setNegTest();

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  rt->addSubtest(rt3);
  rt->addSubtest(rt4);
  rt->addSubtest(rt5);
  rt->addSubtest(rt6);
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doSNP6CopyNumberWorkflow1a()
{
  RT_Test* rt;
  rt = new RT_Test("doSNP6CnWf1a");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp6cn}/GenomeWideSNP_6.spf");
  cmd->addArg("--analysis", "quant-norm.sketch=-1,pm-sum,med-polish,expr.genotype=true.allele-a=true");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA04626_3X.rep3.CEL");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA10851_1X_fosmid_ref.rep1.CEL");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA15510_2X_fosmid.rep2.CEL");
  
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-sum.med-polish.expr.report.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-sum.med-polish.expr.report.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-sum.med-polish.expr.summary.txt");
  check2->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-sum.med-polish.expr.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doSNP6CopyNumberWorkflow1a()
{
  RT_Test* rt;
  rt = new RT_Test("doSNP6CopyNumberWorkflow1a");


  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--analysis", "quant-norm.sketch=-1,pm-sum,med-polish,expr.genotype=true.allele-a=true");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA04626_3X.rep3.CEL");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA10851_1X_fosmid_ref.rep1.CEL");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA15510_2X_fosmid.rep2.CEL");
  
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-sum.med-polish.expr.report.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_genomewidesnp6_testname}/quant-norm.pm-sum.med-polish.expr.report.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-sum.med-polish.expr.summary.txt");
  check2->addArg("--gold", "${gold_probesetsummarize_genomewidesnp6_testname}/quant-norm.pm-sum.med-polish.expr.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doSNP6CopyNumberWorkflow1b()
{
  RT_Test* rt;
  rt = new RT_Test("doSNP6CnWf1b");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp6cn}/GenomeWideSNP_6.spf");
  cmd->addArg("--analysis", "quant-norm.sketch=-1,pm-only,med-polish,expr.genotype=true");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA04626_3X.rep3.CEL");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA10851_1X_fosmid_ref.rep1.CEL");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA15510_2X_fosmid.rep2.CEL");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.expr.report.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.expr.report.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.expr.summary.txt");
  check2->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.expr.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doSNP6CopyNumberWorkflow1b()
{
  RT_Test* rt;
  rt = new RT_Test("doSNP6CopyNumberWorkflow1b");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--analysis", "quant-norm.sketch=-1,pm-only,med-polish,expr.genotype=true");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA04626_3X.rep3.CEL");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA10851_1X_fosmid_ref.rep1.CEL");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA15510_2X_fosmid.rep2.CEL");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.expr.report.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_genomewidesnp6_testname}/quant-norm.pm-only.med-polish.expr.report.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.expr.summary.txt");
  check2->addArg("--gold", "${gold_probesetsummarize_genomewidesnp6_testname}/quant-norm.pm-only.med-polish.expr.summary.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doSNP6CopyNumberWorkflow2()
{
  RT_Test* rt;
  rt = new RT_Test("doSNP6CnWf2");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("--spf-file", "${gold_idata_lib_genomewidesnp6cn}/GenomeWideSNP_6.spf");
  cmd->addArg("--analysis", "quant-norm.sketch=-1,pm-only,med-polish,pca-select");
  cmd->addArg("--meta-probesets", "${gold_idata_lib_genomewidesnp6cn}/GenomeWideSNP_6.na22.dgv-cnvMay07.mps");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA04626_3X.rep3.CEL");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA10851_1X_fosmid_ref.rep1.CEL");
  cmd->addArg("${gold_idata_cel_genomewidesnp6cn}/NA15510_2X_fosmid.rep2.CEL");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.report.txt");
  check2->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.report.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  //MIXED FILE CHECK
  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_mixedfile}");
  check3->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt");
  check3->addArg("--gold", "${gold_idata_psum_testname}/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt");
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--skipLines", "0");
  check3->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doSNP6CopyNumberWorkflow2()
{
  RT_Test* rt;
  rt = new RT_Test("doSNP6CopyNumberWorkflow2");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--analysis", "quant-norm.sketch=-1,pm-only,med-polish,pca-select");
  cmd->addArg("--meta-probesets", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.na22.dgv-cnvMay07.mps");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA04626_3X.rep3.CEL");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA10851_1X_fosmid_ref.rep1.CEL");
  cmd->addArg("${gold_cel_genomewidesnp6}/NA15510_2X_fosmid.rep2.CEL");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_genomewidesnp6_testname}/quant-norm.pm-only.med-polish.pca-select.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_matrix}");
  check2->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.report.txt");
  check2->addArg("--gold", "${gold_probesetsummarize_genomewidesnp6_testname}/quant-norm.pm-only.med-polish.pca-select.report.txt"); 
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--rowSkip", "1"); 
  check2->addArg("--columnSkip", "1");
  check2->addArg("--matchNames", "true");
  check2->addArg("--allowedMismatch", "0");

  //MIXED FILE CHECK
  RT_Check* check3 = rt->newCheck();
  check3->setExe("${apt_check_mixedfile}");
  check3->addArg("--gen", "${out_testname}/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt");
  check3->addArg("--gold", "${gold_probesetsummarize_genomewidesnp6_testname}/quant-norm.pm-only.med-polish.pca-select.pca-select.report.txt");
  check3->addArg("--epsilon", "0.0001");
  check3->addArg("--skipLines", "0");
  check3->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doPlierMMTissueMedianNormTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierMMTissueMedianNormTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "med-norm.lowprecision=true,pm-mm,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd->addArg("-p", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullTissueCels = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fullTissueCels.begin();i < fullTissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/med-norm.pm-mm.plier.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/med-norm.pm-mm.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doPlierMMTissueMedianNormTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPlierMMTissueMedianNormTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "med-norm.lowprecision=true,pm-mm,plier.FixPrecomputed=false.SafetyZero=0.NumericalTolerance=0");
  cmd->addArg("-p", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullTissueCels = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fullTissueCels.begin();i < fullTissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/med-norm.pm-mm.plier.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/med-norm.pm-mm.plier.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doDabgSubsetTest()
{
  RT_Test* rt;
  rt = new RT_Test("doDabgSubsetTest");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "pm-only,dabg");
  cmd->addArg("-c", "${gold_idata_lib_huex10stv2}/HuEx-1_0-st-v2.clf");
  cmd->addArg("-p", "${gold_idata_lib_huex10stv2}/HuEx-1_0-st-v2.pgf");
  cmd->addArg("-b", "${gold_idata_lib_huex10stv2}/antigenomic.bgp");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("-cc-chp-output");
  cmd->addArg("-s", "${gold_idata_lib_huex10stv2}/dabg.subset.txt");
  cmd->addArg("-x", "5 ");

  std::vector<std::string> fullHuexCels = Util::addPrefixSuffix(huexCels, "${gold_idata_cel_huex10stv2}/", "${celsuffixUp}");
  for(std::vector<std::string>::iterator i = fullHuexCels.begin();i < fullHuexCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  string chpFilesSuffix = ".pm-only.dabg.chp";
  const char *chpFiles[] = {"huex_wta_cb_A", "huex_wta_cb_B", "huex_wta_cb_C",
			    "huex_wta_heart_A", "huex_wta_heart_B", "huex_wta_heart_C",
			    "huex_wta_muscle_A", "huex_wta_muscle_B", "huex_wta_muscle_C",
			    "huex_wta_testes_A", "huex_wta_testes_B", "huex_wta_testes_C",
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_idata_psum_testname}/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", chpFilesSuffix);

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/pm-only.dabg.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/pm-only.dabg.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doDabgSubsetTest()
{
  RT_Test* rt;
  rt = new RT_Test("doDabgSubsetTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "pm-only,dabg");
  cmd->addArg("-c", "${gold_lib_huex10stv2}/HuEx-1_0-st-v2.clf");
  cmd->addArg("-p", "${gold_lib_huex10stv2}/HuEx-1_0-st-v2.pgf");
  cmd->addArg("-b", "${gold_lib_huex10stv2}/antigenomic.probe_id-only.bgp");
  cmd->addArg("-o", "${out_testname}");
  cmd->addArg("-cc-chp-output");
  cmd->addArg("-s", "${gold_lib_huex10stv2}/dabg.subset.txt");
  cmd->addArg("-x", "5 ");

  std::vector<std::string> fullHuexCels = Util::addPrefixSuffix(huexCels, "${gold_cel_huex10stv2}/", "${celsuffixUp}");
  for(std::vector<std::string>::iterator i = fullHuexCels.begin();i < fullHuexCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  string chpFilesSuffix = ".pm-only.dabg.chp";
  const char *chpFiles[] = {"huex_wta_cb_A", "huex_wta_cb_B", "huex_wta_cb_C",
			    "huex_wta_heart_A", "huex_wta_heart_B", "huex_wta_heart_C",
			    "huex_wta_muscle_A", "huex_wta_muscle_B", "huex_wta_muscle_C",
			    "huex_wta_testes_A", "huex_wta_testes_B", "huex_wta_testes_C",
			    NULL
  };

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(chpFiles, "${gold_probesetsummarize_huex10stv2_testname}/cc-chp/", chpFilesSuffix);
  gen = Util::addPrefixSuffix(chpFiles, "${out_testname}/cc-chp/", chpFilesSuffix);

  //CALVIN CHP CHECK
  //Generate Check objects for CalvinChpCheck and store them in rt
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_calvinchp}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/pm-only.dabg.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_huex10stv2_testname}/pm-only.dabg.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "true");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* qt_doDabgU133Test()
{
  RT_Test* rt;
  rt = new RT_Test("doDabgU133Test");
  //rt->setIntegrationTestBool(true);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "pm-only,dabg");
  cmd->addArg("-b", "${gold_idata_lib_hgu133plus2}/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp");
  cmd->addArg("-p", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_idata_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullTissueCels = Util::addPrefixSuffix(tissueCels, "${gold_idata_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fullTissueCels.begin();i < fullTissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }
  
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/pm-only.dabg.summary.txt");
  check1->addArg("--gold", "${gold_idata_psum_testname}/pm-only.dabg.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doDabgU133Test()
{
  RT_Test* rt;
  rt = new RT_Test("doDabgU133Test");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_probeset_summarize}");
  cmd->addArg("-a", "pm-only,dabg");
  cmd->addArg("-b", "${gold_lib_hgu133plus2}/pooled-mm-probes.rand-1000-per-bin.probe_id-only.bgp");
  cmd->addArg("-p", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("-c", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.clf");
  cmd->addArg("-x", "5");
  cmd->addArg("-o", "${out_testname}");

  std::vector<std::string> fullTissueCels = Util::addPrefixSuffix(tissueCels, "${gold_cel_hgu133plus2}/", "${celsuffixLower}");
  for(std::vector<std::string>::iterator i = fullTissueCels.begin();i < fullTissueCels.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }
  
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_matrix}");
  check1->addArg("--gen", "${out_testname}/pm-only.dabg.summary.txt");
  check1->addArg("--gold", "${gold_probesetsummarize_hgu133plus2_testname}/pm-only.dabg.summary.txt"); 
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--rowSkip", "1"); 
  check1->addArg("--columnSkip", "1");
  check1->addArg("--matchNames", "false");
  check1->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
