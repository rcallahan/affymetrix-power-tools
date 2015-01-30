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

const char *celFiles[] = {"heart-rep1","heart-rep2","heart-rep3",
			  NULL
};

RT_Test* doSpfTest();
RT_Test* doAgccTest();
RT_Test* doCdfTest();
RT_Test* doPgfTest();
RT_Test* doA5SketchTest();


///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

  RT_Main* rMain;
  rMain = new RT_Main("apt-rt-cel-transformer");
  
  if ( !Fs::dirExists("test-generated") ) {
    Fs::mkdir("test-generated", false);
  }

  rMain->setVar("sdk_top", "../../..");
  rMain->setVar("celsuffix", ".cel");
  rMain->setVar("agcccelsuffix", ".agcc.cel");
  rMain->setVar("gold_lib_hgu133plus2", "${gold_dir}/lib/HG-U133_Plus_2");
  rMain->setVar("gold_cel_hgu133plus2", "${gold_dir}/cel/HG-U133_Plus_2");
  rMain->setVar("gold_celtransformer_hgu133plus2", "${gold_dir}/chipstream/cel-transformer/HG-U133_Plus_2");

  if(argc < 1)
    {
      cout<<"Error: Incorrect number of arguments given"<<endl;
      exit(1);
    }
  //They want us to run all default tests
  else if(argc == 1)
    {
      cout<<"Running All Default Regression Tests"<<endl;
      rMain->addTest(doSpfTest());
      rMain->addTest(doAgccTest());
      rMain->addTest(doCdfTest());
      rMain->addTest(doPgfTest());
      rMain->addTest(doA5SketchTest());
    }
  else
    {
      rMain->parseArgv(argc, argv);

      cout<<"Running Built-In Regression Tests"<<endl;
	 
      if(rMain->getVarVal("doSpfTest") != "")
	{
	  rMain->addTest(doSpfTest());
	}
      if(rMain->getVarVal("doAgccTest") != "")
	{
	  rMain->addTest(doAgccTest());
	}
      if(rMain->getVarVal("doCdfTest") != "")
	{
	  rMain->addTest(doCdfTest());
	}
      if(rMain->getVarVal("doPgfTest") != "")
	{
	  rMain->addTest(doPgfTest());
	}
      if(rMain->getVarVal("doA5SketchTest") != "")
	{
	  rMain->addTest(doA5SketchTest());
	}
    }
  rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doSpfTest()
{
  RT_Test* rt;
  rt = new RT_Test("doSpfTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_transformer}");
  cmd->addArg("-c", "rma-bg,quant-norm");
  cmd->addArg("-o", "${out_dir}");
  cmd->addArg("--spf-file", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.spf");

  std::vector<std::string> fullCelFiles = Util::addPrefixSuffix(celFiles, " ${gold_cel_hgu133plus2}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullCelFiles.begin();i < fullCelFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(celFiles, "${gold_celtransformer_hgu133plus2}/cdf-test/", "${celsuffix}");
  gen = Util::addPrefixSuffix(celFiles, "${out_dir}/", "${celsuffix}");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_cel}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  check->addArg("--epsilon", "0.2");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doAgccTest()
{
  RT_Test* rt;
  rt = new RT_Test("doAgccTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_transformer}");
  cmd->addArg("-c", "rma-bg,quant-norm");
  cmd->addArg("-o", "${out_dir}");
  cmd->addArg("--cdf-file", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");

  std::vector<std::string> fullCelFiles = Util::addPrefixSuffix(celFiles, " ${gold_cel_hgu133plus2}/", "${agcccelsuffix}");
  for(std::vector<std::string>::iterator i = fullCelFiles.begin();i < fullCelFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(celFiles, "${gold_celtransformer_hgu133plus2}/cdf-test/", "${celsuffix}");
  gen = Util::addPrefixSuffix(celFiles, "${out_dir}/", "${agcccelsuffix}");
  
  //Cel Check
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_cel}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  check->addArg("--epsilon", "0.2");
  
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doCdfTest()
{
  RT_Test* rt;
  rt = new RT_Test("doCdfTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_transformer}");
  cmd->addArg("-c", "rma-bg,quant-norm");
  cmd->addArg("-o", "${out_dir}");
  cmd->addArg("--cdf-file", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");

  std::vector<std::string> fullCelFiles = Util::addPrefixSuffix(celFiles, "${gold_cel_hgu133plus2}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullCelFiles.begin();i < fullCelFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  vector<string> gold,gen;
  gold = Util::addPrefixSuffix(celFiles, "${gold_celtransformer_hgu133plus2}/cdf-test/", "${celsuffix}");
  gen = Util::addPrefixSuffix(celFiles, "${out_dir}/", "${celsuffix}");

  //CEL Check
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_cel}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  check->addArg("--epsilon", "0.2");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doPgfTest()
{
  RT_Test* rt;
  rt = new RT_Test("doPgfTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_transformer}");
  cmd->addArg("-c", "rma-bg,quant-norm");
  cmd->addArg("-o", "${out_dir}");
  cmd->addArg("--pgf-file", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.pgf");
  cmd->addArg("--clf-file", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.clf");

  std::vector<std::string> fullCelFiles = Util::addPrefixSuffix(celFiles, " ${gold_cel_hgu133plus2}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullCelFiles.begin();i < fullCelFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  
  std::vector<std::string> gold,gen;
  gold = Util::addPrefixSuffix(celFiles, "${gold_celtransformer_hgu133plus2}/pgf-test/", "${celsuffix}");
  gen = Util::addPrefixSuffix(celFiles, "${out_dir}/", "${celsuffix}");

  //CEL Check
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_cel}");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  check->addArg("--epsilon", "0.2");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doA5SketchTest()
{
  RT_Test* rt;
  rt = new RT_Test("doA5SketchTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_transformer}");
  cmd->addArg("-c", "quant-norm.sketch=1000");
  cmd->addArg("-o", "${out_dir}");
  cmd->addArg("--write-sketch");
  cmd->addArg("--a5-write-sketch");

  std::vector<std::string> fullCelFiles = Util::addPrefixSuffix(celFiles, " ${gold_cel_hgu133plus2}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullCelFiles.begin();i < fullCelFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }

  ///@todo add check for cel file output
  ///@todo add check for txt sketch output file
  ///@todo add check for a5 sketch output file

  return rt;
}
