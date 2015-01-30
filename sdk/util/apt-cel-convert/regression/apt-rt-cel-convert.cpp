////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////
/*
#include "rtest/RT_Args.h"
#include "rtest/RT_Check.h"
#include "rtest/RT_Cmd.h"
#include "rtest/RT_Main.h"
#include "rtest/RT_Test.h"
#include "rtest/RT_Util.h"
*/
#include "rtest/RT_Main.h"
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


const char *xdaCels[] = {
    "heart-rep1",
    "heart-rep2",
    NULL
};

const char *textCels[] = {
    "heart-rep1.text",
    "heart-rep2.text",
    NULL
};

const char *agccCels[] = {
    "heart-rep1.agcc",
    "heart-rep2.agcc",
    NULL
};


RT_Test* doConvertXDAtoAGCC();
RT_Test* doConvertXDAtoTEXT();
RT_Test* doConvertAGCCtoXDA();
RT_Test* doConvertAGCCtoTEXT();
RT_Test* doConvertTEXTtoXDA();
RT_Test* doConvertTEXTtoAGCC();
RT_Test* doConvertRoundRobin();

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if ( !Fs::dirExists("test-generated") ) {
    Fs::mkdir("test-generated", false);
  }

  RT_Main* rMain;
  rMain = new RT_Main("apt-rt-cel-convert");

  //Store names of tests for Help Message
  //Store it in our RT_Main Instance
  std::vector<std::string> testNameList;
  testNameList.push_back("doConvertXDAtoAGCC");
  testNameList.push_back("doConvertXDAtoTEXT");
  testNameList.push_back("doConvertAGCCtoXDA");
  testNameList.push_back("doConvertAGCCtoTEXT");
  testNameList.push_back("doConvertTEXTtoXDA");
  testNameList.push_back("doConvertTEXTtoAGCC");
  testNameList.push_back("doConvertRoundRobin");
  rMain->setTestNameList(testNameList);
  
  //Set Variables
  rMain->setVar("sdktop", "../../..");
  rMain->setVar("gold_cel_hgu133plus2", "${gold_dir}/cel/HG-U133_Plus_2");
  rMain->setVar("celsuffix", ".cel");
  rMain->setVar("out_heartrep1_cel", "${out_dir}/heart-rep1${celsuffix}");
  rMain->setVar("gold_heartrep1_cel", "${gold_cel_hgu133plus2}/heart-rep1${celsuffix}");
  
  if(argc < 1)
    {
      cout<<"Error: Incorrect number of arguments given"<<endl;
      exit(1);
    }
  //They want us to run all default tests
  else if(argc == 1)
    {
      //RUN ALL TESTS OMG
      cout<<"Running all default Regression Tests!"<<endl;
      rMain->addTest(doConvertRoundRobin());
      rMain->addTest(doConvertXDAtoAGCC());
      rMain->addTest(doConvertXDAtoTEXT());
      rMain->addTest(doConvertAGCCtoXDA());
      rMain->addTest(doConvertAGCCtoTEXT());
      rMain->addTest(doConvertTEXTtoXDA());
      rMain->addTest(doConvertTEXTtoAGCC());
    }
  else
    {
      rMain->parseArgv(argc, argv);
 
     
      cout<<"Running Built-In Regression Tests"<<endl;
	 
      if(rMain->getVarVal("doConvertRoundRobin") != "")
	{
	  rMain->addTest(doConvertRoundRobin());
	}
      if(rMain->getVarVal("doConvertXDAtoAGCC") != "")
	{
	  rMain->addTest(doConvertXDAtoAGCC());
	}
      if(rMain->getVarVal("doConvertXDAtoTEXT") != "")
	{
	  rMain->addTest(doConvertXDAtoTEXT());
	}
      if(rMain->getVarVal("doConvertAGCCtoXDA") != "")
	{
	  rMain->addTest(doConvertAGCCtoXDA());
	}
      if(rMain->getVarVal("doConvertAGCCtoTEXT") != "")
	{
	  rMain->addTest(doConvertAGCCtoTEXT());
	}
      if(rMain->getVarVal("doConvertTEXTtoXDA") != "")
	{
	  rMain->addTest(doConvertTEXTtoXDA());
	}
      if(rMain->getVarVal("doConvertTEXTtoAGCC") != "")
	{
	  rMain->addTest(doConvertTEXTtoAGCC());
	}
    }

  rMain->runTests();
}


///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* runTestHarness(const char *cels[], const std::string &inFormat, const std::string &outFormat)
{
  std::string testName = "doConvert_" + inFormat + "_to_" + outFormat;

  RT_Test* rt;
  rt = new RT_Test(testName);
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_convert}");
  cmd->addArg("--out-dir");
  cmd->addArg("${out_dir}");
  cmd->addArg("--format");
  cmd->addArg(outFormat);

  std::vector<std::string> fullCelFiles = Util::addPrefixSuffix(cels, "${gold_cel_hgu133plus2}/", "${celsuffix}");
  for(std::vector<std::string>::iterator i = fullCelFiles.begin();i < fullCelFiles.end();i++)
    {
      std::string stringArg = *i;
      cmd->addArg(stringArg);
    }
  
  std::vector<std::string> gold,gen;
  gen = Util::addPrefixSuffix(cels, "${out_dir}/", "${celsuffix}");
  gold = Util::addPrefixSuffix(cels, "${gold_cel_hgu133plus2}/", "${celsuffix}");
 
  //Cel Check
  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_cel}");
  check->addArg("--epsilon", "0.00001");
  check->setGenFiles(gen);
  check->setGoldFiles(gold);
  //rt->addCheckMult("${apt_check_cel}", gen, gold, checkArgs1); 
  
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doConvertXDAtoAGCC()  
{ 
  RT_Test* rt;
  rt = runTestHarness(xdaCels,"xda", "agcc"); 
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doConvertXDAtoTEXT()  
{ 
  RT_Test* rt;
  rt = runTestHarness(xdaCels,"xda", "text"); 
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doConvertAGCCtoXDA()  
{
  RT_Test* rt;
  rt = runTestHarness(agccCels,"agcc","xda" ); 
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doConvertAGCCtoTEXT() 
{
  RT_Test* rt;
  rt = runTestHarness(agccCels,"agcc","text"); 
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doConvertTEXTtoXDA()  
{ 
  RT_Test* rt;
  rt = runTestHarness(textCels,"text","xda" ); 
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doConvertTEXTtoAGCC() 
{
  RT_Test* rt;
  rt = runTestHarness(textCels,"text","agcc"); 
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doConvertRoundRobin() 
{
  RT_Test* rt;
  rt = new RT_Test("doConvertRoundRobin"); 
  RT_Test* rt1;
  rt1 = new RT_Test("doConvertRoundRobin-part1");
  RT_Test* rt2;
  rt2 = new RT_Test("doConvertRoundRobin-part2");
  RT_Test* rt3;
  rt3 = new RT_Test("doConvertRoundRobin-part3");

  RT_Cmd* cmd1 = rt1->newCmd();
  cmd1->setExe("${apt_cel_convert}");
  cmd1->addArg("--out-dir", "${out_dir}"); 
  cmd1->addArg("--format", "text"); 
  cmd1->addArg("${gold_heartrep1_cel}"); 

  RT_Check* check1 = rt1->newCheck();
  check1->setExe("${apt_check_cel}");
  check1->addArg("--gen", "${out_heartrep1_cel}");
  check1->addArg("--gold", "${gold_heartrep1_cel}");
  check1->addArg("--epsilon", "0.00001");

  RT_Cmd* cmd2 = rt2->newCmd();
  cmd2->setExe("${apt_cel_convert}");
  cmd2->addArg("--in-place");
  cmd2->addArg("--format", "agcc");
  cmd2->addArg("${out_heartrep1_cel}");
  
  RT_Check* check2 = rt2->newCheck();
  check2->setExe("${apt_check_cel}");
  check2->addArg("--gen", "${out_heartrep1_cel}");
  check2->addArg("--gold", "${gold_heartrep1_cel}");
  check2->addArg("--epsilon", "0.00001");

  RT_Cmd* cmd3 = rt3->newCmd();
  cmd3->setExe("${apt_cel_convert}");
  cmd3->addArg("--in-place");
  cmd3->addArg("--format", "xda");
  cmd3->addArg("${out_heartrep1_cel}");

  RT_Check* check3 = rt3->newCheck();
  check3->setExe("${apt_check_cel}");
  check3->addArg("--gen", "${out_heartrep1_cel}");
  check3->addArg("--gold", "${gold_heartrep1_cel}");
  check3->addArg("--epsilon", "0.00001");

  rt->addSubtest(rt1);
  rt->addSubtest(rt2);
  rt->addSubtest(rt3);
  return rt;
}
