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

//
#include "rtest/RT_Args.h"
#include "rtest/RT_Check.h"
#include "rtest/RT_Cmd.h"
#include "rtest/RT_Main.h"
#include "rtest/RT_Test.h"
#include "rtest/RT_Util.h"
//
#include "util/Fs.h"
#include "util/Util.h"
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


const std::string CMD_EXE = "../../../output/amd64-pc-linux/bin/apt-chp-to-txt";
const std::string CALVIN_CHP_CHECK = "../../../output/amd64-pc-linux/bin/apt-check-calvinchp";
const std::string CHP_CHECK = "../../../output/amd64-pc-linux/bin/apt-check-chp";
const std::string MATRIX_CHECK = "../../../output/amd64-pc-linux/bin/apt-check-matrix";
const std::string MIXED_FILE_CHECK = "../../../output/amd64-pc-linux/bin/apt-check-mixedfile";
const std::string CEL_CHECK = "../../../output/amd64-pc-linux/bin/apt-check-cel";
const std::string TEXT_FILE_CHECK = "../../../output/amd64-pc-linux/bin/apt-check-textfile";


RT_Test* doCdfXdaGenoMpam();
RT_Test* doCdfXdaGenoBrlmm();
RT_Test* doCdfXdaGenoDm();
RT_Test* doCdfXdaExpMas5();
RT_Test* doCdfXdaGenoMpamBad();
RT_Test* doTest(const std::string &testName, const std::string &fileName, const std::string &cdf = "");
RT_Test* doBadChipTypeTest(const std::string &testName, const std::string &fileName, const std::string &cdf);
///@todo regression test for "no-header" option
///@todo regression test for "no-body" option
///@todo regression test for "chp-files" option

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

 RT_Main* rMain;
 rMain = new RT_Main("apt-rt-chp-to-txt");

 if ( !Fs::dirExists("test-generated") ) {
   Fs::mkdir("test-generated", false);
 }


 rMain->setVar("sdk_top", "../../..");
 rMain->setVar("gold_lib", "${gold_dir}/lib");

 if(argc < 1)
   {
     cout<<"Error: Incorrect number of arguments given"<<endl;
     exit(1);
   }
 //They want us to run all default tests
 else if(argc == 1)
   {
     rMain->addTest(doCdfXdaGenoMpam());
     rMain->addTest(doCdfXdaGenoBrlmm());
     rMain->addTest(doCdfXdaGenoDm());
     rMain->addTest(doCdfXdaExpMas5());
     rMain->addTest(doCdfXdaGenoMpamBad());
   }
 else
   {
     rMain->parseArgv(argc, argv);

     cout<<"Running Built-In Regression Tests"<<endl;

     if(rMain->getVarVal("doCdfXdaGenoMpam") != "")
       {
	 rMain->addTest(doCdfXdaGenoMpam());
       }
     if(rMain->getVarVal("doCdfXdaGenoBrlmm") != "")
       {
	 rMain->addTest(doCdfXdaGenoBrlmm());
       }
     if(rMain->getVarVal("doCdfXdaGenoDm") != "")
       {
	 rMain->addTest(doCdfXdaGenoDm());
       }
     if(rMain->getVarVal("doCdfXdaExpMas5") != "")
       {
	 rMain->addTest(doCdfXdaExpMas5());
       }
     if(rMain->getVarVal("doCdfXdaGenoMpamBad") != "")
       {
	 rMain->addTest(doCdfXdaGenoMpamBad());
       }
   }
 rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doTest(const std::string &testName, const std::string &fileName, const std::string &cdfFile)
{
  RT_Test* rt;
  rt = new RT_Test(testName);

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_chp_to_txt}");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("--cdf-file", cdfFile);
  cmd->addArg("${gold_dir}/util/chp-to-txt/chp/" + fileName);

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_mixedfile}");
  check->addArg("--gen", "${out_testname}/" + fileName + ".txt");
  check->addArg("--gold", Fs::join("${gold_dir}/util/chp-to-txt",testName,fileName + ".txt"));
  check->addArg("--epsilon", "0.00001");
  check->addArg("--skipLines", "0");
  check->addArg("--allowedMismatch", "0");
  
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doBadChipTypeTest(const std::string &testName, const std::string &fileName, const std::string &cdfFile)
{
  RT_Test* rt;
  rt = new RT_Test(testName);
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_chp_to_txt}");
  cmd->addArg("-v");
  cmd->addArg("-1");
  cmd->addArg("--out-dir", "${out_testname}");
  cmd->addArg("--cdf-file", cdfFile);
  cmd->addArg("${gold_dir}/util/chp-to-txt/chp/" + fileName);

  //Negative Test
  rt->setExpectedRetVal(1);

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doCdfXdaGenoMpam()
{
  RT_Test* rt;
  rt = doTest("doCdfXdaGenoMpam","Mapping10K_Xba131.xda.genotyping.mpam.chp", "${gold_lib}/Mapping10K_Xba131/Mapping10K_Xba131.CDF");
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doCdfXdaGenoBrlmm()
{
  RT_Test* rt;
  rt = doTest("doCdfXdaGenoBrlmm","Mapping250K_Sty.xda.genotyping.brlmm.chp", "${gold_lib}/Mapping250K_Sty/Mapping250K_Sty.cdf");
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doCdfXdaGenoDm()
{
  RT_Test* rt;
  rt = doTest("doCdfXdaGenoDm","Mapping50K_Hind240.xda.genotyping.dm.chp", "${gold_lib}/Mapping50K_Hind240/Mapping50K_Hind240.cdf");
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doCdfXdaExpMas5()
{
  RT_Test* rt;
  rt = doTest("doCdfXdaExpMas5","HG-U133_Plus_2.xda.expression.mas5.chp", "${gold_lib}/HG-U133_Plus_2/HG-U133_Plus_2.cdf");
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doCdfXdaGenoMpamBad()
{
  RT_Test* rt;
  
  //Ensure that chip type mismatch will fail so we can do Negative Test
  Fs::fileCopy("${gold_lib}/Mapping10K_Xba131/Mapping10K_Xba131.CDF","${out_dir}/foo.CDF");
  
  rt = doBadChipTypeTest("doCdfXdaGenoMpamBad","Mapping10K_Xba131.xda.genotyping.mpam.chp", "${out_dir}/foo.CDF");
  return rt;
}
