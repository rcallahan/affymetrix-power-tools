////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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

RT_Test* doTranscriptClusterTest();
RT_Test* doProbesetIdTest();

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

  RT_Main* rMain;
  rMain = new RT_Main("apt-rt-tsv-join");

  if ( !Fs::dirExists("tsv-join-test-generated") ) {
    Fs::mkdir("tsv-join-test-generated", false);
  }

  rMain->setVar("sdk_top", "../../..");
  rMain->setVar("gold_idata_lib_huex", "${gold_dir}/idata/lib/HuEx-1_0-st-v2");
  rMain->setVar("gold_idata_tsvjoin", "${gold_dir}/idata/tsv-join");
  rMain->setVar("tsv_out_dir", "tsv-join-test-generated");

  if(argc < 1)
    {
      cout<<"Error: Incorrect number of arguments given"<<endl;
      exit(1);
    }
  //They want us to run all default tests
  else if(argc == 1)
    {
      //RUN ALL TESTS
      cout<<"Running all default Regression Tests!"<<endl;
      rMain->addTest(doTranscriptClusterTest());
      rMain->addTest(doProbesetIdTest());
    }
  else
    {
      rMain->parseArgv(argc, argv);

      cout<<"Running Built-In Regression Tests"<<endl;
	 
      if(rMain->getVarVal("doTranscriptClusterTest") != "")
	{
	  rMain->addTest(doTranscriptClusterTest());
	}
      if(rMain->getVarVal("doProbesetIdTest") != "")
	{
	  rMain->addTest(doProbesetIdTest());
	}
    }

  rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doTranscriptClusterTest()
{
  RT_Test* rt;
  rt = new RT_Test("doTranscriptClusterTest");
 
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_tsv_join}");
  cmd->addArg("-k", "transcript_cluster_id");
  cmd->addArg("-o", "${tsv_out_dir}/join-transcript-cluster.txt");
  cmd->addArg("${gold_idata_lib_huex}/HuEx-1_0-st-probeset-annot.85000.csv");
  cmd->addArg("${gold_idata_lib_huex}/HuEx-1_0-st-transcript-annot.16000.csv");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_textfile}");
  check->addArg("--gen", "${tsv_out_dir}/join-transcript-cluster.txt");
  check->addArg("--gold", "${gold_idata_tsvjoin}/join-transcript-cluster.txt");
  check->addArg("--skipLines", "34");
 
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doProbesetIdTest()
{
  RT_Test* rt;
  rt = new RT_Test("doProbesetIdTest");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_tsv_join}");
  cmd->addArg("-k", "probeset_id");
  cmd->addArg("-o", "${tsv_out_dir}/join-probeset-id.txt");
  cmd->addArg("${gold_idata_lib_huex}/HuEx-1_0-st-transcript-annot.csv");
  cmd->addArg("${gold_idata_tsvjoin}/quant-norm.pm-gcbg.iter-plier.summary.txt");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_textfile}");
  check->addArg("--gen", "${tsv_out_dir}/join-probeset-id.txt");
  check->addArg("--gold", "${gold_idata_tsvjoin}/join-probeset-id.txt");
  check->addArg("--skipLines", "34");
  
  return rt;
}
