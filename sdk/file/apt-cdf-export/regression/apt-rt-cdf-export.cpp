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

RT_Test* exportCdf();

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

 RT_Main* rMain;
 rMain = new RT_Main("apt-rt-cdf-export");

 if ( !Fs::dirExists("test-generated") ) {
   Fs::mkdir("test-generated", false);
 }

 rMain->setVar("sdk_top", "../../..");

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
     rMain->addTest(exportCdf());
   }
 else
   {
     rMain->parseArgv(argc, argv);

     cout<<"Running Built-In Regression Tests"<<endl;
	 
     if(rMain->getVarVal("exportCdf") != "")
       {
	 rMain->addTest(exportCdf());
       }
   }

 rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* exportCdf()
{
  RT_Test* rt;
  rt = new RT_Test("exportCdf");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cdf_export}");;
  cmd->addArg("-c", "${sdk_top}/rawq/test/data/Test3.CDF");
  cmd->addArg("-ps", "Pae_16SrRNA_s_at");
  cmd->addArg("-ps", "AFFX-Athal-Actin_5_r_at > ${out_dir}/export.Test3.txt");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_textfile}");
  check->addArg("--gen", "${out_dir}/export.Test3.txt");
  check->addArg("--gold", "data/export.Test3.txt");
  check->addArg("--skipLines", "0");

  return rt;
}
