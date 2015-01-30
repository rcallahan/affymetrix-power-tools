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

RT_Test* testNumericDifferences();
RT_Test* testMixedDifferences();

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
 RT_Main* rMain;
 rMain = new RT_Main("apt-rt-matrix-diff");

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
     rMain->addTest(testNumericDifferences());
     rMain->addTest(testMixedDifferences());
   }
 else
   {
     rMain->parseArgv(argc, argv);

     cout<<"Running Built-In Regression Tests"<<endl;

     if(rMain->getVarVal("testNumericDifferences") != "")
       {
	 rMain->addTest(testNumericDifferences());
       }
     if(rMain->getVarVal("testMixedDifferences") != "")
       {
	 rMain->addTest(testMixedDifferences());
       }
   }
 rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* testNumericDifferences()
{
  RT_Test* rt;
  rt = new RT_Test("testMixedDifferences");
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_matrix_diff}");
  cmd->addArg("-l", "1");
  cmd->addArg("-e", "0.1");
  cmd->addArg("data/test-numeric-gold.txt data/test-numeric-generated.txt");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* testMixedDifferences()
{
  RT_Test* rt;
  rt = new RT_Test("testMixedDifferences");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_matrix_diff}");
  cmd->addArg("--mixed");
  cmd->addArg("-l", "2");
  cmd->addArg("-p", "data/test-mixed-gold.txt data/test-mixed-generated.txt");

  return rt;
}
