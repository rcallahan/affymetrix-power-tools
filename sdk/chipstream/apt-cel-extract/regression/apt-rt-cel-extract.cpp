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


RT_Test* doExtractIntensities();
RT_Test* doExtractIntensitiesSubset();
RT_Test* doExtractIntensitiesRaw();
RT_Test* doExtractMultiChannelIntensitiesRaw();
RT_Test* doExtractMultiChannelIntensitiesSubset();

//std::string tissueCels = "${gold_cel_hgu133plus2}/pancrease-rep1.cel ${gold_cel_hgu133plus2}/pancrease-rep2.cel ${gold_cel_hgu133plus2}/pancrease-rep3.cel";

//std::string axiomCels = " ${gold_cel_axiomgwhusnp}/NA18526.CEL ${gold_cel_axiomgwhusnp}/NA18856.CEL ${gold_cel_axiomgwhusnp}/NA18862.CEL";

///@todo it would be good to add a lot more tests. w/ analysis string specifically.
///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
 RT_Main* rMain;
 rMain = new RT_Main("apt-rt-cel-extract");

 if ( !Fs::dirExists("test-generated") ) {
   Fs::mkdir("test-generated", false);
 }

 rMain->setVar("sdk_top", "../../..");
 rMain->setVar("gold_lib_hgu133plus2", "${gold_dir}/lib/HG-U133_Plus_2");
 rMain->setVar("gold_celextract_hgu133plus2", "${gold_dir}/chipstream/cel-extract/HG-U133_Plus_2");
 rMain->setVar("gold_celextract_axiomgwhusnp", "${gold_dir}/chipstream/cel-extract/Axiom_GW_Hu_SNP");
 rMain->setVar("gold_cel_axiomgwhusnp", "${gold_dir}/cel/Axiom_GW_Hu_SNP");
 rMain->setVar("gold_cel_hgu133plus2", "${gold_dir}/cel/HG-U133_Plus_2");
 rMain->setVar("gold_lib_axiomgwhusnp", "${gold_dir}/lib/Axiom_GW_Hu_SNP");
 rMain->setVar("tissuecel1", "${gold_cel_hgu133plus2}/pancrease-rep1.cel");
 rMain->setVar("tissuecel2", "${gold_cel_hgu133plus2}/pancrease-rep2.cel");
 rMain->setVar("tissuecel3", "${gold_cel_hgu133plus2}/pancrease-rep3.cel");
 rMain->setVar("axiomcel1", "${gold_cel_axiomgwhusnp}/NA18526.CEL");
 rMain->setVar("axiomcel2", "${gold_cel_axiomgwhusnp}/NA18856.CEL");
 rMain->setVar("axiomcel3", "${gold_cel_axiomgwhusnp}/NA18862.CEL");

 if(argc < 1)
   {
     cout<<"Error: Incorrect number of arguments given"<<endl;
     exit(1);
   }
 //They want us to run all default tests
 else if(argc == 1)
   {
     rMain->addTest(doExtractIntensities());
     rMain->addTest(doExtractIntensitiesSubset());
     rMain->addTest(doExtractIntensitiesRaw());
     rMain->addTest(doExtractMultiChannelIntensitiesRaw());
     rMain->addTest(doExtractMultiChannelIntensitiesSubset());
   }
 else
   {
     rMain->parseArgv(argc, argv);

     cout<<"Running Built-In Regression Tests"<<endl;
     if(rMain->getVarVal("doExtractIntensities") != "")
       {
	 rMain->addTest(doExtractIntensities());
       }
     if(rMain->getVarVal("doExtractIntensitiesSubset") != "")
       {
	 rMain->addTest(doExtractIntensitiesSubset());
       }
     if(rMain->getVarVal("doExtractIntensitiesRaw") != "")
       {
	 rMain->addTest(doExtractIntensitiesRaw());
       }
     if(rMain->getVarVal("doExtractMultiChannelIntensitiesRaw") != "")
       {
	 rMain->addTest(doExtractMultiChannelIntensitiesRaw());
       }
     if(rMain->getVarVal("doExtractMultiChannelIntensitiesSubset") != "")
       {
	 rMain->addTest(doExtractMultiChannelIntensitiesSubset());
       }
   }
 rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doExtractIntensities()
{
  RT_Test* rt;
  rt = new RT_Test("doExtractIntensities");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_extract}");
  cmd->addArg("-o", "${out_dir}/extract.intensities.txt");
  cmd->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd->addArg("${tissuecel1}");
  cmd->addArg("${tissuecel2}");
  cmd->addArg("${tissuecel3}");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_mixedfile}");
  check->addArg("--gen", "${out_dir}/extract.intensities.txt");
  check->addArg("--gold", "${gold_celextract_hgu133plus2}/extract.intensities.txt");
  check->addArg("--epsilon", "0.1");
  check->addArg("--skipLines", "0");
  check->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doExtractIntensitiesSubset()
{
  RT_Test* rt;
  rt = new RT_Test("doExtractIntensitiesSubset");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_extract}");
  cmd->addArg("--probeset-ids=${gold_celextract_hgu133plus2}/ps.txt");
  cmd->addArg("-o", "${out_dir}/extract.intensities.sub.txt");
  cmd->addArg("-d", "${gold_lib_hgu133plus2}/HG-U133_Plus_2.cdf");
  cmd->addArg("${tissuecel1}");
  cmd->addArg("${tissuecel2}");
  cmd->addArg("${tissuecel3}");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_mixedfile}");
  check->addArg("--gen", "${out_dir}/extract.intensities.sub.txt");
  check->addArg("--gold", "${gold_celextract_hgu133plus2}/extract.intensities.sub.txt");
  check->addArg("--epsilon", "0.1");
  check->addArg("--skipLines", "0");
  check->addArg("--allowedMismatch", "0");
  
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doExtractIntensitiesRaw()
{
  RT_Test* rt;
  rt = new RT_Test("doExtractIntensitiesRaw");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_extract}");
  cmd->addArg("-o", "${out_dir}/extract.intensities.raw.txt");
  cmd->addArg("${tissuecel1}");
  cmd->addArg("${tissuecel2}");
  cmd->addArg("${tissuecel3}");

  RT_Check* check = rt->newCheck();
  check->setExe("${apt_check_mixedfile}");
  check->addArg("--gen", "${out_dir}/extract.intensities.raw.txt");
  check->addArg("--gold", "${gold_celextract_hgu133plus2}/extract.intensities.raw.txt");
  check->addArg("--epsilon", "0.1");
  check->addArg("--skipLines", "0");
  check->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doExtractMultiChannelIntensitiesRaw()
{
  RT_Test* rt;
  rt = new RT_Test("doExtractMultiChannelIntensitiesRaw");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_extract}");
  cmd->addArg("-o", "${out_dir}/extract.mc.intensities.raw.txt");
  cmd->addArg("${axiomcel1}");
  cmd->addArg("${axiomcel2}");
  cmd->addArg("${axiomcel3}");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_mixedfile}");
  check1->addArg("--gen", "${out_dir}/extract.mc.intensities.raw.txt.channel0");
  check1->addArg("--gold", "${gold_celextract_axiomgwhusnp}/extract.mc.intensities.raw.txt.channel0");
  check1->addArg("--epsilon", "0.1");
  check1->addArg("--skipLines", "0");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_dir}/extract.mc.intensities.raw.txt.channel1");
  check2->addArg("--gold", "${gold_celextract_axiomgwhusnp}/extract.mc.intensities.raw.txt.channel1");
  check2->addArg("--epsilon", "0.1");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doExtractMultiChannelIntensitiesSubset()
{
  RT_Test* rt;
  rt = new RT_Test("doExtractMultiChannelIntensitiesSubset");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_cel_extract}");
  cmd->addArg("--probeset-ids=../../../regression-data/data/chipstream/cel-extract/Axiom_GW_Hu_SNP/ps.txt");
  cmd->addArg("-o", "${out_dir}/extract.mc.intensities.sub.txt");
  cmd->addArg("-d", "${gold_lib_axiomgwhusnp}/Axiom_GW_Hu_SNP.s3.cdf");
  cmd->addArg("${axiomcel1}");
  cmd->addArg("${axiomcel2}");
  cmd->addArg("${axiomcel3}");


  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_mixedfile}");
  check1->addArg("--gen", "${out_dir}/extract.mc.intensities.sub.txt.channel0");
  check1->addArg("--gold", "${gold_celextract_axiomgwhusnp}/extract.mc.intensities.sub.txt.channel0");
  check1->addArg("--epsilon", "0.1");
  check1->addArg("--skipLines", "0");
  check1->addArg("--allowedMismatch", "0");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->addArg("--gen", "${out_dir}/extract.mc.intensities.sub.txt.channel1");
  check2->addArg("--gold", "${gold_celextract_axiomgwhusnp}/extract.mc.intensities.sub.txt.channel1");
  check2->addArg("--epsilon", "0.1");
  check2->addArg("--skipLines", "0");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}
