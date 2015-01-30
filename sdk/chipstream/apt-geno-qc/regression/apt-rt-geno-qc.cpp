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


RT_Test* doGenoQc500k();
RT_Test* doGenoQcSnp5();
RT_Test* doGenoQcSnp6();
// the "regression-data/data/lib/Axiom_GW_Hu_SNP/s3" directory is missing.
RT_Test* doGenoQcAxiomCdf();
RT_Test* doGenoQcAxiomSpf();

///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

  RT_Main* rMain;
  rMain = new RT_Main("apt-rt-geno-qc");

  if ( !Fs::dirExists("test-generated" )) {
    Fs::mkdir("test-generated", false);
  }

  rMain->setVar("sdk_top", "../../..");
  rMain->setVar("dmCallsSuffix", ".dm.txt");
  rMain->setVar("gold_lib_mapping250ksty", "${gold_dir}/lib/Mapping250K_Sty");
  rMain->setVar("gold_genoqc_mapping250ksty", "${gold_dir}/chipstream/geno-qc/Mapping250K_Sty");
  rMain->setVar("gold_genoqc_mapping250ksty_testname", "${gold_dir}/chipstream/geno-qc/Mapping250K_Sty/${testname}");
  rMain->setVar("gold_lib_genomewidesnp5", "${gold_dir}/lib/GenomeWideSNP_5");
  rMain->setVar("gold_genoqc_genomewidesnp5_testname", "${gold_dir}/chipstream/geno-qc/GenomeWideSNP_5/${testname}");
  rMain->setVar("gold_lib_genomewidesnp6", "${gold_dir}/lib/GenomeWideSNP_6");
  rMain->setVar("gold_genoqc_genomewidesnp6_testname", "${gold_dir}/chipstream/geno-qc/GenomeWideSNP_6/${testname}");
  rMain->setVar("gold_genoqc_genomewidesnp6", "${gold_dir}/chipstream/geno-qc/GenomeWideSNP_6");
  rMain->setVar("gold_lib_axiomgwhusnp", "${gold_dir}/lib/Axiom_GW_Hu_SNP");
  rMain->setVar("gold_genoqc_axiomgwhusnp_testname", "${gold_dir}/chipstream/geno-qc/Axiom_GW_Hu_SNP/${testname}");
  rMain->setVar("gold_genoqc_axiomgwhusnp", "${gold_dir}/chipstream/geno-qc/Axiom_GW_Hu_SNP");

  if(argc < 1)
    {
      cout<<"Error: Incorrect number of arguments given"<<endl;
      exit(1);
    }
  //They want us to run all default tests
  else if(argc == 1)
    {
      rMain->addTest(doGenoQc500k());
      rMain->addTest(doGenoQcSnp5());
      rMain->addTest(doGenoQcSnp6());
      rMain->addTest(doGenoQcAxiomCdf());
      rMain->addTest(doGenoQcAxiomSpf());
    }
  else
    {
      rMain->parseArgv(argc, argv);

      cout<<"Running Built-In Regression Tests"<<endl;

      if(rMain->getVarVal("doGenoQc500k") != "")
	{
	  rMain->addTest(doGenoQc500k());
	}
      if(rMain->getVarVal("doGenoQcSnp5") != "")
	{
	  rMain->addTest(doGenoQcSnp5());
	}
      if(rMain->getVarVal("doGenoQcSnp6") != "")
	{
	  rMain->addTest(doGenoQcSnp6());
	}
      if(rMain->getVarVal("doGenoQcAxiomCdf") != "")
	{
	  rMain->addTest(doGenoQcAxiomCdf());
	}
      if(rMain->getVarVal("doGenoQcAxiomSpf") != "")
	{
	  rMain->addTest(doGenoQcAxiomSpf());
	}
       
    }
  rMain->runTests();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doGenoQc500k()
{
  RT_Test* rt;
  rt = new RT_Test("doGenoQcTest500k");

  //std::string dmCallsSuffix = ".dm.txt";
  const char *dmCalls[] = 
    {
      "NA06985_B01_Sty_Plate1", "NA06991_B03_Sty_Plate1", "NA06993_B02_Sty_Plate1",
      "NA06994_A11_Sty_Plate1", "NA07000_A10_Sty_Plate1", "NA07019_A09_Sty_Plate1",
      "NA07022_A08_Sty_Plate1", "NA07029_A12_Sty_Plate1", "NA07034_B05_Sty_Plate1",
      "NA07048_B06_Sty_Plate1", "NA07055_B04_Sty_Plate1", "NA07056_A07_Sty_Plate1",
      "NA07345_B10_Sty_Plate1", "NA07348_B12_Sty_Plate1", "NA07357_B11_Sty_Plate1",
      "NA10846_A06_Sty_Plate1", "NA10847_A03_Sty_Plate1", "NA10851_B09_Sty_Plate1",
      "NA10854_C09_Sty_Plate1", "NA10855_C12_Sty_Plate1",
      NULL
    };
  
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_geno_qc}");
  cmd->addArg("--dm-out", "${out_testname}");
  cmd->addArg("--cdf-file", "${gold_lib_mapping250ksty}/Mapping250K_Sty.cdf");
  cmd->addArg("--qca-file", "${gold_lib_mapping250ksty}/Mapping250K_Sty.qca");
  cmd->addArg("--qcc-file", "${gold_lib_mapping250ksty}/Mapping250K_Sty.qcc");
  cmd->addArg("--cel-files", "${gold_genoqc_mapping250ksty}/cel-files.txt");
  cmd->addArg("--out-file", "${out_testname}/geno-qc-sty-test.gqc");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_mixedfile}");
  check1->addArg("--gen", "${out_testname}/geno-qc-sty-test.gqc");
  check1->addArg("--gold", "${gold_genoqc_mapping250ksty}/doGenoQcTest500k/geno-qc-sty-test.gqc");
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--skipLines", "1");
  check1->addArg("--allowedMismatch", "0");

  vector<string> dmTest = Util::addPrefixSuffix(dmCalls, "${out_testname}/", "${dmCallsSuffix}");
  vector<string> dmGold = Util::addPrefixSuffix(dmCalls, "${gold_genoqc_mapping250ksty_testname}/", "${dmCallsSuffix}");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->setGenFiles(dmTest);
  check2->setGoldFiles(dmGold);
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "1");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doGenoQcSnp5()
{
  RT_Test* rt;
  rt = new RT_Test("doGenoQcTestSnp5");

  const char *dmCalls[] = 
    {
      "NA06985_GW5_C", "NA06991_GW5_C", "NA06993_GW5_C",
      "NA06994_GW5_C", "NA07000_GW5_C", "NA07019_GW5_C",     
      "NA07022_GW5_C", "NA07029_GW5_C", "NA07034_GW5_C",
      "NA07048_GW5_C", "NA07055_GW5_C", "NA07056_GW5_C",
      "NA07345_GW5_C", "NA07348_GW5_C", "NA07357_GW5_C",
      "NA10830_GW5_C", "NA10831_GW5_C", "NA10835_GW5_C",
      "NA10838_GW5_C", "NA10839_GW5_C",
      NULL
    };

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_geno_qc}");
  cmd->addArg("--dm-out", "${out_testname}");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.cdf");
  cmd->addArg("--qca-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.qca");
  cmd->addArg("--qcc-file", "${gold_lib_genomewidesnp5}/GenomeWideSNP_5.qcc");
  cmd->addArg("--cel-files", "../../../regression-data/data/chipstream/geno-qc/GenomeWideSNP_5/cel-files.txt");
  cmd->addArg("--out-file", "${out_testname}/geno-qc-snp5-test.gqc");


  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_mixedfile}");
  check1->addArg("--gen", "${out_testname}/geno-qc-snp5-test.gqc");
  check1->addArg("--gold", "${gold_genoqc_genomewidesnp5_testname}/geno-qc-snp5-test.gqc");
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--skipLines", "1");
  check1->addArg("--allowedMismatch", "0");
  
  vector<string> dmTest = Util::addPrefixSuffix(dmCalls, "${out_testname}/", "${dmCallsSuffix}");
  vector<string> dmGold = Util::addPrefixSuffix(dmCalls, "${gold_genoqc_genomewidesnp5_testname}/", "${dmCallsSuffix}");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->setGenFiles(dmTest);
  check2->setGoldFiles(dmGold);
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "1");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doGenoQcSnp6()
{
  RT_Test* rt;
  rt = new RT_Test("doGenoQcTestSnp6");

  const char *dmCalls[] = 
    {
      "/NA06985_GW6_C", "/NA06993_GW6_C", "/NA07000_GW6_C",
      "/NA07019_GW6_C", "/NA07029_GW6_C", "/NA07056_GW6_C",   
      "/NA07345_GW6_C", "/NA10830_GW6_C", "/NA10831_GW6_C",      
      "/NA10839_GW6_C", "/NA10855_GW6_C", "/NA10857_GW6_C",
      "/NA10860_GW6_C", "/NA10863_GW6_C", "/NA11881_GW6_C",
      "/NA11882_GW6_C", "/NA11993_GW6_C", "/NA12003_GW6_C",
      "/NA12004_GW6_C",

      NULL
    };
    
  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_geno_qc}");
  cmd->addArg("--dm-out", "${out_testname}");
  cmd->addArg("--cdf-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.cdf");
  cmd->addArg("--qca-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.r2.qca");
  cmd->addArg("--qcc-file", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.r2.qcc");
  cmd->addArg("--chrX-probes", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.chrXprobes");
  cmd->addArg("--chrY-probes", "${gold_lib_genomewidesnp6}/GenomeWideSNP_6.chrYprobes");
  cmd->addArg("--cel-files", "${gold_genoqc_genomewidesnp6}/fas_cel_files.txt");
  cmd->addArg("--out-file", "${out_testname}/geno-qc-snp6-test.gqc");
  
  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_mixedfile}");
  check1->addArg("--gen", "${out_testname}/geno-qc-snp6-test.gqc");
  check1->addArg("--gold", "${gold_genoqc_genomewidesnp6_testname}/geno-qc-snp6-test.gqc");
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--skipLines", "1");
  check1->addArg("--allowedMismatch", "0");
  
  vector<string> dmTest = Util::addPrefixSuffix(dmCalls, "${out_testname}/", "${dmCallsSuffix}");
  vector<string> dmGold = Util::addPrefixSuffix(dmCalls, "${gold_genoqc_genomewidesnp6_testname}/", "${dmCallsSuffix}");

  RT_Check* check2 = rt->newCheck();
  check2->setExe("${apt_check_mixedfile}");
  check2->setGenFiles(dmTest);
  check2->setGoldFiles(dmGold);
  check2->addArg("--epsilon", "0.0001");
  check2->addArg("--skipLines", "1");
  check2->addArg("--allowedMismatch", "0");

  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doGenoQcAxiomCdf()
{
  RT_Test* rt;
  rt = new RT_Test("doGenoQcAxiomCdf");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_geno_qc}");
  cmd->addArg("--cdf-file", "${gold_lib_axiomgwhusnp}/Axiom_GW_Hu_SNP.s3.cdf");
  cmd->addArg("--qca-file", "${gold_lib_axiomgwhusnp}/Axiom_GW_Hu_SNP.s3.qca");
  cmd->addArg("--qcc-file", "${gold_lib_axiomgwhusnp}/Axiom_GW_Hu_SNP.s3.qcc");
  cmd->addArg("--cel-files", "${gold_genoqc_axiomgwhusnp}/celfiles.txt");
  cmd->addArg("--out-file", "${out_testname}/doGenoQcAxiomCdf.gqc");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_mixedfile}");
  check1->addArg("--gen", "${out_testname}/doGenoQcAxiomCdf.gqc");
  check1->addArg("--gold", "${gold_genoqc_axiomgwhusnp_testname}/doGenoQcAxiomCdf.gqc");
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--skipLines", "1");
  check1->addArg("--allowedMismatch", "0");
  
  return rt;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
RT_Test* doGenoQcAxiomSpf()
{
  RT_Test* rt;
  rt = new RT_Test("doGenoQcAxiomSpf");

  RT_Cmd* cmd = rt->newCmd();
  cmd->setExe("${apt_geno_qc}");
  cmd->addArg("--spf-file", "${gold_lib_axiomgwhusnp}/Axiom_GW_Hu_SNP.s3.spf");
  cmd->addArg("--qca-file", "${gold_lib_axiomgwhusnp}/Axiom_GW_Hu_SNP.s3.qca");
  cmd->addArg("--qcc-file", "${gold_lib_axiomgwhusnp}/Axiom_GW_Hu_SNP.s3.qcc");
  cmd->addArg("--cel-files", "${gold_genoqc_axiomgwhusnp}/celfiles.txt");
  cmd->addArg("--out-file", "${out_testname}/doGenoQcAxiomSpf.gqc");

  RT_Check* check1 = rt->newCheck();
  check1->setExe("${apt_check_mixedfile}");
  check1->addArg("--gen", "${out_testname}/doGenoQcAxiomSpf.gqc");
  check1->addArg("--gold", "${gold_genoqc_axiomgwhusnp_testname}/doGenoQcAxiomSpf.gqc");
  check1->addArg("--epsilon", "0.0001");
  check1->addArg("--skipLines", "1");
  check1->addArg("--allowedMismatch", "0");
  
  return rt;
}
