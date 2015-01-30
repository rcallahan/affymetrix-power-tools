////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

/**
 * @file   ChpCheckTest.cpp
 * @author vliber
 * last change by rsatin on 10/13/09
 *
 * @brief  Testing the ChipCheck functions.
 *
 */

#include "util/AffxString.h"
#include "util/ChpCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
//
#include "util/CPPTest/Setup.h"

using namespace std;

class ChpCheckTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ChpCheckTest );
  CPPUNIT_TEST( testHeaderMap );
  CPPUNIT_TEST( testDataMap );
  CPPUNIT_TEST( testCheckMap );
  CPPUNIT_TEST( testHeaderExp );   
  CPPUNIT_TEST( testDataExp );
  CPPUNIT_TEST( testCheckExp );
  CPPUNIT_TEST( testCheckFracDiffPValue );
  CPPUNIT_TEST( testCheckFracDiffSignal );
  CPPUNIT_TEST( testCheckFracDiffConf );
  CPPUNIT_TEST_SUITE_END();

private:
  /** utility functions */
  int message_error_count(string message, size_t fpos=0);

public: 
  virtual void tearDown();
  /** unit test functions */
  void testHeaderMap();
  void testDataMap();
  void testCheckMap(); 
  void testHeaderExp();
  void testDataExp();
  void testCheckExp();
  void testCheckFracDiffPValue();
  void testCheckFracDiffSignal();
  void testCheckFracDiffConf();
};

CPPUNIT_TEST_SUITE_REGISTRATION( ChpCheckTest );

void ChpCheckTest::tearDown () 
{
  // clean up after failed test cases (pops any dead objects from stack)
  while( Verbose::getParam().m_MsgHandler.size() > 1 )
    Verbose::popMsgHandler();
}

int ChpCheckTest::message_error_count( string message, size_t fpos ) {
  long error_count = 0;
  // parse difference count from output in log file of form:
  // Max diff: 0.00197983 is greater than expected (0.0001) [confidence]
  // 19 of 150 (12.6667%) were different.
  size_t fpos_last = message.find(" were different",fpos);
  if( fpos_last != string::npos ) {
    size_t fpos_start = message.find_last_of("\n", fpos_last);
    if( fpos_last != string::npos ) {
      char *endptr;
      int radix = 10;
      std::string strCount = message.substr( fpos_start+1, fpos_last-fpos_start );
      error_count = strtol( strCount.c_str(), &endptr, radix ); 
	}
  }
  return (int)error_count;
}

void ChpCheckTest::testHeaderMap()
{
  cout<<endl;
  cout<<endl;
  Verbose::out(1, "***ChpCheck testcases***");
  Verbose::out(1, "ChpCheckTest::testHeaderMap");
  //positive 
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmmGold_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_small.chp");
  ChpCheck chpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(chpCheckTest1.check(errorMsg));
  POSITIVE_TEST(chpCheckTest1.check(errorMsg));

  //negative test cases

  //Goal is to receive a message: Error encountered: Error: cols not the same.
  generated.clear();
  gold.clear();
  gold.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmmGold_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad2_small.chp");
  ChpCheck chpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest5.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: rows not the same.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad3_small.chp");
  ChpCheck chpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest6.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: ChipType not the same.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad10_small.chp");
  ChpCheck chpCheckTest8 (generated, gold);
  CPPUNIT_ASSERT(! chpCheckTest8.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: AlgVersion not the same.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad9_small.chp");
  ChpCheck chpCheckTest9 (generated, gold);
  CPPUNIT_ASSERT(! chpCheckTest9.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: ProgID not the same.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad6_small.chp");
  ChpCheck chpCheckTest10 (generated, gold);
  CPPUNIT_ASSERT(! chpCheckTest10.check(errorMsg));

  //first algorithm parameter
  //Goal is to receive a message: 
  //Error: Test missing field: 'program-name' Error: Test missing field: 'program-name'
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad11_small.chp");
  ChpCheck chpCheckTest12 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest12.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest12.check(errorMsg))
    cout << errorMsg << endl;
#endif
  //Error: for field 'program-name' expecting: 'apt-probeset-genotype' got: 'ast-probeset-genotype'
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad12_small.chp");
  ChpCheck chpCheckTest13 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest13.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest13.check(errorMsg))
    cout << errorMsg << endl;
#endif
  //last algorithm parameter
  //Error: Test missing field: 'analysis-text' Error: Test missing field: 'analysis-text' 
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad13_small.chp");
  ChpCheck chpCheckTest14 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest14.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest14.check(errorMsg))
    cout << errorMsg << endl;
#endif
  //Error: for field 'analysis-text' expecting: 'quant-norm.sketch=50000,pm-only,brlmm.transform=ccs.K=4' got: 'quant-norm.sketch=50000,pm-only,brlmm.transform=cct.K=4' 
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad14_small.chp");
  ChpCheck chpCheckTest15 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest15.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest15.check(errorMsg))
    cout << errorMsg << endl;
#endif

 //summary statistics
 //no summary in the file

}

void ChpCheckTest::testDataMap()
{
  cout<<endl;
  Verbose::out(1, "ChpCheckTest::testDataMap");
  //positive 
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmmGold_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_small.chp");
  ChpCheck chpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(chpCheckTest1.check(errorMsg)); 
  POSITIVE_TEST(chpCheckTest1.check(errorMsg));

  //m_Prefix=""
  ChpCheck chpCheckTest2 (generated, gold, 0, string(""));
  CPPUNIT_ASSERT(chpCheckTest2.check(errorMsg));


  //negative test cases

  //first probe set data change
  //gold 536.5 gen 536.563 eps 0.0001 vliber output
  //Max diff: 0.0625 is greater than expected (0.0001) [confidence]
  //1 of 262338 (0.000381188%) were different.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad15.chp: checked 262338 genotype entries.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad15.chp: max confidence diff is: 0.0625
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad15.chp: max pvalue diff is: 0
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad15.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad15_small.chp");
  ChpCheck chpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest6.check(errorMsg));

  
  //25-th probe set data change

  //gold 0.0118848 gen 0.0128613 eps 0.0001 vliber output
  //Max diff: 0.000976563 is greater than expected (0.0001) [confidence]
  //1 of 262338 (0.000381188%) were different.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad16.chp: checked 262338 genotype entries.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad16.chp: max confidence diff is: 0.000976563
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad16.chp: max pvalue diff is: 0
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad16.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad16_small.chp");
  ChpCheck chpCheckTest7 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest7.check(errorMsg));
  
  //gold 714.957 gen 778.957 eps 0.0001 vliber output
  //Max diff: 64 is greater than expected (0.0001) [pvalue]
  //1 of 262338 (0.000381188%) were different.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad17.chp: checked 262338 genotype entries.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad17.chp: max confidence diff is: 0
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad17.chp: max pvalue diff is: 64
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad17.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad17_small.chp");
  ChpCheck chpCheckTest8 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest8.check(errorMsg));
  

  //gold 61.059 gen 66.1181 eps 0.0001 vliber output
  //Max diff: 5.05903 is greater than expected (0.0001) [pvalue]
  //1 of 262338 (0.000381188%) were different.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad18.chp: checked 262338 genotype entries.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad18.chp: max confidence diff is: 0
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad18.chp: max pvalue diff is: 5.05903
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad18.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad18_small.chp");
  ChpCheck chpCheckTest9 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest9.check(errorMsg));
  
  //gold 0.725673 gen 0.788173 eps 0.0001 vliber output
  //Max diff: 0.0625 is greater than expected (0.0001) [pvalue]
  //1 of 262338 (0.000381188%) were different.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad19.chp: checked 262338 genotype entries.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad19.chp: max confidence diff is: 0
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad19.chp: max pvalue diff is: 0.0625
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad19.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad19_small.chp");
  ChpCheck chpCheckTest10 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest10.check(errorMsg));
  
  //gold 3.40282e+038 gen 3.40282e+038 eps 0.0001 vliber output
  //Max diff: 3.24519e+032 is greater than expected (0.0001) [pvalue]
  //1 of 262338 (0.000381188%) were different.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad20.chp: checked 262338 genotype entries.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad20.chp: max confidence diff is: 0
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad20.chp: max pvalue diff is: 3.24519e+032
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad20.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad20_small.chp");
  ChpCheck chpCheckTest11 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest11.check(errorMsg));

  //change confidence and 1 pvalue only
  //gold 0.0118848 gen 0.012869 eps 0.0001
  //gold 714.957 gen 714.962 eps 0.0001
  //Max diff: 0.000984192 is greater than expected (0.0001) [confidence]
  //Max diff: 0.00488281 is greater than expected (0.0001) [pvalue]
  //1 of 262338 (0.000381188%) were different.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad21.chp: checked 262338 genotype entries.
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad21.chp: max confidence diff is: 0.000984192
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad21.chp: max pvalue diff is: 0.00488281
  //input\chp\map\CCL-256.1D_NSP.brlmm_bad21.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad21_small.chp");
  ChpCheck chpCheckTest12 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest12.check(errorMsg));

  //positive m_DiffAllowed == error == 1
  ChpCheck chpCheckTest13 (generated, gold, 1);
  CPPUNIT_ASSERT(chpCheckTest13.check(errorMsg));
  //positive eps=0.0049
  ChpCheck chpCheckTest14 (generated, gold, 0,"apt-", 0.0049);
  CPPUNIT_ASSERT(chpCheckTest14.check(errorMsg));
}

void ChpCheckTest::testCheckMap()
{
  cout<<endl;
  Verbose::out(1, "ChpCheckTest::testCheckMap");
  //positive: binary
  vector<string> generated, gold;
  gold.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmmGold_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_small.chp");
  std::string errorMsg ("");
  ChpCheck chpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(chpCheckTest1.check(errorMsg));
  POSITIVE_TEST(chpCheckTest1.check(errorMsg));
      
  //negative: 
  //Goal is to receive a message: CelCheck::check() - generated and gold vectors must be same size.
  gold.clear();
  generated.clear();
  gold.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmmGold_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256D_NSP.brlmm_small.chp");
  ChpCheck chpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest2.check (errorMsg));
  
  //bad path to the gold file
  //Goal is to receive a message: FATAL ERROR: Failed to get entry for gold (0). File: input1\chp\map\CCL-256.1D_NSP.brlmmGold.chp Error: 
  gold.clear();
  generated.clear();
  gold.push_back("input1/chp/map/CCL-256.1D_NSP.brlmmGold_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm.chp");
  ChpCheck chpCheckTest3 (generated, gold);
  Verbose::setLevel(-1);
  CPPUNIT_ASSERT(!chpCheckTest3.check(errorMsg));
  Verbose::setLevel(3);

  
  //bad path to the generated file
  //Goal is to receive a message: FATAL ERROR: ChpCheck::ChpCheck() - unknown CHP type.
  errorMsg="";
  gold.clear();
  generated.clear();
  gold.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmmGold.chp");
  generated.push_back("input1/chp/map/CCL-256.1D_NSP.brlmm_small.chp");
  ChpCheck chpCheckTest4 (generated, gold);
  Verbose::setLevel(-1);
  CPPUNIT_ASSERT(!chpCheckTest4.check(errorMsg));
  Verbose::setLevel(3);
  
}


void ChpCheckTest::testHeaderExp()
{
  
  cout<<endl;
  Verbose::out(1, "ChpCheckTest::testHeaderExp");
  vector<string> generated, gold;
  std::string errorMsg ("");
 
  //positive 
  generated.clear();
  gold.clear();
  gold.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1Gold_small.CHP");
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP");
  ChpCheck chpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(chpCheckTest1.check(errorMsg));
  POSITIVE_TEST(chpCheckTest1.check(errorMsg));

  //negative test cases
     
  //Goal is to receive a message: Error encountered: Error: cols not the same.
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad2_small.CHP");
  ChpCheck chpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest5.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: rows not the same.
  errorMsg.clear();
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad3_small.CHP");
  ChpCheck chpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest6.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: NumProbeSets not the same.
  errorMsg.clear();
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad4_small.CHP");
  ChpCheck chpCheckTest7 (generated, gold); 
  CPPUNIT_ASSERT(!chpCheckTest7.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: NumProbeSets not the same.
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad4_small_exception.CHP");
  ChpCheck chpCheckTest7a (generated, gold); 
  CPPUNIT_ASSERT(!chpCheckTest7a.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: Assay Type not the same.
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad5_small.CHP");
  ChpCheck chpCheckTest11 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest11.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: ChipType not the same.
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad7_small.CHP");
  ChpCheck chpCheckTest8 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest8.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: AlgVersion not the same. 
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad8_small.CHP");
  ChpCheck chpCheckTest9 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest9.check(errorMsg));

  //Goal is to receive a message: Error encountered: Error: ProgID not the same
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad6_small.CHP");
  ChpCheck chpCheckTest10 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest10.check(errorMsg));
  
  //first algorithm parameter
  //Goal is to receive a message: Error: Test missing field: 'HZ'
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad9_small.CHP");
  ChpCheck chpCheckTest12 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest12.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest12.check(errorMsg))
    cout << errorMsg << endl;
#endif
  //Error: for field 'HZ' expecting: '4' got: 'D' 
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad10_small.CHP");
  ChpCheck chpCheckTest13 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest13.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest13.check(errorMsg))
    cout << errorMsg << endl;
#endif

  //last algorithm parameter
  //Error: Test missing field: 'SFGene' 
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad11_small.CHP");
  ChpCheck chpCheckTest14 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest14.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest14.check(errorMsg))
    cout << errorMsg << endl;
#endif
  //Error: for field 'SFGene' expecting: 'All' got: 'A|l' 
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad12_small.CHP");
  ChpCheck chpCheckTest15 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest15.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest15.check(errorMsg))
    cout << errorMsg << endl;
#endif

  //first summary statistics
  //Error: Test missing field: 'Background' 
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad13_small.CHP");
  ChpCheck chpCheckTest16 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest16.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest16.check(errorMsg))
    cout << errorMsg << endl;
#endif
 //Error: for field 'Background' expecting: 'Avg:47.76,Stdev:1.00,Max:50.7,Min:45.0 ' got: 'AvgJ47.76,Stdev:1.00,Max:50.7,Min:45.0 '
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad14_small.CHP");
  ChpCheck chpCheckTest17 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest17.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest17.check(errorMsg))
    cout << errorMsg << endl;
#endif

  //last summary statistics
  //Error: Test missing field: 'Central-' 
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad15_small.CHP");
  ChpCheck chpCheckTest18 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest18.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest18.check(errorMsg))
    cout << errorMsg << endl;
#endif
  //Error: for field 'Central-' expecting: 'Avg:21536,Count:9' got: 'Avg:21636,Count:9'
  errorMsg = ("");
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad16_small.CHP");
  ChpCheck chpCheckTest19 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest19.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if (! chpCheckTest19.check(errorMsg))
    cout << errorMsg << endl;
#endif

}

void ChpCheckTest::testDataExp()
{
  cout<<endl;
  Verbose::out(1, "ChpCheckTest::testDataExp");
  //positive 
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP");
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1Gold_small.CHP");
  ChpCheck chpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(chpCheckTest1.check(errorMsg));
  POSITIVE_TEST(chpCheckTest1.check(errorMsg));

  //m_Prefix=""
  ChpCheck chpCheckTest2 (generated, gold, 0, string(""));
  CPPUNIT_ASSERT(chpCheckTest2.check(errorMsg));

  //negative
  //signal is different for first probe set
  //gold 62.8927 gen 62.8929//vliber output
  //Max diff: 0.00012207 is greater than expected (0.0001) [signal]
  //1 of 54675 (0.00182899%) were different.
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad19.CHP: checked 54675 expression entries.
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad19.CHP: max signal diff is: 0.00012207
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad19.CHP: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad19_small.CHP");
  ChpCheck chpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest5.check(errorMsg));
  
  //positive m_DiffAllowed == error count==1
  ChpCheck chpCheckTest6 (generated, gold, 1);
  CPPUNIT_ASSERT(chpCheckTest6.check(errorMsg));
  
  //positive eps=0.002
  //gold 62.8927 gen 62.8929 eps 0.002//vliber output
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad19.CHP: checked 54675 expression entries.
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad19.CHP: max signal diff is: 0.00012207
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad19.CHP: chip is same.
  ChpCheck chpCheckTest7 (generated, gold,0,"apt-",0.002);
  CPPUNIT_ASSERT(chpCheckTest7.check(errorMsg));


  //signal is different for first  and second set
  //gold 62.8927 gen 62.8929 eps 0.0001 vliber output
  //gold 79.5101 gen 79.514 eps 0.0001 vliber output
  //Max diff: 0.00390625 is greater than expected (0.0001) [signal]
  //2 of 54675 (0.00365798%) were different.
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad20.CHP: checked 54675 expression entries.
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad20.CHP: max signal diff is: 0.00390625
  //input\chp\exp\100907_CTGF2_U133plus2_t1bad20.CHP: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad20_small.CHP");
  ChpCheck chpCheckTest8 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest8.check(errorMsg));
  
  //positive m_DiffAllowed == error count==2
  ChpCheck chpCheckTest9 (generated, gold, 2);
  CPPUNIT_ASSERT(chpCheckTest9.check(errorMsg));

  //positive eps=0.004
  ChpCheck chpCheckTest10 (generated, gold,0,"apt-",0.004);
  CPPUNIT_ASSERT(chpCheckTest10.check(errorMsg));
	 
}

void ChpCheckTest::testCheckExp()
{
  cout<<endl;
  Verbose::out(1, "ChpCheckTest::testCheckExp");
  //positive: binary
  vector<string> generated, gold;
  gold.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1Gold_small.CHP");
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP");
  std::string errorMsg ("");
  ChpCheck chpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(chpCheckTest1.check(errorMsg));
  POSITIVE_TEST(chpCheckTest1.check(errorMsg));
    
  //negative: 
  //Goal is to receive a message: CelCheck::check() - generated and gold vectors must be same size.
  gold.clear();
  generated.clear();
  gold.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1Gold_small.CHP");
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP");
  generated.push_back(INPUT+"/chp/map/100907_CTGF3_U133plus2_t1_small.CHP");
  ChpCheck chpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!chpCheckTest2.check (errorMsg));
  
  //bad path to the generated file
  //Goal is to receive a message: FATAL ERROR: ChpCheck::ChpCheck() - unknown CHP type.
  gold.clear();
  generated.clear();
  gold.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1Gold_small.CHP");
  generated.push_back("input1/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP");
  ChpCheck chpCheckTest3 (generated, gold);
  Verbose::setLevel(-1);
  CPPUNIT_ASSERT(!chpCheckTest3.check(errorMsg));
  Verbose::setLevel(3);

  //bad path to the gold file 
  //Exception has been thrown to receive an exception uncomment last statement
  //todo vliber run time exception if generated has probe set count > than gold line 289 of ChpCheck.h
  gold.clear();
  generated.clear();
  gold.push_back("input1/chp/exp/100907_CTGF2_U133plus2_t1Gold_small.CHP");
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP");
  ChpCheck chpCheckTest4 (generated, gold);
  //CPPUNIT_ASSERT(!chpCheckTest4.check(errorMsg),Except);
}

void ChpCheckTest::testCheckFracDiffPValue()
{
  // message output streamed to string object
  ostringstream MessageHandler_oss;
  MsgStream msgHandler(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler);

  cout<<endl;
  Verbose::out(1, "ChpCheckTest::testCheckFracDiffPValue");
  //positive: binary
  vector<string> generated, gold;
  gold.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmmGold_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad17_small.chp");
  std::string errorMsg ("");

  int diffAllowed;
  std::string prefix;
  double eps;
  bool bCheckHeaders;
  double frac;

  // gold     generated  difference  fraction
  // 714.957  778.957    64.0000000  0.0821611

  // test case: absent vs present optional argument frac==eps
  errorMsg = "";
  Verbose::out(1, "**testcase1a-old behaviour**");    // 1 differences
  // Max diff: 0.00390625 is greater than expected (0.0001) [signal]
  // 2 of 54675 (0.00365798%) were different.
  ChpCheck chpCheckTest1a ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.1 );
  size_t fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(!chpCheckTest1a.check(errorMsg));
  CPPUNIT_ASSERT( ChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 1 );
  errorMsg = "";
  Verbose::out(1, "**testcase1b-new behaviour**");    // 0 differences
  ChpCheck chpCheckTest1b ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.1, 
           bCheckHeaders=true, frac=0.1 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(chpCheckTest1b.check(errorMsg));
  
  Verbose::popMsgHandler();

  return;
}

void ChpCheckTest::testCheckFracDiffSignal()
{
  // message output streamed to string object
  ostringstream MessageHandler_oss;
  MsgStream msgHandler(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler);

  cout<<endl;
  Verbose::out(1, "ChpCheckTest::testCheckFracDiffSignal");
  //positive: binary
  vector<string> generated, gold;
  gold.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1_small.CHP");
  generated.push_back(INPUT+"/chp/exp/100907_CTGF2_U133plus2_t1bad20_small.CHP");
  std::string errorMsg ("");

  int diffAllowed;
  std::string prefix;
  double eps;
  bool bCheckHeaders;
  double frac;

  // gold     generated  difference  fraction
  // 62.8927  62.8929    0.0002000   0.0000032
  // 79.5101  79.514     0.0039000   0.0000490

  // test case: absent vs present optional argument frac==eps
  errorMsg = "";
  Verbose::out(1, "**testcase1a-old behaviour**");    // 2 differences
  // Max diff: 0.00390625 is greater than expected (0.0001) [signal]
  // 2 of 54675 (0.00365798%) were different.
  ChpCheck chpCheckTest1a ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.0001 );
  size_t fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(!chpCheckTest1a.check(errorMsg));
  CPPUNIT_ASSERT( ChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 2 );
  errorMsg = "";
  Verbose::out(1, "**testcase1b-new behaviour**");    // 0 differences
  ChpCheck chpCheckTest1b ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.0001, 
           bCheckHeaders=true, frac=0.0001 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(chpCheckTest1b.check(errorMsg));
  
  // test case: absent vs present optional argument frac<eps
  errorMsg = "";
  Verbose::out(1, "**testcase2a-old behaviour**");    // 1 differences
  ChpCheck chpCheckTest2a ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.001 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(!chpCheckTest2a.check(errorMsg));
  CPPUNIT_ASSERT( ChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 1 );
  errorMsg = "";
  Verbose::out(1, "**testcase2b-new behaviour**");    // 1 differences
  ChpCheck chpCheckTest2b ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.001, 
           bCheckHeaders=true, frac=0.00004 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(!chpCheckTest2b.check(errorMsg));
  CPPUNIT_ASSERT( ChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 1 );

  // test case: absent vs present optional argument frac>eps
  errorMsg = "";
  Verbose::out(1, "**testcase3a-old behaviour**");    //2 differences
  ChpCheck chpCheckTest3a ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.00001 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(!chpCheckTest3a.check(errorMsg));
  CPPUNIT_ASSERT( ChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 2 );
  errorMsg = "";
  Verbose::out(1, "**testcase3b-new behaviour**");    // 0 differences
  ChpCheck chpCheckTest3b ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.00001, 
           bCheckHeaders=true, frac=0.00006 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(chpCheckTest3b.check(errorMsg));
  
  Verbose::popMsgHandler();

  return;
}

void ChpCheckTest::testCheckFracDiffConf()
{
  // message output streamed to string object
  ostringstream MessageHandler_oss;
  MsgStream msgHandler(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler);

  cout<<endl;
  Verbose::out(1, "ChpCheckTest::testCheckFracDiffConf");
  //positive: binary
  vector<string> generated, gold;
  gold.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmmGold_small.chp");
  generated.push_back(INPUT+"/chp/map/CCL-256.1D_NSP.brlmm_bad15_small.chp");
  std::string errorMsg ("");

  int diffAllowed;
  std::string prefix;
  double eps;
  bool bCheckHeaders;
  double frac;

  // gold	generated	difference	fraction
  // 536.5  536.563     0.0630000   0.0001174

  // test case: absent vs present optional argument frac==eps
  errorMsg = "";
  Verbose::out(1, "**testcase1a-old behaviour**");    // 1 differences
  //Max diff: 0.0625 is greater than expected (0.001) [confidence]
  //1 of 262338 (0.000381188%) were different.
  ChpCheck chpCheckTest1a ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.001 );
  size_t fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(!chpCheckTest1a.check(errorMsg));
  CPPUNIT_ASSERT( ChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 1 );
  errorMsg = "";
  Verbose::out(1, "**testcase1b-new behaviour**");    // 0 differences
  ChpCheck chpCheckTest1b ( generated, gold, 
           diffAllowed=0, prefix="apt-", eps=0.001, 
           bCheckHeaders=true, frac=0.001 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT(chpCheckTest1b.check(errorMsg));
  
  Verbose::popMsgHandler();

  return;
}