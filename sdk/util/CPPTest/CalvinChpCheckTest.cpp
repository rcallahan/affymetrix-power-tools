////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
 * @file   CalvinChpCheckTest.cpp
 * @author vliber
 * last change by rsatin on 10/13/09
 *
 * @brief  Testing the CelCheck functions.
 *
 */

//
#include "util/CalvinChpCheck.h"
#include "util/Util.h"
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


class CalvinChpCheckTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( CalvinChpCheckTest );
  
  CPPUNIT_TEST(testHeader_multi_data_gen);
  CPPUNIT_TEST(testData_multi_data_gen);
  CPPUNIT_TEST(testData_Exp);
  CPPUNIT_TEST(testData_multi_data_copynumber);
  CPPUNIT_TEST(testData_multi_data_copynumber_variation);
  CPPUNIT_TEST(testData_multi_data_dmet1);
  CPPUNIT_TEST(testData_multi_data_dmet2);
  CPPUNIT_TEST(testData_multi_data_dmet3);
  CPPUNIT_TEST(testData_multi_data_frac_diff);
  CPPUNIT_TEST(testData_multi_header_frac_diff);
  CPPUNIT_TEST(testData_multi_header_string_diff);
  
  CPPUNIT_TEST_SUITE_END();

private:
  /** utility functions */
  int message_error_count(string message, size_t fpos=0);

public: 
  virtual void tearDown();
  /** unit test functions */
  void testHeader_multi_data_gen();
  void testData_multi_data_gen();
  void testData_Exp();
  void testData_multi_data_copynumber();
  void testData_multi_data_copynumber_variation();
  void testData_multi_data_dmet1();
  void testData_multi_data_dmet2();
  void testData_multi_data_dmet3();
  void testData_multi_data_frac_diff();
  void testData_multi_header_frac_diff();
  void testData_multi_header_string_diff();
  
};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( CalvinChpCheckTest );

void CalvinChpCheckTest::tearDown () 
{
  // clean up after failed test cases (pops any dead objects from stack)
  while( Verbose::getParam().m_MsgHandler.size() > 1 )
    Verbose::popMsgHandler();
}

int CalvinChpCheckTest::message_error_count( string message, size_t fpos ) {
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

void CalvinChpCheckTest::testHeader_multi_data_gen()
{
  cout<<endl;
  cout<<endl;
  Verbose::out(1, "***CalvinChpCheckTest testcases***");
  Verbose::out(1, "CalvinChpCheckTest::testHeader_multi_data_gen");
  
  CHPData chp1, chp2;
  CHPFileReader reader;
  vector<string> generated, gold;
  std::string errorMsg ("");

  //positive 
  gold.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp");
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_small.chp");  
  CalvinChpCheck CalvinChpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(CalvinChpCheckTest1.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest1.check(errorMsg));  

  //negative
  //probe set count > then actual number of the probe sets in the file  000067a0h 7-a  // 000DE136 (original value)
  // file ate the same. Expesting result is positive on both OS, but it is negative on Linux todo vliber
  gold.clear();
  generated.clear();
  gold.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small_bad.chp");
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_small_bad16.chp");  
  CalvinChpCheck CalvinChpCheckTest1a (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest1a.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest1a.check(errorMsg)); 
  
  //Error encountered: CalvinChpCheck::check() - generated and gold vectors must be same size.
  gold.clear();
  generated.clear();
  gold.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_small.chp");
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_small.chp");
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_small.chp");
  CalvinChpCheck CalvinChpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest2.check(errorMsg));
  
  //first character should be 59
  //Error encountered: Error: AGCC library exception:
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad1_small.chp");
  CalvinChpCheck CalvinChpCheckTest3 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest3.check(errorMsg));
  
  //version value is not correct 
  //cannot produce the mesage: Error: different chp versions
  //App. throw an exception with message: Error encountered: Error: AGCC library exception: todo vliber
  //both files return the same version=1, although gen.version=3
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad2_small.chp");
  CalvinChpCheck CalvinChpCheckTest4_1 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest4_1.check(errorMsg));
  
  //# of data groups is not correct
  //Error encountered: Error: Different group counts. 
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3_small.chp");
  CalvinChpCheck CalvinChpCheckTest5_1 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest5_1.check(errorMsg));
  
  //bad chip type
  //message: Error encountered: Error: Different CHP Types: affymetrix-multi-data-type-analysis, affymetrix-}ulti-data-type-analysis
  //Error encountered: Error: Unable to compare CHPs of type affymetrix-}ulti-data-type-analysis
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad4_small.chp");
  CalvinChpCheck CalvinChpCheckTest6_1 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest6_1.check(errorMsg));
  
  //bad array type
  //Error encountered: Error: Different array types. [GenomeWideSNP_6 != GenomeWydeSNP_6]
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad5_small.chp");
  CalvinChpCheck CalvinChpCheckTest7 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest7.check(errorMsg));

  //bad algorithm name
  //Error encountered: Error: Different algorithm names
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad6_small.chp");
  CalvinChpCheck CalvinChpCheckTest8 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest8.check(errorMsg));

  //bad algorithm version
  //Error encountered: Error: Different algorithm version. 
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad7_small.chp");
  CalvinChpCheck CalvinChpCheckTest9 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest9.check(errorMsg));

  //bad summary field
  //bad summary value 
  //bug in the bool ParameterNameValueTypeMostlySame() parameter order is not corect.
  //gold 26.7443 gen 6846.54 eps 0.0001
  //gold female gen gemale
  //Error: Test missing field: 'call_rate' *** row 00004520h/e from 61 to 62  call_rate to cbll_rate
  //Error: for field 'het_rate' expecting: '26.744295' got: '6846.539551'  **** row 000045c0h/1 from 41 to 45
  //Error: for field 'em-cluster-chrX-het-contrast_gender' expecting: 'female' got: 'gemale'   ****  row 00004d10h/9 from 66 to 67
  errorMsg="";
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8_small.chp");
  CalvinChpCheck CalvinChpCheckTest10 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest10.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if(!CalvinChpCheckTest10.check(errorMsg))
	 cout<<errorMsg<<endl;
#endif


  //bad algorithm field
  //name apt-opt-initial-cblls
  //bad algorithm value
  //>gold NA10831_GW6_C.CEL gen NA20831_GW6_C.CEL
  //gold false gen falsf 
  //gold 2.4 gen 3.4
  //gold 0000027696-1179523014-0000013966-0000021724-0000016941 gen 0000027696-1179523014-0000013966-0000021724-0000016041
  //ignore field apt-block-size has been change *** row 00000f60h/e from 63 to 64  apt-blodk-size
  //Error: Test missing field: 'apt-opt-initial-calls' *** row 00002810h/a from 61 to 62  apt-opt-initial-cblls
  //Error: for field 'apt-opt-cel-9' expecting: 'NA10831_GW6_C.CEL' got: 'NA20831_GW6_C.CEL'  *** row 000016d0h/8 from 31 to 32  NA20831_GW6_C.CEL
  //Error: for field 'apt-opt-output-summaries' expecting: 'false' got: 'falsf' *** row 000028b0h/e from 65 to 66  falsf
  //Error: for field 'quantification-version' expecting: '2.4' got: '3.4'  *** row 00003090h/e from 32 to 33  2.4 to 3.4
  //Error: for field 'apt-opt-cel-guid-19' expecting: '0000027696-1179523014-0000013966-0000021724-0000016941' got: '0000027696-1179523014-0000013966-0000021724-0000016041' *** row 00004460h/a from 39 to 30 
  errorMsg="";
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad9_small.chp");
  CalvinChpCheck CalvinChpCheckTest11 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest11.check(errorMsg));
  // APT-992
#ifndef _WIN32
  if(!CalvinChpCheckTest11.check(errorMsg))
	 cout<<errorMsg<<endl; 
#endif
}
void CalvinChpCheckTest::testData_multi_data_gen()
{
  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_data_gen");

  //negative
  // for all files change probe set count to 10 from original 909622 *** row 000067a0h/8-9 to 0 and 000067a0h/a to A
  //bad rows count
  //Error encountered: Wrong number of genotyping probesets. *** row 000067a0h/a from A to B   original: 10 new:11
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_1_small.chp");
  generated.push_back("input/calvin/map/NA06985_GW6_C.birdseed-v2_bad10_small.chp");
  CalvinChpCheck CalvinChpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest1.check(errorMsg)); 
 
  //bad probesetname
  //Error encountered: Different names. gold: 'SNP_A-2131660' test: 'SNQ_A-2131660  *** row 000067b0h/1 from 50 to 51
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad11_small.chp");
  CalvinChpCheck CalvinChpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest2.check(errorMsg));
  
  //bad call for first probe set 
  //Error encountered: Different calls for snp: 'SNP_A-2131660'. gold: 'SNP_A-2131660' test: 'SNP_A-2131660' *** row 000067c0h/0 from 8 to 9
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad12_small.chp");
  CalvinChpCheck CalvinChpCheckTest3 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest3.check(errorMsg));

  //bad confidence for first probe set
  //gold 0.00210923 gen 0.00470674 eps 0.0001 *** row 000067c0h/2 from 0A to 9A
  //Max diff: 0.00259751 is greater than expected (0.0001) [confidence]
  //1 of 10 (10%) were different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad13_small.chp");
  CalvinChpCheck CalvinChpCheckTest4_1 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest4_1.check(errorMsg));

  //positive eps=diff=0.0026
  //gold 0.00210923 gen 0.00470674 eps 0.0026
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad13.chp: Genotype: checked 10 genotype entries.
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad13.chp: Genotype: max confidence diff is: 0.00259751
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad13.chp: Genotype: chip is equivalent.
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad13_small.chp");
  CalvinChpCheck CalvinChpCheckTest4_2 (generated, gold,0,L"apt-",0.0026);
  CPPUNIT_ASSERT(CalvinChpCheckTest4_2.check(errorMsg));

  //positive maxDiff=1
  //gold 0.00210923 gen 0.00470674 eps 0.0001
  //Max diff: 0.00259751 is greater than expected (0.0001) [confidence]
  //1 of 909622 (0.000109936%) were different.
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad13.chp: Genotype: checked 10 genotype entries.
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad13.chp: Genotype: max confidence diff is: 0.00259751
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad13.chp: Genotype: chip is equivalent.
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad13_small.chp");
  CalvinChpCheck CalvinChpCheckTest4_3 (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest4_3.check(errorMsg));
  
  //bad signal vector for second probe set
  //Error encountered: Param vectors different for snp: 'SNP_A-1967418.  ****row 000067e0h 8-f to 0 row 000067f0h /0 to 0
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad14_small.chp");
  CalvinChpCheck CalvinChpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest5.check(errorMsg));
  
  //bad signals
  //Error encountered: Param vectors different for snp: 'SNP_A-1967418'. ****row 000067e0h 9 from A9 to a1
  //gold 4035.93 gen 4999.85 eps 0.0001
  //Error encountered: Param vectors different for snp: 'SNP_A-1969580'. ****row 00006810h 0 from 7c to 9c
  //Error encountered: Param vectors different for snp: 'SNP_A-4263484'. ****row 00006830h 6 from 08 to 02
  //3 of 10 (30%) were different..
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad15.chp: Genotype: checked 10 genotype entries.
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad15.chp: Genotype: max confidence diff is: 0
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad15.chp: Genotype: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad15_small.chp");
  CalvinChpCheck CalvinChpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest6.check(errorMsg));
  
  //3 of 10 (30%) were different.
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad15.chp: Genotype: checked 10 genotype entries.
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad15.chp: Genotype: max confidence diff is: 0
  //input\calvin\map\NA06985_GW6_C.birdseed-v2_bad15.chp: Genotype: chip is equivalent.
  generated.clear();
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad15_small.chp");
  CalvinChpCheck CalvinChpCheckTest7 (generated, gold,3);
  CPPUNIT_ASSERT(CalvinChpCheckTest7.check(errorMsg));
 
}

void CalvinChpCheckTest::testData_Exp()
{
  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_Exp");
  
  //positive

   //input\calvin\expr\3AJW02021805.default.rma.chp: checked 22283 quantification entries.
   //input\calvin\expr\3AJW02021805.default.rma.chp: max signal diff is: 0
   //input\calvin\expr\3AJW02021805.default.rma.chp: chip is equivalent.
   vector<string> generated, gold;
   std::string errorMsg ("");
   gold.push_back(INPUT+"/calvin/expr/3AJW02021805.default.rma_Gold.chp");
   generated.push_back(INPUT+"/calvin/expr/3AJW02021805.default.rma.chp");
   CalvinChpCheck CalvinChpCheckTest1 (generated, gold);
   CPPUNIT_ASSERT(CalvinChpCheckTest1.check(errorMsg));
   POSITIVE_TEST(CalvinChpCheckTest1.check(errorMsg))

  //negative

  //bad rows count
  //Error encountered: Wrong number of expression probesets. *** row 00004ed0h/a from 0B to 1B   original: 22283 new:22299 ***
  generated.clear();
  generated.push_back(INPUT+"/calvin/expr/3AJW02021805.default.rma_Bad1.chp");
  CalvinChpCheck CalvinChpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest2.check(errorMsg));
  
  //bad probesetname *** row 00004ee0h/9 from 35 to 37 and row 00004f00h/4 from 46 to 56 ***
  //Error encountered: Different names. gold: 'AFFX-BioB-5_at' test: 'AFFX-BioB-7_at'
  //Error encountered: Different names. gold: 'AFFX-BioB-M_at' test: 'AFVX-BioB-M_at'
  //2 of 22283 (0.00897545%) were different.
  //input\calvin\expr\3AJW02021805.default.rma_Bad2.chp: checked 22283 quantification entries.
  //input\calvin\expr\3AJW02021805.default.rma_Bad2.chp: max signal diff is: 0
  //input\calvin\expr\3AJW02021805.default.rma_Bad2.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/expr/3AJW02021805.default.rma_Bad2.chp");
  CalvinChpCheck CalvinChpCheckTest3 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest3.check(errorMsg));
  //difCount=2, but result is the same: negative. Is it correct? todo vliber
  //multiDataExpressionSame() should be if(!success && m_DiffAllowed >= numDiff) 
  CalvinChpCheck CalvinChpCheckTest4 (generated,gold,2);
  CPPUNIT_ASSERT(!CalvinChpCheckTest4.check(errorMsg));
  
  //bad quantification *** row 00004ef0h/c from 0C to 5C and row 00004f10h/f from 00 to 01 ***
  //gold 7.72029 gen 7.73005 eps 0.0001
  //gold 8.81257 gen 8.81282 eps 0.0001
  //Max diff: 0.00976563 is greater than expected (0.0001) [quantification]
  //2 of 22283 (0.00897545%) were different.
  //input\calvin\expr\3AJW02021805.default.rma_Bad3.chp: checked 22283 quantification entries.
  //input\calvin\expr\3AJW02021805.default.rma_Bad3.chp: max signal diff is: 0.00976563
  //input\calvin\expr\3AJW02021805.default.rma_Bad3.chp: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/expr/3AJW02021805.default.rma_Bad3.chp");
  CalvinChpCheck CalvinChpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest5.check(errorMsg));
  //difCount=2, but result is the same: negative. Is it correct? todo vliber
  CalvinChpCheck CalvinChpCheckTest6 (generated, gold,2);
  CPPUNIT_ASSERT(!CalvinChpCheckTest6.check(errorMsg));
  //eps=0.01>maxDiff positive result
  //input\calvin\expr\3AJW02021805.default.rma_Bad3.chp: checked 22283 quantification entries.
  //input\calvin\expr\3AJW02021805.default.rma_Bad3.chp: max signal diff is: 0.00976563
  //input\calvin\expr\3AJW02021805.default.rma_Bad3.chp: chip is equivalent.
  CalvinChpCheck CalvinChpCheckTest7 (generated, gold,0,L"apt-",0.01);
  CPPUNIT_ASSERT(CalvinChpCheckTest7.check(errorMsg));
  
}
void CalvinChpCheckTest::testData_multi_data_copynumber()
{
  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_data_copynumber");

  //positive
  
  //input\calvin\copynumber\NA06985_GW6_C.CN5_small.CNCHP: CopyNumber: checked 30 copynumber entries.
  //input\calvin\copynumber\NA06985_GW6_C.CN5_small.CNCHP: CopyNumber: max Log2Ratio diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_small.CNCHP: CopyNumber: max SmoothSignal diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_small.CNCHP: CopyNumber: max AlleleicDifference diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_small.CNCHP: CopyNumber: chip is equivalent.
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_Gold_small.CNCHP");
  generated.push_back("input/calvin/copynumber/NA06985_GW6_C.CN5_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(CalvinChpCheckTest1.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest1.check(errorMsg));

  //negative

  //bad rows count
  //Error encountered: Wrong number of copynumber probesets. *** row 00006430h/b from CF to CE   original: 30 new:17
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad1_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest2.check(errorMsg));
 
  //bad probesetname
  //Error encountered: Different names. gold: 'CN_473981' test: 'CN_T73981  *** row 000064c0h/1 from 34 to 54
  //chip is different
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad2_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest3 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest3.check(errorMsg));
  
  //bad chromosome
  //*** row 000064c0h/b from 01 to 03 gold: 1 test: 3
  //Error encountered: Different Chromosome for snp: 'CN_473981'. gold: 'CN_473981' test: 'CN_473981'  
  //chip is different
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad3_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest4 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest4.check(errorMsg));

  //bad position
  //*** row 000064c0h/f from 23 to 24 gold: 52771 test: 52772
  //Error encountered: Different Position for snp: 'CN_473981'. gold: 'CN_473981' test: 'CN_473981'  
  //chip is different
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad4_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest5.check(errorMsg));

  //bad CNState
  //*** row 000064d0h/0 from 40 to 39
  //Error encountered: Different CNState for snp: 'CN_473981'. gold: 'CN_473981' test: 'CN_473981'
  //chip is different
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad5_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest6.check(errorMsg));
  
  //bad Log2Patio
  // *** row 000067f0h/2 from 3E to 3F
  //Max diff: 0.530446 is greater than expected (0.0001) [Log2Ratio]
  //1 of 30 (3.33333%) were different.
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad6_small.CNCHP: CopyNumber: checked 30 copynumber entries.
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad6_small.CNCHP: CopyNumber: max Log2Ratio diff is: 0.530446
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad6_small.CNCHP: CopyNumber: max SmoothSignal diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad6_small.CNCHP: CopyNumber: max AlleleicDifference diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad6_small.CNCHP: CopyNumber: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad6_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest7 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest7.check(errorMsg));
   
  //max Log2Ratio diff is: 0.530446
  //chip is equivalent.
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad6_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest7a (generated, gold,0,L"apt-",0.55);
  CPPUNIT_ASSERT(CalvinChpCheckTest7a.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest7a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest7b (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest7b.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest7b.check(errorMsg));
 
  //bad SmoothSignal
  // *** row 000067f0h/6 from 3D to 3E
  //Max diff: 0.157968 is greater than expected (0.0001) [SmoothSignal]
  //1 of 30 (3.33333%) were different.
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad7_small.CNCHP: CopyNumber: checked 30 copynumber entries.
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad7_small.CNCHP: CopyNumber: max Log2Ratio diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad7_small.CNCHP: CopyNumber: max SmoothSignal diff is: 0.157968
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad7_small.CNCHP: CopyNumber: max AlleleicDifference diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad7_small.CNCHP: CopyNumber: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad7_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest8 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest8.check(errorMsg));
  
  //max SmoothSignal diff is: 0.157968
  //chip is equivalent.
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad7_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest8a (generated, gold, 0,L"apt-",0.16);
  CPPUNIT_ASSERT(CalvinChpCheckTest8a.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest8a.check(errorMsg));   

  //Max diff: 0.157968 is greater than expected (0.0001) [SmoothSignal]
  //1 of 30 (3.33333%) were different.
  //chip is equivalent.
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad7_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest8b (generated, gold, 1);
  CPPUNIT_ASSERT(CalvinChpCheckTest8b.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest8b.check(errorMsg));
 
  //bad AlleleDifference
  //*** row 00006800h/1 from 2B to 2C
  //Max diff: 1.19209e-007 is greater than expected (0) [AllelicDifference]
  //1 of 30 (3.33333%) were different.
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad8_small.CNCHP: CopyNumber: checked 30 copynumber entries.
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad8_small.CNCHP: CopyNumber: max Log2Ratio diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad8_small.CNCHP: CopyNumber: max SmoothSignal diff is: 0
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad8_small.CNCHP: CopyNumber: max AlleleicDifference diff is: 1.19209e-007
  //input\calvin\copynumber\NA06985_GW6_C.CN5_bad8_small.CNCHP: CopyNumber: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad8_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest9 (generated, gold,0,L"apt-",0);
  CPPUNIT_ASSERT(!CalvinChpCheckTest9.check(errorMsg));
  
  // max AlleleicDifference diff is: 1.19209e-007
  //chip is equivalent
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad8_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest9a (generated, gold, 0,L"apt-",0.1);
  CPPUNIT_ASSERT(CalvinChpCheckTest9a.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest9a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest9b (generated, gold, 1);
  CPPUNIT_ASSERT(CalvinChpCheckTest9b.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest9b.check(errorMsg));

  //bad LOH
  //*** row 000067f0h/d from 00 to 01
  //Error encountered: Different LOH for snp: 'SNP_A-8575125' gold: 'SNP_A-8575125' test: 'SNP_A-8575125'
  //1 of 30 (3.33333%) were different.
  //chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_bad9_small.CNCHP");
  CalvinChpCheck CalvinChpCheckTest10 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest10.check(errorMsg));
}

void CalvinChpCheckTest::testData_multi_data_copynumber_variation()
{
  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_data_copynumber_variation");;
  
  //positive
  
  //input\calvin\copynumber\GenomeWideSNP_6.agcc.multidata_copynumber_var_Gold.chp: CopyNumberVariation: checked 1204 CNV entries.
  //input\calvin\copynumber\GenomeWideSNP_6.agcc.multidata_copynumber_var_Gold.chp: CopyNumberVariation: max Signal diff is: 0
  //input\calvin\copynumber\GenomeWideSNP_6.agcc.multidata_copynumber_var_Gold.chp: CopyNumberVariation: max Confidence diff is: 0
  //input\calvin\copynumber\GenomeWideSNP_6.agcc.multidata_copynumber_var_Gold.chp: CopyNumberVariation: chip is equivalent.
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var.chp");
  generated.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var_Gold.chp");
  CalvinChpCheck CalvinChpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(CalvinChpCheckTest1.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest1.check(errorMsg));

  //negative

  //bad rows count
  //Error encountered: Error encountered: Wrong number of CNV probesets. *** row 00003a10h/1 from 04 to 03   original: 1204 new:948
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var_bad1.chp");
  CalvinChpCheck CalvinChpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest2.check(errorMsg));
 
  //bad probesetname
  //Error encountered: Different names. gold: 'CNP12' test: 'CN`12  *** row 00003a80h/2 from 50 to 60
  //chip is different
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var_bad2.chp");
  CalvinChpCheck CalvinChpCheckTest3 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest3.check(errorMsg));
  
  //bad signal
  //Error encountered: Different signals. gold: '1.3728' test: '1.4978
  //Max diff: 0.125 is greater than expected (0.0001) [signal]
  //1 of 1204 (0.0830565%) were different.
  //chip is different
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var_bad3.chp");
  CalvinChpCheck CalvinChpCheckTest4 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest4.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest4a (generated, gold,0,L"apt-",0.125);
  CPPUNIT_ASSERT(CalvinChpCheckTest4a.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest4a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest4b (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest4b.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest4b.check(errorMsg));

  //bad call
  //Error encountered: Different calls. gold: '' test: '  *** row 00003a20h/3 from 02 to 03 todo vliber bad output value
  //1 of 1204 (0.0830565%) were different.
  //chip is different
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var_bad4.chp");
  CalvinChpCheck CalvinChpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest5.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest5a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest5a.check(errorMsg));
  POSITIVE_TEST(CalvinChpCheckTest5a.check(errorMsg));
  
  //bad confidence
  //Error encountered: Different confidences. gold: '1' test: '1.125
  //Max diff: 0.125 is greater than expected (0.0001) [confidence]
  //1 of 1204 (0.0830565%) were different.
  //chip is different
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var_bad5.chp");
  CalvinChpCheck CalvinChpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest6.check(errorMsg));

  //positive

  //chip is equivalent
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var_bad5.chp");
  CalvinChpCheck CalvinChpCheckTest6a (generated, gold, 0,L"apt-",0.13);
  CPPUNIT_ASSERT(CalvinChpCheckTest6a.check(errorMsg));
  //chip is equivalent
  generated.clear();
  generated.push_back(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var_bad5.chp");
  CalvinChpCheck CalvinChpCheckTest6b (generated, gold, 1);
  CPPUNIT_ASSERT(CalvinChpCheckTest6b.check(errorMsg));
  
}


void CalvinChpCheckTest::testData_multi_data_dmet1()
{

  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_data_dmet::DmetBiAllelicMultiDataType");

  //positive

  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp");
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet.chp");
  CalvinChpCheck CalvinChpCheckTest1 (generated, gold);
  CPPUNIT_ASSERT(CalvinChpCheckTest1.check(errorMsg));  
  POSITIVE_TEST(CalvinChpCheckTest1.check(errorMsg));
 
  //negative

  //bad rows count
  //Error encountered: Wrong number of genotyping probesets. *** row 00007c20h/f from 6E to 6F  orig 1902 new 1903 
   generated.clear();
   generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad7.chp");
   CalvinChpCheck CalvinChpCheckTest8 (generated, gold);
   CPPUNIT_ASSERT(!CalvinChpCheckTest8.check(errorMsg));

  //bad probesetname
  //Error encountered: Different names. gold: 'AM_10001' test: 'QM_10001  *** row 00007c30h/4 from 41 to 51
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad1.chp");
  CalvinChpCheck CalvinChpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest2.check(errorMsg));
  
  //bad call for first probe set 
  //Error encountered: Different calls for cnv: 'AM_10001'. gold: ' ' test: '' *** 00007c30h/e from FF to EF todo vliber no output values
  //1 of 1902 (0.0525762%) were different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad2.chp");
  CalvinChpCheck CalvinChpCheckTest3 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest3.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest3a(generated,gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest3a.check(errorMsg)); 
  
  //bad confidance for first probe set 
  //Max diff: 8.58993e+009 is greater than expected (0.0001) [confidence] *** 00007c30h/f from C0 to D0.
  //1 of 1902 (0.0525762%) were different.
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad8.chp: DmetBiAllelic: checked 1902 genotype entries.
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad8.chp: DmetBiAllelic: max confidence diff is: 8.58993e+009
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad8.chp: DmetBiAllelic: max signal diff is: 0
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad8.chp: DmetBiAllelic: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad8.chp");
  CalvinChpCheck CalvinChpCheckTest9 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest9.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest9a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest9a.check(errorMsg));

  //bad forced call for first probe set 
  //Error encountered: Different forced calls for cnv: 'AM_10001'. gold: ' ' test: '' *** 00007c40h/3 from FF to EF todo vliber no output values
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad3.chp");
  CalvinChpCheck CalvinChpCheckTest4 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest4.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest4a(generated,gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest4a.check(errorMsg));

  //bad contextA call for first probe set 
  //Error encountered: Different contextA for cnv: 'AM_10001'. gold: ' ' test: '' *** 00007c40h/d from FF to EF todo vliber no output values
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad4.chp");
  CalvinChpCheck CalvinChpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest5.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest5a(generated,gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest5a.check(errorMsg));
  
  //bad contextB call for first probe set 
  //Error encountered: Different contextB for cnv: 'AM_10001'. gold: ' ' test: '' *** 00007c40h/c from FF to EF todo vliber no output values
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad5.chp");
  CalvinChpCheck CalvinChpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest6.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest6a(generated,gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest6a.check(errorMsg));

  //bad signalA call for first probe set 
  //*** 00007c60h/2 from C0 to D0. 
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad6.chp");
  CalvinChpCheck CalvinChpCheckTest7 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest7.check(errorMsg));

  //bad signalB call for first probe set 
  //*** 00007c60h/6 from C0 to A0. 
  //Max diff: 2 is greater than expected (0.0001) [signal]
  //1 of 1902 (0.0525762%) were different.
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad9.chp: DmetBiAllelic: checked 1902 genotype entries.
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad9.chp: DmetBiAllelic: max confidence diff is: 0
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad9.chp: DmetBiAllelic: max signal diff is: 2
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad9.chp: DmetBiAllelic: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad9.chp");
  CalvinChpCheck CalvinChpCheckTest10 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest10.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest10a (generated, gold, 0,L"apt-",2);
  CPPUNIT_ASSERT(CalvinChpCheckTest10a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest10b (generated, gold, 1);
  CPPUNIT_ASSERT(CalvinChpCheckTest10b.check(errorMsg));
}

void CalvinChpCheckTest::testData_multi_data_dmet2()
{

  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_data_dmet::DmetCopyNumberMultiDataType");
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp");

  //negative

  //bad rows count
  //Error encountered: Wrong number of genotyping probesets. *** row 000163c0h/e from 05 to 04  orig 5 new 4 
   generated.clear();
   generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad10.chp");
   CalvinChpCheck CalvinChpCheckTest1 (generated, gold);
   CPPUNIT_ASSERT(!CalvinChpCheckTest1.check(errorMsg));

  //bad probesetname
  //Error encountered: Different names. gold: 'CN_CYP2A6' test: 'CN_3YP2A6  *** row 000163d0h/6 from 43 to 33
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad11.chp");
  CalvinChpCheck CalvinChpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest2.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest2a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest2a.check(errorMsg));
  
  //bad calls
  //Error encountered: Different calls for cnv: 'CN_CYP2A6'. gold: '-2' test: '-18'  *** row 000163d0h/e from EE to FE
  //1 of 5 (20%) were different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad12.chp");
  CalvinChpCheck CalvinChpCheckTest3 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest3.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest3a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest3a.check(errorMsg));

  //bad confidance
  //Max diff: 8.58993e+009 is greater than expected (0.0001) [confidence]  *** row 000163d0h/f from C0 to D0
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad13.chp");
  CalvinChpCheck CalvinChpCheckTest4 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest4.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest4a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest4a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest4b (generated, gold,0,L"apt-",10000000000.0);
  CPPUNIT_ASSERT(CalvinChpCheckTest4b.check(errorMsg));

  //bad forse calls
  //Error encountered: Different forced calls for cnv: 'CN_CYP2A6'. gold: '-2' test: '-82'  *** row 000163e0h/4 from FE to AE
  //1 of 5 (20%) were different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad14.chp");
  CalvinChpCheck CalvinChpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest5.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest5a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest5a.check(errorMsg));

  //bad lower calls
  //Error encountered: Different lower for cnv: 'CN_CYP2A6'. gold: '-2' test: '-2.00098'  *** row 000163e0h/b from 00 to 10
  //1 of 5 (20%) were different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad15.chp");
  CalvinChpCheck CalvinChpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest6.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest6a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest6a.check(errorMsg));

  //bad upper calls
  //Error encountered: Different upper for cnv: 'CN_CYP2A6'. gold: '-2' test: '-2.01563'  *** row 00016400h/e from 00 to 01
  //1 of 5 (20%) were different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad16.chp");
  CalvinChpCheck CalvinChpCheckTest7 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest7.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest7a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest7a.check(errorMsg));

  //Max diff: 0.000976563 is greater than expected (0.0001) [signal]
  //1 of 5 (20%) were different.
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad17.chp: DmetCopyNumber: checked 5 genotype entries.
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad17.chp: DmetCopyNumber: max confidence diff is: 0
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad17.chp: DmetCopyNumber: max signal diff is: 0.000976563
  //input\calvin\dmet\ZQ_D10_gCtrl_1.dmet_bad17.chp: DmetCopyNumber: chip is different.
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad17.chp");
  CalvinChpCheck CalvinChpCheckTest8 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest8.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest8a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest8a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest8b (generated, gold,0,L"apt-",0.001);
  CPPUNIT_ASSERT(CalvinChpCheckTest8b.check(errorMsg));
}

void CalvinChpCheckTest::testData_multi_data_dmet3()
{

  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_data_dmet::DmetMultiAllelicMultiDataType");
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp");

  //negative

  //bad rows count
  //Error encountered: Wrong number of genotyping probesets. *** row 00015d10h/3 from 1D to 1E  orig 29 new 4 
   generated.clear();
   generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad18.chp");
   CalvinChpCheck CalvinChpCheckTest1 (generated, gold);
   CPPUNIT_ASSERT(!CalvinChpCheckTest1.check(errorMsg));

  //bad probesetname
  //Error encountered: Different names. gold: 'AM_10002' test: 'AO_10002  *** row 00015d10h/9 from 4D to 4F
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad19.chp");
  CalvinChpCheck CalvinChpCheckTest2 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest2.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest2a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest2a.check(errorMsg));

  //bad call
  //Error encountered: Different calls for cnv: 'AM_10002'. gold: ' ' test: ''todo vliber not regognizable characters in output
  //1 of 29 (3.44828%) were different.*** row 00015d20h/2 from FF to EF
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad20.chp");
  CalvinChpCheck CalvinChpCheckTest3 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest3.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest3a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest3a.check(errorMsg));

  //bad confidence
  //Max diff: 0.000976563 is greater than expected (0.0001) [confidence]
  //1 of 29 (3.44828%) were different.*** row 00015d20h/5 from 00 to 10
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad21.chp");
  CalvinChpCheck CalvinChpCheckTest4 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest4.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest4a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest4a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest4b (generated, gold,0,L"apt-",0.001);
  CPPUNIT_ASSERT(CalvinChpCheckTest4b.check(errorMsg));

  //bad farce call
  //Error encountered: Different forced calls for cnv: 'AM_10002'. gold: ' ' test: '' todo liber bad character
  //1 of 29 (3.44828%) were different.*** row 00015d20h/7 from FF to EF
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad22.chp");
  CalvinChpCheck CalvinChpCheckTest5 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest5.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest5a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest5a.check(errorMsg));

  //bad allele count
  //Error encountered: Different alleleCount for cnv: 'AM_10002'. gold: '' test: '' todo liber bad character
  //1 of 29 (3.44828%) were different.*** row 00015d20h/8 from 03 to 04
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad23.chp");
  CalvinChpCheck CalvinChpCheckTest6 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest6.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest6a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest6a.check(errorMsg));
  
  //bad signal
  //Max diff: 0.015625 is greater than expected (0.0001) [signal] 
  //1 of 29 (3.44828%) were different.*** row 00015d20h/a from 00 to 01
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad24.chp");
  CalvinChpCheck CalvinChpCheckTest7 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest7.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest7a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest7a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest7b (generated, gold,0,L"apt-",0.016);
  CPPUNIT_ASSERT(CalvinChpCheckTest7b.check(errorMsg));

  //bad signal
  //Max diff: 0.25 is greater than expected (0.0001) [signal]
  //1 of 29 (3.44828%) were different.*** row 00015d20h/e from 00 to 10
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad25.chp");
  CalvinChpCheck CalvinChpCheckTest8 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest8.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest8a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest8a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest8b (generated, gold,0,L"apt-",0.25);
  CPPUNIT_ASSERT(CalvinChpCheckTest8b.check(errorMsg));

  //bad signal
  //Max diff: 0.75 is greater than expected (0.0001) [signal]
  //1 of 29 (3.44828%) were different.*** row 00015d30h/2 from 00 to 30, 
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad26.chp");
  CalvinChpCheck CalvinChpCheckTest9 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest9.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest9a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest9a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest9b (generated, gold,0,L"apt-",0.75);
  CPPUNIT_ASSERT(CalvinChpCheckTest9b.check(errorMsg));
  
  //bad signal
  //Max diff: 1.25 is greater than expected (0.0001) [signal]
  //1 of 29 (3.44828%) were different.
  //*** row 00015d20h/a from 00 to 10, 00015d20h/e from 00 to 40, 00015d30h/2 from 00 to 50
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_dmet_bad24-26.chp");
  CalvinChpCheck CalvinChpCheckTest10 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest10.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest10a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest10a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest10b (generated, gold,0,L"apt-",1.25);
  CPPUNIT_ASSERT(CalvinChpCheckTest10b.check(errorMsg));

  //Different multiple context
  //1 of 29 (3.44828%) were different.
  //*** row 00015d40h/1 from FF to AF, 00015d40h/2 from FF to BF, 00015d40h/3 from FF to CF
  //*** row 00015d40h/4 from FF to DF, 00015d40h/5 from EF to 40, 00015d40h/6 from FF to 9F
  //Error encountered: Different contextA for cnv: 'AM_10002'. gold: ' ' test: '' todo vliber bad output character
  //Error encountered: Different contextB for cnv: 'AM_10002'. gold: ' ' test: '' todo vliber bad output character
  //Error encountered: Different contextC for cnv: 'AM_10002'. gold: ' ' test: '' todo vliber bad output character
  //Error encountered: Different contextD for cnv: 'AM_10002'. gold: ' ' test: '' todo vliber bad output character
  //Error encountered: Different contextE for cnv: 'AM_10002'. gold: ' ' test: '' todo vliber bad output character
  //Error encountered: Different contextF for cnv: 'AM_10002'. gold: ' ' test: '' todo vliber bad output character
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad27.chp");
  CalvinChpCheck CalvinChpCheckTest11 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest11.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest11a (generated, gold,1);
  CPPUNIT_ASSERT(CalvinChpCheckTest11a.check(errorMsg));
  
  // multiple SNP with differen errors
  //*** row 00015d20h/2 from FF to EF, 00015d60h/6 from 10 to 30, 00015da0h/7 from FF to AF
  //Error encountered: Different calls for cnv: 'AM_10002'. gold: ' ' test: ''
  //Error encountered: Different contextA for cnv: 'AM_10131'. gold: ' ' test: ''
  //Max diff: 0.0078125 is greater than expected (0.0001) [signal]
  generated.clear();
  generated.push_back(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_bad28.chp");
  CalvinChpCheck CalvinChpCheckTest12 (generated, gold);
  CPPUNIT_ASSERT(!CalvinChpCheckTest12.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest12a (generated, gold,1);
  CPPUNIT_ASSERT(!CalvinChpCheckTest12a.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest12b (generated, gold,2);
  CPPUNIT_ASSERT(!CalvinChpCheckTest12b.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest12c (generated, gold,3);
  CPPUNIT_ASSERT(CalvinChpCheckTest12c.check(errorMsg));
  CalvinChpCheck CalvinChpCheckTest12d (generated, gold,2,L"apt-",0.008);
  CPPUNIT_ASSERT(CalvinChpCheckTest12d.check(errorMsg));
 }

void CalvinChpCheckTest::testData_multi_data_frac_diff()
{
  // message output streamed to string object
  ostringstream MessageHandler_oss;
  MsgStream msgHandler(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler);

  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_data_frac_diff");
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp");
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp");

  int diffAllowed;
  std::wstring prefix;
  double eps;
  bool bCheckHeaders;
  double fraction;

  // test case: absent vs present optional argument frac==eps
  Verbose::out(1, "**testcase1a-old behaviour**");    // 2 differences
  CalvinChpCheck test1a (generated, gold, diffAllowed=0, prefix=L"apt-", eps=0.001 );
  size_t fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT( !test1a.check(errorMsg) );
  CPPUNIT_ASSERT( CalvinChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 2 );
  Verbose::out(1, "**testcase1b-new behaviour**");    // 0 differences
  CalvinChpCheck test1b (generated, gold, diffAllowed=0, prefix=L"apt-", eps=0.001, bCheckHeaders=true, fraction=0.001 );
  CPPUNIT_ASSERT( test1b.check(errorMsg) );

  // test case: absent vs present optional argument frac>eps
  Verbose::out(1, "**testcase2a-old behaviour**");    // 19 differences
  CalvinChpCheck test2a (generated, gold, diffAllowed=0, prefix=L"apt-", eps=0.0001 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT( !test2a.check(errorMsg) );
  CPPUNIT_ASSERT( CalvinChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 19 );
  Verbose::out(1, "**testcase2b-new behaviour**");    // 0 differences
  CalvinChpCheck test2b (generated, gold, diffAllowed=0, prefix=L"apt-", eps=0.0001, bCheckHeaders=true, fraction=0.001 );
  CPPUNIT_ASSERT( test2b.check(errorMsg) );

  // test case: absent vs present optional argument frac<eps
  Verbose::out(1, "**testcase3a-old behaviour**");    // 19 differences
  CalvinChpCheck test3a (generated, gold, diffAllowed=0, prefix=L"apt-", eps=0.0001 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT( !test3a.check(errorMsg) );
  CPPUNIT_ASSERT( CalvinChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 19 );
  fpos_start = MessageHandler_oss.str().size();
  Verbose::out(1, "**testcase3b-new behaviour**");    // 6 differences
  CalvinChpCheck test3b (generated, gold, diffAllowed=0, prefix=L"apt-", eps=0.0001, bCheckHeaders=true, fraction=0.00001 );
  CPPUNIT_ASSERT( !test3b.check(errorMsg) );
  CPPUNIT_ASSERT( CalvinChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 6 );

  Verbose::popMsgHandler();

  return;
}

void CalvinChpCheckTest::testData_multi_header_frac_diff()
{
  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_data_frac_diff");
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp");
  generated.push_back(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8c_small.chp");

  int diffAllowed;
  std::wstring prefix;
  double eps;
  bool bCheckHeaders;
  double fraction;

  // test case: absent vs present optional argument frac==eps
  errorMsg = "";
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "**testcase1a-old behaviour**");    // 2 differences
#endif
  CalvinChpCheck test1a (generated, gold, diffAllowed=0, prefix=L"affymetrix-", eps=0.001 );
  CPPUNIT_ASSERT( !test1a.check(errorMsg) );
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "errorMsg: " + errorMsg );
#endif
  // Error: for field 'het_rate' expecting: '26.744295' got: '26.705233'
  // Error: for field 'hom_rate' expecting: '73.161270' got: '73.172989'
  CPPUNIT_ASSERT( errorMsg.find("'het_rate'") != string::npos );
  CPPUNIT_ASSERT( errorMsg.find("'hom_rate'") != string::npos );
  errorMsg = "";
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "**testcase1b-new behaviour**");    // 1 differences
#endif
  CalvinChpCheck test1b (generated, gold, diffAllowed=0, prefix=L"affymetrix-", eps=0.001, bCheckHeaders=true, fraction=0.001 );
  CPPUNIT_ASSERT( !test1b.check(errorMsg) );
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "errorMsg: " + errorMsg );
#endif
  // Error: for field 'het_rate' expecting: '26.744295' got: '26.705233'
  CPPUNIT_ASSERT( errorMsg.find("'het_rate'") != string::npos );
  CPPUNIT_ASSERT( errorMsg.find("'hom_rate'") == string::npos );

  // test case: absent vs present optional argument frac>eps
  errorMsg = "";
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "**testcase2a-old behaviour**");    // 2 differences
#endif
  CalvinChpCheck test2a (generated, gold, diffAllowed=0, prefix=L"affymetrix-", eps=0.001 );
  CPPUNIT_ASSERT( !test2a.check(errorMsg) );
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "errorMsg: " + errorMsg );
#endif
  CPPUNIT_ASSERT( errorMsg.find("'het_rate'") != string::npos );
  CPPUNIT_ASSERT( errorMsg.find("'hom_rate'") != string::npos );
  errorMsg = "";
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "**testcase2b-new behaviour**");    // 0 differences
#endif
  CalvinChpCheck test2b (generated, gold, diffAllowed=0, prefix=L"affymetrix-", eps=0.001, bCheckHeaders=true, fraction=0.002 );
  CPPUNIT_ASSERT( test2b.check(errorMsg) );

  // test case: absent vs present optional argument frac<eps
  errorMsg = "";
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "**testcase3a-old behaviour**");    // 1 differences
#endif
  CalvinChpCheck test3a (generated, gold, diffAllowed=0, prefix=L"affymetrix-", eps=0.02 );
  CPPUNIT_ASSERT( !test3a.check(errorMsg) );
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "errorMsg: " + errorMsg );
#endif
  CPPUNIT_ASSERT( errorMsg.find("'het_rate'") != string::npos );
  CPPUNIT_ASSERT( errorMsg.find("'hom_rate'") == string::npos );
  errorMsg = "";
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "**testcase3b-new behaviour**");    // 0 differences
#endif
  CalvinChpCheck test3b (generated, gold, diffAllowed=0, prefix=L"affymetrix-", eps=0.02, bCheckHeaders=true, fraction=0.002 );
  CPPUNIT_ASSERT( test3b.check(errorMsg) );  

  return;
}

void CalvinChpCheckTest::testData_multi_header_string_diff()
{
  // message output streamed to string object
  ostringstream MessageHandler_oss;
  MsgStream msgHandler(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler);

  cout<<endl;
  Verbose::out(1, "CalvinChpCheckTest::testData_multi_header_string_diff");
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/calvin/map/NA06985_AxiomGWASHuSNP1_HeaderDiff1.chp");
  generated.push_back(INPUT+"/calvin/map/NA06985_AxiomGWASHuSNP1_HeaderDiff2.chp");

  int diffAllowed;
  std::wstring prefix;
  double eps;
  bool bCheckHeaders;
  double fraction;

  errorMsg = "";
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "**testcase1a-new behaviour**");    // 2 SNPs, 2 header fields
#endif
  CalvinChpCheck test1a (generated, gold, diffAllowed=0, prefix=L"apt-", eps=0.001 );
  size_t fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT( !test1a.check(errorMsg) );
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "errorMsg: " + errorMsg );
#endif
  CPPUNIT_ASSERT( CalvinChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 2 );
  // Error: for field 'apt-program-name' expecting: 'apt-probeset-genotype' got: 'apt-probeset-genotype.exe' 
  // Error: for field 'apt-command-line' expecting: '../../releases/2009-08-19_apt_head/os-x/apt-probeset-genotype --cc-chp-output ...'
  CPPUNIT_ASSERT( errorMsg.find("'apt-program-name'") != string::npos );
  CPPUNIT_ASSERT( errorMsg.find("'apt-opt-cc-md-chp-out-dir'") != string::npos );

  errorMsg = "";
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "**testcase1b-new behaviour**");    // 0 SNPs, 2 header fields
#endif
  CalvinChpCheck test1b (generated, gold, diffAllowed=0, prefix=L"apt-", eps=0.001, bCheckHeaders=true, fraction=0.01 );
  fpos_start = MessageHandler_oss.str().size();
  CPPUNIT_ASSERT( !test1b.check(errorMsg) );
  // APT-992
#ifndef _WIN32
  Verbose::out(1, "errorMsg: " + errorMsg );
#endif
  CPPUNIT_ASSERT( CalvinChpCheckTest::message_error_count( MessageHandler_oss.str(), fpos_start ) == 0 );
  // Error: for field 'apt-program-name' expecting: 'apt-probeset-genotype' got: 'apt-probeset-genotype.exe' 
  // Error: for field 'apt-command-line' expecting: '../../releases/2009-08-19_apt_head/os-x/apt-probeset-genotype --cc-chp-output ...'
  CPPUNIT_ASSERT( errorMsg.find("'apt-program-name'") != string::npos );
  CPPUNIT_ASSERT( errorMsg.find("'apt-opt-cc-md-chp-out-dir'") != string::npos );

  Verbose::popMsgHandler();

  return;
}
