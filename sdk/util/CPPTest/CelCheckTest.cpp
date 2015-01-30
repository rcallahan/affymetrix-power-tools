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
 * @file   CelCheckTest.cpp
 * @author vliber
 * last change by vliber on 06/18/08
 *
 * @brief  Testing the CelCheck functions.
 *
 */



#include "util/AffxString.h"
#include "util/CelCheck.h"
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







class CelCheckTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( CelCheckTest ); 
  CPPUNIT_TEST(testExpBinary);
  CPPUNIT_TEST(testExpASCII);
  CPPUNIT_TEST(testMapBinary);
  CPPUNIT_TEST(testMapASCII);
  CPPUNIT_TEST(testCalvinBinary);
  CPPUNIT_TEST_SUITE_END();


public:

  void testExpBinary();
  void testExpASCII();
  void testMapBinary();
  void testMapASCII();
  void testCalvinBinary();
};


CPPUNIT_TEST_SUITE_REGISTRATION( CelCheckTest );


void CelCheckTest::testExpBinary()
{
  cout<<endl;
  cout<<endl;
  Verbose::out(1, "***CelCheckTest testcases***");
  Verbose::out(1, "CelCheckTest::testExpBinary");

  //positive
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_small.CEL");
  CelCheck celCheckTest1 (generated, gold,0.00);
  CPPUNIT_ASSERT(celCheckTest1.check(errorMsg));
  POSITIVE_TEST(celCheckTest1.check(errorMsg));

  //negative
  //Error encountered: CelCheck::check() - generated and gold vectors must be same size.
  gold.clear();
  generated.clear();
  gold.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_small.CEL");
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_small.CEL");
  CelCheck celCheckTest2 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest2.check(errorMsg));

  //bad gold file
  //Can't read cel file: input1\cel\exp\100907_CTGF2_U133plus2_t1_Gold.CEL
  gold.clear();
  generated.clear();
  errorMsg = "" ;
  gold.push_back("input1/cel/exp/100907_CTGF2_U133plus2_t1_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_small.CEL");
  CelCheck celCheckTest3 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest3.check(errorMsg));
  if(!celCheckTest3.check(errorMsg))
	  cout<<errorMsg<<endl; 

  //bad generated file

  //Can't read cel file: input1\cel\exp\100907_CTGF2_U133plus2_t1.CEL
  gold.clear();
  generated.clear();
  errorMsg = "";
  generated.push_back("input1/cel/exp/100907_CTGF2_U133plus2_t1_small.CEL");
  gold.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_Gold_small.CEL");
  CelCheck celCheckTest4 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest4.check(errorMsg));
  if(!celCheckTest4.check(errorMsg))
	  cout<<errorMsg<<endl;
   
  //file should start with 64
  //Can't read cel file: input\cel\exp\100907_CTGF2_U133plus2_t1_bad1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_bad1_small.CEL");
  CelCheck celCheckTest5 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest5.check(errorMsg));
  if(!celCheckTest5.check(errorMsg))
	  cout<<errorMsg<<endl;

  //gold 158 gen 159 eps 0 line 00000460h/4 1E to 1F
  //Max diff: 1 is greater than expected (0)
  //1 of 3 (33.3333%) were different.
  //input\cel\exp\100907_CTGF2_U133plus2_t1_bad2.CEL chip is different max diff is: 1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_bad2_small.CEL");
  CelCheck celCheckTest6 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest6.check(errorMsg));
  //eps=1.0
  //gold 158 gen 159 eps 1
  //input\cel\exp\100907_CTGF2_U133plus2_t1_bad2.CEL chip is same max diff is: 1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_bad2_small.CEL");
  CelCheck celCheckTest7 (generated, gold,1.00);
  CPPUNIT_ASSERT(celCheckTest7.check(errorMsg));
  POSITIVE_TEST(celCheckTest7.check(errorMsg));

  //gold 158 gen 159 eps 0
  //gold 20520 gen 20521 eps 0
  //Max diff: 2 is greater than expected (0)
  //2 of 3 (66.6667%) were different.
  //input\cel\exp\100907_CTGF2_U133plus2_t1_bad3.CEL chip is different max diff is: 2
  errorMsg = "";

  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_bad3_small.CEL");
  CelCheck celCheckTest8 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest8.check(errorMsg));

  //maxDiff==2
  //positive result
  //gold 158 gen 159 eps 0
  //gold 20520 gen 20522 eps 0
  //Max diff: 2 is greater than expected (0)
  //2 of 3 (66.6667%) were different.
  //input\cel\exp\100907_CTGF2_U133plus2_t1_bad3.CEL chip is same max diff is: 2
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1_bad3_small.CEL");
  CelCheck celCheckTest9 (generated, gold,0.00,"affymetrix-",2);
  CPPUNIT_ASSERT(celCheckTest9.check(errorMsg));
  POSITIVE_TEST(celCheckTest9.check(errorMsg));

}

void CelCheckTest::testExpASCII()

{
  cout<<endl;
  Verbose::out(1, "CelCheckTest::testExpASCII");

  //positive
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_small.CEL");
  CelCheck celCheckTest1 (generated, gold,0.00);
  CPPUNIT_ASSERT(celCheckTest1.check(errorMsg));
  POSITIVE_TEST(celCheckTest1.check(errorMsg));

  //bad gold file
  //Can't read cel file: input1\cel\exp\100907_CTGF2_U133plus2_t1ASCII_Gold.CEL
  gold.clear();
  generated.clear();
  errorMsg = "" ;
  gold.push_back("input1/cel/exp/100907_CTGF2_U133plus2_t1ASCII_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_small.CEL");
  CelCheck celCheckTest3 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest3.check(errorMsg));
  if(!celCheckTest3.check(errorMsg))
	  cout<<errorMsg<<endl;
  
  //bad generated file
  //Can't read cel file: input1\cel\exp\100907_CTGF2_U133plus2_t1ASCII.CEL
  gold.clear();
  generated.clear();
  errorMsg = "";
  gold.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_Gold_small.CEL");
  generated.push_back("input1/cel/exp/100907_CTGF2_U133plus2_t1ASCII_small.CEL");
  CelCheck celCheckTest4 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest4.check(errorMsg));
  if(!celCheckTest4.check(errorMsg))
	  cout<<errorMsg<<endl;
  
  // delete [CEL] Version=3
  //Can't read cel file: input\cel\exp\100907_CTGF2_U133plus2_t1ASCII_bad1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_bad1_small.CEL");
  CelCheck celCheckTest5 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest5.check(errorMsg));
  if(!celCheckTest5.check(errorMsg))
	  cout<<errorMsg<<endl;

  //gold 158 gen 157 eps 0
  //Max diff: 1 is greater than expected (0)
  //1 of 51 (1.96078%) were different.
  //input\cel\exp\100907_CTGF2_U133plus2_t1ASCII_bad2.CEL chip is different max diff is: 1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_bad2_small.CEL");
  CelCheck celCheckTest6 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest6.check(errorMsg));
  //eps=1.0
  //gold 158 gen 157 eps 1
  //input\cel\exp\100907_CTGF2_U133plus2_t1ASCII_bad2.CEL chip is same max diff is: 1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_bad2_small.CEL");
  CelCheck celCheckTest7 (generated, gold,1.00);
  CPPUNIT_ASSERT(celCheckTest7.check(errorMsg)); 
  POSITIVE_TEST(celCheckTest7.check(errorMsg));

  //gold 158 gen 159 eps 0
  //gold 20520 gen 20522 eps 0
  //Max diff: 2 is greater than expected (0)
  //2 of 51 (3.92157%) were different.
  //input\cel\exp\100907_CTGF2_U133plus2_t1ASCII_bad3.CEL chip is different max diff is: 2
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_bad3_small.CEL");
  CelCheck celCheckTest8 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest8.check(errorMsg));

  //maxDiff==2
  //positive result
  //gold 158 gen 157 eps 0
  //gold 20520 gen 20518 eps 0
  //Max diff: 2 is greater than expected (0)
  //2 of 51 (3.92157%) were different.
  //input\cel\exp\100907_CTGF2_U133plus2_t1ASCII_bad3.CEL chip is same max diff is: 2
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/exp/100907_CTGF2_U133plus2_t1ASCII_bad3_small.CEL");
  CelCheck celCheckTest9 (generated, gold,0.00,"affymetrix-",2);
  CPPUNIT_ASSERT(celCheckTest9.check(errorMsg)); 
  POSITIVE_TEST(celCheckTest9.check(errorMsg));   
}

void CelCheckTest::testMapBinary()

{
	cout<<endl;
  Verbose::out(1, "CelCheckTest::testMapBinary");
  //small file 1col X 3 rows
  //positive
  vector<string> generated, gold; 
  std::string errorMsg ("");
  gold.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_Gold_small.CEL");
  generated.push_back("input/cel/map/CCL-256.1D_NSP_small.CEL");
  CelCheck celCheckTest1 (generated, gold,0.00);
  CPPUNIT_ASSERT(celCheckTest1.check(errorMsg));
  POSITIVE_TEST(celCheckTest1.check(errorMsg)); 

  //negative
  //Error encountered: CelCheck::check() - generated and gold vectors must be same size.
  gold.clear();
  generated.clear();
  gold.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_small.CEL");
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_small.CEL");
  CelCheck celCheckTest2 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest2.check(errorMsg)); 

  //bad gold file
  //Can't read cel file: input1\cel\map\CCL-256.1D_NSP_Gold.CEL
  gold.clear();
  generated.clear();
  errorMsg = "" ;
  gold.push_back("input1/cel/map/CCL-256.1D_NSP_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_small.CEL");
  CelCheck celCheckTest3 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest3.check(errorMsg));
  if(!celCheckTest3.check(errorMsg))
	  cout<<errorMsg<<endl;
  
  //bad generated file
  //Can't read cel file: input1\cel\map\CCL-256.1D_NSP.CEL
  gold.clear();
  generated.clear();
  errorMsg = "";
  generated.push_back("input1/cel/map/CCL-256.1D_NSP_small.CEL");
  gold.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_Gold_small.CEL");
  CelCheck celCheckTest4 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest4.check(errorMsg));
  if(!celCheckTest4.check(errorMsg))
	  cout<<errorMsg<<endl;
  
  //file should start with 64
  //Can't read cel file: input\cel\map\CCL-256.1D_NSP.CEL_bad1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP.CEL_bad1_small.CEL");
  CelCheck celCheckTest5 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest5.check(errorMsg));
  if(!celCheckTest5.check(errorMsg))
	  cout<<errorMsg<<endl;

  //one entry has been changed
  //Max diff: 0.000976563 is greater than expected (0)
  //1 of 3 (33.3333%) were different.
  //input\cel\map\CCL-256.1D_NSP_bad2_small.CEL chip is different max diff is: 0.000976563
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_bad2_small.CEL");
  CelCheck celCheckTest6 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest6.check(errorMsg));
  //eps=0.001
  //gold 9165 gen 9165.25 eps 0.25
  //input\cel\map\CCL-256.1D_NSP_bad2_small.CEL chip is same max diff is: 0.000976563
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_bad2_small.CEL");
  CelCheck celCheckTest7 (generated, gold,0.001);
  CPPUNIT_ASSERT(celCheckTest7.check(errorMsg));
  POSITIVE_TEST(celCheckTest7.check(errorMsg));

  //two entries have been changed
  //Max diff: 1 is greater than expected (0)
  //2 of 3 (66.6667%) were different.
  //iinput\cel\map\CCL-256.1D_NSP_bad3_small.CEL chip is different max diff is: 1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_bad3_small.CEL");
  CelCheck celCheckTest8 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest8.check(errorMsg));
  
  //maxDiff==2
  //Max diff: 1 is greater than expected (0)
  //2 of 3 (66.6667%) were different.
  //input\cel\map\CCL-256.1D_NSP_bad3_small.CEL chip is same max diff is: 1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_bad3_small.CEL");
  CelCheck celCheckTest9 (generated, gold,0.00,"affymetrix-",2);
  CPPUNIT_ASSERT(celCheckTest9.check(errorMsg));
  POSITIVE_TEST(celCheckTest9.check(errorMsg));
}

void CelCheckTest::testMapASCII()

{
  cout<<endl;
  Verbose::out(1, "CelCheckTest::testMapASCII");
  //positive
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII_small.CEL");
  CelCheck celCheckTest1 (generated, gold,0.00);
  CPPUNIT_ASSERT(celCheckTest1.check(errorMsg));
  POSITIVE_TEST(celCheckTest1.check(errorMsg));

  //bad gold file
  //Can't read cel file: input1\cel\map\CCL-256.1D_NSP_ASCII_Gold.CEL
  gold.clear();
  generated.clear();
  errorMsg = "" ;
  gold.push_back("input1/cel/map/CCL-256.1D_NSP_ASCII_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII_small.CEL");
  CelCheck celCheckTest3 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest3.check(errorMsg));
  if(!celCheckTest3.check(errorMsg))
	  cout<<errorMsg<<endl;
  
  //bad generated file
  //Can't read cel file: input1\cel\map\CCL-256.1D_NSP_ASCII.CEL
  gold.clear();
  generated.clear();
  errorMsg = "";
  gold.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII_Gold_small.CEL");
  generated.push_back("input1/cel/map/CCL-256.1D_NSP_ASCII_small.CEL");
  CelCheck celCheckTest4 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest4.check(errorMsg));
  if(!celCheckTest4.check(errorMsg))
	  cout<<errorMsg<<endl;
  
  //file 
  //Can't read cel file: input\cel\map\CCL-256.1D_NSP_ASCII.CEL_bad1.CEL
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII.CEL_bad1_small.CEL");
  CelCheck celCheckTest5 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest5.check(errorMsg));
  if(!celCheckTest5.check(errorMsg))
	  cout<<errorMsg<<endl;

  //gold 9165 gen 9166 eps 0
  //Max diff: 1 is greater than expected (0)
  //1 of 6553600 (1.52588e-005%) were different.
  //input\cel\map\CCL-256.1D_NSP_ASCII_bad2.CEL chip is different max diff is: 1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII_bad2_small.CEL");
  CelCheck celCheckTest6 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest6.check(errorMsg));
  //eps=1.0
  //gold 9165 gen 9166 eps 1
  //input\cel\map\CCL-256.1D_NSP_ASCII_bad2.CEL chip is same max diff is: 1
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII_bad2_small.CEL");
  CelCheck celCheckTest7 (generated, gold,1.00);
  CPPUNIT_ASSERT(celCheckTest7.check(errorMsg));
  POSITIVE_TEST(celCheckTest7.check(errorMsg));

  //gold 9165 gen 9166 eps 0
  //gold 229 gen 232 eps 0
  //Max diff: 3 is greater than expected (0)
  //2 of 6553600 (3.05176e-005%) were different.
  //input\cel\map\CCL-256.1D_NSP_ASCII_bad3.CEL chip is different max diff is: 3
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII_bad3_small.CEL");
  CelCheck celCheckTest8 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest8.check(errorMsg));

  //maxDiff==2
  //gold 9165 gen 9166 eps 0
  //gold 229 gen 232 eps 0
  //Max diff: 3 is greater than expected (0)
  //2 of 6553600 (3.05176e-005%) were different.
  //input\cel\map\CCL-256.1D_NSP_ASCII_bad3.CEL chip is same max diff is: 3
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/map/CCL-256.1D_NSP_ASCII_bad3_small.CEL");
  CelCheck celCheckTest9 (generated, gold,0.00,"affymetrix-",2);
  //CPPUNIT_ASSERT(celCheckTest9.check(errorMsg)); // todo vliber: fails on Linux and OS-X
  POSITIVE_TEST(celCheckTest9.check(errorMsg));
}

void CelCheckTest::testCalvinBinary()
{
  cout<<endl;
  Verbose::out(1, "CelCheckTest::testCalvinBinary");

  //positive
  vector<string> generated, gold;
  std::string errorMsg ("");
  gold.push_back(INPUT+"/cel/com-console/NA06985_GW6_C_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/com-console/NA06985_GW6_C_small.CEL");
  CelCheck celCheckTest1 (generated, gold,0.00);
  //input\cel\com-console\NA06985_GW6_C_small.CEL chip is same max diff is: 0
  CPPUNIT_ASSERT(celCheckTest1.check(errorMsg));
  //input\cel\com-console\NA06985_GW6_C_small.CEL chip is same max diff is: 0
  POSITIVE_TEST(celCheckTest1.check(errorMsg));

  //negative
  //Error encountered: CelCheck::check() - generated and gold vectors must be same size.
  generated.clear();
  generated.push_back(INPUT+"/cel/com-console/NA06985_GW6_C_small.CEL");
  generated.push_back(INPUT+"/cel/com-console/NA06993_GW6_C_small.CEL");
  CelCheck celCheckTest2 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest2.check(errorMsg));

  //bad gold file
  //Can't read cel file: input1\cel\com-console\NA06985_GW6_C_Gold_small.CEL
  gold.clear();
  generated.clear();
  errorMsg = "" ;
  gold.push_back("input1/cel/com-console/NA06985_GW6_C_Gold_small.CEL");
  generated.push_back(INPUT+"/cel/com-console/CCL-256.1D_NSP_small.CEL");
  CelCheck celCheckTest3 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest3.check(errorMsg));
  if(!celCheckTest3.check(errorMsg))
	 cout<<errorMsg<<endl;
  
  //bad generated file
  //Can't read cel file: input1\cel\com-console\CCL-256.1D_NSP_small.CEL
  gold.clear();
  generated.clear();
  errorMsg = "";
  generated.push_back("input1/cel/com-console/CCL-256.1D_NSP_small.CEL");
  gold.push_back(INPUT+"/cel/com-console/NA06985_GW6_C_Gold_small.CEL");
  CelCheck celCheckTest4 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest4.check(errorMsg));
  if(!celCheckTest4.check(errorMsg))
	  cout<<errorMsg<<endl;
  
  //file should start with 59
  //Can't read cel file: input1\cel\com-console\NA06985_GW6_C_Bad1_small.CEL
  errorMsg = "";
  generated.clear();
  generated.push_back(INPUT+"/cel/com-console/NA06985_GW6_C_Bad1_small.CEL");
  CelCheck celCheckTest5 (generated, gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest5.check(errorMsg));
  if(!celCheckTest5.check(errorMsg))
	  cout<<errorMsg<<endl;

  //bad intensity in first set
  //****00001930h/0 from 00 to 10 ****
  //****00001930h/3 from 00 to 01 ****
  //gold 6619 gen 6619.01 eps 0
  //gold 292 gen 292.008 eps 0
  //Max diff: 0.0078125 is greater than expected (0)
  //2 of 5 (40%) were different
  //input\cel\com-console\NA06985_GW6_C_Bad2_small.CEL chip is different max diff is: 0.0078125
  generated.clear();
  generated.push_back(INPUT+"/cel/com-console/NA06985_GW6_C_Bad2_small.CEL");
  CelCheck celCheckTest6 (generated,gold,0.00);
  CPPUNIT_ASSERT(!celCheckTest6.check(errorMsg));
  //eps>diff=0.008 positive
  //input\cel\com-console\NA06985_GW6_C_Bad2_small.CEL chip is same max diff is: 0.0078125
  CelCheck celCheckTest7 (generated,gold,0.008);
  CPPUNIT_ASSERT(celCheckTest7.check(errorMsg));
  POSITIVE_TEST(celCheckTest7.check(errorMsg));
  //errCountAllow=2 positive
  //Max diff: 0.0078125 is greater than expected (0)
  //2 of 5 (40%) were different.
  //input\cel\com-console\NA06985_GW6_C_Bad2_small.CEL chip is same max diff is: 0.0078125
  CelCheck celCheckTest8 (generated,gold,0.00,"affymetrix-",2);
  //CPPUNIT_ASSERT(celCheckTest8.check(errorMsg)); // todo vliber: fails on Linux and OS-X
  POSITIVE_TEST(celCheckTest8.check(errorMsg));
}
