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

// You should have received a copy of the GNU General Public License
// along with this program;if not, write to the
// Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

////////////////////////////////////////////////////////////////

/**

 * @file   CalvinEquivalentTest.cpp
 * @author vliber
 * last change by rsatin on 10/13/09
 *
 @brief  Testing the Calvin::equivalent() function.
 *
 */


#include "util/AffxString.h"
#include "util/Convert.h"
#include "util/Util.h"
#include "util/Verbose.h"
#include "calvin_files/utils/src/Calvin.h"
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

class CalvinEquivalentTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( CalvinEquivalentTest );
  CPPUNIT_TEST(testEquivalentPositive);
  CPPUNIT_TEST(testEquivalentFileHeader);
  CPPUNIT_TEST(testEquivalentGenericDataHeader);
  CPPUNIT_TEST(testEquivalentDataGroup_oneDataSet);
  CPPUNIT_TEST(testEquivalentDataGroup_multipleDataSet);
  CPPUNIT_TEST(testEquivalentFileExistence);
  CPPUNIT_TEST(testHeaderPresentAbsentIgnore);
  CPPUNIT_TEST(testEquivalentFracDiff);
  CPPUNIT_TEST(testEquivalentHeaderFracDiff);
  CPPUNIT_TEST_SUITE_END();

private:
  /** utility functions */
  int message_error_count(string message);

public: 
  virtual void tearDown();
  /** unit test functions */
  void testEquivalentPositive();
  void testEquivalentFileHeader();
  void testEquivalentGenericDataHeader();
  void testEquivalentDataGroup_oneDataSet();
  void testEquivalentDataGroup_multipleDataSet();
  void testEquivalentFileExistence();
  void testHeaderPresentAbsentIgnore();
  void testEquivalentFracDiff();
  void testEquivalentHeaderFracDiff();
  
};
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(CalvinEquivalentTest );

void CalvinEquivalentTest::tearDown () 
{
  // clean up after failed test cases (pops any dead objects from stack)
  while( Verbose::getParam().m_MsgHandler.size() > 1 )
    Verbose::popMsgHandler();
}

int CalvinEquivalentTest::message_error_count( string message ) {
  size_t fpos = 0;
  size_t fpos_last = 0;
  // find output from last test case and end of log
  while( fpos != string::npos ) {
    fpos_last = fpos;
    fpos = message.find("Comparing ",fpos+1);
  }
  fpos = fpos_last;
  // add up number of differences reported
  long error_count = 0;
  while( fpos != string::npos ) {
    fpos = message.find("NumberFailures = ",fpos);
	if( fpos != string::npos ) {
      char *endptr;
      int radix = 10;
	  fpos = fpos + strlen("NumberFailures = ");
      std::string strCount = message.substr( fpos, message.length()-fpos );
      error_count = error_count + strtol( strCount.c_str(), &endptr, radix ); 
      }
    }
  return error_count;
}

void CalvinEquivalentTest::testEquivalentPositive()
{
  std::cout<<std::endl;
  std::cout<<std::endl;
  Verbose::out(1, "***CalvinEquivalentTest testcases***");
  Verbose::out(1, "**CalvinEquivalentTest::testEquivalent**");
  std::set<std::string> setIgnore;
  std::map<std::string, float> mapEpsilon;

  //positive genotyping data
  CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

  //positive copynumber data
  CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_small.CNCHP", INPUT+"/calvin/copynumber/NA06985_GW6_C.CN5_small.CNCHP", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true)==1);
  CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var.chp", INPUT+"/calvin/copynumber/GenomeWideSNP_6.agcc.multidata_copynumber_var.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //expression data
  CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/expr/3AJW02021805.default.rma.chp", INPUT+"/calvin/expr/3AJW02021805.default.rma.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //dmet data  
  CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
}  

void CalvinEquivalentTest::testEquivalentFileHeader()
{
  std::cout<<std::endl; 
  Verbose::out(1, "**CalvinEquivalentTest::testEquivalentFileHeader**");
  std::set<std::string> setIgnore;
  std::map<std::string, float> mapEpsilon;
  
  //first character should be 59(Magic) todo vliber: add message
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad1_small.chp
  //affymetrix_calvin_exceptions::CalvinException caught.        Affymetrix GeneChip Command Console library has thrown an exception of type class affymetrix_calvin_exceptions::CalvinException Message is: 
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad1_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //version value is not correct(gold=1, gen=3) todo vliber: add message
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad2_small.chp
  //affymetrix_calvin_exceptions::CalvinException caught.        Affymetrix GeneChip Command Console library has thrown an exception of type class affymetrix_calvin_exceptions::CalvinException Message is: 
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad2_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //delete one parameter from gen file
  //>Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad2a_small.chp
  //WARNING: File Header does not have the same number of Parameters.
  //File Header Parameter name not found: program-name
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad2a_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //delete one parameter from gold file
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad2a_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp
  //WARNING: File Header does not have the same number of Parameters.
  //File Header Parameter name not found: 
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad2a_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
 
  //data groups count is not correct
  //original=1, change to 2
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3_small.chp
  //Expecting result: Files do not have the same number of Data Groups. Gold = 1, generated = 2
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   
  }

void CalvinEquivalentTest::testEquivalentGenericDataHeader()
{
  std::cout<<std::endl;
  Verbose::out(1, "**CalvinEquivalentTest::testEquivalentGenericDataHeader**");
  std::set<std::string> setIgnore;
  std::map<std::string, float> mapEpsilon;

  //00000010h/9   6D to 7D
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad4_small.chp
  //File Header Parameter Value mismatch for FileTypeIdentifier affymetrix-multi-data-type-analysis != affymetrix-}ulti-data-type-analysis
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad4_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //00000030h/f   2D to 3D
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad4a_small.chp
  //File Header Parameter Value mismatch for FileIdentifier 0000010154-1197306474-1669795230-1795436647-0174193476 != 0000010154=1197306474-1669795230-1795436647-0174193476
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad4a_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //00000070h/8   2D to 3D
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad4b_small.chp
  //File Header Parameter Value mismatch for FileLocale en-US != en=US
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad4b_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
 
  
  //00000160h/a 79 to 78
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad5a_small.chp
  //File Header Parameter name not found: affymetrix-array-type
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad5a_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

  //00000240h/8 65 to 66
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad5b_small.chp
  //File Header Parameter Type mismatch for parameter: affymetrix-array-type
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad5b_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

  //00000180h/8 69 to 79
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad5_small.chp
  //File Header Parameter Value mismatch for affymetrix-array-type GenomeWideSNP_6 != GenomeWydeSNP_6
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad5_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
	
  //new file
  //for field 'het_rate' expecting: '26.744295' got: '6846.539551'  **** row 000045c0h/1 from 41 to 45
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad8a_small.chp
  //File Header Parameter Value is out of spec. for affymetrix-chipsummary-het_rate Difference = 6819.7954101563 Epsilon = 0.0001
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8a_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  //new file
  setIgnore.insert("affymetrix-chipsummary-het_rate");
  //for field 'em-cluster-chrX-het-contrast_gender' expecting: 'female' got: 'gemale'   ****  row 00004d10h/9 from 66 to 67
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad8a_small.chp
  //File Header Parameter Value mismatch for affymetrix-chipsummary-em-cluster-chrX-het-contrast_gender female != gemale 
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8a_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  //new file
  setIgnore.insert("affymetrix-chipsummary-em-cluster-chrX-het-contrast_gender");
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad8a_small.chp
  //Files are equivalent.
  CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8a_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //00004e00h/c  01 to 02
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad8b_small.chp
  //Files do not have the same number of Parent Headers.
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8b_small.chp", setIgnore, mapEpsilon, 0.00f, 1.0, true));
  
  setIgnore.clear();
  //000045c0h/3 from F4 to A4
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad8c_small.chp
  //File Header Parameter Value is out of spec. for affymetrix-chipsummary-het_rate Difference = 0.0390625 Epsilon = 0
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8c_small.chp", setIgnore, mapEpsilon, 0.00f, 1.0, true));
  //Files are equivalent
  CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8c_small.chp", setIgnore, mapEpsilon, 0.04f, 1.0, true));
  mapEpsilon.clear();
  //File Header Parameter Value is out of spec. for affymetrix-chipsummary-het_rate Difference = 0.0390625 Epsilon = 0.0099999997
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8c_small.chp", setIgnore, mapEpsilon, 0.00f, 1.0, true));
  mapEpsilon["affymetrix-chipsummary-het_rate"]=0.04f;
  mapEpsilon["affymetrix-chipsummary-hom_rate"]=0.00f;
  //File Header Parameter Value is out of spec. for affymetrix-chipsummary-hom_rate Difference = 0.01171875 Epsilon = 0
  CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8c_small.chp", setIgnore, mapEpsilon, 0.00f, 1.0, true));
  mapEpsilon["affymetrix-chipsummary-hom_rate"]=0.012f;
  //Files are equivalent
  CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8c_small.chp", setIgnore, mapEpsilon, 0.00f, 1.0, true));
  
  }

void CalvinEquivalentTest::testEquivalentDataGroup_oneDataSet()
{
  std::cout<<std::endl; 
  Verbose::out(1, "**CalvinEquivalentTest::testEquivalentDataGroup_oneDataSet**");
  std::set<std::string> setIgnore;
  std::map<std::string, float> mapEpsilon;
  

   Verbose::out(1, "**Data Group name comparison**");
  //Data Group name mismatch 000066d0h/2 4D to AD
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3a_small.chp
  //Expecting result: File Data Group name mismatch: MultiData != ¡ultiData
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3a_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
   //Data Sets count is not correct 000066c0h/c 01 to 02
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3c_small.chp
   //Expecting result: Data Group MultiData does not have the same number of Data Sets.
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3c_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
   Verbose::out(1, "**Data Set comparison**");
   //Data set name mismatch 000066f0h/0 47 to 57
  //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3b_small.chp
  //Expecting result: File Data Set name mismatch: MultiData.Genotype != MultiData.Wenotype
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3b_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

   //number of Columns is not the same 00006700h/6 06 to 05 
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3d_small.chp
   //Expecting result: Data Set MultiData.Genotype does not have the same number of Columns. Gold = 6, Generated = 5
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3d_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
   //number of Columns is not the same 000067a0h/a 96 to 97 
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3e_small.chp
   //Data Set MultiData.Genotype does not have the same number of Rows. Gold = 150, Generated = 151
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3e_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   
   //00006790h/d   61  to 71
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3f_small.chp
   //File Data Set column name mismatch: MultiData.Genotype.Forced Call != MultiData.Genotype.Forced Cqll
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3f_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   
   //number of Columns is not the same 00006790h/d 61 to 71 
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3g_small.chp
   //Value Type mismatch for Data Set Column MultiData.Genotype.Forced Call, Gold = 1, Generated = 2
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3g_small.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
   Verbose::out(1, "**comparison all data for one row **");
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3h_small.chp
   //Value mismatch for Data Set Column MultiData.Genotype.ProbeSetName at row 1 //000067dh0/8 31 to 32
   //Value mismatch for Data Set Column MultiData.Genotype.Call at row 1 //000067eh0/3 07 to 06
   //000067eh0/5 14 to 54
   //Value is out of spec. for Data Set Column MultiData.Genotype.Confidence MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999794286
   //000067eh0/b 1B to 6B
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal A MaxDifference at row 1 MaxDifference = 0.0024414063 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //000067eh0/f EF to EE
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.0469360352 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value mismatch for Data Set Column MultiData.Genotype.Forced Call at row 1 //000067fh0/0 07 to 05
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3h_small.chp", setIgnore, mapEpsilon, 0.00f, 1.0, true));
   
   Verbose::out(1, "**use epsilon value for comparison**");
   //Value is out of spec. for Data Set Column MultiData.Genotype.Confidence MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999794286
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal A MaxDifference at row 1 MaxDifference = 0.0024414063 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.046875 NumberFailures = 1 Epsilon = 0 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.000f, 1.0, true));
   //Value is out of spec. for Data Set Column MultiData.Genotype.Confidence MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0.0024999999 Correlation = 0.9999794286
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.046875 NumberFailures = 1 Epsilon = 0.0024999999 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.0025f, 1.0, true));
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.046875 NumberFailures = 1 Epsilon = 0.0040000002 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.004f, 1.0, true));
   //Files are equivalent
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.047f, 1.0, true));
   
  
   Verbose::out(1, "**use setIgnore to avoid value comparison**");
   //Value is out of spec. for Data Set Column MultiData.Genotype.Confidence MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999794286
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal A MaxDifference at row 1 MaxDifference = 0.0024414063 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.046875 NumberFailures = 1 Epsilon = 0 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.000f, 1.0, true));
   setIgnore.clear();
   setIgnore.insert("MultiData.Genotype.Confidence");
   setIgnore.insert("MultiData.Genotype.Signal A");
   setIgnore.insert("MultiData.Genotype.Signal B");
   //Files are equivalent
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.000f, 1.0, true));
   setIgnore.clear();
   
   
   Verbose::out(1, "**use map epsilon value for comparison**");
   //Value is out of spec. for Data Set Column MultiData.Genotype.Confidence MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999794286
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal A MaxDifference at row 1 MaxDifference = 0.0024414063 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.046875 NumberFailures = 1 Epsilon = 0 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.000f, 1.0, true));
   mapEpsilon["MultiData.Genotype.Confidence"]=0.004f;
   mapEpsilon["MultiData.Genotype.Signal A"]=0.0024f;
   mapEpsilon["MultiData.Genotype.Signal B"]=0.046f;
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal A MaxDifference at row 1 MaxDifference = 0.0024414063 NumberFailures = 1 Epsilon = 0.0024000001 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.046875 NumberFailures = 1 Epsilon = 0.0460000001 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.000f, 1.0, true));
   mapEpsilon["MultiData.Genotype.Signal A"]=0.00245f;
   mapEpsilon["MultiData.Genotype.Signal B"]=0.046f;
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.046875 NumberFailures = 1 Epsilon = 0.0460000001 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.000f, 1.0, true));
   mapEpsilon["MultiData.Genotype.Signal B"]=0.0469f;
   //Files are equivalent
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.000f, 1.0, true));
   mapEpsilon.clear();

  
   Verbose::out(1, "**correlation**");
   //Value is out of spec. for Data Set Column MultiData.Genotype.Confidence MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0.0024999999 Correlation = 0.9999794286
   //Value is out of spec. for Data Set Column MultiData.Genotype.Signal B MaxDifference at row 1 MaxDifference = 0.046875 NumberFailures = 1 Epsilon = 0.0024999999 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.0025f, 1.0, true));
   //Value is out of spec. for Data Set Column MultiData.Genotype.Confidence MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999794286
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.00f, 0.99999, true));
   //Files are equivalent
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3i_small.chp", setIgnore, mapEpsilon, 0.00f, 0.9998, true));
  
   Verbose::out(1, "**message limit**");
   //000067b0h/0 4E to 5E
   //000067c0h/0 08 to 07
   //000067c0h/d 08 to 07
   //000067d0h/2 53 to 54
   //000067f0h/5 53 to 54
   //00006810h/8 53 to 54
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3j_small.chp
   //Value mismatch for Data Set Column MultiData.Genotype.ProbeSetName at row 0
   //Value mismatch for Data Set Column MultiData.Genotype.ProbeSetName at row 1
   //Value mismatch for Data Set Column MultiData.Genotype.ProbeSetName at row 2
   //Value mismatch for Data Set Column MultiData.Genotype.ProbeSetName at row 3
   //Value mismatch for Data Set Column MultiData.Genotype.Call at row 0
   //Value mismatch for Data Set Column MultiData.Genotype.Forced Call at row 0
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3j_small.chp", setIgnore, mapEpsilon, 0.00f, 1.0, true));
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3j_small.chp
   //Value mismatch for Data Set Column MultiData.Genotype.ProbeSetName at row 0
   //Value mismatch for Data Set Column MultiData.Genotype.ProbeSetName at row 1
   //Message limit reached. Additional messages will be suppressed.
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp", INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad3j_small.chp", setIgnore, mapEpsilon, 0.00f, 1.0, true,2));
   
  }

void CalvinEquivalentTest::testEquivalentDataGroup_multipleDataSet()
{
  std::cout<<std::endl; 
  Verbose::out(1, "CalvinEquivalentTest::testEquivalentDataGroup_multipleDataSet");
  std::set<std::string> setIgnore;
  std::map<std::string, float> mapEpsilon;

  Verbose::out(1, "**Data Group name comparison**");
   //Data Group name mismatch 00007b10h/7 4D to AD
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_30.chp
   //File Data Group name mismatch: MultiData != ¡ultiData
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_30.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
   //Data Sets count is not correct 00007b10h/1 03 to 04
   //Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad3c_small.chp
   //Expecting result: Data Group MultiData does not have the same number of Data Sets.
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_31.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
   Verbose::out(1, "**Data Set comparison**");
   //Data set name mismatch 00007b40h/1 41 to 40
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_32.chp
   //File Data Set name mismatch: MultiData.DmetBiAllelic != MultiData.DmetBi@llelic
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_32a.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //Data set name mismatch 00015b20h/2/2 44 to 43
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_32a.chp
   //File Data Set name mismatch: MultiData.DmetMultiAllelic != MultiData.CmetMultiAllelic
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_32b.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //Data set name mismatch 000162e0h/9 44 to 43
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_32c.chp
   //File Data Set name mismatch: MultiData.DmetCopyNumber != MultiData.CmetCopyNumber
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_32c.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

   //00007b50h/5 08 to 07 
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_33a.chp
   //Data Set MultiData.DmetBiAllelic does not have the same number of Columns. Gold = 8, Generated = 7
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_33a.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //00015b40h/8 11 to 10
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_33b.chp
   //Data Set MultiData.DmetMultiAllelic does not have the same number of Columns. Gold = 17, Generated = 16 
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_33b.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //00016300h/b 07 to 06
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_33c.chp
   //Data Set MultiData.DmetCopyNumber does not have the same number of Columns. Gold = 7, Generated = 6
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_33c.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

   //00007c20h/f 6E to 5E 
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_34a.chp
   //Data Set MultiData.DmetBiAllelic does not have the same number of Rows. Gold = 1902, Generated = 1886
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_34a.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //00015d10h/3 1D to 2D 
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_34b.chp
   //Data Set MultiData.DmetMultiAllelic does not have the same number of Rows. Gold = 29, Generated = 45
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_34b.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //000163c0h/e 05 to 04 
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_34c.chp
   //Data Set MultiData.DmetCopyNumber does not have the same number of Rows. Gold = 5, Generated = 4
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_34c.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

   //00007bb0h/0 63 to 64
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_35a.chp
   //File Data Set column name mismatch: MultiData.DmetBiAllelic.Forced Call != MultiData.DmetBiAllelic.Forded Call
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_35a.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //00015c30h/0 61 to 62
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_35b.chp
   //File Data Set column name mismatch: MultiData.DmetMultiAllelic.Signal D != MultiData.DmetMultiAllelic.Signbl D
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_35b.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //000163a0h/0 4E to 4C
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet_35c.chp
   //File Data Set column name mismatch: MultiData.DmetCopyNumber.CN_force != MultiData.DmetCopyNumber.CNOforce
   //File Data Set column name mismatch: MultiData.DmetCopyNumber.CN_lower != MultiData.DmetCopyNumber.CL_lower
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_35c.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

   //00007c20h/7 01 to 02
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1.dmet_36a.chp
   //Value Type mismatch for Data Set Column MultiData.DmetBiAllelic.Context B, Gold = 1, Generated = 2
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_36a.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //00015cf0h/0 01 to 02
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1.dmet_36b.chp
   //Value Type mismatch for Data Set Column MultiData.DmetMultiAllelic.Context E, Gold = 1, Generated = 2
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_36b.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   //000163c0h/6 06 to 07
   //Comparing ./input/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp and ./input/calvin/dmet/ZQ_D10_gCtrl_1.dmet_36c.chp
   //Value Type mismatch for Data Set Column MultiData.DmetCopyNumber.CN_upper, Gold = 6, Generated = 7
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_36c.chp", setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
   mapEpsilon.clear();
   setIgnore.clear();
   
   Verbose::out(1, "**use epsilon value for comparison**");
   //00007c40h/6 00 to 90
   //00007c60h/8 00 to 02
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal A MaxDifference at row 0 MaxDifference = 0.0087890625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal B MaxDifference at row 1 MaxDifference = 0.0001220703 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //00015d20h/5 00 to 04
   //00015d20h/f 04 to 14
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Confidence MaxDifference at row 0 MaxDifference = 0.0002441406 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999999993
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Signal A MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //00016450h/5 00 to 02
   //00016460h/f 00 to 01
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_estim MaxDifference at row 4 MaxDifference = 0.0009765625 NumberFailures = 1 Epsilon = 0 Correlation = nan
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_upper MaxDifference at row 3 MaxDifference = 0.0001220703 NumberFailures = 1 Epsilon = 0 Correlation = nan
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal A MaxDifference at row 0 MaxDifference = 0.0087890625 NumberFailures = 1 Epsilon = 0.00013 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Confidence MaxDifference at row 0 MaxDifference = 0.0002441406 NumberFailures = 1 Epsilon = 0.00013 Correlation = 0.9999999993
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Signal A MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0.00013 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_estim MaxDifference at row 4 MaxDifference = 0.0009765625 NumberFailures = 1 Epsilon = 0.00013 Correlation = nan
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.00013f, 1.0, true));
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal A MaxDifference at row 0 MaxDifference = 0.0087890625 NumberFailures = 1 Epsilon = 0.00025 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Signal A MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0.00025 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_estim MaxDifference at row 4 MaxDifference = 0.0009765625 NumberFailures = 1 Epsilon = 0.00025 Correlation = nan
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.00025f, 1.0, true));
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal A MaxDifference at row 0 MaxDifference = 0.0087890625 NumberFailures = 1 Epsilon = 0.0040000002 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.004f, 1.0, true));
   //Files are equivalent.
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.009f, 1.0, true));
      
  
   Verbose::out(1, "**use setIgnore to avoid value comparison**");
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal A MaxDifference at row 0 MaxDifference = 0.0087890625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal B MaxDifference at row 1 MaxDifference = 0.0001220703 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Confidence MaxDifference at row 0 MaxDifference = 0.0002441406 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999999993
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Signal A MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_estim MaxDifference at row 4 MaxDifference = 0.0009765625 NumberFailures = 1 Epsilon = 0 Correlation = nan
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_upper MaxDifference at row 3 MaxDifference = 0.0001220703 NumberFailures = 1 Epsilon = 0 Correlation = nan
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   setIgnore.clear();
   setIgnore.insert("MultiData.DmetBiAllelic.Signal A");
   setIgnore.insert("MultiData.DmetBiAllelic.Signal B");
   setIgnore.insert("MultiData.DmetMultiAllelic.Confidence");
   setIgnore.insert("MultiData.DmetMultiAllelic.Signal A");
   setIgnore.insert("MultiData.DmetCopyNumber.CN_estim");
   setIgnore.insert("MultiData.DmetCopyNumber.CN_upper");
   //Files are equivalent.
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   setIgnore.clear();
   
   
   Verbose::out(1, "**use map epsilon value for comparison**");
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal A MaxDifference at row 0 MaxDifference = 0.0087890625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal B MaxDifference at row 1 MaxDifference = 0.0001220703 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Confidence MaxDifference at row 0 MaxDifference = 0.0002441406 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999999993
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Signal A MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_estim MaxDifference at row 4 MaxDifference = 0.0009765625 NumberFailures = 1 Epsilon = 0 Correlation = nan
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_upper MaxDifference at row 3 MaxDifference = 0.0001220703 NumberFailures = 1 Epsilon = 0 Correlation = nan
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   mapEpsilon["MultiData.DmetBiAllelic.Signal A"]=0.009f;
   mapEpsilon["MultiData.DmetBiAllelic.Signal B"]=0.00001f;
   mapEpsilon["MultiData.DmetMultiAllelic.Confidence"]=0.0001f;
   mapEpsilon["MultiData.DmetMultiAllelic.Signal A"]=0.0f;
   mapEpsilon["MultiData.DmetCopyNumber.CN_estim"]=0.0f;
   mapEpsilon["MultiData.DmetCopyNumber.CN_upper"]=0.00013f;
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal B MaxDifference at row 1 MaxDifference = 0.0001220703 NumberFailures = 1 Epsilon = 1e-005 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Confidence MaxDifference at row 0 MaxDifference = 0.0002441406 NumberFailures = 1 Epsilon = 0.0001 Correlation = 0.9999999993
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Signal A MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_estim MaxDifference at row 4 MaxDifference = 0.0009765625 NumberFailures = 1 Epsilon = 0 Correlation = nan
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   mapEpsilon["MultiData.DmetBiAllelic.Signal B"]=0.00013f;
   mapEpsilon["MultiData.DmetMultiAllelic.Signal A"]=0.004f;
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Confidence MaxDifference at row 0 MaxDifference = 0.0002441406 NumberFailures = 1 Epsilon = 0.0001 Correlation = 0.9999999993
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_estim MaxDifference at row 4 MaxDifference = 0.0009765625 NumberFailures = 1 Epsilon = 0 Correlation = nan
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   mapEpsilon["MultiData.DmetMultiAllelic.Confidence"]=0.00025f;
   mapEpsilon["MultiData.DmetCopyNumber.CN_estim"]=0.001f;
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Signal A MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Files are equivalent.
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   mapEpsilon.clear();

  
   Verbose::out(1, "**correlation**");
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal A MaxDifference at row 0 MaxDifference = 0.0087890625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetBiAllelic.Signal B MaxDifference at row 1 MaxDifference = 0.0001220703 NumberFailures = 1 Epsilon = 0 Correlation = 1
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Confidence MaxDifference at row 0 MaxDifference = 0.0002441406 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999999993
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Signal A MaxDifference at row 1 MaxDifference = 0.00390625 NumberFailures = 1 Epsilon = 0 Correlation = 1
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   //Value is out of spec. for Data Set Column MultiData.DmetMultiAllelic.Confidence MaxDifference at row 0 MaxDifference = 0.0002441406 NumberFailures = 1 Epsilon = 0 Correlation = 0.9999999993
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 0.99999999999, true));
   //Files are equivalent.
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_37.chp", setIgnore, mapEpsilon, 0.0000f, 0.99, true));
   Verbose::out(1, "**correlation only for deta set #3**");
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN Call MaxDifference at row 2 MaxDifference = 1 NumberFailures = 3 Epsilon = 0 Correlation = 0.4082482905
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN_force MaxDifference at row 1 MaxDifference = 2 NumberFailures = 4 Epsilon = 0 Correlation = 0.7453559925
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_39_1.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_39_2.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   //Value is out of spec. for Data Set Column MultiData.DmetCopyNumber.CN Call MaxDifference at row 2 MaxDifference = 1 NumberFailures = 3 Epsilon = 0 Correlation = 0.4082482905
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_39_1.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_39_2.chp", setIgnore, mapEpsilon, 0.0000f, 0.74, true));
   //Files are equivalent.
   CPPUNIT_ASSERT(Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_39_1.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_39_2.chp", setIgnore, mapEpsilon, 0.0000f, 0.40, true));
   
   Verbose::out(1, "**message limit**");
   //Value mismatch for Data Set Column MultiData.DmetBiAllelic.ProbeSetName at row 0
   //Value mismatch for Data Set Column MultiData.DmetBiAllelic.ProbeSetName at row 3
   //Value mismatch for Data Set Column MultiData.DmetBiAllelic.ProbeSetName at row 36
   //Value mismatch for Data Set Column MultiData.DmetMultiAllelic.ProbeSetName at row 20
   //Value mismatch for Data Set Column MultiData.DmetMultiAllelic.ProbeSetName at row 22
   //Value mismatch for Data Set Column MultiData.DmetMultiAllelic.ProbeSetName at row 27
   //Value mismatch for Data Set Column MultiData.DmetCopyNumber.ProbeSetName at row 1
   //Value mismatch for Data Set Column MultiData.DmetCopyNumber.ProbeSetName at row 2
   //Value mismatch for Data Set Column MultiData.DmetCopyNumber.ProbeSetName at row 3
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_38.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true));
   //Value mismatch for Data Set Column MultiData.DmetBiAllelic.ProbeSetName at row 0
   //Value mismatch for Data Set Column MultiData.DmetBiAllelic.ProbeSetName at row 3
   //Message limit reached. Additional messages will be suppressed.
   //Value mismatch for Data Set Column MultiData.DmetMultiAllelic.ProbeSetName at row 20
   //Value mismatch for Data Set Column MultiData.DmetMultiAllelic.ProbeSetName at row 22
   //Message limit reached. Additional messages will be suppressed.
   //Value mismatch for Data Set Column MultiData.DmetCopyNumber.ProbeSetName at row 1
   //Value mismatch for Data Set Column MultiData.DmetCopyNumber.ProbeSetName at row 2
   //Message limit reached. Additional messages will be suppressed.
   CPPUNIT_ASSERT(!Calvin::equivalent(INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1_Gold.dmet.chp", INPUT+"/calvin/dmet/ZQ_D10_gCtrl_1.dmet_38.chp", setIgnore, mapEpsilon, 0.0000f, 1.0, true,2));
 
  }

void CalvinEquivalentTest::testEquivalentFileExistence() {

  cout<<endl;
  Verbose::out(1, "**CalvinEquivalentTest::testEquivalentFileExistence**");

  std::map<std::string, float> mapEpsilon;

  std::set<std::string> setIgnore;
  setIgnore.insert("calvin-file");
   
  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  string firstFileName = INPUT+"/calvin/map/no_such_file1.chp";
  string secondFileName = INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";

  Verbose::out(1, "**testcase1-first-file-does-not-exist**");
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001F, 0.99999F, true, 1000, 0.0001F));

  // verify "File not found" error message is reported followed by filename causing error
  size_t fpos_fnf = MessageHandler1_oss.str().find("File not found:");
  CPPUNIT_ASSERT( fpos_fnf != string::npos );
  CPPUNIT_ASSERT( MessageHandler1_oss.str().find("no_such_file1.chp",fpos_fnf) != string::npos );

  Verbose::popMsgHandler();

  // message output streamed to string object
  ostringstream MessageHandler2_oss;
  MsgStream msgHandler2(2, &MessageHandler2_oss);
  Verbose::pushMsgHandler(&msgHandler2);

  firstFileName = INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp";
  secondFileName = INPUT+"/calvin/map/no_such_file2.chp";

  // test case: second file does not exist
  Verbose::out(1, "**testcase2-second-file-does-not-exist**");
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001F, 0.99999F, true, 1000, 0.0001F));

  // verify "File not found" error message is reported followed by filename causing error
  fpos_fnf = MessageHandler2_oss.str().find("File not found:");
  CPPUNIT_ASSERT( fpos_fnf != string::npos );
  CPPUNIT_ASSERT( MessageHandler2_oss.str().find("no_such_file2.chp",fpos_fnf) != string::npos );

  Verbose::popMsgHandler();
}

void CalvinEquivalentTest::testHeaderPresentAbsentIgnore()
{
  std::cout<<std::endl; 
  Verbose::out(1, "**CalvinEquivalentTest::testHeaderPresentAbsentIgnore**");

  string firstFileName = INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp";
  string secondFileName = INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_HdrParamPresAbs.chp";

  // first file has parameter "affymetrix-algorithm-param-apt-program-name"
  // second file has parameter "affymetrix-algorithm-param-apt-program-xxxx"

  std::map<std::string, float> mapEpsilon;
  std::set<std::string> setIgnore;

  //NEGATIVE CASE: one parameter missing from gen file, different parameter missing from gold file (not ignored)
  //>Comparing ./input/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp and ./input/calvin/map/NA06985_GW6_C.birdseed-v2_bad2a_small.chp
  //WARNING: File Header does not have the same number of Parameters.
  //File Header Parameter name not found: program-name
  Verbose::out(1, "**testcase1 - deleted parameters not ignored**");
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  // add first file parameter (missing from second file) to ignore list 
  setIgnore.insert("affymetrix-algorithm-param-apt-program-name");

  Verbose::out(1, "**testcase2a - deleted gold file parameter in ignore list**");
  CPPUNIT_ASSERT(Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));

  Verbose::out(1, "**testcase2b - deleted generated file parameter ignore list**");
  CPPUNIT_ASSERT(!Calvin::equivalent(secondFileName, firstFileName, setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  // add second file parameter (missing from first file) to ignore list 
  setIgnore.insert("affymetrix-algorithm-param-apt-program-xxxx");

   //POSITIVE CASE: delete one parameter from gen file (ignored)
  Verbose::out(1, "**testcase3a - deleted gen and gold file parameters in ignored list**");
  CPPUNIT_ASSERT(Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
  
  //POSITIVE CASE: delete one parameter from gold file (ignored)
  Verbose::out(1, "**testcase3b - deleted gen and gold file parameters in ignored list**");
  CPPUNIT_ASSERT(Calvin::equivalent(secondFileName, firstFileName, setIgnore, mapEpsilon, 0.0001f, 0.9999989, true));
 
  return;
}
void CalvinEquivalentTest::testEquivalentFracDiff()
{
  cout<<endl;
  Verbose::out(1, "**CalvinEquivalentTest::testEquivalentFracDiff**");

  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  std::map<std::string, float> mapEpsilon;

  std::set<std::string> setIgnore;
  setIgnore.insert("calvin-file");
  
  string firstFileName = INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_ValueDiff1_small.chp";
  string secondFileName = INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_ValueDiff2_small.chp";

  // test case: absent vs present optional argument frac==eps
  Verbose::out(1, "**testcase1-old behaviour**");    // 2 differences
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.001F));
  CPPUNIT_ASSERT( message_error_count(MessageHandler1_oss.str()) == 2 );
  Verbose::out(1, "**testcase1a-new behaviour**");   // 0 differences
  CPPUNIT_ASSERT(Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.001F, 1.00001, true, 1000, 0.001F));

  // test case: absent vs present optional argument frac>eps
  Verbose::out(1, "**testcase2-old behaviour**");    // 19 differences
  //CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001F)); //TODO rsatin default dCorrelationCutoff ignores some differences
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001F, 1.00001));
  CPPUNIT_ASSERT( message_error_count(MessageHandler1_oss.str()) == 19 );
  Verbose::out(1, "**testcase2a-new behaviour**");   // 0 differences
  CPPUNIT_ASSERT(Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001F, 1.00001, true, 1000, 0.001F));

  // test case: absent vs present optional argument frac<eps
  Verbose::out(1, "**testcase3-old behaviour**");    // 19 differences
  //CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001F)); //TODO rsatin default dCorrelationCutoff ignores some differences
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001F, 1.00001));
  CPPUNIT_ASSERT( message_error_count(MessageHandler1_oss.str()) == 19 );
  Verbose::out(1, "**testcase3a-new behaviour**");   // 6 differences
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon, 0.0001F, 1.00001, true, 1000, 0.00001F));
  CPPUNIT_ASSERT( message_error_count(MessageHandler1_oss.str()) == 6 );

  Verbose::popMsgHandler();

  return;
}

void CalvinEquivalentTest::testEquivalentHeaderFracDiff()
{
  cout<<endl;
  Verbose::out(1, "CalvinEquivalentTest::testEquivalentHeaderFracDiff");

  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  std::map<std::string, float> mapEpsilon_empty;

  std::map<std::string, float> mapEpsilon_full;
  mapEpsilon_full["affymetrix-chipsummary-het_rate"]=0.0401f;
  mapEpsilon_full["affymetrix-chipsummary-hom_rate"]=0.0101f;

  std::set<std::string> setIgnore;
  setIgnore.insert("calvin-file");
  
  string firstFileName = INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_Gold_small.chp";
  string secondFileName = INPUT+"/calvin/map/NA06985_GW6_C.birdseed-v2_bad8c_small.chp";

  // test case: absent vs present optional argument frac==eps applied to header parameter differences

  Verbose::out(1, "**testcase1-old behaviour**");    // 2 differences
  // File Header Parameter Value is out of spec. for affymetrix-chipsummary-het_rate Difference = 0.0390625 Epsilon = 0.0090100002
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon_empty, 0.00901F));

  Verbose::out(1, "**testcase1a-new behaviour-fEpsilon used**");   // 1 differences
  // File Header Parameter Value is out of spec. for affymetrix-chipsummary-het_rate Difference = 0.0390625 Epsilon = 0.0090100002
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon_empty, 0.00901F, 1.00001, true, 1000, 0.0014F));

  Verbose::out(1, "**testcase1b-new behaviour-fEpsilon used**");   // 0 differences
  CPPUNIT_ASSERT(Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon_empty, 0.00901F, 1.00001, true, 1000, 0.0016F));

  Verbose::out(1, "**testcase1c-new behaviour-mapEpsilon used**");   // 1 differences
  // File Header Parameter Value is out of spec. for affymetrix-chipsummary-hom_rate Difference = 0.01171875 Epsilon = 0.0100999996
  CPPUNIT_ASSERT(!Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon_full, 0.00901F, 1.00001, true, 1000, 0.00015F));

  Verbose::out(1, "**testcase1d-new behaviour-mapEpsilon used**");   // 0 differences
  CPPUNIT_ASSERT(Calvin::equivalent(firstFileName, secondFileName, setIgnore, mapEpsilon_full, 0.00901F, 1.00001, true, 1000, 0.00017F));

  Verbose::popMsgHandler();

  return;
}
