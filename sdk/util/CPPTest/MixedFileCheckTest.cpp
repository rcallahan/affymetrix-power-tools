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





//
#include "util/AffxString.h"
#include "util/Err.h"
#include "util/MixedFileCheck.h"
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

/**
 * @class MixedFileCheckTest
 * @brief cppunit class for testing functions from Util/MixedFileCheck class.
 *last change by rsatin on 09/17/09
 */


class MixedFileCheckTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MixedFileCheckTest );
  CPPUNIT_TEST( test );  
  CPPUNIT_TEST( test_small_diff );
  CPPUNIT_TEST( test_setMaxError );
  CPPUNIT_TEST( test_skipDataLines );
  CPPUNIT_TEST( test_reportLineColumn );
  CPPUNIT_TEST_SUITE_END();

public:
  
  virtual void tearDown();
  void test();
  void test_small_diff();
  void test_setMaxError();
  void test_skipDataLines();
  void test_reportLineColumn();

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MixedFileCheckTest );


void MixedFileCheckTest::tearDown () 
{
  // clean up after failed test cases (pops any dead objects from stack)
  while( Verbose::getParam().m_MsgHandler.size() > 1 )
    Verbose::popMsgHandler();
}


void MixedFileCheckTest::test()
{
   cout<<endl;
   Verbose::out(1, "***MixedFileCheck testcase***");
   Verbose::out(1, "MixedFileCheckTest::test");

  //positive
  MixedFileCheck mixedFileCheck1 (INPUT+"/testMixed1.txt",INPUT+"/testMixed1Gold.txt", 0.0,0,0);
  std::string errorMsg ("");
  CPPUNIT_ASSERT(mixedFileCheck1.check(errorMsg));
  POSITIVE_TEST(mixedFileCheck1.check(errorMsg));

  //negative
  // goal is to receive a message: Unable to open generated file ./input1/test1.txt -OK
  MixedFileCheck mixedFileCheck2 ("./input1/test1.txt",INPUT+"/test1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck2.check(errorMsg)); 
  
  
  //goal is to receive a message: Unable to open gold file input1\test1Gold.txt - OK
  MixedFileCheck mixedFileCheck3 (INPUT+"/test1.txt","input1/test1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck3.check(errorMsg)); 
  
  
  //goal is to receive a message:e The generated file, ./input/test5Bad.txt, has fewer lines than the gold file, ./input/test1Gold.txt -OK
  MixedFileCheck mixedFileCheck4 (INPUT+"/test5Bad.txt",INPUT+"/test1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck4.check (errorMsg));
  

  //goal is to receive a message: The generated file, ./input/test1.txt, has more lines than the gold file, ./input/test1Gold_Bad1.txt -OK
  MixedFileCheck mixedFileCheck5 (INPUT+"/test1.txt",INPUT+"/test1Gold_Bad1.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck5.check (errorMsg));
  

  //goal is to receive a message: Unable to open generated file -OK 
  MixedFileCheck mixedFileCheck6 ("","", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck6.check(errorMsg));
  
  //delete one cell data in the generated file
  //goal is to receive a message: Strings differ: asdfg and 5.3456789 
  //Mismatch reading generated file input/testMixed1Bad.txt:
  //gold line: '3 0 3.000123 asdfg 5.3456789 6.123456789 3'
  //generated line: '3 0 3.000123  5.3456789 6.123456789 3' 
  MixedFileCheck mixedFileCheck7 (INPUT+"/testMixed1Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck7.check(errorMsg));
  

  //differnce in numeric data
  //Numbers differ: 3.234 and 3
  //goal is to receive a message: There were 1 instances where numeric fields differed by more than the accepted tolerance: only 0 are allowed
  MixedFileCheck mixedFileCheck8 (INPUT+"/testMixed2Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck8.check(errorMsg));
  //positive allowMismath=1
  MixedFileCheck mixedFileCheck9 (INPUT+"/testMixed2Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,0,1);
  CPPUNIT_ASSERT(mixedFileCheck9.check(errorMsg));
  //positive eps=diff
  MixedFileCheck mixedFileCheck10 (INPUT+"/testMixed2Bad.txt",INPUT+"/testMixed1Gold.txt", 0.24,0,0);
  CPPUNIT_ASSERT(mixedFileCheck10.check(errorMsg));
  //positive skip first 4 lines (diff on line 4)
  MixedFileCheck mixedFileCheck11 (INPUT+"/testMixed2Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,4,0);
  CPPUNIT_ASSERT(mixedFileCheck11.check(errorMsg));


  //difference in the string data
  //goal is to receive a message:Strings differ: asdfg and aerty
  //Mismatch reading generated file input/testMixed3Bad.txt:
  //gold line: '3 0 3.000123 asdfg 5.3456789 6.123456789 3'
  //generated line: '3 0 3.000123 aerty 5.3456789 6.123456789 3'
  MixedFileCheck mixedFileCheck12 (INPUT+"/testMixed3Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck12.check(errorMsg));
  //positive skip first 5 lines (diff on line 5)
  MixedFileCheck mixedFileCheck13 (INPUT+"/testMixed3Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,5,0);
  CPPUNIT_ASSERT(mixedFileCheck13.check(errorMsg));


  //different fields count
  //goal is to receive a message: Unequal amount of whitespace delimited fileds in files
  //Mismatch reading generated file input/testMixed4Bad.txt:
  //gold line: '2 0 2.000123 3.234 4.3456789 5.123456789 2'
  //generated line: '2 0 2.000123 3.234 4.3456789 5.123456789 '
  MixedFileCheck mixedFileCheck14 (INPUT+"/testMixed4Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck14.check(errorMsg));
  

  //difference: string & double
  //goal is to receive a message: Strings differ: 7.3456789 and rrrrrr
  //Mismatch reading generated file input/testMixed5Bad.txt:
  //gold line: '5 0 5.000123 6.234 7.3456789 8.123456789 5'
  //generated line: '5 0 5.000123 6.234 rrrrrr 8.123456789 5'
  MixedFileCheck mixedFileCheck15 (INPUT+"/testMixed5Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(!mixedFileCheck15.check(errorMsg));
  
  //positive
  //additional empty last columnn in the generated file
  MixedFileCheck mixedFileCheck16 (INPUT+"/testMixed6Bad.txt",INPUT+"/testMixed1Gold.txt", 0.0,0,0);
  CPPUNIT_ASSERT(mixedFileCheck16.check(errorMsg));
}


void MixedFileCheckTest::test_small_diff()
{
  cout<<endl;
  Verbose::out(1, "MixedFileCheckTest::test_small_diff");

  std::string errorMsg ("");
  std::string strCount ("");
  char *endptr;
  int radix = 10;
  long errorCount = 0;

  // test case: absent vs present optional argument frac==eps
  Verbose::out(1, "**testcase1-old behaviour**");
  MixedFileCheck mixedFileCheck1 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.001, 0,0);
  CPPUNIT_ASSERT(!mixedFileCheck1.check(errorMsg));        // 2 differences
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 2 );
  Verbose::out(1, "**testcase1a-new behaviour**");
  MixedFileCheck mixedFileCheck1a (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.001, 0,0,0.001);
  CPPUNIT_ASSERT(mixedFileCheck1a.check(errorMsg));        // 0 differences
  POSITIVE_TEST(mixedFileCheck1a.check(errorMsg));

  // test case: absent vs present optional argument frac>eps
  Verbose::out(1, "**testcase2-old behaviour**");
  MixedFileCheck mixedFileCheck2 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 0,0);
  CPPUNIT_ASSERT(!mixedFileCheck2.check(errorMsg));        // 19 differences
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 19 );  
  Verbose::out(1, "**testcase2a-new behaviour**");
  MixedFileCheck mixedFileCheck2a (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 0,0,0.001);
  CPPUNIT_ASSERT(mixedFileCheck2a.check(errorMsg));        // 0 differences
  POSITIVE_TEST(mixedFileCheck2a.check(errorMsg));

  // test case: absent vs present optional argument frac<eps
  Verbose::out(1, "**testcase3-old behaviour**");
  MixedFileCheck mixedFileCheck3 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 0,0);
  CPPUNIT_ASSERT(!mixedFileCheck3.check(errorMsg));        // 19 differences
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 19 );  
  Verbose::out(1, "**testcase3a-new behaviour**");
  MixedFileCheck mixedFileCheck3a (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 0,0,0.00001);
  CPPUNIT_ASSERT(!mixedFileCheck3a.check(errorMsg));       // 6 differences
  POSITIVE_TEST(mixedFileCheck2a.check(errorMsg));
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 6 );  

  return;
}


void MixedFileCheckTest::test_setMaxError()
{
  cout<<endl;
  Verbose::out(1, "MixedFileCheckTest::test_setMaxError");

  std::string errorMsg ("");
  std::string strCount ("");
  char *endptr;
  int radix = 10;
  long errorCount = 0;
  long errorReportMax = 13;

  Verbose::out(1, "**testcase1-no-limit-specified**");

  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  MixedFileCheck mixedFileCheck1 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 0,0);
  CPPUNIT_ASSERT(!mixedFileCheck1.check(errorMsg));        // 19 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 19 );  

  // validate all differences were reported
  int find_count = 0;
  size_t fpos = 0;
  while( fpos != string::npos ) {
    fpos = MessageHandler1_oss.str().find("Numbers differ at",fpos);
	if(fpos != string::npos) {
      fpos++;
      find_count++;
      }
    }
  CPPUNIT_ASSERT( find_count == errorCount );

  // validate warning not reported
  CPPUNIT_ASSERT( MessageHandler1_oss.str().find("Number of differences exceeds maximum number") == string::npos );

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase2-limit-exceeded**");

  // message output streamed to string object
  ostringstream MessageHandler2_oss;
  MsgStream msgHandler2(2, &MessageHandler2_oss);
  Verbose::pushMsgHandler(&msgHandler2);

  MixedFileCheck mixedFileCheck2 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 0,0);
  mixedFileCheck2.setMaxError(errorReportMax);
  CPPUNIT_ASSERT(!mixedFileCheck2.check(errorMsg));        // 19 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 19 );  

  // validate count of differences reported matches specified limit
  find_count = 0;
  fpos = 0;
  while( fpos != string::npos ) {
    fpos = MessageHandler2_oss.str().find("Numbers differ at",fpos);
	if(fpos != string::npos) {
      fpos++;
      find_count++;
      }
    }
  CPPUNIT_ASSERT( find_count == errorReportMax );

  // validate warning is reported exactly once
  CPPUNIT_ASSERT( (fpos = MessageHandler2_oss.str().find("Number of differences exceeds maximum number"),0) != string::npos );
  CPPUNIT_ASSERT( MessageHandler2_oss.str().find("Number of differences exceeds maximum number",fpos+1) == string::npos );

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase3-limit-not-exceeded**");

  // message output streamed to string object
  ostringstream MessageHandler3_oss;
  MsgStream msgHandler3(2, &MessageHandler3_oss);
  Verbose::pushMsgHandler(&msgHandler3);

  MixedFileCheck mixedFileCheck3 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 0,0,0.00001);
  mixedFileCheck2.setMaxError(errorReportMax);
  CPPUNIT_ASSERT(!mixedFileCheck3.check(errorMsg));       // 6 differences
  POSITIVE_TEST(mixedFileCheck3.check(errorMsg));
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 6 );  

  // validate all differences were reported
  find_count = 0;
  fpos = 0;
  while( fpos != string::npos ) {
    fpos = MessageHandler3_oss.str().find("Numbers differ at",fpos);
	if(fpos != string::npos) {
      fpos++;
      find_count++;
      }
    }
  CPPUNIT_ASSERT( find_count == errorCount );

  // validate warning not reported
  CPPUNIT_ASSERT( MessageHandler3_oss.str().find("Number of differences exceeds maximum number") == string::npos );

  Verbose::popMsgHandler();

}

void MixedFileCheckTest::test_skipDataLines()
{
  cout<<endl;
  Verbose::out(1, "MixedFileCheckTest::test_skipDataLines");

  std::string errorMsg ("");
  std::string strCount ("");
  char *endptr;
  int radix = 10;
  long errorCount = 0;

  Verbose::out(1, "**testcase1-no-data-lines-skipped**");

  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  MixedFileCheck mixedFileCheck1 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 0,0);
  CPPUNIT_ASSERT(!mixedFileCheck1.check(errorMsg));        // 19 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 19 );  

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase2-one-data-lines-skipped**");

  // message output streamed to string object
  ostringstream MessageHandler2_oss;
  MsgStream msgHandler2(2, &MessageHandler2_oss);
  Verbose::pushMsgHandler(&msgHandler2);

  MixedFileCheck mixedFileCheck2 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 1,0);
  CPPUNIT_ASSERT(!mixedFileCheck2.check(errorMsg));        // 19 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 19 );  

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase3-two-data-lines-skipped**");

  // message output streamed to string object
  ostringstream MessageHandler3_oss;
  MsgStream msgHandler3(2, &MessageHandler3_oss);
  Verbose::pushMsgHandler(&msgHandler3);

  MixedFileCheck mixedFileCheck3 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 2,0);
  CPPUNIT_ASSERT(!mixedFileCheck3.check(errorMsg));        // 19 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 16 );  

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase4-25-data-lines-skipped**");

  // message output streamed to string object
  ostringstream MessageHandler4_oss;
  MsgStream msgHandler4(2, &MessageHandler4_oss);
  Verbose::pushMsgHandler(&msgHandler4);

  MixedFileCheck mixedFileCheck4 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 25,0);
  CPPUNIT_ASSERT(!mixedFileCheck4.check(errorMsg));        // 16 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 16 );  

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase5-26-data-lines-skipped**");

  // message output streamed to string object
  ostringstream MessageHandler5_oss;
  MsgStream msgHandler5(2, &MessageHandler5_oss);
  Verbose::pushMsgHandler(&msgHandler5);

  MixedFileCheck mixedFileCheck5 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 26,0);
  CPPUNIT_ASSERT(!mixedFileCheck5.check(errorMsg));        // 15 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 15 );  

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase6-all-data-lines-skipped**");

  // message output streamed to string object
  ostringstream MessageHandler6_oss;
  MsgStream msgHandler6(2, &MessageHandler6_oss);
  Verbose::pushMsgHandler(&msgHandler6);

  MixedFileCheck mixedFileCheck6 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.0001, 9999,0);
  CPPUNIT_ASSERT(mixedFileCheck6.check(errorMsg));        // 0 differences 

  Verbose::popMsgHandler();

}

void MixedFileCheckTest::test_reportLineColumn()
{
  cout<<endl;
  Verbose::out(1, "MixedFileCheckTest::test_reportLineColumn");

  std::string errorMsg ("");
  std::string strCount ("");
  char *endptr;
  int radix = 10;
  long errorCount = 0;

  Verbose::out(1, "**testcase1-gold-fewer-header_lines**");

  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  MixedFileCheck mixedFileCheck1 (INPUT+"/testMixedDiff1.txt",INPUT+"/testMixedDiff2.txt",0.001, 0,0);
  CPPUNIT_ASSERT(!mixedFileCheck1.check(errorMsg));        // 2 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 2 );  

  // validate expected line and column number of each of two differences was reported
  CPPUNIT_ASSERT( MessageHandler1_oss.str().find("Numbers differ at gold line 177 generated line 179 column 1:") != string::npos );
  CPPUNIT_ASSERT( MessageHandler1_oss.str().find("Numbers differ at gold line 215 generated line 217 column 1:") != string::npos );

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase2-gold-more-header_lines**");

  // message output streamed to string object
  ostringstream MessageHandler2_oss;
  MsgStream msgHandler2(2, &MessageHandler2_oss);
  Verbose::pushMsgHandler(&msgHandler2);

  MixedFileCheck mixedFileCheck2 (INPUT+"/testMixedDiff2.txt",INPUT+"/testMixedDiff1.txt",0.001, 0,0);
  CPPUNIT_ASSERT(!mixedFileCheck2.check(errorMsg));        // 2 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 2 );  

  // validate expected ZERO-based line and column number of each of two differences was reported
  CPPUNIT_ASSERT( MessageHandler2_oss.str().find("Numbers differ at gold line 179 generated line 177 column 1:") != string::npos );
  CPPUNIT_ASSERT( MessageHandler2_oss.str().find("Numbers differ at gold line 217 generated line 215 column 1:") != string::npos );

  Verbose::popMsgHandler();

  Verbose::out(1, "**testcase3-gold-same-header_lines**");

  // message output streamed to string object
  ostringstream MessageHandler3_oss;
  MsgStream msgHandler3(2, &MessageHandler3_oss);
  Verbose::pushMsgHandler(&msgHandler3);

  MixedFileCheck mixedFileCheck3 (INPUT+"/testMixedDiff3.txt",INPUT+"/testMixedDiff1.txt",0.001, 0,0);
  CPPUNIT_ASSERT(!mixedFileCheck3.check(errorMsg));        // 2 differences 
  strCount = errorMsg.substr( strlen("There were"), errorMsg.length()-strlen("There were") );
  errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  CPPUNIT_ASSERT( errorCount == 2 );  

  // validate expected ZERO-based line and column number of each of two differences was reported
  CPPUNIT_ASSERT( MessageHandler3_oss.str().find("Numbers differ at gold line 179 generated line 179 column 1:") != string::npos );
  CPPUNIT_ASSERT( MessageHandler3_oss.str().find("Numbers differ at gold line 217 generated line 217 column 1:") != string::npos );

  Verbose::popMsgHandler();

}
