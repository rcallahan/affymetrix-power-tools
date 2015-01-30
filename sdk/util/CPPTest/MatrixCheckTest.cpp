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





#include "util/AffxString.h"
#include "util/MatrixCheck.h"
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
 * @class MatrixCheckTest
 * @brief cppunit class for testing functions from Util/MatrixCheck class.
 *last change by rsatin on 09/17/09
 */


class MatrixCheckTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MatrixCheckTest );
  CPPUNIT_TEST( test );
  CPPUNIT_TEST( test_small_diff );
  CPPUNIT_TEST( test_setPrintMismatch );  
  CPPUNIT_TEST( test_setPrintMismatchMax );  
  CPPUNIT_TEST( test_reportLineColumn );
  CPPUNIT_TEST_SUITE_END();

public:
  
  virtual void tearDown();
  void test();
  void test_small_diff();
  void test_setPrintMismatch();
  void test_setPrintMismatchMax();
  void test_reportLineColumn();

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MatrixCheckTest );


void MatrixCheckTest::tearDown () 
{
  // clean up after failed test cases (pops any dead objects from stack)
  while( Verbose::getParam().m_MsgHandler.size() > 1 )
    Verbose::popMsgHandler();
}


void MatrixCheckTest ::test()
{
   cout<<endl;
   Verbose::out(1, "***MatrixCheck testcase***");
   Verbose::out(1, "MatrixCheckTest::test");
  //positive
  
  MatrixCheck textCheck1 (INPUT+"/test1.txt",INPUT+"/test1Gold.txt", 0.0,1,0,false,0);
  std::string errorMsg ("");
  CPPUNIT_ASSERT(textCheck1.check(errorMsg));
  POSITIVE_TEST(textCheck1.check(errorMsg));
  
  //play with ignored rows and errors count
  //Goal is to receive a message: File: input\test1Bad.txt vs .\input\test1Gold.txt: Expecting no more than 0 found: 4 -OK
  MatrixCheck textCheck2 (INPUT+"/test1Bad.txt",INPUT+"/test1Gold.txt", 0.0,1,0,false,0);
  CPPUNIT_ASSERT(!textCheck2.check(errorMsg));
  if (! textCheck2.check (errorMsg))
    cout << errorMsg << endl;
  //positive result, because diff. count is 4
  MatrixCheck textCheck3 (INPUT+"/test1Bad.txt",INPUT+"/test1Gold.txt", 0.0,1,0,false,4);
  CPPUNIT_ASSERT(textCheck3.check(errorMsg));
  
  //Goal is to receive a message: File: input\test3Bad.txt vs .\input\test1Gold.txt: Expecting no more than 0 found: 3 -OK
  MatrixCheck textCheck4 (INPUT+"/test3Bad.txt",INPUT+"/test1Gold.txt", 0.0,4,0,false,0);
  CPPUNIT_ASSERT(!textCheck4.check(errorMsg));
  if (! textCheck4.check (errorMsg))
    cout << errorMsg << endl;
  //positive result, because diff. count is 3 
  MatrixCheck textCheck5 (INPUT+"/test3Bad.txt", INPUT+"/test1Gold.txt", 0.0,4,0,false,3);
  CPPUNIT_ASSERT(textCheck5.check(errorMsg));
  
  //Goal is to receive a message: File: input\test3Bad.txt vs .\input\test1Gold.txt: Expecting no more than 0 found: 1 -OK
  MatrixCheck textCheck6 (INPUT+"/test3Bad.txt",INPUT+"/test1Gold.txt", 0.0,5,0,false,0);
  CPPUNIT_ASSERT(!textCheck6.check(errorMsg));
  if (! textCheck6.check (errorMsg))
    cout << errorMsg << endl;
  //positive result, because diff. count is 1 
  MatrixCheck textCheck7 (INPUT+"/test3Bad.txt",INPUT+"/test1Gold.txt", 0.0,5,0,false,1);
  CPPUNIT_ASSERT(textCheck7.check(errorMsg));
  //positive result, because diff. count is 0 after skipping 8 rows 
  MatrixCheck textCheck8(INPUT+"/test3Bad.txt",INPUT+"/test1Gold.txt", 0.0,8,0,false,0);
  CPPUNIT_ASSERT(textCheck8.check(errorMsg));
  
  //play with ignored col and errors count
  //Goal is to receive a message: File: input\test3Bad.txt vs .\input\test1Gold.txt: Expecting no more than 0 found: 1 -OK
  MatrixCheck textCheck9 (INPUT+"/test3Bad.txt",INPUT+"/test1Gold.txt", 0.0,1,1,false,0);
  CPPUNIT_ASSERT(!textCheck9.check(errorMsg));
  if (! textCheck9.check (errorMsg))
    cout << errorMsg << endl;
  //positive result, because diff. count is 1 after skipping 1 col.
  MatrixCheck textCheck10 (INPUT+"/test3Bad.txt",INPUT+"/test1Gold.txt", 0.0,1,1,false,1);
  CPPUNIT_ASSERT(textCheck10.check(errorMsg));

  //play with epsilon
  //Goal is to receive a message: File: input\test4Bad.txt vs .\input\test1Gold.txt: Expecting no more than 0 found: 1 -OK
  MatrixCheck textCheck11 (INPUT+"/test4Bad.txt",INPUT+"/test1Gold.txt", 0.0,1,0,false,0);
  CPPUNIT_ASSERT(!textCheck11.check(errorMsg));
  if (! textCheck11.check (errorMsg))
    cout << errorMsg << endl;
  //positive result, because epsion == diff.
  MatrixCheck textCheck12 (INPUT+"/test4Bad.txt",INPUT+"/test1Gold.txt", 0.002,1,0,false,0);
  CPPUNIT_ASSERT(textCheck12.check(errorMsg));

  //play with matchNames
  //Goal is to receive a message: File: input\test2.txt vs .\input\test1Gold.txt: Expecting no more than 0 found: 60 -OK
  MatrixCheck textCheck13 (INPUT+"/test2.txt",INPUT+"/test1Gold.txt", 0.0,1,0,false,0);
  CPPUNIT_ASSERT(!textCheck13.check(errorMsg));
  if (! textCheck13.check (errorMsg))
    cout << errorMsg << endl;
  //positive result, because matchNames.
  MatrixCheck textCheck14 (INPUT+"/test2.txt",INPUT+"/test1Gold.txt", 0.0,1,0,true,0);
  CPPUNIT_ASSERT(textCheck14.check(errorMsg));

  //negative
  // Goal is to receive a message: Can't open file: .\input1\test1.txt to read. -OK
  MatrixCheck textCheck15 ("./input1/test1.txt","./input/test1Gold.txt", 0.0,1,0,false,0);
  CPPUNIT_ASSERT(!textCheck15.check (errorMsg));
  if (! textCheck15.check (errorMsg))
    cout << errorMsg << endl;
  
  // Goal is to receive a message: Can't open file: input1\test1Gold.txt to read. - OK
  MatrixCheck textCheck16 (INPUT+"/test1.txt","input1/test1Gold.txt", 0.0,1,0,false,0);
  CPPUNIT_ASSERT(!textCheck16.check(errorMsg)); 
  if (! textCheck16.check (errorMsg))
    cout << errorMsg << endl;
  
  // Goal is to receive a message: invalid ColSkip, RowSkip, and/or Epsilon. - OK
  MatrixCheck textCheck17 (INPUT+"/test1.txt",INPUT+"/test1Gold.txt", -0.01,1,0,false,0);
  CPPUNIT_ASSERT(!textCheck17.check(errorMsg)); 
  if (! textCheck17.check (errorMsg))
    cout << errorMsg << endl;
  
  //Goal is to receive a message: FATAL ERROR: RowFile::matrixFromFile() - Number of skipCols >= number of cols. -OK
  //MatrixCheck textCheck18 ("input/test1.txt","input/test1Gold.txt", 0.0,1,7,false,0);
  //CPPUNIT_ASSERT(!textCheck18.check(errorMsg));  
}


void MatrixCheckTest::test_small_diff()
{
   cout<<endl;
   Verbose::out(1, "MatrixCheckTest::test_small_diff");
  //positive

  std::string errorMsg ("");
  std::string strCount ("");
  size_t strStart;
  char *endptr;
  int radix = 10;
  long errorCount = 0;

  Verbose::out(1, "**testcase1-old behaviour**");
  MatrixCheck textCheck1a (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.001,1,0,false,0);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck1a.check(errorMsg));
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount == 2 );

  Verbose::out(1, "**testcase1-new behaviour**");
  MatrixCheck textCheck1b (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.001,1,0,false,0,0.001);
  errorMsg.clear();
  CPPUNIT_ASSERT(textCheck1b.check(errorMsg));
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount == 0 );

  Verbose::out(1, "**testcase2-old behaviour**");
  MatrixCheck textCheck2a (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.0001,1,0,false,0);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck2a.check(errorMsg));
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount == 10 );

  Verbose::out(1, "**testcase2-new behaviour**");
  MatrixCheck textCheck2b (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.0001,1,0,false,0,0.001);
  errorMsg.clear();
  CPPUNIT_ASSERT(textCheck2b.check(errorMsg));
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount == 0 );

  Verbose::out(1, "**testcase3-old behaviour**");
  MatrixCheck textCheck3a (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.0001,1,0,false,0);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck3a.check(errorMsg));
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount == 10 );

  Verbose::out(1, "**testcase3-new behaviour**");
  MatrixCheck textCheck3b (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.0001,1,0,false,0,0.00001);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck3b.check(errorMsg));
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount == 4 );

}


void MatrixCheckTest::test_setPrintMismatch()
{
  cout<<endl;
  Verbose::out(1, "MatrixCheckTest::test_setPrintMismatch");

  std::string errorMsg ("");
  std::string strCount ("");
  char *endptr;
  int radix = 10;

  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  Verbose::out(1, "**testcase1-old behaviour**");
  MatrixCheck textCheck1a (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.001,1,0,false,0);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck1a.check(errorMsg));

  // get difference count from error message
  long errorCount = 0;
  size_t strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount>0 );

  // validate NO differences were reported on the console for default case
  CPPUNIT_ASSERT( MessageHandler1_oss.str().find(") Diff:") == string::npos );

  Verbose::popMsgHandler();

  // message output streamed to string object
  ostringstream MessageHandler2_oss;
  MsgStream msgHandler2(2, &MessageHandler2_oss);
  Verbose::pushMsgHandler(&msgHandler2);

  Verbose::out(1, "**testcase1-new behaviour**");
  MatrixCheck textCheck1b (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.001,1,0,false,0);
  textCheck1b.setPrintMismatch(true);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck1b.check(errorMsg));

  // get difference count from error message
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount>0 );

  // validate ALL differences were reported
  int find_count = 0;
  size_t fpos = 0;
  while( fpos != string::npos ) {
    fpos = MessageHandler2_oss.str().find(") Diff:",fpos);
	if(fpos != string::npos) {
      fpos++;
      find_count++;
      }
    }
  CPPUNIT_ASSERT( find_count == errorCount );

  Verbose::popMsgHandler();
}


void MatrixCheckTest::test_setPrintMismatchMax()
{
  cout<<endl;
  Verbose::out(1, "MatrixCheckTest::test_setPrintMismatchMax");

  std::string errorMsg ("");
  std::string strCount ("");
  char *endptr;
  int radix = 10;
  int PrintMismatchMax = 5;

  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  Verbose::out(1, "**testcase1-old behaviour**");
  MatrixCheck textCheck1a (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.0001,1,0,false,0);
  textCheck1a.setPrintMismatch(true);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck1a.check(errorMsg));

  // get difference count from error message
  long errorCount = 0;
  size_t strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount>0 );

  // validate all differences were reported
  int find_count = 0;
  size_t fpos = 0;
  while( fpos != string::npos ) {
    fpos = MessageHandler1_oss.str().find(") Diff:",fpos);
	if(fpos != string::npos) {
      fpos++;
      find_count++;
      }
    }
  CPPUNIT_ASSERT( find_count == errorCount );

  // validate print mismatch limit not displayed
  CPPUNIT_ASSERT( MessageHandler1_oss.str().find("Number of differences exceeds maximum number") == string::npos );

  Verbose::popMsgHandler();

  // message output streamed to string object
  ostringstream MessageHandler2_oss;
  MsgStream msgHandler2(2, &MessageHandler2_oss);
  Verbose::pushMsgHandler(&msgHandler2);

  Verbose::out(1, "**testcase1-new behaviour**");
  MatrixCheck textCheck1b (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.0001,1,0,false,0);
  textCheck1b.setPrintMismatch(true);
  textCheck1b.setPrintMismatchMax(PrintMismatchMax);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck1b.check(errorMsg));

  // get difference count from error message
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount>0 );

  // validate all differences were reported
  find_count = 0;
  fpos = 0;
  while( fpos != string::npos ) {
    fpos = MessageHandler2_oss.str().find(") Diff:",fpos);
	if(fpos != string::npos) {
      fpos++;
      find_count++;
      }
    }
  CPPUNIT_ASSERT( find_count == PrintMismatchMax );

  // validate print mismatch limit displayed
  CPPUNIT_ASSERT( MessageHandler2_oss.str().find("Number of differences exceeds maximum number") > 0 );

  Verbose::popMsgHandler();
}


void MatrixCheckTest::test_reportLineColumn()
{
  cout<<endl;
  Verbose::out(1, "MatrixCheckTest::test_reportLineColumn");

  std::string errorMsg ("");
  std::string strCount ("");
  char *endptr;
  int radix = 10;

  // message output streamed to string object
  ostringstream MessageHandler1_oss;
  MsgStream msgHandler1(2, &MessageHandler1_oss);
  Verbose::pushMsgHandler(&msgHandler1);

  Verbose::out(1, "**testcase1-skip-one-row-zero-column**");
  MatrixCheck textCheck1a (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.001,1,0,false,0);
  textCheck1a.setPrintMismatch(true);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck1a.check(errorMsg));

  // get difference count from error message
  long errorCount = 0;
  size_t strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount==2 );

  // validate expected ZERO-based line and column number of each of two differences was reported
  CPPUNIT_ASSERT( MessageHandler1_oss.str().find("row: 6 col: 2 ") != string::npos );
  CPPUNIT_ASSERT( MessageHandler1_oss.str().find("row: 9 col: 2 ") != string::npos );

  Verbose::popMsgHandler();

  // message output streamed to string object
  ostringstream MessageHandler2_oss;
  MsgStream msgHandler2(2, &MessageHandler2_oss);
  Verbose::pushMsgHandler(&msgHandler2);

  Verbose::out(1, "**testcase2-skip-three-row-one-column**");
  MatrixCheck textCheck2a (INPUT+"/testMatrixDiff1.txt",INPUT+"/testMatrixDiff2.txt", 0.001,3,1,false,0);
  textCheck2a.setPrintMismatch(true);
  errorMsg.clear();
  CPPUNIT_ASSERT(!textCheck2a.check(errorMsg));

  // get difference count from error message
  errorCount = 0;
  strStart = errorMsg.find("found: ");
  if( strStart != string::npos ) {
    strCount = errorMsg.substr( strStart+strlen("found: "), errorMsg.length()-strStart-strlen("found: ") );
    errorCount = strtol( strCount.c_str(), &endptr, radix ); 
  }
  CPPUNIT_ASSERT( errorCount>0 );

  // validate expected ZERO-based line and column number of each of two differences was reported
  //CPPUNIT_ASSERT( MessageHandler2_oss.str().find("row: 6 col: 2 ") != string::npos );   //TODO rsatin row,col changed, is it a bug?
  //CPPUNIT_ASSERT( MessageHandler2_oss.str().find("row: 9 col: 2 ") != string::npos );  //TODO rsatin row,col changed, is it a bug?

}

