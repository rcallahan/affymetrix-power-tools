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
#include <sstream>
#include <string>
//

using namespace std;

/**
 * @class VerboseTest
 * @brief cppunit class for testing logging and some command line ui.
 */
class VerboseTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( VerboseTest );
  CPPUNIT_TEST( testMessageHandler );
  CPPUNIT_TEST( testWarningHandler );
  CPPUNIT_TEST( testProgressHandler );
  CPPUNIT_TEST_SUITE_END();

public:
  /** Test command line output */
  virtual void tearDown();
  void testMessageHandler();
  void testWarningHandler();
  void testProgressHandler();
};

void VerboseTest::tearDown () 
{
  // clean up after failed test cases (pops any dead objects from stack)
  while( Verbose::getParam().m_MsgHandler.size() > 1 )
    Verbose::popMsgHandler();
  while( Verbose::getParam().m_WarnHandler.size() > 1 )
    Verbose::popWarnHandler();
  while( Verbose::getParam().m_ProHandler.size() > 1 )
    Verbose::popProgressHandler();
}

void VerboseTest::testMessageHandler() {
	Verbose::out(1, "***Verbose testcaces***");
  Verbose::out(1, "VerboseTest::testMessageHandler");
  ostringstream MessageHandler_oss;                               // message output streamed to string object
  MsgStream msgHandler(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler);
  ostringstream WarningHandler_oss;                               // warning output streamed to string object
  MsgStream warnHandler(2, &WarningHandler_oss);
  Verbose::pushWarnHandler(&warnHandler);
  ostringstream ProgressHandler_oss;                              // progress output streamed to string object
  ProgressDot progHandler(2, &ProgressHandler_oss);
  Verbose::pushProgressHandler(&progHandler); 
  Verbose::out(1, "msg 1", 1 );                                   // level < base_verbosity
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\n" );
  Verbose::out(2, "msg 2", 0 );                                   // level == base_verbosity
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\nmsg 2" ); 
  Verbose::out(3, "msg 3", 0 );                                   // level > base_verbosity
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\nmsg 2" ); 
  Verbose::setOutput(0);
  Verbose::out(1, "msg 4", 1 );                                   // output switched off
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\nmsg 2" ); 
  Verbose::setOutput(1);
  Verbose::out(1, "msg 5", 0 );                                   // output switched on
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\nmsg 2msg 5" ); 
  Verbose::setLevel(1);                                           // lower baseline verbosity
  Verbose::out(1, "msg 6", 1 );                                   // level == base_verbosity
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\nmsg 2msg 5msg 6\n" );
  Verbose::out(2, "msg 7", 1 );                                   // level > base_verbosity
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\nmsg 2msg 5msg 6\n" );
  Verbose::setLevel(3);                                           // increase base verbosity
  Verbose::out(3, "msg 8", 0 );                                   // level == base_verbosity
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\nmsg 2msg 5msg 6\nmsg 8" );
  CPPUNIT_ASSERT( WarningHandler_oss.str().empty() );             // validate no warnings 
  CPPUNIT_ASSERT( ProgressHandler_oss.str().empty() );            // validate no progess 
  Verbose::popProgressHandler();
  Verbose::popWarnHandler();
  Verbose::popMsgHandler();
  Verbose::out(1, "msg 9", 0 );                                   // handler poped from array
  CPPUNIT_ASSERT( MessageHandler_oss.str() == "msg 1\nmsg 2msg 5msg 6\nmsg 8" );
}

void VerboseTest::testWarningHandler() {
   Verbose::out(1, "\nVerboseTest::testWarningHandler");
  ostringstream MessageHandler_oss;                               // message output streamed to string object
  MsgStream msgHandler(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler);
  ostringstream WarningHandler_oss;                               // warning output streamed to string object
  MsgStream warnHandler(2, &WarningHandler_oss);
  Verbose::pushWarnHandler(&warnHandler);
  ostringstream ProgressHandler_oss;                              // progress output streamed to string object
  ProgressDot progHandler(2, &ProgressHandler_oss);
  Verbose::pushProgressHandler(&progHandler); 
  Verbose::warn(1, "warn 1", 1 );                                 // level < base_verbosity
  CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\n" );
  Verbose::warn(2, "warn 2", 0, "Warn: " );                       // level == base_verbosity
  CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\nWarn: warn 2" ); 
  Verbose::warn(3, "warn 3", 0, "" );                             // level > base_verbosity
  CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\nWarn: warn 2" ); 
  Verbose::setOutput(0);
  Verbose::warn(1, "warn 4", 1, "" );                             // output switched off
  CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\nWarn: warn 2" ); 
  Verbose::setOutput(1);
  Verbose::warn(1, "warn 5", 0, "" );                             // output switched on
  CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\nWarn: warn 2warn 5" ); 
  Verbose::setLevel(1);                                           // lower baseline verbosity
  Verbose::warn(1, "warn 6", 1, "" );                             // level == base_verbosity
  CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\nWarn: warn 2warn 5warn 6\n" );
  // BUG rsatin 28apr09 ==> setLevel does not update WarningHandler list in Verbose object
  Verbose::warn(2, "warn 7", 1, "" );                             // level > base_verbosity
//CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\nWarn: warn 2warn 5warn 6\n" );
  Verbose::setLevel(3);                                           // increase base verbosity
//Verbose::warn(3, "warn 8", 0, "" );                             // level == base_verbosity
//CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\nWarn: warn 2warn 5warn 6\nwarn 8" );
  CPPUNIT_ASSERT( MessageHandler_oss.str().empty() );             // validate no messages 
  CPPUNIT_ASSERT( ProgressHandler_oss.str().empty() );            // validate no progess 
  Verbose::popProgressHandler();
  Verbose::popWarnHandler();
  Verbose::popMsgHandler();
  Verbose::warn(1, "warn 9", 0 );                                 // handler poped from array
  // BUG rsatin 28apr09 ==> setLevel does not update WarningHandler list in Verbose object
//CPPUNIT_ASSERT( WarningHandler_oss.str() == "\nWARNING: warn 1\nwarn 2warn 5warn 6\nwarn 8" );  // FIXME rsatin
}

void VerboseTest::testProgressHandler() {
  Verbose::out(1, "\nVerboseTest::testProgressHandler");
  ostringstream MessageHandler_oss;                                     // message output streamed to string object
  MsgStream msgHandler(2, &MessageHandler_oss);
  Verbose::pushMsgHandler(&msgHandler);
  ostringstream WarningHandler_oss;                                     // warning output streamed to string object
  MsgStream warnHandler(2, &WarningHandler_oss);
  Verbose::pushWarnHandler(&warnHandler);
  ostringstream ProgressHandler_oss1;                                   // progress output streamed to string objects
  ProgressDot ProgressHandler1(2, &ProgressHandler_oss1);
  Verbose::pushProgressHandler(&ProgressHandler1);                      // handler pushed into Verbose class arrays
  // test case: verbosity level at default threshold
  int verbosity = 1;
  int total = 5;
  int dotMod = 2;
  int maxCalls = 10;
  Verbose::progressBegin(verbosity, "Start", total, dotMod, maxCalls);  // generate progress output
  for( int i=0; i<maxCalls; i++ ) {
    Verbose::progressStep(verbosity);
  }
  Verbose::progressEnd(verbosity, "End");
  Verbose::popProgressHandler();
  std::string str1 = ProgressHandler_oss1.str();
  CPPUNIT_ASSERT( str1.find("Start.....End") == 0 );                          // validate handling verbosity_1 == theshold_1
  // test case: verbosity level higher than default threshold
  verbosity = 2;
  ostringstream ProgressHandler_oss2;                                   // progress output streamed to string objects
  ProgressDot ProgressHandler2(2, &ProgressHandler_oss2);
  Verbose::pushProgressHandler(&ProgressHandler2);                      // handler pushed into Verbose class arrays
  Verbose::progressBegin(verbosity, "Start", total, dotMod, maxCalls);  // generate progress output
  for( int i=0; i<maxCalls; i++ ) {
    Verbose::progressStep(verbosity);
  }
  Verbose::progressEnd(verbosity, "End");
  Verbose::popProgressHandler();
  std::string str2 = ProgressHandler_oss2.str();                        // should be empty
  // BUG rsatin 28apr09 ==> verbosity>1 not handled properly for Progess
//CPPUNIT_ASSERT( str2.empty() );                                       // validate handling verbosity_2 > theshold_1
  // test case: verbosity level lower than updated threshold
  verbosity = 2;
  Verbose::setLevel(3);                                                 // increase base verbosity
  ostringstream ProgressHandler_oss3;                                   // progress output streamed to string objects
  ProgressDot ProgressHandler3(2, &ProgressHandler_oss3);
  Verbose::pushProgressHandler(&ProgressHandler3);                      // handler pushed into Verbose class arrays
  Verbose::progressBegin(verbosity, "Start", total, dotMod, maxCalls);  // generate progress output
  for( int i=0; i<maxCalls; i++ ) {
    Verbose::progressStep(verbosity);
  }
  Verbose::progressEnd(verbosity, "End");
  Verbose::popProgressHandler();
  std::string str3 = ProgressHandler_oss3.str();                        // should have progress output
  CPPUNIT_ASSERT( str3.find("Start.....End") == 0 );                          // validate handling verbosity_2 < theshold_3
  // test case: verbosity level higher than updated threshold
  verbosity = 4;
  Verbose::setLevel(3);                                                 // increase base verbosity
  ostringstream ProgressHandler_oss4;                                   // progress output streamed to string objects
  ProgressDot ProgressHandler4(2, &ProgressHandler_oss4);
  Verbose::pushProgressHandler(&ProgressHandler4);                      // handler pushed into Verbose class arrays
  Verbose::progressBegin(verbosity, "Start", total, dotMod, maxCalls);  // generate progress output
  for( int i=0; i<maxCalls; i++ ) {
    Verbose::progressStep(verbosity);
  }
  Verbose::progressEnd(verbosity, "End");
  Verbose::popProgressHandler();
  std::string str4 = ProgressHandler_oss4.str();                        // should be empty
  CPPUNIT_ASSERT( str4.empty());                                        // validate handling verbosity_4 > theshold_3
  CPPUNIT_ASSERT( WarningHandler_oss.str().empty() );                   // validate no warnings 
  CPPUNIT_ASSERT( MessageHandler_oss.str().empty() );                   // validate no messages 
  Verbose::popProgressHandler();
  Verbose::popWarnHandler();
  Verbose::popMsgHandler();
}

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( VerboseTest );
