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



#include "util/Err.h"
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
#include <vector>
//
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif
//
#include "util/CPPTest/Setup.h"

using namespace std;


/**
 * @class ErrTest
 * @brief cppunit class for testing error report functions.
 * last change by vliber on 01/22/08
 */
class ErrTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ErrTest );
  CPPUNIT_TEST( TestError );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void TestError();
  
};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ErrTest );



void ErrTest::TestError()
{
	Verbose::out(1, "***ErrTest testcase***");
    Verbose::out(1, "ErrTest::TestError");
    VerboseErrHandler *h1 = new VerboseErrHandler();
    VerboseErrHandler *h2 = new VerboseErrHandler();
    int size = Err::getParam().m_ErrHandlers.size();
    Err::pushHandler(h1);
    Err::pushHandler(h2);
    CPPUNIT_ASSERT(Err::getParam().m_ErrHandlers.size()==size+2);
    Err::popHandler();
    CPPUNIT_ASSERT(Err::getParam().m_ErrHandlers.size()==size+1);
    Err::popHandler();
    CPPUNIT_ASSERT(Err::getParam().m_ErrHandlers.size()==size);
    Err::setThrowStatus(true);//this order required for the following assert calls.
    POSITIVE_TEST(Err::check(true,"test2"));
    //Goal is to receive a message: FATAL ERROR: test2a -OK
    NEGATIVE_TEST(Err::check(false,"test2a"),Except);
    //Goal is to receive a message: FATAL ERROR: test2b -OK
    NEGATIVE_TEST(Err::errAbort("test2b"),Except);
    delete h1;
    delete h2;
}

