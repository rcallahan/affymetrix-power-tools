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
#include "util/Util.h"
#include "util/Verbose.h"
#include <util/MsgStream.h>
#include <util/OutputMessageStream.h>
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <string>
//

using namespace std;

/**
 * @class OutputMessageStreamTest
 * @brief cppunit class for testing conversion functions.
 */
class OutputMessageStreamTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( OutputMessageStreamTest );
  CPPUNIT_TEST( testostream_str );
  CPPUNIT_TEST( testostream_cout );
  CPPUNIT_TEST( testproperties );
  CPPUNIT_TEST_SUITE_END();

public:
  /** Test reading rows. */
  void testostream_str();
  void testostream_cout();
  void testproperties();
  void setUp();
};

void OutputMessageStreamTest::setUp() {
}

void OutputMessageStreamTest::testostream_str()
{
	cout<<endl;
	Verbose::out(1, "***OutputMessageStreamTest testcases***");
	 Verbose::out(1, "OutputMessageStreamTest::testostream_str");
	std::stringstream str;
	OutputMessageStream out(1, &str);
	out.Write(1, string("test1"));
	string s = str.str();
	CPPUNIT_ASSERT(s == "test1");

	out.Write(2, string("test2"));
	s = str.str();
	CPPUNIT_ASSERT(s == "test1");

	out.Write(1, string("\ntest2"));
	s = str.str();
	CPPUNIT_ASSERT(s == "test1\ntest2");

    //********vliber*********
	out.Write(1, string(" #$_%-|}"));
	CPPUNIT_ASSERT(str.str() == "test1\ntest2 #$_%-|}");
}

void OutputMessageStreamTest::testostream_cout()
{
	Verbose::out(1, "OutputMessageStreamTest::testostream_cout");
	OutputMessageStream out(1);
	out.Write(1, string("test"));	
}

void OutputMessageStreamTest::testproperties()
{
	cout<<endl;
	Verbose::out(1, "OutputMessageStreamTest::testproperties");
	OutputMessageStream out(1);
	CPPUNIT_ASSERT(out.GetLevel() == 1);
	out.SetLevel(2);
	CPPUNIT_ASSERT(out.GetLevel() == 2);
	//********vliber*********
	std::stringstream str;
	OutputMessageStream out1(2, &str);
	out1.SetLevel(1);
	out1.Write(2, string("test"));
	string s = str.str();
	CPPUNIT_ASSERT(s == "");
}

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( OutputMessageStreamTest );
