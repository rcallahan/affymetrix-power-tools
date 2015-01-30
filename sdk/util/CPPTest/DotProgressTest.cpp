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

//
#include "util/DotProgress.h"
#include "util/Util.h"
#include "util/Verbose.h"
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
 * @class DotProgressTest
 * @brief cppunit class for testing conversion functions.
 * last change by vliber on 01/22/08
 */
class DotProgressTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( DotProgressTest );
  CPPUNIT_TEST( testdot );
  CPPUNIT_TEST( testnodots );
  CPPUNIT_TEST_SUITE_END();

public:
  /** Test reading rows. */
  void testdot();
  void testnodots();
  void setUp();
};

void DotProgressTest::setUp() {
  // no setup needed.
}

void DotProgressTest::testnodots()
{
	std::stringstream str;
	DotProgress dot(&str);
	dot.SetStepProperties(1, 10, 10);
	for (int i=0; i<100; i++)
		dot.Step(2);
	string s=str.str();
	CPPUNIT_ASSERT(s.length() == 0);
}

void DotProgressTest::testdot()
{
	std::cout<<endl;
	Verbose::out(1, "***DotProgress testcases***");
	Verbose::out(1, "DotProgressTest::testdot");
	std::stringstream str;
	DotProgress dot(&str);
	dot.SetStepProperties(1, 10, 10);
	for (int i=0; i<100; i++)
		dot.Step(1);
	string s=str.str();
	CPPUNIT_ASSERT(s.length() == 10);
	CPPUNIT_ASSERT(s == "..........");
	//********vliber*************
	//20
	dot.SetStepProperties(2, 10, 10);
	for (int i=0; i<100; i++)
		dot.Step(1);
	s=str.str();
	CPPUNIT_ASSERT(s.length() == 20);
	CPPUNIT_ASSERT(s == "....................");
   
	//add 5
	dot.SetStepProperties(1, 10, 20);
	for (int i=0; i<100; i++)
		dot.Step(1);
	s=str.str();
	cout<<endl;
	CPPUNIT_ASSERT(s.length() == 25);
	CPPUNIT_ASSERT(s == ".........................");
   
	//add 5
	dot.SetStepProperties(1, 20, 20);
	for (int i=0; i<100; i++)
		dot.Step(1);
	s=str.str();
	CPPUNIT_ASSERT(s.length() == 30);
	CPPUNIT_ASSERT(s == "..............................");
    
	//add 10
	dot.SetStepProperties(1, 20, 20);
	for (int i=0; i<200; i++)
		dot.Step(1);
	s=str.str();
	CPPUNIT_ASSERT(s.length() == 40);
	CPPUNIT_ASSERT(s == "........................................");

	// the same
	dot.SetStepProperties(0, 20, 20);
	for (int i=0; i<200; i++)
		dot.Step(1);
	s=str.str();
	CPPUNIT_ASSERT(s.length() == 40);
	CPPUNIT_ASSERT(s == "........................................");
    
	//add 100
	dot.SetStepProperties(1, 20, 0);
	for (int i=0; i<100; i++)
		dot.Step(1);
	s=str.str();
	CPPUNIT_ASSERT(s.length() == 140);
	CPPUNIT_ASSERT(s == "............................................................................................................................................");
}


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DotProgressTest );
