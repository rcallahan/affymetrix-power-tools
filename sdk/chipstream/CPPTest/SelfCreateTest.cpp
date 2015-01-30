////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

/** 
 * @file   SelfCreateTest.cpp
 * @author csugne
 * @date   Tue Nov  8 08:54:40 PST 2005
 * 
 * @brief  Testing the SelfCreate functions.
 * 
 */
#ifndef SELFCREATETEST_H
#define SELFCREATETEST_H

#include "chipstream/SelfCreate.h"
//
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <string>
#include <vector>
//

using namespace std;
/**
 * @class SelfCreateTest
 * @brief cppunit class for testing conversion functions.
 */
class SelfCreateTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SelfCreateTest );
  CPPUNIT_TEST( testFillInNameParam );
  CPPUNIT_TEST_SUITE_END();

public:
  // blank test.
  void testFillInNameParam();
  
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( SelfCreateTest );

void SelfCreateTest::testFillInNameParam() {
  map<string, string> paramMap;
  string name;
  string s("targetName.paramOne=param1.paramTwo=2.0.param3=paramThree");
  string s2("quant-norm.sketch=10000.bioc=true");
  string s3("quant-norm.sketch=10000.bioc=2.0");
  /* Try with float entry. */
  SelfCreate::fillInNameParam(s, name, paramMap);
  CPPUNIT_ASSERT( name == "targetName" );
  CPPUNIT_ASSERT( paramMap["paramOne"] == "param1" );
  CPPUNIT_ASSERT( paramMap["paramTwo"] == "2.0" );
  CPPUNIT_ASSERT( paramMap["param3"] == "paramThree" );
  CPPUNIT_ASSERT( paramMap.find("notThere") == paramMap.end() );
  /* Try with no float entry. */
  SelfCreate::fillInNameParam(s2, name, paramMap);
  CPPUNIT_ASSERT( name == "quant-norm" );
  CPPUNIT_ASSERT( paramMap["sketch"] == "10000" );
  CPPUNIT_ASSERT( paramMap["bioc"] == "true" );
  /* Try with terminating float entry. */
  SelfCreate::fillInNameParam(s3, name, paramMap);
  CPPUNIT_ASSERT( name == "quant-norm" );
  CPPUNIT_ASSERT( paramMap["sketch"] == "10000" );
  CPPUNIT_ASSERT( paramMap["bioc"] == "2.0" );
}

#endif 
