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


#ifndef __PARAMETERNAMEVALUEPAIRTEST_H_
#define __PARAMETERNAMEVALUEPAIRTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class ParameterNameValuePairTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( ParameterNameValuePairTest );

	CPPUNIT_TEST(test_assignment);
	CPPUNIT_TEST(test_equality_to_other_object);
	CPPUNIT_TEST(test_equality_to_string);
	CPPUNIT_TEST(test_accessibility);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void test_assignment();
	void test_equality_to_other_object();
	void test_equality_to_string();
	void test_accessibility();

};

#endif // __PARAMETERNAMEVALUEPAIRTEST_H_
