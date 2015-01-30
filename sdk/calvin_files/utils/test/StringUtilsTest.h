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

#ifndef __STRINGUTILSTEST_H_
#define __STRINGUTILSTEST_H_

#include <cppunit/extensions/HelperMacros.h>

class StringUtilsTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE (StringUtilsTest);

	CPPUNIT_TEST (testCreation);
	CPPUNIT_TEST (testmethod_STLTrimLeft);
	CPPUNIT_TEST (testmethod_STLTrimRight);
	CPPUNIT_TEST (testmethod_STLTrimLeft_wide_version);
	CPPUNIT_TEST (testmethod_STLTrimRight_wide_version);
	CPPUNIT_TEST (testmethod_ConvertWCSToMBS);
	CPPUNIT_TEST (testmethod_ConvertMBSToWCS);
	CPPUNIT_TEST (testmethod_FormatString);
	CPPUNIT_TEST (testmethod_ToString);

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testmethod_STLTrimLeft();
	void testmethod_STLTrimRight();
	void testmethod_STLTrimLeft_wide_version();
	void testmethod_STLTrimRight_wide_version();
	void testmethod_ConvertWCSToMBS();
	void testmethod_ConvertMBSToWCS();
	void testmethod_FormatString();
	void testmethod_ToString();
};

#endif // __STRINGUTILSTEST_H_
