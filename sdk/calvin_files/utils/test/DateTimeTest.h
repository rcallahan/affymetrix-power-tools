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

#ifndef __DATETIMETEST_H_
#define __DATETIMETEST_H_

#include <cppunit/extensions/HelperMacros.h>

class DateTimeTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( DateTimeTest );

	CPPUNIT_TEST ( testCreation );
	CPPUNIT_TEST ( testproperty_Date );
	CPPUNIT_TEST ( testproperty_Time );
	CPPUNIT_TEST ( testmethod_Clear );
	CPPUNIT_TEST ( testmethod_GetCurrentDateTime );
	CPPUNIT_TEST ( testmethod_Parse );
	CPPUNIT_TEST ( testmethod_ParseFail );
	CPPUNIT_TEST ( FormatDateTest );
	CPPUNIT_TEST ( FormatTimeTest );
	CPPUNIT_TEST ( FormatDateTimeTest );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();

	void testCreation();
	void testproperty_Date();
	void testproperty_Time();
	void testmethod_Clear();
	void testmethod_GetCurrentDateTime();
	void testmethod_Parse();
	void testmethod_ParseFail();
	void FormatDateTest();
	void FormatTimeTest();
	void FormatDateTimeTest();

};

#endif // __DATETIMETEST_H_
