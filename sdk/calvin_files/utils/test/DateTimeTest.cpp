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
#include "calvin_files/utils/test/DateTimeTest.h"
//
#include "calvin_files/exception/src/InterpretationException.h"
#include "calvin_files/utils/src/DateTime.h"
//
#include <cstdio>
#include <map>
//

using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( DateTimeTest );

void DateTimeTest::setUp()
{
}

void DateTimeTest::tearDown()
{
}

void DateTimeTest::testCreation()
{
	DateTime dt;
	CPPUNIT_ASSERT( 1 );
}

void DateTimeTest::testproperty_Date()
{
	DateTime dt;
	std::wstring date = L"Jan 15 2005";
	dt.Date(date);
	CPPUNIT_ASSERT( dt.Date() == date );
}

void DateTimeTest::testproperty_Time()
{
	DateTime dt;
	std::wstring time = L"12:10Z";
	dt.Time(time);
	CPPUNIT_ASSERT( dt.Time() == time );
}

void DateTimeTest::testmethod_Clear()
{
	DateTime dt;
	std::wstring time = L"12:10Z";
	dt.Time(time);
	std::wstring date = L"Jan 15 2005";
	dt.Date(date);
	dt.Clear();
	CPPUNIT_ASSERT( dt.Time().compare(L"") == 0);
	CPPUNIT_ASSERT( dt.Date().compare(L"") == 0);
}

void DateTimeTest::testmethod_GetCurrentDateTime()
{
	DateTime dt = DateTime::GetCurrentDateTime();

	int year=-1;
	int month=-1;
	int day=-1;
	int hour=-2;
	int minute=-2;
	int second=-2;
	int nparsed=0;

	nparsed = swscanf(dt.Date().c_str(), L"%d-%d-%d", &year, &month, &day);
	CPPUNIT_ASSERT(nparsed == 3);
	CPPUNIT_ASSERT( year > 2004 && year < 3000 );
	CPPUNIT_ASSERT( month > 0 && month < 13 );
	CPPUNIT_ASSERT( day > 0 && day < 32 );

	nparsed = swscanf(dt.Time().c_str(), L"%d:%d:%d", &hour, &minute, &second);
	CPPUNIT_ASSERT(nparsed == 3);
	CPPUNIT_ASSERT( hour > -1 && hour < 25 );
	CPPUNIT_ASSERT( minute > -1 && minute < 61 );
	CPPUNIT_ASSERT( second > -1 && second < 61 );
}

void DateTimeTest::testmethod_Parse()
{
	DateTime dt1 = DateTime::GetCurrentDateTime();
	std::wstring dateString1 = dt1.ToString();
	DateTime dt2 = DateTime::Parse(dateString1);
	std::wstring dateString2 = dt2.ToString();
	CPPUNIT_ASSERT(dateString1 == dateString2);
	CPPUNIT_ASSERT(dt2.IsUTC() == true);
	std::wstring dateString3(L"2005-12-25T10:13:45");
	DateTime dt3 = DateTime::Parse(dateString3);
	CPPUNIT_ASSERT(dt3.Date() == L"2005-12-25");
	CPPUNIT_ASSERT(dt3.Time() == L"10:13:45");
	CPPUNIT_ASSERT(dt3.IsUTC() == false);
	std::wstring dateString4(L"2005/12/24T10:13:46Z");
	DateTime dt4 = DateTime::Parse(dateString4);
	CPPUNIT_ASSERT(dt4.Date() == L"2005-12-24");
	CPPUNIT_ASSERT(dt4.Time() == L"10:13:46");
	CPPUNIT_ASSERT(dt4.IsUTC() == true);
}

void DateTimeTest::testmethod_ParseFail()
{
	std::wstring dateString(L"This is a bogus time");
	CPPUNIT_ASSERT_THROW(DateTime::Parse(dateString), affymetrix_calvin_exceptions::FormatException);
	// slighly less bogus time
	dateString = L"Jan 5 2005T10:13:45";
	CPPUNIT_ASSERT_THROW(DateTime::Parse(dateString), affymetrix_calvin_exceptions::FormatException);
}

void DateTimeTest::FormatDateTest()
{
	u_int32_t d = 22;
	u_int32_t m = 11;
	u_int32_t y = 1995;
	CPPUNIT_ASSERT(DateTime::FormatDate(y, m, d) == L"1995-11-22");
}

void DateTimeTest::FormatTimeTest()
{
	u_int32_t h = 7;
	u_int32_t m = 6;
	u_int32_t s = 5;
	CPPUNIT_ASSERT(DateTime::FormatTime(h, m, s) == L"07:06:05");

}

void DateTimeTest::FormatDateTimeTest()
{
	u_int32_t d = 22;
	u_int32_t m = 11;
	u_int32_t y = 1995;
	u_int32_t h = 7;
	u_int32_t min = 6;
	u_int32_t s = 5;
	CPPUNIT_ASSERT(DateTime::FormatDateTime(y, m, d, h, min, s, false) == L"1995-11-22T07:06:05");

}
