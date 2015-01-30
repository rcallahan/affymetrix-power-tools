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


#include "calvin_files/utils/test/StringUtilsTest.h"
//
#include "calvin_files/utils/src/StringUtils.h"
//

using namespace std;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( StringUtilsTest );

#define TEST_NON_EXISTANT_FILE "./no_file"
#define TEST_FILE "./FileUtilsTest.cpp"

void StringUtilsTest::setUp()
{
}

void StringUtilsTest::tearDown()
{
}

void StringUtilsTest::testCreation()
{
	StringUtils *utils = new StringUtils;
	CPPUNIT_ASSERT( utils != NULL );
	delete utils;
}

void StringUtilsTest::testmethod_STLTrimLeft()
{
	std::string s1("---Hello---");
	std::string s2("   Goodbye   ");
	StringUtils::STLTrimLeft(s1);
	CPPUNIT_ASSERT(s1 == "---Hello---");
	StringUtils::STLTrimLeft(s2, '-');
	CPPUNIT_ASSERT(s2 == "   Goodbye   ");
	StringUtils::STLTrimLeft(s1, '-');
	CPPUNIT_ASSERT(s1 == "Hello---");
	StringUtils::STLTrimLeft(s2);
	CPPUNIT_ASSERT(s2 == "Goodbye   ");
}

void StringUtilsTest::testmethod_STLTrimRight()
{
	std::string s1("---Hello---");
	std::string s2("   Goodbye   ");
	StringUtils::STLTrimRight(s1);
	CPPUNIT_ASSERT(s1 == "---Hello---");
	StringUtils::STLTrimRight(s2, '-');
	CPPUNIT_ASSERT(s2 == "   Goodbye   ");
	StringUtils::STLTrimRight(s1, '-');
	CPPUNIT_ASSERT(s1 == "---Hello");
	StringUtils::STLTrimRight(s2);
	CPPUNIT_ASSERT(s2 == "   Goodbye");
}

void StringUtilsTest::testmethod_STLTrimLeft_wide_version()
{
	std::wstring s1(L"---Hello---");
	std::wstring s2(L"   Goodbye   ");
	StringUtils::STLTrimLeft(s1);
	CPPUNIT_ASSERT(s1 == L"---Hello---");
	StringUtils::STLTrimLeft(s2, '-');
	CPPUNIT_ASSERT(s2 == L"   Goodbye   ");
	StringUtils::STLTrimLeft(s1, '-');
	CPPUNIT_ASSERT(s1 == L"Hello---");
	StringUtils::STLTrimLeft(s2);
	CPPUNIT_ASSERT(s2 == L"Goodbye   ");
}

void StringUtilsTest::testmethod_STLTrimRight_wide_version()
{
	std::wstring s1(L"---Hello---");
	std::wstring s2(L"   Goodbye   ");
	StringUtils::STLTrimRight(s1);
	CPPUNIT_ASSERT(s1 == L"---Hello---");
	StringUtils::STLTrimRight(s2, '-');
	CPPUNIT_ASSERT(s2 == L"   Goodbye   ");
	StringUtils::STLTrimRight(s1, '-');
	CPPUNIT_ASSERT(s1 == L"---Hello");
	StringUtils::STLTrimRight(s2);
	CPPUNIT_ASSERT(s2 == L"   Goodbye");
}

void StringUtilsTest::testmethod_ConvertWCSToMBS()
{
	std::wstring wide(L"Affymetrix");
	std::string s = StringUtils::ConvertWCSToMBS(wide);
	CPPUNIT_ASSERT(s == "Affymetrix");
}

void StringUtilsTest::testmethod_ConvertMBSToWCS()
{
	std::string wide("Affymetrix");
	std::wstring s = StringUtils::ConvertMBSToWCS(wide);
	CPPUNIT_ASSERT(s == L"Affymetrix");

}

void StringUtilsTest::testmethod_FormatString()
{
	wchar_t buffer[64];
	int count = 64;
	int value = 101;
	FormatString1(buffer, count, L"%d", value);
	CPPUNIT_ASSERT(std::wstring(buffer) == L"101");

	int min=0;
	int max=100;
	FormatString2(buffer, count, L"[%d..%d]", min, max);
	CPPUNIT_ASSERT(std::wstring(buffer) == L"[0..100]");

	int a=0;
	int b=1;
	int c=2;
	FormatString3(buffer, count, L"%d%d%d", a, b, c);
	CPPUNIT_ASSERT(std::wstring(buffer) == L"012");

	int d=3;
	FormatString4(buffer, count, L"%d%d%d%d", a, b, c, d);
	CPPUNIT_ASSERT(std::wstring(buffer) == L"0123");

	int e=4;
	FormatString5(buffer, count, L"%d%d%d%d%d", a, b, c, d, e);
	CPPUNIT_ASSERT(std::wstring(buffer) == L"01234");
}

void StringUtilsTest::testmethod_ToString()
{
	CPPUNIT_ASSERT(StringUtils::ToString(1, 2, '0') == L"01");
	CPPUNIT_ASSERT(StringUtils::ToString(1, 3, '0') == L"001");
	CPPUNIT_ASSERT(StringUtils::ToString(1000, 3, '0') == L"1000");
}
