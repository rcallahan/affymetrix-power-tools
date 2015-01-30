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
#include "calvin_files/data/test/CHPDataTest.h"
//
#include "calvin_files/data/src/CHPData.h"
#include "calvin_files/utils/src/AffyStlCollectionTypes.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPDataTest );

void CHPDataTest::setUp()
{
}

void CHPDataTest::tearDown()
{
}

void CHPDataTest::testCreation()
{
	CHPData data("CHPData", CHP_EXPRESSION_ASSAY_TYPE);
	CPPUNIT_ASSERT(1);
}

void CHPDataTest::FilenameTest()
{
	CHPData data;

	data.SetFilename("CHPData");
	CPPUNIT_ASSERT(data.GetFilename() == "CHPData");
}

void CHPDataTest::SetEntryCountTest()
{
	CHPData data("CHPData", CHP_EXPRESSION_ASSAY_TYPE);

	u_int32_t p1 = 100;
	data.SetEntryCount(p1, false);
	CPPUNIT_ASSERT(1);
}

void CHPDataTest::ArrayTypeTest()
{
	CHPData data;
	CPPUNIT_ASSERT_NO_THROW(data.SetArrayType(L"coolio"));
	CPPUNIT_ASSERT(data.GetArrayType() == L"coolio");
}

void CHPDataTest::AlgParamTest()
{
	CHPData data;
	data.AddAlgParam(L"param1", L"value");
	data.AddAlgParam(L"param1", 1);
	data.AddAlgParam(L"param1", 1.0f);
	CPPUNIT_ASSERT(1);
}

void CHPDataTest::ChipSumTest()
{
	CHPData data;
	data.AddChipSum(L"param1", L"value");
	data.AddChipSum(L"param2", 1);
	data.AddChipSum(L"param3", 123.321f);

	ParameterNameValueType p = data.GetChipSum(L"param1");
	CPPUNIT_ASSERT(p.GetName() == L"param1");
	CPPUNIT_ASSERT(p.GetValueText() == L"value");

	p = data.GetChipSum(L"param2");
	CPPUNIT_ASSERT(p.GetName() == L"param2");
	CPPUNIT_ASSERT(p.GetValueInt32() == 1);

	p = data.GetChipSum(L"param3");
	CPPUNIT_ASSERT(p.GetName() == L"param3");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.GetValueFloat(), 123.321, 0.001f);
}

void CHPDataTest::ProgIdTest()
{
	CHPData data;
	CPPUNIT_ASSERT_NO_THROW(data.SetProgId(L"progId"));
	CPPUNIT_ASSERT(data.GetProgId() == L"progId");
}

void CHPDataTest::AlgVersionTest()
{
	CHPData data;
	CPPUNIT_ASSERT_NO_THROW(data.SetAlgVersion(L"version"));
	CPPUNIT_ASSERT(data.GetAlgVersion() == L"version");
}

void CHPDataTest::AlgNameTest()
{
	CHPData data;
	CPPUNIT_ASSERT_NO_THROW(data.SetAlgName(L"alg name"));
	CPPUNIT_ASSERT(data.GetAlgName() == L"alg name");
}

void CHPDataTest::RowsTest()
{
	CHPData data;
	CPPUNIT_ASSERT_NO_THROW(data.SetRows(234));
	CPPUNIT_ASSERT(data.GetRows() == 234);
}

void CHPDataTest::ColsTest()
{
	CHPData data;
	CPPUNIT_ASSERT_NO_THROW(data.SetCols(8760));
	CPPUNIT_ASSERT(data.GetCols() == 8760);
}
