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
#include "calvin_files/data/test/CHPTilingDataTest.h"
//
#include "calvin_files/data/src/CHPTilingData.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPTilingDataTest );

void CHPTilingDataTest::setUp()
{
}

void CHPTilingDataTest::tearDown()
{
}

void CHPTilingDataTest::test_FileName()
{
	CHPTilingData data;
	data.SetFilename("file");
	CPPUNIT_ASSERT( data.GetFilename() == "file");
}

void CHPTilingDataTest::test_NumberSequences()
{
	CHPTilingData data;
	data.SetNumberSequences(10);
	CPPUNIT_ASSERT(data.GetNumberSequences() == 10);
}

void CHPTilingDataTest::test_AlgName()
{
	CHPTilingData data;
	data.SetAlgName(L"alg");
	CPPUNIT_ASSERT(data.GetAlgName() == L"alg");
}

void CHPTilingDataTest::test_AlgVersion()
{
	CHPTilingData data;
	data.SetAlgVersion(L"1.0");
	CPPUNIT_ASSERT(data.GetAlgVersion() == L"1.0");
}

void CHPTilingDataTest::test_AlgParams()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;
	CHPTilingData data;

	param.SetName(L"n1");
	param.SetValueText(L"v1");
	params.push_back(param);

	param.SetName(L"n2");
	param.SetValueText(L"v2");
	params.push_back(param);

	data.AddAlgParams(params);

	ParameterNameValueTypeList params_out = data.GetAlgParams();
	CPPUNIT_ASSERT(params_out.size() == 2);
	ParameterNameValueTypeList::iterator it = params_out.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"n1");
	CPPUNIT_ASSERT(param.GetValueText() == L"v1");
	++it;
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"n2");
	CPPUNIT_ASSERT(param.GetValueText() == L"v2");

}
