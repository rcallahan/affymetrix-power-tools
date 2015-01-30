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
#include "calvin_files/data/test/DataGroupHeaderTest.h"
//
#include "calvin_files/data/src/ColumnInfo.h"
#include "calvin_files/data/src/DataGroupHeader.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataGroupHeaderTest );

void DataGroupHeaderTest::setUp()
{
	header = new DataGroupHeader(L"dataGroup");
}

void DataGroupHeaderTest::tearDown()
{
	delete header;
}

void DataGroupHeaderTest::testCreation()
{
	DataGroupHeader hdr(L"dataGroup");
	CPPUNIT_ASSERT(1);
}

void DataGroupHeaderTest::NameTest()
{
	std::wstring p1 = L"some_name";
	header->SetName(p1);
	std::wstring p2 = header->GetName();
	CPPUNIT_ASSERT(p1 == p2);
}

void DataGroupHeaderTest::DataSetCntTest()
{
	CPPUNIT_ASSERT(1);
}

void DataGroupHeaderTest::DataSetTest()
{
	CPPUNIT_ASSERT(1);
}

void DataGroupHeaderTest::FindDataSetHeaderTest()
{
	// Create DataSetHeaders
	DataSetHeader dph1;
	dph1.SetName(L"pixel intensity");
	ParameterNameValueType param;
	param.SetName(L"Scanner");
	param.SetValueText(L"M10");
	dph1.AddNameValParam(param);
	dph1.AddUShortColumn(L"");
	dph1.SetRowCnt(1);

	DataSetHeader dph2;
	dph2.SetName(L"grid coordinates");
	param.SetName(L"Corner Pattern");
	param.SetValueText(L"Checkerboard");
	dph2.AddNameValParam(param);
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");

	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");

	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");

	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.AddUShortColumn(L"");
	dph2.SetRowCnt(1);

	header->AddDataSetHdr(dph1);
	header->AddDataSetHdr(dph2);

	const DataSetHeader* dph = header->FindDataSetHeader(L"none");
	CPPUNIT_ASSERT(0 == dph);
	dph = header->FindDataSetHeader(dph1.GetName());
	CPPUNIT_ASSERT(0 != dph);
	CPPUNIT_ASSERT(dph->GetName() == dph1.GetName());
	dph = header->FindDataSetHeader(dph2.GetName());
	CPPUNIT_ASSERT(0 != dph);
	CPPUNIT_ASSERT(dph->GetName() == dph2.GetName());
}
