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
#include "calvin_files/data/test/DataSetHeaderTest.h"
//
#include "calvin_files/data/src/ColumnInfo.h"
#include "calvin_files/data/src/DataSetHeader.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetHeaderTest );

void DataSetHeaderTest::setUp()
{
	header = new DataSetHeader();
}

void DataSetHeaderTest::tearDown()
{
	delete header;
}

void DataSetHeaderTest::testCreation()
{
	DataSetHeader hdr;
	CPPUNIT_ASSERT(1);
}

void DataSetHeaderTest::NameTest()
{
	std::wstring p1 = L"some_name";
	header->SetName(p1);
	std::wstring p2 = header->GetName();
	CPPUNIT_ASSERT(p1 == p2);
}

void DataSetHeaderTest::NameValPairCntTest()
{
	std::wstring p[] = {L"Uma",L"Thurman",L"Angelina",L"Jolie",L"Isabella ",L"Rossellini"};
	ParameterNameValueType t1;
	t1.SetName(p[0]);
	t1.SetValueText(p[1]);
	header->AddNameValParam(t1);

	ParameterNameValueType t2;
	t2.SetName(p[2]);
	t2.SetValueText(p[3]);
	header->AddNameValParam(t2);

	ParameterNameValueType t3;
	t3.SetName(p[4]);
	t3.SetValueText(p[5]);
	header->AddNameValParam(t3);
	CPPUNIT_ASSERT(header->GetNameValParamCnt() == 3);
}

void DataSetHeaderTest::NameValPairTest()
{
	std::wstring p[] = {L"Uma",L"Thurman",L"Angelina",L"Jolie",L"Isabella ",L"Rossellini"};
	ParameterNameValueType t1;
	t1.SetName(p[0]);
	t1.SetValueText(p[1]);
	header->AddNameValParam(t1);

	ParameterNameValueType t2;
	t2.SetName(p[2]);
	t2.SetValueText(p[3]);
	header->AddNameValParam(t2);

	ParameterNameValueType t3;
	t3.SetName(p[4]);
	t3.SetValueText(p[5]);
	header->AddNameValParam(t3);
	ParameterNameValueTypeConstIt begin;
	ParameterNameValueTypeConstIt end;
	header->GetNameValIterators(begin, end);
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetName() == p[0]);
		CPPUNIT_ASSERT(begin->GetValueText() == p[1]);
		begin++;
	}
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetName() == p[2]);
		CPPUNIT_ASSERT(begin->GetValueText() == p[3]);
		begin++;
	}
	if(begin != end)
	{
		CPPUNIT_ASSERT(begin->GetName() == p[4]);
		CPPUNIT_ASSERT(begin->GetValueText() == p[5]);
	}
}

void DataSetHeaderTest::ColumnNameTest()
{
	header->AddByteColumn(L"Fred");
	header->AddIntColumn(L"Wilma");
	header->AddShortColumn(L"");
	CPPUNIT_ASSERT(header->GetColumnInfo(0).GetName() == L"Fred");
	CPPUNIT_ASSERT(header->GetColumnInfo(1).GetName() == L"Wilma");
	CPPUNIT_ASSERT(header->GetColumnInfo(2).GetName() == L"");
}

void DataSetHeaderTest::ColumnValueTypeTest()
{
	header->AddByteColumn(L"");
	header->AddIntColumn(L"");
	header->AddShortColumn(L"");
	CPPUNIT_ASSERT(header->GetColumnInfo(0).GetColumnType() == ByteColType);
	CPPUNIT_ASSERT(header->GetColumnInfo(1).GetColumnType() == IntColType);
	CPPUNIT_ASSERT(header->GetColumnInfo(2).GetColumnType() == ShortColType);
}

void DataSetHeaderTest::ColumnSizeTest()
{
	header->AddByteColumn(L"");
	header->AddIntColumn(L"");
	header->AddShortColumn(L"");
	CPPUNIT_ASSERT(header->GetColumnInfo(0).GetSize() == sizeof(int8_t));
	CPPUNIT_ASSERT(header->GetColumnInfo(1).GetSize() == sizeof(int32_t));
	CPPUNIT_ASSERT(header->GetColumnInfo(2).GetSize() == sizeof(int16_t));
}

void DataSetHeaderTest::RowCntTest()
{
	int rows = 6;
	header->SetRowCnt(rows);
	CPPUNIT_ASSERT(header->GetRowCnt() == rows);
}

void DataSetHeaderTest::AssignmentTest()
{
	header->SetName(L"pixel data");
	header->AddUShortColumn(L"Intensity");
	header->SetRowCnt(101);
	ParameterNameValueType t1;
	t1.SetName(L"Scanner");
	t1.SetValueText(L"M10");
	header->AddNameValParam(t1);
	ParameterNameValueType t2;
	t2.SetName(L"Pixel Size");
	t2.SetValueText(L"0.051");
	header->AddNameValParam(t2);

	DataSetHeader assignee;
	assignee = *header;

	// Check that the assignment worked
	CPPUNIT_ASSERT(assignee.GetName() == header->GetName());
	CPPUNIT_ASSERT(assignee.GetColumnCnt() == header->GetColumnCnt());
	CPPUNIT_ASSERT(assignee.GetNameValParamCnt() == header->GetNameValParamCnt());
	CPPUNIT_ASSERT(assignee.GetRowCnt() == header->GetRowCnt());
	CPPUNIT_ASSERT(assignee.GetColumnInfo(0) == header->GetColumnInfo(0));

	ParameterNameValueTypeConstIt nvpAssigneeBegin, nvpAssigneeEnd;
	assignee.GetNameValIterators(nvpAssigneeBegin, nvpAssigneeEnd);
	ParameterNameValueTypeConstIt nvpHeaderBegin, nvpHeaderEnd;
	assignee.GetNameValIterators(nvpHeaderBegin, nvpHeaderEnd);
	CPPUNIT_ASSERT(*nvpAssigneeBegin == *nvpHeaderBegin);
	++nvpAssigneeBegin;
	++nvpHeaderBegin;
	CPPUNIT_ASSERT(*nvpAssigneeBegin == *nvpHeaderBegin);
	++nvpAssigneeBegin;
	++nvpHeaderBegin;
	CPPUNIT_ASSERT(nvpAssigneeBegin == nvpAssigneeEnd);
	CPPUNIT_ASSERT(nvpHeaderBegin == nvpHeaderEnd);

	// TBD: no test of row offset array.
}


void DataSetHeaderTest::CopyCTorTest()
{
	header->SetName(L"pixel data");
	header->AddUShortColumn(L"Intensity");
	header->SetRowCnt(101);
	ParameterNameValueType t1;
	t1.SetName(L"Scanner");
	t1.SetValueText(L"M10");
	header->AddNameValParam(t1);
	ParameterNameValueType t2;
	t2.SetName(L"Pixel Size");
	t2.SetValueText(L"0.051");
	header->AddNameValParam(t2);

	DataSetHeader assignee(*header);

	// Check that the copy construction worked
	CPPUNIT_ASSERT(assignee.GetName() == header->GetName());
	CPPUNIT_ASSERT(assignee.GetColumnCnt() == header->GetColumnCnt());
	CPPUNIT_ASSERT(assignee.GetNameValParamCnt() == header->GetNameValParamCnt());
	CPPUNIT_ASSERT(assignee.GetRowCnt() == header->GetRowCnt());
	CPPUNIT_ASSERT(assignee.GetColumnInfo(0) == header->GetColumnInfo(0));
	CPPUNIT_ASSERT(assignee.GetColumnCnt() == header->GetColumnCnt());

	ParameterNameValueTypeConstIt nvpAssigneeBegin, nvpAssigneeEnd;
	assignee.GetNameValIterators(nvpAssigneeBegin, nvpAssigneeEnd);
	ParameterNameValueTypeConstIt nvpHeaderBegin, nvpHeaderEnd;
	assignee.GetNameValIterators(nvpHeaderBegin, nvpHeaderEnd);
	CPPUNIT_ASSERT(*nvpAssigneeBegin == *nvpHeaderBegin);
	++nvpAssigneeBegin;
	++nvpHeaderBegin;
	CPPUNIT_ASSERT(*nvpAssigneeBegin == *nvpHeaderBegin);
	++nvpAssigneeBegin;
	++nvpHeaderBegin;
	CPPUNIT_ASSERT(nvpAssigneeBegin == nvpAssigneeEnd);
	CPPUNIT_ASSERT(nvpHeaderBegin == nvpHeaderEnd);

	// TBD: no test of row offset array.
}
void DataSetHeaderTest::DataStartFilePosTest()
{
	CPPUNIT_ASSERT(header->GetDataStartFilePos() == 0);
	CPPUNIT_ASSERT_NO_THROW(header->SetDataStartFilePos(12));
	CPPUNIT_ASSERT(header->GetDataStartFilePos() == 12);
}

void DataSetHeaderTest::HeaderStartFilePosTest()
{
	CPPUNIT_ASSERT(header->GetHeaderStartFilePos() == 0);
	CPPUNIT_ASSERT_NO_THROW(header->SetHeaderStartFilePos(61209));
	CPPUNIT_ASSERT(header->GetHeaderStartFilePos() == 61209);
}
