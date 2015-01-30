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

#include "calvin_files/writers/test/DataSetUpdaterTest.h"
//
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/DataSetUpdater.h"
#include "calvin_files/writers/src/GenericFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetUpdaterTest );

#define TEST_FILE "data_file.datasetupdater_test"

void DataSetUpdaterTest::setUp()
{
}

void DataSetUpdaterTest::CreateReferenceFile()
{
	// Create the data object.
	GenericData data;
	data.Header().SetFilename(TEST_FILE);

	DataGroupHeader dgHdr(L"data");
	DataSetHeader dsHdr;

	dsHdr.SetName(L"test");
	dsHdr.AddAsciiColumn(L"string8", 64);
	dsHdr.AddUnicodeColumn(L"string16", 64);
	dsHdr.AddByteColumn(L"int8");
	dsHdr.AddShortColumn(L"int16");
	dsHdr.AddIntColumn(L"int32");
	dsHdr.AddUByteColumn(L"uint8");
	dsHdr.AddUShortColumn(L"uint16");
	dsHdr.AddUIntColumn(L"uint32");
	dsHdr.AddFloatColumn(L"float");
	dsHdr.SetRowCnt(3);

	dgHdr.AddDataSetHdr(dsHdr);

	data.Header().AddDataGroupHdr(dgHdr);
	data.Header().GetGenericDataHdr()->SetFileTypeId("affymetrix.test.data");



	// Write the data object to the file
	GenericFileWriter writer(&data.Header());
	writer.WriteHeader();

	DataGroupWriter &dataGroupWriter = writer.GetDataGroupWriter(0);
	dataGroupWriter.WriteHeader();

	//int iSet=0;
	DataSetWriterIt beginSet;
	DataSetWriterIt endSet;
	dataGroupWriter.GetDataSetWriters(beginSet, endSet);

	beginSet->WriteHeader();

	int32_t dataSetSz = beginSet->GetDataSetSize();
	int32_t offset = writer.GetFilePos();
	writer.SeekFromCurrentPos(dataSetSz + 1);
	beginSet->UpdateNextDataSetOffset();

	dataGroupWriter.UpdateNextDataGroupPos();
	
	writer.SeekFromBeginPos(offset);


	// Write the data.
	beginSet->Write("first_row", 64);
	beginSet->Write(L"first_row", 64);
	beginSet->Write((int8_t) 1);
	beginSet->Write((int16_t) 1);
	beginSet->Write((int32_t) 1);
	beginSet->Write((u_int8_t) 1);
	beginSet->Write((u_int16_t) 1);
	beginSet->Write((u_int32_t) 1);
	beginSet->Write(1.0f);

	beginSet->Write("second_row", 64);
	beginSet->Write(L"second_row", 64);
	beginSet->Write((int8_t) 2);
	beginSet->Write((int16_t) 2);
	beginSet->Write((int32_t) 2);
	beginSet->Write((u_int8_t) 2);
	beginSet->Write((u_int16_t) 2);
	beginSet->Write((u_int32_t) 2);
	beginSet->Write(2.0f);

	beginSet->Write("third_row", 64);
	beginSet->Write(L"third_row", 64);
	beginSet->Write((int8_t) 3);
	beginSet->Write((int16_t) 3);
	beginSet->Write((int32_t) 3);
	beginSet->Write((u_int8_t) 3);
	beginSet->Write((u_int16_t) 3);
	beginSet->Write((u_int32_t) 3);
	beginSet->Write(3.0f);
}

void DataSetUpdaterTest::tearDown()
{
}

void DataSetUpdaterTest::testUpdateFail()
{
	DataSetUpdater u;
	CPPUNIT_ASSERT_THROW(u.Initialize("no_file"), FileNotFoundException);
	CPPUNIT_ASSERT_THROW(u.Initialize("../../parsers/data/test.file.full_array_file"), InvalidFileTypeException);
}

typedef struct _RowDataType
{
	string str;
	wstring wstr;
	int8_t b;
	int16_t s;
	int32_t i;
	u_int8_t ub;
	u_int16_t us;
	u_int32_t ui;
	float f;
} RowDataType;

void GetRow(int row, RowDataType &d)
{
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename(TEST_FILE);
	reader.Open(data);

	DataSet *set = data.DataSet(0, 0);
	set->Open();

	set->GetData(row, 0, d.str);
	set->GetData(row, 1, d.wstr);
	set->GetData(row, 2, d.b);
	set->GetData(row, 3, d.s);
	set->GetData(row, 4, d.i);
	set->GetData(row, 5, d.ub);
	set->GetData(row, 6, d.us);
	set->GetData(row, 7, d.ui);
	set->GetData(row, 8, d.f);

	set->Close();
	set->Delete();
}

void DataSetUpdaterTest::testUpdateInt8()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	int8_t newval = 12;
	upd.Update(0, 0, 2, 2, newval);
	newval = 22;
	upd.Update(0, 0, 1, 2, newval);

	RowDataType data;
	GetRow(2, data);
	CPPUNIT_ASSERT(data.b == 12);
	GetRow(1, data);
	CPPUNIT_ASSERT(data.b == 22);
}

void DataSetUpdaterTest::testUpdateInt16()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	int16_t newval = 12;
	upd.Update(0, 0, 2, 3, newval);
	newval = 22;
	upd.Update(0, 0, 1, 3, newval);

	RowDataType data;
	GetRow(2, data);
	CPPUNIT_ASSERT(data.s == 12);
	GetRow(1, data);
	CPPUNIT_ASSERT(data.s == 22);
}

void DataSetUpdaterTest::testUpdateInt32()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	int32_t newval = 12;
	upd.Update(0, 0, 2, 4, newval);
	newval = 22;
	upd.Update(0, 0, 1, 4, newval);

	RowDataType data;
	GetRow(2, data);
	CPPUNIT_ASSERT(data.i == 12);
	GetRow(1, data);
	CPPUNIT_ASSERT(data.i == 22);
}

void DataSetUpdaterTest::testUpdateUInt8()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	u_int8_t newval = 12;
	int row=2;
	int col=5;
	upd.Update(0, 0, row, col, newval);

	RowDataType data;
	GetRow(row, data);
	CPPUNIT_ASSERT(data.ub == 12);
}

void DataSetUpdaterTest::testUpdateUInt16()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	u_int16_t newval = 12;
	int row=2;
	int col=6;
	upd.Update(0, 0, row, col, newval);

	RowDataType data;
	GetRow(row, data);
	CPPUNIT_ASSERT(data.us == 12);
}

void DataSetUpdaterTest::testUpdateUInt32()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	u_int32_t newval = 12;
	int row=2;
	int col=7;
	upd.Update(0, 0, row, col, newval);

	RowDataType data;
	GetRow(row, data);
	CPPUNIT_ASSERT(data.ui == 12);
}

void DataSetUpdaterTest::testUpdateFloat()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	upd.Update(0, 0, 2, 8, 12.0f);
	upd.Update(0, 0, 0, 8, 22.0f);

	RowDataType data;
	GetRow(0, data);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.f, 22.0f, 0.000001f);
	GetRow(2, data);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.f, 12.0f, 0.000001f);
}

void DataSetUpdaterTest::testUpdateString8()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	upd.Update(0,0,1,0,"newval");
	upd.Update(0,0,0,0,"another_value");
	upd.Update(0,0,2,0,"1");
	RowDataType data;
	GetRow(1, data);
	CPPUNIT_ASSERT(data.str == "newval");
	GetRow(0, data);
	CPPUNIT_ASSERT(data.str == "another_value");
	GetRow(2, data);
	CPPUNIT_ASSERT(data.str == "1");
}

void DataSetUpdaterTest::testUpdateString16()
{
	CreateReferenceFile();
	DataSetUpdater upd;
	upd.Initialize(TEST_FILE);
	upd.Update(0,0,1,1,L"newval");
	upd.Update(0,0,0,1,L"one");
	upd.Update(0,0,2,1,L"two");
	RowDataType data;
	GetRow(0, data);
	CPPUNIT_ASSERT(data.wstr == L"one");
	GetRow(1, data);
	CPPUNIT_ASSERT(data.wstr == L"newval");
	GetRow(2, data);
	CPPUNIT_ASSERT(data.wstr == L"two");
}
