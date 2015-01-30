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
#include "calvin_files/data/test/DataSetTest.h"
//
#include "calvin_files/data/src/DataException.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;

#define TEST_DATA_DAT_FILE "../data/test.file.data_dat"

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetTest );

DataSetTest::DataSetTest()
{
	data = 0;
	dataPlane = 0;
#ifdef _MSC_VER
	fileHandle = INVALID_HANDLE_VALUE;
	fileMapHandle = NULL;
#endif
}

DataSetTest::~DataSetTest()
{

}

void DataSetTest::setUp()
{
	data = new GenericData;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);
	reader.ReadHeader(*data);
	DataGroupHdrIt dchBegin, dchEnd;
	data->Header().GetDataGroupIts(dchBegin, dchEnd);
	CPPUNIT_ASSERT(dchBegin != dchEnd);
	DataSetHdrIt dphBegin, dphEnd;
	dchBegin->GetDataSetIterators(dphBegin, dphEnd);
	void* handle = MapFile(name);
	dataPlane = new DataSet(data->Header().GetFilename(), *dphBegin, handle);
}

void DataSetTest::tearDown()
{
	UnmapFile();
	if (dataPlane)
		dataPlane->Delete();
	delete data;
}

void DataSetTest::testCreation()
{
	DataSetHeader header;
	DataSet* dc = new DataSet(TEST_DATA_DAT_FILE, header,0);
	CPPUNIT_ASSERT(dc != 0);
	dc->Delete();
}

void DataSetTest::testmethod_Delete()
{
	DataSetHeader header;
	DataSet* dc = new DataSet(TEST_DATA_DAT_FILE, header,0);
	CPPUNIT_ASSERT(dc != 0);
	CPPUNIT_ASSERT_NO_THROW(dc->Delete());
	// TBD: How to assert that Delete destroyed dc?
}

void DataSetTest::testmethod_Open()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	// TBD: How to assert that the DataGroup is open?
}

void DataSetTest::testmethod_Close()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
	// TBD: How to assert that the DataGroup is closed?
}

void DataSetTest::testproperty_BytesPerRow()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(dataPlane->BytesPerRow() == 2);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testproperty_Rows()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(dataPlane->Rows() == 100);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testproperty_Cols()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(dataPlane->Cols() == 1);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testproperty_Header()
{
	// Get the DataSetHeader
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	const DataSetHeader& header = dataPlane->Header();

	// Get the DataSetHeader from the FileHeader to check against
	DataGroupHdrIt dchBegin, dchEnd;
	data->Header().GetDataGroupIts(dchBegin, dchEnd);
	CPPUNIT_ASSERT(dchBegin != dchEnd);
	DataSetHdrIt dphBegin, dphEnd;
	dchBegin->GetDataSetIterators(dphBegin, dphEnd);

	CPPUNIT_ASSERT(header.GetName() == dphBegin->GetName());
	CPPUNIT_ASSERT(header.GetDataSize() == dphBegin->GetDataSize());
	CPPUNIT_ASSERT(header.GetRowSize() == dphBegin->GetRowSize());
	CPPUNIT_ASSERT(header.GetNameValParamCnt() == dphBegin->GetNameValParamCnt());
	CPPUNIT_ASSERT(header.GetColumnInfo(0) == dphBegin->GetColumnInfo(0));
	CPPUNIT_ASSERT(header.GetColumnCnt() == dphBegin->GetColumnCnt());
	CPPUNIT_ASSERT(header.GetRowCnt() == dphBegin->GetRowCnt());
}

void DataSetTest::testmethod_ReadUInt16_Single()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 0, value));
		CPPUNIT_ASSERT(value == expected);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_ReadUInt16_Vector()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	Uint16Vector values;
	CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(0, 0, rows, values));

	// Check the values
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		CPPUNIT_ASSERT(values[row] == expected);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_ReadUInt16_Vector_TooManyRows()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	Uint16Vector values;
	CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(0, 0, rows*2, values));
	CPPUNIT_ASSERT(values.size() == rows);

	// Check the values
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		CPPUNIT_ASSERT(values[row] == expected);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_ReadUInt16_Vector_InChunks()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	int32_t chunks = 7;
	int32_t rowsPerChunk = rows/chunks;

	for (int32_t chunk = 0; chunk < chunks; ++chunk)
	{
		int32_t startRow = chunk*rowsPerChunk;
		Uint16Vector values;
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(0, startRow, rowsPerChunk, values));

		int32_t rowsInVector = rowsPerChunk;
		if (dataPlane->Rows() < startRow+rowsPerChunk)
			rowsInVector = dataPlane->Rows() - startRow;

		CPPUNIT_ASSERT(values.size() == rowsInVector);

		// Check the values
		for( int32_t row = startRow; row < rowsInVector+startRow; ++row )
		{
			u_int16_t expected = (u_int16_t)(row*10+row);
			CPPUNIT_ASSERT(values[row-startRow] == expected);
		}
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_ReadUInt16_Raw()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	u_int16_t* values = new u_int16_t[rows];
	CPPUNIT_ASSERT_NO_THROW(dataPlane->GetDataRaw(0, 0, rows, values));

	// Check the values
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		CPPUNIT_ASSERT(values[row] == expected);
	}

	delete [] values;

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_ReadUInt16_Raw_TooManyRows()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	u_int16_t* values = new u_int16_t[rows*2];
	u_int32_t rowsRead = 0;
	CPPUNIT_ASSERT_NO_THROW(rowsRead = dataPlane->GetDataRaw(0, 0, rows*2, values));
	CPPUNIT_ASSERT(rowsRead == rows);

	// Check the values
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		CPPUNIT_ASSERT(values[row] == expected);
	}

	delete [] values;

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_ReadUInt16_Raw_InChunks()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	int32_t chunks = 7;
	int32_t rowsPerChunk = rows/chunks;

	for (int32_t chunk = 0; chunk < chunks; ++chunk)
	{
		int32_t startRow = chunk*rowsPerChunk;
		u_int16_t* values = new u_int16_t[rows*2];
		u_int32_t rowsRead = 0;
		CPPUNIT_ASSERT_NO_THROW(rowsRead = dataPlane->GetDataRaw(0, startRow, rowsPerChunk, values));

		int32_t rowsInArray = rowsPerChunk;
		if (dataPlane->Rows() < startRow+rowsPerChunk)
			rowsInArray = dataPlane->Rows() - startRow;

		CPPUNIT_ASSERT(rowsRead == rowsInArray);

		// Check the values
		for( int32_t row = startRow; row < rowsInArray+startRow; ++row )
		{
			u_int16_t expected = (u_int16_t)(row*10+row);
			CPPUNIT_ASSERT(values[row-startRow] == expected);
		}

		delete [] values;
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_CheckRowColumnAndType_DataSetNotOpenException()
{
	CPPUNIT_ASSERT_THROW(dataPlane->CheckRowColumnAndType(0, 0, ByteColType), affymetrix_calvin_exceptions::DataSetNotOpenException);
}

void DataSetTest::testmethod_CheckRowColumnAndType_ColumnOutOfBoundsException()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_THROW(dataPlane->CheckRowColumnAndType(0, 1, UShortColType), affymetrix_calvin_exceptions::ColumnIndexOutOfBoundsException);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_CheckRowColumnAndType_RowOutofBoundsException()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_THROW(dataPlane->CheckRowColumnAndType(100, 0, UShortColType), affymetrix_calvin_exceptions::RowIndexOutOfBoundsException);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_CheckRowColumnAndType_UnexpectedColumnTypeException()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_THROW(dataPlane->CheckRowColumnAndType(99, 0, UIntColType), affymetrix_calvin_exceptions::UnexpectedColumnTypeException);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest::testmethod_ReadUInt16_Single_DataSetNotOpenException()
{ 
	u_int16_t value;
	CPPUNIT_ASSERT_THROW(dataPlane->GetData(0, 0, value), affymetrix_calvin_exceptions::DataSetNotOpenException);
}

void DataSetTest::testmethod_ReadUInt16_Vector_DataSetNotOpenException()
{ 
	Uint16Vector values;
	CPPUNIT_ASSERT_THROW(dataPlane->GetData(0, 0, 1, values), affymetrix_calvin_exceptions::DataSetNotOpenException);
}

void DataSetTest::testmethod_ReadUInt16_Raw_DataSetNotOpenException()
{ 
	u_int16_t* values = new u_int16_t[1];
	CPPUNIT_ASSERT_THROW(dataPlane->GetDataRaw(0, 0, 1, values), affymetrix_calvin_exceptions::DataSetNotOpenException);
	delete [] values;
}

void* DataSetTest::MapFile(const std::string& fileName)
{
#ifdef _MSC_VER
	if (fileHandle == INVALID_HANDLE_VALUE)
	{
		// Create the file.
		fileHandle = CreateFile(fileName.c_str(), GENERIC_READ, FILE_SHARE_READ,
				NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
		if (fileHandle == INVALID_HANDLE_VALUE)
		{
			return 0;
		}
	}
	if (fileMapHandle == NULL)
	{
		// Use the current size of the file.
		DWORD dwSizeHigh = 0;
		DWORD dwSizeLow  = 0;
		fileMapHandle = CreateFileMapping(fileHandle, NULL, PAGE_READONLY, dwSizeHigh, dwSizeLow, NULL);
		if (fileMapHandle == NULL)
			return 0;
	}
	return fileMapHandle;
#else
	return 0;
#endif
}

void DataSetTest::UnmapFile()
{
#ifdef _MSC_VER
	if (fileHandle != INVALID_HANDLE_VALUE)
	{
		if (fileMapHandle != NULL)
		{
			CloseHandle(fileMapHandle);
			fileMapHandle = NULL;
		}
		CloseHandle (fileHandle);
		fileHandle = INVALID_HANDLE_VALUE;
	}
#endif
}
