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
#include "calvin_files/data/test/DataSetTest_FStream.h"
//
#include "calvin_files/data/src/DataException.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include "util/Fs.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;

#define TEST_DATA_DAT_FILE "../data/test.file.data_dat"

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetTest_FStream );

DataSetTest_FStream::DataSetTest_FStream()
{
	data = 0;
	dataPlane = 0;
}

DataSetTest_FStream::~DataSetTest_FStream()
{

}

void DataSetTest_FStream::setUp()
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
        Fs::aptOpen(ifs, name, std::ios::in | std::ios::binary);
	dataPlane = new DataSet(data->Header().GetFilename(), *dphBegin, ifs);
}

void DataSetTest_FStream::tearDown()
{
	ifs.close();
	if (dataPlane)
		dataPlane->Delete();
	delete data;
}

void DataSetTest_FStream::testCreation()
{
	DataSetHeader header;
	DataSet* dc = new DataSet(TEST_DATA_DAT_FILE, header,ifs);
	CPPUNIT_ASSERT(dc != 0);
	dc->Delete();
}

void DataSetTest_FStream::testmethod_Delete()
{
	DataSetHeader header;
	DataSet* dc = new DataSet(TEST_DATA_DAT_FILE, header,ifs);
	CPPUNIT_ASSERT(dc != 0);
	CPPUNIT_ASSERT_NO_THROW(dc->Delete());
	// TBD: How to assert that Delete destroyed dc?
}

void DataSetTest_FStream::testmethod_Open()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	// TBD: How to assert that the DataGroup is open?
}

void DataSetTest_FStream::testmethod_Close()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
	// TBD: How to assert that the DataGroup is closed?
}

void DataSetTest_FStream::testproperty_BytesPerRow()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(dataPlane->BytesPerRow() == 2);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest_FStream::testproperty_Rows()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(dataPlane->Rows() == 100);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest_FStream::testproperty_Cols()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT(dataPlane->Cols() == 1);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest_FStream::testproperty_Header()
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

void DataSetTest_FStream::testmethod_ReadUInt16_Single()
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

void DataSetTest_FStream::testmethod_ReadUInt16_Vector()
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

void DataSetTest_FStream::testmethod_ReadUInt16_Vector_TooManyRows()
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

void DataSetTest_FStream::testmethod_ReadUInt16_Vector_InChunks()
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

void DataSetTest_FStream::testmethod_ReadUInt16_Raw()
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

void DataSetTest_FStream::testmethod_ReadUInt16_Raw_TooManyRows()
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

void DataSetTest_FStream::testmethod_ReadUInt16_Raw_InChunks()
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

void DataSetTest_FStream::testmethod_CheckRowColumnAndType_DataSetNotOpenException()
{
	CPPUNIT_ASSERT_THROW(dataPlane->CheckRowColumnAndType(0, 0, ByteColType), affymetrix_calvin_exceptions::DataSetNotOpenException);
}

void DataSetTest_FStream::testmethod_CheckRowColumnAndType_ColumnOutOfBoundsException()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_THROW(dataPlane->CheckRowColumnAndType(0, 1, UShortColType), affymetrix_calvin_exceptions::ColumnIndexOutOfBoundsException);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest_FStream::testmethod_CheckRowColumnAndType_RowOutofBoundsException()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_THROW(dataPlane->CheckRowColumnAndType(100, 0, UShortColType), affymetrix_calvin_exceptions::RowIndexOutOfBoundsException);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest_FStream::testmethod_CheckRowColumnAndType_UnexpectedColumnTypeException()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_THROW(dataPlane->CheckRowColumnAndType(99, 0, UIntColType), affymetrix_calvin_exceptions::UnexpectedColumnTypeException);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest_FStream::testmethod_ReadUInt16_Single_DataSetNotOpenException()
{ 
	u_int16_t value;
	CPPUNIT_ASSERT_THROW(dataPlane->GetData(0, 0, value), affymetrix_calvin_exceptions::DataSetNotOpenException);
}

void DataSetTest_FStream::testmethod_ReadUInt16_Vector_DataSetNotOpenException()
{ 
	Uint16Vector values;
	CPPUNIT_ASSERT_THROW(dataPlane->GetData(0, 0, 1, values), affymetrix_calvin_exceptions::DataSetNotOpenException);
}

void DataSetTest_FStream::testmethod_ReadUInt16_Raw_DataSetNotOpenException()
{ 
	u_int16_t* values = new u_int16_t[1];
	CPPUNIT_ASSERT_THROW(dataPlane->GetDataRaw(0, 0, 1, values), affymetrix_calvin_exceptions::DataSetNotOpenException);
	delete [] values;
}
