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
#include "calvin_files/data/test/DataGroupTest.h"
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

const std::string TEST_DATA_DAT_FILE = "../data/test.file.data_dat";

CPPUNIT_TEST_SUITE_REGISTRATION( DataGroupTest );

DataGroupTest::DataGroupTest()
{
	data = 0;
	dc = 0;
#ifdef _MSC_VER
	fileHandle = INVALID_HANDLE_VALUE;
	fileMapHandle = NULL;
#endif
}

DataGroupTest::~DataGroupTest()
{

}

void DataGroupTest::setUp()
{
	data = new GenericData;
	GenericFileReader reader;
	reader.SetFilename(TEST_DATA_DAT_FILE);
	reader.ReadHeader(*data);
	DataGroupHdrIt dchBegin, dchEnd;
	data->Header().GetDataGroupIts(dchBegin, dchEnd);
	void* handle = MapFile(TEST_DATA_DAT_FILE);
	dc = new DataGroup(TEST_DATA_DAT_FILE, *dchBegin, handle);
}

void DataGroupTest::tearDown()
{
	delete data;
	delete dc;
	UnmapFile();
}

void DataGroupTest::CreationTest()
{
	CPPUNIT_ASSERT(dc != 0);
}

void DataGroupTest::HeaderTest()
{
	CPPUNIT_ASSERT(dc->Header().GetName() == L"First Data Cube");
}

void DataGroupTest::DataSetByIndexTest()
{
	DataSet* dataPlane = dc->DataSet(0);
	CPPUNIT_ASSERT(dataPlane);

	// Check the DataSet data
	dataPlane->Open();
	CPPUNIT_ASSERT(dataPlane->BytesPerRow() == 2);
	CPPUNIT_ASSERT(dataPlane->Rows() == 100);
	CPPUNIT_ASSERT(dataPlane->Cols() == 1);

	int32_t rows = dataPlane->Rows();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 0, value));
		CPPUNIT_ASSERT(value == expected);
	}

	dataPlane->Close();	// optional - Delete will call this
	dataPlane->Delete();
}

void DataGroupTest::DataSetByNameTest()
{
	DataSet* dataPlane = dc->DataSet(L"acquired data");
	CPPUNIT_ASSERT(dataPlane);

	// Check the DataSet data
	dataPlane->Open();
	CPPUNIT_ASSERT(dataPlane->BytesPerRow() == 2);
	CPPUNIT_ASSERT(dataPlane->Rows() == 100);
	CPPUNIT_ASSERT(dataPlane->Cols() == 1);

	int32_t rows = dataPlane->Rows();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 0, value));
		CPPUNIT_ASSERT(value == expected);
	}

	dataPlane->Close();	// optional - Delete will call this
	dataPlane->Delete();
}

void* DataGroupTest::MapFile(const std::string& fileName)
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

void DataGroupTest::UnmapFile()
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
