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
#include "calvin_files/data/test/DataGroupTest_FStreamLoadEntireDataSet.h"
//
#include "calvin_files/data/src/DataException.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include "util/Fs.h"
//
#include <cstring>
#include <fstream>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;

const std::string TEST_DATA_DAT_FILE = "../data/test.file.data_dat";

CPPUNIT_TEST_SUITE_REGISTRATION( DataGroupTest_FStreamLoadEntireDataSet );

DataGroupTest_FStreamLoadEntireDataSet::DataGroupTest_FStreamLoadEntireDataSet()
{
	data = 0;
	dc = 0;
}

DataGroupTest_FStreamLoadEntireDataSet::~DataGroupTest_FStreamLoadEntireDataSet()
{

}

void DataGroupTest_FStreamLoadEntireDataSet::setUp()
{
	data = new GenericData;
	GenericFileReader reader;
	reader.SetFilename(TEST_DATA_DAT_FILE);
	reader.ReadHeader(*data);
	DataGroupHdrIt dchBegin, dchEnd;
	data->Header().GetDataGroupIts(dchBegin, dchEnd);

        Fs::aptOpen(ifs, TEST_DATA_DAT_FILE, std::ios::in | std::ios::binary);
	dc = new DataGroup(TEST_DATA_DAT_FILE, *dchBegin, ifs, true);
}

void DataGroupTest_FStreamLoadEntireDataSet::tearDown()
{
	ifs.close();
	delete data;
	delete dc;
}

void DataGroupTest_FStreamLoadEntireDataSet::CreationTest()
{
	CPPUNIT_ASSERT(dc != 0);
}

void DataGroupTest_FStreamLoadEntireDataSet::HeaderTest()
{
	CPPUNIT_ASSERT(dc->Header().GetName() == L"First Data Cube");
}

void DataGroupTest_FStreamLoadEntireDataSet::DataSetByIndexTest()
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

void DataGroupTest_FStreamLoadEntireDataSet::DataSetByNameTest()
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
