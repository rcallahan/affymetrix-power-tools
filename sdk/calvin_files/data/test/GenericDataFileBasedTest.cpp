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
#include "calvin_files/data/test/GenericDataFileBasedTest.h"
//
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;

const std::string TEST_DATA_DAT_FILE = "../data/test.file.data_dat";

CPPUNIT_TEST_SUITE_REGISTRATION( GenericDataFileBasedTest );

GenericDataFileBasedTest::GenericDataFileBasedTest()
{

}

GenericDataFileBasedTest::~GenericDataFileBasedTest()
{

}

void GenericDataFileBasedTest::setUp()
{
}

void GenericDataFileBasedTest::tearDown()
{
}

void GenericDataFileBasedTest::CreationTest()
{
	GenericData data;
	CPPUNIT_ASSERT(1);	// zzzz
}

/*
 * Test retrieving a DataGroup configured to read data using memory-mapping
 */
void GenericDataFileBasedTest::DataGroupMMTest()
{
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename(TEST_DATA_DAT_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data, GenericFileReader::ReadNoDataGroupHeader));

	u_int32_t addr = data.Header().GetFirstDataGroupFilePos();
	DataGroup dc = data.DataGroup(addr);

	// Check that we can get a DataSet and use it
	DataSet* dataPlane = dc.DataSet(0);
	CPPUNIT_ASSERT(dataPlane);

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

/*
 * Test retrieving a DataGroup configured to read data using small memory footprint fstream access.
 */
void GenericDataFileBasedTest::DataGroupFSreamTest()
{
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename(TEST_DATA_DAT_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data, GenericFileReader::ReadNoDataGroupHeader));
	data.UseMemoryMapping(false);
	data.LoadEntireDataSetHint(false);

	u_int32_t addr = data.Header().GetFirstDataGroupFilePos();
	DataGroup dc = data.DataGroup(addr);

	// Check that we can get a DataSet and use it
	DataSet* dataPlane = dc.DataSet(0);
	CPPUNIT_ASSERT(dataPlane);

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

void GenericDataFileBasedTest::DataGroupFStreamLoadEntireDataSetTest()
{
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename(TEST_DATA_DAT_FILE);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data, GenericFileReader::ReadNoDataGroupHeader));
	data.UseMemoryMapping(false);
	data.LoadEntireDataSetHint(true);

	u_int32_t addr = data.Header().GetFirstDataGroupFilePos();
	DataGroup dc = data.DataGroup(addr);

	// Check that we can get a DataSet and use it
	DataSet* dataPlane = dc.DataSet(0);
	CPPUNIT_ASSERT(dataPlane);

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
