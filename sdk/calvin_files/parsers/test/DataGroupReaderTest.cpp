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
#include "calvin_files/parsers/test/DataGroupReaderTest.h"
//
#include "calvin_files/parsers/src/DataGroupHeaderReader.h"
#include "calvin_files/parsers/src/DataGroupReader.h"
#include "calvin_files/parsers/src/FileHeaderReader.h"
//
#include "util/Fs.h"
//

#define TEST_DATA_DAT_FILE "../data/test.file.data_dat"

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataGroupReaderTest );

void DataGroupReaderTest::setUp()
{
  
  Fs::aptOpen(is,TEST_DATA_DAT_FILE, std::ios::in | std::ios::binary);

	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	DataGroupHeaderReader dchReader;
	dchReader.ReadAll(is, fh, fhReader.GetDataGroupCnt());
}

void DataGroupReaderTest::tearDown()
{
	is.close();
	fh.Clear();
}

void DataGroupReaderTest::CreationTest()
{
	DataGroupHeader& dch = fh.GetDataGroup(0);
	DataGroupReader reader(is, dch);
	CPPUNIT_ASSERT(1);
}

void DataGroupReaderTest::GetDataGroupNameTest()
{
	DataGroupHeader& dch = fh.GetDataGroup(0);
	DataGroupReader reader(is, dch);

	CPPUNIT_ASSERT(reader.GetDataGroupName() == L"First Data Cube");
}

void DataGroupReaderTest::GetDataSetReaderByIndexTest()
{
	DataGroupHeader& dch = fh.GetDataGroup(0);
	DataGroupReader reader(is, dch);

	DataSetReader dpReader = reader.GetDataSetReader(0);
	CPPUNIT_ASSERT(1);	// we got here

	// Now check some values
	int32_t rows = dch.GetDataSet(0).GetRowCnt();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dpReader.Read(value));
		CPPUNIT_ASSERT(value == expected);
	}
}

void DataGroupReaderTest::GetDataSetReaderByNameTest()
{
	DataGroupHeader& dch = fh.GetDataGroup(0);
	DataGroupReader reader(is, dch);

	DataSetReader dpReader = reader.GetDataSetReader(L"acquired data");
	CPPUNIT_ASSERT(1);	// we got here

	// Now check some values
	int32_t rows = dch.GetDataSet(0).GetRowCnt();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dpReader.Read(value));
		CPPUNIT_ASSERT(value == expected);
	}
}
