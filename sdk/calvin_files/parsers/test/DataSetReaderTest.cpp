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
#include "calvin_files/parsers/test/DataSetReaderTest.h"
//
#include "calvin_files/parsers/src/DataGroupHeaderReader.h"
#include "calvin_files/parsers/src/DataSetReader.h"
#include "calvin_files/parsers/src/FileHeaderReader.h"
//
#include "util/Fs.h"
//

#define TEST_DATA_DAT_FILE "../data/test.file.data_dat"

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetReaderTest );

void DataSetReaderTest::setUp()
{
  Fs::aptOpen(is, TEST_DATA_DAT_FILE, std::ios::in | std::ios::binary);

	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	DataGroupHeaderReader dchReader;
	dchReader.ReadAll(is, fh, fhReader.GetDataGroupCnt());
}

void DataSetReaderTest::tearDown()
{
	is.close();
	fh.Clear();
}

void DataSetReaderTest::CreationTest()
{
	DataSetHeader& dph = fh.GetDataGroup(0).GetDataSet(0);
	DataSetReader reader(is, dph);
	CPPUNIT_ASSERT(1);
}

void DataSetReaderTest::ReadDataSetTest()
{
	DataSetHeader& dph = fh.GetDataGroup(0).GetDataSet(0);
	DataSetReader reader(is, dph);

	CPPUNIT_ASSERT(reader.GetDataSetName() == L"acquired data");

	int32_t rows = dph.GetRowCnt();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(reader.Read(value));
		CPPUNIT_ASSERT(value == expected);
	}
}
