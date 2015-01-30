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
#include "calvin_files/parsers/test/DataGroupHeaderReaderTest.h"
//
#include "calvin_files/parsers/src/DataGroupHeaderReader.h"
#include "calvin_files/parsers/src/FileHeaderReader.h"
//
#include "util/Fs.h"
//

const std::string TEST_DATA_DAT_FILE = "../data/test.file.data_dat";

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataGroupHeaderReaderTest );

void DataGroupHeaderReaderTest::setUp()
{
  Fs::aptOpen(is, TEST_DATA_DAT_FILE.c_str(), std::ios::in | std::ios::binary);
}

void DataGroupHeaderReaderTest::tearDown()
{
	is.close();
	fh.Clear();
}

void DataGroupHeaderReaderTest::CreationTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	DataGroupHeaderReader dchReader;
	CPPUNIT_ASSERT(1);
}

void DataGroupHeaderReaderTest::ReadHeaderTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	// File position is at the start of the first DataGroupHeader
	DataGroupHeader dch;
	DataGroupHeaderReader dchReader;
	CPPUNIT_ASSERT_NO_THROW(dchReader.ReadHeader(is, dch));

	CPPUNIT_ASSERT(dch.GetName() == L"First Data Cube");
	CPPUNIT_ASSERT(dch.GetDataSetCnt() == 0);	// doesn't read in the DataSetHeader
	CPPUNIT_ASSERT(dch.GetDataSetPos() == 0x2d1);
	CPPUNIT_ASSERT(dch.GetNextGroupPos() == 0x468);
}

void DataGroupHeaderReaderTest::ReadOneDataGroupHeaderTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	// File position is at the start of the first DataGroupHeader
	DataGroupHeader dch;
	DataGroupHeaderReader dchReader;
	CPPUNIT_ASSERT_NO_THROW(dchReader.Read(is, dch));

	CPPUNIT_ASSERT(dch.GetName() == L"First Data Cube");
	CPPUNIT_ASSERT(dch.GetDataSetCnt() == 1);
	CPPUNIT_ASSERT(dch.GetDataSetPos() == 0x2d1);
	CPPUNIT_ASSERT(dch.GetNextGroupPos() == 0x468);

	// Check that the DataSetHeader information was read.
	DataSetHdrIt dphBegin, dphEnd;
	CPPUNIT_ASSERT_NO_THROW(dch.GetDataSetIterators(dphBegin,dphEnd));
	CPPUNIT_ASSERT(dphBegin != dphEnd);
	CPPUNIT_ASSERT(dphBegin->GetName() == L"acquired data");
	CPPUNIT_ASSERT(dphBegin->GetDataStartFilePos() != 0);
	CPPUNIT_ASSERT(dphBegin->GetNameValParamCnt() == 2);
	CPPUNIT_ASSERT(dphBegin->GetRowSize() == sizeof(u_int16_t));
	CPPUNIT_ASSERT(dphBegin->GetDataSize() == sizeof(u_int16_t)*100);
	CPPUNIT_ASSERT(dphBegin->GetRowCnt() == 100);
	CPPUNIT_ASSERT(dphBegin->GetColumnCnt() == 1);
	CPPUNIT_ASSERT(dphBegin->GetDataStartFilePos() == 0x3A0);

	// Check the DataSetHeader name value pairs
	ParameterNameValueTypeConstIt nvpBegin, nvpEnd;
	dphBegin->GetNameValIterators(nvpBegin, nvpEnd);
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == L"Scanner");
	CPPUNIT_ASSERT(nvpBegin->GetValueText() == L"M10");
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == L"Pixel Size");
	CPPUNIT_ASSERT(nvpBegin->GetValueFloat() == 0.051f);
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin == nvpEnd);

	// Check the data dataSet columns
	CPPUNIT_ASSERT(dphBegin->GetColumnInfo(0).GetColumnType() == UShortColType);
	CPPUNIT_ASSERT(dphBegin->GetColumnInfo(0).GetSize() == sizeof(u_int16_t));

	// Check that there are not any more DataSetHeaders
	++dphBegin;
	CPPUNIT_ASSERT(dphBegin == dphEnd);
}

void DataGroupHeaderReaderTest::ReadMinimumInfoForOneDataGroupHeaderTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	// File position is at the start of the first DataGroupHeader
	DataGroupHeader dch;
	DataGroupHeaderReader dchReader;
	CPPUNIT_ASSERT_NO_THROW(dchReader.ReadMinimumInfo(is, dch));

	CPPUNIT_ASSERT(dch.GetName() == L"First Data Cube");
	CPPUNIT_ASSERT(dch.GetDataSetCnt() == 1);
	CPPUNIT_ASSERT(dch.GetDataSetPos() == 0x2d1);
	CPPUNIT_ASSERT(dch.GetNextGroupPos() == 0x468);

	// Check that the minimum DataSetHeader information was read.
	DataSetHdrIt dphBegin, dphEnd;
	CPPUNIT_ASSERT_NO_THROW(dch.GetDataSetIterators(dphBegin,dphEnd));
	CPPUNIT_ASSERT(dphBegin->GetName() == L"acquired data");
	CPPUNIT_ASSERT(dphBegin->GetDataStartFilePos() != 0);
	CPPUNIT_ASSERT(dphBegin->GetNameValParamCnt() == 0);
	CPPUNIT_ASSERT(dphBegin->GetRowSize() == 0);
	CPPUNIT_ASSERT(dphBegin->GetDataSize() == 0);
	CPPUNIT_ASSERT(dphBegin->GetRowCnt() == 0);
	CPPUNIT_ASSERT(dphBegin->GetColumnCnt() == 0);
	CPPUNIT_ASSERT(dphBegin->GetDataStartFilePos() == 0x3A0);

	// Check that there are not any more DataSetHeaders
	++dphBegin;
	CPPUNIT_ASSERT(dphBegin == dphEnd);
}

void DataGroupHeaderReaderTest::ReadAllTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	// File position is at the start of the first DataGroupHeader
	DataGroupHeaderReader dchReader;
	CPPUNIT_ASSERT_NO_THROW(dchReader.ReadAll(is, fh, fhReader.GetDataGroupCnt()));

	// Check some header information
	CPPUNIT_ASSERT(fh.GetDataGroupCnt() == 1);
	CPPUNIT_ASSERT(fh.GetNumDataGroups() == 1);

	DataGroupHdrIt dchBegin, dchEnd;
	CPPUNIT_ASSERT_NO_THROW(fh.GetDataGroupIts(dchBegin, dchEnd));
	CPPUNIT_ASSERT(dchBegin->GetName() == L"First Data Cube");
	CPPUNIT_ASSERT(dchBegin->GetDataSetCnt() == 1);
	CPPUNIT_ASSERT(dchBegin->GetDataSetPos() == 0x2d1);
	CPPUNIT_ASSERT(dchBegin->GetNextGroupPos() == 0x468);

	// Check that the DataSetHeader information was read.
	DataSetHdrIt dphBegin, dphEnd;
	CPPUNIT_ASSERT_NO_THROW(dchBegin->GetDataSetIterators(dphBegin,dphEnd));
	CPPUNIT_ASSERT(dphBegin != dphEnd);
	CPPUNIT_ASSERT(dphBegin->GetName() == L"acquired data");
	CPPUNIT_ASSERT(dphBegin->GetDataStartFilePos() != 0);
	CPPUNIT_ASSERT(dphBegin->GetNameValParamCnt() == 2);
	CPPUNIT_ASSERT(dphBegin->GetRowSize() == sizeof(u_int16_t));
	CPPUNIT_ASSERT(dphBegin->GetDataSize() == sizeof(u_int16_t)*100);
	CPPUNIT_ASSERT(dphBegin->GetRowCnt() == 100);
	CPPUNIT_ASSERT(dphBegin->GetColumnCnt() == 1);
	CPPUNIT_ASSERT(dphBegin->GetDataStartFilePos() == 0x3A0);

	// Check the DataSetHeader name value pairs
	ParameterNameValueTypeConstIt nvpBegin, nvpEnd;
	dphBegin->GetNameValIterators(nvpBegin, nvpEnd);
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == L"Scanner");
	CPPUNIT_ASSERT(nvpBegin->GetValueText() == L"M10");
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == L"Pixel Size");
	CPPUNIT_ASSERT(nvpBegin->GetValueFloat() == 0.051f);
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin == nvpEnd);

	// Check the data dataSet columns
	CPPUNIT_ASSERT(dphBegin->GetColumnInfo(0).GetColumnType() == UShortColType);
	CPPUNIT_ASSERT(dphBegin->GetColumnInfo(0).GetSize() == sizeof(u_int16_t));

	// Check that there are not any more DataSetHeaders
	++dphBegin;
	CPPUNIT_ASSERT(dphBegin == dphEnd);

	// Check that there are not any more DataGroupHeaders
	++dchBegin;
	CPPUNIT_ASSERT(dchBegin == dchEnd);
}

void DataGroupHeaderReaderTest::ReadAllMinimumInfoTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	// File position is at the start of the first DataGroupHeader
	DataGroupHeaderReader dchReader;
	CPPUNIT_ASSERT_NO_THROW(dchReader.ReadAllMinimumInfo(is, fh, fhReader.GetDataGroupCnt()));

	// Check some header information
	CPPUNIT_ASSERT(fh.GetDataGroupCnt() == 1);
	CPPUNIT_ASSERT(fh.GetNumDataGroups() == 1);

	DataGroupHdrIt dchBegin, dchEnd;
	CPPUNIT_ASSERT_NO_THROW(fh.GetDataGroupIts(dchBegin, dchEnd));
	CPPUNIT_ASSERT(dchBegin->GetName() == L"First Data Cube");
	CPPUNIT_ASSERT(dchBegin->GetDataSetCnt() == 1);
	CPPUNIT_ASSERT(dchBegin->GetDataSetPos() == 0x2d1);
	CPPUNIT_ASSERT(dchBegin->GetNextGroupPos() == 0x468);

	// Check that the minimum DataSetHeader information was read.
	DataSetHdrIt dphBegin, dphEnd;
	CPPUNIT_ASSERT_NO_THROW(dchBegin->GetDataSetIterators(dphBegin,dphEnd));
	CPPUNIT_ASSERT(dphBegin->GetName() == L"acquired data");
	CPPUNIT_ASSERT(dphBegin->GetDataStartFilePos() != 0);
	CPPUNIT_ASSERT(dphBegin->GetNameValParamCnt() == 0);
	CPPUNIT_ASSERT(dphBegin->GetRowSize() == 0);
	CPPUNIT_ASSERT(dphBegin->GetDataSize() == 0);
	CPPUNIT_ASSERT(dphBegin->GetRowCnt() == 0);
	CPPUNIT_ASSERT(dphBegin->GetColumnCnt() == 0);
	CPPUNIT_ASSERT(dphBegin->GetDataStartFilePos() == 0x3A0);

	// Check that there are not any more DataSetHeaders
	++dphBegin;
	CPPUNIT_ASSERT(dphBegin == dphEnd);

	// Check that there are not any more DataGroupHeaders
	++dchBegin;
	CPPUNIT_ASSERT(dchBegin == dchEnd);
}
