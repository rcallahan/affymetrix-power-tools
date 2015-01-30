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
#include "calvin_files/parsers/test/DataSetHeaderReaderTest.h"
//
#include "calvin_files/parsers/src/DataGroupHeaderReader.h"
#include "calvin_files/parsers/src/FileHeaderReader.h"
//
#include "util/Fs.h"
//

const std::string TEST_DATA_DAT_FILE = "../data/test.file.data_dat";

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetHeaderReaderTest );

void DataSetHeaderReaderTest::setUp()
{
  Fs::aptOpen(is, TEST_DATA_DAT_FILE.c_str(), std::ios::in | std::ios::binary);
}

void DataSetHeaderReaderTest::tearDown()
{
	is.close();
	fh.Clear();
}

void DataSetHeaderReaderTest::CreationTest()
{
	DataSetHeaderReader dphReader;
	CPPUNIT_ASSERT(1);
}

void DataSetHeaderReaderTest::ReadOneDataSetHeaderTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	DataGroupHeader dch;
	DataGroupHeaderReader dchReader;
	dchReader.ReadHeader(is, dch);

	// File position is at the start of the first DataSetHeader in the first DataGroup
	DataSetHeader dph;
	DataSetHeaderReader dphReader;
	u_int32_t nextDataSet = 0;
	CPPUNIT_ASSERT_NO_THROW(nextDataSet = dphReader.Read(is, dph));

	// Check the DataSetHeader information was read
	CPPUNIT_ASSERT(dph.GetName() == L"acquired data");
	CPPUNIT_ASSERT(dph.GetDataStartFilePos() != 0);
	CPPUNIT_ASSERT(dph.GetNameValParamCnt() == 2);
	CPPUNIT_ASSERT(dph.GetRowSize() == sizeof(u_int16_t));
	CPPUNIT_ASSERT(dph.GetDataSize() == sizeof(u_int16_t)*100);
	CPPUNIT_ASSERT(dph.GetRowCnt() == 100);
	CPPUNIT_ASSERT(dph.GetColumnCnt() == 1);
	CPPUNIT_ASSERT(dph.GetDataStartFilePos() == 0x3A0);

	// Check the data dataSet name value pairs
	ParameterNameValueTypeConstIt nvpBegin, nvpEnd;
	dph.GetNameValIterators(nvpBegin, nvpEnd);
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == L"Scanner");
	CPPUNIT_ASSERT(nvpBegin->GetValueText() == L"M10");
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == L"Pixel Size");
	CPPUNIT_ASSERT(nvpBegin->GetValueFloat() == 0.051f);
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin == nvpEnd);
}

void DataSetHeaderReaderTest::ReadMinimumInfoForOneDataSetHeaderTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	DataGroupHeader dch;
	DataGroupHeaderReader dchReader;
	dchReader.ReadHeader(is, dch);

	// File position is at the start of the first DataSetHeader in the first DataGroup
	DataSetHeader dph;
	DataSetHeaderReader dphReader;
	u_int32_t nextDataSet = 0;
	CPPUNIT_ASSERT_NO_THROW(nextDataSet = dphReader.ReadMinimumInfo(is, dph));

	// Check that the minimum DataSetHeader information was read.
	CPPUNIT_ASSERT(dph.GetName() == L"acquired data");
	CPPUNIT_ASSERT(dph.GetDataStartFilePos() != 0);
	CPPUNIT_ASSERT(dph.GetNameValParamCnt() == 0);
	CPPUNIT_ASSERT(dph.GetRowSize() == 0);
	CPPUNIT_ASSERT(dph.GetDataSize() == 0);
	CPPUNIT_ASSERT(dph.GetRowCnt() == 0);
	CPPUNIT_ASSERT(dph.GetColumnCnt() == 0);
	CPPUNIT_ASSERT(dph.GetDataStartFilePos() == 0x3A0);
}

void DataSetHeaderReaderTest::ReadAllTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	DataGroupHeader dch;
	DataGroupHeaderReader dchReader;
	dchReader.ReadHeader(is, dch);

	int32_t dpCount = 1; // fh.GetDataGroup(0).GetDataSetCnt(); is not ready at this time

	// File position is at the start of the first DataSetHeader in the first DataGroup
	DataSetHeaderReader dphReader;
	CPPUNIT_ASSERT_NO_THROW(dphReader.ReadAll(is, dch, dpCount));

	// Check that the right number of DataSets were read
	CPPUNIT_ASSERT(dch.GetDataSetCnt() == 1);

	// Check the DataSetHeader information was read
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
	CPPUNIT_ASSERT(dphBegin->GetColumnInfo(0).GetName() == L"Pixel");
	CPPUNIT_ASSERT(dphBegin->GetColumnInfo(0).GetColumnType() == UShortColType);
	CPPUNIT_ASSERT(dphBegin->GetColumnInfo(0).GetSize() == sizeof(u_int16_t));

	// Check that there are not any more DataSetHeaders
	++dphBegin;
	CPPUNIT_ASSERT(dphBegin == dphEnd);
}

void DataSetHeaderReaderTest::ReadAllMinimumInfoTest()
{
	FileHeaderReader fhReader(is, fh);
	fhReader.Read();

	DataGroupHeader dch;
	DataGroupHeaderReader dchReader;
	dchReader.ReadHeader(is, dch);

	int32_t dpCount = 1; // fh.GetDataGroup(0).GetDataSetCnt(); is not ready at this time

	// File position is at the start of the first DataSetHeader in the first DataGroup
	DataSetHeaderReader dphReader;
	CPPUNIT_ASSERT_NO_THROW(dphReader.ReadAllMinimumInfo(is, dch, dpCount));

	// Check that the right number of DataSets were read
	CPPUNIT_ASSERT(dch.GetDataSetCnt() == 1);

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
