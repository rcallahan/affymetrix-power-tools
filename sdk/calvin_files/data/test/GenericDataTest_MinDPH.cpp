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
//#include <math.h>
#include "calvin_files/data/test/GenericDataTest_MinDPH.h"
//
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( GenericDataTest_MinDPH );

#define TEST_DATA_DAT_FILE "../data/test.file.data_dat"

GenericDataTest_MinDPH::GenericDataTest_MinDPH()
{
	data = 0;
}

GenericDataTest_MinDPH::~GenericDataTest_MinDPH()
{

}

void GenericDataTest_MinDPH::setUp()
{
	// Create generic data header
	data = new GenericData;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);
	reader.ReadHeader(*data, GenericFileReader::ReadMinDataGroupHeader);
}

void GenericDataTest_MinDPH::tearDown()
{
	delete data;
}

void GenericDataTest_MinDPH::CreationTest()
{
	CPPUNIT_ASSERT(data);
}

void GenericDataTest_MinDPH::DataSetByIndexTest()
{
	// Check that a MinDP will be read-in completely before being accessed through GenericData.
	DataSet* dp = data->DataSet(0,0);
	CPPUNIT_ASSERT(dp);
	const DataSetHeader& dph = dp->Header();

	CPPUNIT_ASSERT(dph.GetName() == L"acquired data");
	CPPUNIT_ASSERT(dph.GetDataStartFilePos() != 0);
	CPPUNIT_ASSERT(dph.GetNameValParamCnt() == 2);
	CPPUNIT_ASSERT(dph.GetRowSize() == sizeof(u_int16_t));
	CPPUNIT_ASSERT(dph.GetDataSize() == sizeof(u_int16_t)*100);
	CPPUNIT_ASSERT(dph.GetRowCnt() == 100);
	CPPUNIT_ASSERT(dph.GetColumnCnt() == 1);
	CPPUNIT_ASSERT(dph.GetDataStartFilePos() == 0x3a0);

	// Check the data group name value pairs
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

	// Check the data set columns
	CPPUNIT_ASSERT(dph.GetColumnInfo(0).GetColumnType() == UShortColType);
	CPPUNIT_ASSERT(dph.GetColumnInfo(0).GetSize() == sizeof(u_int16_t));

	dp->Delete();
}

void GenericDataTest_MinDPH::DataSetByNameTest()
{
	// Check that a MinDP will be read-in completely before being accessed through GenericData.
	DataSet* dp = data->DataSet(L"First Data Cube", L"acquired data");
	CPPUNIT_ASSERT(dp);
	const DataSetHeader& dph = dp->Header();

	CPPUNIT_ASSERT(dph.GetName() == L"acquired data");
	CPPUNIT_ASSERT(dph.GetDataStartFilePos() != 0);
	CPPUNIT_ASSERT(dph.GetNameValParamCnt() == 2);
	CPPUNIT_ASSERT(dph.GetRowSize() == sizeof(u_int16_t));
	CPPUNIT_ASSERT(dph.GetDataSize() == sizeof(u_int16_t)*100);
	CPPUNIT_ASSERT(dph.GetRowCnt() == 100);
	CPPUNIT_ASSERT(dph.GetColumnCnt() == 1);
	CPPUNIT_ASSERT(dph.GetDataStartFilePos() == 0x3a0);

	// Check the data group name value pairs
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

	// Check the data set columns
	CPPUNIT_ASSERT(dph.GetColumnInfo(0).GetColumnType() == UShortColType);
	CPPUNIT_ASSERT(dph.GetColumnInfo(0).GetSize() == sizeof(u_int16_t));

	dp->Delete();
}

void GenericDataTest_MinDPH::FileIdentifierTest()
{
	AffymetrixGuidType id = data->FileIdentifier();
	CPPUNIT_ASSERT(id == "test-dat-guid");
}

void GenericDataTest_MinDPH::ArrayFileIdentifierTest()
{
	AffymetrixGuidType id = data->ArrayFileIdentifier();
	CPPUNIT_ASSERT(id == "test-array-guid");
}

void GenericDataTest_MinDPH::HeaderTest()
{
	FileHeader& fhdr = data->Header();
	CPPUNIT_ASSERT(fhdr.GetGenericDataHdr()->GetFileId() == "test-dat-guid");
}

void GenericDataTest_MinDPH::DataGroupNamesTest()
{
	std::vector<std::wstring> names;
	data->DataGroupNames(names);
	CPPUNIT_ASSERT(names.size() == 1);
	CPPUNIT_ASSERT(names[0] == L"First Data Cube");

}

void GenericDataTest_MinDPH::DataSetCntTest()
{
	CPPUNIT_ASSERT_THROW(data->DataSetCnt(L"none"), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT(data->DataSetCnt(L"First Data Cube") == 1);
	CPPUNIT_ASSERT_THROW(data->DataSetCnt(1), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT(data->DataSetCnt(0) == 1);

}

void GenericDataTest_MinDPH::DataSetNamesTest()
{
	// By name
	std::vector<std::wstring> names;
	CPPUNIT_ASSERT_THROW(data->DataSetNames(L"none", names), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_NO_THROW(data->DataSetNames(L"First Data Cube", names));
	CPPUNIT_ASSERT(names.size() == 1);
	CPPUNIT_ASSERT(names[0] == L"acquired data");

	// By index
	CPPUNIT_ASSERT_THROW(data->DataSetNames(1, names), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_NO_THROW(data->DataSetNames(0, names));
	CPPUNIT_ASSERT(names.size() == 1);
	CPPUNIT_ASSERT(names[0] == L"acquired data");
}
