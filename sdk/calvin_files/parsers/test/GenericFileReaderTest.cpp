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
#include "calvin_files/parsers/test/GenericFileReaderTest.h"
//
#include "calvin_files/array/src/ArrayId.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_parameter;

CPPUNIT_TEST_SUITE_REGISTRATION( GenericFileReaderTest );

#define TEST_DATA_FILE_HEADER_ONLY "../data/test.file.data_header_only"
#define TEST_NON_EXISTANT_FILE "../data/file_does_not_exist";
#define TEST_DATA_FILE_INVALID "../data/test.file.array_header_only"
#define TEST_DATA_INVALID_VERSION "../data/test.file.data_invalid_version"
#define TEST_DATA_FULL_FILE "../data/test.file.full_data_file"
#define TEST_DATA_DAT_FILE "../data/test.file.data_dat"


void GenericFileReaderTest::setUp()
{
}

void GenericFileReaderTest::tearDown()
{
}

void GenericFileReaderTest::testCreation()
{
	GenericFileReader reader;
	CPPUNIT_ASSERT(1);
}

void GenericFileReaderTest::testproperty_FileName()
{
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);
	CPPUNIT_ASSERT( reader.GetFilename() == name );
}

void GenericFileReaderTest::testmethod_ReadHeaderMinDP_when_file_exists()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data, GenericFileReader::ReadMinDataGroupHeader));
	CPPUNIT_ASSERT(data.Header().GetFilename() == TEST_DATA_DAT_FILE);

	Check_TEST_DATA_DAT_FILE_GenericHeader(data);

	// Check the DataGroupHeader information.
	CPPUNIT_ASSERT(data.Header().GetDataGroupCnt() == 1);
	CPPUNIT_ASSERT(data.Header().GetNumDataGroups() == 1);

	DataGroupHdrIt dchBegin, dchEnd;
	CPPUNIT_ASSERT_NO_THROW(data.Header().GetDataGroupIts(dchBegin, dchEnd));
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

void GenericFileReaderTest::testmethod_ReadHeaderMinDP_when_file_does_not_exist()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_NON_EXISTANT_FILE;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_THROW(reader.ReadHeader(data, GenericFileReader::ReadMinDataGroupHeader), FileNotFoundException);
}

void GenericFileReaderTest::testmethod_ReadHeaderMinDP_when_file_has_wrong_version()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_INVALID_VERSION;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_THROW(reader.ReadHeader(data, GenericFileReader::ReadMinDataGroupHeader), InvalidVersionException);
}

void GenericFileReaderTest::testmethod_ReadHeaderMinDP_when_file_is_not_valid()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_FILE_INVALID;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_THROW(reader.ReadHeader(data, GenericFileReader::ReadMinDataGroupHeader),InvalidFileTypeException);
}

void GenericFileReaderTest::testmethod_ReadHeader_when_file_has_wrong_version()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_INVALID_VERSION;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_THROW(reader.ReadHeader(data), InvalidVersionException);
}

void GenericFileReaderTest::testmethod_ReadHeader_when_file_is_not_valid()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_FILE_INVALID;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_THROW(reader.ReadHeader(data), InvalidFileTypeException);
}

void GenericFileReaderTest::testmethod_ReadHeader_when_file_does_not_exist()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_NON_EXISTANT_FILE;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_THROW(reader.ReadHeader(data), FileNotFoundException);
}

void GenericFileReaderTest::testmethod_ReadHeader_when_file_exists()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data));
	CPPUNIT_ASSERT(data.Header().GetFilename() == TEST_DATA_DAT_FILE);

	Check_TEST_DATA_DAT_FILE_GenericHeader(data);

	// Check the DataGroupHeader information.
	CPPUNIT_ASSERT(data.Header().GetDataGroupCnt() == 1);

	DataGroupHdrIt dchBegin, dchEnd;
	CPPUNIT_ASSERT_NO_THROW(data.Header().GetDataGroupIts(dchBegin, dchEnd));
	CPPUNIT_ASSERT(dchBegin->GetName() == L"First Data Cube");
	CPPUNIT_ASSERT(dchBegin->GetDataSetCnt() == 1);
	CPPUNIT_ASSERT(dchBegin->GetDataSetPos() == 0x2d1);
	CPPUNIT_ASSERT(dchBegin->GetNextGroupPos() == 0x468);

	// Check that the minimum DataSetHeader information was read.
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

void GenericFileReaderTest::Check_TEST_DATA_DAT_FILE_GenericHeader(GenericData& data)
{
	// Check the file generic header
	GenericDataHeader* gdh = data.Header().GetGenericDataHdr();
	CPPUNIT_ASSERT(gdh->GetFileTypeId() == SCAN_ACQUISITION_DATA_TYPE);
	CPPUNIT_ASSERT(gdh->GetFileId() == "test-dat-guid");
	CPPUNIT_ASSERT(data.FileIdentifier() == "test-dat-guid");
	CPPUNIT_ASSERT(gdh->GetFileCreationTime() == L"2004-07-04T11:12:13Z");
	CPPUNIT_ASSERT(gdh->GetLocale() == L"en-US");

	ParameterNameValueTypeIt nvpBegin, nvpEnd;
	gdh->GetNameValIterators(nvpBegin, nvpEnd);
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == ARRAY_TYPE_PARAM_NAME);
	CPPUNIT_ASSERT(nvpBegin->GetValueText() == L"Hg-U133A");
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == ARRAY_BARCODE_PARAM_NAME);
	CPPUNIT_ASSERT(nvpBegin->GetValueText() == L"Barcode");
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == L"Parameter1");
	CPPUNIT_ASSERT(nvpBegin->GetValueText() == L"Value1");
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin == nvpEnd);

	// Check the parent header 
	CPPUNIT_ASSERT(gdh->GetParentCnt() == 1);
	GenDataHdrVectorIt parentBegin, parentEnd;
	gdh->GetParentIterators(parentBegin, parentEnd);
	CPPUNIT_ASSERT(parentBegin != parentEnd);
	CPPUNIT_ASSERT(parentBegin->GetFileTypeId() == ARRAY_TYPE_IDENTIFIER);
	CPPUNIT_ASSERT(parentBegin->GetFileId() == "test-array-guid");
	CPPUNIT_ASSERT(data.ArrayFileIdentifier() == "test-array-guid");
	CPPUNIT_ASSERT(parentBegin->GetFileCreationTime() == L"2004-07-01T13:14:15Z");
	CPPUNIT_ASSERT(parentBegin->GetLocale() == L"en-US");
	
	parentBegin->GetNameValIterators(nvpBegin, nvpEnd);
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == ARRAY_TYPE_PARAM_NAME);
	CPPUNIT_ASSERT(nvpBegin->GetValueText() == L"Hg-U133A");
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin != nvpEnd);
	CPPUNIT_ASSERT(nvpBegin->GetName() == ARRAY_LOT_PARAM_NAME);
	CPPUNIT_ASSERT(nvpBegin->GetValueText() == L"Thanks alot");
	++nvpBegin;
	CPPUNIT_ASSERT(nvpBegin == nvpEnd);

	++parentBegin;
	CPPUNIT_ASSERT(parentBegin == parentEnd);
}

void GenericFileReaderTest::OpenTest()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);

	CPPUNIT_ASSERT_NO_THROW(reader.Open(data));
	reader.Close();
}

void GenericFileReaderTest::GetDataGroupCntTest()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);

	reader.Open(data);
	CPPUNIT_ASSERT(reader.GetDataGroupCnt() == 1);
	reader.Close();

}

void GenericFileReaderTest::GetDataGroupReaderByIndexTest()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);

	reader.Open(data);
	DataGroupReader dcReader = reader.GetDataGroupReader(0);
	DataSetReader dpReader = dcReader.GetDataSetReader(0);

	int32_t rows = data.Header().GetDataGroup(0).GetDataSet(0).GetRowCnt();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dpReader.Read(value));
		CPPUNIT_ASSERT(value == expected);
	}
	reader.Close();
}

void GenericFileReaderTest::GetDataGroupReaderByNameTest()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);

	reader.Open(data);
	DataGroupReader dcReader = reader.GetDataGroupReader(L"First Data Cube");
	DataSetReader dpReader = dcReader.GetDataSetReader(0);

	int32_t rows = data.Header().GetDataGroup(0).GetDataSet(0).GetRowCnt();
	for( int32_t row = 0; row < rows; ++row )
	{
		u_int16_t expected = (u_int16_t)(row*10+row);
		u_int16_t value;
		CPPUNIT_ASSERT_NO_THROW(dpReader.Read(value));
		CPPUNIT_ASSERT(value == expected);
	}
	reader.Close();
}

void GenericFileReaderTest::CloseNoOpenTest()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);

	CPPUNIT_ASSERT_NO_THROW(reader.Close());
}

void GenericFileReaderTest::ReadHeaderNoDataGroupHeaderTest()
{
	GenericData data;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE;
	reader.SetFilename(name);
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data, GenericFileReader::ReadNoDataGroupHeader));
	CPPUNIT_ASSERT(data.Header().GetFilename() == TEST_DATA_DAT_FILE);

	Check_TEST_DATA_DAT_FILE_GenericHeader(data);

	// Check the DataGroupHeader information.
	CPPUNIT_ASSERT(data.Header().GetDataGroupCnt() == 0);
	CPPUNIT_ASSERT(data.Header().GetNumDataGroups() == 1);


}

void GenericFileReaderTest::ReadHeaderOfAFileWithMultipleDataGroups()
{
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename("../data/small_CDF_file");
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data, GenericFileReader::ReadAllHeaders));

	// Get the data group names
	std::vector<std::wstring> names;
	data.DataGroupNames(names);
	u_int32_t count = data.DataGroupCnt();

	wchar_t name[100];

	// Check the names
	for (u_int32_t i = 0; i < count; ++i)
	{
		std::wstring expectedName;
		if (i == 0)
			expectedName = L"Probe Set Names";
		else
		{
			FormatString1(name, 100, L"biob_%d", i-1);
			expectedName = name;
		}

		CPPUNIT_ASSERT(names[i] == expectedName);
	}
}

void GenericFileReaderTest::ReadHeaderOfFileWithReserveStringParameters()
{
	GenericData data;
	GenericFileReader reader;
	reader.SetFilename("../data/small_dat_file_with_reserved_string_parameters");
	CPPUNIT_ASSERT_NO_THROW(reader.ReadHeader(data, GenericFileReader::ReadAllHeaders));

	// Check array id because it has reserved space and is stored as an ASCII string.
	CPPUNIT_ASSERT(data.ArrayIdentifier() == "smellsliketeenspirit");

	// Get the raw parameter to check the length.
	bool found = false;
	int32_t parentHeaderCount = data.Header().GetGenericDataHdr()->GetParentCnt();
	for (int32_t index = 0; index < parentHeaderCount; ++index)
	{
		GenericDataHeader genericHdr = data.Header().GetGenericDataHdr()->GetParent(index);
		if (genericHdr.GetFileTypeId() == ARRAY_TYPE_IDENTIFIER)
		{
			ParameterNameValueType nvt;
			CPPUNIT_ASSERT(genericHdr.FindNameValParam(ARRAY_ID_PARAM_NAME, nvt));
			CPPUNIT_ASSERT(nvt.GetMIMEValue().Size() == AFFY_GUID_LEN);
			found = true;
			break;
		}
	}

	CPPUNIT_ASSERT(found);

	// Check reserved wchar string.
	ParameterNameValueType result;
	CPPUNIT_ASSERT(data.Header().GetGenericDataHdr()->FindNameValParam(L"fixedlen", result));
	CPPUNIT_ASSERT(result.GetValueText() == L"twenty-five");
	CPPUNIT_ASSERT(result.GetMIMEValue().Size() == 25*sizeof(u_int16_t));
}
