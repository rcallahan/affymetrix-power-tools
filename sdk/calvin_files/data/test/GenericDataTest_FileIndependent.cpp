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
#include "calvin_files/data/test/GenericDataTest_FileIndependent.h"
//
#include "calvin_files/array/src/ArrayId.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/utils/src/DateTime.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;

CPPUNIT_TEST_SUITE_REGISTRATION( GenericDataTest_FileIndependent );

GenericDataTest_FileIndependent::GenericDataTest_FileIndependent()
{
	data = 0;
	header = 0;
	parent = 0;
	dch = 0;
	dphPI = 0;
	dphGrid = 0;
}

GenericDataTest_FileIndependent::~GenericDataTest_FileIndependent()
{

}

void GenericDataTest_FileIndependent::setUp()
{
	// Create generic data header
	data = new GenericData;

	data->Header().SetFilename("../data/test.file.data_dat");

	header = new GenericDataHeader;
	header->SetFileCreationTime(L"20040823T17:06:00Z");
//	header->SetFileCreationTime(DateTime::GetCurrentDateTime().ToString().c_str()); // change time to wstring?
	header->SetFileTypeId(INTENSITY_DATA_TYPE);
//	header->SetFileId(AffymetrixGuid::GenerateNewGuid());
	header->SetFileId("someuniquedatfileid");
	header->SetLocale(L"en-US");

	// Create parent array file header
	parent = new GenericDataHeader;
	parent->SetFileTypeId(ARRAY_TYPE_IDENTIFIER);
	parent->SetLocale(L"en-US");
//	parent->SetFileId(AffymetrixGuid::GenerateNewGuid());
	parent->SetFileId("someuniquearrayfileid");
	parent->SetFileCreationTime(L"20031225T18:23:00Z");

	ParameterNameValueType nvt;
	nvt.SetName(ARRAY_ID_PARAM_NAME);
	nvt.SetValueAscii("arrayidis17");
	parent->AddNameValParam(nvt);

	header->AddParent(*parent);

	// Add GenericDataHeader to the FileHeader
	data->Header().SetGenericDataHdr(*header);

	// Create DataGroupHeaders
	dch = new DataGroupHeader;
	dch->SetName(L"Default");

	// Create DataSetHeaders
	dphPI = new DataSetHeader;
	dphPI->SetName(L"pixel intensity");
	ParameterNameValueType param;
	param.SetName(L"Scanner");
	param.SetValueText(L"M10");
	dphPI->AddNameValParam(param);
	dphPI->AddUShortColumn(L"Intensity");
	dphPI->SetRowCnt(1);

	dphGrid = new DataSetHeader;
	dphGrid->SetName(L"grid coordinates");
	param.SetName(L"Corner Pattern");
	param.SetValueText(L"Checkerboard");
	dphGrid->AddNameValParam(param);
	dphGrid->AddUShortColumn(L"GridULX");
	dphGrid->AddUShortColumn(L"GridULY");
	dphGrid->AddUShortColumn(L"GridURX");
	dphGrid->AddUShortColumn(L"GridURY");

	dphGrid->AddUShortColumn(L"GridLRX");
	dphGrid->AddUShortColumn(L"GridLRY");
	dphGrid->AddUShortColumn(L"GridLLX");
	dphGrid->AddUShortColumn(L"GridLLY");
	dphGrid->SetRowCnt(1);

	// Add DataGroupHeaders
	dch->AddDataSetHdr(*dphPI);
	dch->AddDataSetHdr(*dphGrid);

	data->Header().AddDataGroupHdr(*dch);
}

void GenericDataTest_FileIndependent::tearDown()
{
	delete data;
	delete header;
	delete parent;
	delete dch;
	delete dphPI;
	delete dphGrid;
}

void GenericDataTest_FileIndependent::testCreation()
{
	GenericData Generic;
	CPPUNIT_ASSERT(1);
}

void GenericDataTest_FileIndependent::testproperty_Header()
{
	FileHeader& fhdr = data->Header();
	CPPUNIT_ASSERT(fhdr.GetGenericDataHdr()->GetFileId() == header->GetFileId());
}

void GenericDataTest_FileIndependent::testmethod_DataSet_ByIndex()
{
	DataSet* dp = data->DataSet(0,0);
	CPPUNIT_ASSERT(0 != dp);
	dp->Delete();
	dp = data->DataSet(0,1);
	CPPUNIT_ASSERT(0 != dp);
	dp->Delete();
	CPPUNIT_ASSERT_THROW(data->DataSet(1,0), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_THROW(data->DataSet(0,2), affymetrix_calvin_exceptions::DataSetNotFoundException);
	// TODO: Add more asserts to make sure that I am getting back the right data set.
}

void GenericDataTest_FileIndependent::testmethod_DataSet_ByName()
{
	GenericData d;
	CPPUNIT_ASSERT_THROW(data->DataSet(L"none", dphPI->GetName()), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_THROW(data->DataSet(L"Default", L"none"), affymetrix_calvin_exceptions::DataSetNotFoundException);
	DataSet* dp = data->DataSet(dch->GetName(), dphPI->GetName());
	CPPUNIT_ASSERT(0 != dp);
	dp->Delete();
	dp = data->DataSet(dch->GetName(), dphGrid->GetName());
	CPPUNIT_ASSERT(0 != dp);
	dp->Delete();
	// TODO: Add more asserts to make sure that I am getting back the right data set.
}

void GenericDataTest_FileIndependent::testmethod_FindDataSetHeader_ByName()
{
	const DataSetHeader* dph = data->FindDataSetHeader(0, dphPI->GetName());
	CPPUNIT_ASSERT(0 == dph);
	dph = data->FindDataSetHeader(dch, L"none");
	CPPUNIT_ASSERT(0 == dph);
	dph = data->FindDataSetHeader(dch, dphPI->GetName());
	CPPUNIT_ASSERT(0 != dph);
	CPPUNIT_ASSERT(dph->GetName() == dphPI->GetName());
	dph = data->FindDataSetHeader(dch, dphGrid->GetName());
	CPPUNIT_ASSERT(0 != dph);
	CPPUNIT_ASSERT(dph->GetName() == dphGrid->GetName());
}

void GenericDataTest_FileIndependent::testmethod_FindDataSetHeader_ByIndex()
{
	const DataSetHeader* dph = data->FindDataSetHeader(0, 0);
	CPPUNIT_ASSERT(0 == dph);
	dph = data->FindDataSetHeader(dch, 2);
	CPPUNIT_ASSERT(0 == dph);
	dph = data->FindDataSetHeader(dch, 0);
	CPPUNIT_ASSERT(0 != dph);
	CPPUNIT_ASSERT(dph->GetName() == dphPI->GetName());
	dph = data->FindDataSetHeader(dch, 1);
	CPPUNIT_ASSERT(0 != dph);
	CPPUNIT_ASSERT(dph->GetName() == dphGrid->GetName());
}

void GenericDataTest_FileIndependent::testmethod_DataGroupNames()
{
	std::vector<std::wstring> names;
	data->DataGroupNames(names);
	CPPUNIT_ASSERT(names.size() == 1);
	CPPUNIT_ASSERT(names[0] == dch->GetName());
}

void GenericDataTest_FileIndependent::testproperty_DataGroupCnt()
{
	int32_t count = data->DataGroupCnt();
	CPPUNIT_ASSERT(count == 1);
}

void GenericDataTest_FileIndependent::testproperty_ArrayFileIdentifier()
{
	AffymetrixGuidType id = data->ArrayFileIdentifier();
	CPPUNIT_ASSERT(id == parent->GetFileId());
}

void GenericDataTest_FileIndependent::testproperty_FileIdentifier()
{
	AffymetrixGuidType id = data->FileIdentifier();
	CPPUNIT_ASSERT(id == header->GetFileId());
}

void GenericDataTest_FileIndependent::testproperty_ArrayIdentifier()
{
	AffymetrixGuidType id = data->ArrayIdentifier();
	CPPUNIT_ASSERT(id == "arrayidis17");
}

void GenericDataTest_FileIndependent::testmethod_Clear()
{
	int32_t count = data->DataGroupCnt();
	CPPUNIT_ASSERT(count == 1);
	data->Clear();
	count = data->DataGroupCnt();
	CPPUNIT_ASSERT(count == 0);
}

void GenericDataTest_FileIndependent::testmethod_DataSetCnt()
{
	CPPUNIT_ASSERT_THROW(data->DataSetCnt(L"none"), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT(data->DataSetCnt(L"Default") == 2);
	CPPUNIT_ASSERT_THROW(data->DataSetCnt(1), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT(data->DataSetCnt(0) == 2);
}

void GenericDataTest_FileIndependent::testmethod_DataSetNames()
{
	// By name
	std::vector<std::wstring> names;
	CPPUNIT_ASSERT_THROW(data->DataSetNames(L"none", names), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_NO_THROW(data->DataSetNames(L"Default", names));
	CPPUNIT_ASSERT(names.size() == 2);
	CPPUNIT_ASSERT(names[0] == dphPI->GetName());
	CPPUNIT_ASSERT(names[1] == dphGrid->GetName());

	// By index
	CPPUNIT_ASSERT_THROW(data->DataSetNames(1, names), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_NO_THROW(data->DataSetNames(0, names));
	CPPUNIT_ASSERT(names.size() == 2);
	CPPUNIT_ASSERT(names[0] == dphPI->GetName());
	CPPUNIT_ASSERT(names[1] == dphGrid->GetName());
}

void GenericDataTest_FileIndependent::testmethod_FindDataGroupHeader_ByName()
{
	const DataGroupHeader* dch = data->FindDataGroupHeader(L"none");
	CPPUNIT_ASSERT(0 == dch);
	dch = data->FindDataGroupHeader(L"default");
	CPPUNIT_ASSERT(0 == dch);
	dch = data->FindDataGroupHeader(L"Default");
	CPPUNIT_ASSERT(0 != dch);
	CPPUNIT_ASSERT(dch->GetName() == L"Default");
}

void GenericDataTest_FileIndependent::testmethod_FindDataGroupHeader_ByIndex()
{
	const DataGroupHeader* dch = data->FindDataGroupHeader(1);
	CPPUNIT_ASSERT(0 == dch);
	dch = data->FindDataGroupHeader(0);
	CPPUNIT_ASSERT(0 != dch);
	CPPUNIT_ASSERT(dch->GetName() == L"Default");
}
