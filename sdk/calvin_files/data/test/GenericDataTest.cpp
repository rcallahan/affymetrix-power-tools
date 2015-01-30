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
#include "calvin_files/data/test/GenericDataTest.h"
//
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
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

CPPUNIT_TEST_SUITE_REGISTRATION( GenericDataTest );


const std::string SMALL_CEL_FILE = "../data/small_cel_file";


GenericDataTest::GenericDataTest()
{
	data = 0;
}

GenericDataTest::~GenericDataTest()
{

}

void GenericDataTest::setUp()
{
	data = new GenericData;
	GenericFileReader reader;
	reader.SetFilename(SMALL_CEL_FILE);
	reader.Open(*data);
}

void GenericDataTest::tearDown()
{
	data->Close();
	delete data;
}

void GenericDataTest::CreationTest()
{
	CPPUNIT_ASSERT(data);
}


void GenericDataTest::FileIdentifierTest()
{
	AffymetrixGuidType id = data->FileIdentifier();
	CPPUNIT_ASSERT(id.length() > 0);
}

void GenericDataTest::ArrayFileIdentifierTest()
{
	AffymetrixGuidType id = data->ArrayFileIdentifier();
	CPPUNIT_ASSERT(id.length() == 0);
}

void GenericDataTest::HeaderTest()
{
	FileHeader& fhdr = data->Header();
	CPPUNIT_ASSERT(fhdr.GetFilename() == SMALL_CEL_FILE);
}

void GenericDataTest::DataGroupCntTest()
{
	CPPUNIT_ASSERT(data->DataGroupCnt() == 1);
}

void GenericDataTest::DataGroupNamesTest()
{
	WStringVector names;
	CPPUNIT_ASSERT_NO_THROW(data->DataGroupNames(names));
	CPPUNIT_ASSERT(names.size() == 1);
	CPPUNIT_ASSERT(names[0] == L"Default Group");
}

void GenericDataTest::DataSetCntByIndexTest()
{
	CPPUNIT_ASSERT(data->DataSetCnt(0) == 5);
}

void GenericDataTest::DataSetCntByIndexErrorTest()
{
	CPPUNIT_ASSERT_THROW(data->DataSetCnt(1), affymetrix_calvin_exceptions::DataGroupNotFoundException);
}

void GenericDataTest::DataSetCntByNameTest()
{
	CPPUNIT_ASSERT(data->DataSetCnt(L"Default Group") == 5);
}

void GenericDataTest::DataSetCntByNameErrorTest()
{
	CPPUNIT_ASSERT_THROW(data->DataSetCnt(L"Nonsense"), affymetrix_calvin_exceptions::DataGroupNotFoundException);
}

void GenericDataTest::DataSetNamesByGroupIndexTest()
{
	WStringVector names;
	CPPUNIT_ASSERT_NO_THROW(data->DataSetNames(0, names));
	CPPUNIT_ASSERT(names.size() == 5);
	CPPUNIT_ASSERT(names[0] == L"Intensity");
	CPPUNIT_ASSERT(names[1] == L"StdDev");
	CPPUNIT_ASSERT(names[2] == L"Pixel");
	CPPUNIT_ASSERT(names[3] == L"Outlier");
	CPPUNIT_ASSERT(names[4] == L"Mask");
}

void GenericDataTest::DataSetNamesByGroupIndexErrorTest()
{
	WStringVector names;
	CPPUNIT_ASSERT_THROW(data->DataSetNames(1, names), affymetrix_calvin_exceptions::DataGroupNotFoundException);
}

void GenericDataTest::DataSetNamesByGroupNameTest()
{
	WStringVector names;
	CPPUNIT_ASSERT_NO_THROW(data->DataSetNames(L"Default Group", names));
	CPPUNIT_ASSERT(names.size() == 5);
	CPPUNIT_ASSERT(names[0] == L"Intensity");
	CPPUNIT_ASSERT(names[1] == L"StdDev");
	CPPUNIT_ASSERT(names[2] == L"Pixel");
	CPPUNIT_ASSERT(names[3] == L"Outlier");
	CPPUNIT_ASSERT(names[4] == L"Mask");
}

void GenericDataTest::DataSetNamesByGroupNameErrorTest()
{
	WStringVector names;
	CPPUNIT_ASSERT_THROW(data->DataSetNames(L"Nonsense", names), affymetrix_calvin_exceptions::DataGroupNotFoundException);
}

void GenericDataTest::DataSetByIndexUsingMMTest()
{
	data->UseMemoryMapping(true);	// Default
	DataSet* ds = data->DataSet(0, 0);
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"Intensity");
	ds->Delete();
}

void GenericDataTest::DataSetByNameUsingMMTest()
{
	data->UseMemoryMapping(true);	// Default
	DataSet* ds = data->DataSet(L"Default Group", L"Intensity");
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"Intensity");
	ds->Delete();
}

void GenericDataTest::DataSetByIndexUsingMMErrorTest()
{
	data->UseMemoryMapping(true);	// Default
	CPPUNIT_ASSERT_THROW(data->DataSet(1, 0), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_THROW(data->DataSet(0, 5), affymetrix_calvin_exceptions::DataSetNotFoundException);
}

void GenericDataTest::DataSetByNameUsingMMErrorTest()
{
	data->UseMemoryMapping(true);	// Default
	CPPUNIT_ASSERT_THROW(data->DataSet(L"Nonsense Group", L"Intensity"), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_THROW(data->DataSet(L"Default Group", L"Nonsense"), affymetrix_calvin_exceptions::DataSetNotFoundException);
}

void GenericDataTest::DataSetByIndexUsingFStreamTest()
{
	data->UseMemoryMapping(false);
	DataSet* ds = data->DataSet(0, 0);
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"Intensity");
	ds->Delete();
}

void GenericDataTest::DataSetByNameUsingFStreamTest()
{
	data->UseMemoryMapping(false);
	DataSet* ds = data->DataSet(L"Default Group", L"Intensity");
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"Intensity");
	ds->Delete();
}

void GenericDataTest::DataSetByIndexUsingFStreamWithLoadEntireDataSetTest()
{
	data->UseMemoryMapping(false);
	data->LoadEntireDataSetHint(true);
	DataSet* ds = data->DataSet(0, 0);
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"Intensity");
	ds->Delete();
}

void GenericDataTest::DataSetByNameUsingFStreamWithLoadEntireDataSetTest()
{
	data->UseMemoryMapping(false);
	data->LoadEntireDataSetHint(true);
	DataSet* ds = data->DataSet(L"Default Group", L"Intensity");
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"Intensity");
	ds->Delete();
}

void GenericDataTest::DataSetByIndexUsingFStreamErrorTest()
{
	data->UseMemoryMapping(false);
	CPPUNIT_ASSERT_THROW(data->DataSet(1, 0), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_THROW(data->DataSet(0, 5), affymetrix_calvin_exceptions::DataSetNotFoundException);
}

void GenericDataTest::DataSetByNameUsingFStreamErrorTest()
{
	data->UseMemoryMapping(false);
	CPPUNIT_ASSERT_THROW(data->DataSet(L"Nonsense Group", L"Intensity"), affymetrix_calvin_exceptions::DataGroupNotFoundException);
	CPPUNIT_ASSERT_THROW(data->DataSet(L"Default Group", L"Nonsense"), affymetrix_calvin_exceptions::DataSetNotFoundException);
}

void GenericDataTest::DataSetByIndexUsingMMAndFStreamTest()
{
	data->UseMemoryMapping(false);
	DataSet* ds = data->DataSet(0, 0);
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"Intensity");
	ds->Delete();
	data->UseMemoryMapping(true);
	ds = data->DataSet(0, 1);
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"StdDev");
	ds->Delete();
}

void GenericDataTest::DataSetByNameUsingMMAndFStreamTest()
{
	data->UseMemoryMapping(false);
	DataSet* ds = data->DataSet(L"Default Group", L"Intensity");
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"Intensity");
	ds->Delete();
	data->UseMemoryMapping(true);
	ds = data->DataSet(L"Default Group", L"StdDev");
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Header().GetName() == L"StdDev");
	ds->Delete();
}

void GenericDataTest::DataGroupUsingMMTest()
{
	data->UseMemoryMapping(true);
	DataGroup dg = data->DataGroup(0x495);
	CPPUNIT_ASSERT(dg.Header().GetName() == L"Default Group");
}

void GenericDataTest::DataGroupUsingFStreamTest()
{
	data->UseMemoryMapping(false);
	DataGroup dg = data->DataGroup(0x495);
	CPPUNIT_ASSERT(dg.Header().GetName() == L"Default Group");
}

void GenericDataTest::DataGroupUsingFStreamLoadEntireTest()
{
	data->UseMemoryMapping(false);
	data->LoadEntireDataSetHint(true);
	DataGroup dg = data->DataGroup(0x495);
	CPPUNIT_ASSERT(dg.Header().GetName() == L"Default Group");
}

void GenericDataTest::DataGroupUsingMMAndFStreamTest()
{
	data->UseMemoryMapping(true);
	DataGroup dg1 = data->DataGroup(0x495);
	CPPUNIT_ASSERT(dg1.Header().GetName() == L"Default Group");
	data->UseMemoryMapping(false);
	DataGroup dg2 = data->DataGroup(0x495);
	CPPUNIT_ASSERT(dg2.Header().GetName() == L"Default Group");
}

void GenericDataTest::FindDataGroupHeaderByNameTest()
{
	const DataGroupHeader* dgh = data->FindDataGroupHeader(L"none");
	CPPUNIT_ASSERT(0 == dgh);
	dgh = data->FindDataGroupHeader(L"default Group");
	CPPUNIT_ASSERT(0 == dgh);
	dgh = data->FindDataGroupHeader(L"Default Group");
	CPPUNIT_ASSERT(0 != dgh);
	CPPUNIT_ASSERT(dgh->GetName() == L"Default Group");
}

void GenericDataTest::FindDataGroupHeaderByIndexTest()
{
	const DataGroupHeader* dch = data->FindDataGroupHeader(1);
	CPPUNIT_ASSERT(0 == dch);
	dch = data->FindDataGroupHeader(0);
	CPPUNIT_ASSERT(0 != dch);
	CPPUNIT_ASSERT(dch->GetName() == L"Default Group");
}

void GenericDataTest::FindDataSetHeaderByNameTest()
{
	DataGroupHeader* dgh = data->FindDataGroupHeader(L"Default Group");
	CPPUNIT_ASSERT(dgh);
	const DataSetHeader* dsh = data->FindDataSetHeader(0, L"none");
	CPPUNIT_ASSERT(0 == dsh);
	dsh = data->FindDataSetHeader(dgh, L"intensity");
	CPPUNIT_ASSERT(0 == dsh);
	dsh = data->FindDataSetHeader(dgh, L"Intensity");
	CPPUNIT_ASSERT(0 != dsh);
	CPPUNIT_ASSERT(dsh->GetName() == L"Intensity");
	dsh = data->FindDataSetHeader(dgh, L"StdDev");
	CPPUNIT_ASSERT(0 != dsh);
	CPPUNIT_ASSERT(dsh->GetName() == L"StdDev");
}

void GenericDataTest::FindDataSetHeaderByIndexTest()
{
	DataGroupHeader* dgh = data->FindDataGroupHeader(L"Default Group");
	CPPUNIT_ASSERT(dgh);
	const DataSetHeader* dsh = data->FindDataSetHeader(0, 0);
	CPPUNIT_ASSERT(0 == dsh);
	dsh = data->FindDataSetHeader(dgh, 5);
	CPPUNIT_ASSERT(0 == dsh);
	dsh = data->FindDataSetHeader(dgh, 0);
	CPPUNIT_ASSERT(0 != dsh);
	CPPUNIT_ASSERT(dsh->GetName() == L"Intensity");
	dsh= data->FindDataSetHeader(dgh, 1);
	CPPUNIT_ASSERT(0 != dsh);
	CPPUNIT_ASSERT(dsh->GetName() == L"StdDev");
}

void GenericDataTest::ClearTest()
{
	int32_t count = data->DataGroupCnt();
	CPPUNIT_ASSERT(count == 1);
	data->Clear();
	count = data->DataGroupCnt();
	CPPUNIT_ASSERT(count == 0);
}

