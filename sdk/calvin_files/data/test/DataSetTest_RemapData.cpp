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
#include "calvin_files/data/test/DataSetTest_RemapData.h"
//
#include "calvin_files/data/src/DataException.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;

// NOTE:
// Requires the existence of a large DAT file in affy/sdk/calvin_files/data/data.
// This test is not included in the build because the file is not in CVS.
#define TEST_DATA_LARGE_DAT "../data/large_DAT_file"

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetTest_RemapData );

DataSetTest_RemapData::DataSetTest_RemapData()
{
	data = 0;
	dataPlane = 0;
}

DataSetTest_RemapData::~DataSetTest_RemapData()
{

}

void DataSetTest_RemapData::setUp()
{
	data = new GenericData;
	GenericFileReader reader;
	std::string name = TEST_DATA_LARGE_DAT;
	reader.SetFilename(name);
	try
	{
		reader.ReadHeader(*data);
		dataPlane = data->DataSet(0,0);
	}
	catch(affymetrix_calvin_exceptions::CalvinException&)
	{
		dataPlane = 0;
	}
}

void DataSetTest_RemapData::tearDown()
{
	if (dataPlane)
		dataPlane->Delete();
	delete data;
}

void DataSetTest_RemapData::CreationTest()
{
	// Check that the DataSet was created.
	CPPUNIT_ASSERT(dataPlane);
}

void DataSetTest_RemapData::RemapTest()
{
	CPPUNIT_ASSERT(dataPlane);
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	int32_t cols = dataPlane->Cols();

	Uint16Vector v;

	// Get data from the first view
	dataPlane->GetData(0, 0, 10, v);

	// Check the values.
	u_int16_t inten = 10;
	for (int32_t row = 0; row < 10; ++row, ++inten)
	{
		if (inten > 46000)
		{
			inten = 0;
		}
		CPPUNIT_ASSERT(inten == v[row]);
	}

	// Get data from another view - force a remap
	dataPlane->GetData(0, 150000000, 10, v);

	// Check the values.
	inten = 10;
	for (int32_t row = 0; row < 150000010; ++row, ++inten)
	{
		if (inten > 46000)
		{
			inten = 0;
		}

		if (row >= 150000000 && row < 150000010)
			CPPUNIT_ASSERT(inten == v[row-150000000]);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());

}

void DataSetTest_RemapData::RemapBackwardTest()
{
	CPPUNIT_ASSERT(dataPlane);
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	int32_t cols = dataPlane->Cols();

	Uint16Vector v;
	Uint16Vector vExpected1;
	Uint16Vector vExpected2;
	Uint16Vector vExpected3;
	Uint16Vector vExpected4;

	vExpected1.resize(10);
	vExpected2.resize(10);
	vExpected3.resize(10);
	vExpected4.resize(10);

	// Compute expected results
	u_int16_t inten = 10;
	for (int32_t row = 0; row < 150000010; ++row, ++inten)
	{
		if (inten > 46000)
		{
			inten = 0;
		}

		if ( row >= 0 && row < 10)
			vExpected1[row] = inten;

		if (row >= 149999980 && row < 149999990)
			vExpected2[row-149999980] = inten;

		if (row >= 149999990 && row < 150000000)
			vExpected3[row-149999990] = inten;

		if (row >= 150000000 && row < 150000010)
			vExpected4[row-150000000] = inten;
	}

	// Jump to 150,000,000 pixels into the file
	// The map should start 150,000,000 and extend 200 MB from there
	dataPlane->GetData(0, 150000000, 10, v);

	// Check the results

	for (int32_t row=0; row<10; ++row)
	{
		CPPUNIT_ASSERT(vExpected4[row] == v[row]);
	}

	// Jump forward in the file.  This should cause a remap.
	dataPlane->GetData(0, 149999990, 10, v);

	for (int32_t row=0; row<10; ++row)
	{
		if (vExpected3[row] == v[row])
			CPPUNIT_ASSERT(vExpected3[row] == v[row]);
	}

	// Jump forward in the file, again.  No remap.
	dataPlane->GetData(0, 149999980, 10, v);

	for (int32_t row=0; row<10; ++row)
	{
		CPPUNIT_ASSERT(vExpected2[row] == v[row]);
	}

	// Get data from the first view
	dataPlane->GetData(0, 0, 10, v);

	for (int32_t row=0; row<10; ++row)
	{
		CPPUNIT_ASSERT(vExpected1[row] == v[row]);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}

void DataSetTest_RemapData::RemapSizeTest()
{
	// A bug was found where a very large amount of data
	// would be mapped when data at the end of the file
	// was requested and a remap was required.

	CPPUNIT_ASSERT(dataPlane);
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();

	u_int16_t inten = 10;

	CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(rows-1, 0, inten));
}