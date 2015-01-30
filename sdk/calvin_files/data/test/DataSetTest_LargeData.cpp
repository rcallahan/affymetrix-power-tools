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
#include "calvin_files/data/test/DataSetTest_LargeData.h"
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
// Requires the exsitance of a large DAT file in affy/sdk/calvin_files/data/data.
// This test is not included in the build because the file is not in CVS.
#define TEST_DATA_LARGE_DAT "../data/large_DAT_file"

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetTest_LargeData );

DataSetTest_LargeData::DataSetTest_LargeData()
{
	data = 0;
	dataPlane = 0;
}

DataSetTest_LargeData::~DataSetTest_LargeData()
{

}

void DataSetTest_LargeData::setUp()
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

void DataSetTest_LargeData::tearDown()
{
	if (dataPlane)
		dataPlane->Delete();
	delete data;
}

void DataSetTest_LargeData::CreationTest()
{
	// Check that the DataSet was created.
	CPPUNIT_ASSERT(dataPlane);
}

// A bug was discovered where requesting
// more than 200MB in GetDataRaw on a
// DataSet with one column (like a DAT)
// caused a crash.
void DataSetTest_LargeData::GetDataRaw300MBTest()
{
	CPPUNIT_ASSERT(dataPlane);
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	int32_t cols = dataPlane->Cols();

	int32_t threeHundredMB = 150000000; // 2 bytes per pixel = 300 MB
	u_int16_t* values = new u_int16_t[threeHundredMB];

	CPPUNIT_ASSERT_NO_THROW(dataPlane->GetDataRaw(0, 0, threeHundredMB, values));

	// Check the values
	u_int16_t inten = 10;
	for (int32_t row = 0; row < 150000000; ++row, ++inten)
	{
		if (inten > 46000)
		{
			inten = 0;
		}

		// Spot check
		if (row % 99 == 0)
			CPPUNIT_ASSERT(inten == values[row]);			

		// Check at the 200MB boundary.
		if (row >= 99999990 && row < 100000010)
			CPPUNIT_ASSERT(inten == values[row]);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());

}

// A similar problem to the GetDataRaw
// bug above exists in GetData that returns
// a vector.
void DataSetTest_LargeData::GetData300MBTest()
{
	CPPUNIT_ASSERT(dataPlane);
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	int32_t cols = dataPlane->Cols();

	int32_t threeHundredMB = 150000000; // 2 bytes per pixel = 300MB
	Uint16Vector v;

	CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(0, 0, threeHundredMB, v));

		// Check the values
	u_int16_t inten = 10;
	for (int32_t row = 0; row < 150000000; ++row, ++inten)
	{
		if (inten > 46000)
		{
			inten = 0;
		}

		// Spot check
		if (row % 99 == 0)
			CPPUNIT_ASSERT(inten == v[row]);			

		// Check at the 200MB boundary.
		if (row >= 99999990 && row < 100000010)
			CPPUNIT_ASSERT(inten == v[row]);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());

}