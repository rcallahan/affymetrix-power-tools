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
#include "calvin_files/data/test/DataSetTest_DATGridSet.h"
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

#define TEST_DATA_DAT_FILE_WITH_GRIDS "../data/test.file.data_dat_with_grid"

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetTest_DATGridSet );

DataSetTest_DATGridSet::DataSetTest_DATGridSet()
{
	data = 0;
	dataPlane = 0;
}

DataSetTest_DATGridSet::~DataSetTest_DATGridSet()
{

}

void DataSetTest_DATGridSet::setUp()
{
	data = new GenericData;
	GenericFileReader reader;
	std::string name = TEST_DATA_DAT_FILE_WITH_GRIDS;
	reader.SetFilename(name);
	reader.ReadHeader(*data);
	try
	{
		dataPlane = data->DataSet(0,1);	// grids is the second DataSet of the first DataGroup.
	}
	catch(affymetrix_calvin_exceptions::CalvinException&)
	{
		dataPlane = 0;
	}
}

void DataSetTest_DATGridSet::tearDown()
{
	if (dataPlane)
		dataPlane->Delete();
	delete data;
}

void DataSetTest_DATGridSet::testCreation()
{
	// Check that the DataSet was created.
	CPPUNIT_ASSERT(dataPlane);
}

void DataSetTest_DATGridSet::testSingleElementAccess()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	int32_t cols = dataPlane->Cols();
	for (int32_t row = 0; row < rows; ++row)
	{
		for (int32_t col = 0; col < cols; ++col)
		{
			// Compute the expected value
			float expected = (float)(row*100 + col/2); // corner = col/2

			float value;
			CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, col, value));
			CPPUNIT_ASSERT(value == expected);
		}
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());

}

