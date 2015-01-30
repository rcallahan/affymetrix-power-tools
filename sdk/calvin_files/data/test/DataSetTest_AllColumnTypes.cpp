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
#include "calvin_files/data/test/DataSetTest_AllColumnTypes.h"
//
#include "calvin_files/data/src/DataException.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstring>
#include <string>
//
#include <stdio.h>

using namespace std;
using namespace affymetrix_calvin_io;

#define TEST_DATA_ALL_COLUMN_TYPES "../data/test.file.data_all_column_types"

#ifdef _MSC_VER
#pragma warning(disable: 4996) // ignore deprecated functions warning
#endif

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetTest_AllColumnTypes );

DataSetTest_AllColumnTypes::DataSetTest_AllColumnTypes()
{
	data = 0;
	dataPlane = 0;
}

DataSetTest_AllColumnTypes::~DataSetTest_AllColumnTypes()
{

}

void DataSetTest_AllColumnTypes::setUp()
{
	// Open one DataSet.
	data = new GenericData;
	GenericFileReader reader;
	std::string name = TEST_DATA_ALL_COLUMN_TYPES;
	reader.SetFilename(name);
	reader.ReadHeader(*data);
	try
	{
		dataPlane = data->DataSet(0,0);	// there is only one DataSet
	}
	catch(affymetrix_calvin_exceptions::CalvinException&)
	{
		dataPlane = 0;
	}
}

void DataSetTest_AllColumnTypes::tearDown()
{
	if (dataPlane)
		dataPlane->Delete();
	delete data;
}

void DataSetTest_AllColumnTypes::testCreation()
{
	// Check that the DataSet was created.
	CPPUNIT_ASSERT(dataPlane);
}

void DataSetTest_AllColumnTypes::testOpenDataSet()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);
	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());
}


void DataSetTest_AllColumnTypes::testAccessColumnValues()
{
	CPPUNIT_ASSERT(dataPlane->Open() == true);

	int32_t rows = dataPlane->Rows();
	//int32_t cols = dataPlane->Cols();

	// types to receive values.
	int8_t b;
	u_int8_t ub;
	int16_t s;
	u_int16_t us;
	int32_t i;
	u_int32_t ui;
	std::string str;
	std::wstring wstr;
	float f;

	for (int32_t row = 0; row < rows; ++row)
	{
		// expected values
		char expected_str[10];
		wchar_t expected_wstr[15];
		int8_t expected_b = 1+10*row;
		u_int8_t expected_ub = 2+10*row;
		sprintf(expected_str, "%d", 3+10*row);
		int16_t expected_s = 4+10*row;
		u_int16_t expected_us = 5+10*row;
		int32_t expected_i = 6+10*row;
		u_int32_t expected_ui = 7+10*row;
		FormatString1(expected_wstr, 15, L"%d", 8+10*row);
		float expected_f = (float)(9+10*row);

		// Get and check the values
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 0, b));
		CPPUNIT_ASSERT(b == expected_b);
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 1, ub));
		CPPUNIT_ASSERT(ub == expected_ub);
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 2, str));
		CPPUNIT_ASSERT(str == expected_str);
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 3, s));
		CPPUNIT_ASSERT(s == expected_s);
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 4, us));
		CPPUNIT_ASSERT(us == expected_us);
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 5, i));
		CPPUNIT_ASSERT(i == expected_i);
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 6, ui));
		CPPUNIT_ASSERT(ui == expected_ui);
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 7, wstr));
		CPPUNIT_ASSERT(wstr == expected_wstr);
		CPPUNIT_ASSERT_NO_THROW(dataPlane->GetData(row, 8, f));
		CPPUNIT_ASSERT(f == expected_f);
	}

	CPPUNIT_ASSERT_NO_THROW(dataPlane->Close());

}

