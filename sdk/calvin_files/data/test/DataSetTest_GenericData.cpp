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
#include "calvin_files/data/test/DataSetTest_GenericData.h"
//
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//

using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetTest_GenericData );

DataSetTest_GenericData::DataSetTest_GenericData()
{
}

DataSetTest_GenericData::~DataSetTest_GenericData()
{
}

void DataSetTest_GenericData::setUp()
{
}

void DataSetTest_GenericData::tearDown()
{
}

void DataSetTest_GenericData::ColsFromDataSetWithZeroRows()
{
	GenericData gd;
	GenericFileReader reader;
	reader.SetFilename("../data/small_cel_file_with_dataset_of_zero_rows");
	reader.ReadHeader(gd);

	// Open a DataSet that has 0 rows.
	DataSet* ds = gd.DataSet(0,3);
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Cols() == 2);
	CPPUNIT_ASSERT(ds->Rows() == 0);

	ds->Delete();

	// Open another DataSet that has 0 rows.
	ds = gd.DataSet(0,4);
	CPPUNIT_ASSERT(ds);
	CPPUNIT_ASSERT(ds->Cols() == 2);
	CPPUNIT_ASSERT(ds->Rows() == 0);

	ds->Delete();
}
