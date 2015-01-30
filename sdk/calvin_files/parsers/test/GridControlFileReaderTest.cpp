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

#include "calvin_files/parsers/test/GridControlFileReaderTest.h"
//
#include "calvin_files/parsers/src/GridControlFileReader.h"
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_grid_control;
using namespace affymetrix_calvin_utilities;

#define TEST_FILE_NAME "../data/test.file.grc"

CPPUNIT_TEST_SUITE_REGISTRATION( GridControlFileReaderTest );

void GridControlFileReaderTest::setUp()
{
}

void GridControlFileReaderTest::tearDown()
{
}

void GridControlFileReaderTest::testCreation()
{
	GridControlFileReader *reader = new GridControlFileReader;
	CPPUNIT_ASSERT ( reader != NULL );
	delete reader;
}

void GridControlFileReaderTest::testRead()
{
	GridControlFileReader reader;
	GridControlData data;

	reader.Read(TEST_FILE_NAME, data);

	CPPUNIT_ASSERT(data.GetRows() == 5);
	CPPUNIT_ASSERT(data.GetColumns() == 5);
	CPPUNIT_ASSERT(data.GetNumB1Probes() == 3);
	CPPUNIT_ASSERT(data.GetNumB2Probes() == 2);
	CPPUNIT_ASSERT(data.GetNumNSProbes() == 1);

	CPPUNIT_ASSERT(data.GetB1(0).x == 0);
	CPPUNIT_ASSERT(data.GetB1(0).y == 1);
	CPPUNIT_ASSERT(data.GetB1(1).x == 2);
	CPPUNIT_ASSERT(data.GetB1(1).y == 3);
	CPPUNIT_ASSERT(data.GetB1(2).x == 4);
	CPPUNIT_ASSERT(data.GetB1(2).y == 5);

	CPPUNIT_ASSERT(data.GetB2(0).x == 10);
	CPPUNIT_ASSERT(data.GetB2(0).y == 11);
	CPPUNIT_ASSERT(data.GetB2(1).x == 12);
	CPPUNIT_ASSERT(data.GetB2(1).y == 13);

	CPPUNIT_ASSERT(data.GetNS(0).x == 10);
	CPPUNIT_ASSERT(data.GetNS(0).y == 20);
}

void GridControlFileReaderTest::testRead_no_file()
{
	GridControlFileReader reader;
	GridControlData data;

	CPPUNIT_ASSERT_THROW(reader.Read("no_file.grc", data), affymetrix_calvin_exceptions::FileNotFoundException);
}

void GridControlFileReaderTest::testRead_invalid_file_version()
{
	GridControlFileReader reader;
	GridControlData data;

	CPPUNIT_ASSERT_THROW(reader.Read("../data/test.file.data_invalid_version", data), affymetrix_calvin_exceptions::InvalidVersionException);
}

void GridControlFileReaderTest::testRead_invalid_file_type()
{
	GridControlFileReader reader;
	GridControlData data;

	CPPUNIT_ASSERT_THROW(reader.Read("../data/test.file.full_array_file", data), affymetrix_calvin_exceptions::InvalidFileTypeException);
}
