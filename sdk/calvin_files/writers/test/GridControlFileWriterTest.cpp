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

#include "calvin_files/writers/test/GridControlFileWriterTest.h"
//
#include "calvin_files/parsers/src/GridControlFileReader.h"
#include "calvin_files/writers/src/GridControlFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_grid_control;
using namespace affymetrix_calvin_utilities;

#define TEST_FILE_NAME "./test.new.grc"

CPPUNIT_TEST_SUITE_REGISTRATION( GridControlFileWriterTest );

void GridControlFileWriterTest::setUp()
{
}

void GridControlFileWriterTest::tearDown()
{
}

void GridControlFileWriterTest::testCreation()
{
	GridControlFileWriter *writer = new GridControlFileWriter;
	CPPUNIT_ASSERT ( writer != NULL );
	delete writer;
}

void GridControlFileWriterTest::WriteTest()
{
	GridControlData data;
	FeatureCoordinate coord;
	GridControlFileReader reader;
	GridControlData dataRead;
	GridControlFileWriter writer;

	data.SetRows(5);
	data.SetColumns(5);
	data.ResizeB1(3);
	data.ResizeB2(2);
	data.ResizeNS(1);

	coord.x = 0;
	coord.y = 1;
	data.SetB1(0, coord);
	coord.x = 2;
	coord.y = 3;
	data.SetB1(1, coord);
	coord.x = 4;
	coord.y = 5;
	data.SetB1(2, coord);

	coord.x = 10;
	coord.y = 11;
	data.SetB2(0, coord);
	coord.x = 12;
	coord.y = 13;
	data.SetB2(1, coord);

	coord.x = 10;
	coord.y = 20;
	data.SetNS(0, coord);

	CPPUNIT_ASSERT(writer.Write(TEST_FILE_NAME, data) == true);



	reader.Read(TEST_FILE_NAME, dataRead);

	CPPUNIT_ASSERT(dataRead.GetRows() == data.GetRows());
	CPPUNIT_ASSERT(dataRead.GetColumns() == data.GetColumns());
	CPPUNIT_ASSERT(dataRead.GetNumB1Probes() == data.GetNumB1Probes());
	CPPUNIT_ASSERT(dataRead.GetNumB2Probes() == data.GetNumB2Probes());
	CPPUNIT_ASSERT(dataRead.GetNumNSProbes() == data.GetNumNSProbes());

	CPPUNIT_ASSERT(dataRead.GetB1(0).x == data.GetB1(0).x);
	CPPUNIT_ASSERT(dataRead.GetB1(0).y == data.GetB1(0).y);
	CPPUNIT_ASSERT(dataRead.GetB1(1).x == data.GetB1(1).x);
	CPPUNIT_ASSERT(dataRead.GetB1(1).y == data.GetB1(1).y);
	CPPUNIT_ASSERT(dataRead.GetB1(2).x == data.GetB1(2).x);
	CPPUNIT_ASSERT(dataRead.GetB1(2).y == data.GetB1(2).y);

	CPPUNIT_ASSERT(dataRead.GetB2(0).x == data.GetB2(0).x);
	CPPUNIT_ASSERT(dataRead.GetB2(0).y == data.GetB2(0).y);
	CPPUNIT_ASSERT(dataRead.GetB2(1).x == data.GetB2(1).x);
	CPPUNIT_ASSERT(dataRead.GetB2(1).y == data.GetB2(1).y);

	CPPUNIT_ASSERT(dataRead.GetNS(0).x == data.GetNS(0).x);
	CPPUNIT_ASSERT(dataRead.GetNS(0).y == data.GetNS(0).y);

}
