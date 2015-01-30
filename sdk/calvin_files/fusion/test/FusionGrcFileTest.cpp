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


#include "calvin_files/fusion/test/FusionGrcFileTest.h"
//
#include "calvin_files/fusion/src/FusionGrcFileReader.h"
//

using namespace affymetrix_fusion_io;
using namespace affymetrix_grid_control;
using namespace std;

#define INVALID_FILE "../../parsers/data/test.file.full_array_file"
#define CALVIN_FILE "../../parsers/data/test.file.grc"
#define GCOS_FILE "../../../file/CPPTest/data/test.grc"

CPPUNIT_TEST_SUITE_REGISTRATION( FusionGrcFileReaderTest );

void FusionGrcFileReaderTest::setUp()
{
}

void FusionGrcFileReaderTest::tearDown()
{
}

void FusionGrcFileReaderTest::testCreation()
{
	FusionGrcFileReader *reader = new FusionGrcFileReader;
	CPPUNIT_ASSERT(reader != NULL);
	delete reader;
}

void FusionGrcFileReaderTest::testmethod_Read_invalid_file()
{
	FusionGrcFileReader reader;
	GridControlData data;
	//CPPUNIT_ASSERT( reader.Read(INVALID_FILE, data) == false);
}

void FusionGrcFileReaderTest::testmethod_Read_missing_file()
{
	FusionGrcFileReader reader;
	GridControlData data;
	CPPUNIT_ASSERT( reader.Read("no_file", data) == false);
}

void FusionGrcFileReaderTest::testmethod_Read_calvin_grc_file()
{
	FusionGrcFileReader reader;
	GridControlData data;
	CPPUNIT_ASSERT( reader.Read(CALVIN_FILE, data) == true);

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

void FusionGrcFileReaderTest::testmethod_Read_gcos_grc_file()
{
	FusionGrcFileReader reader;
	GridControlData data;
	CPPUNIT_ASSERT( reader.Read(GCOS_FILE, data) == true);

	CPPUNIT_ASSERT( data.GetRows() == 2) ;
	CPPUNIT_ASSERT( data.GetColumns() == 2) ;

	CPPUNIT_ASSERT( data.GetNumB1Probes() == 1) ;
	CPPUNIT_ASSERT( data.GetB1(0).x == 1);
	CPPUNIT_ASSERT( data.GetB1(0).y == 0);

	CPPUNIT_ASSERT( data.GetNumB2Probes() == 2) ;
	CPPUNIT_ASSERT( data.GetB2(0).x == 0);
	CPPUNIT_ASSERT( data.GetB2(0).y == 0);
	CPPUNIT_ASSERT( data.GetB2(1).x == 0);
	CPPUNIT_ASSERT( data.GetB2(1).y == 1);

	CPPUNIT_ASSERT( data.GetNumNSProbes() == 1) ;
	CPPUNIT_ASSERT( data.GetNS(0).x == 1);
	CPPUNIT_ASSERT( data.GetNS(0).y == 1);

}
