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
#include "calvin_files/writers/test/CalvinCHPQuantificationFileUpdaterTest.h"
//
#include "calvin_files/parsers/src/CHPQuantificationFileReader.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationFileUpdater.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;

CPPUNIT_TEST_SUITE_REGISTRATION( CalvinCHPQuantificationFileUpdaterTest );

#define TEST_FILE "data_file.quantificationupdater_test"

void CalvinCHPQuantificationFileUpdaterTest::setUp()
{
}

void CalvinCHPQuantificationFileUpdaterTest::CreateReferenceFile()
{
	CHPQuantificationData data(TEST_FILE);
	data.SetEntryCount(4, 10);
	CHPQuantificationFileWriter *writer = new CHPQuantificationFileWriter(data);
	ProbeSetQuantificationData e;

	writer->SeekToDataSet();
	e.name = "1";
	e.quantification = 10.0f;
	writer->WriteEntry(e);
	e.name = "2";
	e.quantification = 20.0f;
	writer->WriteEntry(e);
	e.name = "3";
	e.quantification = 30.0f;
	writer->WriteEntry(e);
	e.name = "4";
	e.quantification = 40.0f;
	writer->WriteEntry(e);

	delete writer;
}

void CalvinCHPQuantificationFileUpdaterTest::tearDown()
{
}

void CalvinCHPQuantificationFileUpdaterTest::testQuantification()
{
	CreateReferenceFile();

	CalvinCHPQuantificationFileUpdater upd;
	upd.Initialize(TEST_FILE);
	upd.UpdateQuantification(0, 123.0f);
	upd.UpdateQuantification(2, 222.0f);

	CHPQuantificationData data;
	CHPQuantificationFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(TEST_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
	CPPUNIT_ASSERT(data.GetEntryCount() == 4);

	ProbeSetQuantificationData e;
	data.GetQuantificationEntry(0, e);
	CPPUNIT_ASSERT(e.name == "1");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 123.0f, 0.000001f);
	data.GetQuantificationEntry(1, e);
	CPPUNIT_ASSERT(e.name == "2");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 20.0f, 0.000001f);
	data.GetQuantificationEntry(2, e);
	CPPUNIT_ASSERT(e.name == "3");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 222.0f, 0.000001f);
	data.GetQuantificationEntry(3, e);
	CPPUNIT_ASSERT(e.name == "4");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 40.0f, 0.000001f);
}
