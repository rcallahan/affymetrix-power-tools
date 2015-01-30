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
#include "calvin_files/writers/test/CalvinCHPQuantificationDetectionFileUpdaterTest.h"
//
#include "calvin_files/parsers/src/CHPQuantificationDetectionFileReader.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationDetectionFileUpdater.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationDetectionFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;

CPPUNIT_TEST_SUITE_REGISTRATION( CalvinCHPQuantificationDetectionFileUpdaterTest );

#define TEST_FILE "data_file.quantification_detection_updater_test"

void CalvinCHPQuantificationDetectionFileUpdaterTest::setUp()
{
}

void CalvinCHPQuantificationDetectionFileUpdaterTest::CreateReferenceFile()
{
	CHPQuantificationDetectionData data(TEST_FILE);
	data.SetEntryCount(4, 10);
	CHPQuantificationDetectionFileWriter *writer = new CHPQuantificationDetectionFileWriter(data);
	ProbeSetQuantificationDetectionData e;

	writer->SeekToDataSet();
	e.name = "1";
	e.quantification = 10.0f;
    e.pvalue = 0.1f;
	writer->WriteEntry(e);
	e.name = "2";
	e.quantification = 20.0f;
    e.pvalue = 0.2f;
	writer->WriteEntry(e);
	e.name = "3";
	e.quantification = 30.0f;
    e.pvalue = 0.3f;
	writer->WriteEntry(e);
	e.name = "4";
	e.quantification = 40.0f;
    e.pvalue = 0.4f;
	writer->WriteEntry(e);

	delete writer;
}

void CalvinCHPQuantificationDetectionFileUpdaterTest::tearDown()
{
}

void CalvinCHPQuantificationDetectionFileUpdaterTest::testQuantificationDetection()
{
	CreateReferenceFile();

	CalvinCHPQuantificationDetectionFileUpdater upd;
	upd.Initialize(TEST_FILE);
	upd.UpdateQuantification(0, 123.0f);
	upd.UpdateQuantification(2, 222.0f);
    upd.UpdateDetection(1, 0.22f);

	CHPQuantificationDetectionData data;
	CHPQuantificationDetectionFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(TEST_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
	CPPUNIT_ASSERT(data.GetEntryCount() == 4);

	ProbeSetQuantificationDetectionData e;
	data.GetQuantificationDetectionEntry(0, e);
	CPPUNIT_ASSERT(e.name == "1");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 123.0f, 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.pvalue, 0.1f, 0.000001f);
	data.GetQuantificationDetectionEntry(1, e);
	CPPUNIT_ASSERT(e.name == "2");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 20.0f, 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.pvalue, 0.22f, 0.000001f);
	data.GetQuantificationDetectionEntry(2, e);
	CPPUNIT_ASSERT(e.name == "3");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 222.0f, 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.pvalue, 0.3f, 0.000001f);
	data.GetQuantificationDetectionEntry(3, e);
	CPPUNIT_ASSERT(e.name == "4");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 40.0f, 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.pvalue, 0.4f, 0.000001f);
}
