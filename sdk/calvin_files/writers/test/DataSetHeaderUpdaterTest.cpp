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
#include "calvin_files/writers/test/DataSetHeaderUpdaterTest.h"
//
#include "calvin_files/parsers/src/CHPQuantificationDetectionFileReader.h"
#include "calvin_files/parsers/src/CHPQuantificationFileReader.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationDetectionFileWriter.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationFileWriter.h"
#include "calvin_files/writers/src/DataSetHeaderUpdater.h"
//
#include "util/Fs.h"
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;

CPPUNIT_TEST_SUITE_REGISTRATION( DataSetHeaderUpdaterTest );

#define TEST_FILE "data_file.quantification_updater_test"

void DataSetHeaderUpdaterTest::setUp()
{
}

void DataSetHeaderUpdaterTest::CreateReferenceQDFile()
{
	CHPQuantificationDetectionData data(TEST_FILE);
	data.SetEntryCount(4, 10);

    ParameterNameValueType p;
    p.SetName(L"affymetrix-dataset-md5");
    p.SetValueAscii("", 128);
    data.GetDataSetHeader().AddNameValParam(p);
    p.SetName(L"other-param");
    p.SetValueAscii("", 128);
    data.GetDataSetHeader().AddNameValParam(p);

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

void DataSetHeaderUpdaterTest::CreateReferenceQFile()
{
	CHPQuantificationData data(TEST_FILE);
	data.SetEntryCount(4, 10);

    ParameterNameValueType p;
    p.SetName(L"affymetrix-dataset-md5");
    p.SetValueAscii("", 128);
    data.GetDataSetHeader().AddNameValParam(p);
    p.SetName(L"other-param");
    p.SetValueAscii("", 128);
    data.GetDataSetHeader().AddNameValParam(p);

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

void DataSetHeaderUpdaterTest::tearDown()
{
}

void DataSetHeaderUpdaterTest::testUpdateQDFile()
{
	CreateReferenceQDFile();

	CHPQuantificationDetectionData data;
	CHPQuantificationDetectionFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(TEST_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
    DataSetHeader hdr = data.GetDataSetHeader();
    data.Clear();

   	std::ofstream os;
        Fs::aptOpen(os, TEST_FILE, std::ios::out|std::ios::binary|std::ios::in);
    DataSetHeaderUpdater upd(os, hdr);
    
    ParameterNameValueType nvt;
    nvt.SetName(L"affymetrix-dataset-md5");
    nvt.SetValueAscii("abc", 128);
    CPPUNIT_ASSERT(upd.UpdateParameter(nvt) == true);
    nvt.SetName(L"other-param");
    nvt.SetValueAscii("xyz", 128);
    CPPUNIT_ASSERT(upd.UpdateParameter(nvt) == true);
    os.close();


    nvt.SetName(L"");
    nvt.SetValueAscii("");

	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
    hdr = data.GetDataSetHeader();
    CPPUNIT_ASSERT(hdr.FindNameValParam(L"other-param", nvt) == true);
    CPPUNIT_ASSERT(nvt.GetValueAscii() == "xyz");
    CPPUNIT_ASSERT(hdr.FindNameValParam(L"affymetrix-dataset-md5", nvt) == true);
    CPPUNIT_ASSERT(nvt.GetValueAscii() == "abc");


}

void DataSetHeaderUpdaterTest::testUpdateQFile()
{
	CreateReferenceQFile();

	CHPQuantificationData data;
	CHPQuantificationFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(TEST_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
    DataSetHeader hdr = data.GetDataSetHeader();
    data.Clear();

   	std::ofstream os;
        Fs::aptOpen(os, TEST_FILE, std::ios::out|std::ios::binary|std::ios::in);
    DataSetHeaderUpdater upd(os, hdr);
    
    ParameterNameValueType nvt;
    nvt.SetName(L"affymetrix-dataset-md5");
    nvt.SetValueAscii("abc", 128);
    CPPUNIT_ASSERT(upd.UpdateParameter(nvt) == true);
    nvt.SetName(L"other-param");
    nvt.SetValueAscii("xyz", 128);
    CPPUNIT_ASSERT(upd.UpdateParameter(nvt) == true);
    os.close();


    nvt.SetName(L"");
    nvt.SetValueAscii("");

	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
    hdr = data.GetDataSetHeader();
    CPPUNIT_ASSERT(hdr.FindNameValParam(L"other-param", nvt) == true);
    CPPUNIT_ASSERT(nvt.GetValueAscii() == "xyz");
    CPPUNIT_ASSERT(hdr.FindNameValParam(L"affymetrix-dataset-md5", nvt) == true);
    CPPUNIT_ASSERT(nvt.GetValueAscii() == "abc");


}
