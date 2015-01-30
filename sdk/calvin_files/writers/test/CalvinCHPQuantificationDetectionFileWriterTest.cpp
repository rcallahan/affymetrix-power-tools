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
#include "calvin_files/writers/test/CalvinCHPQuantificationDetectionFileWriterTest.h"
//
#include "calvin_files/data/src/ProbeSetQuantificationDetectionData.h"
#include "calvin_files/parsers/src/CHPQuantificationDetectionFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationDetectionFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPQuantificationDetectionFileWriterTest );

void CHPQuantificationDetectionFileWriterTest::setUp() {}

void CHPQuantificationDetectionFileWriterTest::tearDown(){}

void CHPQuantificationDetectionFileWriterTest::testCreation()
{
	CHPQuantificationDetectionData fHdr("CHP_quantification_detection_file_empty");
	CHPQuantificationDetectionFileWriter* w = new CHPQuantificationDetectionFileWriter(fHdr);
	CPPUNIT_ASSERT(1);
	delete w;
}

void CHPQuantificationDetectionFileWriterTest::WriteTest()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPQuantificationDetectionData data("CHP_quantification_detection_file");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(2, 10);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPQuantificationDetectionFileWriter *writer = new CHPQuantificationDetectionFileWriter(data);
	affymetrix_calvin_data::ProbeSetQuantificationDetectionData e;

	writer->SeekToDataSet();
	e.name = "abc";
	e.quantification = 10.0f;
    e.pvalue = 0.1f;
	writer->WriteEntry(e);
	e.name = "xyz";
	e.quantification = 20.0f;
    e.pvalue = 0.2f;
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);




	CHPQuantificationDetectionData data2;
	CHPQuantificationDetectionFileReader reader;
	reader.SetFilename("CHP_quantification_detection_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount() == 2);

	ParameterNameValueTypeList p = data2.GetAlgParams();
	ParameterNameValueTypeList::iterator it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	p = data2.GetSummaryParams();
	it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");


	data2.GetQuantificationDetectionEntry(0, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 10.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.pvalue, 0.10f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "abc");
    data2.GetQuantificationDetectionEntry(1, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 20.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.pvalue, 0.20f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "xyz");

}

void CHPQuantificationDetectionFileWriterTest::WriteTestId()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPQuantificationDetectionData data("CHP_quantification_detection_file_id");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(2);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPQuantificationDetectionFileWriter *writer = new CHPQuantificationDetectionFileWriter(data);
	affymetrix_calvin_data::ProbeSetQuantificationDetectionData e;

	writer->SeekToDataSet();
	e.id = 10;
	e.quantification = 10.0f;
    e.pvalue = 0.1f;
	writer->WriteEntry(e);
	e.id = 20;
	e.quantification = 20.0f;
    e.pvalue = 0.2f;
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);


	CHPQuantificationDetectionData data2;
	CHPQuantificationDetectionFileReader reader;
	reader.SetFilename("CHP_quantification_detection_file_id");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount() == 2);

	ParameterNameValueTypeList p = data2.GetAlgParams();
	ParameterNameValueTypeList::iterator it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	p = data2.GetSummaryParams();
	it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");


	data2.GetQuantificationDetectionEntry(0, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 10.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.pvalue, 0.1f, 0.0001f);
	CPPUNIT_ASSERT(e.id == 10);
	CPPUNIT_ASSERT(e.name == "");
	data2.GetQuantificationDetectionEntry(1, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 20.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.pvalue, 0.2f, 0.0001f);
	CPPUNIT_ASSERT(e.id == 20);
	CPPUNIT_ASSERT(e.name == "");

}
