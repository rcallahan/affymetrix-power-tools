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

#include "calvin_files/writers/test/CalvinCHPTilingFileWriterTest.h"
//
#include "calvin_files/data/src/CHPTilingEntry.h"
#include "calvin_files/parsers/src/CHPTilingFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/CalvinCHPTilingFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPTilingFileWriterTest );

void CHPTilingFileWriterTest::setUp() {}

void CHPTilingFileWriterTest::tearDown(){}

void CHPTilingFileWriterTest::testCreation()
{
	CHPTilingData fHdr("CHP_tiling_file_empty");
	CHPTilingFileWriter* w = new CHPTilingFileWriter(fHdr);
	CPPUNIT_ASSERT(1);
	delete w;
}

void CHPTilingFileWriterTest::WriteTest()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPTilingData data("CHP_tiling_file");

	data.SetAlgName(L"tile");
	data.SetAlgVersion(L"1.0");
	data.SetNumberSequences(2);
	param.SetName(L"p1");
	param.SetValueText(L"v1");
	params.push_back(param);
	data.AddAlgParams(params);

	TilingSequenceData seq;
	seq.name = L"n1";
	seq.groupName = L"g1";
	seq.version = L"v1";
	param.SetName(L"seq1_p1");
	param.SetValueText(L"seq1_v1");
	seq.parameters.push_back(param);
	data.AddTilingSequenceData(2, seq);

	seq.name = L"n2";
	seq.groupName = L"g2";
	seq.version = L"v2";
	param.SetName(L"seq2_p1");
	param.SetValueText(L"seq2_v1");
	seq.parameters.clear();
	seq.parameters.push_back(param);
	data.AddTilingSequenceData(3, seq);

	CHPTilingFileWriter *writer = new CHPTilingFileWriter(data);
	CHPTilingEntry e;

	writer->SeekToDataSet(0);
	e.position =10;
	e.value = 10.0;
	writer->WriteTilingEntry(e);
	e.position = 20;
	e.value = 20.0;
	writer->WriteTilingEntry(e);

	writer->SeekToDataSet(1);
	e.position =11;
	e.value = 11.0;
	writer->WriteTilingEntry(e);
	e.position = 21;
	e.value = 21.0;
	writer->WriteTilingEntry(e);
	e.position = 31;
	e.value = 31.0;
	writer->WriteTilingEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);




	CHPTilingData data2;
	CHPTilingFileReader reader;
	reader.SetFilename("CHP_tiling_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"tile");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetNumberSequences() == 2);

	params = data2.GetAlgParams();
	ParameterNameValueTypeList::iterator it = params.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"v1");
	CPPUNIT_ASSERT(data2.GetNumberSequences() == 2);

	data2.OpenTilingSequenceDataSet(0);
	seq = data2.GetTilingSequenceData();
	CPPUNIT_ASSERT(seq.name == L"n1");
	CPPUNIT_ASSERT(seq.groupName == L"g1");
	CPPUNIT_ASSERT(seq.version == L"v1");
	param = *seq.parameters.begin();
	CPPUNIT_ASSERT(param.GetName() == L"seq1_p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"seq1_v1");
	CPPUNIT_ASSERT(data2.GetTilingSequenceEntryCount(0) == 2);

	data2.GetTilingSequenceEntry(0, e);
	CPPUNIT_ASSERT(e.position == 10);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 10.0f, 0.00001f);
	data2.GetTilingSequenceEntry(1, e);
	CPPUNIT_ASSERT(e.position == 20);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 20.0f, 0.00001f);

	data2.OpenTilingSequenceDataSet(1);
	seq = data2.GetTilingSequenceData();
	CPPUNIT_ASSERT(seq.name == L"n2");
	CPPUNIT_ASSERT(seq.groupName == L"g2");
	CPPUNIT_ASSERT(seq.version == L"v2");
	param = *seq.parameters.begin();
	CPPUNIT_ASSERT(param.GetName() == L"seq2_p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"seq2_v1");
	CPPUNIT_ASSERT(data2.GetTilingSequenceEntryCount(1) == 3);

	data2.GetTilingSequenceEntry(0, e);
	CPPUNIT_ASSERT(e.position == 11);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 11.0f, 0.00001f);
	data2.GetTilingSequenceEntry(1, e);
	CPPUNIT_ASSERT(e.position == 21);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 21.0f, 0.00001f);
	data2.GetTilingSequenceEntry(2, e);
	CPPUNIT_ASSERT(e.position == 31);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 31.0f, 0.00001f);

}
