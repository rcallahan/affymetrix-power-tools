////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
#include "calvin_files/writers/test/CalvinCHPMultiDataFileUpdaterTest.h"
//
#include "calvin_files/parsers/src/CHPMultiDataFileReader.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileUpdater.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;

CPPUNIT_TEST_SUITE_REGISTRATION( CalvinCHPMultiDataFileUpdaterTest );

#define TEST1_FILE "data_file.MultiDataupdater_test1"
#define TEST2_FILE "data_file.MultiDataupdater_test2"
#define TEST3_FILE "data_file.MultiDataupdater_test3"

void CalvinCHPMultiDataFileUpdaterTest::setUp()
{
}

void CalvinCHPMultiDataFileUpdaterTest::CreateReferenceFile1()
{
	CHPMultiDataData data(TEST1_FILE);
	data.SetEntryCount(GenotypeMultiDataType, 4, 10);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	ProbeSetMultiDataGenotypeData e;

	writer->SeekToDataSet(GenotypeMultiDataType);
	e.name = "1";
    e.call = 1;
    e.confidence = 10.0f;
	writer->WriteEntry(e);
	e.name = "2";
    e.call = 2;
    e.confidence = 20.0f;
	writer->WriteEntry(e);
	e.name = "3";
    e.call = 3;
    e.confidence = 30.0f;
	writer->WriteEntry(e);
	e.name = "4";
    e.call = 4;
    e.confidence = 40.0f;
	writer->WriteEntry(e);

	delete writer;
}

void CalvinCHPMultiDataFileUpdaterTest::CreateReferenceFile2()
{
	CHPMultiDataData data(TEST2_FILE);
    vector<ColumnInfo> cols;
    ParameterNameValueType nv;

    ByteColumn bcol(L"byte");
    cols.push_back(bcol);

    UByteColumn ubcol(L"ubyte");
    cols.push_back(ubcol);

    ShortColumn scol(L"short");
    cols.push_back(scol);

    UShortColumn uscol(L"ushort");
    cols.push_back(uscol);

    IntColumn icol(L"int");
    cols.push_back(icol);

    UIntColumn uicol(L"uint");
    cols.push_back(uicol);

    FloatColumn fcol(L"float");
    cols.push_back(fcol);

    ASCIIColumn acol(L"ascii", 7);
    cols.push_back(acol);

    UnicodeColumn tcol(L"text", 10);
    cols.push_back(tcol);


	ProbeSetMultiDataGenotypeData e;
	data.SetEntryCount(GenotypeMultiDataType, 4, 10, cols);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);

    nv.SetName(L"byte");
    nv.SetValueInt8(8);
    e.metrics.push_back(nv);
    nv.SetName(L"ubyte");
    nv.SetValueUInt8(8);
    e.metrics.push_back(nv);

    nv.SetName(L"short");
    nv.SetValueInt16(16);
    e.metrics.push_back(nv);
    nv.SetName(L"ushort");
    nv.SetValueUInt16(16);
    e.metrics.push_back(nv);

    nv.SetName(L"int");
    nv.SetValueInt32(32);
    e.metrics.push_back(nv);
    nv.SetName(L"uint");
    nv.SetValueUInt32(32);
    e.metrics.push_back(nv);

    nv.SetName(L"float");
    nv.SetValueFloat(44.0f);
    e.metrics.push_back(nv);

    nv.SetName(L"ascii");
    nv.SetValueAscii("ascii");
    e.metrics.push_back(nv);

    nv.SetName(L"text");
    nv.SetValueText(L"text");
    e.metrics.push_back(nv);


	writer->SeekToDataSet(GenotypeMultiDataType);
	e.name = "1";
    e.call = 1;
    e.confidence = 10.0f;
	writer->WriteEntry(e);
	e.name = "2";
    e.call = 2;
    e.confidence = 20.0f;
	writer->WriteEntry(e);
	e.name = "3";
    e.call = 3;
    e.confidence = 30.0f;
	writer->WriteEntry(e);
	e.name = "4";
    e.call = 4;
    e.confidence = 40.0f;
	writer->WriteEntry(e);

	delete writer;
}

void CalvinCHPMultiDataFileUpdaterTest::CreateReferenceFile3()
{
	CHPMultiDataData data(TEST3_FILE);
    vector<ColumnInfo> cols;
    ParameterNameValueType nv;

    ByteColumn bcol(L"byte");
    cols.push_back(bcol);

    UByteColumn ubcol(L"ubyte");
    cols.push_back(ubcol);

    ShortColumn scol(L"short");
    cols.push_back(scol);

    UShortColumn uscol(L"ushort");
    cols.push_back(uscol);

    IntColumn icol(L"int");
    cols.push_back(icol);

    UIntColumn uicol(L"uint");
    cols.push_back(uicol);

    FloatColumn fcol(L"float");
    cols.push_back(fcol);

    ASCIIColumn acol(L"ascii", 7);
    cols.push_back(acol);

    UnicodeColumn tcol(L"text", 10);
    cols.push_back(tcol);


	ProbeSetMultiDataCopyNumberData e;
	ProbeSetMultiDataCytoRegionData c;
    ProbeSetMultiDataCopyNumberVariationRegionData v;
	data.SetEntryCount(CopyNumberMultiDataType, 4, 10, cols);
	data.SetEntryCount(CytoMultiDataType, 2, 10);
    data.SetEntryCount(CopyNumberVariationMultiDataType, 2, 10);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);

    nv.SetName(L"byte");
    nv.SetValueInt8(8);
    e.metrics.push_back(nv);
    nv.SetName(L"ubyte");
    nv.SetValueUInt8(8);
    e.metrics.push_back(nv);

    nv.SetName(L"short");
    nv.SetValueInt16(16);
    e.metrics.push_back(nv);
    nv.SetName(L"ushort");
    nv.SetValueUInt16(16);
    e.metrics.push_back(nv);

    nv.SetName(L"int");
    nv.SetValueInt32(32);
    e.metrics.push_back(nv);
    nv.SetName(L"uint");
    nv.SetValueUInt32(32);
    e.metrics.push_back(nv);

    nv.SetName(L"float");
    nv.SetValueFloat(44.0f);
    e.metrics.push_back(nv);

    nv.SetName(L"ascii");
    nv.SetValueAscii("ascii");
    e.metrics.push_back(nv);

    nv.SetName(L"text");
    nv.SetValueText(L"text");
    e.metrics.push_back(nv);


	writer->SeekToDataSet(CopyNumberMultiDataType);
	e.name = "1";
    e.chr = 1;
    e.position = 10;
	writer->WriteEntry(e);
	e.name = "2";
    e.chr = 2;
    e.position = 20;
	writer->WriteEntry(e);
	e.name = "3";
    e.chr = 3;
    e.position = 30;
	writer->WriteEntry(e);
	e.name = "4";
    e.chr = 4;
    e.position = 40;
	writer->WriteEntry(e);

	writer->SeekToDataSet(CytoMultiDataType);
	c.name = "1";
    c.chr= 1;
    c.startPosition = 1;
    c.stopPosition = 2;
    c.call = 1;
    c.confidenceScore = 10.0f;
	writer->WriteEntry(c);
	c.name = "2";
    c.chr= 2;
    c.startPosition = 2;
    c.stopPosition = 3;
    c.call = 2;
    c.confidenceScore = 20.0f;
	writer->WriteEntry(c);

    writer->SeekToDataSet(CopyNumberVariationMultiDataType);
	v.name = "1";
    v.signal = 1;
    v.call = 1;
    v.confidenceScore = 10.0f;
	writer->WriteEntry(v);
	v.name = "2";
    v.signal = 2;
    v.call = 2;
    v.confidenceScore = 20.0f;
	writer->WriteEntry(v);

	delete writer;
}

void CalvinCHPMultiDataFileUpdaterTest::tearDown()
{
}

void CalvinCHPMultiDataFileUpdaterTest::testMultiData1()
{
	CreateReferenceFile1();

    vector<ColumnInfo> cols;
	CalvinCHPMultiDataFileUpdater upd;
    ProbeSetMultiDataGenotypeData e;
	upd.Initialize(TEST1_FILE);

    e.call = 11;
    e.confidence = 111.0f;
	upd.UpdateMultiData(GenotypeMultiDataType, 0, e);

    e.call = 22;
    e.confidence = 222.0f;
	upd.UpdateMultiData(GenotypeMultiDataType, 2, e);

	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(TEST1_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
	CPPUNIT_ASSERT(data.GetEntryCount(GenotypeMultiDataType) == 4);

    data.GetGenotypeEntry(GenotypeMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.name == "1");
	CPPUNIT_ASSERT(e.call == 11);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 111.0f, 0.000001f);
	data.GetGenotypeEntry(GenotypeMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.name == "2");
	CPPUNIT_ASSERT(e.call == 2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 20.0f, 0.000001f);
	data.GetGenotypeEntry(GenotypeMultiDataType, 2, e);
	CPPUNIT_ASSERT(e.name == "3");
	CPPUNIT_ASSERT(e.call == 22);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 222.0f, 0.000001f);
	data.GetGenotypeEntry(GenotypeMultiDataType, 3, e);
	CPPUNIT_ASSERT(e.name == "4");
	CPPUNIT_ASSERT(e.call == 4);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 40.0f, 0.000001f);
}

void CalvinCHPMultiDataFileUpdaterTest::testMultiData2()
{
	CreateReferenceFile2();

    vector<ColumnInfo> cols;
    ByteColumn bcol(L"byte");
    cols.push_back(bcol);

    UByteColumn ubcol(L"ubyte");
    cols.push_back(ubcol);

    ShortColumn scol(L"short");
    cols.push_back(scol);

    UShortColumn uscol(L"ushort");
    cols.push_back(uscol);

    IntColumn icol(L"int");
    cols.push_back(icol);

    UIntColumn uicol(L"uint");
    cols.push_back(uicol);

    FloatColumn fcol(L"float");
    cols.push_back(fcol);

    ASCIIColumn acol(L"ascii", 7);
    cols.push_back(acol);

    UnicodeColumn tcol(L"text", 10);
    cols.push_back(tcol);

	CalvinCHPMultiDataFileUpdater upd;
    ProbeSetMultiDataGenotypeData e;
	upd.Initialize(TEST2_FILE);

    e.metrics.resize(9);

    e.call = 11;
    e.confidence = 111.0f;
    e.metrics[0].SetValueInt8(9);
    e.metrics[1].SetValueUInt8(10);
    e.metrics[2].SetValueInt16(17);
    e.metrics[3].SetValueUInt16(18);
    e.metrics[4].SetValueInt32(33);
    e.metrics[5].SetValueUInt32(34);
    e.metrics[6].SetValueFloat(55.0f);
    e.metrics[7].SetValueAscii("text");
    e.metrics[8].SetValueText(L"ascii");
	upd.UpdateMultiData(GenotypeMultiDataType, 0, e, cols);

    e.call = 22;
    e.confidence = 222.0f;
    e.metrics[0].SetValueInt8(10);
    e.metrics[1].SetValueUInt8(11);
    e.metrics[2].SetValueInt16(18);
    e.metrics[3].SetValueUInt16(19);
    e.metrics[4].SetValueInt32(34);
    e.metrics[5].SetValueUInt32(35);
    e.metrics[6].SetValueFloat(66.0f);
    e.metrics[7].SetValueAscii("text2");
    e.metrics[8].SetValueText(L"ascii2");
	upd.UpdateMultiData(GenotypeMultiDataType, 2, e, cols);

	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(TEST2_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
	CPPUNIT_ASSERT(data.GetEntryCount(GenotypeMultiDataType) == 4);

    data.GetGenotypeEntry(GenotypeMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.name == "1");
	CPPUNIT_ASSERT(e.call == 11);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 111.0f, 0.000001f);
    CPPUNIT_ASSERT(e.metrics[0].GetValueInt8() == 9);
    CPPUNIT_ASSERT(e.metrics[1].GetValueUInt8() == 10);
    CPPUNIT_ASSERT(e.metrics[2].GetValueInt16() == 17);
    CPPUNIT_ASSERT(e.metrics[3].GetValueUInt16() == 18);
    CPPUNIT_ASSERT(e.metrics[4].GetValueInt32() == 33);
    CPPUNIT_ASSERT(e.metrics[5].GetValueUInt32() == 34);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.metrics[6].GetValueFloat(), 55.0f, 0.00001f);
    CPPUNIT_ASSERT(e.metrics[7].GetValueAscii() == "text");
    CPPUNIT_ASSERT(e.metrics[8].GetValueText() == L"ascii");

    data.GetGenotypeEntry(GenotypeMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.name == "2");
	CPPUNIT_ASSERT(e.call == 2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 20.0f, 0.000001f);
    CPPUNIT_ASSERT(e.metrics[0].GetValueInt8() == 8);
    CPPUNIT_ASSERT(e.metrics[1].GetValueUInt8() == 8);
    CPPUNIT_ASSERT(e.metrics[2].GetValueInt16() == 16);
    CPPUNIT_ASSERT(e.metrics[3].GetValueUInt16() == 16);
    CPPUNIT_ASSERT(e.metrics[4].GetValueInt32() == 32);
    CPPUNIT_ASSERT(e.metrics[5].GetValueUInt32() == 32);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.metrics[6].GetValueFloat(), 44.0f, 0.00001f);
    CPPUNIT_ASSERT(e.metrics[7].GetValueAscii() == "ascii");
    CPPUNIT_ASSERT(e.metrics[8].GetValueText() == L"text");

    data.GetGenotypeEntry(GenotypeMultiDataType, 2, e);
	CPPUNIT_ASSERT(e.name == "3");
	CPPUNIT_ASSERT(e.call == 22);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 222.0f, 0.000001f);
    CPPUNIT_ASSERT(e.metrics[0].GetValueInt8() == 10);
    CPPUNIT_ASSERT(e.metrics[1].GetValueUInt8() == 11);
    CPPUNIT_ASSERT(e.metrics[2].GetValueInt16() == 18);
    CPPUNIT_ASSERT(e.metrics[3].GetValueUInt16() == 19);
    CPPUNIT_ASSERT(e.metrics[4].GetValueInt32() == 34);
    CPPUNIT_ASSERT(e.metrics[5].GetValueUInt32() == 35);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.metrics[6].GetValueFloat(), 66.0f, 0.00001f);
    CPPUNIT_ASSERT(e.metrics[7].GetValueAscii() == "text2");
    CPPUNIT_ASSERT(e.metrics[8].GetValueText() == L"ascii2");

    data.GetGenotypeEntry(GenotypeMultiDataType, 3, e);
	CPPUNIT_ASSERT(e.name == "4");
	CPPUNIT_ASSERT(e.call == 4);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 40.0f, 0.000001f);
    CPPUNIT_ASSERT(e.metrics[0].GetValueInt8() == 8);
    CPPUNIT_ASSERT(e.metrics[1].GetValueUInt8() == 8);
    CPPUNIT_ASSERT(e.metrics[2].GetValueInt16() == 16);
    CPPUNIT_ASSERT(e.metrics[3].GetValueUInt16() == 16);
    CPPUNIT_ASSERT(e.metrics[4].GetValueInt32() == 32);
    CPPUNIT_ASSERT(e.metrics[5].GetValueUInt32() == 32);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.metrics[6].GetValueFloat(), 44.0f, 0.00001f);
    CPPUNIT_ASSERT(e.metrics[7].GetValueAscii() == "ascii");
    CPPUNIT_ASSERT(e.metrics[8].GetValueText() == L"text");
}

void CalvinCHPMultiDataFileUpdaterTest::testMultiData3()
{
	CreateReferenceFile3();

    vector<ColumnInfo> cols;
    ByteColumn bcol(L"byte");
    cols.push_back(bcol);

    UByteColumn ubcol(L"ubyte");
    cols.push_back(ubcol);

    ShortColumn scol(L"short");
    cols.push_back(scol);

    UShortColumn uscol(L"ushort");
    cols.push_back(uscol);

    IntColumn icol(L"int");
    cols.push_back(icol);

    UIntColumn uicol(L"uint");
    cols.push_back(uicol);

    FloatColumn fcol(L"float");
    cols.push_back(fcol);

    ASCIIColumn acol(L"ascii", 7);
    cols.push_back(acol);

    UnicodeColumn tcol(L"text", 10);
    cols.push_back(tcol);

	CalvinCHPMultiDataFileUpdater upd;
    ProbeSetMultiDataCopyNumberData e;
	ProbeSetMultiDataCytoRegionData c;
    ProbeSetMultiDataCopyNumberVariationRegionData v;
	upd.Initialize(TEST3_FILE);

    e.metrics.resize(9);

    e.chr = 11;
    e.position = 111;
    e.metrics[0].SetValueInt8(9);
    e.metrics[1].SetValueUInt8(10);
    e.metrics[2].SetValueInt16(17);
    e.metrics[3].SetValueUInt16(18);
    e.metrics[4].SetValueInt32(33);
    e.metrics[5].SetValueUInt32(34);
    e.metrics[6].SetValueFloat(55.0f);
    e.metrics[7].SetValueAscii("text");
    e.metrics[8].SetValueText(L"ascii");
	upd.UpdateMultiData(CopyNumberMultiDataType, 0, e, cols);

    e.chr = 22;
    e.position = 222;
    e.metrics[0].SetValueInt8(10);
    e.metrics[1].SetValueUInt8(11);
    e.metrics[2].SetValueInt16(18);
    e.metrics[3].SetValueUInt16(19);
    e.metrics[4].SetValueInt32(34);
    e.metrics[5].SetValueUInt32(35);
    e.metrics[6].SetValueFloat(66.0f);
    e.metrics[7].SetValueAscii("text2");
    e.metrics[8].SetValueText(L"ascii2");
	upd.UpdateMultiData(CopyNumberMultiDataType, 2, e, cols);

    vector<ColumnInfo> cycols;
    c.call = 11;
    c.confidenceScore = 111.0f;
    c.chr = 11;
    c.startPosition = 11;
    c.stopPosition = 22;
    upd.UpdateMultiData(CytoMultiDataType, 0, c, cycols);

    c.call = 2;
    c.confidenceScore = 20.0f;
    c.chr = 22;
    c.startPosition = 22;
    c.stopPosition = 33;
    upd.UpdateMultiData(CytoMultiDataType, 1, c, cycols);

    vector<ColumnInfo> cnvcols;
    v.call = 11;
    v.confidenceScore = 111.0f;
    v.signal = 111.0f;
    upd.UpdateMultiData(CopyNumberVariationMultiDataType, 0, v, cnvcols);

    v.call = 22;
    v.confidenceScore = 222.0f;
    v.signal = 222.0f;
    upd.UpdateMultiData(CopyNumberVariationMultiDataType, 1,v, cnvcols);

	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename(TEST3_FILE));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));
	CPPUNIT_ASSERT(data.GetEntryCount(CopyNumberMultiDataType) == 4);
	CPPUNIT_ASSERT(data.GetEntryCount(CytoMultiDataType) == 2);
    CPPUNIT_ASSERT(data.GetEntryCount(CopyNumberVariationMultiDataType) == 2);

    data.GetCopyNumberEntry(CopyNumberMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.name == "1");
	CPPUNIT_ASSERT(e.chr == 11);
    CPPUNIT_ASSERT(e.position == 111);
    CPPUNIT_ASSERT(e.metrics[0].GetValueInt8() == 9);
    CPPUNIT_ASSERT(e.metrics[1].GetValueUInt8() == 10);
    CPPUNIT_ASSERT(e.metrics[2].GetValueInt16() == 17);
    CPPUNIT_ASSERT(e.metrics[3].GetValueUInt16() == 18);
    CPPUNIT_ASSERT(e.metrics[4].GetValueInt32() == 33);
    CPPUNIT_ASSERT(e.metrics[5].GetValueUInt32() == 34);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.metrics[6].GetValueFloat(), 55.0f, 0.00001f);
    CPPUNIT_ASSERT(e.metrics[7].GetValueAscii() == "text");
    CPPUNIT_ASSERT(e.metrics[8].GetValueText() == L"ascii");

    data.GetCopyNumberEntry(CopyNumberMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.name == "2");
	CPPUNIT_ASSERT(e.chr == 2);
	CPPUNIT_ASSERT(e.position == 20);
    CPPUNIT_ASSERT(e.metrics[0].GetValueInt8() == 8);
    CPPUNIT_ASSERT(e.metrics[1].GetValueUInt8() == 8);
    CPPUNIT_ASSERT(e.metrics[2].GetValueInt16() == 16);
    CPPUNIT_ASSERT(e.metrics[3].GetValueUInt16() == 16);
    CPPUNIT_ASSERT(e.metrics[4].GetValueInt32() == 32);
    CPPUNIT_ASSERT(e.metrics[5].GetValueUInt32() == 32);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.metrics[6].GetValueFloat(), 44.0f, 0.00001f);
    CPPUNIT_ASSERT(e.metrics[7].GetValueAscii() == "ascii");
    CPPUNIT_ASSERT(e.metrics[8].GetValueText() == L"text");

    data.GetCopyNumberEntry(CopyNumberMultiDataType, 2, e);
	CPPUNIT_ASSERT(e.name == "3");
	CPPUNIT_ASSERT(e.chr == 22);
	CPPUNIT_ASSERT(e.position == 222);
    CPPUNIT_ASSERT(e.metrics[0].GetValueInt8() == 10);
    CPPUNIT_ASSERT(e.metrics[1].GetValueUInt8() == 11);
    CPPUNIT_ASSERT(e.metrics[2].GetValueInt16() == 18);
    CPPUNIT_ASSERT(e.metrics[3].GetValueUInt16() == 19);
    CPPUNIT_ASSERT(e.metrics[4].GetValueInt32() == 34);
    CPPUNIT_ASSERT(e.metrics[5].GetValueUInt32() == 35);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.metrics[6].GetValueFloat(), 66.0f, 0.00001f);
    CPPUNIT_ASSERT(e.metrics[7].GetValueAscii() == "text2");
    CPPUNIT_ASSERT(e.metrics[8].GetValueText() == L"ascii2");

    data.GetCopyNumberEntry(CopyNumberMultiDataType, 3, e);
	CPPUNIT_ASSERT(e.name == "4");
	CPPUNIT_ASSERT(e.chr == 4);
	CPPUNIT_ASSERT(e.position == 40);
    CPPUNIT_ASSERT(e.metrics[0].GetValueInt8() == 8);
    CPPUNIT_ASSERT(e.metrics[1].GetValueUInt8() == 8);
    CPPUNIT_ASSERT(e.metrics[2].GetValueInt16() == 16);
    CPPUNIT_ASSERT(e.metrics[3].GetValueUInt16() == 16);
    CPPUNIT_ASSERT(e.metrics[4].GetValueInt32() == 32);
    CPPUNIT_ASSERT(e.metrics[5].GetValueUInt32() == 32);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.metrics[6].GetValueFloat(), 44.0f, 0.00001f);
    CPPUNIT_ASSERT(e.metrics[7].GetValueAscii() == "ascii");
    CPPUNIT_ASSERT(e.metrics[8].GetValueText() == L"text");


    data.GetCytoEntry(CytoMultiDataType, 0, c);
	CPPUNIT_ASSERT(c.name == "1");
	CPPUNIT_ASSERT(c.call == 11);
	CPPUNIT_ASSERT(c.chr == 11);
    CPPUNIT_ASSERT(c.startPosition == 11);
    CPPUNIT_ASSERT(c.stopPosition == 22);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(c.confidenceScore, 111.0f, 0.00001f);

    data.GetCytoEntry(CytoMultiDataType, 1, c);
	CPPUNIT_ASSERT(c.name == "2");
	CPPUNIT_ASSERT(c.chr == 22);
    CPPUNIT_ASSERT(c.startPosition == 22);
    CPPUNIT_ASSERT(c.stopPosition == 33);
	CPPUNIT_ASSERT(c.call == 2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(c.confidenceScore, 20.0f, 0.00001f);

    data.GetCopyNumberVariationEntry(CopyNumberVariationMultiDataType, 0, v);
	CPPUNIT_ASSERT(v.name == "1");
	CPPUNIT_ASSERT(v.call == 11);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v.confidenceScore, 111.0f, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v.signal, 111.0f, 0.00001f);

    data.GetCopyNumberVariationEntry(CopyNumberVariationMultiDataType, 1, v);
	CPPUNIT_ASSERT(v.name == "2");
	CPPUNIT_ASSERT(v.call == 22);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v.confidenceScore, 222.0f, 0.00002f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v.signal, 222.0f, 0.00002f);

}
