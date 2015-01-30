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
#include "calvin_files/writers/test/CalvinCHPMultiDataFileWriterTest.h"
//
#include "calvin_files/data/src/ProbeSetMultiDataData.h"
#include "calvin_files/parsers/src/CHPMultiDataFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
//
#include <stdio.h>

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPMultiDataFileWriterTest );

void CHPMultiDataFileWriterTest::setUp() {}

void CHPMultiDataFileWriterTest::tearDown(){}

void CHPMultiDataFileWriterTest::testCreation()
{
	CHPMultiDataData fHdr("CHP_MultiData_file_empty");
	CHPMultiDataFileWriter* w = new CHPMultiDataFileWriter(fHdr);
	CPPUNIT_ASSERT(1);
	delete w;
}

void CHPMultiDataFileWriterTest::WriteTestGeno()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(GenotypeMultiDataType, 2, 10);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataGenotypeData e;

	writer->SeekToDataSet(GenotypeMultiDataType);
	e.name = "abc";
	e.confidence = 10.0f;
	writer->WriteEntry(e);
	e.name = "xyz";
	e.confidence = 20.0f;
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

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


    data2.GetGenotypeEntry(GenotypeMultiDataType, 0, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "abc");
	data2.GetGenotypeEntry(GenotypeMultiDataType, 1, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidence, 20.0f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "xyz");

}

void CHPMultiDataFileWriterTest::WriteTestCN()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file_cn");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(CopyNumberMultiDataType, 2, 10);
    data.SetEntryCount(CytoMultiDataType, 2, 10);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData e;
	affymetrix_calvin_data::ProbeSetMultiDataCytoRegionData c;

	writer->SeekToDataSet(CopyNumberMultiDataType);
	e.name = "abc";
    e.chr = 10;
    e.position = 11;
	writer->WriteEntry(e);
	e.name = "xyz";
    e.chr = 20;
    e.position = 21;
	writer->WriteEntry(e);

	writer->SeekToDataSet(CytoMultiDataType);
	c.name = "abc";
    c.call = 10;
    c.confidenceScore = 11.0f;
    c.chr = 10;
    c.startPosition = 10;
    c.stopPosition = 11;
	writer->WriteEntry(c);
	c.name = "xyz";
    c.call = 20;
    c.confidenceScore = 21.0f;
    c.chr = 20;
    c.startPosition = 20;
    c.stopPosition = 21;
	writer->WriteEntry(c);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file_cn");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(CopyNumberMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

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

  data2.GetCopyNumberEntry(CopyNumberMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.chr == 10);
	CPPUNIT_ASSERT(e.position == 11);
	CPPUNIT_ASSERT(e.name == "abc");
	data2.GetCopyNumberEntry(CopyNumberMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.chr == 20);
	CPPUNIT_ASSERT(e.position == 21);
	CPPUNIT_ASSERT(e.name == "xyz");

  data2.GetCytoEntry(CytoMultiDataType, 0, c);
	CPPUNIT_ASSERT(c.call == 10);
	CPPUNIT_ASSERT(c.chr == 10);
	CPPUNIT_ASSERT(c.startPosition == 10);
	CPPUNIT_ASSERT(c.stopPosition == 11);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(c.confidenceScore, 11.0f, 0.0001f);
	CPPUNIT_ASSERT(c.name == "abc");
	data2.GetCytoEntry(CytoMultiDataType, 1, c);
	CPPUNIT_ASSERT(c.call == 20);
	CPPUNIT_ASSERT(c.chr == 20);
	CPPUNIT_ASSERT(c.startPosition == 20);
	CPPUNIT_ASSERT(c.stopPosition == 21);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(c.confidenceScore, 21.0f, 0.0001f);
	CPPUNIT_ASSERT(c.name == "xyz");
}

void CHPMultiDataFileWriterTest::WriteTestDMET()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("DMET_CHP_MultiData");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(DmetBiAllelicMultiDataType, 2, 10);
	data.SetEntryCount(DmetMultiAllelicMultiDataType, 2, 10);
	data.SetEntryCount(DmetCopyNumberMultiDataType, 2, 10);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);

	//write the DMET genotype data set
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::DmetBiAllelicData e1;
	writer->SeekToDataSet(DmetBiAllelicMultiDataType);
	e1.name = "abc";
	e1.call = 10;
	e1.confidence = 10.0f;
	e1.force = 11;
	e1.signalA = 22.0f;
	e1.signalB = 34.0f;
	e1.contextB = 77;
	writer->WriteEntry(e1);
	e1.name = "xyz";
	e1.call = 11;
	e1.confidence = 44.0f;
	e1.force = 33;
	e1.signalA = 22.6f;
	e1.signalB = 34.2f;
	e1.contextA = 99;
	e1.contextB = 98;
	writer->WriteEntry(e1);

	//write the DMET multi-allelic data set
	affymetrix_calvin_data::DmetMultiAllelicData e2;
	writer->SeekToDataSet(DmetMultiAllelicMultiDataType);
	e2.name = "abc";
	e2.call = 10;
	e2.confidence = 10.0f;
	e2.force = 11;
	e2.alleleCount = 77;
	e2.signalA = 22.0f;
	e2.signalB = 34.0f;
	e2.signalC = 34.0f;
	e2.signalD = 34.0f;
	e2.signalE = 34.0f;
	e2.signalF = 34.0f;
	e2.contextA = 17;
	e2.contextB = 77;
	e2.contextC = 77;
	e2.contextD = 77;
	e2.contextE = 77;
	e2.contextF = 77;
	writer->WriteEntry(e2);
	e2.name = "xyz";
	e2.call = 11;
	e2.confidence = 44.0f;
	e2.force = 33;
	e2.alleleCount = 66;
	e2.signalA = 22.0f;
	e2.signalB = 34.0f;
	e2.signalC = 34.0f;
	e2.signalD = 34.0f;
	e2.signalE = 34.0f;
	e2.signalF = 34.0f;
	e2.contextA = 17;
	e2.contextB = 77;
	e2.contextC = 77;
	e2.contextD = 77;
	e2.contextE = 77;
	e2.contextF = 77;
	writer->WriteEntry(e2);

	//write the DMET copy number data set
	affymetrix_calvin_data::DmetCopyNumberData e3;
	writer->SeekToDataSet(DmetCopyNumberMultiDataType);
	e3.name = "abc";
	e3.call = 10;
	e3.confidence = 10.0f;
	e3.force = 11;
	e3.estimate = 22.0f;
	e3.lower = 34.0f;
	e3.upper = 34.0f;
	writer->WriteEntry(e3);
	e3.name = "xyz";
	e3.call = 11;
	e3.confidence = 44.0f;
	e3.force = 33;
	e3.estimate = 22.0f;
	e3.lower = 34.0f;
	e3.upper = 34.0f;
	writer->WriteEntry(e3);

	delete writer;

	CPPUNIT_ASSERT(1);

  /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("DMET_CHP_MultiData");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(DmetBiAllelicMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(DmetMultiAllelicMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(DmetCopyNumberMultiDataType) == 2);

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

  data2.GetEntry(DmetBiAllelicMultiDataType, 0, e1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e1.confidence, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(e1.name == "abc");
	data2.GetEntry(DmetBiAllelicMultiDataType, 1, e1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e1.confidence, 44.0f, 0.0001f);
	CPPUNIT_ASSERT(e1.name == "xyz");

	data2.GetEntry(DmetMultiAllelicMultiDataType, 0, e2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e2.confidence, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(e2.name == "abc");
	data2.GetEntry(DmetMultiAllelicMultiDataType, 1, e2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e2.confidence, 44.0f, 0.0001f);
	CPPUNIT_ASSERT(e2.name == "xyz");

	data2.GetEntry(DmetCopyNumberMultiDataType, 0, e3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e3.confidence, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(e3.name == "abc");
	data2.GetEntry(DmetCopyNumberMultiDataType, 1, e3);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e3.confidence, 44.0f, 0.0001f);
	CPPUNIT_ASSERT(e3.name == "xyz");
}

void CHPMultiDataFileWriterTest::WriteTestCNV()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file_cnv");

	data.SetAlgName(L"canary");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(CopyNumberVariationMultiDataType, 2, 10);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);

	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataCopyNumberVariationRegionData e;

	writer->SeekToDataSet(CopyNumberVariationMultiDataType);
	e.name = "abc";
    e.signal = 100;
    e.call = 1;
    e.confidenceScore = 1.1;
	writer->WriteEntry(e);
	e.name = "xyz";
    e.signal = 200;
    e.call = 2;
    e.confidenceScore = 2.2;
	writer->WriteEntry(e);
	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file_cnv");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"canary");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
    CPPUNIT_ASSERT(data2.GetEntryCount(CopyNumberVariationMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(CopyNumberMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

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

    data2.GetCopyNumberVariationEntry(CopyNumberVariationMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.call == 1);
	CPPUNIT_ASSERT(e.name == "abc");
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidenceScore, 1.1, 0.0001f);
	data2.GetCopyNumberVariationEntry(CopyNumberVariationMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.call == 2);
	CPPUNIT_ASSERT(e.name == "xyz");
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidenceScore, 2.2, 0.00001f);
   
}

void CHPMultiDataFileWriterTest::WriteTestExp()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(ExpressionMultiDataType, 2, 10);

	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataExpressionData e;

	writer->SeekToDataSet(ExpressionMultiDataType);
	e.name = "abc";
    e.quantification = 10.0f;
	writer->WriteEntry(e);
	e.name = "xyz";
	e.quantification = 20.0f;
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);

    data2.GetExpressionEntry(ExpressionMultiDataType, 0, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "abc");
	data2.GetExpressionEntry(ExpressionMultiDataType, 1, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.quantification, 20.0f, 0.0001f);
	CPPUNIT_ASSERT(e.name == "xyz");

}

void CHPMultiDataFileWriterTest::WriteTestAll()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_file");

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(ExpressionMultiDataType, 1, 10);
	data.SetEntryCount(GenotypeMultiDataType, 2, 10, cols);

	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::ProbeSetMultiDataExpressionData ex;
    affymetrix_calvin_data::ProbeSetMultiDataGenotypeData gn;

    // expression
	writer->SeekToDataSet(ExpressionMultiDataType);
	ex.name = "ex1";
    ex.quantification = 10.0f;
	writer->WriteEntry(ex);

    // genotype
    writer->SeekToDataSet(GenotypeMultiDataType);
	gn.name = "gn1";
	gn.call = 1;
    gn.confidence = 11.0f;
    gn.metrics.clear();
    param.SetValueInt32(1);
    gn.metrics.push_back(param);
    param.SetValueFloat(2.0f);
    gn.metrics.push_back(param);
	writer->WriteEntry(gn);

	gn.name = "gn2";
	gn.call = 2;
    gn.confidence = 22.0f;;
    gn.metrics.clear();
    param.SetValueInt32(2);
    gn.metrics.push_back(param);
    param.SetValueFloat(3.0f);
    gn.metrics.push_back(param);
	writer->WriteEntry(gn);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 1);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 2);

    // expression
    CPPUNIT_ASSERT(data2.GetNumMetricColumns(ExpressionMultiDataType) == 0);
    data2.GetExpressionEntry(ExpressionMultiDataType, 0, ex);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(ex.name == "ex1");

    // genotype
    CPPUNIT_ASSERT(data2.GetNumMetricColumns(GenotypeMultiDataType) == 2);
    CPPUNIT_ASSERT(data2.GetMetricColumnName(GenotypeMultiDataType, 0) == L"int");
    CPPUNIT_ASSERT(data2.GetMetricColumnName(GenotypeMultiDataType, 1) == L"float");

    data2.GetGenotypeEntry(GenotypeMultiDataType, 0, gn);
	CPPUNIT_ASSERT(gn.call == 1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, 11.0f, 0.0001f);
	CPPUNIT_ASSERT(gn.name == "gn1");
    CPPUNIT_ASSERT(gn.metrics.size() == 2);
    param = gn.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 1);
    param = gn.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 2.0f, 0.00001f);

    data2.GetGenotypeEntry(GenotypeMultiDataType, 1, gn);
	CPPUNIT_ASSERT(gn.call == 2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, 22.0f, 0.0001f);
	CPPUNIT_ASSERT(gn.name == "gn2");
    CPPUNIT_ASSERT(gn.metrics.size() == 2);
    param = gn.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 2);
    param = gn.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 3.0f, 0.00001f);
}

void CHPMultiDataFileWriterTest::WriteTestChrSummary()
{
    affymetrix_calvin_data::ChromosomeMultiDataSummaryData p;

	CHPMultiDataData data("CHP_MultiData_chr_sum_file");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(ChromosomeSummaryMultiDataType, 3, 5);
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
    writer->SeekToDataSet(ChromosomeSummaryMultiDataType);

    float fval = 1.0f;
    p.chr=1;
	p.display = "1";
	p.startIndex=(u_int32_t)fval++;
	p.markerCount=(u_int32_t)fval++;
    p.minSignal=fval++;
    p.maxSignal=fval++;
    p.medianCnState=fval++;
    p.homFrequency=fval++;
    p.hetFrequency=fval++;
    //p.mosaicism=fval++;
    writer->WriteEntry(p);
    p.chr=255;
	p.display="255";
	p.startIndex=(u_int32_t)fval++;
	p.markerCount=(u_int32_t)fval++;
    p.minSignal=fval++;
    p.maxSignal=fval++;
    p.medianCnState=fval++;
    p.homFrequency=fval++;
    p.hetFrequency=fval++;
    //p.mosaicism=fval++;
    writer->WriteEntry(p);
    p.chr=3;
	p.display="3";
	p.startIndex=(u_int32_t)fval++;
	p.markerCount=(u_int32_t)fval++;
    p.minSignal=fval++;
    p.maxSignal=fval++;
    p.medianCnState=fval++;
    p.homFrequency=fval++;
    p.hetFrequency=fval++;
    //p.mosaicism=fval++;
    writer->WriteEntry(p);
    delete writer;

    fval = 1.0f;
	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_chr_sum_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(ChromosomeSummaryMultiDataType) == 3);

    data2.GetChromosomeSummaryEntry(ChromosomeSummaryMultiDataType, 0, p);
   	CPPUNIT_ASSERT(p.chr == 1);
   	CPPUNIT_ASSERT(p.display == "1");
	CPPUNIT_ASSERT(p.startIndex == (u_int32_t)fval++);
	CPPUNIT_ASSERT(p.markerCount == (u_int32_t)fval++);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.minSignal, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.maxSignal, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.medianCnState, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.homFrequency, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.hetFrequency, fval++, 0.0001f);
	//CPPUNIT_ASSERT_DOUBLES_EQUAL(p.mosaicism, fval++, 0.0001f);

    data2.GetChromosomeSummaryEntry(ChromosomeSummaryMultiDataType, 1, p);
   	CPPUNIT_ASSERT(p.chr == 255);
   	CPPUNIT_ASSERT(p.display == "255");
	CPPUNIT_ASSERT(p.startIndex == (u_int32_t)fval++);
	CPPUNIT_ASSERT(p.markerCount == (u_int32_t)fval++);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.minSignal, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.maxSignal, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.medianCnState, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.homFrequency, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.hetFrequency, fval++, 0.0001f);
	//CPPUNIT_ASSERT_DOUBLES_EQUAL(p.mosaicism, fval++, 0.0001f);

    data2.GetChromosomeSummaryEntry(ChromosomeSummaryMultiDataType, 2, p);
   	CPPUNIT_ASSERT(p.chr == 3);
   	CPPUNIT_ASSERT(p.display == "3");
	CPPUNIT_ASSERT(p.startIndex == (u_int32_t)fval++);
	CPPUNIT_ASSERT(p.markerCount == (u_int32_t)fval++);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.minSignal, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.maxSignal, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.medianCnState, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.homFrequency, fval++, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p.hetFrequency, fval++, 0.0001f);
	//CPPUNIT_ASSERT_DOUBLES_EQUAL(p.mosaicism, fval++, 0.0001f);

}

void CHPMultiDataFileWriterTest::WriteTestChrSegment()
{
    affymetrix_calvin_data::ChromosomeSegmentData p;

	CHPMultiDataData data("CHP_MultiData_chr_seg_file");

	ParameterNameValueType param;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");

    const MultiDataType dataTypes[] = {
        SegmentCNMultiDataType,
        SegmentLOHMultiDataType,
        SegmentCNNeutralLOHMultiDataType,
        SegmentNormalDiploidMultiDataType,
        SegmentMosaicismMultiDataType,
        SegmentNoCallMultiDataType
    };
    int nTypes = sizeof(dataTypes) / sizeof(MultiDataType);

    for (int itype=0; itype<nTypes; itype++)
    {
        if (dataTypes[itype] == SegmentNormalDiploidMultiDataType)
            data.SetEntryCount(SegmentNormalDiploidMultiDataType, 3+itype, 10, cols);
        else
            data.SetEntryCount(dataTypes[itype], 3+itype, 10);
    }
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);

    int ival=0;
    //char buf[10];
    for (int itype=0; itype<nTypes; itype++)
    {
    	writer->SeekToDataSet(dataTypes[itype]);
        for (int i=0; i<3+itype; i++)
        {
            p.chr = ival % 200;
            ++ival;
            p.startPosition = ival++;
            p.stopPosition = ival++;
            p.markerCount = ival++;
			p.meanMarkerDistance = ival++;
            p.segmentId = ival++;
            p.metrics.clear();
	        if (dataTypes[itype] == SegmentNormalDiploidMultiDataType)
            {
                param.SetValueInt32(ival++);
                p.metrics.push_back(param);
                param.SetValueFloat(ival++);
                p.metrics.push_back(param);
            }
            writer->WriteEntry(p);
            p.metrics.clear();
        }
    }
    delete writer;


	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_chr_seg_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");

    ival = 0;
    for (int itype=0; itype<nTypes; itype++)
    {
    	CPPUNIT_ASSERT(data2.GetEntryCount(dataTypes[itype]) == 3+itype);
    
        if (dataTypes[itype] == SegmentNormalDiploidMultiDataType)
        {
            CPPUNIT_ASSERT(data2.GetMetricColumnName(SegmentNormalDiploidMultiDataType, 0) == L"int");
            CPPUNIT_ASSERT(data2.GetMetricColumnName(SegmentNormalDiploidMultiDataType, 1) == L"float");
        }
        for (int i=0; i<3+itype; i++)
        {
            data2.GetChromosomeSegmentEntry(dataTypes[itype], i, p);
   	        CPPUNIT_ASSERT(p.chr == ival%200);
            ++ival;
	        CPPUNIT_ASSERT(p.startPosition == ival++);
            CPPUNIT_ASSERT(p.stopPosition == ival++);
            CPPUNIT_ASSERT(p.markerCount == ival++);
			CPPUNIT_ASSERT(p.meanMarkerDistance == ival++);
            CPPUNIT_ASSERT(p.segmentId == ival++);
            if (dataTypes[itype] == SegmentNormalDiploidMultiDataType)
            {
                param = p.metrics[0];
                CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
                CPPUNIT_ASSERT(param.GetValueInt32() == ival++);
                param = p.metrics[1];
                CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), ival++, 0.00001f);
            }
        }
    }

}

void CHPMultiDataFileWriterTest::WriteTestChrSegmentEx()
{
    affymetrix_calvin_data::ChromosomeSegmentDataEx p;

	CHPMultiDataData data("CHP_MultiData_chr_seg_file_ex");

	ParameterNameValueType param;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");

    const MultiDataType dataTypes[] = {
        /*SegmentCNMultiDataType,
        SegmentLOHMultiDataType,
        SegmentCNNeutralLOHMultiDataType,
        SegmentNormalDiploidMultiDataType,
        SegmentMosaicismMultiDataType,
        SegmentNoCallMultiDataType,*/
        SegmentGenotypeConcordanceMultiDataType,
        SegmentGenotypeDiscordanceMultiDataType,
        SegmentCNLossLOHConcordanceMultiDataType,
        SegmentCNNeutralLOHConcordanceMultiDataType,
        SegmentHeteroUPDMultiDataType,
        SegmentIsoUPDMultiDataType,
        SegmentDenovoCopyNumberMultiDataType,
        SegmentHemizygousParentOfOriginMultiDataType
    };
    int nTypes = sizeof(dataTypes) / sizeof(MultiDataType);

    for (int itype=0; itype<nTypes; itype++)
    {
        if (dataTypes[itype] == SegmentMosaicismMultiDataType)
            data.SetEntryCount(SegmentMosaicismMultiDataType, 3+itype, 10, cols);
        else
            data.SetEntryCount(dataTypes[itype], 3+itype, 10);
    }
	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);

    int ival=0;
    //char buf[10];
    for (int itype=0; itype<nTypes; itype++)
    {
    	writer->SeekToDataSet(dataTypes[itype]);
        for (int i=0; i<3+itype; i++)
        {
			p.referenceSampleKey = ival++;
			p.familialSampleKey = ival++;
            p.chr = ival % 200;
            ++ival;
            p.call = ival % 200;
            ++ival;
            p.confidence = ival++;
            p.startPosition = ival++;
            p.stopPosition = ival++;
            p.markerCount = ival++;
			p.homozygosity = ival++;
			p.heterozygosity = ival++;
            p.segmentId = ival++;
            p.metrics.clear();
	        if (dataTypes[itype] == SegmentMosaicismMultiDataType)
            {
                param.SetValueInt32(ival++);
                p.metrics.push_back(param);
                param.SetValueFloat(ival++);
                p.metrics.push_back(param);
            }
            writer->WriteEntry(p);
            p.metrics.clear();
        }
    }
    delete writer;


	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_chr_seg_file_ex");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");

    ival = 0;
    for (int itype=0; itype<nTypes; itype++)
    {
    	CPPUNIT_ASSERT(data2.GetEntryCount(dataTypes[itype]) == 3+itype);
    
        if (dataTypes[itype] == SegmentMosaicismMultiDataType)
        {
            CPPUNIT_ASSERT(data2.GetMetricColumnName(SegmentMosaicismMultiDataType, 0) == L"int");
            CPPUNIT_ASSERT(data2.GetMetricColumnName(SegmentMosaicismMultiDataType, 1) == L"float");
        }
        for (int i=0; i<3+itype; i++)
        {
            data2.GetChromosomeSegmentEntry(dataTypes[itype], i, p);
			CPPUNIT_ASSERT(p.referenceSampleKey == ival++);
			CPPUNIT_ASSERT(p.familialSampleKey == ival++);
   	        CPPUNIT_ASSERT(p.chr == ival%200);
            ++ival;
	        CPPUNIT_ASSERT(p.call == ival%200);
            ++ival;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(p.confidence, ival++, 0.0001f);
	        CPPUNIT_ASSERT(p.startPosition == ival++);
            CPPUNIT_ASSERT(p.stopPosition == ival++);
            CPPUNIT_ASSERT(p.markerCount == ival++);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(p.homozygosity, ival++, 0.0001f);
			CPPUNIT_ASSERT_DOUBLES_EQUAL(p.heterozygosity, ival++, 0.0001f);
            CPPUNIT_ASSERT(p.segmentId == ival++);
            if (dataTypes[itype] == SegmentMosaicismMultiDataType)
            {
                param = p.metrics[0];
                CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
                CPPUNIT_ASSERT(param.GetValueInt32() == ival++);
                param = p.metrics[1];
                CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), ival++, 0.00001f);
            }
        }
    }

}

void CHPMultiDataFileWriterTest::WriteTestFamilial()
{
    affymetrix_calvin_data::FamilialSample s;
    affymetrix_calvin_data::FamilialSegmentOverlap o;

	CHPMultiDataData data("CHP_MultiData_familial_file");

	ParameterNameValueType param;

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");

    data.SetEntryCount(FamilialSegmentOverlapsMultiDataType, 3, 10, 10, 10);
    data.SetEntryCount(FamilialSamplesMultiDataType, 3, 10, 10, 10, 10);

    CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);

    int ival=0;
    char buf[10];

	writer->SeekToDataSet(FamilialSamplesMultiDataType);
    for (int i=0; i<3; i++)
    {
        s.sampleKey = ival++;
        sprintf(buf,"%i",ival++);
	    s.arrID = buf;
        sprintf(buf,"%i",ival++);
	    s.chpID = buf;
        //wsprintf(wbuf,L"%i",ival++);
	    s.chpFilename = L"wbuf";
        sprintf(buf,"%i",ival++);
	    s.role = buf;
        s.roleValidity = (ival%2);
        ++ival;
	    s.roleConfidence = ival++;
        writer->WriteEntry(s);
    }
	writer->SeekToDataSet(FamilialSegmentOverlapsMultiDataType);
    for (int i=0; i<3; i++)
    {
        sprintf(buf,"%i",ival++);
        o.segmentType = buf;
        o.referenceSampleKey = ival++;
        sprintf(buf,"%i",ival++);
        o.referenceSegmentID = buf;
        o.familialSampleKey = ival++;
        sprintf(buf,"%i",ival++);
        o.familialSegmentID = buf;
        writer->WriteEntry(o);
    }
    delete writer;


	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_familial_file");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");

	CPPUNIT_ASSERT(data2.GetEntryCount(FamilialSamplesMultiDataType) == 3);
	CPPUNIT_ASSERT(data2.GetEntryCount(FamilialSegmentOverlapsMultiDataType) == 3);
    
    ival = 0;
    for (int i=0; i<3; i++)
    {
        data2.GetFamilialSampleEntry(FamilialSamplesMultiDataType, i, s);
        CPPUNIT_ASSERT(s.sampleKey == ival++);
        sprintf(buf,"%i",ival++);
	    CPPUNIT_ASSERT(s.arrID == buf);
        sprintf(buf,"%i",ival++);
	    CPPUNIT_ASSERT(s.chpID == buf);
        //wsprintf(wbuf,L"%i",ival++);
	    CPPUNIT_ASSERT(s.chpFilename == L"wbuf");
        sprintf(buf,"%i",ival++);
	    CPPUNIT_ASSERT(s.role == buf);
        CPPUNIT_ASSERT(s.roleValidity == (ival%2));
        ++ival;
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(s.roleConfidence, ival++, 0.0001f);
    }
    for (int i=0; i<3; i++)
    {
        data2.GetFamilialSegmentOverlapEntry(FamilialSegmentOverlapsMultiDataType, i, o);
        sprintf(buf,"%i",ival++);
        CPPUNIT_ASSERT(o.segmentType == buf);
        CPPUNIT_ASSERT(o.referenceSampleKey == ival++);
        sprintf(buf,"%i",ival++);
        CPPUNIT_ASSERT(o.referenceSegmentID == buf);
        CPPUNIT_ASSERT(o.familialSampleKey == ival++);
        sprintf(buf,"%i",ival++);
        CPPUNIT_ASSERT(o.familialSegmentID == buf);
    }

}

void CHPMultiDataFileWriterTest::WriteTestAllelePeaks()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;
    affymetrix_calvin_data::AllelePeaks s;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	CHPMultiDataData data("CHP_MultiData_allele_peaks");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(AllelePeaksMultiDataType, 2, 10, cols);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::AllelePeaks e;

	writer->SeekToDataSet(AllelePeaksMultiDataType);
	e.name = "abc";
    e.chr = 10;
    e.position = 11;
	e.peaks.clear();
    param.SetValueInt32(1);
    e.peaks.push_back(param);
    param.SetValueFloat(2.0f);
    e.peaks.push_back(param);
	writer->WriteEntry(e);

	e.name = "xyz";
    e.chr = 20;
    e.position = 21;
	e.peaks.clear();
    param.SetValueInt32(2);
    e.peaks.push_back(param);
    param.SetValueFloat(3.0f);
    e.peaks.push_back(param);
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_allele_peaks");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetNumMetricColumns(AllelePeaksMultiDataType) == 2);
    CPPUNIT_ASSERT(data2.GetMetricColumnName(AllelePeaksMultiDataType, 0) == L"int");
    CPPUNIT_ASSERT(data2.GetMetricColumnName(AllelePeaksMultiDataType, 1) == L"float");

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(AllelePeaksMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

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

	data2.GetEntry(AllelePeaksMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.chr == 10);
	CPPUNIT_ASSERT(e.position == 11);
	CPPUNIT_ASSERT(e.name == "abc");
	CPPUNIT_ASSERT(e.peaks.size() == 2);
    param = e.peaks[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 1);
    param = e.peaks[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 2.0f, 0.00001f);


	data2.GetEntry(AllelePeaksMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.chr == 20);
	CPPUNIT_ASSERT(e.position == 21);
	CPPUNIT_ASSERT(e.name == "xyz");
	CPPUNIT_ASSERT(e.peaks.size() == 2);
    param = e.peaks[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 2);
    param = e.peaks[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 3.0f, 0.00001f);

}

void CHPMultiDataFileWriterTest::WriteTestMarkerABSignals()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;
  //affymetrix_calvin_data::MarkerABSignals s;

    vector<ColumnInfo> cols;
    IntColumn icol(L"int");
    cols.push_back(icol);
    FloatColumn fcol(L"float");
    cols.push_back(fcol);

	CHPMultiDataData data("CHP_MultiData_marker_signals");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(MarkerABSignalsMultiDataType, 2, 10, cols);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::MarkerABSignals e;

	writer->SeekToDataSet(MarkerABSignalsMultiDataType);
	e.index = 1;
	e.metrics.clear();
    param.SetValueInt32(1);
    e.metrics.push_back(param);
    param.SetValueFloat(2.0f);
    e.metrics.push_back(param);
	writer->WriteEntry(e);

	e.index = 2;
	e.metrics.clear();
    param.SetValueInt32(2);
    e.metrics.push_back(param);
    param.SetValueFloat(3.0f);
    e.metrics.push_back(param);
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);

    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_marker_signals");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetNumMetricColumns(MarkerABSignalsMultiDataType) == 2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(MarkerABSignalsMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

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

	data2.GetEntry(MarkerABSignalsMultiDataType, 0, e);	
	CPPUNIT_ASSERT(e.index == 1);
	CPPUNIT_ASSERT(e.metrics.size() == 2);
    param = e.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 1);
    param = e.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 2.0f, 0.00001f);

	data2.GetEntry(MarkerABSignalsMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.index == 2);
	CPPUNIT_ASSERT(e.metrics.size() == 2);
    param = e.metrics[0];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::Int32Type);
    CPPUNIT_ASSERT(param.GetValueInt32() == 2);
    param = e.metrics[1];
    CPPUNIT_ASSERT(param.GetParameterType() == ParameterNameValueType::FloatType);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.GetValueFloat(), 3.0f, 0.00001f);
}


void CHPMultiDataFileWriterTest::WriteTestCytoGeno()
{
	ParameterNameValueTypeList params;
	ParameterNameValueType param;

	CHPMultiDataData data("CHP_MultiData_cytogenotype_data");

	data.SetAlgName(L"sig");
	data.SetAlgVersion(L"1.0");
	data.SetArrayType(L"test3");
	data.SetEntryCount(CytoGenotypeCallMultiDataType, 2, 10);

	param.SetName(L"an1");
	param.SetValueText(L"av1");
	params.push_back(param);
	data.AddAlgParams(params);

	params.clear();
	param.SetName(L"sn1");
	param.SetValueText(L"sv1");
	params.push_back(param);
	data.AddSummaryParams(params);


	CHPMultiDataFileWriter *writer = new CHPMultiDataFileWriter(data);
	affymetrix_calvin_data::CytoGenotypeCallData e;

	writer->SeekToDataSet(CytoGenotypeCallMultiDataType);
	e.index = 0;
	e.call = 1;
	e.aSignal = 1.0f;
	e.bSignal = 1.0f;
	e.confidence = 1.0f;
	e.signalStrength = 1.0f;
	e.forcedCall = 1;
	e.contrast = 1.0f;
	writer->WriteEntry(e);

	e.index = 1;
	e.call = 2;
	e.aSignal = 2.0f;
	e.bSignal = 2.0f;
	e.confidence = 2.0f;
	e.signalStrength = 2.0f;
	e.forcedCall = 2;
	e.contrast = 2.0f;
	writer->WriteEntry(e);

	delete writer;

	CPPUNIT_ASSERT(1);


    /////////////////////////////

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("CHP_MultiData_cytogenotype_data");
	reader.Read(data2);

	CPPUNIT_ASSERT(data2.GetAlgName() == L"sig");
	CPPUNIT_ASSERT(data2.GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(data2.GetArrayType() == L"test3");
	CPPUNIT_ASSERT(data2.GetEntryCount(CytoGenotypeCallMultiDataType) == 2);
	CPPUNIT_ASSERT(data2.GetEntryCount(ExpressionMultiDataType) == 0);

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


    data2.GetEntry(CytoGenotypeCallMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.index == 0);
	CPPUNIT_ASSERT(e.call == 1);
	CPPUNIT_ASSERT(e.confidence == 1);
	CPPUNIT_ASSERT(e.forcedCall == 1);
	CPPUNIT_ASSERT(e.aSignal == 1.0f);
	CPPUNIT_ASSERT(e.bSignal == 1.0f);
	CPPUNIT_ASSERT(e.signalStrength == 1.0f);
	CPPUNIT_ASSERT(e.contrast == 1.0f);
	data2.GetEntry(CytoGenotypeCallMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.index == 1);
	CPPUNIT_ASSERT(e.call == 2);
	CPPUNIT_ASSERT(e.confidence == 2);
	CPPUNIT_ASSERT(e.forcedCall == 2);
	CPPUNIT_ASSERT(e.aSignal == 2.0f);
	CPPUNIT_ASSERT(e.bSignal == 2.0f);
	CPPUNIT_ASSERT(e.signalStrength == 2.0f);
	CPPUNIT_ASSERT(e.contrast == 2.0f);

}


