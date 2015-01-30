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
#include "calvin_files/parsers/test/CHPMultiDataFileReaderTest.h"
//
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/parsers/src/CHPMultiDataFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
//
#include <cmath>
#include <cstring>
#include <string>
//

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_exceptions;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_parameter;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPMultiDataFileReaderTest );

void CHPMultiDataFileReaderTest::setUp()
{
}

void CHPMultiDataFileReaderTest::tearDown()
{
}

void CHPMultiDataFileReaderTest::testCreation()
{
	CHPMultiDataFileReader reader;
	CPPUNIT_ASSERT(1);
}

void CHPMultiDataFileReaderTest::testReadCN()
{
   	affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData e;
	ParameterNameValueType param;
	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("../data/CHP_MultiData_file_cn");
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
}

void CHPMultiDataFileReaderTest::testReadCNV()
{
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberVariationRegionData e;
	ParameterNameValueType param;

    CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("../data/CHP_MultiData_file_cnv");
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
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidenceScore, 1.1f, 0.0001f);
    CPPUNIT_ASSERT(e.signal == 100);
	data2.GetCopyNumberVariationEntry(CopyNumberVariationMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.call == 2);
	CPPUNIT_ASSERT(e.name == "xyz");
    CPPUNIT_ASSERT_DOUBLES_EQUAL(e.confidenceScore, 2.2, 0.0001f);
    CPPUNIT_ASSERT(e.signal == 200);
}

void CHPMultiDataFileReaderTest::testRead()
{
    affymetrix_calvin_data::ProbeSetMultiDataExpressionData ex;
    affymetrix_calvin_data::ProbeSetMultiDataGenotypeData gn;
	ParameterNameValueType param;

	CHPMultiDataData data2;
	CHPMultiDataFileReader reader;
	reader.SetFilename("../data/CHP_MultiData_file");
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

void CHPMultiDataFileReaderTest::testReadCyto()
{
	/*
	const char *file = "\\\\affymetrix.com\\shares\\santaclara\\Transfer\\Walt\\HapMap-As_NA18526_C09_01_NN_20081218.Cytogenetics.cyto2.cychp";

	CHPMultiDataData data;
	CHPMultiDataFileReader reader;
	reader.SetFilename(file);
	reader.Read(data);

	int n;
	MarkerABSignals signals;
	AllelePeaks peaks;

	n = data.GetEntryCount(AllelePeaksMultiDataType);
	for (int i=0; i<n; i++)
	{
		data.GetEntry(AllelePeaksMultiDataType, i, peaks);
	}
	n = data.GetEntryCount(MarkerABSignalsMultiDataType);
	for (int i=0; i<n; i++)
	{
		data.GetEntry(MarkerABSignalsMultiDataType, i, signals);
	}
	*/
}
