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
#include "calvin_files/fusion/test/FusionCHPMultiDataDataTest.h"
//
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
//

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_data;

CPPUNIT_TEST_SUITE_REGISTRATION( FusionCHPMultiDataDataTest );

void FusionCHPMultiDataDataTest::setUp()
{
}

void FusionCHPMultiDataDataTest::tearDown()
{
}

void FusionCHPMultiDataDataTest::testCreation()
{
	CPPUNIT_ASSERT(1);
}

void FusionCHPMultiDataDataTest::testReadNonGenotype()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_expression_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPMultiDataData *genoChp = FusionCHPMultiDataData::FromBase(chp); 
	CPPUNIT_ASSERT(genoChp == NULL);
	delete chp;
}

void FusionCHPMultiDataDataTest::testRead()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_MultiData_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPMultiDataData *genoChp = FusionCHPMultiDataData::FromBase(chp); 
	CPPUNIT_ASSERT(genoChp != NULL);


    affymetrix_calvin_data::ProbeSetMultiDataExpressionData ex;
    affymetrix_calvin_data::ProbeSetMultiDataGenotypeData gn;
	ParameterNameValueType param;

	CPPUNIT_ASSERT(genoChp->GetAlgName() == L"sig");
	CPPUNIT_ASSERT(genoChp->GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(genoChp->GetArrayType() == L"test3");
	CPPUNIT_ASSERT(genoChp->GetEntryCount(ExpressionMultiDataType) == 1);
	CPPUNIT_ASSERT(genoChp->GetEntryCount(GenotypeMultiDataType) == 2);
	CPPUNIT_ASSERT(genoChp->GetEntryCount(ExpressionControlMultiDataType) == 3);
	CPPUNIT_ASSERT(genoChp->GetEntryCount(GenotypeControlMultiDataType) == 1);

    // expression
    CPPUNIT_ASSERT(genoChp->GetNumMetricColumns(ExpressionMultiDataType) == 0);
    genoChp->GetExpressionEntry(ExpressionMultiDataType, 0, ex);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, 10.0f, 0.0001f);
	CPPUNIT_ASSERT(ex.name == "ex1");

    // genotype
    CPPUNIT_ASSERT(genoChp->GetNumMetricColumns(GenotypeMultiDataType) == 2);
    CPPUNIT_ASSERT(genoChp->GetMetricColumnName(GenotypeMultiDataType, 0) == L"int");
    CPPUNIT_ASSERT(genoChp->GetMetricColumnName(GenotypeMultiDataType, 1) == L"float");

    genoChp->GetGenotypeEntry(GenotypeMultiDataType, 0, gn);
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

    genoChp->GetGenotypeEntry(GenotypeMultiDataType, 1, gn);
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


    // genotype control
    CPPUNIT_ASSERT(genoChp->GetNumMetricColumns(GenotypeControlMultiDataType) == 0);
    genoChp->GetGenotypeEntry(GenotypeControlMultiDataType, 0, gn);
	CPPUNIT_ASSERT(gn.call == 2);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(gn.confidence, 22.0f, 0.0001f);
	CPPUNIT_ASSERT(gn.metrics.size() == 0);
	CPPUNIT_ASSERT(gn.name == "gc1");



    // expression control
    CPPUNIT_ASSERT(genoChp->GetNumMetricColumns(ExpressionControlMultiDataType) == 0);
	genoChp->GetExpressionEntry(ExpressionControlMultiDataType, 0, ex);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, 20.0f, 0.0001f);
	CPPUNIT_ASSERT(ex.name == "ec1");
	genoChp->GetExpressionEntry(ExpressionControlMultiDataType, 1, ex);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, 30.0f, 0.0001f);
	CPPUNIT_ASSERT(ex.name == "ec2");
	genoChp->GetExpressionEntry(ExpressionControlMultiDataType, 2, ex);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ex.quantification, 40.0f, 0.0001f);
	CPPUNIT_ASSERT(ex.name == "ec3");

	delete genoChp;
}

void FusionCHPMultiDataDataTest::testCN()
{
    FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_MultiData_file_cn");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPMultiDataData *cnChp = FusionCHPMultiDataData::FromBase(chp); 
	CPPUNIT_ASSERT(cnChp != NULL);

    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData e;
	ParameterNameValueType param;

	CPPUNIT_ASSERT(cnChp->GetAlgName() == L"sig");
	CPPUNIT_ASSERT(cnChp->GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(cnChp->GetArrayType() == L"test3");
	CPPUNIT_ASSERT(cnChp->GetEntryCount(CopyNumberMultiDataType) == 2);
	CPPUNIT_ASSERT(cnChp->GetEntryCount(GenotypeMultiDataType) == 0);
	CPPUNIT_ASSERT(cnChp->GetEntryCount(ExpressionMultiDataType) == 0);
	CPPUNIT_ASSERT(cnChp->GetEntryCount(ExpressionControlMultiDataType) == 0);
	CPPUNIT_ASSERT(cnChp->GetEntryCount(GenotypeControlMultiDataType) == 0);

	ParameterNameValueTypeList p = cnChp->GetAlgParams();
	ParameterNameValueTypeList::iterator it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	p = cnChp->GetSummaryParams();
	it = p.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");

    cnChp->GetCopyNumberEntry(CopyNumberMultiDataType, 0, e);
	CPPUNIT_ASSERT(e.chr == 10);
	CPPUNIT_ASSERT(e.position == 11);
	CPPUNIT_ASSERT(e.name == "abc");
	cnChp->GetCopyNumberEntry(CopyNumberMultiDataType, 1, e);
	CPPUNIT_ASSERT(e.chr == 20);
	CPPUNIT_ASSERT(e.position == 21);
	CPPUNIT_ASSERT(e.name == "xyz");

	delete cnChp;
}

void FusionCHPMultiDataDataTest::testReadCyto()
{
	return;

	const char *file = "\\\\affymetrix.com\\shares\\santaclara\\Transfer\\Walt\\HapMap-As_NA18526_C09_01_NN_20081218.Cytogenetics.cyto2.cychp";

	FusionCHPData *chp = FusionCHPDataReg::Read(file);
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPMultiDataData *cnChp = FusionCHPMultiDataData::FromBase(chp); 
	CPPUNIT_ASSERT(cnChp != NULL);

	int n;
	MarkerABSignals signals;
	AllelePeaks peaks;
	n = cnChp->GetEntryCount(AllelePeaksMultiDataType);
	for (int i=0; i<n; i++)
	{
		cnChp->GetAllelePeakEntry(AllelePeaksMultiDataType, i, peaks);
	}
	n = cnChp->GetEntryCount(MarkerABSignalsMultiDataType);
	for (int i=433; i<n; i++)
	{
		cnChp->GetMarkerABSignalsEntry(MarkerABSignalsMultiDataType, i, signals);
	}
	delete cnChp;
}
