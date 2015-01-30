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
#include "calvin_files/fusion/test/FusionCHPQuantificationDetectionDataTest.h"
//
#include "calvin_files/fusion/src/FusionCHPQuantificationDetectionData.h"
//

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( FusionCHPQuantificationDetectionDataTest );

void FusionCHPQuantificationDetectionDataTest::setUp()
{
}

void FusionCHPQuantificationDetectionDataTest::tearDown()
{
}

void FusionCHPQuantificationDetectionDataTest::testCreation()
{
	CPPUNIT_ASSERT(1);
}

void FusionCHPQuantificationDetectionDataTest::TestFileId()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_quantification_detection_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPQuantificationDetectionData *sigChp = FusionCHPQuantificationDetectionData::FromBase(chp); 
	CPPUNIT_ASSERT(sigChp != NULL);
	CPPUNIT_ASSERT(sigChp->FileId() == "0000065535-1152158444-0000019912-0000001869-0000011538");
	CPPUNIT_ASSERT(sigChp->GetGenericData()->FileIdentifier() == "0000065535-1152158444-0000019912-0000001869-0000011538");
}

void FusionCHPQuantificationDetectionDataTest::testReadNonSignal()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_expression_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPQuantificationDetectionData *sigChp = FusionCHPQuantificationDetectionData::FromBase(chp); 
	CPPUNIT_ASSERT(sigChp == NULL);
	delete chp;
}

void FusionCHPQuantificationDetectionDataTest::testRead()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_quantification_detection_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPQuantificationDetectionData *sigChp = FusionCHPQuantificationDetectionData::FromBase(chp); 
	CPPUNIT_ASSERT(sigChp != NULL);

	CPPUNIT_ASSERT(sigChp->GetAlgName() == L"sig");
	CPPUNIT_ASSERT(sigChp->GetAlgVersion() == L"1.0");
	CPPUNIT_ASSERT(sigChp->GetArrayType() == L"test3");
	CPPUNIT_ASSERT(sigChp->GetEntryCount() == 2);
	
	ParameterNameValueTypeList params = sigChp->GetAlgParams();
	CPPUNIT_ASSERT(params.size() == 1);
	ParameterNameValueTypeList::iterator it=params.begin();
	ParameterNameValueType &param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"an1");
	CPPUNIT_ASSERT(param.GetValueText() == L"av1");

	params = sigChp->GetSummaryParams();
	CPPUNIT_ASSERT(params.size() == 1);
	it=params.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"sn1");
	CPPUNIT_ASSERT(param.GetValueText() == L"sv1");

	affymetrix_calvin_data::ProbeSetQuantificationDetectionData e;
	sigChp->GetQuantificationDetectionEntry(0, e);
	CPPUNIT_ASSERT(e.name == "abc");
    CPPUNIT_ASSERT(e.quantification == 10.0f);
    CPPUNIT_ASSERT(e.pvalue == 0.10f);
	sigChp->GetQuantificationDetectionEntry(1, e);
	CPPUNIT_ASSERT(e.name == "xyz");
	CPPUNIT_ASSERT(e.quantification == 20.0f);
    CPPUNIT_ASSERT(e.pvalue == 0.20f);

	
	delete sigChp;
}
