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
#include "calvin_files/fusion/test/FusionCHPQuantificationDataTest.h"
//
#include "calvin_files/fusion/src/FusionCHPQuantificationData.h"
//

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( FusionCHPQuantificationDataTest );

void FusionCHPQuantificationDataTest::setUp()
{
}

void FusionCHPQuantificationDataTest::tearDown()
{
}

void FusionCHPQuantificationDataTest::testCreation()
{
	CPPUNIT_ASSERT(1);
}

void FusionCHPQuantificationDataTest::TestFileId()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_quantification_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPQuantificationData *sigChp = FusionCHPQuantificationData::FromBase(chp); 
	CPPUNIT_ASSERT(sigChp != NULL);
	CPPUNIT_ASSERT(sigChp->FileId() == "0000065535-1152133546-0000000153-0000003902-0000014604");
	CPPUNIT_ASSERT(sigChp->GetGenericData()->FileIdentifier() == "0000065535-1152133546-0000000153-0000003902-0000014604");
}

void FusionCHPQuantificationDataTest::testReadNonSignal()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_expression_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPQuantificationData *sigChp = FusionCHPQuantificationData::FromBase(chp); 
	CPPUNIT_ASSERT(sigChp == NULL);
	delete chp;
}

void FusionCHPQuantificationDataTest::testRead()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_quantification_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPQuantificationData *sigChp = FusionCHPQuantificationData::FromBase(chp); 
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

	affymetrix_calvin_data::ProbeSetQuantificationData e;
	sigChp->GetQuantificationEntry(0, e);
	CPPUNIT_ASSERT(e.name == "abc");
    CPPUNIT_ASSERT(e.quantification == 10.0f);
	sigChp->GetQuantificationEntry(1, e);
	CPPUNIT_ASSERT(e.name == "xyz");
	CPPUNIT_ASSERT(e.quantification == 20.0f);

	
	delete sigChp;
}
