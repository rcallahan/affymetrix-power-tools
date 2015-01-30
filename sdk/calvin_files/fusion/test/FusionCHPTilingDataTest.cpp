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
#include "calvin_files/fusion/test/FusionCHPTilingDataTest.h"
//
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "calvin_files/fusion/src/FusionCHPTilingData.h"
//

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( FusionCHPTilingDataTest );

void FusionCHPTilingDataTest::setUp()
{
}

void FusionCHPTilingDataTest::tearDown()
{
}

void FusionCHPTilingDataTest::testCreation()
{
	CPPUNIT_ASSERT(1);
}

void FusionCHPTilingDataTest::TestFileId()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_tiling_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPTilingData *tileChp = FusionCHPTilingData::FromBase(chp); 
	CPPUNIT_ASSERT(tileChp != NULL);
	CPPUNIT_ASSERT(tileChp->FileId() == "0000039321-1131034775-0000032391-0000005436-0000004827");
	CPPUNIT_ASSERT(tileChp->GetGenericData()->FileIdentifier() == "0000039321-1131034775-0000032391-0000005436-0000004827");
}

void FusionCHPTilingDataTest::testReadNonTiling()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_expression_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPTilingData *tileChp = FusionCHPTilingData::FromBase(chp); 
	CPPUNIT_ASSERT(tileChp == NULL);
	delete chp;
}

void FusionCHPTilingDataTest::testRead()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../../parsers/data/CHP_tiling_file");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPTilingData *tileChp = FusionCHPTilingData::FromBase(chp); 
	CPPUNIT_ASSERT(tileChp != NULL);

	CPPUNIT_ASSERT(tileChp->GetNumberSequences() == 2);
	CPPUNIT_ASSERT(tileChp->GetAlgName() == L"tile");
	CPPUNIT_ASSERT(tileChp->GetAlgVersion() == L"1.0");

	ParameterNameValueTypeList params = tileChp->GetAlgParams();
	CPPUNIT_ASSERT(params.size() == 1);
	ParameterNameValueTypeList::iterator it=params.begin();
	ParameterNameValueType &param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"v1");

	const double eps=0.00001;
	CHPTilingEntry e;
	TilingSequenceData seq;

	tileChp->OpenTilingSequenceDataSet(0);
	seq = tileChp->GetTilingSequenceData();
	
	CPPUNIT_ASSERT(seq.name == L"n1");
	CPPUNIT_ASSERT(seq.groupName == L"g1");
	CPPUNIT_ASSERT(seq.version == L"v1");
	CPPUNIT_ASSERT(seq.parameters.size() == 1);
	it = seq.parameters.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"seq1_p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"seq1_v1");

	CPPUNIT_ASSERT(tileChp->GetTilingSequenceEntryCount(0) == 2);
	tileChp->GetTilingSequenceEntry(0, e);
	CPPUNIT_ASSERT(e.position == 10);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 10.0f, eps);

	tileChp->GetTilingSequenceEntry(1, e);
	CPPUNIT_ASSERT(e.position == 20);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 20.0f, eps);


	tileChp->OpenTilingSequenceDataSet(1);
	seq = tileChp->GetTilingSequenceData();
	
	CPPUNIT_ASSERT(seq.name == L"n2");
	CPPUNIT_ASSERT(seq.groupName == L"g2");
	CPPUNIT_ASSERT(seq.version == L"v2");
	CPPUNIT_ASSERT(seq.parameters.size() == 1);
	it = seq.parameters.begin();
	param = *it;
	CPPUNIT_ASSERT(param.GetName() == L"seq2_p1");
	CPPUNIT_ASSERT(param.GetValueText() == L"seq2_v1");

	CPPUNIT_ASSERT(tileChp->GetTilingSequenceEntryCount(1) == 3);
	tileChp->GetTilingSequenceEntry(0, e);
	CPPUNIT_ASSERT(e.position == 11);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 11.0f, eps);

	tileChp->GetTilingSequenceEntry(1, e);
	CPPUNIT_ASSERT(e.position == 21);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 21.0f, eps);

	tileChp->GetTilingSequenceEntry(2, e);
	CPPUNIT_ASSERT(e.position == 31);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.value, 31.0f, eps);
	
	delete tileChp;
}
