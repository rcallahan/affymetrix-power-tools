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
#include "calvin_files/fusion/test/FusionCHPGenericDataTest.h"
//
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "calvin_files/fusion/src/FusionCHPGenericData.h"
//

using namespace std;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( FusionCHPGenericDataTest );

void FusionCHPGenericDataTest::setUp()
{
}

void FusionCHPGenericDataTest::tearDown()
{
}

void FusionCHPGenericDataTest::TestFileId()
{
	FusionCHPData *chp = FusionCHPDataReg::Read("../data/small_cel_file_partial_datheader");
	CPPUNIT_ASSERT(chp != NULL);
	FusionCHPGenericData *tileChp = FusionCHPGenericData::FromBase(chp); 
	CPPUNIT_ASSERT(tileChp != NULL);
	CPPUNIT_ASSERT(tileChp->FileId() == "0000065535-1147280233-0000005844-0000011008-0000006224");
	CPPUNIT_ASSERT(tileChp->GetGenericData()->FileIdentifier() == "0000065535-1147280233-0000005844-0000011008-0000006224");

}

void FusionCHPGenericDataTest::testRead()
{
	FusionCHPData *chp;
	FusionCHPGenericData *genchp;

	chp = FusionCHPDataReg::Read("../../parsers/data/no_file_exists");
	genchp = FusionCHPGenericData::FromBase(chp);
	CPPUNIT_ASSERT(chp == NULL);
	CPPUNIT_ASSERT(genchp == NULL);
	delete chp;

	chp = FusionCHPDataReg::Read("../../parsers/data/CHP_genotype_file");
	genchp = FusionCHPGenericData::FromBase(chp); 
	CPPUNIT_ASSERT(chp != NULL);
	CPPUNIT_ASSERT(genchp == NULL);
	delete chp;

	chp = FusionCHPDataReg::Read("../../parsers/data/CHP_expression_file");
	genchp = FusionCHPGenericData::FromBase(chp); 
	CPPUNIT_ASSERT(chp != NULL);
	CPPUNIT_ASSERT(genchp == NULL);
	delete chp;

	chp = FusionCHPDataReg::Read("../../parsers/data/CHP_tiling_file");
	genchp = FusionCHPGenericData::FromBase(chp); 
	CPPUNIT_ASSERT(chp != NULL);
	CPPUNIT_ASSERT(genchp == NULL);
	delete chp;

	chp = FusionCHPDataReg::Read("../../parsers/data/CHP_reseq_file");
	genchp = FusionCHPGenericData::FromBase(chp); 
	CPPUNIT_ASSERT(chp != NULL);
	CPPUNIT_ASSERT(genchp == NULL);
	delete chp;

	chp = FusionCHPDataReg::Read("../../parsers/data/CHP_universal_file");
	genchp = FusionCHPGenericData::FromBase(chp); 
	CPPUNIT_ASSERT(chp != NULL);
	CPPUNIT_ASSERT(genchp == NULL);
	delete chp;

	chp = FusionCHPDataReg::Read("../../../file/CPPTest/data/test.exp.abs.CHP");
	genchp = FusionCHPGenericData::FromBase(chp); 
	CPPUNIT_ASSERT(chp != NULL);
	CPPUNIT_ASSERT(genchp == NULL);
	delete chp;
	

}
