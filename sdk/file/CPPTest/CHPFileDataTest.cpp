////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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


#include "file/CPPTest/CHPFileDataTest.h"
//
#include "file/CHPFileData.h"
//
#include <cmath>
//

CPPUNIT_TEST_SUITE_REGISTRATION( CCHPFileDataTest );


#define EXP_COMP_CHP_FILE "./data/test.exp.comp.CHP"
#define EXP_ABS_CHP_FILE "./data/test.exp.abs.CHP"
#define NO_FILE "./data/no_file.CHP"
#define TAG_V11_FILE "./data/tag_v11.CHP"
#define TAG_XDA_FILE "./data/tag_XDA.CHP"
#define RESEQ_V1_FILE "./data/reseq_xda_v1.chp"
#define RESEQ_V2_FILE "./data/reseq_xda_v2.chp"
#define RESEQ_OLD_FILE "./data/reseq_v13.chp"

using namespace affxchp;

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

void CCHPFileDataTest::setUp()
{
}

void CCHPFileDataTest::tearDown()
{
}

void CCHPFileDataTest::testproperty_Header_Dimensions()
{
	affxchp::CCHPFileHeader header;
	header.SetRows(1);
	header.SetCols(2);
	CPPUNIT_ASSERT( header.GetRows() == 1);
	CPPUNIT_ASSERT( header.GetCols() == 2);
}

void CCHPFileDataTest::testproperty_Header_NumProbeSets()
{
	affxchp::CCHPFileHeader header;
	header.SetNumProbeSets(10);
	CPPUNIT_ASSERT( header.GetNumProbeSets() == 10);
}

void CCHPFileDataTest::testproperty_Header_AssayType()
{
	affxchp::CCHPFileHeader header;
	header.SetAssayType(affxchp::CCHPFileHeader::Expression);
	CPPUNIT_ASSERT( header.GetAssayType() == affxchp::CCHPFileHeader::Expression);
}

void CCHPFileDataTest::testproperty_Header_AlgName()
{
	affxchp::CCHPFileHeader header;
	header.SetAlgName("alg");
	CPPUNIT_ASSERT( header.GetAlgName() == "alg");
}

void CCHPFileDataTest::testproperty_Header_AlgVersion()
{
	affxchp::CCHPFileHeader header;
	header.SetAlgVersion("1.0");
	CPPUNIT_ASSERT( header.GetAlgVersion() == "1.0");
}

void CCHPFileDataTest::testproperty_Header_AlgParams()
{
	affxchp::CCHPFileHeader header;

	TagValuePairTypeList& params = header.AlgorithmParameters();
	TagValuePairType pair;
	pair.Tag = "tag1";
	pair.Value = "value1";
	params.push_back(pair);
	pair.Tag = "tag2";
	pair.Value = "value2";
	params.push_back(pair);
	CPPUNIT_ASSERT(header.GetAlgorithmParameter("tag1") == "value1");
	CPPUNIT_ASSERT(header.GetAlgorithmParameter("tag2") == "value2");
	int count;
	TagValuePairTypeList::iterator iter;
	for (count = 0, iter = params.begin(); iter != params.end(); iter++, count++)
	{
		switch (count)
		{
		case 0:
		CPPUNIT_ASSERT( iter->Tag == "tag1");
		CPPUNIT_ASSERT( iter->Value == "value1");
		break;
		case 1:
		CPPUNIT_ASSERT( iter->Tag == "tag2");
		CPPUNIT_ASSERT( iter->Value == "value2");
		break;
		}
	}
}

void CCHPFileDataTest::testproperty_Header_SummaryParams()
{
	affxchp::CCHPFileHeader header;

	TagValuePairTypeList& params = header.SummaryParameters();
	TagValuePairType pair;
	pair.Tag = "tag1";
	pair.Value = "value1";
	params.push_back(pair);
	pair.Tag = "tag2";
	pair.Value = "value2";
	params.push_back(pair);
	CPPUNIT_ASSERT(header.GetSummaryParameter("tag1") == "value1");
	CPPUNIT_ASSERT(header.GetSummaryParameter("tag2") == "value2");
	int count;
	TagValuePairTypeList::iterator iter;
	for (count = 0, iter = params.begin(); iter != params.end(); iter++, count++)
	{
		switch (count)
		{
		case 0:
		CPPUNIT_ASSERT( iter->Tag == "tag1");
		CPPUNIT_ASSERT( iter->Value == "value1");
		break;
		case 1:
		CPPUNIT_ASSERT( iter->Tag == "tag2");
		CPPUNIT_ASSERT( iter->Value == "value2");
		break;
		}
	}
}

void CCHPFileDataTest::testproperty_Header_ChipType()
{
	affxchp::CCHPFileHeader header;
	header.SetChipType("chip type");
	CPPUNIT_ASSERT( header.GetChipType() == "chip type");
}

void CCHPFileDataTest::testproperty_Header_ParentCellFile()
{
	affxchp::CCHPFileHeader header;
	header.SetParentCellFile("test.cel");
	CPPUNIT_ASSERT( header.GetParentCellFile() == "test.cel");
}

void CCHPFileDataTest::testproperty_Header_ProgID()
{
	affxchp::CCHPFileHeader header;
	header.SetProgID("prog_id");
	CPPUNIT_ASSERT( header.GetProgID() == "prog_id");
}

void CCHPFileDataTest::testproperty_Header_BackgroundZoneInfo()
{
	affxchp::CCHPFileHeader header;

	BackgroundZoneInfo& info = header.GetBackgroundZoneInfo();
	info.number_zones = 2;
	info.smooth_factor = 2.3f;
	BackgroundZoneTypeList& zones = header.GetBackgroundZones();
	BackgroundZoneType zone;
	zone.centerx = 1;
	zone.centery = 2;
	zone.background = 2.0f;
	zones.push_back(zone);
	zone.centerx = 3;
	zone.centery = 4;
	zone.background = 3.0f;
	zones.push_back(zone);
	CPPUNIT_ASSERT(CompareFloats(header.GetBackgroundZone(1,2).background, 2.0f));
	CPPUNIT_ASSERT(CompareFloats(header.GetBackgroundZone(3,4).background, 3.0f));
}

void CCHPFileDataTest::testproperty_Paths()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(EXP_COMP_CHP_FILE);
	CPPUNIT_ASSERT(chp.GetFileName() == EXP_COMP_CHP_FILE);
}

void CCHPFileDataTest::testmethod_ReadHeader_fail()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(chp.ReadHeader() == false);
}

void CCHPFileDataTest::testmethod_Read_fail()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(chp.Read() == false);
}

void CCHPFileDataTest::testmethod_Exists()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(NO_FILE);
	CPPUNIT_ASSERT(chp.Exists() == false);
	chp.SetFileName(EXP_COMP_CHP_FILE);
	CPPUNIT_ASSERT(chp.Exists() == true);
	chp.SetFileName(EXP_ABS_CHP_FILE);
	CPPUNIT_ASSERT(chp.Exists() == true);
}

void CCHPFileDataTest::testmethod_IsXDACompatibleFile()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(EXP_COMP_CHP_FILE);
	CPPUNIT_ASSERT(chp.IsXDACompatibleFile() == true);
	chp.SetFileName(EXP_ABS_CHP_FILE);
	CPPUNIT_ASSERT(chp.IsXDACompatibleFile() == true);
}

void CCHPFileDataTest::testmethod_Read_Exp_Comp()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(EXP_COMP_CHP_FILE);
	CPPUNIT_ASSERT( chp.Read() == true );
	CPPUNIT_ASSERT( chp.GetFileName() == EXP_COMP_CHP_FILE);

	// Accessors for header information.
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 120);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 120);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 8);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "TestExon");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Plier");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("p1") == "1.1.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("p2") == "2.1.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cp1") == "1.2.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cp2") == "2.2.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "TestExon.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "test_id");

	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 5);
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZoneInfo().smooth_factor, 3.8f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,0).background, 3.2f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(88,33).background, 8.9f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,60).background, 13.1f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(60,0).background, 5.0f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(120,120).background, 10.8f));
	int count;
	BackgroundZoneTypeList& zones = chp.GetHeader().GetBackgroundZones();
	BackgroundZoneTypeList::iterator iter;
	for (count = 0, iter = zones.begin(); iter != zones.end(); iter++, count++)
	{
		switch (count)
		{
		case 0:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 3.2f));
		break;
		case 1:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 88.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 33.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 8.9f));
		break;
		case 2:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 60.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 13.1f));
		break;
		case 3:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 60.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 5.0f));
		break;
		case 4:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 10.8f));
		break;
		}
	}

	for (int i = 0; i < chp.GetHeader().GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->DetectionPValue, (float)(0.05 - (i / 1000.0))));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->Signal, (float)(1.1 + i)));
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumPairs == 3 + i);
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumUsedPairs == 2 + i);
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->Detection == (i % 4));
		switch (chp.GetExpressionResults(i)->Detection)
		{
		case (ABS_PRESENT_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "P");
		break;
		case (ABS_MARGINAL_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "M");
		break;
		case (ABS_ABSENT_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "A");
		break;
		case (ABS_NO_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "No Call");
		break;
		default:
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "");
		break;
		}
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->m_HasCompResults == true);
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->ChangePValue, (float)(0.04 - (i / 1000.0))));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->SignalLogRatio, (float)(1.1 + i)));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->SignalLogRatioLow, (float)(-1.1 + i)));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->SignalLogRatioHigh, (float)(10.1 + i)));
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumCommonPairs == 2 + i);
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->Change == (i % 6 + 1));
		switch (chp.GetExpressionResults(i)->Change)
		{
		case (COMP_INCREASE_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetChangeString() == "I");
		break;
		case (COMP_DECREASE_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetChangeString() == "D");
		break;
		case (COMP_MOD_INCREASE_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetChangeString() == "MI");
		break;
		case (COMP_MOD_DECREASE_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetChangeString() == "MD");
		break;
		case (COMP_NO_CHANGE_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetChangeString() == "NC");
		break;
		case (COMP_NO_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetChangeString() == "No Call");
		break;
		default:
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetChangeString() == "");
		break;
		}
	}
}

void CCHPFileDataTest::testmethod_Read_Exp_Abs()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(EXP_ABS_CHP_FILE);
	CPPUNIT_ASSERT( chp.Read() == true );
	CPPUNIT_ASSERT( chp.GetFileName() == EXP_ABS_CHP_FILE);

	// Accessors for header information.
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 120);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 120);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 8);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "TestExon");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Plier");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("p1") == "1.1.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("p2") == "2.1.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cp1") == "1.2.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cp2") == "2.2.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "TestExon.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "test_id");

	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 5);
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZoneInfo().smooth_factor, 3.8f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,0).background, 3.2f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(88,33).background, 8.9f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,60).background, 13.1f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(60,0).background, 5.0f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(120,120).background, 10.8f));
	int count;
	BackgroundZoneTypeList& zones = chp.GetHeader().GetBackgroundZones();
	BackgroundZoneTypeList::iterator iter;
	for (count = 0, iter = zones.begin(); iter != zones.end(); iter++, count++)
	{
		switch (count)
		{
		case 0:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 3.2f));
		break;
		case 1:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 88.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 33.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 8.9f));
		break;
		case 2:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 60.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 13.1f));
		break;
		case 3:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 60.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 5.0f));
		break;
		case 4:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 10.8f));
		break;
		}
	}

	for (int i = 0; i < chp.GetHeader().GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->DetectionPValue, (float)(0.05 - (i / 1000.0))));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->Signal, (float)(1.1 + i)));
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumPairs == 3 + i);
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumUsedPairs == 2 + i);
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->Detection == (i % 4));
		switch (chp.GetExpressionResults(i)->Detection)
		{
		case (ABS_PRESENT_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "P");
		break;
		case (ABS_MARGINAL_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "M");
		break;
		case (ABS_ABSENT_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "A");
		break;
		case (ABS_NO_CALL):
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "No Call");
		break;
		default:
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->GetDetectionString() == "");
		break;
		}
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->m_HasCompResults == false);
	}
}

void CCHPFileDataTest::testmethod_ReadHeader_Exp_Comp()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(EXP_COMP_CHP_FILE);
	CPPUNIT_ASSERT( chp.Read() == true );
	CPPUNIT_ASSERT( chp.GetFileName() == EXP_COMP_CHP_FILE);

	// Accessors for header information.
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 120);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 120);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 8);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "TestExon");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Plier");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("p1") == "1.1.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("p2") == "2.1.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cp1") == "1.2.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cp2") == "2.2.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "TestExon.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "test_id");

	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 5);
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZoneInfo().smooth_factor, 3.8f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,0).background, 3.2f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(88,33).background, 8.9f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,60).background, 13.1f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(60,0).background, 5.0f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(120,120).background, 10.8f));
	int count;
	BackgroundZoneTypeList& zones = chp.GetHeader().GetBackgroundZones();
	BackgroundZoneTypeList::iterator iter;
	for (count = 0, iter = zones.begin(); iter != zones.end(); iter++, count++)
	{
		switch (count)
		{
		case 0:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 3.2f));
		break;
		case 1:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 88.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 33.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 8.9f));
		break;
		case 2:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 60.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 13.1f));
		break;
		case 3:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 60.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 5.0f));
		break;
		case 4:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 10.8f));
		break;
		}
	}
}

void CCHPFileDataTest::testmethod_ReadHeader_Exp_Abs()
{
	affxchp::CCHPFileData chp;
	chp.SetFileName(EXP_ABS_CHP_FILE);
	CPPUNIT_ASSERT( chp.Read() == true );
	CPPUNIT_ASSERT( chp.GetFileName() == EXP_ABS_CHP_FILE);

	// Accessors for header information.
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 120);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 120);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 8);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "TestExon");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Plier");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("p1") == "1.1.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("p2") == "2.1.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cp1") == "1.2.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cp2") == "2.2.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "TestExon.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "test_id");

	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 5);
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZoneInfo().smooth_factor, 3.8f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,0).background, 3.2f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(88,33).background, 8.9f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,60).background, 13.1f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(60,0).background, 5.0f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(120,120).background, 10.8f));
	int count;
	BackgroundZoneTypeList& zones = chp.GetHeader().GetBackgroundZones();
	BackgroundZoneTypeList::iterator iter;
	for (count = 0, iter = zones.begin(); iter != zones.end(); iter++, count++)
	{
		switch (count)
		{
		case 0:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 3.2f));
		break;
		case 1:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 88.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 33.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 8.9f));
		break;
		case 2:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 60.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 13.1f));
		break;
		case 3:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 60.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 0.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 5.0f));
		break;
		case 4:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 10.8f));
		break;
		}
	}
}

void CCHPFileDataTest::testmethod_ReadHeader_TagV11()
{
	CCHPFileData chp;
	chp.SetFileName(TAG_V11_FILE);
	CPPUNIT_ASSERT(chp.ReadHeader() == true);

	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 11 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 105);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 105);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 2050);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "GenFlex");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Hybridization");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "4.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "c:\\genechip\\testdata\\universal\\1109-A-05a.CEL");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "GeneChipAnalysis.HybBaseCall.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 0);
	CPPUNIT_ASSERT( chp.GetHeader().SummaryParameters().size() == 0);
	CPPUNIT_ASSERT( chp.GetHeader().AlgorithmParameters().size() == 0);
}

void CCHPFileDataTest::testmethod_Read_TagV11()
{
	CCHPFileData chp;
	chp.SetFileName(TAG_V11_FILE);
	CPPUNIT_ASSERT(chp.Read() == true);
	
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 11 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 105);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 105);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 2050);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "GenFlex");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Hybridization");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "4.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "c:\\genechip\\testdata\\universal\\1109-A-05a.CEL");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "GeneChipAnalysis.HybBaseCall.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 0);
	CPPUNIT_ASSERT( chp.GetHeader().SummaryParameters().size() == 0);
	CPPUNIT_ASSERT( chp.GetHeader().AlgorithmParameters().size() == 0);

	CPPUNIT_ASSERT( CompareFloats(chp.GetUniversalResults(0)->GetBackground(), 114.595f) == true);
	CPPUNIT_ASSERT( CompareFloats(chp.GetUniversalResults(1)->GetBackground(), 118.9f) == true);
	CPPUNIT_ASSERT( CompareFloats(chp.GetUniversalResults(2049)->GetBackground(), 114.854f) == true);
}

void CCHPFileDataTest::testmethod_ReadHeader_TagXDA()
{
	CCHPFileData chp;
	chp.SetFileName(TAG_XDA_FILE);
	CPPUNIT_ASSERT(chp.ReadHeader() == true);

	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 105);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 105);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 2050);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "GenFlex");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Hybridization");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "5.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "C:\\Program Files\\Affymetrix\\GeneChip\\Affy_Data\\Data\\1109-A-05a.CEL");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "GDMTAnalysis.HybBaseCall.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 0);
	CPPUNIT_ASSERT( chp.GetHeader().SummaryParameters().size() == 1);
	TagValuePairTypeList::iterator it = chp.GetHeader().SummaryParameters().begin();
	TagValuePairType param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "WaveLength");
	CPPUNIT_ASSERT( param.Value == "");
	CPPUNIT_ASSERT( chp.GetHeader().AlgorithmParameters().size() == 0);
}

void CCHPFileDataTest::testmethod_Read_TagXDA()
{
	CCHPFileData chp;
	chp.SetFileName(TAG_XDA_FILE);
	CPPUNIT_ASSERT(chp.Read() == true);
	
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 105);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 105);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 2050);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "GenFlex");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Hybridization");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "5.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "C:\\Program Files\\Affymetrix\\GeneChip\\Affy_Data\\Data\\1109-A-05a.CEL");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "GDMTAnalysis.HybBaseCall.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 0);
	CPPUNIT_ASSERT( chp.GetHeader().SummaryParameters().size() == 1);
	TagValuePairTypeList::iterator it = chp.GetHeader().SummaryParameters().begin();
	TagValuePairType param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "WaveLength");
	CPPUNIT_ASSERT( param.Value == "");
	CPPUNIT_ASSERT( chp.GetHeader().AlgorithmParameters().size() == 0);

	CPPUNIT_ASSERT( CompareFloats(chp.GetUniversalResults(0)->GetBackground(), 114.59585f) == true);
	CPPUNIT_ASSERT( CompareFloats(chp.GetUniversalResults(1)->GetBackground(), 118.9f) == true);
	CPPUNIT_ASSERT( CompareFloats(chp.GetUniversalResults(2049)->GetBackground(), 114.85418f) == true);
}

void CCHPFileDataTest::testmethod_Read_ReseqXDA_v1()
{
	CCHPFileData chp;
	chp.SetFileName(RESEQ_V1_FILE);
	CPPUNIT_ASSERT(chp.Read() == true);

	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 488);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 639);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 2);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "DCNtagIQr510989");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "CustomSeq");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "2");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "C:\\Data\\GCOS\\Data\\5303_DCN_01.CEL");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "GDMTAnalysis.VDABaseCall.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 0);
	CPPUNIT_ASSERT( chp.GetHeader().SummaryParameters().size() == 0);
	CPPUNIT_ASSERT( chp.GetHeader().AlgorithmParameters().size() == 10);
	TagValuePairTypeList::iterator it = chp.GetHeader().AlgorithmParameters().begin();
	TagValuePairType param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "NoSignal");
	CPPUNIT_ASSERT( param.Value == "1.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "WeakSignal");
	CPPUNIT_ASSERT( param.Value == "20.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "AberrantSNR2");
	CPPUNIT_ASSERT( param.Value == "20.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "StrandLLR");
	CPPUNIT_ASSERT( param.Value == "0.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "TotalLLR");
	CPPUNIT_ASSERT( param.Value == "75.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "PerfectCallThreshold");
	CPPUNIT_ASSERT( param.Value == "2.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "ModelType");
	CPPUNIT_ASSERT( param.Value == "0");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "FinalMaxHet");
	CPPUNIT_ASSERT( param.Value == "0.900000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "NeighborhoodRule");
	CPPUNIT_ASSERT( param.Value == "0.500000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "SampleReliability");
	CPPUNIT_ASSERT( param.Value == "0.750000");

	CResequencingResults *p = chp.GetResequencingResults();
	CPPUNIT_ASSERT(p != NULL);

	int n = p->GetCalledBasesSize();
	CPPUNIT_ASSERT(n == 27670);
	CPPUNIT_ASSERT(p->GetCalledBase(0) == 'c');
	CPPUNIT_ASSERT(p->GetCalledBase(1) == 'c');
	CPPUNIT_ASSERT(p->GetCalledBase(2) == 'c');
	CPPUNIT_ASSERT(p->GetCalledBase(3) == 'a');
	CPPUNIT_ASSERT(p->GetCalledBase(4) == 'g');
	CPPUNIT_ASSERT(p->GetCalledBase(5) == 'n');
	CPPUNIT_ASSERT(p->GetCalledBase(6) == 'c');

	CPPUNIT_ASSERT(p->GetCalledBase(n-1) == 'n');
	CPPUNIT_ASSERT(p->GetCalledBase(n-44) == 'a');

	const double eps = 1e-5;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(0), 66.465454, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(1), 120.09094, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(2), 153.29333, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(3), 120.06390, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(n-1), 0.0, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(n-6), 75.680908, eps);

	CPPUNIT_ASSERT(p->GetForceCallsSize() == 0);
	CPPUNIT_ASSERT(p->GetOrigCallsSize() == 0);
}

void CCHPFileDataTest::testmethod_Read_ReseqXDA_v2()
{
	CCHPFileData chp;
	chp.SetFileName(RESEQ_V2_FILE);
	CPPUNIT_ASSERT(chp.Read() == true);

	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 2 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 960);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 1008);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 2);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "Tristezar520098");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Reseq2");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "2.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "C:\\GeneChip\\Affy_Data\\Data\\Hyb01004 CTV-T36 Expt 1926.CEL");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "GDMTAnalysis.Reseq2BaseCall.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 0);
	CPPUNIT_ASSERT( chp.GetHeader().SummaryParameters().size() == 0);
	CPPUNIT_ASSERT( chp.GetHeader().AlgorithmParameters().size() == 1);
	TagValuePairTypeList::iterator it = chp.GetHeader().AlgorithmParameters().begin();
	TagValuePairType param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "QualityScore");
	CPPUNIT_ASSERT( param.Value == "3");

	CResequencingResults *p = chp.GetResequencingResults();
	CPPUNIT_ASSERT(p != NULL);

	int n = p->GetCalledBasesSize();
	CPPUNIT_ASSERT(n == 111830);
	CPPUNIT_ASSERT(p->GetCalledBase(0) == 'a');
	CPPUNIT_ASSERT(p->GetCalledBase(1) == 'c');
	CPPUNIT_ASSERT(p->GetCalledBase(2) == 'a');
	CPPUNIT_ASSERT(p->GetCalledBase(3) == 'g');
	CPPUNIT_ASSERT(p->GetCalledBase(4) == 'c');
	CPPUNIT_ASSERT(p->GetCalledBase(5) == 'g');
	CPPUNIT_ASSERT(p->GetCalledBase(6) == 'a');

	CPPUNIT_ASSERT(p->GetCalledBase(n-1) == 'g');
	CPPUNIT_ASSERT(p->GetCalledBase(n-2) == 'g');
	CPPUNIT_ASSERT(p->GetCalledBase(n-3) == 't');

	const double eps = 1e-5;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(0), 18.859673, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(1), 14.588547, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(2), 15.489632, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(3), 25.993202, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(n-1), 32.353516, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(n-2), 16.441238, eps);

	CPPUNIT_ASSERT(p->GetOrigCallsSize() == 0);

	n = p->GetForceCallsSize();
	CPPUNIT_ASSERT(n == 61065);

	ForceCallType force;

	force = p->GetForceCall(0);
	CPPUNIT_ASSERT(force.position == 24);
	CPPUNIT_ASSERT(force.call == 'n');
	CPPUNIT_ASSERT(force.reason == QUALITY_SCORE_THR_FORCE_CALL);

	force = p->GetForceCall(1);
	CPPUNIT_ASSERT(force.position == 31);
	CPPUNIT_ASSERT(force.call == 'n');
	CPPUNIT_ASSERT(force.reason == QUALITY_SCORE_THR_FORCE_CALL);

	force = p->GetForceCall(2);
	CPPUNIT_ASSERT(force.position == 39);
	CPPUNIT_ASSERT(force.call == 'n');
	CPPUNIT_ASSERT(force.reason == QUALITY_SCORE_THR_FORCE_CALL);

	force = p->GetForceCall(n-1);
	CPPUNIT_ASSERT(force.position == 111795);
	CPPUNIT_ASSERT(force.call == 'n');
	CPPUNIT_ASSERT(force.reason == QUALITY_SCORE_THR_FORCE_CALL);
}

void CCHPFileDataTest::testmethod_Read_Reseq_old_file()
{
	CCHPFileData chp;
	chp.SetFileName(RESEQ_OLD_FILE);
	CPPUNIT_ASSERT(chp.Read() == true);

	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 13 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 488);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 639);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 0);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "DCNtagIQr510989");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "CustomSeq");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "S:\\GDAS2\\Test_Files\\Resequence\\DCN Validation Data\\5303_DCN_01.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "GDMTAnalysis.VDABaseCall.1");
	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 0);
	CPPUNIT_ASSERT( chp.GetHeader().SummaryParameters().size() == 0);
	CPPUNIT_ASSERT( chp.GetHeader().AlgorithmParameters().size() == 10);
	TagValuePairTypeList::iterator it = chp.GetHeader().AlgorithmParameters().begin();
	TagValuePairType param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "SampleReliability");
	CPPUNIT_ASSERT( param.Value == "0.750000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "NeighborhoodRule");
	CPPUNIT_ASSERT( param.Value == "0.500000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "FinalMaxHet");
	CPPUNIT_ASSERT( param.Value == "0.900000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "ModelType");
	CPPUNIT_ASSERT( param.Value == "0");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "PerfectCallThreshold");
	CPPUNIT_ASSERT( param.Value == "2.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "TotalLLR");
	CPPUNIT_ASSERT( param.Value == "75.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "StrandLLR");
	CPPUNIT_ASSERT( param.Value == "0.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "AberrantSNR2");
	CPPUNIT_ASSERT( param.Value == "20.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "WeakSignal");
	CPPUNIT_ASSERT( param.Value == "20.000000");
	++it;
	param = (TagValuePairType) *it;
	CPPUNIT_ASSERT( param.Tag == "NoSignal");
	CPPUNIT_ASSERT( param.Value == "1.000000");


	CResequencingResults *p = chp.GetResequencingResults();
	CPPUNIT_ASSERT(p != NULL);

	int n = p->GetCalledBasesSize();
	CPPUNIT_ASSERT(n == 29270);
	CPPUNIT_ASSERT(p->GetCalledBase(0) == 'c');
	CPPUNIT_ASSERT(p->GetCalledBase(1) == 'c');
	CPPUNIT_ASSERT(p->GetCalledBase(2) == 'c');
	CPPUNIT_ASSERT(p->GetCalledBase(3) == 'a');
	CPPUNIT_ASSERT(p->GetCalledBase(4) == 'g');
	CPPUNIT_ASSERT(p->GetCalledBase(5) == 't');
	CPPUNIT_ASSERT(p->GetCalledBase(6) == 'c');

	CPPUNIT_ASSERT(p->GetCalledBase(n-1) == 'n');
	CPPUNIT_ASSERT(p->GetCalledBase(n-44) == 'a');


	const double eps = 1e-5;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(0), 90.058533, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(1), 130.95447, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(2), 166.22278, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(3), 143.23492, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(n-1), 0.0, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(n-2), 0.0, eps);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(p->GetScore(n-6), 75.936035, eps);

	CPPUNIT_ASSERT(p->GetForceCallsSize() == 0);
	CPPUNIT_ASSERT(p->GetOrigCallsSize() == 0);

}
