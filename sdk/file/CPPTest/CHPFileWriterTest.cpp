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

#include "file/CPPTest/CHPFileWriterTest.h"
//
#include "file/CDFFileData.h"
#include "file/CHPFileWriter.h"
//
#include <cmath>
//

// Delimiter character in DAT header 
#define DELIMCHAR 0x14

CPPUNIT_TEST_SUITE_REGISTRATION( CCHPFileWriterTest );

#ifdef _MSC_VER
#define EXP_COMP_OUT_FILE ".\\data\\test.exp.comp.out.CHP"
#define EXP_ABS_OUT_FILE ".\\data\\test.exp.abs.out.CHP"
#else
#define EXP_COMP_OUT_FILE "./data/test.exp.comp.out.CHP"
#define EXP_ABS_OUT_FILE "./data/test.exp.abs.out.CHP"
#endif

#define GENO_OUT_FILE "./data/test.geno.out.CHP"
#define UNIV_OUT_FILE "./data/test.univ.out.CHP"
#define RESEQ_OUT_FILE "./data/test.reseq.out.CHP"

using namespace affxchp;

static bool CompareFloats(float f1, float f2)
{
	const float EPS = 0.0000001f;
	return (fabs(f1-f2) < EPS);
}

void CCHPFileWriterTest::setUp()
{
}

void CCHPFileWriterTest::tearDown()
{
}


void CCHPFileWriterTest::testfunction_Write_Exp_Abs_NoBuffer()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.SetAlgorithmName("Plier");
	chp.SetAlgorithmVersion("1.0");
	chp.SetParentCelFileName("TestChip.cel");
	chp.SetProgID("prog_id");
	chp.AddAlgorithmParameter("param1", "1.1.0");
	chp.AddAlgorithmParameter("param2", "2.1.0");
	chp.AddChipSummaryParameter("cparam1", "1.2.0");
	chp.AddChipSummaryParameter("cparam2", "2.2.0");
	chp.AddBackgroundInfo(3, 2.8f);
	chp.AddBackgroundZone(0, 0, 3.2f);
	chp.AddBackgroundZone(60, 80, 8.0f);
	chp.AddBackgroundZone(120, 120, 10.8f);
	
	chp.InitializeForWriting(150, 150, 5, "TestChip", affxcdf::ExpressionProbeSetType, false);

	chp.SetFileName(EXP_ABS_OUT_FILE);
	CPPUNIT_ASSERT(chp.CreateNewFile() == true);
	CPPUNIT_ASSERT(chp.SaveHeader() == true);

	unsigned char det[] = {ABS_PRESENT_CALL, ABS_MARGINAL_CALL, ABS_ABSENT_CALL, ABS_NO_CALL, ABS_PRESENT_CALL};

	for (int i = 0; i < 5; i++)
	{
		CExpressionProbeSetResults entry;
		entry.DetectionPValue = (float)(0.03 - (i / 1000.0));
		entry.Signal = (float) (2.1 + i);
		entry.NumPairs = 4 + i;
		entry.NumUsedPairs = 3 + i;
		entry.m_HasCompResults = false;
		entry.Detection = det[i];

		chp.SaveExpressionEntry(&entry);
	}
	CPPUNIT_ASSERT(chp.Close());


	chp.Clear();
	chp.SetFileName(EXP_ABS_OUT_FILE);
	CPPUNIT_ASSERT( chp.Read() == true);
	CPPUNIT_ASSERT( chp.GetFileName() == EXP_ABS_OUT_FILE);

	// Accessors for header information.
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 150);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 150);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 5);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "TestChip");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Plier");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("param1") == "1.1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("param2") == "2.1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cparam1") == "1.2.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cparam2") == "2.2.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "TestChip.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "prog_id");

	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 3);
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZoneInfo().smooth_factor, 2.8f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,0).background, 3.2f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(60,80).background, 8.0f));
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
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 60.f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 80.f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 8.0f));
		break;
		case 2:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 10.8f));
		break;
		}
	}

	for (int i = 0; i < chp.GetHeader().GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->DetectionPValue, (float)(0.03 - (i / 1000.0))));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->Signal, (float)(2.1 + i)));
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumPairs == 4 + i);
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumUsedPairs == 3 + i);
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

void CCHPFileWriterTest::testfunction_Write_Genotyping_NoBuffer()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.InitializeForWriting(150, 150, 2, "TestChip", affxcdf::GenotypingProbeSetType, false);
	chp.SetAlgorithmName("Geno");
	chp.SetAlgorithmVersion("1.0");
	chp.SetParentCelFileName("TestChip.cel");
	chp.SetProgID("prog_id");
	chp.AddAlgorithmParameter("param1", "1.1.0");
	chp.AddChipSummaryParameter("cparam1", "1.2.0");

	chp.SetFileName(GENO_OUT_FILE);
	CPPUNIT_ASSERT(chp.CreateNewFile() == true);
	CPPUNIT_ASSERT(chp.SaveHeader() == true);

	affxchp::CGenotypeProbeSetResults entry;
	entry.AlleleCall = ALLELE_A_CALL;
	entry.Confidence = 1.0f;
	entry.pvalue_AA = 0.1f;
	entry.pvalue_AB = 0.2f;
	entry.pvalue_BB = 0.3f;
	entry.pvalue_NoCall = 0.4f;
	chp.SaveMappingEntry(&entry);

	entry.AlleleCall = ALLELE_AB_CALL;
	entry.Confidence = 2.0f;
	entry.pvalue_AA = 0.11f;
	entry.pvalue_AB = 0.12f;
	entry.pvalue_BB = 0.13f;
	entry.pvalue_NoCall = 0.14f;
	chp.SaveMappingEntry(&entry);
	CPPUNIT_ASSERT(chp.Close());

	chp.Clear();


	// Read the results and compare.
	affxchpwriter::CCHPFileWriter chp2;
	chp2.SetFileName(GENO_OUT_FILE);
	CPPUNIT_ASSERT( chp2.Read() == true);
	CPPUNIT_ASSERT( chp2.GetFileName() == GENO_OUT_FILE);

	CPPUNIT_ASSERT( chp2.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp2.GetHeader().GetCols() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetRows() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetNumProbeSets() == 2);
	CPPUNIT_ASSERT( chp2.GetHeader().GetChipType() == "TestChip");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgName() == "GenotypingGeno");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgorithmParameter("param1") == "1.1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetSummaryParameter("cparam1") == "1.2.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetParentCellFile() == "TestChip.cel");
	CPPUNIT_ASSERT( chp2.GetHeader().GetProgID() == "prog_id");

	affxchp::CGenotypeProbeSetResults *pEntry;

	pEntry = chp2.GetGenotypingResults(0);
	CPPUNIT_ASSERT( pEntry->AlleleCall == ALLELE_A_CALL);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->Confidence, 1.0f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_AA, 0.1f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_AB, 0.2f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_BB, 0.3f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_NoCall, 0.4f, 0.001f);

	pEntry = chp2.GetGenotypingResults(1);
	CPPUNIT_ASSERT( pEntry->AlleleCall == ALLELE_AB_CALL);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->Confidence, 2.0f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_AA, 0.11f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_AB, 0.12f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_BB, 0.13f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_NoCall, 0.14f, 0.001f);
}

void CCHPFileWriterTest::testfunction_Write_Universal_NoBuffer()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.InitializeForWriting(150, 150, 2, "TestChip", affxcdf::TagProbeSetType, false);
	chp.SetAlgorithmName("Tag");
	chp.SetAlgorithmVersion("1.0");
	chp.SetParentCelFileName("TestChip.cel");
	chp.SetProgID("prog_id");
	chp.AddAlgorithmParameter("param1", "1.1.0");
	chp.AddChipSummaryParameter("cparam1", "1.2.0");

	chp.SetFileName(UNIV_OUT_FILE);
	CPPUNIT_ASSERT(chp.CreateNewFile() == true);
	CPPUNIT_ASSERT(chp.SaveHeader() == true);

	affxchp::CUniversalProbeSetResults entry;
	entry.SetBackground(1.0f);
	chp.SaveUniversalEntry(&entry);

	entry.SetBackground(2.0f);
	chp.SaveUniversalEntry(&entry);
	CPPUNIT_ASSERT(chp.Close());

	chp.Clear();

	// Read the results and compare.
	affxchpwriter::CCHPFileWriter chp2;
	chp2.SetFileName(UNIV_OUT_FILE);
	CPPUNIT_ASSERT( chp2.Read() == true);
	CPPUNIT_ASSERT( chp2.GetFileName() == UNIV_OUT_FILE);

	CPPUNIT_ASSERT( chp2.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp2.GetHeader().GetCols() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetRows() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetNumProbeSets() == 2);
	CPPUNIT_ASSERT( chp2.GetHeader().GetChipType() == "TestChip");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgName() == "Tag");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgorithmParameter("param1") == "1.1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetSummaryParameter("cparam1") == "1.2.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetParentCellFile() == "TestChip.cel");
	CPPUNIT_ASSERT( chp2.GetHeader().GetProgID() == "prog_id");

	affxchp::CUniversalProbeSetResults *pEntry;

	pEntry = chp2.GetUniversalResults(0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( pEntry->GetBackground(), 1.0f, 0.0001f);
	pEntry = chp2.GetUniversalResults(1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( pEntry->GetBackground(), 2.0f, 0.0001f);
}

void CCHPFileWriterTest::testproperty_Header_AlgName()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.SetAlgorithmName("alg");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "alg");
}

void CCHPFileWriterTest::testproperty_Header_AlgVersion()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.SetAlgorithmVersion("1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1.0");
}

void CCHPFileWriterTest::testproperty_Header_AlgParams()
{
	affxchpwriter::CCHPFileWriter chp;

	chp.AddAlgorithmParameter("tag1", "value1");
	chp.AddAlgorithmParameter("tag2", "value2");
	CPPUNIT_ASSERT(chp.GetHeader().GetAlgorithmParameter("tag1") == "value1");
	CPPUNIT_ASSERT(chp.GetHeader().GetAlgorithmParameter("tag2") == "value2");
	int count;
	TagValuePairTypeList::iterator iter;
	TagValuePairTypeList& params = chp.GetHeader().AlgorithmParameters();
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

void CCHPFileWriterTest::testproperty_Header_SummaryParams()
{
	affxchpwriter::CCHPFileWriter chp;

	chp.AddChipSummaryParameter("tag1", "value1");
	chp.AddChipSummaryParameter("tag2", "value2");
	CPPUNIT_ASSERT(chp.GetHeader().GetSummaryParameter("tag1") == "value1");
	CPPUNIT_ASSERT(chp.GetHeader().GetSummaryParameter("tag2") == "value2");
	int count;
	TagValuePairTypeList::iterator iter;
	TagValuePairTypeList& params = chp.GetHeader().SummaryParameters();
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

void CCHPFileWriterTest::testproperty_Header_ParentCellFile()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.SetParentCelFileName("test.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "test.cel");
}

void CCHPFileWriterTest::testproperty_Header_ProgID()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.SetProgID("prog_id");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "prog_id");
}

void CCHPFileWriterTest::testproperty_Header_BackgroundZoneInfo()
{
	affxchpwriter::CCHPFileWriter chp;

	chp.AddBackgroundInfo(2, 2.3f);
	chp.AddBackgroundZone(1, 2, 2.0f);
	chp.AddBackgroundZone(3, 4, 3.0f);
	CPPUNIT_ASSERT(CompareFloats(chp.GetHeader().GetBackgroundZone(1,2).background, 2.0f));
	CPPUNIT_ASSERT(CompareFloats(chp.GetHeader().GetBackgroundZone(3,4).background, 3.0f));
}

void CCHPFileWriterTest::testfunction_Write_Exp_Comp()
{
	affxchpwriter::CCHPFileWriter chp;
	setData(chp, true);
	chp.SetFileName(EXP_COMP_OUT_FILE);
	CPPUNIT_ASSERT(chp.CreateNewFile() == true);
	CPPUNIT_ASSERT(chp.Save() == true);

	chp.Clear();
	chp.SetFileName(EXP_COMP_OUT_FILE);
	CPPUNIT_ASSERT( chp.Read() == true);
	CPPUNIT_ASSERT( chp.GetFileName() == EXP_COMP_OUT_FILE);

	// Accessors for header information.
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 150);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 150);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 5);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "TestChip");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Plier");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("param1") == "1.1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("param2") == "2.1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cparam1") == "1.2.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cparam2") == "2.2.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "TestChip.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "prog_id");

	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 3);
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZoneInfo().smooth_factor, 2.8f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,0).background, 3.2f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(60,80).background, 8.0f));
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
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 60.f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 80.f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 8.0f));
		break;
		case 2:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 10.8f));
		break;
		}
	}

	for (int i = 0; i < chp.GetHeader().GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->DetectionPValue, (float)(0.03 - (i / 1000.0))));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->Signal, (float)(2.1 + i)));
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumPairs == 4 + i);
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumUsedPairs == 3 + i);
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
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->ChangePValue, (float)(0.02 - (i / 1000.0))));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->SignalLogRatio, (float)(3.3 + i)));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->SignalLogRatioLow, (float)(-2.1 + i)));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->SignalLogRatioHigh, (float)(13.9 + i)));
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

void CCHPFileWriterTest::testfunction_Write_Exp_Abs()
{
	affxchpwriter::CCHPFileWriter chp;
	setData(chp, false);
	chp.SetFileName(EXP_ABS_OUT_FILE);
	CPPUNIT_ASSERT(chp.CreateNewFile() == true);
	CPPUNIT_ASSERT(chp.Save() == true);

	chp.Clear();
	chp.SetFileName(EXP_ABS_OUT_FILE);
	CPPUNIT_ASSERT( chp.Read() == true);
	CPPUNIT_ASSERT( chp.GetFileName() == EXP_ABS_OUT_FILE);

	// Accessors for header information.
	CPPUNIT_ASSERT( chp.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp.GetHeader().GetCols() == 150);
	CPPUNIT_ASSERT( chp.GetHeader().GetRows() == 150);
	CPPUNIT_ASSERT( chp.GetHeader().GetNumProbeSets() == 5);
	CPPUNIT_ASSERT( chp.GetHeader().GetChipType() == "TestChip");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgName() == "Plier");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("param1") == "1.1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgorithmParameter("param2") == "2.1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cparam1") == "1.2.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetSummaryParameter("cparam2") == "2.2.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp.GetHeader().GetParentCellFile() == "TestChip.cel");
	CPPUNIT_ASSERT( chp.GetHeader().GetProgID() == "prog_id");

	CPPUNIT_ASSERT( chp.GetHeader().GetBackgroundZoneInfo().number_zones == 3);
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZoneInfo().smooth_factor, 2.8f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(0,0).background, 3.2f));
	CPPUNIT_ASSERT( CompareFloats(chp.GetHeader().GetBackgroundZone(60,80).background, 8.0f));
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
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 60.f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 80.f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 8.0f));
		break;
		case 2:
		CPPUNIT_ASSERT( CompareFloats(iter->centerx , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->centery , 120.0f));
		CPPUNIT_ASSERT( CompareFloats(iter->background , 10.8f));
		break;
		}
	}

	for (int i = 0; i < chp.GetHeader().GetNumProbeSets(); i++)
	{
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->DetectionPValue, (float)(0.03 - (i / 1000.0))));
		CPPUNIT_ASSERT(CompareFloats(chp.GetExpressionResults(i)->Signal, (float)(2.1 + i)));
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumPairs == 4 + i);
		CPPUNIT_ASSERT(chp.GetExpressionResults(i)->NumUsedPairs == 3 + i);
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

void CCHPFileWriterTest::setData(CCHPFileWriter& chp, bool bComparison)
{
	chp.SetAlgorithmName("Plier");
	chp.SetAlgorithmVersion("1.0");
	chp.SetParentCelFileName("TestChip.cel");
	chp.SetProgID("prog_id");
	chp.AddAlgorithmParameter("param1", "1.1.0");
	chp.AddAlgorithmParameter("param2", "2.1.0");
	chp.AddChipSummaryParameter("cparam1", "1.2.0");
	chp.AddChipSummaryParameter("cparam2", "2.2.0");
	chp.AddBackgroundInfo(3, 2.8f);
	chp.AddBackgroundZone(0, 0, 3.2f);
	chp.AddBackgroundZone(60, 80, 8.0f);
	chp.AddBackgroundZone(120, 120, 10.8f);
	
	chp.InitializeForWriting(150, 150, 5, "TestChip", affxcdf::ExpressionProbeSetType);
	
	for (int i = 0; i < 5; i++)
	{
		chp.GetExpressionResults(i)->DetectionPValue = (float)(0.03 - (i / 1000.0));
		chp.GetExpressionResults(i)->Signal = (float) (2.1 + i);
		chp.GetExpressionResults(i)->NumPairs = 4 + i;
		chp.GetExpressionResults(i)->NumUsedPairs = 3 + i;
		if (bComparison)
		{ 
			chp.GetExpressionResults(i)->m_HasCompResults = true;
			chp.GetExpressionResults(i)->ChangePValue = (float)(0.02 - (i / 1000.0));
			chp.GetExpressionResults(i)->SignalLogRatio = (float)(3.3 + i);
			chp.GetExpressionResults(i)->SignalLogRatioLow = (float)(-2.1 + i);
			chp.GetExpressionResults(i)->SignalLogRatioHigh = (float)(13.9 + i);
			chp.GetExpressionResults(i)->NumCommonPairs = 2 + i;
		}
		else
		{
			chp.GetExpressionResults(i)->m_HasCompResults = false;
		}
	}
	chp.GetExpressionResults(0)->Detection = ABS_PRESENT_CALL;
	chp.GetExpressionResults(1)->Detection = ABS_MARGINAL_CALL;
	chp.GetExpressionResults(2)->Detection = ABS_ABSENT_CALL;
	chp.GetExpressionResults(3)->Detection = ABS_NO_CALL;
	chp.GetExpressionResults(4)->Detection = ABS_PRESENT_CALL;
	if (bComparison)
	{
		chp.GetExpressionResults(0)->Change = COMP_INCREASE_CALL;
		chp.GetExpressionResults(1)->Change = COMP_DECREASE_CALL;
		chp.GetExpressionResults(2)->Change = COMP_MOD_INCREASE_CALL;
		chp.GetExpressionResults(3)->Change = COMP_MOD_DECREASE_CALL;
		chp.GetExpressionResults(4)->Change = COMP_NO_CHANGE_CALL;
	}
}

void CCHPFileWriterTest::testassignment_CResequencingResults()
{
	CResequencingResults a;
	CResequencingResults b;

	a.ResizeCalledBases(2);
	a.SetCalledBase(0, 'a');
	a.SetCalledBase(1, 'c');

	a.ResizeScores(2);
	a.SetScore(0, 1.0f);
	a.SetScore(1, 2.0f);

	ForceCallType f;
	a.ResizeForceCalls(2);
	f.position = 0;
	f.call = 'g';
	f.reason = 'b';
	a.SetForceCall(0, f);
	f.position = 1;
	f.call = 't';
	f.reason = 'c';
	a.SetForceCall(1, f);

	BaseCallType o;
	a.ResizeOrigCalls(2);
	o.position = 0;
	o.call = 't';
	a.SetOrigCall(0, o);
	o.position = 1;
	o.call = 'g';
	a.SetOrigCall(1, o);

	b = a;
	a.Clear();

	CPPUNIT_ASSERT(b.GetCalledBasesSize() == 2);
	CPPUNIT_ASSERT(b.GetScoresSize() == 2);
	CPPUNIT_ASSERT(b.GetForceCallsSize() == 2);
	CPPUNIT_ASSERT(b.GetOrigCallsSize() == 2);

	CPPUNIT_ASSERT(b.GetCalledBase(0) == 'a');
	CPPUNIT_ASSERT(b.GetCalledBase(1) == 'c');

	CPPUNIT_ASSERT_DOUBLES_EQUAL(b.GetScore(0), 1.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(b.GetScore(1), 2.0f, 0.0001f);

	f = b.GetForceCall(0);
	CPPUNIT_ASSERT(f.position == 0);
	CPPUNIT_ASSERT(f.call == 'g');
	CPPUNIT_ASSERT(f.reason == 'b');

	f = b.GetForceCall(1);
	CPPUNIT_ASSERT(f.position == 1);
	CPPUNIT_ASSERT(f.call == 't');
	CPPUNIT_ASSERT(f.reason == 'c');

	o = b.GetOrigCall(0);
	CPPUNIT_ASSERT(o.position == 0);
	CPPUNIT_ASSERT(o.call == 't');

	o = b.GetOrigCall(1);
	CPPUNIT_ASSERT(o.position == 1);
	CPPUNIT_ASSERT(o.call == 'g');

}
void CCHPFileWriterTest::testfunction_Write_Genotyping()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.InitializeForWriting(150, 150, 2, "TestChip", affxcdf::GenotypingProbeSetType);
	chp.SetAlgorithmName("Geno");
	chp.SetAlgorithmVersion("1.0");
	chp.SetParentCelFileName("TestChip.cel");
	chp.SetProgID("prog_id");
	chp.AddAlgorithmParameter("param1", "1.1.0");
	chp.AddChipSummaryParameter("cparam1", "1.2.0");


	affxchp::CGenotypeProbeSetResults entry;
	entry.AlleleCall = ALLELE_A_CALL;
	entry.Confidence = 1.0f;
	entry.pvalue_AA = 0.1f;
	entry.pvalue_AB = 0.2f;
	entry.pvalue_BB = 0.3f;
	entry.pvalue_NoCall = 0.4f;
	chp.SetMappingEntry(0, &entry);

	entry.AlleleCall = ALLELE_AB_CALL;
	entry.Confidence = 2.0f;
	entry.pvalue_AA = 0.11f;
	entry.pvalue_AB = 0.12f;
	entry.pvalue_BB = 0.13f;
	entry.pvalue_NoCall = 0.14f;
	chp.SetMappingEntry(1, &entry);

	chp.SetFileName(GENO_OUT_FILE);
	CPPUNIT_ASSERT(chp.CreateNewFile() == true);
	CPPUNIT_ASSERT(chp.Save() == true);
	chp.Clear();


	// Read the results and compare.
	affxchpwriter::CCHPFileWriter chp2;
	chp2.SetFileName(GENO_OUT_FILE);
	CPPUNIT_ASSERT( chp2.Read() == true);
	CPPUNIT_ASSERT( chp2.GetFileName() == GENO_OUT_FILE);

	CPPUNIT_ASSERT( chp2.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp2.GetHeader().GetCols() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetRows() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetNumProbeSets() == 2);
	CPPUNIT_ASSERT( chp2.GetHeader().GetChipType() == "TestChip");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgName() == "GenotypingGeno");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgorithmParameter("param1") == "1.1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetSummaryParameter("cparam1") == "1.2.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetParentCellFile() == "TestChip.cel");
	CPPUNIT_ASSERT( chp2.GetHeader().GetProgID() == "prog_id");

	affxchp::CGenotypeProbeSetResults *pEntry;

	pEntry = chp2.GetGenotypingResults(0);
	CPPUNIT_ASSERT( pEntry->AlleleCall == ALLELE_A_CALL);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->Confidence, 1.0f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_AA, 0.1f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_AB, 0.2f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_BB, 0.3f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_NoCall, 0.4f, 0.001f);

	pEntry = chp2.GetGenotypingResults(1);
	CPPUNIT_ASSERT( pEntry->AlleleCall == ALLELE_AB_CALL);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->Confidence, 2.0f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_AA, 0.11f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_AB, 0.12f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_BB, 0.13f, 0.001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(pEntry->pvalue_NoCall, 0.14f, 0.001f);

}

void CCHPFileWriterTest::testfunction_Write_Universal()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.InitializeForWriting(150, 150, 2, "TestChip", affxcdf::TagProbeSetType);
	chp.SetAlgorithmName("Tag");
	chp.SetAlgorithmVersion("1.0");
	chp.SetParentCelFileName("TestChip.cel");
	chp.SetProgID("prog_id");
	chp.AddAlgorithmParameter("param1", "1.1.0");
	chp.AddChipSummaryParameter("cparam1", "1.2.0");


	affxchp::CUniversalProbeSetResults entry;
	entry.SetBackground(1.0f);
	chp.SetUniversalEntry(0, &entry);

	entry.SetBackground(2.0f);
	chp.SetUniversalEntry(1, &entry);

	chp.SetFileName(UNIV_OUT_FILE);
	CPPUNIT_ASSERT(chp.CreateNewFile() == true);
	CPPUNIT_ASSERT(chp.Save() == true);
	chp.Clear();




	// Read the results and compare.
	affxchpwriter::CCHPFileWriter chp2;
	chp2.SetFileName(UNIV_OUT_FILE);
	CPPUNIT_ASSERT( chp2.Read() == true);
	CPPUNIT_ASSERT( chp2.GetFileName() == UNIV_OUT_FILE);

	CPPUNIT_ASSERT( chp2.GetHeader().GetVersionNumber() == 1 );
	CPPUNIT_ASSERT( chp2.GetHeader().GetCols() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetRows() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetNumProbeSets() == 2);
	CPPUNIT_ASSERT( chp2.GetHeader().GetChipType() == "TestChip");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgName() == "Tag");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgorithmParameter("param1") == "1.1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetSummaryParameter("cparam1") == "1.2.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetParentCellFile() == "TestChip.cel");
	CPPUNIT_ASSERT( chp2.GetHeader().GetProgID() == "prog_id");

	affxchp::CUniversalProbeSetResults *pEntry;

	pEntry = chp2.GetUniversalResults(0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( pEntry->GetBackground(), 1.0f, 0.0001f);
	pEntry = chp2.GetUniversalResults(1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL( pEntry->GetBackground(), 2.0f, 0.0001f);

}

void CCHPFileWriterTest::testfunction_Write_Resequencing()
{
	affxchpwriter::CCHPFileWriter chp;
	chp.InitializeForWriting(150, 150, 2, "TestChip", affxcdf::ResequencingProbeSetType);
	chp.SetAlgorithmName("Tag");
	chp.SetAlgorithmVersion("1.0");
	chp.SetParentCelFileName("TestChip.cel");
	chp.SetProgID("prog_id");
	chp.AddAlgorithmParameter("param1", "1.1.0");
	chp.AddChipSummaryParameter("cparam1", "1.2.0");



	CResequencingResults a;
	a.ResizeCalledBases(2);
	a.SetCalledBase(0, 'a');
	a.SetCalledBase(1, 'c');

	a.ResizeScores(2);
	a.SetScore(0, 1.0f);
	a.SetScore(1, 2.0f);

	ForceCallType f;
	a.ResizeForceCalls(2);
	f.position = 0;
	f.call = 'g';
	f.reason = 'b';
	a.SetForceCall(0, f);
	f.position = 1;
	f.call = 't';
	f.reason = 'c';
	a.SetForceCall(1, f);

	BaseCallType o;
	a.ResizeOrigCalls(2);
	o.position = 0;
	o.call = 't';
	a.SetOrigCall(0, o);
	o.position = 1;
	o.call = 'g';
	a.SetOrigCall(1, o);

	chp.SetResequencingResults(&a);

	chp.SetFileName(RESEQ_OUT_FILE);
	CPPUNIT_ASSERT(chp.CreateNewFile() == true);
	CPPUNIT_ASSERT(chp.Save() == true);





	// Read the results and compare.
	affxchpwriter::CCHPFileWriter chp2;
	chp2.SetFileName(RESEQ_OUT_FILE);
	CPPUNIT_ASSERT( chp2.Read() == true);
	CPPUNIT_ASSERT( chp2.GetFileName() == RESEQ_OUT_FILE);

	CPPUNIT_ASSERT( chp2.GetHeader().GetVersionNumber() == 2 );
	CPPUNIT_ASSERT( chp2.GetHeader().GetCols() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetRows() == 150);
	CPPUNIT_ASSERT( chp2.GetHeader().GetNumProbeSets() == 2);
	CPPUNIT_ASSERT( chp2.GetHeader().GetChipType() == "TestChip");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgName() == "Tag");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgorithmParameter("param1") == "1.1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetSummaryParameter("cparam1") == "1.2.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetAlgVersion() == "1.0");
	CPPUNIT_ASSERT( chp2.GetHeader().GetParentCellFile() == "TestChip.cel");
	CPPUNIT_ASSERT( chp2.GetHeader().GetProgID() == "prog_id");

	CResequencingResults *b = chp2.GetResequencingResults();

	CPPUNIT_ASSERT(b->GetCalledBasesSize() == 2);
	CPPUNIT_ASSERT(b->GetScoresSize() == 2);
	CPPUNIT_ASSERT(b->GetForceCallsSize() == 2);
	CPPUNIT_ASSERT(b->GetOrigCallsSize() == 2);

	CPPUNIT_ASSERT(b->GetCalledBase(0) == 'a');
	CPPUNIT_ASSERT(b->GetCalledBase(1) == 'c');

	CPPUNIT_ASSERT_DOUBLES_EQUAL(b->GetScore(0), 1.0f, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(b->GetScore(1), 2.0f, 0.0001f);

	f = b->GetForceCall(0);
	CPPUNIT_ASSERT(f.position == 0);
	CPPUNIT_ASSERT(f.call == 'g');
	CPPUNIT_ASSERT(f.reason == 'b');

	f = b->GetForceCall(1);
	CPPUNIT_ASSERT(f.position == 1);
	CPPUNIT_ASSERT(f.call == 't');
	CPPUNIT_ASSERT(f.reason == 'c');

	o = b->GetOrigCall(0);
	CPPUNIT_ASSERT(o.position == 0);
	CPPUNIT_ASSERT(o.call == 't');

	o = b->GetOrigCall(1);
	CPPUNIT_ASSERT(o.position == 1);
	CPPUNIT_ASSERT(o.call == 'g');
}
