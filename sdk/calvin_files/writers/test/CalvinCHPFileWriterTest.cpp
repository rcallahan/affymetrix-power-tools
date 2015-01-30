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

//
#include "calvin_files/writers/test/CalvinCHPFileWriterTest.h"
//
#include "calvin_files/data/src/CHPBackgroundZone.h"
#include "calvin_files/data/src/CHPExpressionEntry.h"
#include "calvin_files/data/src/CHPGenotypeEntry.h"
#include "calvin_files/data/src/CHPReseqEntry.h"
#include "calvin_files/data/src/CHPUniversalEntry.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/writers/src/CalvinCHPFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CHPFileWriterTest );

void CHPFileWriterTest::setUp() {}

void CHPFileWriterTest::tearDown(){}

void CHPFileWriterTest::testCreation()
{
	CHPData fHdr("CHP_file", CHP_EXPRESSION_ASSAY_TYPE);
	CHPFileWriter* w = new CHPFileWriter(fHdr);
	CPPUNIT_ASSERT(1);
	delete w;
}

void CHPFileWriterTest::WriteExpressionTest()
{
	CHPData data("CHP_expression_file", CHP_EXPRESSION_ASSAY_TYPE);
	data.SetEntryCount(2, 16, true);
	data.SetBackgroundZoneCnt(2);
	writer = new CHPFileWriter(data);

	CHPExpressionEntry e1("probe set 1",10,11.0f,17.8f,6,5,true,2,56.9f,45.89f,42.0f,47.0f,2);
	CHPExpressionEntry e2("probe set 2",10,1.0f,7.8f,6,5,true,2,5.9f,4.89f,2.0f,7.0f,2);
	writer->SeekToDataSet();
	writer->WriteExpressionEntry(e1);
	writer->WriteExpressionEntry(e2);

	CHPBackgroundZone bgz1(11.0f,17.8f,56.9f,45.89f);
	CHPBackgroundZone bgz2(1.0f,7.8f,6.9f,5.89f);
	writer->SeekToBgSet();
	writer->WriteBackgroundZone(bgz1);
	writer->WriteBackgroundZone(bgz2);
	CPPUNIT_ASSERT(1);

	delete writer;


	CHPFileReader reader;
	reader.SetFilename("CHP_expression_file");
	CHPData data2;
	reader.Read(data2);

	CHPExpressionEntry e;
	data2.GetEntry(0, e);
	CPPUNIT_ASSERT(e.GetProbeSetName() == e1.GetProbeSetName());
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetSignal(), e1.GetSignal(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetDetectionPValue(), e1.GetDetectionPValue(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetChangePValue(), e1.GetChangePValue(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetSigLogRatio(), e1.GetSigLogRatio(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetSigLogRatioLo(), e1.GetSigLogRatioLo(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetSigLogRatioHi(), e1.GetSigLogRatioHi(), 0.000001f);
	CPPUNIT_ASSERT(e.GetDetection() ==e1.GetDetection());
	CPPUNIT_ASSERT(e.GetNumPairs() ==e1.GetNumPairs());
	CPPUNIT_ASSERT(e.GetNumPairsUsed() ==e1.GetNumPairsUsed());
	CPPUNIT_ASSERT(e.GetHasComparisonData() ==e1.GetHasComparisonData());
	CPPUNIT_ASSERT(e.GetChange() ==e1.GetChange());
	CPPUNIT_ASSERT(e.GetCommonPairs() ==e1.GetCommonPairs());

	data2.GetEntry(1, e);
	CPPUNIT_ASSERT(e.GetProbeSetName() == e2.GetProbeSetName());
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetSignal(), e2.GetSignal(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetDetectionPValue(), e2.GetDetectionPValue(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetChangePValue(), e2.GetChangePValue(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetSigLogRatio(), e2.GetSigLogRatio(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetSigLogRatioLo(), e2.GetSigLogRatioLo(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetSigLogRatioHi(), e2.GetSigLogRatioHi(), 0.000001f);
	CPPUNIT_ASSERT(e.GetDetection() ==e2.GetDetection());
	CPPUNIT_ASSERT(e.GetNumPairs() ==e2.GetNumPairs());
	CPPUNIT_ASSERT(e.GetNumPairsUsed() ==e2.GetNumPairsUsed());
	CPPUNIT_ASSERT(e.GetHasComparisonData() ==e2.GetHasComparisonData());
	CPPUNIT_ASSERT(e.GetChange() ==e2.GetChange());
	CPPUNIT_ASSERT(e.GetCommonPairs() ==e2.GetCommonPairs());


	CHPBackgroundZone bgz;
	data2.GetBackgroundZone(0, bgz);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetBackground(), bgz1.GetBackground(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterX(), bgz1.GetCenterX(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterY(), bgz1.GetCenterY(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetSmoothFactor(), bgz1.GetSmoothFactor(), 0.000001f);

	data2.GetBackgroundZone(1, bgz);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetBackground(), bgz2.GetBackground(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterX(), bgz2.GetCenterX(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterY(), bgz2.GetCenterY(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetSmoothFactor(), bgz2.GetSmoothFactor(), 0.000001f);

}

void CHPFileWriterTest::WriteGenotypeTest()
{
	CHPData data("CHP_genotype_file", CHP_GENOTYPING_ASSAY_TYPE);
	data.SetEntryCount(2, 16);
	data.SetBackgroundZoneCnt(2);
	writer = new CHPFileWriter(data);

	CHPGenotypeEntry e1("probe set 1",1,11.0f,17.8f,6.0f,5.0f,2.0f,56.9f,45.89f);
	CHPGenotypeEntry e2("probe set 2",2,1.0f,7.8f,6.0f,5.0f,2.0f,5.9f,4.8f);
	writer->SeekToDataSet();
	writer->WriteGenotypeEntry(e1);
	writer->WriteGenotypeEntry(e2);

	CHPBackgroundZone bgz1(11.0f,17.8f,56.9f,45.89f);
	CHPBackgroundZone bgz2(1.0f,7.8f,6.9f,5.89f);
	writer->SeekToBgSet();
	writer->WriteBackgroundZone(bgz1);
	writer->WriteBackgroundZone(bgz2);
	CPPUNIT_ASSERT(1);

	delete writer;

	CHPFileReader reader;
	reader.SetFilename("CHP_genotype_file");
	CHPData data2;
	reader.Read(data2);

	CHPGenotypeEntry e;
	data2.GetEntry(0, e);
	CPPUNIT_ASSERT(e.GetProbeSetName() == e1.GetProbeSetName());
	CPPUNIT_ASSERT(e.GetCall() == e1.GetCall());
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetConfidence(), e1.GetConfidence(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetRAS1(), e1.GetRAS1(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetRAS2(), e1.GetRAS2(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetAACall(), e1.GetAACall(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetABCall(), e1.GetABCall(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetBBCall(), e1.GetBBCall(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetNoCall(), e1.GetNoCall(), 0.00001f);


	data2.GetEntry(1, e);
	CPPUNIT_ASSERT(e.GetProbeSetName() == e2.GetProbeSetName());
	CPPUNIT_ASSERT(e.GetCall() == e2.GetCall());
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetConfidence(), e2.GetConfidence(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetRAS1(), e2.GetRAS1(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetRAS2(), e2.GetRAS2(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetAACall(), e2.GetAACall(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetABCall(), e2.GetABCall(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetBBCall(), e2.GetBBCall(), 0.00001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetNoCall(), e2.GetNoCall(), 0.00001f);


	CHPBackgroundZone bgz;
	data2.GetBackgroundZone(0, bgz);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetBackground(), bgz1.GetBackground(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterX(), bgz1.GetCenterX(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterY(), bgz1.GetCenterY(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetSmoothFactor(), bgz1.GetSmoothFactor(), 0.000001f);

	data2.GetBackgroundZone(1, bgz);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetBackground(), bgz2.GetBackground(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterX(), bgz2.GetCenterX(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterY(), bgz2.GetCenterY(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetSmoothFactor(), bgz2.GetSmoothFactor(), 0.000001f);

}

void CHPFileWriterTest::WriteUniversalTest()
{
	CHPData data("CHP_universal_file", CHP_UNIVERSAL_ASSAY_TYPE);
	data.SetEntryCount(2, 8);
	data.SetBackgroundZoneCnt(2);
	writer = new CHPFileWriter(data);

	CHPUniversalEntry e1(11.0f);
	CHPUniversalEntry e2(1.0f);
	writer->SeekToDataSet();
	writer->WriteUniversalEntry(e1);
	writer->WriteUniversalEntry(e2);

	CHPBackgroundZone bgz1(11.0f,17.8f,56.9f,45.89f);
	CHPBackgroundZone bgz2(1.0,7.8f,6.9f,5.89f);
	writer->SeekToBgSet();
	writer->WriteBackgroundZone(bgz1);
	writer->WriteBackgroundZone(bgz2);
	CPPUNIT_ASSERT(1);

	delete writer;




	CHPFileReader reader;
	reader.SetFilename("CHP_universal_file");
	CHPData data2;
	reader.Read(data2);

	CHPUniversalEntry e;

	data2.GetEntry(0, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetBackground(), e1.GetBackground(), 0.000001f);
	data2.GetEntry(1, e);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.GetBackground(), e2.GetBackground(), 0.000001f);
	
	CHPBackgroundZone bgz;
	data2.GetBackgroundZone(0, bgz);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetBackground(), bgz1.GetBackground(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterX(), bgz1.GetCenterX(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterY(), bgz1.GetCenterY(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetSmoothFactor(), bgz1.GetSmoothFactor(), 0.000001f);

	data2.GetBackgroundZone(1, bgz);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetBackground(), bgz2.GetBackground(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterX(), bgz2.GetCenterX(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetCenterY(), bgz2.GetCenterY(), 0.000001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(bgz.GetSmoothFactor(), bgz2.GetSmoothFactor(), 0.000001f);
}

void CHPFileWriterTest::WriteReseqTest()
{
	CHPData data("CHP_reseq_file", CHP_RESEQUENCING_ASSAY_TYPE);
	data.SetEntryCount(5, 10);
	data.SetBackgroundZoneCnt(1);
	data.SetForceCnt(2);
	data.SetOrigCnt(3);
	writer = new CHPFileWriter(data);
	writer->SeekToDataSet();
	
	CHPReseqEntry e;
	e.call = 'a';
	e.score = 1.0f;
	writer->WriteReseqEntry(e);
	e.call = 'c';
	e.score = 2.0f;
	writer->WriteReseqEntry(e);
	e.call = 'g';
	e.score = 3.0f;
	writer->WriteReseqEntry(e);
	e.call = 't';
	e.score = 4.0f;
	writer->WriteReseqEntry(e);
	e.call = 'n';
	e.score = 5.0f;
	writer->WriteReseqEntry(e);

	writer->SeekToForceSet();
	CHPReseqForceCall force;
	force.position = 1;
	force.call = 'a';
	force.reason = CC_SATURATION_LEVEL_FORCE_CALL;
	writer->WriteForceCall(force);
	force.position = 2;
	force.call = 'c';
	force.reason = CC_WEAK_SIGNAL_THR_FORCE_CALL;
	writer->WriteForceCall(force);

	writer->SeekToOrigCallSet();
	CHPReseqOrigCall orig;
	orig.position = 3;
	orig.call = 't';
	writer->WriteOrigCall(orig);
	orig.position = 4;
	orig.call = 'a';
	writer->WriteOrigCall(orig);
	orig.position = 5;
	orig.call = 'g';
	writer->WriteOrigCall(orig);

	CPPUNIT_ASSERT(1);

	delete writer;




	CHPFileReader reader;
	reader.SetFilename("CHP_reseq_file");
	CHPData data2;
	reader.Read(data2);

	data2.GetEntry(0, e);
	CPPUNIT_ASSERT(e.call == 'a');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 1.0f, 0.0001f);

	data2.GetEntry(1, e);
	CPPUNIT_ASSERT(e.call == 'c');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 2.0f, 0.0001f);

	data2.GetEntry(2, e);
	CPPUNIT_ASSERT(e.call == 'g');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 3.0f, 0.0001f);

	data2.GetEntry(3, e);
	CPPUNIT_ASSERT(e.call == 't');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 4.0f, 0.0001f);

	data2.GetEntry(4, e);
	CPPUNIT_ASSERT(e.call == 'n');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 5.0f, 0.0001f);


	data2.GetForceCall(0, force);
	CPPUNIT_ASSERT(force.call == 'a');
	CPPUNIT_ASSERT(force.reason == CC_SATURATION_LEVEL_FORCE_CALL);
	CPPUNIT_ASSERT(force.position == 1);

	data2.GetForceCall(1, force);
	CPPUNIT_ASSERT(force.call == 'c');
	CPPUNIT_ASSERT(force.reason == CC_WEAK_SIGNAL_THR_FORCE_CALL);
	CPPUNIT_ASSERT(force.position == 2);


	data2.GetOrigCall(0, orig);
	CPPUNIT_ASSERT(orig.call == 't');
	CPPUNIT_ASSERT(orig.position == 3);

	data2.GetOrigCall(1, orig);
	CPPUNIT_ASSERT(orig.call == 'a');
	CPPUNIT_ASSERT(orig.position == 4);

	data2.GetOrigCall(2, orig);
	CPPUNIT_ASSERT(orig.call == 'g');
	CPPUNIT_ASSERT(orig.position == 5);

}
