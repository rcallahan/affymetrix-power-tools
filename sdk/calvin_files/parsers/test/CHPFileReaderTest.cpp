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
#include "calvin_files/parsers/test/CHPFileReaderTest.h"
//
#include "calvin_files/data/src/CHPBackgroundZone.h"
#include "calvin_files/data/src/CHPExpressionEntry.h"
#include "calvin_files/data/src/CHPReseqEntry.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/AffymetrixParameterConsts.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
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

CPPUNIT_TEST_SUITE_REGISTRATION( CHPFileReaderTest );

void CHPFileReaderTest::setUp()
{
}

void CHPFileReaderTest::tearDown()
{
}

void CHPFileReaderTest::CreationTest()
{
	CHPFileReader reader;
	CPPUNIT_ASSERT(1);
}

void CHPFileReaderTest::ReadCHPExpressionFileTest()
{
	CHPData data;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename("../data/CHP_expression_file"));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == "../data/CHP_expression_file");
	CPPUNIT_ASSERT(data.GetEntryCount() == 10);
	CPPUNIT_ASSERT(data.GetBackgroundZoneCnt() == 2);

	CHPExpressionEntry entry1;
	data.GetEntry(0, entry1);
	CPPUNIT_ASSERT(entry1.GetProbeSetName() == "probe set 1");
	CPPUNIT_ASSERT(entry1.GetDetection() == 10);
	CPPUNIT_ASSERT(entry1.GetDetectionPValue() == 11.0f);
	float signal = entry1.GetSignal();
	CPPUNIT_ASSERT(signal == 17.8f);
	CPPUNIT_ASSERT(entry1.GetNumPairs() == 6);
	CPPUNIT_ASSERT(entry1.GetNumPairsUsed() == 5);
	CPPUNIT_ASSERT(entry1.GetChange() == 2);
	float chgPVal = entry1.GetChangePValue();
	CPPUNIT_ASSERT(chgPVal == 56.9f);
	CPPUNIT_ASSERT(entry1.GetSigLogRatio() == 45.89f);
	CPPUNIT_ASSERT(entry1.GetSigLogRatioLo() == 42.0f);
	CPPUNIT_ASSERT(entry1.GetSigLogRatioHi() == 47.0f);
	CPPUNIT_ASSERT(entry1.GetCommonPairs() == 2);

	CHPExpressionEntry entry2;
	data.GetEntry(1, entry2);
	CPPUNIT_ASSERT(entry2.GetProbeSetName() == "probe set 2");
	CPPUNIT_ASSERT(entry2.GetDetection() == 10);
	CPPUNIT_ASSERT(entry2.GetDetectionPValue() == 1.0f);
	CPPUNIT_ASSERT(entry2.GetSignal() == 7.8f);
	CPPUNIT_ASSERT(entry2.GetNumPairs() == 6);
	CPPUNIT_ASSERT(entry2.GetNumPairsUsed() == 5);
	CPPUNIT_ASSERT(entry2.GetChange() == 2);
	CPPUNIT_ASSERT(entry2.GetChangePValue() == 5.9f);
	CPPUNIT_ASSERT(entry2.GetSigLogRatio() == 4.89f);
	CPPUNIT_ASSERT(entry2.GetSigLogRatioLo() == 2.0f);
	CPPUNIT_ASSERT(entry2.GetSigLogRatioHi() == 7.0f);
	CPPUNIT_ASSERT(entry2.GetCommonPairs() == 2);

	CHPBackgroundZone zone1;
	data.GetBackgroundZone(0, zone1);
	CPPUNIT_ASSERT(zone1.GetCenterX() == 11.0f);
	CPPUNIT_ASSERT(zone1.GetCenterY() == 17.8f);
	CPPUNIT_ASSERT(zone1.GetBackground() == 56.9f);
	CPPUNIT_ASSERT(zone1.GetSmoothFactor() == 45.89f);

	CHPBackgroundZone zone2;
	data.GetBackgroundZone(1, zone2);
	CPPUNIT_ASSERT(zone2.GetCenterX() == 1.0f);
	CPPUNIT_ASSERT(zone2.GetCenterY() == 7.8f);
	CPPUNIT_ASSERT(zone2.GetBackground() == 6.9f);
	CPPUNIT_ASSERT(zone2.GetSmoothFactor() == 5.89f);
}

void CHPFileReaderTest::ReadCHPGenotypingFileTest()
{
	CHPData data;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename("../data/CHP_genotype_file"));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == "../data/CHP_genotype_file");
	CPPUNIT_ASSERT(data.GetEntryCount() == 2);
	CPPUNIT_ASSERT(data.GetBackgroundZoneCnt() == 2);

	CHPGenotypeEntry entry1;
	data.GetEntry(0, entry1);
	CPPUNIT_ASSERT(entry1.GetProbeSetName() == "probe set 1");
	CPPUNIT_ASSERT(entry1.GetCall() == 1);
	CPPUNIT_ASSERT(entry1.GetConfidence() == 11.0f);
	float ras1 = entry1.GetRAS1();
	CPPUNIT_ASSERT(ras1 == 17.8f);
	CPPUNIT_ASSERT(entry1.GetRAS2() == 6.0f);
	CPPUNIT_ASSERT(entry1.GetAACall() == 5.0f);
	CPPUNIT_ASSERT(entry1.GetABCall() == 2.0f);
	CPPUNIT_ASSERT(entry1.GetBBCall() == 56.9f);
	CPPUNIT_ASSERT(entry1.GetNoCall() == 45.89f);

	CHPGenotypeEntry entry2;
	data.GetEntry(1, entry2);
	CPPUNIT_ASSERT(entry2.GetProbeSetName() == "probe set 2");
	CPPUNIT_ASSERT(entry2.GetCall() == 2);
	CPPUNIT_ASSERT(entry2.GetConfidence() == 1.0f);
	CPPUNIT_ASSERT(entry2.GetRAS1() == 7.8f);
	CPPUNIT_ASSERT(entry2.GetRAS2() == 6.0f);
	CPPUNIT_ASSERT(entry2.GetAACall() == 5.0f);
	CPPUNIT_ASSERT(entry2.GetABCall() == 2.0f);
	CPPUNIT_ASSERT(entry2.GetBBCall() == 5.9f);
	CPPUNIT_ASSERT(entry2.GetNoCall() == 4.8f);

	CHPBackgroundZone zone1;
	data.GetBackgroundZone(0, zone1);
	CPPUNIT_ASSERT(zone1.GetCenterX() == 11.0f);
	CPPUNIT_ASSERT(zone1.GetCenterY() == 17.8f);
	CPPUNIT_ASSERT(zone1.GetBackground() == 56.9f);
	CPPUNIT_ASSERT(zone1.GetSmoothFactor() == 45.89f);

	CHPBackgroundZone zone2;
	data.GetBackgroundZone(1, zone2);
	CPPUNIT_ASSERT(zone2.GetCenterX() == 1.0f);
	CPPUNIT_ASSERT(zone2.GetCenterY() == 7.8f);
	CPPUNIT_ASSERT(zone2.GetBackground() == 6.9f);
	CPPUNIT_ASSERT(zone2.GetSmoothFactor() == 5.89f);
}

void CHPFileReaderTest::ReadCHPUniversalFileTest()
{
	CHPData data;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename("../data/CHP_universal_file"));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == "../data/CHP_universal_file");
	CPPUNIT_ASSERT(data.GetEntryCount() == 2);
	CPPUNIT_ASSERT(data.GetBackgroundZoneCnt() == 2);

	CHPUniversalEntry entry1;
	data.GetEntry(0, entry1);
	CPPUNIT_ASSERT(entry1.GetBackground() == 11.0f);

	CHPUniversalEntry entry2;
	data.GetEntry(1, entry2);
	CPPUNIT_ASSERT(entry2.GetBackground() == 1.0f);

	CHPBackgroundZone zone1;
	data.GetBackgroundZone(0, zone1);
	CPPUNIT_ASSERT(zone1.GetCenterX() == 11.0f);
	CPPUNIT_ASSERT(zone1.GetCenterY() == 17.8f);
	CPPUNIT_ASSERT(zone1.GetBackground() == 56.9f);
	CPPUNIT_ASSERT(zone1.GetSmoothFactor() == 45.89f);

	CHPBackgroundZone zone2;
	data.GetBackgroundZone(1, zone2);
	CPPUNIT_ASSERT(zone2.GetCenterX() == 1.0f);
	CPPUNIT_ASSERT(zone2.GetCenterY() == 7.8f);
	CPPUNIT_ASSERT(zone2.GetBackground() == 6.9f);
	CPPUNIT_ASSERT(zone2.GetSmoothFactor() == 5.89f);
}

void CHPFileReaderTest::ReadCHPReseqFileTest()
{
	CHPData data;
	CHPFileReader reader;
	CPPUNIT_ASSERT_NO_THROW(reader.SetFilename("../data/CHP_reseq_file"));
	CPPUNIT_ASSERT_NO_THROW(reader.Read(data));

	CPPUNIT_ASSERT(data.GetFilename() == "../data/CHP_reseq_file");

	CPPUNIT_ASSERT(data.GetEntryCount() == 5);
	CPPUNIT_ASSERT(data.GetBackgroundZoneCnt() == 1);
	CPPUNIT_ASSERT(data.GetForceCnt() == 2);
	CPPUNIT_ASSERT(data.GetOrigCnt() == 3);

	const double eps = 1e-5;
	CHPReseqEntry e;
	data.GetEntry(0, e);
	CPPUNIT_ASSERT(e.call == 'a');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 1.0f, eps);

	data.GetEntry(1, e);
	CPPUNIT_ASSERT(e.call == 'c');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 2.0f, eps);

	data.GetEntry(2, e);
	CPPUNIT_ASSERT(e.call == 'g');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 3.0f, eps);

	data.GetEntry(3, e);
	CPPUNIT_ASSERT(e.call == 't');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 4.0f, eps);

	data.GetEntry(4, e);
	CPPUNIT_ASSERT(e.call == 'n');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(e.score, 5.0f, eps);

	CHPReseqForceCall force;
	data.GetForceCall(0, force);
	CPPUNIT_ASSERT(force.position == 1);
	CPPUNIT_ASSERT(force.call == 'a');
	CPPUNIT_ASSERT(force.reason == CC_SATURATION_LEVEL_FORCE_CALL);

	data.GetForceCall(1, force);
	CPPUNIT_ASSERT(force.position == 2);
	CPPUNIT_ASSERT(force.call == 'c');
	CPPUNIT_ASSERT(force.reason == CC_WEAK_SIGNAL_THR_FORCE_CALL);

	CHPReseqOrigCall orig;
	data.GetOrigCall(0, orig);
	CPPUNIT_ASSERT(orig.position == 3);
	CPPUNIT_ASSERT(orig.call == 't');

	data.GetOrigCall(1, orig);
	CPPUNIT_ASSERT(orig.position == 4);
	CPPUNIT_ASSERT(orig.call == 'a');

	data.GetOrigCall(2, orig);
	CPPUNIT_ASSERT(orig.position == 5);
	CPPUNIT_ASSERT(orig.call == 'g');

}
