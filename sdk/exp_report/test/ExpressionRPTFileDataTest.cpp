////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#include "exp_report/test/ExpressionRPTFileDataTest.h"
//
#include "exp_report/src/ExpressionRPTFileData.h"
//

using namespace ExpressionReport;
using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( ExpressionRPTFileDataTest );

void ExpressionRPTFileDataTest::setUp()
{
}

void ExpressionRPTFileDataTest::tearDown()
{
}

void ExpressionRPTFileDataTest::testRead()
{
	ExpressionRPTFileData reader;
	ExpressionReportData data;
	string file = "../data/test.rpt";
	reader.FileName() = file;
	CPPUNIT_ASSERT(reader.Read(data) == true);

	CPPUNIT_ASSERT(data.Date() == "08:28AM 05/23/2005");
	CPPUNIT_ASSERT(data.CHPFileName() == "hta_test_data.B12.CHP");
	CPPUNIT_ASSERT(data.ArrayType() == "U133AAofAv2");
	CPPUNIT_ASSERT(data.AlgName() == "Statistical");
	CPPUNIT_ASSERT(data.ProbePairThreshold() == 8);
	CPPUNIT_ASSERT(data.AntiSenseControls() == true);

	NameValuePairList::iterator it = data.AlgParams().begin();
	CPPUNIT_ASSERT(it->name == "Alpha1");
	CPPUNIT_ASSERT(it->value == "0.05");
	++it;
	CPPUNIT_ASSERT(it->name == "Alpha2");
	CPPUNIT_ASSERT(it->value == "0.065");
	++it;
	CPPUNIT_ASSERT(it->name == "Tau");
	CPPUNIT_ASSERT(it->value == "0.015");
	++it;
	CPPUNIT_ASSERT(it->name == "Noise (RawQ)");
	CPPUNIT_ASSERT(it->value == "2.370");
	++it;
	CPPUNIT_ASSERT(it->name == "Scale Factor (SF)");
	CPPUNIT_ASSERT(it->value == "0.879");
	++it;
	CPPUNIT_ASSERT(it->name == "TGT Value");
	CPPUNIT_ASSERT(it->value == "500");
	++it;
	CPPUNIT_ASSERT(it->name == "Norm Factor (NF)");
	CPPUNIT_ASSERT(it->value == "1.000");
	++it;
	CPPUNIT_ASSERT(it == data.AlgParams().end());

	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.BackgroundStats().avg, 333.23, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.BackgroundStats().std, 7.71, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.BackgroundStats().min, 314.60, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.BackgroundStats().max, 358.90, 0.0001f);

	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.NoiseStats().avg, 7.00, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.NoiseStats().std, 0.58, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.NoiseStats().min, 5.70, 0.0001f);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(data.NoiseStats().max, 9.40, 0.0001f);

	NameAvgCountList::iterator cit = data.ControlStats().begin();
	CPPUNIT_ASSERT(cit->name == "Corner+");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(cit->avg, 464, 0.01f);
	CPPUNIT_ASSERT(cit->count == 512);
	++cit;
	CPPUNIT_ASSERT(cit->name == "Corner-");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(cit->avg, 12208, 0.01f);
	CPPUNIT_ASSERT(cit->count == 512);
	++cit;
	CPPUNIT_ASSERT(cit->name == "Central-");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(cit->avg, 16353, 0.01f);
	CPPUNIT_ASSERT(cit->count == 9);
	++cit;
	CPPUNIT_ASSERT(cit->name == "Central+");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(cit->avg, 6353, 0.01f);
	CPPUNIT_ASSERT(cit->count == 9);
	++cit;

	CPPUNIT_ASSERT(data.ProbeSetResults().NumSets() == 22944);

	CPPUNIT_ASSERT(data.ProbeSetResults().PresentCalls().Count() == 14015);
	float exp = 14015*1013.5f;
	float act = data.ProbeSetResults().PresentCalls().Signal();
	CPPUNIT_ASSERT_DOUBLES_EQUAL(act, exp, 0.01f);

	CPPUNIT_ASSERT(data.ProbeSetResults().MarginalCalls().Count() == 280);
	exp = 280*191.1f;
	act = data.ProbeSetResults().MarginalCalls().Signal();
	CPPUNIT_ASSERT_DOUBLES_EQUAL(act, exp, 0.01f);

	CPPUNIT_ASSERT(data.ProbeSetResults().AbsentCalls().Count() == 8649);
	exp = 8649*43.8f;
	act = data.ProbeSetResults().AbsentCalls().Signal();
	CPPUNIT_ASSERT_DOUBLES_EQUAL(act, exp, 0.01f);

	act = (data.ProbeSetResults().AbsentCalls().Signal() + data.ProbeSetResults().PresentCalls().Signal() + data.ProbeSetResults().MarginalCalls().Signal()) /
		(data.ProbeSetResults().AbsentCalls().Count() + data.ProbeSetResults().PresentCalls().Count() + data.ProbeSetResults().MarginalCalls().Count());
	exp = 637.9f;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(act, exp, 0.1f);



	ExpressionControlResultList::iterator sit = data.SpikeStats().begin();

	CPPUNIT_ASSERT(sit->GetName() == "AFFX-Dap");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sit->GetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET), 13.2f, 0.1f);
	CPPUNIT_ASSERT(sit->GetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionPresent);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sit->GetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET), 0.7f, 0.1f);
	CPPUNIT_ASSERT(sit->GetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET) == ReportDataAccessor::DetectionMarginal);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sit->GetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET), 2.5f, 0.1f);
	CPPUNIT_ASSERT(sit->GetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sit->GetThreeFiveRatio(), 0.19f, 0.001f);

	++sit;
	CPPUNIT_ASSERT(sit->GetName() == "AFFX-Lys");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sit->GetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET), 8.9f, 0.1f);
	CPPUNIT_ASSERT(sit->GetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionNoCall);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sit->GetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET), 21.6f, 0.1f);
	CPPUNIT_ASSERT(sit->GetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(sit->GetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET), 8.3f, 0.1f);
	CPPUNIT_ASSERT(sit->GetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(sit->GetThreeFiveRatio(), 0.93f, 0.001f);

	++sit;
	CPPUNIT_ASSERT(sit == data.SpikeStats().end());


	ExpressionControlResultList::iterator hit = data.HousekeepingStats().begin();

	CPPUNIT_ASSERT(hit->GetName() == "AFFX-HUMISGF3A/M97935");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(hit->GetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET), 764.2f, 0.1f);
	CPPUNIT_ASSERT(hit->GetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionPresent);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(hit->GetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET), 1479.7f, 0.1f);
	CPPUNIT_ASSERT(hit->GetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET) == ReportDataAccessor::DetectionPresent);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(hit->GetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET), 3005.5f, 0.1f);
	CPPUNIT_ASSERT(hit->GetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionPresent);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(hit->GetThreeFiveRatio(), 3.93f, 0.001f);

	++hit;
	CPPUNIT_ASSERT(hit->GetName() == "AFFX-HUMRGE/M10098");
	CPPUNIT_ASSERT_DOUBLES_EQUAL(hit->GetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET), 75.2f, 0.1f);
	CPPUNIT_ASSERT(hit->GetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(hit->GetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET), 99.3f, 0.1f);
	CPPUNIT_ASSERT(hit->GetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET) == ReportDataAccessor::DetectionMarginal);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(hit->GetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET), 115.5f, 0.1f);
	CPPUNIT_ASSERT(hit->GetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionPresent);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(hit->GetThreeFiveRatio(), 1.54f, 0.001f);

	++hit;
	CPPUNIT_ASSERT(hit == data.HousekeepingStats().end());
}
