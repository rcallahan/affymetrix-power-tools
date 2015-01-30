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

#include "exp_report/test/ExpressionProbeSetReporterTest.h"
//
#include "exp_report/src/ExpressionProbeSetReporter.h"
//

using namespace ExpressionReport;

CPPUNIT_TEST_SUITE_REGISTRATION( ExpressionProbeSetReporterTest );

////////////////////////////////////////////////////////////////////////////////////////////////

void ExpressionProbeSetReporterTest::setUp()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////

void ExpressionProbeSetReporterTest::tearDown()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////

class TestReportDataAccessor : public ReportDataAccessor
{
public:
	int GetNumProbeSets() { return 10; }
	int GetNumPairs(int index)
	{
		return index+5; // 5,6,7,8, 9,10,11,12,13,14
	}
	DetectionCall GetDetection(int index) {
		if (index>=9)
			return DetectionPresent;
		else if (index>=0 && index<=7)
			return DetectionAbsent;
		else
			return DetectionMarginal;
	}
	bool IsAntiSense(int index) { return (index>1); }
	float GetSignal(int index) { return (float)(index+1); }

	DetectionCall GetBaselineDetection(int index) { return DetectionNoCall; }
	bool HasComparisonData() { return false; }
	ChangeCall GetChange(int index) { return ChangeNoCall; }
	float GetSignalLogRatio(int index) { return 0.0f; }
        bool IsLogScale() { return false; }
    std::vector<float> GetIntensities(QCProbeSetType qctype) { std::vector<float> v; v.resize(2); v[0] = 1.0f; v[1] = 2.0f; return v; }
    std::string GetProbeSetName(int index) { return "affx1"; }
};

////////////////////////////////////////////////////////////////////////////////////////////////

void ExpressionProbeSetReporterTest::testDetection()
{
	ExpressionProbeSetReporter stats;
	TestReportDataAccessor data;
    std::list<int> probeSetIndicies;
	CPPUNIT_ASSERT(stats.Run(true, 9, &data, NULL, probeSetIndicies) == true);

	CPPUNIT_ASSERT(stats.Results().SpikeStats().size() == 0);
	CPPUNIT_ASSERT(stats.Results().HousekeepingStats().size() == 0);
	CPPUNIT_ASSERT(stats.Results().ProbeSetResults().NumSets() == 6);
	CPPUNIT_ASSERT(stats.Results().ProbeSetResults().PresentCalls().Count() == 1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stats.Results().ProbeSetResults().PresentCalls().Signal(), 10.0, 0.0001);
	CPPUNIT_ASSERT(stats.Results().ProbeSetResults().MarginalCalls().Count() == 1);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stats.Results().ProbeSetResults().MarginalCalls().Signal(), 9.0, 0.0001);
	CPPUNIT_ASSERT(stats.Results().ProbeSetResults().AbsentCalls().Count() == 4);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stats.Results().ProbeSetResults().AbsentCalls().Signal(), 26.0, 0.0001);

	stats.Results().Clear();
	CPPUNIT_ASSERT(stats.Results().ProbeSetResults().NumSets() == 0);
	CPPUNIT_ASSERT(stats.Results().ProbeSetResults().PresentCalls().Count() == 0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stats.Results().ProbeSetResults().PresentCalls().Signal(), 0.0, 0.0001);
	CPPUNIT_ASSERT(stats.Results().ProbeSetResults().MarginalCalls().Count() == 0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stats.Results().ProbeSetResults().MarginalCalls().Signal(), 0.0, 0.0001);
	CPPUNIT_ASSERT(stats.Results().ProbeSetResults().AbsentCalls().Count() == 0);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(stats.Results().ProbeSetResults().AbsentCalls().Signal(), 0.0, 0.0001);
}

////////////////////////////////////////////////////////////////////////////////////////////////

void ExpressionProbeSetReporterTest::testControls()
{
	ExpressionProbeSetReporter stats;
	TestReportDataAccessor data;
    std::list<int> probeSetIndicies;

	ExpressionControls cs;
	ExpressionControl c;
	cs.ProbeArrayType() = "test3";
	c.Name() = "s1";
	c.SetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET, 1);
	c.SetProbeSetIndex(ExpressionControl::MIDDLE_PROBE_SET, 2);
	c.SetProbeSetIndex(ExpressionControl::FIVE_PRIME_PROBE_SET, 3);
	cs.SpikeControls().push_back(c);

	c.Clear();
	c.Name() = "s2";
	c.SetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET, 3);
	cs.SpikeControls().push_back(c);

	CPPUNIT_ASSERT(stats.Run(true, 9, &data, &cs, probeSetIndicies) == true);

	CPPUNIT_ASSERT(stats.Results().SpikeStats().size() == 2);

	ExpressionControlResult r;
	ExpressionControlResultList::const_iterator it = stats.Results().SpikeStats().begin();
	r = *it;
	CPPUNIT_ASSERT(r.GetName() == "s1");
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::THREE_PRIME_PROBE_SET) == true);
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::MIDDLE_PROBE_SET) == true);
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == true);

	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionAbsent);

	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET), 2, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET), 3, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET), 4, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetThreeFiveRatio(), 0.5, 0.001);

	++it;
	r = *it;
	CPPUNIT_ASSERT(r.GetName() == "s2");
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::THREE_PRIME_PROBE_SET) == true);
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::MIDDLE_PROBE_SET) == false);
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == false);

	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET) == ReportDataAccessor::DetectionNoCall);
	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionNoCall);

	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET), 4, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET), -1, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET), -1, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetThreeFiveRatio(), -1.0, 0.001);

	CPPUNIT_ASSERT(stats.Results().HousekeepingStats().size() == 0);




	cs.Clear();
	cs.ProbeArrayType() = "test3";
	c.Name() = "s1";
	c.SetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET, 1);
	c.SetProbeSetIndex(ExpressionControl::MIDDLE_PROBE_SET, 2);
	c.SetProbeSetIndex(ExpressionControl::FIVE_PRIME_PROBE_SET, 3);
	cs.HousekeepingControls().push_back(c);

	c.Clear();
	c.Name() = "s2";
	c.SetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET, 3);
	cs.HousekeepingControls().push_back(c);

	stats.Results().Clear();
	CPPUNIT_ASSERT(stats.Run(true, 9, &data, &cs, probeSetIndicies) == true);

	CPPUNIT_ASSERT(stats.Results().HousekeepingStats().size() == 2);

	it = stats.Results().HousekeepingStats().begin();
	r = *it;
	CPPUNIT_ASSERT(r.GetName() == "s1");
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::THREE_PRIME_PROBE_SET) == true);
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::MIDDLE_PROBE_SET) == true);
	CPPUNIT_ASSERT(r.HasControlResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == true);

	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET) == ReportDataAccessor::DetectionAbsent);
	CPPUNIT_ASSERT(r.GetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == ReportDataAccessor::DetectionAbsent);

	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET), 2, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET), 3, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET), 4, 0.001);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(r.GetThreeFiveRatio(), 0.5, 0.001);

	CPPUNIT_ASSERT(stats.Results().SpikeStats().size() == 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////

void ExpressionProbeSetReporterTest::testControlProbeSet()
{
    
	ExpressionProbeSetReporter stats;
	TestReportDataAccessor data;
	ExpressionControls cs;
    std::list<int> probeSetIndicies;

    CPPUNIT_ASSERT(stats.Run(true, 9, &data, &cs, probeSetIndicies) == true);

    CPPUNIT_ASSERT(stats.Results().ControlStats().size() == 4);

    NameAvgCountList::iterator it;
    it = stats.Results().ControlStats().begin();

    CPPUNIT_ASSERT_DOUBLES_EQUAL(it->avg, 1.5, 0.001);
    CPPUNIT_ASSERT(it->count == 2);

    ++it;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(it->avg, 1.5, 0.001);
    CPPUNIT_ASSERT(it->count == 2);

    ++it;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(it->avg, 1.5, 0.001);
    CPPUNIT_ASSERT(it->count == 2);

    ++it;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(it->avg, 1.5, 0.001);
    CPPUNIT_ASSERT(it->count == 2);
}

////////////////////////////////////////////////////////////////////////////////////////////////

void ExpressionProbeSetReporterTest::testProbeSetSignals()
{
	ExpressionProbeSetReporter stats;
	TestReportDataAccessor data;
	ExpressionControls cs;
    std::list<int> probeSetIndicies;
    probeSetIndicies.push_back(0);

    CPPUNIT_ASSERT(stats.Run(true, 9, &data, &cs, probeSetIndicies) == true);

    CPPUNIT_ASSERT(stats.Results().ProbeSetValues().size() == 1);
    NameFloatValuePair param = *stats.Results().ProbeSetValues().begin();
    CPPUNIT_ASSERT(param.name == "affx1");
    CPPUNIT_ASSERT_DOUBLES_EQUAL(param.value, 1.0, 0.001);
}

////////////////////////////////////////////////////////////////////////////////////////////////

class TestReportDataAccessor2 : public ReportDataAccessor
{
public:
	int GetNumProbeSets() { return 14; }
	int GetNumPairs(int index) { return 20; }
	DetectionCall GetDetection(int index) {
		switch (index)
		{
		case 0:
		case 1:
		case 3:
		case 6:
		case 7:
		case 9:
		case 13:
			return DetectionPresent;
			break;

		case 2:
		case 4:
		case 5:
		case 8:
		case 10:
		case 11:
		case 12:
			return DetectionAbsent;
			break;

		default:
			return DetectionMarginal;
			break;
		}
	}
	bool IsAntiSense(int index) { return true; }
	float GetSignal(int index) { return (float)(index+1); }

	bool HasComparisonData() { return true; }
	DetectionCall GetBaselineDetection(int index) {
		switch (index)
		{
		case 0:
		case 3:
		case 4:
		case 6:
		case 9:
		case 10:
		case 13:
			return DetectionPresent;
			break;

		case 1:
		case 2:
		case 5:
		case 7:
		case 8:
		case 11:
		case 12:
			return DetectionAbsent;
			break;

		default:
			return DetectionMarginal;
			break;
		}
	}

	ChangeCall GetChange(int index)
	{
		if (index < 3)
			return ChangeIncrease;
		else if (index >= 3 && index < 6)
			return ChangeDecrease;
		else if (index >= 6 && index < 9)
			return ChangeModerateIncrease;
		else if (index >= 9 && index < 12)
			return ChangeModerateDecrease;
		else if (index >= 12 && index < 14)
			return ChangeNoChange;
		else
			return ChangeNoCall;
	}
	float GetSignalLogRatio(int index) {
		float val = (index+1)/2.0f;
		if (GetChange(index) == ChangeIncrease || GetChange(index) == ChangeModerateIncrease)
			return val;
		else
			return -val;
	}
        bool IsLogScale() { return false; }
    std::vector<float> GetIntensities(QCProbeSetType qctype) { std::vector<float> v; v.resize(2); v[0] = 1.0f; v[1] = 2.0f; return v; }
    std::string GetProbeSetName(int index) { return ""; }
};

////////////////////////////////////////////////////////////////////////////////////////////////

void ExpressionProbeSetReporterTest::testChange()
{
	ExpressionProbeSetReporter stats;
	TestReportDataAccessor2 data;
    std::list<int> probeSetIndicies;
	CPPUNIT_ASSERT(stats.Run(true, 9, &data, NULL, probeSetIndicies) == true);

	CPPUNIT_ASSERT(stats.Results().IncreaseStats().DetectionPresentCount() == 2);
	CPPUNIT_ASSERT(stats.Results().IncreaseStats().DetectionAbsentCount() == 2);
	CPPUNIT_ASSERT(stats.Results().IncreaseStats().DetectionChangeCount() == 2);
	CPPUNIT_ASSERT(stats.Results().IncreaseStats().ChangeCount() == 3);
	CPPUNIT_ASSERT(stats.Results().IncreaseStats().ModerateCount() == 3);

	CPPUNIT_ASSERT(stats.Results().IncreaseStats().FoldChangeCount(0) == 1);
	CPPUNIT_ASSERT(stats.Results().IncreaseStats().FoldChangeCount(1) == 2);
	CPPUNIT_ASSERT(stats.Results().IncreaseStats().FoldChangeCount(2) == 0);
	CPPUNIT_ASSERT(stats.Results().IncreaseStats().FoldChangeCount(3) == 1);
	CPPUNIT_ASSERT(stats.Results().IncreaseStats().FoldChangeCount(4) == 2);

	CPPUNIT_ASSERT(stats.Results().DecreaseStats().DetectionPresentCount() == 2);
	CPPUNIT_ASSERT(stats.Results().DecreaseStats().DetectionAbsentCount() == 2);
	CPPUNIT_ASSERT(stats.Results().DecreaseStats().DetectionChangeCount() == 2);
	CPPUNIT_ASSERT(stats.Results().DecreaseStats().ChangeCount() == 3);
	CPPUNIT_ASSERT(stats.Results().DecreaseStats().ModerateCount() == 3);

	CPPUNIT_ASSERT(stats.Results().DecreaseStats().FoldChangeCount(0) == 0);
	CPPUNIT_ASSERT(stats.Results().DecreaseStats().FoldChangeCount(1) == 0);
	CPPUNIT_ASSERT(stats.Results().DecreaseStats().FoldChangeCount(2) == 2);
	CPPUNIT_ASSERT(stats.Results().DecreaseStats().FoldChangeCount(3) == 1);
	CPPUNIT_ASSERT(stats.Results().DecreaseStats().FoldChangeCount(4) == 3);

	CPPUNIT_ASSERT(stats.Results().NoChangeStats().ChangeCount() == 2);
	CPPUNIT_ASSERT(stats.Results().NoChangeStats().DetectionAbsentCount() == 1);
	CPPUNIT_ASSERT(stats.Results().NoChangeStats().DetectionChangeCount() == 0);
	CPPUNIT_ASSERT(stats.Results().NoChangeStats().DetectionPresentCount() == 1);
	CPPUNIT_ASSERT(stats.Results().NoChangeStats().FoldChangeCount(0) == 0);
	CPPUNIT_ASSERT(stats.Results().NoChangeStats().FoldChangeCount(1) == 0);
	CPPUNIT_ASSERT(stats.Results().NoChangeStats().FoldChangeCount(2) == 0);
	CPPUNIT_ASSERT(stats.Results().NoChangeStats().FoldChangeCount(3) == 0);
	CPPUNIT_ASSERT(stats.Results().NoChangeStats().FoldChangeCount(4) == 0);

    
    CPPUNIT_ASSERT(stats.Results().ProbeSetValues().size() == 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////
