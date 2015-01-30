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

#include "exp_report/src/ExpressionProbeSetReporter.h"
//
#include "file/CDFFileData.h"
//
#include <cmath>
#include <cstdio>
//

#ifndef max
#define max(a, b)  (((a) > (b)) ? (a) : (b)) 
#endif

using namespace ExpressionReport;


/*
 * Clear the members.
 */
ExpressionProbeSetReporter::ExpressionProbeSetReporter()
{
	Clear();
}

/*
 * Clear the members.
 */
void ExpressionProbeSetReporter::Clear()
{
	results.Clear();
	dataAccessor=NULL;
}

/*
 * Compute the probe set statistics.
 */
bool ExpressionProbeSetReporter::Run(bool anti, int probePairThr, ReportDataAccessor *data, ExpressionControls *controls, const std::list<int> &probeSetIndicies, bool bDifference, bool includeAllProbesets)
{
	Clear();
	results.AntiSenseControls() = anti;
	results.ProbePairThreshold() = probePairThr;
	dataAccessor = data;
	ComputeDetectionProbeSetStats(includeAllProbesets);
	ComputeChangeProbeSetStats(includeAllProbesets);
    ComputeQCControlProbeSets();
    AddProbeSetSignals(probeSetIndicies);
	if (controls != NULL)
	{
		ComputeControlStats(controls->SpikeControls(), results.SpikeStats(), bDifference);
		ComputeControlStats(controls->HousekeepingControls(), results.HousekeepingStats(), bDifference);
	}
	return true;
}

/*
 * Add specific probe set signal values to the list.
 */
void ExpressionProbeSetReporter::AddProbeSetSignals(const std::list<int> &probeSetIndicies)
{
    NameFloatValuePair val;
    for (std::list<int>::const_iterator it=probeSetIndicies.begin(); it!=probeSetIndicies.end(); ++it)
    {
        int iSet = *it;
        val.name = dataAccessor->GetProbeSetName(iSet);
        val.value = dataAccessor->GetSignal(iSet);
        results.ProbeSetValues().push_back(val);
    }
}

/*
 * Determine if probe set should be included.
 */
bool ExpressionProbeSetReporter::IncludeProbeSet(int index)
{
	if (dataAccessor->GetNumPairs(index) == 0)
		return false;

	if (dataAccessor->GetNumPairs(index) >= results.ProbePairThreshold() && dataAccessor->IsAntiSense(index) == results.AntiSenseControls())
		return true;

	return false;
}

/*
 * Compute the average of the vector.
 */
static float Average(std::vector<float> &v)
{
    double sum = 0.0;
    for (std::vector<float>::iterator it=v.begin(); it!=v.end(); it++)
        sum += *it;
    return (float)(sum / v.size());
}

/*
 * Compute the central and corner controls.
 */
void ExpressionProbeSetReporter::ComputeQCControlProbeSets()
{
    std::vector<float> intensities;
    NameAvgCount result;

    intensities = dataAccessor->GetIntensities(ReportDataAccessor::CheckerboardPositiveProbeSetType);
    if (intensities.size() > 0)
    {
        result.avg = Average(intensities);
        result.count = (int) intensities.size();
        result.name = "Raw Corner+";
        results.ControlStats().push_back(result);
    }

    intensities = dataAccessor->GetIntensities(ReportDataAccessor::CheckerboardNegativeProbeSetType);
    if (intensities.size() > 0)
    {
        result.avg = Average(intensities);
        result.count = (int) intensities.size();
        result.name = "Raw Corner-";
        results.ControlStats().push_back(result);
    }

    intensities = dataAccessor->GetIntensities(ReportDataAccessor::CentralCrossPositiveProbeSetType);
    if (intensities.size() > 0)
    {
        result.avg = Average(intensities);
        result.count = (int) intensities.size();
        result.name = "Raw Central+";
        results.ControlStats().push_back(result);
    }

    intensities = dataAccessor->GetIntensities(ReportDataAccessor::CentralCrossNegativeProbeSetType);
    if (intensities.size() > 0)
    {
        result.avg = Average(intensities);
        result.count = (int) intensities.size();
        result.name = "Raw Central-";
        results.ControlStats().push_back(result);
    }    
}

/*
 * Check if comparison data exists and if so then
 * count the number of each change call.
 */
void ExpressionProbeSetReporter::ComputeChangeProbeSetStats(bool includeAllProbesets)
{
	if (dataAccessor->HasComparisonData() == false)
		return;

	// Determine the statistics for the probe sets.
	int nSets = dataAccessor->GetNumProbeSets();
	for (int iSet=0; iSet<nSets; iSet++)
	{
		if (includeAllProbesets == false && IncludeProbeSet(iSet) == false)
			continue;

		// Count the number of change calls.
		ReportDataAccessor::ChangeCall change = dataAccessor->GetChange(iSet);
		ReportDataAccessor::DetectionCall edet = dataAccessor->GetDetection(iSet);
		ReportDataAccessor::DetectionCall bdet = dataAccessor->GetBaselineDetection(iSet);
		switch (change)
		{
		case ReportDataAccessor::ChangeIncrease:
			results.IncreaseStats().IncrementChangeCount();
			if (edet == ReportDataAccessor::DetectionPresent && bdet == ReportDataAccessor::DetectionPresent)
				results.IncreaseStats().IncrementDetectionPresentCount();
			else if (edet == ReportDataAccessor::DetectionPresent && bdet != ReportDataAccessor::DetectionPresent)
				results.IncreaseStats().IncrementDetectionChangeCount();
			else if (edet != ReportDataAccessor::DetectionPresent && bdet != ReportDataAccessor::DetectionPresent)
				results.IncreaseStats().IncrementDetectionAbsentCount();
			break;

		case ReportDataAccessor::ChangeDecrease:
			results.DecreaseStats().IncrementChangeCount();
			if (bdet == ReportDataAccessor::DetectionPresent && edet == ReportDataAccessor::DetectionPresent)
				results.DecreaseStats().IncrementDetectionPresentCount();
			else if (bdet == ReportDataAccessor::DetectionPresent && edet != ReportDataAccessor::DetectionPresent)
				results.DecreaseStats().IncrementDetectionChangeCount();
			else if (edet != ReportDataAccessor::DetectionPresent && bdet != ReportDataAccessor::DetectionPresent)
				results.DecreaseStats().IncrementDetectionAbsentCount();
			break;

		case ReportDataAccessor::ChangeModerateIncrease:
			results.IncreaseStats().IncrementModerateCount();
			if (edet == ReportDataAccessor::DetectionPresent && bdet == ReportDataAccessor::DetectionPresent)
				results.IncreaseStats().IncrementDetectionPresentCount();
			else if (edet == ReportDataAccessor::DetectionPresent && bdet != ReportDataAccessor::DetectionPresent)
				results.IncreaseStats().IncrementDetectionChangeCount();
			else if (edet != ReportDataAccessor::DetectionPresent && bdet != ReportDataAccessor::DetectionPresent)
				results.IncreaseStats().IncrementDetectionAbsentCount();
			break;

		case ReportDataAccessor::ChangeModerateDecrease:
			results.DecreaseStats().IncrementModerateCount();
			if (bdet == ReportDataAccessor::DetectionPresent && edet == ReportDataAccessor::DetectionPresent)
				results.DecreaseStats().IncrementDetectionPresentCount();
			else if (bdet == ReportDataAccessor::DetectionPresent && edet != ReportDataAccessor::DetectionPresent)
				results.DecreaseStats().IncrementDetectionChangeCount();
			else if (edet != ReportDataAccessor::DetectionPresent && bdet != ReportDataAccessor::DetectionPresent)
				results.DecreaseStats().IncrementDetectionAbsentCount();
			break;

		case ReportDataAccessor::ChangeNoChange:
			results.NoChangeStats().IncrementChangeCount();
			if (edet != ReportDataAccessor::DetectionPresent && bdet != ReportDataAccessor::DetectionPresent)
				results.NoChangeStats().IncrementDetectionAbsentCount();
			else if (edet == ReportDataAccessor::DetectionPresent && bdet == ReportDataAccessor::DetectionPresent)
				results.NoChangeStats().IncrementDetectionPresentCount();
			break;

		default:
			break;
		}

		// Count the number of SLR values within each bin range.
		float lower;
		float upper;
		float slr = dataAccessor->GetSignalLogRatio(iSet);
		float abs_slr = fabs(slr);
		switch (change)
		{
		case ReportDataAccessor::ChangeIncrease:
		case ReportDataAccessor::ChangeModerateIncrease:
			for (int ibin=0; ibin<NUMBER_FOLD_CHANGE_BINS && slr>=0; ibin++)
			{
				ChangeStats::GetBinRange(ibin, lower, upper);
				if (abs_slr >= lower && abs_slr < upper)
				{
					results.IncreaseStats().IncrementFoldChangeCount(ibin);
					break;
				}
			}
			break;

		case ReportDataAccessor::ChangeDecrease:
		case ReportDataAccessor::ChangeModerateDecrease:
			for (int ibin=0; ibin<NUMBER_FOLD_CHANGE_BINS && slr<=0; ibin++)
			{
				ChangeStats::GetBinRange(ibin, lower, upper);
				if (abs_slr >= lower && abs_slr < upper)
				{
					results.DecreaseStats().IncrementFoldChangeCount(ibin);
					break;
				}
			}
			break;
		default:
			break;
		}
	}
}

/*
 * Count the number of each call and compute the total signal for each call type.
 */
void ExpressionProbeSetReporter::ComputeDetectionProbeSetStats(bool includeAllProbesets)
{
	// Determine the statistics for the probe sets.
	int nSets = dataAccessor->GetNumProbeSets();
	for (int iSet=0; iSet<nSets; iSet++)
	{
		if (includeAllProbesets == false && IncludeProbeSet(iSet) == false)
			continue;

		// Add the call and signal to the results.
		results.ProbeSetResults().AddSet();
		float fIntensityValue = dataAccessor->GetSignal(iSet);
		ReportDataAccessor::DetectionCall det = dataAccessor->GetDetection(iSet);
		if (det == ReportDataAccessor::DetectionPresent)
		{
			results.ProbeSetResults().PresentCalls().IncrementCount();
			results.ProbeSetResults().PresentCalls().AddSignal(fIntensityValue);
		}
		else if (det == ReportDataAccessor::DetectionAbsent)
		{
			results.ProbeSetResults().AbsentCalls().IncrementCount();
			results.ProbeSetResults().AbsentCalls().AddSignal(fIntensityValue);
		}
		else if (det == ReportDataAccessor::DetectionMarginal)
		{
			results.ProbeSetResults().MarginalCalls().IncrementCount();
			results.ProbeSetResults().MarginalCalls().AddSignal(fIntensityValue);
		}
	}
}

/*
 * Compute the control statistics.
 */
void ExpressionProbeSetReporter::ComputeControlStats(ExpressionControlList &controls, ExpressionControlResultList &controlResults, bool bDifference)
{
	ExpressionControlList::iterator it;
	for (it=controls.begin(); it!=controls.end(); it++)
	{
		ExpressionControl &control = *it;
        float val = 0.0f;
		ExpressionControlResult controlResult;
		controlResult.SetName(control.Name());

		// 3'
		if (control.HasValue(ExpressionControl::THREE_PRIME_PROBE_SET))
		{
            controlResult.SetControlSignalResult(ExpressionControl::THREE_PRIME_PROBE_SET,
                dataAccessor->GetSignal(control.GetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET)));
            controlResult.SetControlDetectionResult(ExpressionControl::THREE_PRIME_PROBE_SET,
                dataAccessor->GetDetection(control.GetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET)));
		}

		// 5'
		if (control.HasValue(ExpressionControl::FIVE_PRIME_PROBE_SET))
		{
            controlResult.SetControlSignalResult(ExpressionControl::FIVE_PRIME_PROBE_SET,
                dataAccessor->GetSignal(control.GetProbeSetIndex(ExpressionControl::FIVE_PRIME_PROBE_SET)));
            controlResult.SetControlDetectionResult(ExpressionControl::FIVE_PRIME_PROBE_SET,
                dataAccessor->GetDetection(control.GetProbeSetIndex(ExpressionControl::FIVE_PRIME_PROBE_SET)));
        }

		// M
		if (control.HasValue(ExpressionControl::MIDDLE_PROBE_SET))
		{
            controlResult.SetControlSignalResult(ExpressionControl::MIDDLE_PROBE_SET,
                dataAccessor->GetSignal(control.GetProbeSetIndex(ExpressionControl::MIDDLE_PROBE_SET)));
            controlResult.SetControlDetectionResult(ExpressionControl::MIDDLE_PROBE_SET,
                dataAccessor->GetDetection(control.GetProbeSetIndex(ExpressionControl::MIDDLE_PROBE_SET)));
        }

		// 3'/5'
		if (control.HasValue(ExpressionControl::THREE_PRIME_PROBE_SET) &&
			control.HasValue(ExpressionControl::FIVE_PRIME_PROBE_SET))
		{
			if (bDifference == false)
			{
				val = dataAccessor->GetSignal(control.GetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET)) /
					max(1.0f, dataAccessor->GetSignal(control.GetProbeSetIndex(ExpressionControl::FIVE_PRIME_PROBE_SET)));
			}
			else
			{
				val = dataAccessor->GetSignal(control.GetProbeSetIndex(ExpressionControl::THREE_PRIME_PROBE_SET)) -
						dataAccessor->GetSignal(control.GetProbeSetIndex(ExpressionControl::FIVE_PRIME_PROBE_SET));
			}
			controlResult.SetThreeFiveRatio(val);
		}

        // Add the control to the list.
		controlResults.push_back(controlResult);
	}
}
