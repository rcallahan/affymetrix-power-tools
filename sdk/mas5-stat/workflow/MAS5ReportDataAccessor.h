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

#ifndef _MAS5ReportDataAccessor_HEADER_
#define _MAS5ReportDataAccessor_HEADER_

/*! \file MAS5ReportDataAccessor.h Defines a class to access the MAS5 results. */

#include "mas5-stat/src/ExpressionAlgorithmImplementation.h"
//
#include "exp_report/src/ReportDataAccessor.h"
//

/*! A class to access the MAS5 results. */
class MAS5ReportDataAccessor : public ExpressionReport::ReportDataAccessor
{
private:
	/*! The expression results. */
	CExpressionAlgorithmImplementation *exp;

public:
	/*! Sets the expression results. */
	MAS5ReportDataAccessor(CExpressionAlgorithmImplementation *data) { exp = data; }
	
	/*! Gets the number of expression probe sets. */
	int GetNumProbeSets() { return exp->GetNumResults(); }

	/*! Gets the number of used probe pairs in a probe set.
	 * @param index The probe set index.
	 * @return The number of probe pairs.
	 */
	int GetNumPairs(int index)
	{
		AbsStatExpressionProbeSetResultType *result = exp->GetAbsStatResult(index);
		return result->NumUsedPairs;
	}

	/*! Checks if a probe set is targeting an anti-sense target.
	 * @param index The probe set index.
	 * @return True if the probe set is designed to interrogate an anti-sense target.
	 */
	bool IsAntiSense(int index)
	{
		std::string name = exp->GetProbeSetName(index);
		std::string ext = name.substr(name.length()-2, 2);
		if (ext == "AT" || ext == "at")
			return true;
		return false;
	}

	/*! Gets the probe set signal value
	 * @param index The probe set index.
	 * @return The signal value.
	 */
	float GetSignal(int index)
	{
		AbsStatExpressionProbeSetResultType *result = exp->GetAbsStatResult(index);
		return result->Signal;
	}

	/*! Gets the probe set detection value
	 * @param index The probe set index.
	 * @return The detection value.
	 */
	DetectionCall GetDetection(int index)
	{
		AbsStatExpressionProbeSetResultType *result = exp->GetAbsStatResult(index);
		return (DetectionCall) result->Detection;
	}

	/*! Checks if the data has comparison results.
	 * @return True if comparison data exists.
	 */
	bool HasComparisonData() { return exp->DoesCompDataExists(); }

	/*! Gets the probe set detection value for the baseline file.
	 * @param index The probe set index.
	 * @return The detection value.
	 */
	DetectionCall GetBaselineDetection(int index)
	{
		AbsStatExpressionProbeSetResultType *result = exp->GetBaselineAbsStatResult(index);
		return (DetectionCall) result->Detection;
	}

	/*! Gets the change call.
	 * @param index The probe set index.
	 * @return The change call.
	 */
	ChangeCall GetChange(int index)
	{
		CompStatExpressionProbeSetResultType *result = exp->GetCompStatResult(index);
		return (ChangeCall) result->Change;
	}

	/*! Gets the signal log ratio.
	 * @param index The probe set index.
	 * @return The SLR value.
	 */
	float GetSignalLogRatio(int index)
	{
		CompStatExpressionProbeSetResultType *result = exp->GetCompStatResult(index);
		return result->SignalLogRatio;
	}

	/*! Determines if the direction is anti-sense.
	 * @return True if anti-sense.
	 */
	bool IsAntiSense()
	{
		int antiSenseCount = 0;
		int n = exp->GetNumResults();
		for (int i=0; i<n; i++)
		{
			if (IsAntiSense(i) == true)
				++antiSenseCount;
		}
		return (antiSenseCount > n/2);
	}

    /*! Don't return anything here. The control (central/corner) probe setse
     * get reported in the MAS5 algorithm.
     * @param qctype The type of QC probe set.
     * @return The list of intensities.
     */
    std::vector<float> GetIntensities(QCProbeSetType qctype)
    {
        std::vector<float> v;
        return v;
    }

    /*! Gets the name of a probe set at the given index position.
     * @param index The probe set index.
     * @return The probe set name.
     */
    std::string GetProbeSetName(int index) { return ""; }
	
	bool IsLogScale()
	{
		return false;
	}
};

#endif
