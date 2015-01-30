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

#ifndef _ReportDataAccessor_HEADER_
#define _ReportDataAccessor_HEADER_

/*! \file ReportDataAccessor.h Defines an interface to access data for the expression report. */

#include <vector>

namespace ExpressionReport
{

/*! An interface to access data for the expression report. */
class ReportDataAccessor
{
public:
	/*! The values of the detection call. */
	typedef enum _DetectionCall
	{
		/*! Probe set is detected. */
		DetectionPresent = 0,

		/*! Probe set is marginally detected. */
		DetectionMarginal,

		/*! Probe set is not detected. */
		DetectionAbsent,

		/*! No call */
		DetectionNoCall

	} DetectionCall;

	/*! The values of the change call */
	typedef enum _ChangeCall
	{
		/*! Increase call for expression comparison analysis */
		ChangeIncrease = 1,

		/*! Decrease call for expression comparison analysis */
		ChangeDecrease,

		/*! Moderate increase call for expression comparison analysis */
		ChangeModerateIncrease,
		
		/*! Moderate decrease call for expression comparison analysis */
		ChangeModerateDecrease,

		/*! No change call for expression comparison analysis */
		ChangeNoChange,

		/*! No call call for expression comparison analysis */
		ChangeNoCall

	} ChangeCall;

        
    /*! The type of QC probe set */
    typedef enum _QCProbeSetType
    {
	    /*! Probes used for the checker board patterns for antisense arrays. */
	    CheckerboardNegativeProbeSetType,

	    /*! Probes used for the checker board patterns for sense arrays. */
	    CheckerboardPositiveProbeSetType,

	    /*! Central cross probes for antisense arrays. */
	    CentralCrossNegativeProbeSetType,

	    /*! Central cross probes for sense arrays. */
	    CentralCrossPositiveProbeSetType,

    } QCProbeSetType;

	/* Constructor */
	ReportDataAccessor() {};

	/*! Destructor */
	virtual ~ReportDataAccessor() {};

	/*! Gets the number of expression probe sets. */
	virtual int GetNumProbeSets() = 0;

	/*! Gets the number of used probe pairs in a probe set.
	 * @param index The probe set index.
	 * @return The number of probe pairs.
	 */
	virtual int GetNumPairs(int index) = 0;

	/*! Checks if a probe set is targeting an anti-sense target.
	 * @param index The probe set index.
	 * @return True if the probe set is designed to interrogate an anti-sense target.
	 */
	virtual bool IsAntiSense(int index) = 0;

	/*! Gets the probe set signal value
	 * @param index The probe set index.
	 * @return The signal value.
	 */
	virtual float GetSignal(int index) = 0;

	/*! Gets the probe set detection value
	 * @param index The probe set index.
	 * @return The detection value.
	 */
	virtual DetectionCall GetDetection(int index) = 0;

	/*! Checks if the data has comparison results.
	 * @return True if comparison data exists.
	 */
	virtual bool HasComparisonData() = 0;

	/*! Gets the probe set detection value for the baseline file.
	 * @param index The probe set index.
	 * @return The detection value.
	 */
	virtual DetectionCall GetBaselineDetection(int index) = 0;

	/*! Gets the change call.
	 * @param index The probe set index.
	 * @return The change call.
	 */
	virtual ChangeCall GetChange(int index) = 0;

	/*! Gets the signal log ratio.
	 * @param index The probe set index.
	 * @return The SLR value.
	 */
	virtual float GetSignalLogRatio(int index) = 0;

    /*! Gets the list of intensities for a QC probe set.
     * @param qctype The type of QC probe set.
     * @return The list of intensities.
     */
    virtual std::vector<float> GetIntensities(QCProbeSetType qctype) = 0;

    /*! Gets the name of a probe set at the given index position.
     * @param index The probe set index.
     * @return The probe set name.
     */
    virtual std::string GetProbeSetName(int index) = 0;

	/*! Check whether the data is reported using log (base 2) scale.
	 */
	virtual bool IsLogScale() = 0;
};

}

#endif
