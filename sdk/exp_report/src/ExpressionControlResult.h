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

#ifndef _ExpressionControlResult_HEADER_
#define _ExpressionControlResult_HEADER_

/*! \file ExpressionControlResult.h Defines classes to store results for 3/m/5 prime controls. */

#include "exp_report/src/ExpressionReportControls.h"
#include "exp_report/src/ReportDataAccessor.h"
//
#include <cstring>
#include <list>
#include <string>
//

namespace ExpressionReport
{

/*! The results for a control. */
class ExpressionControlResult
{
	/*! The name of the control. */	
	private: std::string name;

	/*! Get the name of the control. */
	public: std::string GetName() const { return name; }

	/*! Set the name of the control. */
	public: void SetName(const std::string &n) { name=n; }

	/*! The control signal values. */
	private: float signalValues[NUM_SETS_PER_CONTROL];

    /*! The control detection values. */
    private: ReportDataAccessor::DetectionCall detectionValues[NUM_SETS_PER_CONTROL];

	/*! Three/five prime ratio. */
	private: float threeFiveRatio;

    /*! Returns the control signal result.
     * @param resultType The type of result
     * @return The result value
	 */
	public: float GetControlSignalResult(ExpressionControl::ControlValueType resultType) const
    {
        return signalValues[(int)resultType];
    }

    /*! Sets the control signal result.
     * @param resultType The type of result
     * @param val The value of the result
	 */
    public: void SetControlSignalResult(ExpressionControl::ControlValueType resultType, float val)
    {
        signalValues[(int)resultType] = val;
    }

    /*! Returns the control detection result.
     * @param resultType The type of result
     * @return The result value
	 */
	public: ReportDataAccessor::DetectionCall GetControlDetectionResult(ExpressionControl::ControlValueType resultType) const
    {
        return detectionValues[(int)resultType];
    }

    /*! Sets the control detection result
     * @param resultType The type of result
     * @param val The value of the result
     */
	public: void SetControlDetectionResult(ExpressionControl::ControlValueType resultType, ReportDataAccessor::DetectionCall val)
    {
        detectionValues[(int)resultType] = val;
    }

	/*! Three/five prime ratio. */
	public: float GetThreeFiveRatio() const { return threeFiveRatio; }

	/*! Three/five prime ratio. */
	public: void SetThreeFiveRatio(float val) { threeFiveRatio = val; }

    /*! Flag indicating if control result exists.
    * @param resultType The type of result
    * @return Flag indicating if result exists
	*/
    public: bool HasControlResult(ExpressionControl::ControlValueType resultType) const
    {
        return (signalValues[(int)resultType] >= 0);
    }

	/*! Initializes the class members. */
	public: ExpressionControlResult()
	{
		Clear();
	}

	/*! Initialize the class members. */
	public: void Clear()
	{
		threeFiveRatio = -1.0f;
		for (int i = 0; i < NUM_SETS_PER_CONTROL; i++)
		{
            signalValues[i] = -1.0f;
			detectionValues[i] = ReportDataAccessor::DetectionNoCall;
		}
	}
};

/*! A list of control results. */
typedef std::list<ExpressionControlResult> ExpressionControlResultList;

}

#endif
