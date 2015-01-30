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

#ifndef _Report_HEADER_
#define _Report_HEADER_

/*! \file ExpressionProbeSetReporter.h Defines classes to report on expression probe sets. */

#include "exp_report/src/ExpressionReportControls.h"
#include "exp_report/src/ExpressionReportData.h"
#include "exp_report/src/ReportDataAccessor.h"
//

namespace ExpressionReport
{

/*! Provides reporting of expression probe sets. */
class ExpressionProbeSetReporter
{
public:
	/*! Construtor */
	ExpressionProbeSetReporter();

	/*! Run the report on the input data.
	 * @param anti True if anti-sense probe sets should be reported, false if sense probe sets are to be reported.
	 * @param probePairThr The probe pair threshold.
	 * @param data Accessor to the expression data.
	 * @param controls The controls to collect information for.
     * @param probeSetIndicies A list of probe sets values to include in the report.
	 * @param includeAllProbesets Include all probe sets in the report
	 * @return True if successful.
	 */
	bool Run(bool anti, int probePairThr, ReportDataAccessor *data, ExpressionControls *controls, const std::list<int> &probeSetIndicies, bool bDifference=false, bool includeAllProbesets=false);

	/*! The report results. */
	ExpressionReportData &Results() { return results; }

	/*! Data accessor */
	ReportDataAccessor *DataAccessor() { return dataAccessor; }

	/*! Clears the members. */
	void Clear();

private:
	/*! Compute the detection probe set statistics. */
	void ComputeDetectionProbeSetStats(bool includeAllProbesets);

	/*! Compute the change probe set statistics. */
	void ComputeChangeProbeSetStats(bool includeAllProbesets);

    /*! Compute the corner/central probe sets statistics. */
    void ComputeQCControlProbeSets();

    /*! Add specific probe set signal values to the list.
     * @param probeSetIndicies A list of probe sets values to include in the report.
     */ 
    void AddProbeSetSignals(const std::list<int> &probeSetIndicies);

	/*! Determine if the probe set should be included in the stats.
	 * @return True if the probe set should be included.
	 */
	bool IncludeProbeSet(int index);

	/*! Compute the control statistics. */
	void ComputeControlStats(ExpressionControlList &controls, ExpressionControlResultList &controlResults, bool bDifference);

	/*! The report results. */
	ExpressionReportData results;

	/*! Data accessor */
	ReportDataAccessor *dataAccessor;
};

}

#endif
