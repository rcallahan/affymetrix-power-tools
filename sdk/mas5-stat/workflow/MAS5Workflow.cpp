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

#include "mas5-stat/workflow/MAS5Workflow.h"
//
#include "mas5-stat/workflow/MAS5CHPUtils.h"
#include "mas5-stat/workflow/MAS5ReportDataAccessor.h"
//
#include "exp_report/src/ExpressionProbeSetReporter.h"
//

using namespace ExpressionReport;

#define DEFAULT_PROBE_PAIR_THRESHOLD 1

/*
 * Clear the members.
 */
MAS5Workflow::MAS5Workflow()
{
	Clear();
}

/*
 * Clear the members.
 */
MAS5Workflow::~MAS5Workflow()
{
	Clear();
}

/*
 * Clear the members.
 */
void MAS5Workflow::Clear()
{
	errorMsg = "";
	controls.Clear();
	exp.Clear();
    programName.clear();
    programCompany.clear();
    programId.clear();
    probePairThreshold = DEFAULT_PROBE_PAIR_THRESHOLD;
	saveToLegacyFile = false;
}

/*
 * Run the algorithm and report and save the results to a CHP file
 */
bool MAS5Workflow::Execute(const std::string &celFile, const std::string &baselineCelFile, const std::string &cdfFile, const std::string &chpFile)
{
	errorMsg = "";
	if (exp.RunStat(celFile.c_str(), baselineCelFile.c_str(), cdfFile.c_str()) == false)
	{
		errorMsg = exp.GetError();
		return false;
	}

	// Run the report
	MAS5ReportDataAccessor dataAccessor(&exp);
	ExpressionProbeSetReporter report;
    const std::list<int> probeSetIndicies;
	bool reportPass = report.Run(dataAccessor.IsAntiSense(), probePairThreshold, &dataAccessor, &controls, probeSetIndicies, false);
	if (reportPass == false)
	{
		errorMsg = "Failed to run the report";
		return false;
	}

	// Save the file.
	MAS5CHPUtils utils;
    utils.SetProgramInformation(programName, programCompany, programId);
	if (utils.SaveCHPFile(chpFile, exp, report, reportPass, saveToLegacyFile) == false)
	{
		errorMsg = "Unable to write the CHP file.";
		return false;
	}
	return true;
}
