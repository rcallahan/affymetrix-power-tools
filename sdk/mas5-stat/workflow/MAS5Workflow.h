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

#ifndef _MAS5Workflow_HEADER_
#define _MAS5Workflow_HEADER_

/*! \file MAS5Workflow.h Defines a class to provide the workflow of calling the MAS5 algorithm and saving the results to a CHP file. */

#include "mas5-stat/src/ExpressionAlgorithmImplementation.h"
//
#include "exp_report/src/ExpressionReportControls.h"
//
#include <cstring>
#include <string>
//

/*! Provide the workflow of calling the MAS5 algorithm and saving the results to a CHP file. */
class MAS5Workflow
{
	/*! An error message. */
	std::string errorMsg;

	/*! The MAS5 algorithm implmentation object. */
	CExpressionAlgorithmImplementation exp;

	/*! The probe pair threshold for the report. */
	int probePairThreshold;

	/*! The expression controls for the report. */
	ExpressionReport::ExpressionControls controls;

	/*! Indicates if the data should be saved to the old format CHP file. */
	bool saveToLegacyFile;

    /*! The program name. */
    std::wstring programName;
    
    /*! The company name. */
    std::wstring programCompany;
    
    /*! The ID of the program. */
    std::wstring programId;

public:
	/*! Constructor */
	MAS5Workflow();

	/*! Destructor */
	~MAS5Workflow();

	/*! Clears them members. */
	void Clear();

	/*! Get the last error. */
	std::string GetError() { return errorMsg; }

	/*! Gets the algorithm parameters. */
	CExpStatAlgSettings &AlgParameters() { return exp.GetParameters(); }

	/*! The probe pair threshold for the report. */
	int &ProbePairThreshold() { return probePairThreshold; }

	/*! Indicates if the data should be saved to the old format CHP file. */
	bool &SaveToLegacyFile() { return saveToLegacyFile; }

	/*! The expression controls for the report. */
	ExpressionReport::ExpressionControls &ReportControls() { return controls; }

	/*! Executes the algorithm, report and saves the results to a CHP file.
	 * @param celFile The input CEL file.
	 * @param baselineCelFile The baseline file, empty if no comparison analysis is to be performed.
	 * @param cdfFile The CDF file name (full path).
	 * @param chpFile The name of the output CHP file.
	 * @return True if success.
	 */
	bool Execute(const std::string &celFile, const std::string &baselineCelFile, const std::string &cdfFile, const std::string &chpFile);

    /*! The program name. */
    std::wstring &ProgramName() { return programName; }
    
    /*! The company name. */
    std::wstring &ProgramCompany() { return programCompany; }
    
    /*! The ID of the program. */
    std::wstring &ProgramId() { return programId; }

};

#endif
