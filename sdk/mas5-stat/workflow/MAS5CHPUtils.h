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

#ifndef _MAS5CHPUtils_HEADER_
#define _MAS5CHPUtils_HEADER_

/*! \file MAS5CHPUtils.h Defines a class to assist with writing a MAS5 CHP file. */

#include "mas5-stat/src/ExpressionAlgorithmImplementation.h"
//
#include "calvin_files/data/src/CHPData.h"
#include "exp_report/src/ExpressionProbeSetReporter.h"
#include "file/CHPFileWriter.h"
//
#include <cstring>
#include <string>
//

/*! Provides utilities for writing MAS5 CHP files. */
class MAS5CHPUtils
{
private:

	/*! The command console CHP data object. */
	affymetrix_calvin_io::CHPData *cc_data;

	/*! The GCOS CHP data object. */
	affxchpwriter::CCHPFileWriter *gcos_data;

	/*! The algorithm results. */
	CExpressionAlgorithmImplementation *mas5_data;

	/*! The report results. */
	ExpressionReport::ExpressionProbeSetReporter *report_data;

    /*! The program name. */
    std::wstring programName;
    
    /*! The company name. */
    std::wstring programCompany;
    
    /*! The ID of the program. */
    std::wstring programId;

	/*! Add a name/value pair to the algorithm parameter section.
	 * @param name The parameter name.
	 * @param value The parameter value.
	 */
	void AddAlgParam(const std::wstring &name, const std::wstring &value);

	/*! Add a name/value pair to the algorithm parameter section.
	 * @param name The parameter name.
	 * @param value The parameter value.
	 */
	void AddAlgParam(const std::wstring &name, int value);

	/*! Add a name/value pair to the algorithm parameter section.
	 * @param name The parameter name.
	 * @param value The parameter value.
	 */
	void AddAlgParam(const std::wstring &name, float value);

	/*! Add a name/value pair to the chip summary section.
	 * @param name The parameter name.
	 * @param value The parameter value.
	 */
	void AddChipSum(const std::wstring &name, const std::wstring &value);

	/*! Add a name/value pair to the chip summary section.
	 * @param name The parameter name.
	 * @param value The parameter value.
	 */
	void AddChipSum(const std::wstring &name, float value);

	/*! Add a name/value pair to the chip summary section.
	 * @param name The parameter name.
	 * @param value The parameter value.
	 */
	void AddChipSum(const std::wstring &name, int value);

	/*! Store program information. */
	void StoreProgramInfo();

	/*! Store the algorithm parameters to the chp object. */
	void StoreAlgParams();

	/*! Store the summary parameters to the chp object.
     * @param separateEntries True if the average, stdv, min and max values are to be stored in separate parameters.
     */
	void StoreSummaryParams(bool separateEntries);

	/*! Store the report results to the chp object. */
	void StoreReport();

	/*! Stores the report parameters. */
	void StoreReportParameters();

	/*! Store the report detection stats. */
	void StoreReportDetectionStats();

	/*! Store the report control stats. */
	void StoreReportControlStats();

	/*! Store the report change stats. */
	void StoreReportChangeStats();

	/*! Store a control statistic.
	 * @param controlType A string describing the type of control - Spike or housekeeping.
	 * @param result The control results.
	 */
	void StoreReportControlStat(const wchar_t *controlType,
                              const ExpressionReport::ExpressionControlResult &result);

	/*! Store the parent headers to the command console file. */
	void StoreParentHeaders();

	/*! Save the data to a Command Console CHP file.
	 * @param schpFile The name of the output CHP file.
	 * @return True if successful.
	 */
	bool SaveCommandConsoleCHPFile(const std::string &schpFile);

	/*! Save the data to a GCOS CHP file.
	 * @param schpFile The name of the output CHP file.
	 * @return True if successful.
	 */
	bool SaveGCOSCHPFile(const std::string &schpFile);

public:

	/*! Constructor */
    MAS5CHPUtils() { Clear(); }

    /*! Set the program informaion.
     * @param name The name of the program
     * @param co The company or institution creating the program.
     * @param id The identifier of the program.
     */
    void SetProgramInformation(const std::wstring &name, const std::wstring &co, const std::wstring &id); 

	/*! Clears the members. */
	void Clear();

	/*! Save the data to a CHP file.
	 * @param schpFile The name of the output CHP file.
	 * @param exp The MAS5 results object.
	 * @param report The report results.
	 * @param reportPass True if the report has results.
	 * @param saveToLegacyFile True is saving to older CHP file format.
	 * @return True if successful.
	 */
	bool SaveCHPFile(const std::string &chpFile, CExpressionAlgorithmImplementation &exp, ExpressionReport::ExpressionProbeSetReporter &report, bool reportPass, bool saveToLegacyFile);
};

#endif
