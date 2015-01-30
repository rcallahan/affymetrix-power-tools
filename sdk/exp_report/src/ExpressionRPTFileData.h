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

#ifndef _ExpressionRPTFileData_HEADER_
#define _ExpressionRPTFileData_HEADER_

/*! \file ExpressionRPTFileData.h This file provides reading capaibilities for expression RPT files.
 */

#include "exp_report/src/ExpressionReportData.h"
//
#include <cstring>
#include <fstream>
#include <string>
//

//////////////////////////////////////////////////////////////////////

namespace ExpressionReport
{

/*! This class provides reading capabilities for expression RPT files. */
class ExpressionRPTFileData
{
public:
	/*! Constructor */
	ExpressionRPTFileData() { fileName = ""; }

	/*! Destructor */
	~ExpressionRPTFileData() {}

protected:
	/*! The name of the Expression RPT file */
	std::string fileName;

public:
	/*! The name of the RPT file. */
	std::string &FileName() { return fileName; }

	/*! Reads the contents of the file.
	 * @param data The report data.
	 * @return True if successful.
	 */
	bool Read(ExpressionReportData &data);

	/*! Checks for the existance of a file.
	 * @return True if the file exists.
	 */
	bool Exists();
};

}

#endif
