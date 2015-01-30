////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

#ifndef _ExpressionControlsParameterExtraction_HEADER_
#define _ExpressionControlsParameterExtraction_HEADER_

/*! \file ExpressionControlsParameterExtraction.h Defines a class to read the controls for the from an XML file. */

#include "exp_report/src/ExpressionReportControls.h"
//

/*! Defines a class to read the ExpressionControls parameters from an XML file. */
class ExpressionControlsParameterExtraction
{
public:
	/*! Extract the parameters from the AGCC format XML parameter file.
     * @param fileName The name of the parameter file. 
     * @param probePairThreshold The probe pair threshold for the report.
     * @param controls The controls for the report.
     * @return True if the parameters were extracted
     */
	static bool ExtractParameters(const char *fileName, int &probePairThreshold, ExpressionReport::ExpressionControls &controls);

};

#endif
