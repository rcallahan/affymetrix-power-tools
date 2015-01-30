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

#ifndef _MAS5ParameterExtraction_HEADER_
#define _MAS5ParameterExtraction_HEADER_

/*! \file MAS5ParameterExtraction.h Defines a class to read the MAS5 parameters from an XML file. */

#include "mas5-stat/src/ExpStatAlgSettings.h"
//

/*! Defines a class to read the MAS5 parameters from an XML file. */
class MAS5ParameterExtraction
{
public:
	/*! Extract the parameters from the AGCC format XML parameter file.
     * @param fileName The name of the parameter file. 
     * @param parameters The MAS5 algorithm parameters.
     * @return True if the parameters were extracted
     */
	static bool ExtractParameters(const char *fileName, CExpStatAlgSettings &parameters);

};

#endif
