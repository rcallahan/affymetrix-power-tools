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

#ifndef _ExpressionProbeSetFileExtraction_HEADER_
#define _ExpressionProbeSetFileExtraction_HEADER_

/*! \file ExpressionProbeSetFileExtraction.h Defines a class to read the probe set file. */

#include <cstring>
#include <map>
#include <string>
//

typedef std::map<std::string, std::string> ProbeSetFileEntryMap;

/*! This class is used to read the a list of probe sets from a parameter (XML) file. */
class ExpressionProbeSetFileExtraction
{
public:
	/*! Extract the parameters from the AGCC format XML parameter file.
     * @param fileName The name of the parameter file. 
     * @param probeSets A map of probe set id to probe set name.
     * @return True if the parameters were extracted
     */
	static bool ExtractParameters(const char *fileName, ProbeSetFileEntryMap &probeSets);

};

#endif
