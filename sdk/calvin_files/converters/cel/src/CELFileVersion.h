////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////


#ifndef _CELFileVersion_HEADER_
#define _CELFileVersion_HEADER_

/*! \file CELFileVersion.h Defines version types and functions to determine versions. */

namespace affymetrix_cel_converter
{
/*! The list of available CEL file types/versions. */
typedef enum _CELFileVersion
{
	/*! An unknown or unsupported version. */
	Unknown_Version,

	/*! This is the version 3 or ASCII format supported by GCOS. */
	GCOS_Version3,		

	/*! This is the version 4 or binary or XDA format supported by GCOS. */
	GCOS_Version4,

	/*! This is the version 1 Calvin format. */
	Calvin_Version1

} CELFileVersionType;


/*! Provides for information about CEL file versions. */
class CELFileVersion
{
public:

	/*! Determines the version number of the file.
	* @param fileName The file name of the CEL file.
	* @return The file version.
	*/
	static CELFileVersionType DetermineCELFileVersion(const char *fileName);
};

}

#endif

