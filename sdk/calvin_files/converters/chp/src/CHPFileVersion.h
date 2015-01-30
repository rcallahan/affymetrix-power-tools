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


#ifndef _CHPFileVersion_HEADER_
#define _CHPFileVersion_HEADER_

/*! \file CHPFileVersion.h Defines version types and functions to determine versions. */

namespace affymetrix_chp_converter
{
/*! The list of available CHP file types/versions. */
typedef enum _CHPFileVersion
{
	/*! An unknown or unsupported version. */
	Unknown_Version,

	/*! This is the XDA version format supported by GCOS. */
	GCOS_XDA_Version,

	/*! This is the version 1 Calvin format. */
	Calvin_Version1,

	GCOS_Mas5_Version

} CHPFileVersionType;


/*! Provides for information about CHP file versions. */
class CHPFileVersion
{
public:

	/*! Determines the version number of the file.
	* @param fileName The file name of the CHP file.
	* @return The file version.
	*/
	static CHPFileVersionType DetermineCHPFileVersion(const char *fileName);
};

}

#endif

