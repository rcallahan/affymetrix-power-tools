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


#ifndef _CHPFileConverter_HEADER_
#define _CHPFileConverter_HEADER_

/*! file CHPFileConverter.h Defines a class for converting CHP files between GCOS and Calvin format. */

#include "calvin_files/converters/chp/src/CHPFileConverterErrorCode.h"
#include "calvin_files/converters/chp/src/CHPFileVersion.h"
//
#include <cstring>
#include <string>
//

namespace affymetrix_chp_converter
{

/*! This class will convert a CHP file between GCOS and Calvin formats. */
class CHPFileConverter
{
public:

	/*! Constructor */
	CHPFileConverter();

	/*! Destructor */
	~CHPFileConverter();

	/*! Gets the error code for the conversion process.
	 * @return The error code associated with the last failure.
	 */
	CHPFileConverterErrorCode ErrorCode() const { return errorCode; }

	/*! Convert the file to the specified format.
	 * @param fileName The full path name of the CHP file.
	 * @param libPath The full path to the library files directory.
	 * @param fileVersion The file version to convert to.
	 * @return True if successful.
	 */
	bool ConvertFile(const char *fileName, const char *libPath, CHPFileVersionType fileVersion);

	bool ConvertFile(const char *fileName, 
		const char *libPath, 
		const char *celFile, 
		const char *maskPath, 
		CHPFileVersionType fileVersion,
		bool removeBackup);

private:

	/*! Clears the members. */
	void Clear();

	/*! The error code. */
	CHPFileConverterErrorCode errorCode;
};

}

#endif
