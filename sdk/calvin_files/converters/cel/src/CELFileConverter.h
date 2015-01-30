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


#ifndef _CELFileConverter_HEADER_
#define _CELFileConverter_HEADER_

/*! file CELFileConverter.h Defines a class for converting CEL files between GCOS and Calvin format. */

#include "calvin_files/converters/cel/src/CELFileConversionOptions.h"
#include "calvin_files/converters/cel/src/CELFileConverterErrorCode.h"
#include "calvin_files/converters/cel/src/CELFileVersion.h"
//
#include <cstring>
#include <string>
//

namespace affymetrix_cel_converter
{

/*! This class will convert a CEL file between GCOS and Calvin formats. */
class CELFileConverter
{
public:

	/*! Constructor */
	CELFileConverter();

	/*! Destructor */
	~CELFileConverter();

	/*! Gets the error code for the conversion process.
	 * @return The error code associated with the last failure.
	 */
	CELFileConverterErrorCode ErrorCode() const { return errorCode; }

	/*! Convert the file to the specified format.
	 * @param fileName The full path name of the CEL file.
	 * @param fileVersion The file version to convert to.
	 * @return True if successful.
	 */
	bool ConvertFile(const char *fileName, CELFileVersionType fileVersion, CELFileConversionOptions *options = NULL);

	/*! Convert the file to the specified format.
	 * @param fileName The full path name of the input CEL file.
	 * @param newFileName The full path name of the new CEL file to be written.
	 * @param fileVersion The file version to convert to.
	 * @return True if successful.
	 */
	bool ConvertFile(const char *fileName, const char *newFileName, CELFileVersionType fileVersion, CELFileConversionOptions *options = NULL);

private:

	/*! Clears the members. */
	void Clear();

	/*! The error code. */
	CELFileConverterErrorCode errorCode;

    /*! Basic Checks. */
    bool Checks(const char *fileName, CELFileVersionType fileVersion);

};

}

#endif
