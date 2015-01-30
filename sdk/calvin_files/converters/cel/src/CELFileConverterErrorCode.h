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


#ifndef _CELFileConverterErrorCode_HEADER_
#define _CELFileConverterErrorCode_HEADER_

/*! file CELFileConverterErrorCode.h Defines error codes for converting CEL files. */

#include <cstring>
#include <string>
//

namespace affymetrix_cel_converter
{

/*! Error codes for the conversion process. */
typedef enum _CELFileConverterErrorCode
{
	/*! No error. */
	NoConversionError,

	/*! Invalid inputs. */
	InvalidConversionInputs,

	/*! The file does not exist. */
	FileDoesNotExist,

	/*! Unable to open the input CEL file. */
	UnableToOpenCelFile,

	/*! The format is not valid. */
	InvalidCelFileFormat,

	/*! Unable to rename the input file before the conversion process. */
	UnableToRenameInputFile,

	/*! Unable to write the output file. */
	UnableToWriteTheFile,

	/*! Unable to mix options and no format change */
	UnableToMixConversionOptionsAndNoFormatChange,

    /*! Unable to copy the file */
    UnableToCopyFile,

    /*! Unable to copy the file */
    UnableToParseDatHeader

} CELFileConverterErrorCode;

/*! Gets a string message associated with the error code.
 * @param code The error code.
 * @return A string describing the error.
 */
std::string CELFileConverterErrorMessage(CELFileConverterErrorCode code);

}

#endif
