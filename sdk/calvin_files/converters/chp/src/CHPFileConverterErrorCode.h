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


#ifndef _CHPFileConverterErrorCode_HEADER_
#define _CHPFileConverterErrorCode_HEADER_

/*! file CHPFileConverterErrorCode.h Defines error codes for converting CHP files. */

#include <cstring>
#include <string>
//

namespace affymetrix_chp_converter
{

/*! Error codes for the conversion process. */
typedef enum _CHPFileConverterErrorCode
{
	/*! No error. */
	NoConversionError,

	/*! Invalid inputs. */
	InvalidConversionInputs,

	/*! The file does not exist. */
	FileDoesNotExist,

	/*! Unable to open the input CHP file. */
	UnableToOpenChpFile,

	/*! The format is not valid. */
	InvalidChpFileFormat,

	/*! Unable to rename the input file before the conversion process. */
	UnableToRenameInputFile,

	/*! Unable to read the probe set names from the PSI file or CDF file. */
	UnableToLoadProbeSetNames,

	/*! Unable to write the output file. */
	UnableToWriteTheFile,

	/*! Unable to open the parent cel file. */
	UnableToOpenParentCelFile,

	/*! Unable to read the CDF file. */
	UnableToReadCdfFile,

	/*! Unable to convert a CHP based on the algorithm. */
	InvalidAlgorithmType,

	/*! Unable to convert a CHP based on the assay type. */
	InvalidAssayType

} CHPFileConverterErrorCode;

/*! Gets a string message associated with the error code.
 * @param code The error code.
 * @return A string describing the error.
 */
std::string CHPFileConverterErrorMessage(CHPFileConverterErrorCode code);

}

#endif
