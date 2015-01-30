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


#include "calvin_files/converters/cel/src/CELFileConverterErrorCode.h"
//

/*
 * Return an error message for each error type.
 */
std::string affymetrix_cel_converter::CELFileConverterErrorMessage(CELFileConverterErrorCode code)
{
	std::string error="";
	switch (code)
	{
	case NoConversionError:
		error = "";
		break;

	case InvalidConversionInputs:
		error = "The inputs are invalid.";
		break;

	case FileDoesNotExist:
		error = "The input file does not exist.";
		break;

	case UnableToOpenCelFile:
		error = "Unable to open the input CEL file.";
		break;

	case InvalidCelFileFormat:
		error = "The input CEL file is not a CEL compatible with the converter.";
		break;

	case UnableToRenameInputFile:
		error = "Unable to rename the input file before the conversion process.";
		break;

	case UnableToWriteTheFile:
		error = "Unable to write the output file.";
		break;

    case UnableToMixConversionOptionsAndNoFormatChange:
		error = "Unable to mix conversion options with no file format change.";
		break;

    case UnableToCopyFile:
		error = "Unable to copy file.";
		break;

    case UnableToParseDatHeader:
		error = "Unable to parse DAT header.";
		break;

	default:
		error = "";
		break;
	}
	return error;

}

