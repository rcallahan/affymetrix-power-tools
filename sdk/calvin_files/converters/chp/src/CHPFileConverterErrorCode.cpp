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


#include "calvin_files/converters/chp/src/CHPFileConverterErrorCode.h"
//

/*
 * Return an error message for each error type.
 */
std::string affymetrix_chp_converter::CHPFileConverterErrorMessage(CHPFileConverterErrorCode code)
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

	case UnableToOpenChpFile:
		error = "Unable to open the input CHP file.";
		break;

	case InvalidChpFileFormat:
		error = "The input CHP file is not a CHP compatible with the converter.";
		break;

	case UnableToRenameInputFile:
		error = "Unable to rename the input file before the conversion process.";
		break;

	case UnableToLoadProbeSetNames:
		error = "Unable to read the probe set names from the PSI or CDF file.";
		break;

	case UnableToWriteTheFile:
		error = "Unable to write the output file.";
		break;

	case UnableToOpenParentCelFile:
		error = "Unable to open the parent CEL file which is required for non-XDA CHP file conversion.";
		break;

	case UnableToReadCdfFile:
		error = "Unable to read the CDF file.";
		break;

	case InvalidAlgorithmType:
		error = "The input CHP algorithm type is not supported by the converter.";
		break;

	case InvalidAssayType:
		error = "The input CHP assay type is not supported by the converter.";
		break;

	default:
		error = "";
		break;
	}
	return error;

}

