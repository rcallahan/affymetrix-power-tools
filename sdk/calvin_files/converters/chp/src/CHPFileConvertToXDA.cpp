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


#include "calvin_files/converters/chp/src/CHPFileConvertToXDA.h"
//
#include "calvin_files/converters/chp/src/CHPConversionUtilities.h"
#include "calvin_files/data/src/CHPData.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include "file/CHPFileWriter.h"
//
#include <cstring>
#include <iostream>
#include <string>
//

using namespace std;
using namespace affxchp;
using namespace affxchpwriter;
using namespace affymetrix_chp_converter;
using namespace affymetrix_calvin_io;

/*
 * Initialize the class.
 */
CHPFileConvertToXDA::CHPFileConvertToXDA()
{
	Clear();
	extraParameters = NULL;
}

/*
 * Clear the class
 */
void CHPFileConvertToXDA::Clear()
{
	errorCode = NoConversionError;
}

/*
 * Deallocate any memory.
 */
CHPFileConvertToXDA::~CHPFileConvertToXDA()
{
	Clear();
}

/*
 * Convert the file.
 */
bool CHPFileConvertToXDA::ConvertCalvinFile(const char *fileName, const char *newFile)
{
	Clear();

	// Read the input file.
	CHPData inFile;
	CHPFileReader reader;
	reader.SetFilename(fileName);
	try
	{
		reader.Read(inFile);
	}
	catch (...)
	{
		errorCode = UnableToOpenChpFile;
		return false;
	}

	// Copy the input file object to the output file object.
	CCHPFileWriter outFile;
	outFile.SetFileName(newFile);
	CHPConversionUtilities::Copy(inFile, outFile);

	if (parentFileName.empty() == false)
	{
		// override the parent file name
		// Attempt to remove the path.
		string::size_type posBS = parentFileName.rfind('\\');
		string::size_type posS = parentFileName.rfind('/');
		if (posBS != -1)
		{
			outFile.SetParentCelFileName(parentFileName.substr(posBS+1).c_str());
		}
		else if (posS != -1)
		{
			outFile.SetParentCelFileName(parentFileName.substr(posS+1).c_str());
		}
		else
		{
			outFile.SetParentCelFileName(parentFileName.c_str());
		}
	}

	// Override chip type if it was set in extra parameters
	std::string chipType = GetChipTypeFromExtraParameters();
	if (chipType.length() > 0)
	{
		outFile.GetHeader().SetChipType(chipType.c_str());
	}

	// Write the output file.
	if (outFile.CreateNewFile() == false)
	{
		errorCode = UnableToWriteTheFile;
		return false;
	}
	if (outFile.Save() == false)
	{
		errorCode = UnableToWriteTheFile;
		return false;
	}

	return true;
}

std::string CHPFileConvertToXDA::GetChipTypeFromExtraParameters()
{
	std::string arrayType;
	if (extraParameters != NULL && extraParameters->size() > 0)
	{
		for (ParameterNameValueTypeVector::const_iterator ii = extraParameters->begin(); ii != extraParameters->end(); ++ii)
		{
			if (ii->GetName() == ARRAY_TYPE_PARAM_NAME)
			{
				arrayType = StringUtils::ConvertWCSToMBS(ii->GetValueText());
				break;
			}
		}
	}
	return arrayType;
}

