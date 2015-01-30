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


#include "calvin_files/converters/cel/src/CELFileConvertToXDA.h"
//
#include "calvin_files/converters/cel/src/CELConversionUtilities.h"
#include "calvin_files/parameter/src/CELAlgorithmParameterNames.h"
#include "calvin_files/parsers/src/CelFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include "file/CELFileWriter.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affxcel;
using namespace affymetrix_calvin_io;
using namespace affymetrix_cel_converter;
using namespace affymetrix_calvin_utilities;

/*
 * Initialize the class.
 */
CELFileConvertToXDA::CELFileConvertToXDA()
{
	Clear();
	extraParameters = NULL;
}

/*
 * Clear the class
 */
void CELFileConvertToXDA::Clear()
{
	errorCode = NoConversionError;
}

/*
 * Deallocate any memory.
 */
CELFileConvertToXDA::~CELFileConvertToXDA()
{
	Clear();
}

/*
 * Convert the file.
 */
bool CELFileConvertToXDA::ConvertASCIIFile(const char *fileName, const char *newFile, CELFileConversionOptions *options)
{
	Clear();

	// Read the input file.
	CCELFileData inFile;
	inFile.SetFileName(fileName);
	if (inFile.Read() == false)
	{
		errorCode = UnableToOpenCelFile;
		return false;
	}


	// Copy the input file object to the output file object.
	CCELFileWriter outFile;
	outFile.SetFileName(newFile);
	outFile.SetFileFormat(CCELFileData::XDA_BCEL);
	if(!CELConversionUtilities::Copy(inFile, outFile, options))
    {
        return false;
    }


	// Write the output file.
	if (outFile.WriteXDABCel() == false)
	{
		errorCode = UnableToWriteTheFile;
		return false;
	}

	return true;
}

/*
 * Convert the file.
 */
bool CELFileConvertToXDA::ConvertCalvinFile(const char *fileName, const char *newFile, CELFileConversionOptions *options)
{
	Clear();

	// Read the input file.
	CelFileData inFile;
	CelFileReader reader;
	reader.SetFilename(fileName);
	try
	{
		reader.Read(inFile);
	}
	catch (...)
	{
		errorCode = UnableToOpenCelFile;
		return false;
	}

	// Copy the input file object to the output file object.
	CCELFileWriter outFile;
	outFile.SetFileName(newFile);
	outFile.SetFileFormat(CCELFileData::XDA_BCEL);
	if(!CELConversionUtilities::Copy(inFile, outFile, options))
    {
        return false;
    }

	// Override chip type if it was set in extra parameters
	std::string chipType = GetChipTypeFromExtraParameters();
	if (chipType.length() > 0)
	{
		outFile.SetChipType(chipType.c_str());
	}

	// Write the output file.
	if (outFile.WriteXDABCel() == false)
	{
		errorCode = UnableToWriteTheFile;
		return false;
	}

	return true;
}

std::string CELFileConvertToXDA::GetChipTypeFromExtraParameters()
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
