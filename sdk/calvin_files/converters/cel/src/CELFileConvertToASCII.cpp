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


#include "calvin_files/converters/cel/src/CELFileConvertToASCII.h"
//
#include "calvin_files/converters/cel/src/CELConversionUtilities.h"
#include "calvin_files/data/src/CELData.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/CELAlgorithmParameterNames.h"
#include "calvin_files/parsers/src/CelFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
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
CELFileConvertToASCII::CELFileConvertToASCII()
{
	Clear();
}

/*
 * Clear the class
 */
void CELFileConvertToASCII::Clear()
{
	errorCode = NoConversionError;
}

/*
 * Deallocate any memory.
 */
CELFileConvertToASCII::~CELFileConvertToASCII()
{
	Clear();
}

/*
 * Convert the file.
 */
bool CELFileConvertToASCII::ConvertXDAFile(const char *fileName, const char *newFile, CELFileConversionOptions *options)
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
	outFile.SetFileFormat(CCELFileData::TEXT_CEL);
	CELConversionUtilities::Copy(inFile, outFile, options);


	// Write the output file.
	if (outFile.WriteTextCel() == false)
	{
		errorCode = UnableToWriteTheFile;
		return false;
	}

	return true;
}

/*
 * Convert the file.
 */
bool CELFileConvertToASCII::ConvertCalvinFile(const char *fileName, const char *newFile, CELFileConversionOptions *options)
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
	outFile.SetFileFormat(CCELFileData::TEXT_CEL);
	CELConversionUtilities::Copy(inFile, outFile, options);
	inFile.Clear();

	// Write the output file.
	if (outFile.WriteTextCel() == false)
	{
		errorCode = UnableToWriteTheFile;
		return false;
	}

	return true;
}
