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

//
#include "calvin_files/converters/chp/src/CHPFileConverter.h"
//
#include "calvin_files/converters/chp/src/CHPFileConvertToCalvin.h"
#include "calvin_files/converters/chp/src/CHPFileConvertToXDA.h"
#include "calvin_files/converters/chp/src/CHPFileVersion.h"
#include "calvin_files/converters/chp/src/CmdLine.h"
#include "calvin_files/converters/utils/src/ConverterFileUtils.h"
//
#include <cstring>
#include <iostream>
#include <string>
//

using namespace std;
using namespace affymetrix_chp_converter;


/*
 * Initialize the class.
 */
CHPFileConverter::CHPFileConverter()
{
	Clear();
}

/*
 * Clears the members.
 */
void CHPFileConverter::Clear()
{
	errorCode = NoConversionError;
}

/*
 * Deallocate any memory.
 */
CHPFileConverter::~CHPFileConverter()
{
}

bool CHPFileConverter::ConvertFile(const char *fileName, const char *libPath, CHPFileVersionType fileVersion)
{
	return ConvertFile(fileName, libPath, NULL, NULL, fileVersion, true);
}
/*
 * Convert the file.
 */
bool CHPFileConverter::ConvertFile(const char *fileName, 
																	 const char *libPath, 
																	 const char *celFile, 
																	 const char *maskPath, 
																	 CHPFileVersionType fileVersion, 
																	 bool removeBackup)
{
	Clear();

	// Check if the file exists
	if (ConverterFileExists(fileName) == false)
	{
		errorCode = FileDoesNotExist;
		return false;
	}

	// Determine the format
	CHPFileVersionType inVersion = CHPFileVersion::DetermineCHPFileVersion(fileName);

	// If unknown the return a false.
	if (inVersion == Unknown_Version)
	{
		errorCode = InvalidChpFileFormat;
		return false;
	}

	// If same then don't do anything, the file is already of that format.
	if (inVersion == fileVersion)
	{
		return true;
	}

	// Move the input file name to a backup.
	string bakSuffix = ".bak";
	string bakFile = fileName + bakSuffix;
	if (ConverterMoveFile(fileName, bakFile.c_str(), true) == false)
	{
		errorCode = UnableToRenameInputFile;
		return false;
	}

	// Convert the file.
	bool status = false;
	if (fileVersion == GCOS_XDA_Version && inVersion == Calvin_Version1)
	{
		CHPFileConvertToXDA outFile;
		status = outFile.ConvertCalvinFile(bakFile.c_str(), fileName);
		errorCode = outFile.ErrorCode();
	}
	else if (fileVersion == Calvin_Version1 && inVersion == GCOS_XDA_Version)
	{
		CHPFileConvertToCalvin outFile;
		status = outFile.ConvertXDAFile(bakFile.c_str(), libPath, fileName);
		errorCode = outFile.ErrorCode();
	}
	else if(fileVersion == Calvin_Version1 && inVersion == GCOS_Mas5_Version && celFile != NULL)
	{
		CHPFileConvertToCalvin outFile;
		status = outFile.ConvertMas5File(bakFile.c_str(), celFile, libPath, maskPath, fileName);
		errorCode = outFile.ErrorCode();
	}
	else
	{
		errorCode = InvalidConversionInputs;
		status = false;
	}

	// Check the status and revert if needed.
	if (status == false)
	{
		ConverterMoveFile(bakFile.c_str(), fileName, true);
	}
	else if(removeBackup)
	{
		ConverterRemoveFile(bakFile.c_str());
	}
	return status;
}
