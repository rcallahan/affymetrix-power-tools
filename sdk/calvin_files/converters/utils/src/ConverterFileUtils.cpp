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

#ifdef WIN32
#include "windows.h"
#endif

#include "calvin_files/converters/utils/src/ConverterFileUtils.h"
//
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/array/src/ArrayId.h"
//
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
//

/*
 * Move the file without checks.
 */
static bool MoveWithoutChecks(const char *oldFile, const char *newFile)
{
#ifdef WIN32
	return (MoveFile(oldFile, newFile) == TRUE);
#else
	return (rename(oldFile, newFile) == 0);
#endif
}

static bool CopyWithoutChecks(const char *in, const char *out)
{
    bool result = true;

    std::ifstream is (in,std::ios::binary);
    std::ofstream os (out,std::ios::binary);
    if(!is.good() || !os.good())
        result = false;

    os << is.rdbuf();
    if(!is.good() || !os.good())
        result = false;

    is.close();
    os.close();
    if(!is.good() || !os.good())
        result = false;

    return true;
}

/*
 * Check if the new file needs to be deleted first.
 * Then move the file.
 */
bool ConverterMoveFile(const char *oldFile, const char *newFile, bool overwrite)
{
	// No overwrite of the new file.
	if (overwrite == false)
	{
		// Dont move if the file exists.
		if (ConverterFileExists(newFile) == true)
			return false;

		// Move the file.
		return MoveWithoutChecks(oldFile, newFile);
	}

	// Overwrite if exists.
	else
	{
		// Remove the new file first.
		if (ConverterRemoveFile(newFile) == false)
			return false;

		// Move the file.
		return MoveWithoutChecks(oldFile, newFile);
	}
}

/*
 * Check if the new file needs to be deleted first.
 * Then copy the file.
 */
bool ConverterCopyFile(const char *oldFile, const char *newFile, bool overwrite)
{
	// No overwrite of the new file.
	if (overwrite == false)
	{
		// Dont move if the file exists.
		if (ConverterFileExists(newFile) == true)
			return false;

		// Move the file.
		return CopyWithoutChecks(oldFile, newFile);
	}

	// Overwrite if exists.
	else
	{
		// Remove the new file first.
		if (ConverterRemoveFile(newFile) == false)
			return false;

		// Move the file.
		return CopyWithoutChecks(oldFile, newFile);
	}
}

/*
 * Check if exists, if so then remove it.
 */
bool ConverterRemoveFile(const char *fileName)
{
	if (ConverterFileExists(fileName) == true)
	{
#ifdef WIN32
		if (DeleteFile(fileName) == false)
			return false;
#else
		if (remove(fileName) != 0)
			return false;
#endif
	}
	return true;
}

/*
 * Check if the file exists.
 */
bool ConverterFileExists(const char *fileName)
{
  return (affymetrix_calvin_utilities::FileUtils::Exists(fileName));
}

void ConverterCreateAndAddParentArrayGenericDataHeader(affymetrix_calvin_io::GenericDataHeader& hdr, std::string arrayID, std::wstring arrayBarcode)
{
	// Create a new parent GenericDataHeader and add to the current GenericDataHeader
	affymetrix_calvin_io::GenericDataHeader parentGDH;
	parentGDH.SetFileTypeId(ARRAY_TYPE_IDENTIFIER);

	ParameterNameValueType nvt;
	nvt.SetName(ARRAY_ID_PARAM_NAME);
	nvt.SetValueAscii(arrayID.c_str());
	parentGDH.AddNameValParam(nvt);

	nvt.SetName(ARRAY_BARCODE_PARAM_NAME);
	nvt.SetValueText(arrayBarcode.c_str());
	parentGDH.AddNameValParam(nvt);

	hdr.AddParent(parentGDH);
}
