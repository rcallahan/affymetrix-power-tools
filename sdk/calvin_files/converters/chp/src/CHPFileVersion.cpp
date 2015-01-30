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


#include "calvin_files/converters/chp/src/CHPFileVersion.h"
//
#include "calvin_files/data/src/CHPData.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
//
#include "file/CHPFileData.h"
//

using namespace affxchp;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_chp_converter;

/*
 * Check if the input file is a Calvin file.
 */
static bool IsCalvinCHPFile(const char *fileName)
{
	CHPData data;
	CHPFileReader reader;
	reader.SetFilename(fileName);
	try
	{
		reader.Read(data);
		data.Clear();
		return true;
	}
	catch (affymetrix_calvin_exceptions::CalvinException)
	{
	}
	catch (...)
	{
	}
	return false;
}


/*
 * Open the file as a Calvin file first then try as a GCOS file.
 */
CHPFileVersionType CHPFileVersion::DetermineCHPFileVersion(const char *fileName)
{
	if (IsCalvinCHPFile(fileName) == true)
	{
		return Calvin_Version1;
	}

	CCHPFileData chp;
	chp.SetFileName(fileName);
	if (chp.IsMas5File() == true)
	{
		return GCOS_Mas5_Version;
	}
	if (chp.IsXDACompatibleFile() == true)
	{
		return GCOS_XDA_Version;
	}

	return Unknown_Version;
}
