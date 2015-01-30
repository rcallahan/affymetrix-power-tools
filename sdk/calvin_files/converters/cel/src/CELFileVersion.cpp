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


#include "calvin_files/converters/cel/src/CELFileVersion.h"
//
#include "calvin_files/data/src/CELData.h"
#include "calvin_files/parsers/src/CelFileReader.h"
//
#include "file/CELFileData.h"
//

using namespace affxcel;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_cel_converter;

/*
 * Check if the input file is a Calvin file.
 */
static bool IsCalvinCELFile(const char *fileName)
{
	CelFileData data;
	CelFileReader reader;
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
CELFileVersionType CELFileVersion::DetermineCELFileVersion(const char *fileName)
{
	if (IsCalvinCELFile(fileName) == true)
	{
		return Calvin_Version1;
	}

	CCELFileData cel;
	cel.SetFileName(fileName);
	if (cel.IsXDACompatibleFile() == true)
	{
		return GCOS_Version4;
	}

	if (cel.IsVersion3CompatibleFile() == true)
	{
		return GCOS_Version3;
	}

	return Unknown_Version;
}
