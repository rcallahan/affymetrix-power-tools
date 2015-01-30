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


#include "calvin_files/converters/chp/comparer/CHPComparer.h"
//
#include "calvin_files/converters/chp/comparer/CHPCompareCalvintoCalvin.h"
#include "calvin_files/converters/chp/comparer/CHPCompareGCOStoCalvin.h"
#include "calvin_files/converters/chp/comparer/CHPCompareGCOStoGCOS.h"
#include "calvin_files/converters/chp/src/CHPFileVersion.h"
//
#include <cmath>
//

using namespace affymetrix_comparer;
using namespace affymetrix_chp_converter;
using namespace std;


/*
 * Clear the members.
 */
void CHPComparer::Clear()
{
	differences = "";
}

/*
 * Compare the versions.
 */
bool CHPComparer::CompareFiles(const char *file1, const char *file2)
{
	bool status=false;

	CHPFileVersionType ver1 = CHPFileVersion::DetermineCHPFileVersion(file1);
	CHPFileVersionType ver2 = CHPFileVersion::DetermineCHPFileVersion(file2);

	// GCOS to Calvin comparison.
	if (((ver1 == GCOS_XDA_Version || ver1 == GCOS_Mas5_Version) && ver2 == Calvin_Version1) ||
		(ver1 == Calvin_Version1 && (ver2 == GCOS_XDA_Version || ver2 == GCOS_Mas5_Version)))
	{
		const char *gcosFile = (ver1 == GCOS_XDA_Version ? file1 : file2 );
		const char *calvinFile = (ver1 == Calvin_Version1 ? file1 : file2 );
		CHPCompareGCOStoCalvin compare;
		status = compare.CompareFiles(gcosFile, calvinFile);
		differences = compare.Differences();
	}

	// GCOS to GCOS comparison.
	else if ((ver1 == GCOS_XDA_Version || ver1 == GCOS_Mas5_Version) && (ver2 == GCOS_XDA_Version || ver2 == GCOS_Mas5_Version))
	{
		CHPCompareGCOStoGCOS compare;
		status = compare.CompareFiles(file1, file2);
		differences = compare.Differences();
	}

	// Calvin to Calvin comparison.
	else if (ver1 == Calvin_Version1 && ver2 == Calvin_Version1)
	{
		CHPCompareCalvintoCalvin compare;
		status = compare.CompareFiles(file1, file2);
		differences = compare.Differences();
	}
	else
	{
		differences = "file comparisons not supported.";
		status = false;
	}

	return status;
}
