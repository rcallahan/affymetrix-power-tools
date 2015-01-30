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


#include "calvin_files/converters/chp/comparer/CHPCompareUtils.h"
//
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cmath>
//

using namespace affymetrix_calvin_utilities;
using namespace affymetrix_comparer;
using namespace std;

const string ret = "\r\n";

/*
 * Check if two floats are about the same.
 */
void CHPCompareUtils::CompareFloats(float a, float b, string &diffs, const char *msg)
{
	const float eps = 0.000001f;
	if (fabs(a-b) > eps)
	{
		diffs += msg + ret;
	}		
}

/*
 * Check integers.
 */
void CHPCompareUtils::CompareInts(int a, int b, string &diffs, const char *msg)
{
	if (a != b)
	{
		diffs += msg + ret;
	}
}

/*
 * Check strings.
 */
void CHPCompareUtils::CompareStrings(const string &a, const string &b, string &diffs, const char *msg)
{
	if (a != b)
	{
		diffs += msg + ret;
	}
}

/*
 * Check strings.
 */
void CHPCompareUtils::CompareStrings(const wstring &a, const wstring &b, string &diffs, const char *msg)
{
	if (a != b)
	{
		diffs += msg + ret;
	}
}

/*
 * Check strings.
 */
void CHPCompareUtils::CompareStrings(const string &a, const wstring &b, string &diffs, const char *msg)
{
	if (a != StringUtils::ConvertWCSToMBS(b))
	{
		diffs += msg + ret;
	}
}
