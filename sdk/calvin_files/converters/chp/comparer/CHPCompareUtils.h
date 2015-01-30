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


#ifndef _CHPCompareUtils_HEADER_
#define _CHPCompareUtils_HEADER_

/*! \file CHPCompareUtils.h Defines utilities for the CHP comparison classes. */

#include <cstring>
#include <string>
//

namespace affymetrix_comparer
{

/*! Utilities for CHP file comparisons. */
class CHPCompareUtils
{
public:
	/*! Compares two floats.
	 * @param a The first value to compare.
	 * @param b The second value to compare.
	 * @param diffs A string to hold a message about the differences.
	 * @param msg The message to append to the differences.
	 */
	static void CompareFloats(float a, float b, std::string &diffs, const char *msg);

	/*! Compares two integers.
	 * @param a The first value to compare.
	 * @param b The second value to compare.
	 * @param diffs A string to hold a message about the differences.
	 * @param msg The message to append to the differences.
	 */
	static void CompareInts(int a, int b, std::string &diffs, const char *msg);

	/*! Compares two strings.
	 * @param a The first value to compare.
	 * @param b The second value to compare.
	 * @param diffs A string to hold a message about the differences.
	 * @param msg The message to append to the differences.
	 */
	static void CompareStrings(const std::string &a, const std::string &b, std::string &diffs, const char *msg);

	/*! Compares two strings.
	 * @param a The first value to compare.
	 * @param b The second value to compare.
	 * @param diffs A string to hold a message about the differences.
	 * @param msg The message to append to the differences.
	 */
	static void CompareStrings(const std::wstring &a, const std::wstring &b, std::string &diffs, const char *msg);

	/*! Compares two strings.
	 * @param a The first value to compare.
	 * @param b The second value to compare.
	 * @param diffs A string to hold a message about the differences.
	 * @param msg The message to append to the differences.
	 */
	static void CompareStrings(const std::string &a, const std::wstring &b, std::string &diffs, const char *msg);

};

}

#endif

