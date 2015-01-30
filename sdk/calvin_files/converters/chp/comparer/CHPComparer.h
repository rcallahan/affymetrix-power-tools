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


#ifndef _CHPComparer_HEADER_
#define _CHPComparer_HEADER_

/*! \file CHPComparer.h Defines classes to compare GCOS and Calvin files. */

#include <cstring>
#include <string>
//

namespace affymetrix_comparer
{

/*! Provides for information about CHP file versions. */
class CHPComparer
{
private:
	/*! The differences between the files. */
	std::string differences;

public:

	/*! Clears the members. */
	void Clear();

	/*! Determines the version number of the file.
	* @param file1 The first file to compare.
	* @param file2 The second file to compare.
	* @return True if identical.
	*/
	bool CompareFiles(const char *file1, const char *file2);

	/*! The differences between the files.
	 * @return A string representation of the differences.
	 */
	std::string Differences() const { return differences; }
};

}

#endif

