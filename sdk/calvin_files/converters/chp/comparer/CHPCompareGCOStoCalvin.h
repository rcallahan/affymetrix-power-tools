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


#ifndef _CHPCompareGCOStoCalvin_HEADER_
#define _CHPCompareGCOStoCalvin_HEADER_

/*! \file CHPCompareGCOStoCalvin.h Defines classes to compare GCOS and Calvin files. */

#include "calvin_files/data/src/CHPData.h"
//
#include "file/CHPFileData.h"
//
#include <cstring>
#include <string>
//

namespace affymetrix_comparer
{

/*! Provides for information about CHP file versions. */
class CHPCompareGCOStoCalvin
{
private:
	/*! The differences between the files. */
	std::string differences;

	/*! Compare the file header. */
	void CompareHeader();

	/*! Compare the header, not including alg parameters. */
	void CompareNonParameterHeader();

	/*! Compare the algorithm parameters. */
	void CompareAlgParams();

	/*! Compare the chip summary. */
	void CompareChipSummary();

	/*! Compare the background zone information. */
	void CompareBackground();

	/*! Compare the data. */
	void CompareData();

	/*! Compare the expression data. */
	void CompareExpression();

	/*! Compare the genotyping data. */
	void CompareGenotyping();

	/*! Compare the universal data. */
	void CompareUniversal();

	/*! Compare the resequencing data. */
	void CompareResequencing();

	/*! Compare the files. */
	void CompareFiles();

	/*! Open the files.
	 * @param gcosFile The file name of the gcos CHP file.
	 * @param calvinFile The file name of the calvin CHP file.
	 * @return True if successful.
	 */
	bool OpenFiles(const char *gcosFile, const char *calvinFile);

	/*! The GCOS file object. */
	affxchp::CCHPFileData gcos;

	/*! The calvin file object. */
	affymetrix_calvin_io::CHPData calvin;

public:

	/*! Clears the members. */
	void Clear();

	/*! Determines the version number of the file.
	* @param gcosFile The file name of the gcos CHP file.
	* @param calvinFile The file name of the calvin CHP file.
	* @return True if identical.
	*/
	bool CompareFiles(const char *gcosFile, const char *calvinFile);

	/*! The differences between the files.
	 * @return A string representation of the differences.
	 */
	std::string Differences() const { return differences; }
};

}

#endif

