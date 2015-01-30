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


#ifndef _CELConversionUtilities_HEADER_
#define _CELConversionUtilities_HEADER_

/*! file CELConversionUtilities.h This file contains utilities needed for converting CEL file formats. */

#include "calvin_files/converters/cel/src/CELFileConversionOptions.h"
#include "calvin_files/data/src/CELData.h"
//
#include "file/CELFileData.h"
//
#include <cstring>
#include <string>
//

namespace affymetrix_cel_converter
{

/*! This class will convert a CEL file from GCOS to Calvin format. */
class CELConversionUtilities
{
public:

	/*! Gets the grid coordinates from the calvin data object.
	 * @param inFile The Calvin CEL file object.
	 * @return The grid coordinates.
	 */
	static GridCoordinatesType GetGrid(affymetrix_calvin_io::CelFileData &inFile);

	/*! Gets the Dat header from the calvin data object.
	 * @param inFile The Calvin CEL file object.
	 * @return The Dat header string.
	*/
	static std::string GetDatHeader(affymetrix_calvin_io::CelFileData &inFile);

	/*! Copies the calvin object to the GCOS object.
	 * @param inFile The Calvin CEL file object.
	 * @param outFile The GCOS CEL file object.
	*/
	static bool Copy(affymetrix_calvin_io::CelFileData &inFile, affxcel::CCELFileData &outFile, CELFileConversionOptions *options = NULL);

	/*! Copies a GCOS CEL object to another GCOS CEL object.
	 * @param inFile The GCOS CEL file object to copy.
	 * @param outFile The resulting GCOS CEL file object.
	 */
	static bool Copy(affxcel::CCELFileData &inFile, affxcel::CCELFileData &outFile, CELFileConversionOptions *options = NULL);

	/*! Sets the cell margin in the outFile.  It converts the parameter is required.
	 *	@param inFile The GCOS CEL file holding the CellMargin value.
	 *	@param outFile The GCOS CEL file object to which the cell margin is set.
	 *	@param defaultValue The cell margin value to write if the value is not found in inFile.
	 */
	static void SetCellMargin(affymetrix_calvin_io::CelFileData &inFile, affxcel::CCELFileData &outFile, int defaultValue);

	/*! Finds the algorithm parameter value in the inFile and sets the value to the outFile algorithm parameters.
	 *	If the value is not found in inFile, then the default value is written to the outFile.
	 *	@param name The name of the parameter to find in the inFile.
	 *	@param defaultValue The value to use if the parameter is not found in inFile.
	 *	@param inFile The Calvin CEL file ojbect.
	 *	@param outFile The resulting GCOS CEL object.
	 */
	static void FindAndAddAlgorithmParameter(std::string name, std::string defaultValue, affymetrix_calvin_io::CelFileData &inFile, affxcel::CCELFileData &outFile);


};

}

#endif
