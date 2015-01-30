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


#ifndef _CHPConversionUtilities_HEADER_
#define _CHPConversionUtilities_HEADER_

/*! file CHPConversionUtilities.h This file contains utilities needed for converting CHP file formats. */

#include "calvin_files/data/src/CHPData.h"
//
#include "file/CHPFileWriter.h"
//
#include <cstring>
#include <string>
//

namespace affymetrix_chp_converter
{

/*! This class will copy Calvin data to GCOS CHP objects. */
class CHPConversionUtilities
{
public:

	/*! Copies the calvin object to the GCOS object.
	 * @param inFile The Calvin CHP file object.
	 * @param outFile The GCOS CHP file object.
	*/
	static void Copy(affymetrix_calvin_io::CHPData &inFile, affxchpwriter::CCHPFileWriter &outFile);

	/*! Sets the algorithm name in the CCHPFileWriter object.
	 *	It ensures that the algorithm name value is appropriate before setting it.
	 *	@param outFile The GCOS CHP file object.
	 *	@param algName The algorithm name.
	 */
	static void SetAlgName(affxchpwriter::CCHPFileWriter& outFile, std::string algName);

};

}

#endif
