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


#ifndef _CHPFileConvertToXDA_HEADER_
#define _CHPFileConvertToXDA_HEADER_


/*! file CHPFileConvertToXDA.h This file contains a class to convert a CHP file to XDA format. */

#include "calvin_files/converters/chp/src/CHPFileConverterErrorCode.h"
#include <calvin_files/data/src/CHPData.h>
//

namespace affymetrix_chp_converter
{

/*! This class will convert a CHP file to XDA format. */
class CHPFileConvertToXDA
{
public:

	/*! Constructor */
	CHPFileConvertToXDA();

	/*! Destructor */
	~CHPFileConvertToXDA();

	/*! Set extra parameters that are to be written to the file.
	 *	@param params A collection of ParameterNameValueTypes to add to the file.
	 */
	void SetExtraParameters(const affymetrix_calvin_parameter::ParameterNameValueTypeVector* params) { extraParameters = params; }

	/*! Gets the error code for the conversion process.
	 * @return The error code associated with the last failure.
	 */
	CHPFileConverterErrorCode ErrorCode() const { return errorCode; }

	/*! The full path and name of the parent file.
	 *	@param str The full path and name of the parent file.
	 */
	void SetParentFileName(const char* str) { parentFileName = str; }

	/*! Convert the Calvin file to XDA format.
	 * @param fileName The full path name of the CHP file.
	 * @param newFile The name of the new CHP file.
	 * @return True if successful.
	 */
	bool ConvertCalvinFile(const char *fileName, const char *newFile);

private:

	/*! Extra parameters. Lifetime is not managed by this class. */
	const affymetrix_calvin_parameter::ParameterNameValueTypeVector* extraParameters;

	/*! The error code. */
	CHPFileConverterErrorCode errorCode;

	/*! The parent file name. */
	std::string parentFileName;

	/*! Clear the class members. */
	void Clear();

	/*! Get the array type from extraParameters
	 * @return The array type from extraParameter.  If the array type is not in the extraParameters collection
	 * an empty string is returned.
	 */
	std::string GetChipTypeFromExtraParameters();

};

}

#endif
