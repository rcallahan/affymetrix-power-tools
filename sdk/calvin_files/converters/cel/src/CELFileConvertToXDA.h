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


#ifndef _CELFileConvertToXDA_HEADER_
#define _CELFileConvertToXDA_HEADER_

/*! file CELFileConvertToXDA.h This file contains a class to convert a CEL file to XDA format. */

#include "calvin_files/converters/cel/src/CELFileConversionOptions.h"
#include "calvin_files/converters/cel/src/CELFileConverterErrorCode.h"
#include "calvin_files/data/src/CELData.h"
//

namespace affymetrix_cel_converter
{

/*! This class will convert a CEL file to XDA format. */
class CELFileConvertToXDA
{
public:

	/*! Constructor */
	CELFileConvertToXDA();

	/*! Destructor */
	~CELFileConvertToXDA();

	/*! Set extra parameters that are to be written to the file.
	 *	@param params A collection of ParameterNameValueTypes to add to the file.
	 */
	void SetExtraParameters(const affymetrix_calvin_parameter::ParameterNameValueTypeVector* params) { extraParameters = params; }

	/*! Gets the error code for the conversion process.
	 * @return The error code associated with the last failure.
	 */
	CELFileConverterErrorCode ErrorCode() const { return errorCode; }

	/*! Convert the Calvin file to XDA format.
	 * @param fileName The full path name of the CEL file.
	 * @param newFile The name of the new CEL file.
	 * @return True if successful.
	 */
	bool ConvertCalvinFile(const char *fileName, const char *newFile, CELFileConversionOptions *options = NULL);

	/*! Convert the ASCII file to XDA format.
	 * @param fileName The full path name of the CEL file.
	 * @param newFile The name of the new CEL file.
	 * @return True if successful.
	 */
	bool ConvertASCIIFile(const char *fileName, const char *newFile, CELFileConversionOptions *options = NULL);

private:

	/*! Extra parameters. Lifetime is not managed by this class. */
	const affymetrix_calvin_parameter::ParameterNameValueTypeVector* extraParameters;

	/*! The error code. */
	CELFileConverterErrorCode errorCode;

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
