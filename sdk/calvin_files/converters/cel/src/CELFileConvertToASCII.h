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


#ifndef _CELFileConvertToASCII_HEADER_
#define _CELFileConvertToASCII_HEADER_

/*! file CELFileConvertToASCII.h This file contains a class to convert a CEL file to ASCII format. */

#include "calvin_files/converters/cel/src/CELFileConversionOptions.h"
#include "calvin_files/converters/cel/src/CELFileConverterErrorCode.h"
//

namespace affymetrix_cel_converter
{

/*! This class will convert a CEL file to ASCII format. */
class CELFileConvertToASCII
{
public:

	/*! Constructor */
	CELFileConvertToASCII();

	/*! Destructor */
	~CELFileConvertToASCII();

	/*! Gets the error code for the conversion process.
	 * @return The error code associated with the last failure.
	 */
	CELFileConverterErrorCode ErrorCode() const { return errorCode; }

	/*! Convert the Calvin file to ASCII format.
	 * @param fileName The full path name of the CEL file.
	 * @param newFile The name of the new CEL file.
	 * @return True if successful.
	 */
	bool ConvertCalvinFile(const char *fileName, const char *newFile, CELFileConversionOptions *options = NULL);

	/*! Convert the XDA file to ASCII format.
	 * @param fileName The full path name of the CEL file.
	 * @param newFile The name of the new CEL file.
	 * @return True if successful.
	 */
	bool ConvertXDAFile(const char *fileName, const char *newFile, CELFileConversionOptions *options = NULL);

private:

	/*! The error code. */
	CELFileConverterErrorCode errorCode;

	/*! Clear the class members. */
	void Clear();

};

}

#endif
