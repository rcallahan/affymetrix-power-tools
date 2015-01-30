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


#ifndef _CELFileConvertToCalvin_HEADER_
#define _CELFileConvertToCalvin_HEADER_

/*! file CELFileConvertToCalvin.h This file contains a class to convert a CEL file to Calvin format. */

#include "calvin_files/converters/cel/src/CELFileConversionOptions.h"
#include "calvin_files/converters/cel/src/CELFileConverterErrorCode.h"
#include "calvin_files/data/src/CELData.h"
//
#include "file/CELFileData.h"
//

namespace affymetrix_cel_converter
{

/*! This class will convert a CEL file to Calvin format. */
class CELFileConvertToCalvin
{
public:

	/*! Constructor */
	CELFileConvertToCalvin();

	/*! Destructor */
	~CELFileConvertToCalvin();

	/*! Set the physical array id.
	 *	@param str The physical array id.
	 */
	void SetArrayID(const char *str) { arrayID = str; }

	/*! Set the array barcode.
	 *	@param str The array barcode.
	 */
	void SetArrayBarcode(const wchar_t *str) { arrayBarcode = str; }

	/*! The full path and name of the parent file.
	 *	@param str The full path and name of the parent file.
	 */
	void SetParentFileName(const char* str) { parentFileName = str; }

	/*! Set extra parameters that are to be written to the file.
	 *	@param params A collection of ParameterNameValueTypes to add to the file.
	 */
	void SetExtraParameters(const affymetrix_calvin_parameter::ParameterNameValueTypeVector* params) { extraParameters = params; }

	/*! Gets the error code for the conversion process.
	 * @return The error code associated with the last failure.
	 */
	CELFileConverterErrorCode ErrorCode() const { return errorCode; }

	/*! Convert the file to Calvin format.
	 * @param fileName The full path name of the CEL file.
	 * @param newFile The name of the new CEL file.
	 * @return True if successful.
	 */
	bool ConvertASCIIFile(const char *fileName, const char *newFile, CELFileConversionOptions *options = NULL);

	/*! Convert the file to Calvin format.
	 * @param fileName The full path name of the CEL file.
	 * @param newFile The name of the new CEL file.
	 * @return True if successful.
	 */
	bool ConvertXDAFile(const char *fileName, const char *newFile, CELFileConversionOptions *options = NULL, const char* datFileName = NULL );

protected:

	/*! Fills the GenericDataHeader with the parent file header which it gets either from the parent file
	 *	or generates based on the GCOS CEL file.
	 *	@param hdr The GenericDataHeader that will be filled with the parent file header information.
	 *	@param inFile The GCOS CEL file.
	 */
	bool FillParentGenericDataHeader(affymetrix_calvin_io::GenericDataHeader& hdr, affxcel::CCELFileData& inFile);

	/*! Gets the GenericDataHeader from the parent Command Console file.
	 *	The arrayID member is ignored.
	 *	@param hdr The GenericDataHeader that will be filled with the parent file information.
	 *	@return True if the parent file GenericDataHeader information was found.
	 */
	bool GetParentGenericDataHeader(affymetrix_calvin_io::GenericDataHeader& hdr);

	/*! Generates a parent GenericDataHeader based on the GCOS CEL file.
	 *	If the arrayID member is set, it will be added to the parent header.
	 *	@param hdr The GenericDataHeader that will be filled with the parent information.
	 */
	bool GenerateParentGenericDataHeaderFromGCOS(affymetrix_calvin_io::GenericDataHeader& hdr, affxcel::CCELFileData& inFile);

	/*! Adds any extra parameters from the extraParameter collection to the file header.
	 *	@param data The target CelFileData.
	 */
	void AddExtraParameters(affymetrix_calvin_io::CelFileData& data);

	/*! Gets the scanner type from the extra parameters.
	 *	@return Scanner type
	 */
	std::wstring GetScannerTypeFromExtraParameters();

private:

	/*! The physical array ID. */
	std::string arrayID;

	/*! The array barcode */
	std::wstring arrayBarcode;

	/*! The parent file name. */
	std::string parentFileName;

	/*! Extra parameters. Lifetime is not managed by this class. */
	const affymetrix_calvin_parameter::ParameterNameValueTypeVector* extraParameters;

	/*! The error code. */
	CELFileConverterErrorCode errorCode;

	/*! Clear the class members. */
	void Clear();

    /*! Options for conversion. */
    CELFileConversionOptions *m_Options;

};

}

#endif
