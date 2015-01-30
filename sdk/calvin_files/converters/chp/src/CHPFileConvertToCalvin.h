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


#ifndef _CHPFileConvertToCalvin_HEADER_
#define _CHPFileConvertToCalvin_HEADER_


/*! file CHPFileConvertToCalvin.h This file contains a class to convert a CHP file to Calvin format. */

//
#include "calvin_files/converters/chp/src/CHPFileConverterErrorCode.h"
#include "calvin_files/data/src/GenericDataHeader.h"
#include <calvin_files/data/src/CHPData.h>
//
#include "mas5-stat/src/ExpressionAlgorithmImplementation.h"
//


namespace affymetrix_chp_converter
{

/*! This class will convert a CHP file to Calvin format. */
class CHPFileConvertToCalvin
{
public:

	/*! Constructor */
	CHPFileConvertToCalvin();

	/*! Destructor */
	~CHPFileConvertToCalvin();

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
	CHPFileConverterErrorCode ErrorCode() const { return errorCode; }

	/*! Convert the XDA file to Calvin format.
	 * @param fileName The full path name of the CHP file.
	 * @param newFile The name of the new CHP file.
	 * @return True if successful.
	 */
	bool ConvertXDAFile(const char *fileName, const char *libPath, const char *newFile);

	/*! Convert a non-XDA file to Calvin format.
	 *	@param fileName The full name of the CHP file.
	 *	@param celFile The full name of the CEL file.
	 *	@param libPath The path of the library folder.
	 *	@param maskFile The full name of the mask file.  If there is no mask file, it can be NULL.
	 *	@param newFile The name of the new CHP file.
	 *	@return True if successful.
	 */
	bool ConvertMas5File(const char *fileName, const char *celFile, const char *libPath, const char* maskFile, const char *newFile);

	/*! Validate the conversion of an XDA file to Calvin format.
	 * @param fileName The full path name of the CHP file.
	 * @param newFile The name of the new CHP file.
	 * @return True if the validation was successful. If false is returned, call ErrorCode.
	 */
	bool ValidateXDAFile(const char *fileName, const char *libPath, const char *newFile);

	/*! Validate the conversion of a non-XDA file to Calvin format.
	 *	@param fileName The full name of the CHP file.
	 *	@param celFile The full name of the CEL file.
	 *	@param libPath The path of the library folder.
	 *	@param maskFile The full name of the mask file.  If there is no mask file, it can be NULL.
	 *	@param newFile The name of the new CHP file.
	 *	@return True if the validation was successful.  If false is returned, call ErrorCode.
	 */
	bool ValidateMas5File(const char *fileName, const char *celFile, const char *libPath, const char* maskFile, const char *newFile);

	/*! Compute the zone background for an expression array.
	 *	@param celData An open FusionCELData object.
	 *	@param libPath The path of the library folder.
	 *	@param zonesInfo Receives the results of the background computation.
	 *	@return True if the background was successully computed, otherwise false.
	 */
	bool ComputeBackgroundZoneInfo(FusionCELData& celData, const char* libPath, AllZonesInfoType& zonesInfo) const;

	/*! Compute the zone background for an expression array.
	 *	@param celData An open FusionCELData object.
	 *	@param libPath The path of the library folder.
	 *	@param maskPath The full name of the mask file.  If there is no mask file, it can be NULL.
	 *	@param zonesInfo Receives the results of the background computation.
	 *	@return True if the background was successully computed, otherwise false.
	 */
	bool ComputeBackgroundZoneInfo(FusionCELData& celData, const char* libPath, const char* maskPath, AllZonesInfoType& zonesInfo) const;

protected:

	/*! Fills the GenericDataHeader with the parent file header which it gets either from the parent file
	 *	or generates based on the GCOS CEL file.
	 *	@param hdr The GenericDataHeader that will be filled with the parent file header information.
	 *	@param inFile The GCOS CEL file.
	 */
	void FillParentGenericDataHeader(affymetrix_calvin_io::GenericDataHeader& hdr);

	/*! Gets the GenericDataHeader from the parent Command Console file.
	 *	The arrayID member is ignored.
	 *	@param hdr The GenericDataHeader that will be filled with the parent file information.
	 *	@return True if the parent file GenericDataHeader information was found.
	 */
	bool GetParentGenericDataHeader(affymetrix_calvin_io::GenericDataHeader& hdr);

	/*! Adds any extra parameters from the extraParameter collection to the file header.
	 *	@param data The target CHPData.
	 */
	void AddExtraParameters(affymetrix_calvin_io::CHPData& data);

	/*! Checks that conversion is supported based on the algorithm name.
	 *	@param algorithmName The algorithm name. ExpressionCall, ExpressionStat etc.
	 *	@return True if conversion of a CHP based on the indicated algorithm is supported.
	 */
	bool CheckAlgorithmType(const std::string& algorithmName);

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
	CHPFileConverterErrorCode errorCode;

	/*! Clear the class members. */
	void Clear();
};

}

#endif
