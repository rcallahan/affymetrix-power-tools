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


#ifndef _ConverterFileUtils_HEADER_
#define _ConverterFileUtils_HEADER_

/*! \file ConverterFileUtils.h This file contains file system utility functions. */

#include "calvin_files/data/src/GenericDataHeader.h"
//

/*! Moves a file.
 * @param oldFile The name of the file to move.
 * @param newFile The name of the new file.
 * @param overwrite A flag indicating if an existing file should be deleted prior to the move.
 * @return True if the file was moved.
 */
bool ConverterMoveFile(const char *oldFile, const char *newFile, bool overwrite=false);

/*! Copies a file.
 * @param oldFile The name of the file to copy.
 * @param newFile The name of the new file.
 * @param overwrite A flag indicating if an existing file should be deleted prior to the copy.
 * @return True if the file was copied.
 */
bool ConverterCopyFile(const char *oldFile, const char *newFile, bool overwrite=false);

/*! Deletes a file.
 * @param fileName The name of the file to delete.
 * @return True if the file was deleted.
 */
bool ConverterRemoveFile(const char *fileName);

/*! Checks if a file exists.
 * @param fileName The name of the file to delete.
 * @return True if exists.
 */
bool ConverterFileExists(const char *fileName);

/*! Creates and adds the array set GenericDataHeader as a parent to the hdr parameter.
 *	The arrayID member will be set as the affymetrix-array-id value.
 *	@param hdr The GenericDataHeader that will receive the array set parent header.
 *	@param arrayID. The array ID. (physcial array ID).
 *	@param arrayBarcode.  The array barcode.
 */
void ConverterCreateAndAddParentArrayGenericDataHeader(affymetrix_calvin_io::GenericDataHeader& hdr, std::string arrayID, std::wstring arrayBarcode);


#endif

