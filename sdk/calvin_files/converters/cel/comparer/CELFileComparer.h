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


#ifndef _CELFileComparer_HEADER_
#define _CELFileComparer_HEADER_

/*! \file CELFileComparer.h This file contains a class to compare a GCOS CEL file and Calvin format CEL. */

#include "calvin_files/fusion/src/FusionCELData.h"
//
#include "file/CELFileData.h"
//
#include <cstring>
#include <fstream>
#include <list>
#include <string>
//

using namespace std;

namespace affymetrix_comparers
{

/*! This class will convert a CEL file from GCOS to Calvin format. */
class CELFileComparer
{
public:
	/*! Error codes for the conversion process. */
	typedef enum _CompareErrorCodes
	{
		NoError,				/*! No error. */
		CalvinFileAsInputError,	/*! The input file is a Calvin CEL file. */
		UnableToOpenCelFile,	/*! Unable to open the input CEL file. */
		UnableToWriteTheFile,	/*! Unable to write the output file. */
		ComparisonFailed		/*! Comparison between Calvin and GCOS Cel files failed. */
	} CompareErrorCodes;

	/*! Constructor */
	CELFileComparer();
	CELFileComparer(bool newWriteDetailedDiffFile, float newTolerance, bool newIgnoreMask, bool newCelSummary, bool newIgnoreHeaders);

	/*! Destructor */
	~CELFileComparer();

	/*! Sets the name of the GCOS CEL file.
	 * @param name The full path name of the CEL file.
	 */
	void SetFileNameGCOS(const char *name) { fileNameGCOS = name; }

	/*! Sets the name of the Calvin CEL file.
	 * @param name The full path name of the CEL file.
	 */
	void SetFileNameCalvin(const char *name) { fileNameCalvin = name; }

	/*! Gets the name of the GCOS CEL file.
	 * @return The full path name of the CEL file.
	 */
	std::string GetFileNameGCOS() const { return fileNameGCOS; }

	/*! Gets the name of the GCOS CEL file.
	 * @return The full path name of the CEL file.
	 */
	std::string GetFileNameCalvin() const { return fileNameCalvin; }

	/*! Gets the error code for the conversion process.
	 * @return The error code associated with the last failure.
	 */
	CompareErrorCodes GetErrorCode() const { return errorCode; }

	/*! Convert the file to the Calvin format.
	 * @return True if successful.
	 */
	bool CompareFiles();

	/*! Get the output file name.
	 * @return The output file name.
	 */
	std::string GetOutputFileName() const { return outputFile; }

private:
	/*! The GCOS CEL file object. */
	affxcel::CCELFileData gcos;

	/*! The Calvin CEL file object. */
	affymetrix_fusion_io::FusionCELData fusion;

	ofstream streamCompare;

	/*! The GCOS CEL file name. */
	std::string fileNameGCOS;

	/*! The Calvin CEL file name. */
	std::string fileNameCalvin;

	/*! The name of the output CEL file. */
	std::string outputFile;

	/*! The error code. */
	CompareErrorCodes errorCode;
	
	/*! Reads the CEL file into a data object.
	 * @return True if successful.
	 */
	bool ReadGCOSFile();

	/*! Reads the CEL file into a data object.
	 * @return True if successful.
	 */
	bool ReadCalvinFile();

	/*! Closes the files. */
	void CloseFile();

	/*! Creates a Calvin file with the expression probe sets results.
	 * @return True if successful.
	 */
	bool Compare();

	/*! Helper function.  Trim excess spaces from the end of a string.	*/
	void trim(std::string& str);

	/*! Helper function.  Compare floating point values.	*/
	bool CompareFloats(float, float);

	/*! Detailed Diff file flag						*/
	bool detailedDiffFile;

	/*! Ignore Mask flag					*/
	bool ignoreMask;

	/*! celSummary flag						*/
	bool celSummary;

	/*! Float comparison tolerance					*/
	float float_tolerance;

	/*! Ignore Headers flag					*/
	bool ignoreHeaders;
};

}

#endif
