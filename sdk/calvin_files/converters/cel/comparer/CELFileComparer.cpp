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


#include "calvin_files/converters/cel/comparer/CELFileComparer.h"
//
#include "calvin_files/parsers/src/CelFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include "util/Fs.h"
//
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>
//

using namespace std;
using namespace affxcel;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affymetrix_comparers;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

/*
 * Initialize the class.
 */
CELFileComparer::CELFileComparer(bool newDetailedDiffFile, float newTolerance, bool newIgnoreMask, bool newCelSummary, bool newIgnoreHeaders)
{
	detailedDiffFile = newDetailedDiffFile;
	float_tolerance = newTolerance;
	ignoreMask = newIgnoreMask;
	ignoreHeaders = newIgnoreHeaders;
	celSummary = newCelSummary;
}

/*
 * Deallocate any memory.
 */
CELFileComparer::~CELFileComparer()
{
	CloseFile();
}

/*
 * Closes the open files.
 */
void CELFileComparer::CloseFile()
{
	gcos.Clear();
	streamCompare.close();
}

/*
 * Read the CEL file into a data object.
 */
bool CELFileComparer::ReadGCOSFile()
{
	gcos.SetFileName(fileNameGCOS.c_str());
	bool status = gcos.Read();
	if (status == false)
		errorCode = UnableToOpenCelFile;
	return status;
}

/*
 * Read the CEL file into a data object.
 */
bool CELFileComparer::ReadCalvinFile()
{
	//fusion.SetFileName();
	fusion.ReadEx(fileNameCalvin.c_str(),FusionCELData::CEL_ALL);
	return true;
}

/*
 * Convert the file.
 *
 * First check the validity of the file.
 *	Is it a calvin file? 
 *
 * Read the GCOS file.
 *
 */
bool CELFileComparer::CompareFiles()
{
	errorCode = NoError;

	// Read the GCOS file.
	if (ReadGCOSFile() == false)
		return false;

	// Read the GCOS file.
	ReadCalvinFile();

	// Clear the output file list
	outputFile = fileNameCalvin + ".comparison";

	// Create a Calvin CEL file.
	bool status = Compare();

	// Close the file and return the status.
	CloseFile();

	return status;
}

/*
 * Create a Calvin CEL file.
 *
 * Issues:
 *	DatHeader is missing.
 *	Grid is missing.
 *	Intensities should be float's not int's - do we store them in the GCOS 1 digit precision?
 *	There is a bug in writing the outlier and mask sections when no data exists.
 *
 */
bool CELFileComparer::Compare()
{
	bool compareResults = true;

        Fs::aptOpen(streamCompare, outputFile);
	time_t ltime;
	time(&ltime );

#ifdef WIN32
	wchar_t buf[26];
	_wctime_s(buf, 26, &ltime);
	streamCompare << "Time stamp: " << StringUtils::ConvertWCSToMBS(wstring(buf)) << endl << endl;
#endif

	///---------------------
	/// Commented out by Breck (20060712) for the header varaible was not used. The original purpose was to 
	/// compare the gocs header to the calvin gcos dat header that was inserted. Not sure if the 
	/// gcos header ever made it into the calvin header. If it does, uncomment this code and add 
	///the call to the calvin gcos data header and compare.
	///CCELFileHeaderData &gcosHeader = gcos.GetHeader();
	///---------------------

	CelFileReader cr;
	CelFileData data;
	cr.SetFilename(fileNameCalvin.c_str());
	cr.Read(data);

	if (!ignoreHeaders)
	{
		if(fusion.GetAlg().compare(StringUtils::ConvertMBSToWCS(gcos.GetAlg())) != 0)
		{
			//Algorithm names need not match, do not flag as a failure
			streamCompare << "Algorithm Names mismatch: GCOS = " << gcos.GetAlg() << "\tCalvin = " << StringUtils::ConvertWCSToMBS(fusion.GetAlg()) << endl << endl;
		}

		if(int32_t(gcos.GetCols()) != fusion.GetCols())
		{
			compareResults = false;
			streamCompare << "Columns: GCOS = " << gcos.GetCols() << "\tCalvin = " << fusion.GetCols() << endl << endl;
		}

		if(int32_t(gcos.GetRows()) != fusion.GetRows())
		{
			compareResults = false;
			streamCompare << "Rows: GCOS = " << gcos.GetRows() << "\tCalvin = " << fusion.GetRows() << endl << endl;
		}

		if(fusion.GetChipType().compare(StringUtils::ConvertMBSToWCS(gcos.GetChipType())) != 0)
		{
			compareResults = false;
			streamCompare << "Chip Type: GCOS = " << gcos.GetChipType() << "\tCalvin = " << StringUtils::ConvertWCSToMBS(fusion.GetChipType()) << endl << endl;
		}

		if(fusion.GetCellMargin() != gcos.GetCellMargin())
		{
			compareResults = false;
			streamCompare << "CellMargin: GCOS = " << gcos.GetCellMargin() << "\tCalvin = " << fusion.GetCellMargin() << endl << endl;
		}

		if(fusion.GetNumCells() != gcos.GetNumCells())
		{
			compareResults = false;
			streamCompare << "NumCells: GCOS = " << gcos.GetNumCells() << "\tCalvin = " << fusion.GetNumCells() << endl << endl;
		}
	}

	if (!ignoreMask)
	{
		if(fusion.GetNumMasked() != gcos.GetNumMasked())
		{
			{
				compareResults = false;
				streamCompare << "NumMasked: GCOS = " << gcos.GetNumMasked() << "\tCalvin = " << fusion.GetNumMasked() << endl << endl;
			}
		}
	}

	//If comparing CelSummary data, compare mask total to celSummary total
	if (celSummary)
	{
		std::string dgName("Default Group");
		std::string dsName("CelSummaryReport");

		//Data Comparison
		DataSet* ds1=NULL;
		bool celSummaryDsExists = true;

		try
		{
			ds1 = fusion.GetGenericData()->DataSet(StringUtils::ConvertMBSToWCS(dgName), StringUtils::ConvertMBSToWCS(dsName));
		}
		catch (affymetrix_calvin_exceptions::DataSetNotFoundException)
		{
			streamCompare << "CelSummaryReport data set was not found in the CEL file.  Use -cs option only if comparing CEL files with this dataset."  << endl << endl;
			compareResults = false;
			celSummaryDsExists = false;
		}

		if (celSummaryDsExists)
		{
			ds1->Open();

			int celSumRows = ds1->Rows();

			if(celSumRows != gcos.GetNumMasked())
			{
				{
					compareResults = false;
					streamCompare << "Mismatch in number of CelSummary points vs. GCOS masked points."  << endl << endl;
					streamCompare << "NumMasked: GCOS = " << gcos.GetNumMasked() << "\t NumCelSummaryReport = " << celSumRows << endl << endl;
				}
			}

			// For each value in the CelSummary section, compare GCOS to see that the corresponding location is masked.  If not, it fails.
			// If all the corresponding locations are masked, and the counts match, the result is a pass.

			for (int i=0; i < celSumRows; i++)
			{
				int32_t xLocation;
				int32_t yLocation;

				try
				{
					ds1->GetData(i, 0, xLocation);
					ds1->GetData(i, 1, yLocation);
				}
				catch (...)
				{
					compareResults = false;
					streamCompare << "Error reading Calvin CelSummary data group."  << endl << endl;
					break;
				}

				//If location isn't masked in GCOS, it fails
				if (!gcos.IsMasked(xLocation, yLocation))
				{
					compareResults = false;
					streamCompare << "CelSummary point not masked in GCOS: (" << xLocation << ", " << yLocation << ")" << endl;
				}

			}

			//Add endl to make comparison file easier to read
			if (!compareResults)
				streamCompare << endl;
		}
	}

	if(fusion.GetNumOutliers() != gcos.GetNumOutliers())
	{
		compareResults = false;
		streamCompare << "NumOutliers: GCOS = " << gcos.GetNumOutliers() << "\tCalvin = " << fusion.GetNumOutliers() << endl << endl;
	}


	//if(fusion.GetVersion() != gcos.GetVersion())
	////////{
	////////	compareResults = false;
	////////	streamCompare << "Version: GCOS = " << gcos.GetVersion() << "\tCalvin = " << fusion.GetVersion() << endl << endl;
	////////}

	/// ------------------------
	/// Commented out by Breck (20060712)- This was not doing anything and am not sure why it is here. I do not remember what the original purpose was for.
	/// It may have been for a feature I was trying to get into the converters so that we knew what GCOS version it was converted from so that if 
	/// we had converted to Calvin and then back to gcos, we would know the proper format for gcos as well as data allowed. Not sure if the converter
	/// ever did it.
	//ParameterNameValueType value, ;

	//data.GetFileHeader()->GetGenericDataHdr()->FindNameValParam(L"ConvertedSourceVersion",value); 
	//ParameterNameValueType::ParameterType type = value.GetParameterType();
	//data.GetFileHeader()->GetGenericDataHdr()->FindNameValParam(L"ConverterVersion",value);
	///-------------------------


	/* Compare Algorithm parameters
	 *
	 * Step through each Calvin algorithm parameter, if it exists in GCOS,
	 * determine the data type from Calvin, cast each value, and compare the
	 * results.
	 *
	 * Keep a count of the number of matching parameters, this number is then
	 * compared to the number of GCOS parameters, offset by one to
	 * account for the missing parameter AlgVersion.
	 * 
	 */

	//Compare Headers only if ignoreHeaders is false

	if (!ignoreHeaders)
	{
		//Get Calvin Algorithm Parameters
		FusionTagValuePairTypeList& fusionList = fusion.GetParameters();
		std::list<FusionTagValuePairType>::iterator fi = fusionList.begin();

		//Expected number of matching keys
		int gcosParamCount = gcos.GetNumberAlgorithmParameters() - 1;

		for (int j=0; j < fusion.GetNumberAlgorithmParameters(); j++)
		{
			std::wstring key = fi->Tag;
			std::wstring value = fi->Value;

			std::string gcosKey = StringUtils::ConvertWCSToMBS(key);
			
			//Exception to map IgnoreShiftRowOutliers
			if (key.compare(L"IgnoreShiftRowOutliers") == 0)
				gcosKey = "IgnoreOutliersInShiftRows";

			//Exception to map PoolWidthExtenstion (GCOS) to ExtendPoolWidth (Calvin)
			if (key.compare(L"ExtendPoolWidth") == 0)
				gcosKey = "PoolWidthExtenstion";

			//Exception to map PoolHeightExtenstion (GCOS) to ExtendPoolHeight (Calvin)
			if (key.compare(L"ExtendPoolHeight") == 0)
				gcosKey = "PoolHeightExtension";

			std::string gcosValue = gcos.GetAlgorithmParameter(gcosKey.c_str());

			//If GCOS parameter exists, compare values, decrement count
			if (gcosValue.compare("") != 0)
			{
				gcosParamCount--;

				//Based on type of Calvin data, make comparison

				//ASCII
				if (fi->DetailedType().GetParameterType() == ParameterNameValueType::AsciiType)
				{			
					//Compare
					if (gcosValue.compare(StringUtils::ConvertWCSToMBS(value)) != 0)
					{
						streamCompare << "Mismatched ASCII Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
						compareResults = false;
					}
				}
				//Float
				else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::FloatType)
				{
					float gcosFloat = (float) atof(gcosValue.c_str());
					float calvinFloat = (float) atof( (StringUtils::ConvertWCSToMBS(value)).c_str());

					//Compare floats (exactly)
					if ( gcosFloat != calvinFloat )
					{
						streamCompare << "Mismatched Float Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
						compareResults = false;
					}
				}
				//Int16
				else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::Int16Type)
				{
					int16_t gcosInt16 = (int16_t) ( atoi(gcosValue.c_str()) );
					int16_t calvinInt16 = (int16_t) (atoi( StringUtils::ConvertWCSToMBS(value).c_str() ) );

					//Compare
					if (gcosInt16 != calvinInt16)
					{
						streamCompare << "Mismatched Int16 Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
						compareResults = false;
					}
				}
				//Int32
				else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::Int32Type)
				{
					int32_t gcosInt32 = (int32_t) ( atoi(gcosValue.c_str()) );
					int32_t calvinInt32 = (int32_t) (atoi( StringUtils::ConvertWCSToMBS(value).c_str() ) );

					//Compare
					if (gcosInt32 != calvinInt32)
					{
						//For PoolHeightExtension, PoolWidthExtension parameters, a mismatch is expected, so no error should be thrown
						if ( (gcosKey.compare("PoolWidthExtenstion") != 0) &&  (gcosKey.compare("PoolHeightExtension") != 0) )
						{
							streamCompare << "Mismatched Int32 Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
							compareResults = false;
						}
						else
						{
							streamCompare << "Mismatched Int32 Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
							streamCompare << "Mismatch for PoolHeightExtension, PoolHeightExtension not flagged as errors, mismatches are expected." << endl << endl;
						}
					}
				}
				//Int8
				else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::Int8Type)
				{
					int8_t gcosInt8 = (int8_t) ( atoi(gcosValue.c_str()) );
					int8_t calvinInt8 = (int8_t) (atoi( StringUtils::ConvertWCSToMBS(value).c_str() ) );

					//Compare
					if (gcosInt8 != calvinInt8)
					{
						streamCompare << "Mismatched Int8 Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
						compareResults = false;
					}					
				}
				//Text
				else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::TextType)
				{
					//Compare as string
					//If key is AlgVersion or FeatureExtraction, as mismatch is expected.  Log the difference, but no error
					if (gcosValue.compare(StringUtils::ConvertWCSToMBS(value)) != 0)
					{
						if ( (gcosKey.compare("AlgVersion") != 0) &&  (gcosKey.compare("FeatureExtraction") != 0) )
						{
							streamCompare << "Mismatched Text Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
							compareResults = false;
						}
						else
						{
							streamCompare << "Mismatched Text Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl;
							streamCompare << "Mismatch for AlgVersion, FeatureExtraction not flagged as errors, mismatches are expected." << endl << endl;
						}
					}					
				}
				//UInt8
				else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::UInt8Type)
				{
					uint8_t gcosUInt8 = (uint8_t) ( atoi(gcosValue.c_str()) );
					uint8_t calvinUInt8 = (uint8_t) (atoi( StringUtils::ConvertWCSToMBS(value).c_str() ) );

					//Compare
					if (gcosUInt8 != calvinUInt8)
					{
						streamCompare << "Mismatched UInt8 Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
						compareResults = false;
					}						
				}
				//UInt16
				else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::UInt16Type)
				{
					uint16_t gcosUInt16 = (uint16_t) ( atoi(gcosValue.c_str()) );
					uint16_t calvinUInt16 = (uint16_t) (atoi( StringUtils::ConvertWCSToMBS(value).c_str() ) );

					//Compare
					if (gcosUInt16 != calvinUInt16)
					{
						streamCompare << "Mismatched UInt16 Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
						compareResults = false;
					}				
				}
				//UInt32
				else if (fi->DetailedType().GetParameterType() == ParameterNameValueType::UInt32Type)
				{
					uint32_t gcosUInt32 = (uint32_t) ( atoi(gcosValue.c_str()) );
					uint32_t calvinUInt32 = (uint32_t) (atoi( StringUtils::ConvertWCSToMBS(value).c_str() ) );

					//Compare
					if (gcosUInt32 != calvinUInt32)
					{
						streamCompare << "Mismatched UInt32 Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
						compareResults = false;
					}					
				}
				//Default, compare as string
				else
				{
					//Compare
					if (gcosValue.compare(StringUtils::ConvertWCSToMBS(value)) != 0)
					{
						streamCompare << "Mismatched ASCII Parameter value by key(" << gcosKey << "): GCOS = " << gcosValue << "\tCalvin = " << StringUtils::ConvertWCSToMBS(value) << endl << endl;
						compareResults = false;
					}					
				}


			}

			fi++;
		}

		// Check matching parameter count, count should have decremented to 0 if all parameters were found
		//
		// If totals do not match, iterate through GCOS parameters and list those that are not present
		if (gcosParamCount > 0)
		{
			for (int i = 0; i < gcos.GetNumberAlgorithmParameters(); i++)
			{
				std::wstring calvinKey = StringUtils::ConvertMBSToWCS(gcos.GetAlgorithmParameterTag(i));
				
				//If Calvin does not contain the parameter, report it
				//Ignore AlgVersion, IgnoreShiftRowOutliers, PoolWidthExtenstion, PoolHeightExtenstion, NumDATSubgrids
				if ( (calvinKey.compare(StringUtils::ConvertMBSToWCS("AlgVersion")) != 0)	&& (calvinKey.compare(StringUtils::ConvertMBSToWCS("IgnoreOutliersInShiftRows")) != 0)
					&& (calvinKey.compare(StringUtils::ConvertMBSToWCS("PoolWidthExtenstion")) != 0) && (calvinKey.compare(StringUtils::ConvertMBSToWCS("PoolHeightExtension")) != 0)
					&& (calvinKey.compare(StringUtils::ConvertMBSToWCS("NumDATSubgrids")) != 0))
				{
						if (fusion.GetAlgorithmParameter(calvinKey.c_str()).compare(StringUtils::ConvertMBSToWCS("")) == 0)
						{
							streamCompare << "GCOS parameter " << StringUtils::ConvertWCSToMBS(calvinKey) << " not found in Calvin Cel file." << endl;
							compareResults = false;
						}
				}
			}
		}

		// check to see if we can get the grid coords from the parameter array directly
	#ifdef WIN32
		if(fusion.GetNumberAlgorithmParameters())
		{
			GridCoordinatesType gcosGridType = gcos.GetGridCorners();
			affymetrix_fusion_io::FGridCoords fusionGridType;
			wstring value = fusion.GetAlgorithmParameter(L"GridULX");
			fusionGridType.upperleft.x = (float)_wtof(value.c_str());
			value = fusion.GetAlgorithmParameter(L"GridULY");
			fusionGridType.upperleft.y = (float)_wtof(value.c_str());
			value = fusion.GetAlgorithmParameter(L"GridURX");
			fusionGridType.upperright.x = (float)_wtof(value.c_str());
			value = fusion.GetAlgorithmParameter(L"GridURY");
			fusionGridType.upperright.y = (float)_wtof(value.c_str());
			value = fusion.GetAlgorithmParameter(L"GridLRX");
			fusionGridType.lowerright.x = (float)_wtof(value.c_str());
			value = fusion.GetAlgorithmParameter(L"GridLRY");
			fusionGridType.lowerright.y = (float)_wtof(value.c_str());
			value = fusion.GetAlgorithmParameter(L"GridLLX");
			fusionGridType.lowerleft.x = (float)_wtof(value.c_str());
			value = fusion.GetAlgorithmParameter(L"GridLLY");
			fusionGridType.lowerleft.y = (float)_wtof(value.c_str());

			if((int32_t)fusionGridType.upperleft.x != gcosGridType.upperleft.x   ||
			(int32_t)fusionGridType.upperleft.y != gcosGridType.upperleft.y   ||
			(int32_t)fusionGridType.lowerleft.x != gcosGridType.lowerleft.x   ||
			(int32_t)fusionGridType.lowerleft.y != gcosGridType.lowerleft.y   ||
			(int32_t)fusionGridType.upperright.x != gcosGridType.upperright.x ||
			(int32_t)fusionGridType.upperright.y != gcosGridType.upperright.y ||
			(int32_t)fusionGridType.lowerright.x != gcosGridType.lowerright.x ||
			(int32_t)fusionGridType.lowerright.y != gcosGridType.lowerright.y)
			{
				streamCompare << "GridCoordinatesType:" << endl;
				streamCompare << "\tGCOS lowerleft.x	= " << gcosGridType.lowerleft.x << "\t\tCalvin lowerleft.x		= " << fusionGridType.lowerleft.x  << endl;
				streamCompare << "\tGCOS lowerleft.y	= " << gcosGridType.lowerleft.y << "\t\tCalvin lowerleft.y		= " << fusionGridType.lowerleft.y  << endl;
				streamCompare << "\tGCOS lowerright.x	= " << gcosGridType.lowerright.x << "\t\tCalvin lowerright.x	= " << fusionGridType.lowerright.x << endl;
				streamCompare << "\tGCOS lowerright.y	= " << gcosGridType.lowerright.y << "\t\tCalvin lowerright.y	= " << fusionGridType.lowerright.y << endl;
				streamCompare << "\tGCOS upperleft.x	= " << gcosGridType.upperleft.x << "\t\tCalvin upperleft.x		= " << fusionGridType.upperleft.x  << endl;
				streamCompare << "\tGCOS upperrleft.y	= " << gcosGridType.upperleft.y << "\t\tCalvin upperrleft.y		= " << fusionGridType.upperleft.y  << endl;
				streamCompare << "\tGCOS upperright.x	= " << gcosGridType.upperright.x << "\t\tCalvin upperright.x	= " << fusionGridType.upperright.x << endl;
				streamCompare << "\tGCOS upperright.y	= " << gcosGridType.upperright.y << "\t\tCalvin upperright.y	= " << fusionGridType.upperright.y << endl << endl;
			}
		}
	#endif

		if(StringUtils::ConvertMBSToWCS(gcos.GetDatHeader()).compare(fusion.GetDatHeader()) != 0)
		{
			//DAT headers need not match, do not flag as a failure
			streamCompare << "Dat Header String: GCOS = " << gcos.GetDatHeader() << endl << "Calvin = " << StringUtils::ConvertWCSToMBS(fusion.GetDatHeader()) << endl << endl;
		}

		/*
		 * Compare DAT header strings between GCOS and Calvin files.
		 * 
		 * Header strings are comprised of multiple parts, so the string must be broken down and compared
		 * piece by piece, with exceptions made for those segments that will not match due to designed in
		 * differences in the formats (e.g. filenames, algorithm names may not match).
		 */

		std::string gcosDATHeader = gcos.GetDatHeader();
		std::string calvinDATHeader = StringUtils::ConvertWCSToMBS(fusion.GetDatHeader());

		//Try/Catch block around all comparisons; any exceptions indicate failure

		try
		{
			//Compare min pixel values
			int gcosMinPixel = atoi(gcosDATHeader.substr(gcosDATHeader.find_first_of("[") + 1, gcosDATHeader.find_first_of("..")).c_str());
			int calvinMinPixel = atoi(calvinDATHeader.substr(calvinDATHeader.find_first_of("[") + 1, calvinDATHeader.find_first_of("..")).c_str());

			if(gcosMinPixel != calvinMinPixel)
			{
				compareResults = false;
				streamCompare << "Minimum Pixel values mismatch: GCOS = " << gcosMinPixel << endl << "Calvin = " << calvinMinPixel << endl << endl;
			}

			//Compare max pixel values
			int gcosMaxPixel = atoi(gcosDATHeader.substr(gcosDATHeader.find("..") + 2, gcosDATHeader.find("]")).c_str());
			int calvinMaxPixel = atoi(calvinDATHeader.substr(calvinDATHeader.find("..") + 2, calvinDATHeader.find("]")).c_str());

			if(gcosMaxPixel != calvinMaxPixel)
			{
				compareResults = false;
				streamCompare << "Maximum Pixel values mismatch: GCOS = " << gcosMaxPixel << endl << "Calvin = " << calvinMaxPixel << endl << endl;
			}
			
			//Ignore DAT file names, as they will not necessarily match

			//Remove string start up to columns
			gcosDATHeader = gcosDATHeader.substr(gcosDATHeader.find("CLS="), gcosDATHeader.length() - gcosDATHeader.find("CLS=") );
			calvinDATHeader = calvinDATHeader.substr(calvinDATHeader.find("CLS="), calvinDATHeader.length() - calvinDATHeader.find("CLS="));

			//Compare columns values
			int gcosColumns = atoi(gcosDATHeader.substr(4, gcosDATHeader.find("RWS=") - 4).c_str());
			int calvinColumns = atoi(calvinDATHeader.substr(4, calvinDATHeader.find("RWS=") - 4).c_str());

			if(gcosColumns != calvinColumns)
			{
				compareResults = false;
				streamCompare << "Columns values mismatch: GCOS = " << gcosColumns << endl << "Calvin = " << calvinColumns << endl << endl;
			}

			//Remove string start up to rows
			gcosDATHeader = gcosDATHeader.substr(gcosDATHeader.find("RWS="), gcosDATHeader.length() - gcosDATHeader.find("RWS=") );
			calvinDATHeader = calvinDATHeader.substr(calvinDATHeader.find("RWS="), calvinDATHeader.length() - calvinDATHeader.find("RWS="));

			//Compare rows values
			int gcosRows = atoi(gcosDATHeader.substr(4, gcosDATHeader.find("XIN=") - 4).c_str());
			int calvinRows = atoi(calvinDATHeader.substr(4, calvinDATHeader.find("XIN=") - 4).c_str());

			if(gcosRows != calvinRows)
			{
				compareResults = false;
				streamCompare << "Rows values mismatch: GCOS = " << gcosRows << endl << "Calvin = " << calvinRows << endl << endl;
			}

			//Remove string start up to XIN (x-pixel size)
			gcosDATHeader = gcosDATHeader.substr(gcosDATHeader.find("XIN="), gcosDATHeader.length() - gcosDATHeader.find("XIN=") );
			calvinDATHeader = calvinDATHeader.substr(calvinDATHeader.find("XIN="), calvinDATHeader.length() - calvinDATHeader.find("XIN="));

			//Compare Xin values
			int gcosXin = atoi(gcosDATHeader.substr(4, gcosDATHeader.find("YIN=") - 4).c_str());
			int calvinXin = atoi(calvinDATHeader.substr(4, calvinDATHeader.find("YIN=") - 4).c_str());

			if(gcosColumns != calvinColumns)
			{
				compareResults = false;
				streamCompare << "XIN Pixel size values mismatch: GCOS = " << gcosXin << endl << "calvinXin = " << calvinXin << endl << endl;
			}

			//Remove string start up to YIN (y-pixel size)
			gcosDATHeader = gcosDATHeader.substr(gcosDATHeader.find("YIN="), gcosDATHeader.length() - gcosDATHeader.find("YIN=") );
			calvinDATHeader = calvinDATHeader.substr(calvinDATHeader.find("YIN="), calvinDATHeader.length() - calvinDATHeader.find("YIN="));

			//Compare Yin values
			int gcosYin = atoi(gcosDATHeader.substr(4, gcosDATHeader.find("VE=") - 4).c_str());
			int calvinYin = atoi(calvinDATHeader.substr(4, calvinDATHeader.find("VE=") - 4).c_str());

			if(gcosYin != calvinYin)
			{
				compareResults = false;
				streamCompare << "YIN Pixel size values mismatch: GCOS = " << gcosYin << endl << "Calvin = " << calvinYin << endl << endl;
			}

			//Remove string start up to VE (scan speed)
			//These need not match
			gcosDATHeader = gcosDATHeader.substr(gcosDATHeader.find("VE="), gcosDATHeader.length() - gcosDATHeader.find("VE=") );
			calvinDATHeader = calvinDATHeader.substr(calvinDATHeader.find("VE="), calvinDATHeader.length() - calvinDATHeader.find("VE="));

			//Compare VE values
			int gcosVE= atoi(gcosDATHeader.substr(3, gcosDATHeader.find("VE=") - 4).c_str());
			int calvinVE = atoi(calvinDATHeader.substr(3, calvinDATHeader.find("VE=") - 4).c_str());

			if(gcosVE != calvinVE)
			{
				//Values need not match, do not flag as a failure
				streamCompare << "VE (scan speed) values mismatch: " << endl << "GCOS VE= " << gcosVE << endl << "Calvin VE = " << calvinVE << endl << endl;
			}

			//Trim VE data from headers
			gcosDATHeader = gcosDATHeader.substr(gcosDATHeader.find(" "), gcosDATHeader.length() - gcosDATHeader.find(" ") );
			calvinDATHeader = calvinDATHeader.substr(calvinDATHeader.find(" "), calvinDATHeader.length() - calvinDATHeader.find(" "));
			trim(gcosDATHeader);
			trim(calvinDATHeader);

			//Compare remaining header strings, should match once code is complete
			if (gcosDATHeader.compare(calvinDATHeader) != 0)
			{
				//HACK: Don't flag as a failure yet - change when code is complete, uncommenting line below:
				//
				//compareResults = false;
				//
				streamCompare << "Dat Headers don't match" << endl;
				streamCompare << "GCOS DAT Header   : " << gcosDATHeader << endl;
				streamCompare << "Calvin DAT Header : " << calvinDATHeader << endl << endl;
			}
		}
		catch (...)
		{
				streamCompare << "Dat Headers don't match" << endl;
				streamCompare << "GCOS DAT Header   : " << gcosDATHeader << endl;
				streamCompare << "Calvin DAT Header : " << calvinDATHeader << endl << endl;
		}
	}

	int32_t diffIntensities = 0;
	int32_t diffStdev = 0;
	int32_t diffPixels = 0;
	int32_t diffMasks = 0;
	int32_t diffOutliers = 0;

	int count  = gcos.GetNumCells();
	int x = 0, y = 0;
	CELFileEntryType gcosEntry;
	FusionCELFileEntryType fusionEntry;
	for(int i = 0; i < count; ++i)
	{
		x = gcos.IndexToX(i);
		y = gcos.IndexToY(i);

		if(gcos.XYToIndex(x,y) != fusion.XYToIndex(x,y) || gcos.XYToIndex(x,y) != i || fusion.XYToIndex(x,y) != i)
		{
			compareResults = false;

			if (detailedDiffFile)
				streamCompare << "XYToIndex(Index: " << i << "\t" << x << "," << y << "): GCOS = " << gcos.XYToIndex(x,y) << "\tCalvin = " << fusion.XYToIndex(x,y) << endl << endl;
		}

		if( !CompareFloats(gcos.GetStdv(i), fusion.GetStdv(i)) )
		{
			++diffStdev;

			compareResults = false;
			
			if (detailedDiffFile)
				streamCompare << "Stdv by index(" << i << "): GCOS = " << gcos.GetStdv(i) << "\tCalvin = " << fusion.GetStdv(i) << endl << endl;
		}

		if ( !CompareFloats(gcos.GetStdv(x,y), fusion.GetStdv(x,y)) )
		{
			compareResults = false;
			
			if (detailedDiffFile)
				streamCompare << "Stdv by x,y(" << x << "," << y << "): GCOS = " << gcos.GetStdv(x,y) << "\tCalvin = " << fusion.GetStdv(x,y) << endl << endl;
		}

		if(gcos.GetPixels(i) != fusion.GetPixels(i))
		{
			++diffPixels;

			compareResults = false;
			
			if (detailedDiffFile)
				streamCompare << "Pixels by index(" << i << "): GCOS = " << gcos.GetPixels(i) << "\tCalvin = " << fusion.GetPixels(i) << endl << endl;
		}
		if(gcos.GetPixels(x,y) != fusion.GetPixels(x,y))
		{
			compareResults = false;
			
			if (detailedDiffFile)
				streamCompare << "Pixels by x,y(" << x << "," << y << "): GCOS = " << gcos.GetPixels(x,y) << "\tCalvin = " << fusion.GetPixels(x,y) << endl << endl;
		}

		if(!CompareFloats(gcos.GetIntensity(i), fusion.GetIntensity(i)))
		{
			++diffIntensities;

			compareResults = false;
			
			if (detailedDiffFile)
				streamCompare << "Intensity by index(" << i << "): GCOS = " << gcos.GetIntensity(i) << "\tCalvin = " << fusion.GetIntensity(i) << endl << endl;
		}

		if (!CompareFloats(gcos.GetIntensity(x,y), fusion.GetIntensity(x,y)))
		{
			compareResults = false;
			
			if (detailedDiffFile)
				streamCompare << "Intensity by x,y(" << x << "," << y << "): GCOS = " << gcos.GetIntensity(x,y) << "\tCalvin = " << fusion.GetIntensity(x,y) << endl << endl;
		}

		gcos.GetEntry(i,gcosEntry);
		fusion.GetEntry(i,fusionEntry);

		if(!CompareFloats(fusionEntry.Intensity, gcosEntry.Intensity) || fusionEntry.Pixels != gcosEntry.Pixels || !CompareFloats(fusionEntry.Stdv, gcosEntry.Stdv) )
		{
			compareResults = false;
			if (detailedDiffFile)
			{	
				streamCompare << "Entry by index(" << i << ")" << endl;
				streamCompare << "\tIntensity:	GCOS = " << gcosEntry.Intensity << "\tCalvin = " << fusionEntry.Intensity << endl;
				streamCompare << "\tPixels:     GCOS = " << gcosEntry.Pixels << "\tCalvin = " << fusionEntry.Pixels << endl;
				streamCompare << "\tStdv:       GCOS = " << gcosEntry.Stdv << "\tCalvin = " << fusionEntry.Stdv << endl;
			}
		}

		gcos.GetEntry(x,y,gcosEntry);
		fusion.GetEntry(x,y,fusionEntry);

		if(!CompareFloats(fusionEntry.Intensity, gcosEntry.Intensity) || fusionEntry.Pixels != gcosEntry.Pixels || !CompareFloats(fusionEntry.Stdv, gcosEntry.Stdv))
		{
			compareResults = false;
			if (detailedDiffFile)
			{	
				streamCompare << "Entry by x,y(" << x << "," << y << ")" << endl;
				streamCompare << "\tIntensity:	GCOS = " << gcosEntry.Intensity << "\tCalvin = " << fusionEntry.Intensity << endl;
				streamCompare << "\tPixels:     GCOS = " << gcosEntry.Pixels << "\tCalvin = " << fusionEntry.Pixels << endl;
				streamCompare << "\tStdv:       GCOS = " << gcosEntry.Stdv << "\tCalvin = " << fusionEntry.Stdv << endl << endl;
			}
		}

		if (!ignoreMask)
		{
			if(fusion.IsMasked(i) != gcos.IsMasked(i))
			{
				++diffMasks;

				compareResults = false;
				if (detailedDiffFile)
				{
					if (fusion.IsMasked(i))
					{
						streamCompare << "Is Masked by index(" << i << "), Calvin == true" << endl << endl;
					}
					else
					{
						streamCompare << "Is Masked by index(" << i << "), Calvin == false" << endl << endl;
					}
					if (gcos.IsMasked(i))
					{
						streamCompare << "Is Masked by index(" << i << "), GCOS == true" << endl << endl;
					}
					else
					{
						streamCompare << "Is Masked by index(" << i << "), GCOS == false" << endl << endl;
					}
				}
			}

			if(fusion.IsMasked(x,y) != gcos.IsMasked(x,y))
			{
				compareResults = false;

				if (detailedDiffFile)
				{
					if (fusion.IsMasked(x,y))
					{
						streamCompare << "Is Masked by x,y(" << x << "," << y << "), Calvin == true" << endl << endl;
					}
					else
					{
						streamCompare << "Is Masked by x,y(" << x << "," << y << "), Calvin == false" << endl << endl;
					}
					if (gcos.IsMasked(x,y))
					{
						streamCompare << "Is Masked by x,y(" << x << "," << y << "), GCOS == true" << endl << endl;
					}
					else
					{
						streamCompare << "Is Masked by x,y(" << x << "," << y << "), GCOS == false" << endl << endl;
					}
				}
			}
		}

		if(fusion.IsOutlier(i) != gcos.IsOutlier(i))
		{
			++diffOutliers;

			compareResults = false;
			
			if (detailedDiffFile)
			{
				if (fusion.IsOutlier(i))
				{
					streamCompare << "Is Outlier by index(" << i << "), Calvin == true" << endl << endl;
				}
				else
				{
					streamCompare << "Is Outlier by index(" << i << "), Calvin == false" << endl << endl;
				}
				if (gcos.IsOutlier(i))
				{
					streamCompare << "Is Outlier by index(" << i << "), GCOS == true" << endl << endl;
				}
				else
				{
					streamCompare << "Is Outlier by index(" << i << "), GCOS == false" << endl << endl;
				}
			}

		}

		if(fusion.IsOutlier(x,y) != gcos.IsOutlier(x,y))
		{
			compareResults = false;
			
			if (detailedDiffFile)
			{
				if (fusion.IsOutlier(x,y))
				{
					streamCompare << "Is Outlier by x,y(" << x << "," << y << "), Calvin == true" << endl << endl;
				}
				else
				{
					streamCompare << "Is Outlier by x,y(" << x << "," << y << "), Calvin == false" << endl << endl;
				}
				if (gcos.IsOutlier(x,y))
				{
					streamCompare << "Is Outlier by x,y(" << x << "," << y << "), GCOS == true" << endl << endl;
				}
				else
				{
					streamCompare << "Is Outlier by x,y(" << x << "," << y << "), GCOS == false" << endl << endl;
				}
			}
		}

	}

	// Write a summary report
	streamCompare << endl << "Data Comparison Report" << endl << endl;
	streamCompare << "Compared " << count << " cells." << endl;
	streamCompare << "Intensity differences: " << diffIntensities << endl;
	streamCompare << "Stdev differences: " << diffStdev << endl;
	streamCompare << "Num Pixels differences: " << diffPixels << endl;

	if (!ignoreMask)
		streamCompare << "Mask differences: " << diffMasks << endl;

	streamCompare << "Outlier differences: " << diffOutliers << endl;

	streamCompare.close();

	if(compareResults)
	{
#ifdef WIN32
		DeleteFile(outputFile.c_str());
#else
		remove(outputFile.c_str());
#endif
	}
	else
	{
		errorCode = ComparisonFailed;
	}

	return compareResults;
}

/*
 * Helper function, trims excess white space from a std::string value
 */
void CELFileComparer::trim(std::string& str)
{
  string::size_type pos = str.find_last_not_of(' ');
  if(pos != string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if(pos != string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}

/*
 * Helper function, compares float values to within a tolerance defined by a global constant
 */
bool CELFileComparer::CompareFloats(float value1, float value2)
{
#ifdef WIN32
	float result = abs(value1 - value2);
#else
	float result = (value1 > value2 ? value1 - value2 : value2 - value1);
#endif

	return (result <= float_tolerance);

}
