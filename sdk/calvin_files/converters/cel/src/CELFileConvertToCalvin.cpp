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

//
#include "calvin_files/converters/cel/src/CELFileConvertToCalvin.h"
//
#include "calvin_files/converters/cel/src/CELConversionUtilities.h"
#include "calvin_files/converters/utils/src/ConverterFileUtils.h"
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/fusion-dat/src/GCOSParameterNames.h"
#include "calvin_files/parameter/src/CELAlgorithmParameterNames.h"
#include "calvin_files/parsers/src/DATFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCelFileWriter.h"
//
#include "file/CELFileWriter.h"
//
#include <cstring>
#include <iostream>
#include <string.h>
#include <string>
//

using namespace std;
using namespace affxcel;
using namespace affymetrix_calvin_io;
using namespace affymetrix_cel_converter;
using namespace affymetrix_calvin_utilities;

/*
 * Initialize the class.
 */
CELFileConvertToCalvin::CELFileConvertToCalvin()
{
	Clear();
	extraParameters = NULL;
    m_Options = NULL;
}

/*
 * Clear the class
 */
void CELFileConvertToCalvin::Clear()
{
	errorCode = NoConversionError;
}

/*
 * Deallocate any memory.
 */
CELFileConvertToCalvin::~CELFileConvertToCalvin()
{
	Clear();
}

/*
 * Convert the file, the process is the same as the XDA process.
 */
bool CELFileConvertToCalvin::ConvertASCIIFile(const char *fileName, const char *newFile, CELFileConversionOptions *options)
{
	return ConvertXDAFile(fileName, newFile, options, NULL);
}

/*
 * Convert the file.
 */

bool CELFileConvertToCalvin::ConvertXDAFile(const char *fileName, const char *newFile, CELFileConversionOptions *options, const char * datFileName)
{
	Clear();

    // Set Options for use by member methods
    m_Options = options;

	// Read the input file.
	CCELFileData inFile;
	inFile.SetFileName(fileName);
	if (inFile.Read() == false)
	{
		errorCode = UnableToOpenCelFile;
		return false;
	}

	// Create the calvin data object.
	CelFileData data(newFile);

	// Add the parameters.
	CCELFileHeaderData &gcosHeader = inFile.GetHeader();
	int nparams = gcosHeader.GetNumberAlgorithmParameters();
	ParameterNameValueType param;
	for (int iparam=0; iparam<nparams; iparam++)
	{
		string name = gcosHeader.GetAlgorithmParameterTag(iparam);
		param.SetName(StringUtils::ConvertMBSToWCS(name.c_str()));
		param.SetValueAscii(gcosHeader.GetAlgorithmParameter(name.c_str()));
		data.AddAlgorithmParameter(param);
	}

	// Add the grid to the parameters
	GridCoordinatesType grid = inFile.GetGridCorners();
	param.SetName(GRIDULX_PARAM_NAME);
	param.SetValueFloat((float)grid.upperleft.x);
	data.AddAlgorithmParameter(param);

	param.SetName(GRIDULY_PARAM_NAME);
	param.SetValueFloat((float)grid.upperleft.y);
	data.AddAlgorithmParameter(param);

	param.SetName(GRIDURX_PARAM_NAME);
	param.SetValueFloat((float)grid.upperright.x);
	data.AddAlgorithmParameter(param);

	param.SetName(GRIDURY_PARAM_NAME);
	param.SetValueFloat((float)grid.upperright.y);
	data.AddAlgorithmParameter(param);

	param.SetName(GRIDLRX_PARAM_NAME);
	param.SetValueFloat((float)grid.lowerright.x);
	data.AddAlgorithmParameter(param);

	param.SetName(GRIDLRY_PARAM_NAME);
	param.SetValueFloat((float)grid.lowerright.y);
	data.AddAlgorithmParameter(param);

	param.SetName(GRIDLLX_PARAM_NAME);
	param.SetValueFloat((float)grid.lowerleft.x);
	data.AddAlgorithmParameter(param);

	param.SetName(GRIDLLY_PARAM_NAME);
	param.SetValueFloat((float)grid.lowerleft.y);
	data.AddAlgorithmParameter(param);


	// Add the algorithm name, array type and dimensions.
	data.SetAlgorithmName(StringUtils::ConvertMBSToWCS(gcosHeader.GetAlg()));
    if((m_Options != NULL)&&(m_Options->m_ChipType != NULL)) 
    {
	    data.SetArrayType(StringUtils::ConvertMBSToWCS(m_Options->m_ChipType));
    }
    else
    {
	    data.SetArrayType(StringUtils::ConvertMBSToWCS(gcosHeader.GetChipType()));
    }
	data.SetCols(gcosHeader.GetCols());
	data.SetRows(gcosHeader.GetRows());
	data.SetIntensityCount(gcosHeader.GetCells());
	data.SetStdDevCount(gcosHeader.GetCells());
	data.SetPixelCount(gcosHeader.GetCells());
	data.SetOutlierCount(gcosHeader.GetOutliers());
	data.SetMaskCount(gcosHeader.GetMasked());

	// Allow extra parameters to override parameters from the file.
	AddExtraParameters(data);

	// Add the parent header
	// from ConvertASCIIFile
	if ( datFileName == NULL )
	{
		GenericDataHeader hdr;
		if(!FillParentGenericDataHeader(hdr, inFile))
		{
			return false;
		}
		data.GetFileHeader()->GetGenericDataHdr()->AddParent(hdr);

	}
	else // datFileName is not null
	{
		if (strlen(datFileName) == 0)
		{
			GenericDataHeader hdr;
			if(!FillParentGenericDataHeader(hdr, inFile))
			{
				return false;
			}

			// need to change it to multi-scan-acquisition
			hdr.SetFileTypeId("affymetrix-calvin-multi-scan-acquisition");
			data.GetFileHeader()->GetGenericDataHdr()->AddParent(hdr);

		}
		else
		{
			affymetrix_calvin_io::DATData* dat = new DATData;

			// Read DAT header
			DATFileReader reader;
			reader.SetFilename(datFileName);

			try
			{
				reader.Read(*dat);
			}
			catch(affymetrix_calvin_exceptions::CalvinException&)
			{
		//		CloseDAT();
		//		return false;
			}

			// cell file data
			data.GetFileHeader()->GetGenericDataHdr()->AddParent(*dat->GetFileHeader()->GetGenericDataHdr());
		}
	}

	// Create a writer.
	try
	{
		CelFileWriter writer(data);
		int ncells = gcosHeader.GetCells();
		FloatVector intensities;
		intensities.resize(ncells);
		for (int icel=0; icel<ncells; icel++)
		{
			intensities[icel] = inFile.GetIntensity(icel);
		}
		writer.WriteIntensities(intensities);
		intensities.clear();

		FloatVector stdv;
		stdv.resize(ncells);
		for (int icel=0; icel<ncells; icel++)
		{
			stdv[icel] = inFile.GetStdv(icel);
		}
		writer.WriteStdDevs(stdv);
		stdv.clear();

		Int16Vector pix;
		pix.resize(ncells);
		for (int icel=0; icel<ncells; icel++)
		{
			pix[icel] = inFile.GetPixels(icel);
		}
		writer.WritePixels(pix);
		pix.clear();

		int ioutlier=0;
		int noutliers = gcosHeader.GetOutliers();
		XYCoordVector coords(noutliers);
		for (int icel=0; icel<ncells; icel++)
		{
			if (inFile.IsOutlier(icel) == true)
			{
				coords[ioutlier].xCoord = inFile.IndexToX(icel);
				coords[ioutlier].yCoord = inFile.IndexToY(icel);
				++ioutlier;
			}
		}
		writer.WriteOutlierCoords(coords);
		coords.clear();

		int imasked=0;
		int nmasked = gcosHeader.GetMasked();
		coords.resize(nmasked);
		for (int icel=0; icel<ncells; icel++)
		{
			if (inFile.IsMasked(icel) == true)
			{
				coords[imasked].xCoord = inFile.IndexToX(icel);
				coords[imasked].yCoord = inFile.IndexToY(icel);
				++imasked;
			}
		}
		writer.WriteMaskCoords(coords);
	}
	catch (...)
	{
		errorCode = UnableToWriteTheFile;
		return false;
	}
	return true;
}

/*
 * Fills the GenericDataHeader with either the GenericDataHeader from the parent
 * Command Console file or based on the GCOS CEL file.
 */
bool CELFileConvertToCalvin::FillParentGenericDataHeader(GenericDataHeader& hdr, CCELFileData& inFile)
{
	if (GetParentGenericDataHeader(hdr) == false)
	{
		return GenerateParentGenericDataHeaderFromGCOS(hdr, inFile);
	}
    return true;
}

/*
 * Attempts to read the GenericDataHeader from the parent Command Console file.
 */
bool CELFileConvertToCalvin::GetParentGenericDataHeader(GenericDataHeader& hdr)
{
	if (this->parentFileName.length() == 0)
		return false;
	
	// Read DAT header
	DATData dat;
	DATFileReader reader;
	reader.SetFilename(this->parentFileName.c_str());
	try
	{
		reader.Read(dat);
	}
	catch(affymetrix_calvin_exceptions::CalvinException&)
	{
		return false;
	}

	// Copy GenericDataHeader
	hdr = *dat.GetFileHeader()->GetGenericDataHdr();
	return true;
}

/*
 * Generates the parent GenericDataHeader based on the GCOS CEL file.
 */
bool CELFileConvertToCalvin::GenerateParentGenericDataHeaderFromGCOS(GenericDataHeader& hdr, CCELFileData& inFile)
{
	ParameterNameValueType param;
	CCELFileHeaderData &gcosHeader = inFile.GetHeader();

	hdr.SetFileId("");
	hdr.SetFileTypeId(SCAN_ACQUISITION_DATA_TYPE);
	param.SetName(DAT_HEADER_PARAM_NAME);
	wstring datHeader;
    if((m_Options != NULL)&&(m_Options->m_DATFileName != NULL)) 
    {
        string dat = CELFileConversionOptions::newDatName(gcosHeader.GetDatHeader(), m_Options->m_DATFileName);
        if(dat == "") 
        {
            errorCode = UnableToParseDatHeader;
            return false;
        }
        else
        {
            datHeader = StringUtils::ConvertMBSToWCS(dat);
        }
    }
    else
    {
        datHeader = StringUtils::ConvertMBSToWCS(gcosHeader.GetDatHeader());
    }
	param.SetValueText(datHeader);
	hdr.AddNameValParam(param);
	param.SetName(ARRAY_TYPE_PARAM_NAME);
    if((m_Options != NULL)&&(m_Options->m_ChipType != NULL)) 
    {
	    param.SetValueText(StringUtils::ConvertMBSToWCS(m_Options->m_ChipType));
    }
    else
    {
	    param.SetValueText(StringUtils::ConvertMBSToWCS(gcosHeader.GetChipType()));
    }
	hdr.AddNameValParam(param);

	// Extract date from the DAT header string.
	wstring::size_type posColon = datHeader.find(L':');
	if (posColon != wstring::npos)
	{
		// Extract pixel cols and pixel rows from the DAT header string.
		int32_t rows = 0, cols = 0;
		wstring colonString = datHeader.substr(posColon+1);
		if (swscanf(colonString.c_str(), L"CLS=%d RWS=%d XIN=%*d YIN=%*d VE=%*d ", &rows, &cols) == 2)
		{
			param.SetName(ROWS_PARAM_NAME);
			param.SetValueInt32(rows);
			hdr.AddNameValParam(param);
			param.SetName(COLS_PARAM_NAME);
			param.SetValueInt32(cols);
			hdr.AddNameValParam(param);
		}

		wstring::size_type posDate = colonString.find(L"VE=");
		if (posDate != wstring::npos)
		{
			// Move 17 spaces to get to the start of date string from the start of VE=
			// "VE=%-2d "								= 6 spaces
			// "TMP=%-3d" or "       "	= 7 spaces
			// "%-4d" or "%%-4.1f"			= 4 spaces
			posDate += 17;

			wstring dateString = colonString.substr(posDate);
			// Extract scan date and time, scanner id and scanner type from the DAT header string.
			u_int32_t m = 0, d = 0, y = 0, h = 0, min = 0, s = 0;
			if (swscanf(dateString.c_str(), L"%d/%d/%d %d:%d:%d", &m, &d, &y, &h, &min, &s) == 6)
			{
				if (y < 100)
				{
					if (y < 90)	// before affys time
						y += 2000;
					else
						y += 1900;
				}
				param.SetName(SCAN_DATE_PARAM_NAME);
				param.SetValueText(DateTime::FormatDateTime(y, m, d, h, min, s, false));
				hdr.AddNameValParam(param);

				wstring::size_type posEnd = dateString.find(L'\x14');
				if (posEnd != wstring::npos)
				{
					// Limit the string to be considered
					wstring scannerString = dateString.substr(0, posEnd);

					// Jump over date
					wstring::size_type posStart = scannerString.find(L" ");
					if (posStart != wstring::npos)
					{
						scannerString = scannerString.substr(posStart+1);

						// Jump over time
						posStart = scannerString.find(L" ");
						if (posStart != wstring::npos)
						{
							scannerString = scannerString.substr(posStart+1);

							// The first character should be the Scanner Id; it could be an empty string though
							if (scannerString[0] != L' ')
							{
								// Find the end of the scanner ID.
								posStart = scannerString.find(L"  ");
								if (posStart != wstring::npos)
								{
									wstring scannerID = scannerString.substr(0, posStart);
									//wchar_t scannerID[10];
									//swscanf(scannerString.c_str(), L"%10S", scannerID);

									param.SetName(SCANNER_ID_PARAM_NAME);
									param.SetValueText(scannerID);
									hdr.AddNameValParam(param);
								
									// Jump over the scanner id
									scannerString = scannerString.substr(scannerID.length());
								}
							}

							// Jump over 2 spaces between the scanner id and scanner type
							scannerString = scannerString.substr(2);

							// The first character should be the Scanner Type; it could be an empty string though
							if (scannerString[0] != L' ')
							{
								// Find the end of the scanner type.
								posStart = scannerString.find(L"  ");
								if (posStart != wstring::npos)
								{
									wstring scannerType = scannerString.substr(0, posStart);
									//wchar_t scannerType[10];
									//swscanf(scannerString.c_str(), L"%10S", scannerType);

									param.SetName(SCANNER_TYPE_PARAM_NAME);
									param.SetValueText(scannerType);
									hdr.AddNameValParam(param);
								}
							}

							/// TODO: Check that this works
							// If the param name is not scanner type then, it was not added to the hdr.
							if (param.GetName() != SCANNER_TYPE_PARAM_NAME)
							{
								wstring scannerType = GetScannerTypeFromExtraParameters();
								if (scannerType.length() > 0)
								{
									param.SetName(SCANNER_TYPE_PARAM_NAME);
									param.SetValueText(scannerType);
									hdr.AddNameValParam(param);
								}
							}
							/// TODO: End check
						}
					}
				}
			}
		}
	}

	wstring::size_type posChipInfo = datHeader.find(L'\x14');
	if (posChipInfo != wstring::npos)
	{
		// Extract arc radius, laser spot size, pixel size and image orientation from the DAT header string.
		wstring chipInfoString = datHeader.substr(posChipInfo+1);

		wstring arcRadiusString;
		wstring laserSpotSizeString;
		wstring pixelSizeString;
		wstring orientationString;
		for (int32_t i = 0; i <=10; ++i)
		{
			wstring::size_type nextPos = chipInfoString.find(L'\x14');

			wstring itemString;
			if (nextPos != wstring::npos)
				itemString = chipInfoString.substr(0, nextPos);
			else
				itemString = chipInfoString;

			StringUtils::STLTrimLeft(itemString);
			StringUtils::STLTrimRight(itemString);

			if (i == 6 && itemString.length() != 0)	// filter
			{
				int32_t filter = 0;
				if (swscanf(itemString.c_str(), L"%d", &filter) == 1)
				{
					param.SetName(FILTER_PARAM_NAME);
					param.SetValueInt32(filter);
					hdr.AddNameValParam(param);
				}
			}
			else if (i == 7 && itemString.length() != 0)	// arc radius
			{
				float arcRadius = 0.f;
				if (swscanf(itemString.c_str(), L"%f", &arcRadius) == 1)
				{
					if (arcRadius != 0.f)
					{
						param.SetName(ARCRADIUS_PARAM_NAME);
						param.SetValueFloat(arcRadius);
						hdr.AddNameValParam(param);
					}
				}
			}
			else if (i == 8 && itemString.length() != 0)	// laser spot size
			{
				float spotSize = 0.f;
				if (swscanf(itemString.c_str(), L"%f", &spotSize) == 1)
				{
					if (spotSize != 0.f)
					{
						param.SetName(LASER_SPOTSIZE_PARAM_NAME);
						param.SetValueFloat(spotSize);
						hdr.AddNameValParam(param);
					}
				}
			}
			else if (i == 9 && itemString.length() != 0)	// pixel size
			{
				float pixelSize = 0.f;
				if (swscanf(itemString.c_str(), L"%f", &pixelSize) == 1)
				{
					if (pixelSize != 0.f)
					{
						param.SetName(PIXEL_SIZE_PARAM_NAME);
						param.SetValueFloat(pixelSize);
						hdr.AddNameValParam(param);
					}
				}
			}
			else if (i == 10 && itemString.length() != 0)	// orientation
			{
				u_int32_t orientation = 0;
				if (swscanf(itemString.c_str(), L"%d", &orientation) == 1)
				{
					param.SetName(ORIENTATION_PARAM_NAME);
					param.SetValueUInt8((u_int8_t)orientation);
					hdr.AddNameValParam(param);
				}
			}

			if (nextPos != wstring::npos)
				chipInfoString = chipInfoString.substr(nextPos+1);
		}
	}

	if (arrayID.empty() == false || arrayBarcode.empty() == false)
		// Add ArraySet GenericDataHeader
		ConverterCreateAndAddParentArrayGenericDataHeader(hdr, arrayID, arrayBarcode);

    return true;
}

/*
 * Adds any extra parameters from the extraParameter collection to the file header.
 */
void CELFileConvertToCalvin::AddExtraParameters(CelFileData& data)
{
	if (this->extraParameters == NULL)
		return;

	GenericDataHeader* gdh = data.GetFileHeader()->GetGenericDataHdr();

	// Copy the parameters to the file.
	for (ParameterNameValueTypeVector::const_iterator ii = extraParameters->begin(); ii != extraParameters->end(); ++ii)
	{
		if (ii->GetName() == SCANNER_TYPE_PARAM_NAME)
		{
			continue;
		}

		gdh->AddNameValParam(*ii);
	}
}

/*
 * TODO Check that this works.
 */
std::wstring CELFileConvertToCalvin::GetScannerTypeFromExtraParameters()
{
	if (this->extraParameters == NULL)
		return L"";

	for (ParameterNameValueTypeVector::const_iterator ii = extraParameters->begin(); ii != extraParameters->end(); ++ii)
	{
		if (ii->GetName() == SCANNER_TYPE_PARAM_NAME)
		{
			return ii->GetValueText();
		}
	}
	return L"";
}

