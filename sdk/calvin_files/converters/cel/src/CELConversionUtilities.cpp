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


#include "calvin_files/converters/cel/src/CELConversionUtilities.h"
//
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/parameter/src/CELAlgorithmParameterNames.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstring>
#include <string.h>
#include <string>
//

using namespace std;
using namespace affxcel;
using namespace affymetrix_calvin_io;
using namespace affymetrix_cel_converter;
using namespace affymetrix_calvin_utilities;

#ifdef WIN32
#pragma warning(disable:4996) // Don't show deprecated messages.
#endif

/*
 * Get the grid from the parameters.
 */
GridCoordinatesType CELConversionUtilities::GetGrid(CelFileData &inFile)
{
	GridCoordinatesType grid;
	ParameterNameValueType nvt;

	grid.lowerleft.x = grid.lowerleft.y = 0;
	grid.upperleft.x = grid.upperleft.y = 0;
	grid.lowerright.x = grid.lowerright.y = 0;
	grid.upperright.x = grid.upperright.y = 0;

	if (inFile.FindAlgorithmParameter(GRIDULX_PARAM_NAME, nvt) && nvt.GetParameterType() == ParameterNameValueType::FloatType)
		grid.upperleft.x = (int) nvt.GetValueFloat();

	if (inFile.FindAlgorithmParameter(GRIDULY_PARAM_NAME, nvt) && nvt.GetParameterType() == ParameterNameValueType::FloatType)
		grid.upperleft.y = (int) nvt.GetValueFloat();

	if (inFile.FindAlgorithmParameter(GRIDURX_PARAM_NAME, nvt) && nvt.GetParameterType() == ParameterNameValueType::FloatType)
		grid.upperright.x = (int) nvt.GetValueFloat();

	if (inFile.FindAlgorithmParameter(GRIDURY_PARAM_NAME, nvt) && nvt.GetParameterType() == ParameterNameValueType::FloatType)
		grid.upperright.y = (int) nvt.GetValueFloat();

	if (inFile.FindAlgorithmParameter(GRIDLRX_PARAM_NAME, nvt) && nvt.GetParameterType() == ParameterNameValueType::FloatType)
		grid.lowerright.x = (int) nvt.GetValueFloat();

	if (inFile.FindAlgorithmParameter(GRIDLRY_PARAM_NAME, nvt) && nvt.GetParameterType() == ParameterNameValueType::FloatType)
		grid.lowerright.y = (int) nvt.GetValueFloat();

	if (inFile.FindAlgorithmParameter(GRIDLLX_PARAM_NAME, nvt) && nvt.GetParameterType() == ParameterNameValueType::FloatType)
		grid.lowerleft.x = (int) nvt.GetValueFloat();

	if (inFile.FindAlgorithmParameter(GRIDLLY_PARAM_NAME, nvt) && nvt.GetParameterType() == ParameterNameValueType::FloatType)
		grid.lowerleft.y = (int) nvt.GetValueFloat();

	return grid;
}

/*
 * Get the DAT header from the DAT file header object.
 */
string CELConversionUtilities::GetDatHeader(CelFileData &inFile)
{
	std::wstring datHeader;

	GenDataHdrVectorIt begin, end; 
	inFile.GetFileHeader()->GetGenericDataHdr()->GetParentIterators(begin, end);

	// Find the DAT generic header
	for (GenDataHdrVectorIt ii = begin; ii != end; ++ii)
	{
		std::string s = ii->GetFileTypeId();
		if (ii->GetFileTypeId() == SCAN_ACQUISITION_DATA_TYPE)
		{

			// found the right header, now look for the parameter
			ParameterNameValueType nvt;
			if (ii->FindNameValParam(DAT_HEADER_PARAM_NAME, nvt))
			{
				if (nvt.GetParameterType() == ParameterNameValueType::TextType)
					datHeader = nvt.GetValueText();
			}
			else if (ii->FindNameValParam(PARTIAL_DAT_HEADER_PARAM_NAME, nvt))
			{
				if (nvt.GetParameterType() == ParameterNameValueType::TextType)
				{
					std::wstring partialDatHeader = nvt.GetValueText();

					u_int16_t min = 0;
					u_int16_t max = 0;

					// Find the max and min parameters and append to the string.
					if (ii->FindNameValParam(MAX_PIXEL_INTENSITY_PARAM_NAME, nvt))
					{
						if (nvt.GetParameterType() == ParameterNameValueType::UInt16Type)
							max = nvt.GetValueUInt16();
					}

					if (ii->FindNameValParam(MIN_PIXEL_INTENSITY_PARAM_NAME, nvt))
					{
						if (nvt.GetParameterType() == ParameterNameValueType::UInt16Type)
							min = nvt.GetValueUInt16();
					}

					wchar_t buf[30];
					FormatString2(buf, 30, L"[%d..%d]", min, max);
					datHeader = buf;
					datHeader += partialDatHeader;
				}
			}
			break;
		}
	}
	return StringUtils::ConvertWCSToMBS(datHeader);
}

/*
 * Copy the file contents.
 */
bool CELConversionUtilities::Copy(CelFileData &inFile, CCELFileData &outFile, CELFileConversionOptions *options)
{
	// Set the dimensions
	outFile.SetDimensions(inFile.GetRows(), inFile.GetCols());

	// Set the intensities, stdv and pixels
	FloatVector values;
	int n = inFile.GetNumCells();
	inFile.GetIntensities(0, n, values);
	for (int i=0; i<n; i++)
	{
		outFile.SetIntensity(i, values[i]);
	}
	values.clear();
	inFile.GetStdev(0, n, values);
	for (int i=0; i<n; i++)
	{
		outFile.SetStdv(i, values[i]);
	}
	Int16Vector pixels;
	inFile.GetNumPixels(0, n, pixels);
	for (int i=0; i<n; i++)
	{
		outFile.SetPixels(i, pixels[i]);
	}

	// Set the outliers and masks.
	XYCoordVector xy;
	inFile.GetOutlierCoords(xy);
	n = (int) xy.size();
	for (int i=0; i<n; i++)
	{
		outFile.SetOutlier(xy[i].xCoord, xy[i].yCoord, true);
	}
	xy.clear();
	inFile.GetMaskedCoords(xy);
	n = (int) xy.size();
	for (int i=0; i<n; i++)
	{
		outFile.SetMask(xy[i].xCoord, xy[i].yCoord, true);
	}


	// Set the header.
	outFile.SetGridCorners(GetGrid(inFile));
    if((options != NULL)&&(options->m_ChipType != NULL)) 
    {
	    outFile.SetChipType(options->m_ChipType);
    } 
    else 
    {
	    outFile.SetChipType(StringUtils::ConvertWCSToMBS(inFile.GetArrayType()).c_str());
    }

	// Make sure the algorithm name is "Percentile" so it can import into GCOS.
	outFile.SetAlgorithmName("Percentile");
	//outFile.SetAlgorithmName(StringUtils::ConvertWCSToMBS(inFile.GetAlgorithmName()).c_str());

	// Add the first four algorithm parameters required by GCOS.
	FindAndAddAlgorithmParameter("Percentile", "75", inFile, outFile);
	FindAndAddAlgorithmParameter("CellMargin", "2", inFile, outFile);
	FindAndAddAlgorithmParameter("OutlierHigh", "1.500", inFile, outFile);
	FindAndAddAlgorithmParameter("OutlierLow", "1.004", inFile, outFile);

	// Set the CellMargin
	SetCellMargin(inFile, outFile, 2);

	ParameterNameValueTypeVector algParams;
	inFile.GetAlgorithmParameters(algParams);
	n = (int) algParams.size();
	string tag;
	string value;
	for (int i=0; i<n; i++)
	{
		tag = StringUtils::ConvertWCSToMBS(algParams[i].GetName());
		value = StringUtils::ConvertWCSToMBS(algParams[i].ToString());

		// Filter the first four required GCOS algorithm parameters.
		if (tag == "Percentile" || tag == "CellMargin" || tag == "OutlierHigh" || tag == "OutlierLow")
		{
			continue;
		}

		// Don't include the grid parameters added by the converter.
		if (tag.length() > strlen("Grid") && strncmp(tag.c_str(), "Grid", strlen("Grid")) == 0)
			continue;

		// Add the parameter;
		outFile.AddAlgorithmParameter(tag.c_str(), value.c_str());
	}
    if((options != NULL)&&(options->m_DATFileName != NULL)) 
    {
        std::string datHeader = CELFileConversionOptions::newDatName(GetDatHeader(inFile).c_str(), options->m_DATFileName);
        if(datHeader == "") 
        {
            return false;
        }
        else
        {
	        outFile.GetHeader().SetDatHeader(datHeader.c_str());
        }
    }
    else
    {
	    outFile.GetHeader().SetDatHeader(GetDatHeader(inFile).c_str());
    }
	outFile.GetHeader().SetHeader(outFile.GetHeader().GetHeader().c_str());

    return true;
}

/*
 * Copy the file contents.
 */
bool CELConversionUtilities::Copy(CCELFileData &inFile, CCELFileData &outFile, CELFileConversionOptions *options)
{
	// Set and dimensions
	outFile.SetDimensions(inFile.GetRows(), inFile.GetCols());

	// Set the entries
	int n = inFile.GetNumCells();
	for (int i=0; i<n; i++)
	{
        // Using GetEntry and SetCellEntry
        // causes problems on OS-X PPC
		outFile.SetIntensity(i,inFile.GetIntensity(i));
		outFile.SetStdv(i,inFile.GetStdv(i));
		outFile.SetPixels(i,inFile.GetPixels(i));
		if (inFile.IsOutlier(i) == true)
			outFile.SetOutlier(i, true);
		if (inFile.IsMasked(i) == true)
			outFile.SetMask(i, true);
	}

	// Set the header.
	outFile.SetGridCorners(inFile.GetGridCorners());
    if((options != NULL)&&(options->m_ChipType != NULL)) 
    {
	    outFile.SetChipType(options->m_ChipType);
    } 
    else
    {
	    outFile.SetChipType(inFile.GetChipType().c_str());
    }
	outFile.SetAlgorithmName(inFile.GetAlg().c_str());
	n = inFile.GetNumberAlgorithmParameters();
	string tag;
	string value;
	for (int i=0; i<n; i++)
	{
		tag = inFile.GetAlgorithmParameterTag(i);
		value = inFile.GetAlgorithmParameter(tag.c_str());
		outFile.AddAlgorithmParameter(tag.c_str(), value.c_str());
	}
    if((options != NULL)&&(options->m_DATFileName != NULL)) 
    {
        std::string datHeader = CELFileConversionOptions::newDatName(inFile.GetHeader().GetDatHeader().c_str(), options->m_DATFileName).c_str();
        if(datHeader == "") 
        {
            return false;
        }
        else
        {
	        outFile.GetHeader().SetDatHeader(datHeader.c_str());
        }
    }
    else
    {
	    outFile.GetHeader().SetDatHeader(inFile.GetHeader().GetDatHeader().c_str());
    }
	outFile.GetHeader().SetHeader(inFile.GetHeader().GetHeader().c_str());

    return true;
}

void CELConversionUtilities::SetCellMargin(CelFileData &inFile, CCELFileData &outFile, int defaultValue)
{
	ParameterNameValueType nvt;
	if (inFile.FindAlgorithmParameter(L"CellMargin", nvt))
	{
		if (nvt.GetParameterType() == ParameterNameValueType::Int32Type)
		{
			outFile.GetHeader().SetMargin(nvt.GetValueInt32());
		}
		else if(nvt.GetParameterType() == ParameterNameValueType::AsciiType)
		{
			outFile.GetHeader().SetMargin(atoi(nvt.GetValueAscii().c_str()));
		}
	}
	else
	{
		outFile.GetHeader().SetMargin(defaultValue);
	}
}

void CELConversionUtilities::FindAndAddAlgorithmParameter(string name, string defaultValue, CelFileData &inFile, CCELFileData &outFile)
{
	ParameterNameValueType nvt;
	if (inFile.FindAlgorithmParameter(StringUtils::ConvertMBSToWCS(name), nvt))
	{
		outFile.AddAlgorithmParameter(name.c_str(), StringUtils::ConvertWCSToMBS(nvt.ToString()).c_str());
	}
	else
	{
		outFile.AddAlgorithmParameter(name.c_str(), defaultValue.c_str());
	}
}
