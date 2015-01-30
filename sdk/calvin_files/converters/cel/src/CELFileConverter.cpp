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


#include "calvin_files/converters/cel/src/CELFileConverter.h"
//
#include "calvin_files/converters/cel/src/CELFileConvertToASCII.h"
#include "calvin_files/converters/cel/src/CELFileConvertToCalvin.h"
#include "calvin_files/converters/cel/src/CELFileConvertToXDA.h"
#include "calvin_files/converters/cel/src/CELFileVersion.h"
#include "calvin_files/converters/utils/src/ConverterFileUtils.h"
//
#include <cstring>
#include <iostream>
#include <string>
//

using namespace std;
using namespace affymetrix_cel_converter;

/*
 * Initialize the class.
 */
CELFileConverter::CELFileConverter()
{
	Clear();
}

/*
 * Clears the members.
 */
void CELFileConverter::Clear()
{
	errorCode = NoConversionError;
}

/*
 * Deallocate any memory.
 */
CELFileConverter::~CELFileConverter()
{
}

/*
 * Basic checks before starting.
 */
bool CELFileConverter::Checks(const char *fileName, CELFileVersionType fileVersion) 
{
	Clear();

	// Check if the file exists
	if (ConverterFileExists(fileName) == false)
	{
		errorCode = FileDoesNotExist;
		return false;
	}

	// Determine the format
	CELFileVersionType inVersion = CELFileVersion::DetermineCELFileVersion(fileName);

	// If unknown the return a false.
	if (inVersion == Unknown_Version)
	{
		errorCode = InvalidCelFileFormat;
		return false;
	}

    return true;

}

/*
 * Convert the file (inplace).
 */
bool CELFileConverter::ConvertFile(const char *fileName, CELFileVersionType fileVersion, CELFileConversionOptions *options)
{
    // Basic checks
    if(!Checks(fileName,fileVersion))
        return false;

	// Determine the format
	CELFileVersionType inVersion = CELFileVersion::DetermineCELFileVersion(fileName);

	// If same then done.
	if (inVersion == fileVersion)
	{
        return true;
    }
       
	// Move the input file name to a backup.
	string bakSuffix = ".bak";
	string bakFile = fileName + bakSuffix;
	if (ConverterMoveFile(fileName, bakFile.c_str(), true) == false)
	{
		errorCode = UnableToRenameInputFile;
		return false;
	}


	// Check the status and revert if needed.
	if (ConvertFile(bakFile.c_str(),fileName, fileVersion))
	{
		ConverterRemoveFile(bakFile.c_str());
	    return true;
	}
	else
	{
		ConverterMoveFile(bakFile.c_str(), fileName, true);
	    return false;
	}
}

/*
 * Convert the file.
 */
bool CELFileConverter::ConvertFile(const char *fileName, const char *newFileName, CELFileVersionType fileVersion, CELFileConversionOptions *options)
{

    // Basic checks
    if(!Checks(fileName,fileVersion))
        return false;

	// Determine the format
	CELFileVersionType inVersion = CELFileVersion::DetermineCELFileVersion(fileName);

	// If same then just copy.
	if (inVersion == fileVersion)
	{
        if(options != NULL) 
        {
            errorCode = UnableToMixConversionOptionsAndNoFormatChange;
            return false;
        }
        else
        {
		    if(ConverterCopyFile(fileName, newFileName, true))
            {
                return true;
            }
            else
            {
                errorCode = UnableToCopyFile;
                return false;
            }
        }
	}

	// Convert the file.
	bool status = false;
	if (fileVersion == GCOS_Version3)
	{
		CELFileConvertToASCII outFile;
		if (inVersion == Calvin_Version1)
		{
			status = outFile.ConvertCalvinFile(fileName, newFileName, options);
		}
		else if (inVersion == GCOS_Version4)
		{
			status = outFile.ConvertXDAFile(fileName, newFileName, options);
		}
		errorCode = outFile.ErrorCode();
	}
	else if (fileVersion == GCOS_Version4)
	{
		CELFileConvertToXDA outFile;
		if (inVersion == Calvin_Version1)
		{
			status = outFile.ConvertCalvinFile(fileName, newFileName, options);
		}
		else if (inVersion == GCOS_Version3)
		{
			status = outFile.ConvertASCIIFile(fileName, newFileName, options);
		}
		errorCode = outFile.ErrorCode();
	}
	else if (fileVersion == Calvin_Version1)
	{
		CELFileConvertToCalvin outFile;
		if (inVersion == GCOS_Version3)
		{
			status = outFile.ConvertASCIIFile(fileName, newFileName, options);
		}
		else if (inVersion == GCOS_Version4)
		{
			status = outFile.ConvertXDAFile(fileName, newFileName, options);
		}
		errorCode = outFile.ErrorCode();
	}
	else
	{
		errorCode = InvalidConversionInputs;
		status = false;
	}

	// Check the status and revert if needed.
	if (status == false)
	{
		ConverterRemoveFile(newFileName);
	}

	return status;
}


/*
bool CELFileConverter::CreateFile()
{
	// Get the number of probe sets.
	CCELFileHeaderData &gcosHeader = gcosCel.GetHeader();

	// Create a Calvin CEL object.
	outputFile = fileName.substr(0, fileName.length()-4) + ".calvin.CEL";
	CelFileData data(outputFile);

	// Add the header information.
	int nparams = gcosHeader.GetNumberAlgorithmParameters();
	for (int iparam=0; iparam<nparams; iparam++)
	{
		ParameterNameValueType param;
		string name = gcosHeader.GetAlgorithmParameterTag(iparam);
		param.SetName(StringUtils::ConvertMBSToWCS(name.c_str()));
		param.SetValueText(StringUtils::ConvertMBSToWCS(gcosHeader.GetAlgorithmParameter(name.c_str())));
		data.AddAlgorithmParameter(param);
	}
	data.SetAlgorithmName(StringUtils::ConvertMBSToWCS(gcosHeader.GetAlg()));
	data.SetArrayType(StringUtils::ConvertMBSToWCS(gcosHeader.GetChipType()));
	data.SetCols(gcosHeader.GetCols());
	data.SetRows(gcosHeader.GetRows());
	data.SetIntensityCount(gcosHeader.GetCells());
	data.SetStdDevCount(gcosHeader.GetCells());
	data.SetPixelCount(gcosHeader.GetCells());
	data.SetOutlierCount(gcosHeader.GetOutliers());
	data.SetMaskCount(gcosHeader.GetMasked());

	// Create a writer.
	try
	{
		CelFileWriter writer(data);
		int ncells = gcosHeader.GetCells();
		FloatVector intensities;
		intensities.resize(ncells);
		for (int icel=0; icel<ncells; icel++)
		{
			intensities[icel] = gcosCel.GetIntensity(icel);
		}
		writer.WriteIntensities(intensities);
		intensities.clear();

		FloatVector stdv;
		stdv.resize(ncells);
		for (int icel=0; icel<ncells; icel++)
		{
			stdv[icel] = gcosCel.GetStdv(icel);
		}
		writer.WriteStdDevs(stdv);
		stdv.clear();

		Int16Vector pix;
		pix.resize(ncells);
		for (int icel=0; icel<ncells; icel++)
		{
			pix[icel] = gcosCel.GetPixels(icel);
		}
		writer.WritePixels(pix);
		pix.clear();

		int ioutlier=0;
		int noutliers = gcosHeader.GetOutliers();
		XYCoordVector coords(noutliers);
		for (int icel=0; icel<ncells; icel++)
		{
			if (gcosCel.IsOutlier(icel) == true)
			{
				coords[ioutlier].xCoord = gcosCel.IndexToX(icel);
				coords[ioutlier].yCoord = gcosCel.IndexToY(icel);
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
			if (gcosCel.IsMasked(icel) == true)
			{
				coords[imasked].xCoord = gcosCel.IndexToX(icel);
				coords[imasked].yCoord = gcosCel.IndexToY(icel);
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
*/

