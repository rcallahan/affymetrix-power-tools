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
#include "calvin_files/converters/chp/src/CHPFileConvertToCalvin.h"
//
#include "calvin_files/converters/utils/src/ConverterFileUtils.h"
#include "calvin_files/data/src/CHPData.h"
#include "calvin_files/parsers/src/CelFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPFileWriter.h"
//
#include "file/CDFFileData.h"
#include "file/CHPFileData.h"
#include "file/PSIFileData.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cstring>
#include <iostream>
#include <string.h>
#include <string>
#include <stdio.h>
#include <vector>
//

using namespace std;
using namespace affxchp;
using namespace affxpsi;
using namespace affxcdf;
using namespace affymetrix_chp_converter;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;

/* Constants
 */
const int32_t HZ = 4;
const int32_t VZ = 4;
const float BG = 2.0f;

/*
 * Initialize the class.
 */
CHPFileConvertToCalvin::CHPFileConvertToCalvin()
{
	Clear();
	extraParameters = NULL;
}

/*
 * Clear the class
 */
void CHPFileConvertToCalvin::Clear()
{
	errorCode = NoConversionError;
}

/*
 * Deallocate any memory.
 */
CHPFileConvertToCalvin::~CHPFileConvertToCalvin()
{
	Clear();
}

/*
 * Get the assay type based on the assay type in the CHP file.
 */
static bool GetAssayType(string &assayType, CCHPFileData &inFile)
{
	switch (inFile.GetHeader().GetAssayType())
	{
	case CCHPFileHeader::Expression:
		assayType = CHP_EXPRESSION_ASSAY_TYPE;
		break;
	case CCHPFileHeader::Genotyping:
		assayType = CHP_GENOTYPING_ASSAY_TYPE;
		break;
	case CCHPFileHeader::Resequencing:
		assayType = CHP_RESEQUENCING_ASSAY_TYPE;
		break;
	case CCHPFileHeader::Universal:
		assayType = CHP_UNIVERSAL_ASSAY_TYPE;
		break;
	default:
		return false;
		break;
	}
	return true;
}


static bool CheckAssayType(CCHPFileData &inFile)
{
	if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping)
	{
		return false;
	}
	return true;
}

/*
 * Read the probe set names from either the PSI or CDF file.
 * Note: resequencing arrays do not have probe set names.
 */
static bool ReadProbeSetNames(CCHPFileData &inFile, const char *libPath, vector<string> &psNames)
{
	if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Resequencing)
		return true;

	if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression ||
		inFile.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping)
	{
		string libraryPath = libPath;
		string psiFile = libPath + inFile.GetHeader().GetChipType() + ".PSI";
		CPSIFileData psi;
		psi.SetFileName(psiFile.c_str());
		if (psi.Read() == true)
		{
			int n=psi.GetProbeSetCount();
			psNames.resize(n);
			for (int i=0; i<n; i++)
			{
				psNames[i] = psi.GetProbeSetName(i);
			}
		}
	}

	// Try reading the CDF file if PSI failed and not resequencing.
	if (psNames.size() == 0)
	{
		string libraryPath = libPath;
		string cdfFile = libPath + inFile.GetHeader().GetChipType() + ".CDF";
		CCDFFileData cdf;
		cdf.SetFileName(cdfFile.c_str());
		if (cdf.Read() == true)
		{
			int n = cdf.GetHeader().GetNumProbeSets();
			psNames.resize(n);
			for (int i=0; i<n; i++)
			{
				psNames[i] = cdf.GetProbeSetName(i);
			}
		}
	}

	return (psNames.size() > 0);
}

/*
 * Write the expression data.
 */
static void WriteExpressionData(CCHPFileData &inFile, vector<string> &psNames, CHPFileWriter &writer)
{
	writer.SeekToDataSet();
	int nps = inFile.GetHeader().GetNumProbeSets();
	for (int i=0; i<nps; i++)
	{
		CExpressionProbeSetResults *results = inFile.GetExpressionResults(i);
		CHPExpressionEntry e(
			psNames[i],
			results->Detection,
			results->DetectionPValue,
			results->Signal,
			results->NumPairs,
			results->NumUsedPairs,
			results->m_HasCompResults,
			results->Change,
			results->ChangePValue,
			results->SignalLogRatio,
			results->SignalLogRatioLow,
			results->SignalLogRatioHigh,
			results->NumCommonPairs);
		writer.WriteExpressionEntry(e);
	}
}

/*
 * Write the backgroun data
 */
static void WriteBackground(CCHPFileData &inFile, CHPFileWriter &writer)
{
	writer.SeekToBgSet();
	CCHPFileHeader &gcosHeader = inFile.GetHeader();
	for (BackgroundZoneTypeList::iterator it=gcosHeader.GetBackgroundZones().begin(); it!=gcosHeader.GetBackgroundZones().end(); ++it)
	{
		CHPBackgroundZone bz(it->centerx, it->centery, it->background, gcosHeader.GetBackgroundZoneInfo().smooth_factor);
		writer.WriteBackgroundZone(bz);
	}
}

static void WriteBackground(const AllZonesInfoType& zonesInfo, CHPFileWriter &writer)
{
	writer.SeekToBgSet();
	int cnt = zonesInfo.number_zones;
	for (int i = 0; i < cnt; i++)
	{
		CHPBackgroundZone bz(zonesInfo.pZones[i].center.x, 
			zonesInfo.pZones[i].center.y, 
			zonesInfo.pZones[i].background, 
			zonesInfo.smooth_factor);
		writer.WriteBackgroundZone(bz);
	}
}

/*
 * Write the genotyping data.
 */
static void WriteGenotypingData(CCHPFileData &inFile, vector<string> &psNames, CHPFileWriter &writer)
{
	writer.SeekToDataSet();
	int nps = inFile.GetHeader().GetNumProbeSets();
	for (int i=0; i<nps; i++)
	{
		CGenotypeProbeSetResults *results = inFile.GetGenotypingResults(i);
		CHPGenotypeEntry e(
			psNames[i],
			results->AlleleCall,
			results->Confidence,
			results->RAS1,
			results->RAS2,
			results->pvalue_AA,
			results->pvalue_AB,
			results->pvalue_BB,
			results->pvalue_NoCall);
		writer.WriteGenotypeEntry(e);
	}
}

/*
 * Write the universal data.
 */
static void WriteUniversalData(CCHPFileData &inFile, CHPFileWriter &writer)
{
	writer.SeekToDataSet();
	int nps = inFile.GetHeader().GetNumProbeSets();
	for (int i=0; i<nps; i++)
	{
		CUniversalProbeSetResults *results = inFile.GetUniversalResults(i);
		CHPUniversalEntry e(
			results->GetBackground());
		writer.WriteUniversalEntry(e);
	}
}

/*
 * Write the resequencing data.
 */
static void WriteResequencingData(CCHPFileData &inFile, CHPFileWriter &writer)
{
	writer.SeekToDataSet();
	CHPReseqEntry entry;
	CResequencingResults *results = inFile.GetResequencingResults();
	int n = results->GetCalledBasesSize();
	for (int i=0; i<n; i++)
	{
		entry.call = results->GetCalledBase(i);
		entry.score = results->GetScore(i);
		writer.WriteReseqEntry(entry);
	}

	writer.SeekToForceSet();
	n = results->GetForceCallsSize();
	CHPReseqForceCall force;
	for (int i=0; i<n; i++)
	{
		force.position = results->GetForceCall(i).position;
		force.call = results->GetForceCall(i).call;
		force.reason = results->GetForceCall(i).reason;
		writer.WriteForceCall(force);
	}

	writer.SeekToOrigCallSet();
	n = results->GetOrigCallsSize();
	CHPReseqOrigCall orig;
	for (int i=0; i<n; i++)
	{
		orig.position = results->GetOrigCall(i).position;
		orig.call = results->GetOrigCall(i).call;
		writer.WriteOrigCall(orig);
	}
}

/*
 * Adds zone algorithm parameters to the insertedParams list
 */
static void AddZoneAlgParameters(TagValuePairTypeList& insertedParams)
{
	char buf[10];
	TagValuePairType param;
	param.Tag = "HZ";
	sprintf(buf, "%d", HZ);
	param.Value = buf;
	insertedParams.push_back(param);
	param.Tag = "VZ";
	sprintf(buf, "%d", VZ);
	param.Value = buf;
	insertedParams.push_back(param);
	param.Tag = "BG";
	sprintf(buf, "%d", (int32_t)BG);
	param.Value = buf;
	insertedParams.push_back(param);
}

/*
 * Set the header of the CHP object.
 */
static void SetHeader(CCHPFileData &inFile, CHPData &data, TagValuePairTypeList* insertedAlgParams)
{
	CCHPFileHeader &gcosHeader = inFile.GetHeader();
	data.SetAlgName(StringUtils::ConvertMBSToWCS(gcosHeader.GetAlgName()));
	data.SetAlgVersion(StringUtils::ConvertMBSToWCS(gcosHeader.GetAlgVersion()));
	data.SetArrayType(StringUtils::ConvertMBSToWCS(gcosHeader.GetChipType()));
	data.SetCols(gcosHeader.GetCols());
	data.SetRows(gcosHeader.GetRows());
	data.SetParentCell(StringUtils::ConvertMBSToWCS(gcosHeader.GetParentCellFile()));
	data.SetProgId(StringUtils::ConvertMBSToWCS(gcosHeader.GetProgID()));
	if (insertedAlgParams != NULL)
	{
		for (TagValuePairTypeList::iterator it = insertedAlgParams->begin(); it != insertedAlgParams->end(); ++it)
		{
			data.AddAlgParam(StringUtils::ConvertMBSToWCS(it->Tag), StringUtils::ConvertMBSToWCS(it->Value));
		}
	}
	TagValuePairTypeList &params = gcosHeader.AlgorithmParameters();
	for (TagValuePairTypeList::iterator it = params.begin(); it != params.end(); ++it)
	{
		data.AddAlgParam(StringUtils::ConvertMBSToWCS(it->Tag), StringUtils::ConvertMBSToWCS(it->Value));
	}

	TagValuePairTypeList &summary = gcosHeader.SummaryParameters();
	for (TagValuePairTypeList::iterator it = summary.begin(); it != summary.end(); ++it)
	{
		data.AddChipSum(StringUtils::ConvertMBSToWCS(it->Tag), StringUtils::ConvertMBSToWCS(it->Value));
	}
}

bool CHPFileConvertToCalvin::ComputeBackgroundZoneInfo(FusionCELData& celData, 
																											 const char* libPath,
																											 AllZonesInfoType& zonesInfo) const
{
	return ComputeBackgroundZoneInfo(celData, libPath, NULL, zonesInfo);
}

bool CHPFileConvertToCalvin::ComputeBackgroundZoneInfo(FusionCELData& celData, 
																											 const char* libPath, 
																											 const char* maskPath,
																											 AllZonesInfoType& zonesInfo) const
{
	CExpressionAlgorithmImplementation algorithm;

	CExpStatAlgSettings settings = algorithm.GetParameters();
	settings.Epsilon = 0.5;
	if(maskPath != NULL)
	{
		settings.ScaleMaskFile = string(maskPath);
	}
	settings.NumberHorZones = HZ;
	settings.NumberVertZones = VZ;
	settings.NumberBGCells = BG;

	algorithm.SetLibPath(libPath);
	if(algorithm.ReadCdfFile(celData))
	{
		zonesInfo.number_zones = 16;
		zonesInfo.smooth_factor = 100;
		vector<float> featureIntensity;
		algorithm.ComputeBackGroundZones(&celData, zonesInfo, featureIntensity);
		return true;
	}
	else
	{
		return false;
	}
}

/*
 * Convert the file.
 */
bool CHPFileConvertToCalvin::ConvertXDAFile(const char *fileName, const char *libPath, const char *newFile)
{
	Clear();

	// Read the input file.
	CCHPFileData inFile;
	inFile.SetFileName(fileName);
	if (inFile.Read() == false)
	{
		errorCode = UnableToOpenChpFile;
		return false;
	}

	// Determine the assay type.
	string assayType;
	if (GetAssayType(assayType, inFile) == false)
	{
		errorCode = InvalidChpFileFormat;
		return false;
	}

	// Check that we can convert the assay type to Calvin
	if (CheckAssayType(inFile) == false)
	{
		errorCode = InvalidAssayType;
		return false;
	}

	// Check the algorithm type
	if (CheckAlgorithmType(inFile.GetHeader().GetAlgName()) == false)
	{
		errorCode = InvalidAlgorithmType;
		return false;
	}

	// Create the new CHP object.
	CHPData data(newFile, assayType);

	// Set the data
	SetHeader(inFile, data, NULL);

	// Set the parent name
	if (parentFileName.empty() == false)
	{
		// Attempt to remove the path.
		string::size_type posBS = parentFileName.rfind('\\');
		string::size_type posS = parentFileName.rfind('/');
		if (posBS != -1)
		{
			data.SetParentCell(StringUtils::ConvertMBSToWCS(parentFileName.substr(posBS+1)));
		}
		else if (posS != -1)
		{
			data.SetParentCell(StringUtils::ConvertMBSToWCS(parentFileName.substr(posS+1)));
		}
		else
		{
			data.SetParentCell(StringUtils::ConvertMBSToWCS(parentFileName));
		}
	}

	// Allow extra parameters to override parameters from the file.
	AddExtraParameters(data);

	// Add the parent header
	GenericDataHeader hdr;
	FillParentGenericDataHeader(hdr);
	data.GetFileHeader()->GetGenericDataHdr()->AddParent(hdr);

	// Read the probe set names.
	// First try the PSI file for expression and genotyping arrays.
	vector<string> psNames;
	if (ReadProbeSetNames(inFile, libPath, psNames) == false)
	{
		errorCode = UnableToLoadProbeSetNames;
		return false;
	}

	// Determine the maximum probe set name length
	int maxln = 0;
	int psCount = (int) psNames.size();
	for (int ips=0; ips<psCount; ips++)
		maxln = max(maxln, (int) psNames[ips].size());

	// Set the number of entries.
	if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression ||
		inFile.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping ||
		inFile.GetHeader().GetAssayType() == CCHPFileHeader::Universal)
	{
		int nps = inFile.GetHeader().GetNumProbeSets();
		bool hasComp = false;
		if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression)
		{
			CExpressionProbeSetResults *results = inFile.GetExpressionResults(0);
			hasComp = results->m_HasCompResults;
		}
		data.SetEntryCount(nps, maxln, hasComp);
		data.SetBackgroundZoneCnt((int)inFile.GetHeader().GetBackgroundZones().size());
	}
	else if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Resequencing)
	{
		CResequencingResults *reseq = inFile.GetResequencingResults();
		data.SetEntryCount(reseq->GetCalledBasesSize(), maxln);
		data.SetBackgroundZoneCnt((int)inFile.GetHeader().GetBackgroundZones().size());
		data.SetForceCnt(reseq->GetForceCallsSize());
		data.SetOrigCnt(reseq->GetOrigCallsSize());
	}

	// Write the file
	try
	{

		CHPFileWriter writer(data);
		if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression)
		{
			WriteExpressionData(inFile, psNames, writer);
			WriteBackground(inFile, writer);
		}
		else if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping)
		{
			WriteGenotypingData(inFile, psNames, writer);
			WriteBackground(inFile, writer);
		}

		else if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Universal)
		{
			WriteUniversalData(inFile, writer);
			WriteBackground(inFile, writer);
		}

		else if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Resequencing)
		{
			WriteResequencingData(inFile, writer);
			WriteBackground(inFile, writer);
		}
	}
	catch (...)
	{
		errorCode = UnableToWriteTheFile;
		return false;
	}
	return true;
}

bool CHPFileConvertToCalvin::ConvertMas5File(const char *fileName, const char *celFile, const char *libPath, const char* maskFile, const char *newFile)
{
	Clear();

	// Read the input file.
	CCHPFileData inFile;
	inFile.SetFileName(fileName);
	if (inFile.Read() == false)
	{
		errorCode = UnableToOpenChpFile;
		return false;
	}

	// Determine the assay type.
	string assayType;
	if (GetAssayType(assayType, inFile) == false)
	{
		errorCode = InvalidChpFileFormat;
		return false;
	}

	// Check that we can convert the assay type to Calvin
	if (CheckAssayType(inFile) == false)
	{
		errorCode = InvalidAssayType;
		return false;
	}

	// Check the algorithm type
	if (CheckAlgorithmType(inFile.GetHeader().GetAlgName()) == false)
	{
		errorCode = InvalidAlgorithmType;
		return false;
	}

	// Create the new CHP object.
	CHPData data(newFile, assayType);

	AddExtraParameters(data);

	// Add zone algorithm parameters that apply because this was converted from a non-XDA format
	TagValuePairTypeList insertedParams;
	AddZoneAlgParameters(insertedParams);

	// Set the data
	SetHeader(inFile, data, &insertedParams);

	// Set the parent name
	if (parentFileName.empty() == false)
	{
		// Attempt to remove the path.
		string::size_type posBS = parentFileName.rfind('\\');
		string::size_type posS = parentFileName.rfind('/');
		if (posBS != -1)
		{
			data.SetParentCell(StringUtils::ConvertMBSToWCS(parentFileName.substr(posBS+1)));
		}
		else if (posS != -1)
		{
			data.SetParentCell(StringUtils::ConvertMBSToWCS(parentFileName.substr(posS+1)));
		}
		else
		{
			data.SetParentCell(StringUtils::ConvertMBSToWCS(parentFileName));
		}
	}

	// Add the parent header
	GenericDataHeader hdr;
	FillParentGenericDataHeader(hdr);
	data.GetFileHeader()->GetGenericDataHdr()->AddParent(hdr);

	// Read the probe set names.
	// First try the PSI file for expression and genotyping arrays.
	vector<string> psNames;
	if (ReadProbeSetNames(inFile, libPath, psNames) == false)
	{
		errorCode = UnableToLoadProbeSetNames;
		return false;
	}

	// Determine the maximum probe set name length
	int maxln = 0;
	int psCount = (int) psNames.size();
	for (int ips=0; ips<psCount; ips++)
		maxln = max(maxln, (int) psNames[ips].size());

	// Compute the background
	AllZonesInfoType zonesInfo;
	zonesInfo.number_zones = 0;
	zonesInfo.pZones = NULL;
	zonesInfo.smooth_factor = 0.f;

	if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression)
	{
		// Read the parent CEL file.
		FusionCELData celData;
		celData.SetFileName(celFile);
		if(celFile != NULL && strlen(celFile) > 0 && celData.Exists())
		{
			try
			{
				if(celData.Read() == false)
				{
					errorCode = UnableToOpenParentCelFile;
					return false;
				}
			}
			catch(affymetrix_calvin_exceptions::CalvinException&)
			{
				errorCode = UnableToOpenParentCelFile;
				return false;
			}
		}
		else
		{
			errorCode = UnableToOpenParentCelFile;
			return false;
		}

		// Zones are only computed for expression.
		if (this->ComputeBackgroundZoneInfo(celData, libPath, maskFile, zonesInfo) == false)
		{
			errorCode = UnableToReadCdfFile;
			return false;
		}
	}

	// Set the number of entries.
	if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression ||
		inFile.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping ||
		inFile.GetHeader().GetAssayType() == CCHPFileHeader::Universal)
	{
		int nps = inFile.GetHeader().GetNumProbeSets();
		bool hasComp = false;
		if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression)
		{
			CExpressionProbeSetResults *results = inFile.GetExpressionResults(0);
			hasComp = results->m_HasCompResults;
		}
		data.SetEntryCount(nps, maxln, hasComp);

		data.SetBackgroundZoneCnt(zonesInfo.number_zones);
	}
	else if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Resequencing)
	{
		CResequencingResults *reseq = inFile.GetResequencingResults();
		data.SetEntryCount(reseq->GetCalledBasesSize(), maxln);

		data.SetBackgroundZoneCnt(zonesInfo.number_zones);

		data.SetForceCnt(reseq->GetForceCallsSize());
		data.SetOrigCnt(reseq->GetOrigCallsSize());
	}

	// Write the file
	try
	{

		CHPFileWriter writer(data);
		if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression)
		{
			WriteExpressionData(inFile, psNames, writer);
			WriteBackground(zonesInfo, writer);
		}
		else if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping)
		{
			WriteGenotypingData(inFile, psNames, writer);
			//WriteBackground(zonesInfo, writer);
		}

		else if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Universal)
		{
			WriteUniversalData(inFile, writer);
			//WriteBackground(zonesInfo, writer);
		}

		else if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Resequencing)
		{
			WriteResequencingData(inFile, writer);
			//WriteBackground(zonesInfo, writer);
		}
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
 * Command Console file or from a array set GenericDataHeader created on the fly to
 * include the array ID.
 */
void CHPFileConvertToCalvin::FillParentGenericDataHeader(GenericDataHeader& hdr)
{
	if (GetParentGenericDataHeader(hdr) == false)
	{
		if (arrayID.empty() == false || arrayBarcode.empty() == false)
			ConverterCreateAndAddParentArrayGenericDataHeader(hdr, arrayID, arrayBarcode);
	}
}

/*
 * Attempts to read the GenericDataHeader from the parent Command Console file.
 */
bool CHPFileConvertToCalvin::GetParentGenericDataHeader(GenericDataHeader& hdr)
{
	if (this->parentFileName.length() == 0)
		return false;
	
	// Read CEL header
	CelFileData cel;
	CelFileReader reader;
	reader.SetFilename(this->parentFileName.c_str());
	try
	{
		reader.Read(cel);
	}
	catch(affymetrix_calvin_exceptions::CalvinException&)
	{
		return false;
	}

	// Copy GenericDataHeader
	hdr = *cel.GetFileHeader()->GetGenericDataHdr();
	return true;
}

/*
 * Adds any extra parameters from the extraParameter collection to the file header.
 */
void CHPFileConvertToCalvin::AddExtraParameters(CHPData& data)
{
	if (this->extraParameters == NULL)
		return;

	GenericDataHeader* gdh = data.GetFileHeader()->GetGenericDataHdr();

	// Copy the parameters to the file.
	for (ParameterNameValueTypeVector::const_iterator ii = extraParameters->begin(); ii != extraParameters->end(); ++ii)
	{
		gdh->AddNameValParam(*ii);
	}
}

bool CHPFileConvertToCalvin::CheckAlgorithmType(const std::string& algorithmName)
{
	// Disallow conversion of CHP's based on the empirical algorithm.
#ifndef WIN32
#define stricmp strcasecmp
#endif
	if (stricmp(algorithmName.c_str(), "ExpressionCall") == 0)
	{
		return false;
	}
	return true;
#undef stricmp
}

bool CHPFileConvertToCalvin::ValidateXDAFile(const char *fileName, const char *libPath, const char *newFile)
{
	Clear();

	// Read the input file.
	CCHPFileData inFile;
	inFile.SetFileName(fileName);
	if (inFile.ReadHeader() == false)
	{
		errorCode = UnableToOpenChpFile;
		return false;
	}

	// Check that we can convert the assay type to Calvin
	if (CheckAssayType(inFile) == false)
	{
		errorCode = InvalidAssayType;
		return false;
	}

	// First try the PSI file for expression and genotyping arrays.
	vector<string> psNames;
	if (ReadProbeSetNames(inFile, libPath, psNames) == false)
	{
		errorCode = UnableToLoadProbeSetNames;
		return false;
	}

	return true;
}

bool CHPFileConvertToCalvin::ValidateMas5File(const char *fileName, const char *celFile, const char *libPath, const char* maskFile, const char *newFile)
{
	Clear();

	// Read the input file.
	CCHPFileData inFile;
	inFile.SetFileName(fileName);
	if (inFile.ReadHeader() == false)
	{
		errorCode = UnableToOpenChpFile;
		return false;
	}

		// Check that we can convert the assay type to Calvin
	if (CheckAssayType(inFile) == false)
	{
		errorCode = InvalidAssayType;
		return false;
	}

	// First try the PSI file for expression and genotyping arrays.
	vector<string> psNames;
	if (ReadProbeSetNames(inFile, libPath, psNames) == false)
	{
		errorCode = UnableToLoadProbeSetNames;
		return false;
	}

	if (inFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression)
	{
		// Read the parent CEL file.
		FusionCELData celData;
		celData.SetFileName(celFile);
		if(celFile == NULL || strlen(celFile) == 0 || celData.Exists() == false)
		{
			errorCode = UnableToOpenParentCelFile;
			return false;
		}

		if (celData.ReadHeader() == false)
		{
			errorCode = UnableToOpenParentCelFile;
			return false;
		}

		// Read the CDF file.
		// From CExpressionAlgorithmImplementation::ReadCdfFile
		FusionCDFData m_Cdf;
		string cdfFile = Fs::join(libPath ,StringUtils::ConvertWCSToMBS(celData.GetChipType()) + ".CDF");
		m_Cdf.SetFileName(cdfFile.c_str());
		if (m_Cdf.Exists() == false)
		{
			errorCode = UnableToReadCdfFile;
			return false;
		}
	}
	return true;
}
