////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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

/*! 

\page txt2cnvchp_description Writing genotype results to a "multi-data" CHP file.

This program is provided as an example of using the Fusion SDK to write genotype results to an Affymetrix
"multi-data" CHP file. The multi-data CHP file was designed to store a variety of Affymetrix results,
including genotypes and copy number. The Genotyping Console software currently uses the "multi-data" CHP
file for storing genotype results from the BRLMM-P and Birdseed algorithms as well as copy number results
from its copy number algorithms.

The program is implemented as a command line application. The arguments specify items such as the TSV file
containing the genotype results, algorithm parameters, summary statistics, probe array type, parent CEL file, 
and other items that are stored in the CHP file header.

 */

#ifdef _MSC_VER
#include <windows.h>
#endif
//
#include "calvin_files/converters/chp/src/CmdLine.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileBufferWriter.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
#include "calvin_files/writers/src/GenericDataHeaderUpdater.h"
//
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
//

#ifdef _MSC_VER
#pragma warning(disable: 4996)
#endif

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;


static void GetCNVCHPParentHeaderFromCEL(string celFile, CHPMultiDataData *data)
{
	// Store the CEL header
	if (celFile.length() > 0 && FileUtils::Exists(celFile.c_str()) == true)
	{
		FusionCELData cel;
		cel.SetFileName(celFile.c_str());
        cel.ReadHeader();
	    GenericData *gdata = cel.GetGenericData();
	    if (gdata != NULL)
			data->GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr()); 
	    cel.Close();
	}
}

static void GetCNVCHPParentHeaderFromCNCHP(string cnvchpFile, CHPMultiDataData *data)
{
	// Store the CNCHP header
	if (cnvchpFile.length() > 0 && FileUtils::Exists(cnvchpFile.c_str()) == true)
	{
		FusionCHPData *chp = FusionCHPDataReg::ReadHeader(cnvchpFile);
        FusionCHPMultiDataData *dchp = FusionCHPMultiDataData::FromBase(chp);
        if (dchp != NULL)
		{
			GenericData *gdata = dchp->GetGenericData();
			if (gdata != NULL)
			{
				GenericDataHeader pHeader = *gdata->Header().GetGenericDataHdr();
				if(pHeader.GetParentCnt() > 0)
				{
					data->GetFileHeader()->GetGenericDataHdr()->AddParent(pHeader.GetParent(0));
				}
				else
				{
					throw 0;
				}
			}
		}
        delete chp;
	}
}

/*! Create a "multi-data" CHP file with just the header information. The remainder of the file
 * will be created at a later time using the buffer writer technique.
 * The CHP file will contain only "genotyping" results.
 * @param execId The execution identifier. This identifier is used to identify the batch run that created the CHP files.
 * @param celFile The full path to the parent CEL file. The header of the CEL file is copied to the CHP file.
 * @param outFile The name of the output CHP file.
 * @param extraColNames The names of the extra data columns. Should not include probe set name, call and confidence columns.
 * @param extraColTypes The types (float, int, ubyte) of the extra columns.
 * @param numEntries The number of rows (entries) of results to store in the CHP file.
 * @param maxProbeSetNameLength The maximum length of the probe set names.
 * @param algName The name of the algorithm used to create the results.
 * @param algVersion The algorithm version.
 * @param chipType the chip type, also known as the probe array type.
 * @param programName The name of the program used to create the CHP file.
 * @param programVersion The version of the program.
 * @param programCompany The company or institution who developed the CHP creating software.
 * @param paramNames A list of parameter names to store in the CHP file header.
 * @param paramValues A list of parameter values to store in the CHP file header.
 * @param sumNames A list of summary statistic names to store in the CHP file header.
 * @param sumValues A list of summary statistic values to store in the CHP file header.
 * @param execGUID A AffymetrixGUID type for holding the execGUID. This is so each batch of files processed all have the same exec guid.
*/
static void CreateFileWithHeader
(
	const string &execId,
	const string &celFile,
	const string &outFile,
	const vector<string>& extraColNames,
	const vector<string>& extraColTypes,
	unsigned long numEntries,
	int maxProbeSetNameLength,
	const string &algName,
	const string &algVersion,
	const string &chipType,
	const string &programName,
	const string &programVersion,
	const string &programCompany,
	const string &mapFile,
	const string &regionsFile,
	const vector<string>& paramNames,
	const vector<string>& paramValues,
	const vector<string>& sumNames,
	const vector<string>& sumValues,
	const vector<string>& extraNames,
	const vector<string>& extraValues,
	const MultiDataType &chpDataType,
	const ParameterNameValueTypeList& datasetParams //BO,
//BO	const AffymetrixGuidType& execGUID
)
{
	// Create the vector of extra columns. The sample code here supports only float, 32 bit integers and 8 bit unsigned integers.
	vector<ColumnInfo> extraColumns;
	int ncols = (int)extraColNames.size();
	for (int icol=0; icol<ncols; icol++)
	{
		if (extraColTypes[icol] == "float")
		{
			FloatColumn fcol(StringUtils::ConvertMBSToWCS(extraColNames[icol]));
			extraColumns.push_back(fcol);
		}
		else if (extraColTypes[icol] == "int")
		{
			IntColumn intcol(StringUtils::ConvertMBSToWCS(extraColNames[icol]));
			extraColumns.push_back(intcol);
		}
		else if (extraColTypes[icol] == "ubyte")
		{
			UByteColumn ubcol(StringUtils::ConvertMBSToWCS(extraColNames[icol]));
			extraColumns.push_back(ubcol);
		}
		else
		{
			throw string("Unsupported column type: ") + extraColTypes[icol];
		}
	}

	// Create the data object
	CHPMultiDataData *data = new CHPMultiDataData(outFile);
    data->SetEntryCount(chpDataType, numEntries, maxProbeSetNameLength, extraColumns);
	data->SetAlgName(StringUtils::ConvertMBSToWCS(algName));
	data->SetAlgVersion(StringUtils::ConvertMBSToWCS(algVersion));
	data->SetArrayType(StringUtils::ConvertMBSToWCS(chipType));

	if(chpDataType == CopyNumberVariationMultiDataType)
	{
		string cnvchpFile = celFile;
		GetCNVCHPParentHeaderFromCNCHP(cnvchpFile, data);
	}
	else
	{
		GetCNVCHPParentHeaderFromCEL(celFile, data);
	}

	// Add algorithm parameters to list.
    ParameterNameValueTypeList params;
    ParameterNameValueType param;
    
	if (programName.empty() == false)
	{
		param.SetName(L"program-name");
		param.SetValueText(StringUtils::ConvertMBSToWCS(programName));
		data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
	}

	if (programVersion.empty() == false)
	{
		param.SetName(L"program-version");
		param.SetValueText(StringUtils::ConvertMBSToWCS(programVersion));
		data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
	}

	if (programCompany.empty() == false)
	{
		param.SetName(L"program-company");
		param.SetValueText(StringUtils::ConvertMBSToWCS(programCompany));
		data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
	}

	if(chpDataType == CopyNumberVariationMultiDataType)
	{
		//map file
		param.SetName(L"affymetrix-param-map-file");
        param.SetValueAscii(mapFile);
		data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);

		//regions file
		param.SetName(L"affymetrix-param-regions-file");
        param.SetValueAscii(regionsFile);
		data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);

		//apt exec guid
		/*param.SetName(L"affymetrix-algorithm-apt-exec-guid");
		param.SetValueAscii(AffymetrixGuid::GenerateNewGuid());
		data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);*/
	}

	int nparams = (int) extraNames.size();
	for(int iparam=0; iparam<nparams; iparam++)
	{
		param.SetName(StringUtils::ConvertMBSToWCS(extraNames[iparam]));
        param.SetValueAscii(extraValues[iparam]);
        data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
	}

	nparams = (int) paramNames.size();
	param.SetName(L"exec-guid");
//BO - Just use the execID passed in on the cmdline
//	if(chpDataType == CopyNumberVariationMultiDataType)
//	{
//		param.SetValueAscii(execGUID);
//		param.SetValueAscii(AffymetrixGuid::GenerateNewGuid());
//	}
//	else
//	{
		param.SetValueAscii(execId);
//	}

	params.push_back(param);
	for(int iparam=0; iparam<nparams; iparam++)
	{
		param.SetName(StringUtils::ConvertMBSToWCS(paramNames[iparam]));
        param.SetValueAscii(paramValues[iparam]);
        params.push_back(param);
	}
	if (params.empty() == false)
		data->AddAlgParams(params);

	params.clear();
	nparams = (int) sumNames.size();
	for(int iparam=0; iparam<nparams; iparam++)
	{
		param.SetName(StringUtils::ConvertMBSToWCS(sumNames[iparam]));
        param.SetValueAscii(sumValues[iparam]);
        params.push_back(param);
	}
	if (params.empty() == false)
		data->AddSummaryParams(params);
	
	//Add dataset parameters specific to the multi data chip type
	if(datasetParams.empty() == false)
	{
		DataSetHeader *dsh = data->GetDataSetHeader(chpDataType);
		for (ParameterNameValueTypeList::const_iterator it=datasetParams.begin(); it!=datasetParams.end(); it++)
            dsh->AddNameValParam(*it);
	}

	// Creating the writer object will create the file with the header information.
	CHPMultiDataFileWriter writer(*data);
}

/*! Convert a call to a code
 * @param call The call, or can be a code
 * @return The code
 */
static int CallToCode(const string &call)
{
	if (isdigit(call[0]) != 0)
	{
		return atoi(call.c_str());
	}
	else
	{
		if (call == "AA" || call == "aa")
			return SNP_AA_CALL;
		else if (call == "AB" || call == "ab")
			return SNP_AB_CALL;
		else if (call == "BB" || call == "bb")
			return SNP_BB_CALL;
		else
			return SNP_NO_CALL;
	}
}

/*! Parse a line into a list.
 * @param line The line of data
 * @return A list of strings from the line (separated by tabs).
 */
static list<string> ParseLine(char *line)
{
	list<string> tokens;
	const char *sep = "\t";
	char *token = strtok(line, sep);
	while (token != NULL)
	{
		tokens.push_back(token);
		token = strtok(NULL, sep);
	}
	return tokens;	
}

/*! Adds the genotype results to the CHP file. This function uses the buffer writer object
 * to write the data to the CHP file.
 * @param inFile The name of the input TSV file.
 * @param outFile The name of the output CHP file.
 * @param maxProbeSetNameLength The maximum length of the probe set names.
 * @param extraColTypes The types (float, int, ubyte) of the extra columns.
 */
static void AddFileBody(const string &inFile, const string &outFile, int maxProbeSetNameLength, const vector<string>& extraColTypes)
{
	// Create the buffer writer object
	CHPMultiDataFileBufferWriter genotypeEntryBufferWriter;
	vector<string> fileNames;
	fileNames.push_back(outFile);
    vector<MultiDataType> dataTypes;
    dataTypes.push_back(GenotypeMultiDataType);
    map<MultiDataType, int> maxLen;
    maxLen[GenotypeMultiDataType] = maxProbeSetNameLength;
	genotypeEntryBufferWriter.Initialize(&fileNames, dataTypes, maxLen);

	// Open the file
	ifstream instr(inFile.c_str(), ios::in);
	if (!instr)
		throw string("Unable to read the input file");

	// First line with header
	char line[1024];
	instr.getline(line, 1024);

	// Get the data
	maxProbeSetNameLength = 0;

	ProbeSetMultiDataGenotypeData entry;
	while (instr.getline(line, 1024))
	{
		list<string> cols = ParseLine(line);
		list<string>::const_iterator colIt = cols.begin();

		entry.name = *colIt;
		++colIt;
		entry.call = CallToCode(*colIt);
		++colIt;
		entry.confidence = (float)atof(colIt->c_str());
		++colIt;

		entry.metrics.resize(cols.size() - 3);
		int colIdx = 0;
		for (; colIt!=cols.end(); colIt++, colIdx++)
		{
			if (extraColTypes[colIdx] == "float")
				entry.metrics[colIdx].SetValueFloat((float)atof(colIt->c_str()));
			else if (extraColTypes[colIdx] == "int")
				entry.metrics[colIdx].SetValueInt32(atoi(colIt->c_str()));
			else if (extraColTypes[colIdx] == "ubyte")
				entry.metrics[colIdx].SetValueUInt8(CallToCode(colIt->c_str()));
		}
		genotypeEntryBufferWriter.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, 0, entry);
	}

	//Flush the buffer. The destructor will close the file.
	genotypeEntryBufferWriter.FlushBuffer();
}

/*! Adds the copy number results to the CHP file. This function uses the buffer writer object
 * to write the data to the CNCHP file.
 * @param inFile The name of the input TSV file.
 * @param outFile The name of the output CNCHP file.
 * @param maxProbeSetNameLength The maximum length of the probe set names.
 * @param extraColTypes The types (float, int, ubyte) of the extra columns.
 */
static void AddFileBodyCNCHP(const string &inFile, const string &outFile, int maxProbeSetNameLength, const vector<string>& extraColTypes)
{
	// Create the buffer writer object
	CHPMultiDataFileBufferWriter copyNumberEntryBufferWriter;
	vector<string> fileNames;
	fileNames.push_back(outFile);
    vector<MultiDataType> dataTypes;
    dataTypes.push_back(CopyNumberMultiDataType);
    map<MultiDataType, int> maxLen;
    maxLen[CopyNumberMultiDataType] = maxProbeSetNameLength;
	copyNumberEntryBufferWriter.Initialize(&fileNames, dataTypes, maxLen);
	
	// Open the file
	ifstream instr(inFile.c_str(), ios::in);
	if (!instr)
		throw string("Unable to read the input file");

	// First line with header
	char line[1024];
	instr.getline(line, 1024);

	// Get the data
	maxProbeSetNameLength = 0;

	ProbeSetMultiDataCopyNumberData entry;
	while (instr.getline(line, 1024))
	{
		list<string> cols = ParseLine(line);
		list<string>::const_iterator colIt = cols.begin();

		entry.name = *colIt;
		++colIt;
		entry.chr = (u_int8_t)atof(colIt->c_str());
		++colIt;
		entry.position = (u_int32_t)atof(colIt->c_str());
		++colIt;

		entry.metrics.resize(cols.size() - 3);
		int colIdx = 0;
		for (; colIt!=cols.end(); colIt++, colIdx++)
		{
			if (extraColTypes[colIdx] == "float")
				entry.metrics[colIdx].SetValueFloat((float)atof(colIt->c_str()));
			else if (extraColTypes[colIdx] == "int")
				entry.metrics[colIdx].SetValueInt32(atoi(colIt->c_str()));
			else if (extraColTypes[colIdx] == "ubyte")
				entry.metrics[colIdx].SetValueUInt8(CallToCode(colIt->c_str()));

		}
 		copyNumberEntryBufferWriter.WriteMultiDataCopyNumberEntry(CopyNumberMultiDataType, 0, entry);
	}

	//Flush the buffer. The destructor will close the file.
	copyNumberEntryBufferWriter.FlushBuffer();
}

static float GetNumFromCalculation(long start, long end, time_t time, int row)
{
	long num = end - start;
	if(num <= 0)
		num = 100;
	return (time / (float)num)/++row;
}

static float GetSignal(float signal)
{
	int nSignal = signal;
	signal = signal - nSignal;
	nSignal = nSignal % 3;
	signal = nSignal + signal ;
	if(signal < 0.25)
		signal = 0.25;
	if(signal > 3)
		signal = 3;
	return signal;
}

static int GetCall(float call)
{
	int nCall = call;
	nCall = nCall % 5;
	if(nCall < 0)
		nCall = 0;
	if(nCall > 4)
		nCall = 4;
	return nCall;
}

static float GetAllele(float allele)
{
	int nAllele = allele;
	allele = allele - nAllele;
	if(allele == 0.0)
		allele = (nAllele % 2) * 1.0;
	if(allele < 0)
		allele = 0;
	if(allele > 1)
		allele = 1;
	return allele;
}

static float GetConfident(float confident)
{
	int nConfident = confident;
	confident = confident - nConfident;
	confident = confident / 1.8;
	if(confident < 0)
		confident = 0;
	if(confident > 0.5)
		confident = 0.5;
	return confident;
}

/*! Adds the copy number results to the CHP file. This function uses the buffer writer object
 * to write the data to the CNVCHP file.
 * @param inFile The name of the input TSV file.
 * @param outFile The name of the output CNCHP file.
 * @param maxProbeSetNameLength The maximum length of the probe set names.
 * @param extraColTypes The types (float, int, ubyte) of the extra columns.
 */
static void AddFileBodyCNVCHP(const string &inFile, const string &outFile, int maxProbeSetNameLength, const vector<string>& extraColTypes)
{
	// Create the buffer writer object
	CHPMultiDataFileBufferWriter copyNumberEntryBufferWriter;
	vector<string> fileNames;
	fileNames.push_back(outFile);
    vector<MultiDataType> dataTypes;
    dataTypes.push_back(CopyNumberVariationMultiDataType);
    map<MultiDataType, int> maxLen;
    maxLen[CopyNumberVariationMultiDataType] = maxProbeSetNameLength;
	copyNumberEntryBufferWriter.Initialize(&fileNames, dataTypes, maxLen);
	
	time_t seconds = time(NULL);
	int clk = clock();
	if(clk == 0) clk = 1;
	seconds = (seconds/1000) * clk;

	// Open the file
	ifstream instr(inFile.c_str(), ios::in);
	if (!instr)
		throw string("Unable to read the input file");

	// First line with header
	char line[5000];
	instr.getline(line, 5000);

	// Get the data
	maxProbeSetNameLength = 0;

	ProbeSetMultiDataCopyNumberVariationRegionData entry;
	int row = 0;
	while (instr.getline(line, 5000))
	{
		list<string> cols = ParseLine(line);
		list<string>::const_iterator colIt = cols.begin();
		entry.name = *colIt;
		++colIt;
		++colIt;
		long start = atol(colIt->c_str());
		++colIt;
		long end = atol(colIt->c_str());
		float randomNum = GetNumFromCalculation(start, end, seconds, row);
		entry.signal = GetSignal(randomNum);
		entry.call = (u_int8_t)GetCall(randomNum);
		entry.confidenceScore = GetConfident(randomNum);
		copyNumberEntryBufferWriter.WriteMultiDataCopyNumberVariationRegionEntry(CopyNumberVariationMultiDataType, 0, entry);
		row++;
	}

	//Flush the buffer. The destructor will close the file.
	copyNumberEntryBufferWriter.FlushBuffer();
}

/*! Reads the genotype results from a simple TSV text file.
 * @param inFile The name of the input TSV file.
 * @param maxProbeSetNameLength The maximum length of the probe set names.
 */
static unsigned long ReadData(const string &inFile, int &maxProbeSetNameLength)
{
	// Open the file
	ifstream instr(inFile.c_str(), ios::in);
	if (!instr)
		throw string("Unable to read the input file");

	// First line with header
	char line[1024];
	instr.getline(line, 1024);

	// Get the data
	maxProbeSetNameLength = 0;

	// count keeps track of the number of rows. This now handles what used to be the data list size
	// given we do not read and store the data in the list anymore, we needed this value for later
	// see the main function for how count is used from the return of this fuction.
	unsigned long count = 0;
	while (instr.getline(line, 1024))
	{
		list<string> cols = ParseLine(line);
		maxProbeSetNameLength = max(maxProbeSetNameLength, (int) cols.begin()->size());
		++count;
	}
	return count;
}

/*
 * Store the start index and count of probe sets for the given chromosome.
 */
static void SetChromosomeProbeSetIndexInformation(ParameterNameValueTypeList &params, int chr, unsigned long start, unsigned long count)
{
	ostringstream str;
	str << (int) chr;
	wstring schr = StringUtils::ConvertMBSToWCS(str.str());
	ParameterNameValueType param;
	param.SetName(schr + L":start");
	param.SetValueInt32(start);
	params.push_back(param);
	param.SetName(schr + L":count");
	param.SetValueInt32(count);
	params.push_back(param);
	param.SetName(schr + L":display");
	param.SetValueAscii(ChromosomeToString(chr));
	params.push_back(param);
}

/*! Reads the copy number results from a simple TSV text file.
 * @param inFile The name of the input TSV file.
 * @param maxProbeSetNameLength The maximum length of the probe set names.
 * @param datasetParams The list of parameters that form the part of the dataset header.
 */
static unsigned long ReadCNData(const string &inFile, int &maxProbeSetNameLength, ParameterNameValueTypeList &datasetParams)
{
	// Open the file
	ifstream instr(inFile.c_str(), ios::in);
	if (!instr)
		throw string("Unable to read the input file");

	// First line with header
	char line[1024];
	instr.getline(line, 1024);

	// Get the data
	maxProbeSetNameLength = 0;

	// count keeps track of the number of rows. This now handles what used to be the data list size
	// given we do not read and store the data in the list anymore, we needed this value for later
	// see the main function for how count is used from the return of this fuction.
	unsigned long count = 0;
	int lastChr = -1;
	unsigned long chrRowCount = 0;
	unsigned long start = 0;
	while (instr.getline(line, 1024))
	{
		list<string> cols = ParseLine(line);
		list<string>::iterator colsIt = cols.begin();
		maxProbeSetNameLength = max(maxProbeSetNameLength, (int) colsIt->size());
		++colsIt;
		int chr = atoi(colsIt->c_str());
		// If there is a new chromosome, save the details of the previous one.
		if(lastChr != chr)
		{
			if(lastChr != -1)
				SetChromosomeProbeSetIndexInformation(datasetParams, lastChr, start, chrRowCount);
			
			lastChr = chr;
			start = count;
			chrRowCount = 0;
		}
		++chrRowCount;
		++count;
	}
	SetChromosomeProbeSetIndexInformation(datasetParams, lastChr, start, chrRowCount);

	return count;
}

/*! Reads the copy number variation results from a simple TSV text file.
 * @param inFile The name of the input TSV file.
 * @param maxProbeSetNameLength The maximum length of the probe set names.
 * @param datasetParams The list of parameters that form the part of the dataset header.
 */
static unsigned long ReadCNVData(const string &inFile, int &maxProbeSetNameLength, ParameterNameValueTypeList &datasetParams)
{
	// Open the file
	ifstream instr(inFile.c_str(), ios::in);
	if (!instr)
		throw string("Unable to read the input file");

	// First line with header
	char line[5000];
	instr.getline(line, 5000);

	// Get the data
	maxProbeSetNameLength = 0;

	// count keeps track of the number of rows. This now handles what used to be the data list size
	// given we do not read and store the data in the list anymore, we needed this value for later
	// see the main function for how count is used from the return of this fuction.
	unsigned long count = 0;
	int lastChr = -1;
	unsigned long chrRowCount = 0;
	unsigned long start = 0;
	while (instr.getline(line, 5000))
	{
		list<string> cols = ParseLine(line);
		list<string>::iterator colsIt = cols.begin();
		maxProbeSetNameLength = max(maxProbeSetNameLength, (int) colsIt->size());
		++count;
	}

	return count;
}

/*! Show a help message on the use of the command line parameters. */
static void show_help()
{
	const char *tab = "\t";
    cout << "The cnvchp files will be created in the same folder as the source cnchp files." << endl;
    cout << "The following details the parameters" << endl;

cout << "Must be set as shown; the files are in the txt2cnvchp folder in cvs" << endl;
cout << "-execid 1234 	This is the execution guid assigned to all of the cnvchps.  Can be any string valid for a guid" << endl;
cout << "-cel \\\\10.80.199.147\\share\\DataFromTest\\b475CN		This is the folder of your cnchp or cel files (I have not tested using cel files as a source)" << endl; 	 
cout << "-in C:\\cnvMapFile\\CNP_regions.tsv 			This must be this file.   Change to the folder that is correct for you" << endl;
cout << "-arrayType GenomeWideSNP_6 				Must be this text" << endl;
cout << "-chpType CopyNumberVariation 				Must be this text	" << endl;			
cout << "-regionsFile C:\\cnvMapFile\\cnp_map.b36.txt This must be this file.  Change to the folder that is correct for you" << endl;
cout << "-mapFile C:\\cnvMapFile\\CNP_regions.tsv 	This must be this file.  Change to the folder that is correct for you" << endl;


cout << "These can be any text" << endl;
cout << "-out C:\\CnChpdata1k.chp" << endl; 
cout << "-algName TXT2MCHP" << endl;
cout << "-algVersion version1.0 " << endl;
cout << "-programName GTC3.0" << endl; 
cout << "-programVersion 0.1" << endl; 
cout << "-programCompany Affymetrix " << endl;
	
}

static void ReadFile(const string &celFile)
{
	// Store the CEL header
	if (celFile.length() > 0 && FileUtils::Exists(celFile.c_str()) == true)
	{
		FusionCELData cel;
		cel.SetFileName(celFile.c_str());
        cel.ReadHeader();
	    GenericData *gdata = cel.GetGenericData();
	    if (gdata != NULL)
		{
			*gdata->Header().GetGenericDataHdr(); 
		}
	    cel.Close();
	}
}

static string GetOutputFile(string cnchpFile)
{	
	string chchpFileToLower = cnchpFile;
	int pos = cnchpFile.rfind('\\');
	string filePath = cnchpFile.substr(0, pos);
	string fileName = cnchpFile.substr(pos + 1, cnchpFile.length() - pos);
	pos = fileName.find('.');
	fileName = fileName.substr(0, pos);
	return filePath + "\\" + fileName + ".CNVCHP";
}

static void GetFiles(string *dirPath, string searchExtension, vector<string>& list)
{
	string path = dirPath->c_str();
	dirPath->append("\\").append(searchExtension);
	list.clear();
    WIN32_FIND_DATA fd;
	HANDLE findHandle = ::FindFirstFileA(dirPath->c_str(), &fd);
    if(findHandle == INVALID_HANDLE_VALUE)
        return;
    do
    {
        if(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
        {
            continue;
        }
		list.push_back(path + "\\" + fd.cFileName);
    } while(::FindNextFileA(findHandle, &fd) != 0);

    ::FindClose(findHandle);
}

/*! The main entrance to the program. */
int main(int argc, char * argv[])
{
	try
	{

		// Parse the command line arguments.
		CCmdLine cmdLine;
		if (cmdLine.SplitLine(argc, argv) < 1)
		{
			show_help();
			exit(-1);
		}

		string mapFile = "";
		string regionsFile = "";

		string execId = cmdLine.GetArgument("-execid", 0);
		string celFile = cmdLine.GetSafeArgument("-cel", 0, "");
		string inFile = cmdLine.GetArgument("-in", 0);
		string chpType = cmdLine.GetArgument("-chpType", 0);
		string arrayType = cmdLine.GetArgument("-arrayType", 0);
		string algName = cmdLine.GetArgument("-algName", 0);
		string algVersion = cmdLine.GetArgument("-algVersion", 0);
		string programName = cmdLine.GetArgument("-programName", 0);
		string programVersion = cmdLine.GetArgument("-programVersion", 0);
		string programCompany = cmdLine.GetArgument("-programCompany", 0);

		if(chpType == "CopyNumberVariation")
		{
            mapFile = "c:\\cnv\\cnvmap.hg18.bed";
			regionsFile = "c:\\cnv\\cnvregions.hg18.regions";
		}

		vector<string> colNames;
		int n = cmdLine.GetArgumentCount("-colNames");
		for (int i=0; i<n; i++)
			colNames.push_back(cmdLine.GetArgument("-colNames", i));

		vector<string> colTypes;
		n = cmdLine.GetArgumentCount("-colTypes");
		for (int i=0; i<n; i++)
			colTypes.push_back(cmdLine.GetArgument("-colTypes", i));

		vector<string> paramNames;
		n = cmdLine.GetArgumentCount("-paramNames");
		for (int i=0; i<n; i++)
			paramNames.push_back(cmdLine.GetArgument("-paramNames", i));

		vector<string> paramValues;
		n = cmdLine.GetArgumentCount("-paramValues");
		for (int i=0; i<n; i++)
			paramValues.push_back(cmdLine.GetArgument("-paramValues", i));

		vector<string> sumNames;
		n = cmdLine.GetArgumentCount("-sumNames");
		for (int i=0; i<n; i++)
			sumNames.push_back(cmdLine.GetArgument("-sumNames", i));

		vector<string> sumValues;
		n = cmdLine.GetArgumentCount("-sumValues");
		for (int i=0; i<n; i++)
			sumValues.push_back(cmdLine.GetArgument("-sumValues", i));

		vector<string> extraNames;
		n = cmdLine.GetArgumentCount("-extraNames");
		for (int i=0; i<n; i++)
			extraNames.push_back(cmdLine.GetArgument("-extraNames", i));

		vector<string> extraValues;
		n = cmdLine.GetArgumentCount("-extraValues");
		for (int i=0; i<n; i++)
			extraValues.push_back(cmdLine.GetArgument("-extraValues", i));

		vector<string> inputFiles;
		if(chpType == "CopyNumberVariation")
		{
			string cnchpDir = celFile;
			GetFiles(&cnchpDir, "*.cnchp", inputFiles);
		}
		else
		{
			inputFiles.push_back(celFile);
		}

		for(int iFile = 0; iFile < inputFiles.size(); iFile++)
		{
			celFile = inputFiles[iFile];
			string outFile = "";
			if(chpType == "CopyNumberVariation")
				outFile = GetOutputFile(celFile);
			else
				outFile = cmdLine.GetArgument("-out", 0);

			// Read the TSV input file
			int maxProbeSetNameLength = 0;
			MultiDataType chpDataType ;
			if(chpType == "Genotype" ||  chpType == "")
				chpDataType = GenotypeMultiDataType;
			else if(chpType == "CopyNumber")
				chpDataType = CopyNumberMultiDataType;
			else if(chpType == "CopyNumberVariation")
				chpDataType = CopyNumberVariationMultiDataType;

			// data collection is no longer used, so we have ReadData return the row count so that we can pass it into
			// CreateFileWithHeader call
			unsigned long rowCount =0;
			ParameterNameValueTypeList datasetParams;

			if(chpDataType == GenotypeMultiDataType)
				rowCount = ReadData(inFile, maxProbeSetNameLength);
			else if(chpDataType == CopyNumberMultiDataType)
				rowCount = ReadCNData(inFile, maxProbeSetNameLength, datasetParams);
			else if(chpDataType == CopyNumberVariationMultiDataType)
				rowCount = ReadCNVData(inFile, maxProbeSetNameLength, datasetParams);
			
			// Create the CHP file.
			// Creates the header
			CreateFileWithHeader(execId, celFile, outFile, colNames, colTypes, rowCount, maxProbeSetNameLength,
				algName, algVersion, arrayType, programName, programVersion, programCompany, mapFile, regionsFile, paramNames, paramValues,
				sumNames, sumValues, extraNames, extraValues, chpDataType, datasetParams);//BO, execGUID);

			//// add the body (data)
			if(chpDataType == GenotypeMultiDataType)
				AddFileBody(inFile,outFile, maxProbeSetNameLength, colTypes);
			else if(chpDataType == CopyNumberMultiDataType)
				AddFileBodyCNCHP(inFile,outFile, maxProbeSetNameLength, colTypes);
			else if(chpDataType == CopyNumberVariationMultiDataType)
				AddFileBodyCNVCHP(inFile,outFile, maxProbeSetNameLength, colTypes);
		}
	}
	catch (string s)
	{
		cout << s << endl;
		show_help();
	}
	catch (int e)
	{
		cout << "Invalid argument" << endl;
		show_help();
	}
	catch (...)
	{
		cout << "Unknown error" << endl;
		show_help();
	}
	return 0;
}
