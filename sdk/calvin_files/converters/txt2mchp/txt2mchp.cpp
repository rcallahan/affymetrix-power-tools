////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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

\page txt2mchp_description Writing genotype results to a "multi-data" CHP file.

This program is provided as an example of using the Fusion SDK to write genotype results to an Affymetrix
"multi-data" CHP file. The multi-data CHP file was designed to store a variety of Affymetrix results,
including genotypes and copy number. The Genotyping Console software currently uses the "multi-data" CHP
file for storing genotype results from the BRLMM-P and Birdseed algorithms as well as copy number results
from its copy number algorithms.

The program is implemented as a command line application. The arguments specify items such as the TSV file
containing the genotype results, algorithm parameters, summary statistics, probe array type, parent CEL file, 
and other items that are stored in the CHP file header.

 */
#include "calvin_files/converters/utils/src/CmdLine.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileBufferWriter.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"

#include <vector>
#include <list>
#include <iostream>
#include <fstream>

#ifdef _MSC_VER
#pragma warning(disable: 4996)
#endif

using namespace std;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;

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
	const vector<string>& paramNames,
	const vector<string>& paramValues,
	const vector<string>& sumNames,
	const vector<string>& sumValues,
	const vector<string>& extraNames,
	const vector<string>& extraValues
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
    data->SetEntryCount(GenotypeMultiDataType, numEntries, maxProbeSetNameLength, extraColumns);
	data->SetAlgName(StringUtils::ConvertMBSToWCS(algName));
	data->SetAlgVersion(StringUtils::ConvertMBSToWCS(algVersion));
	data->SetArrayType(StringUtils::ConvertMBSToWCS(chipType));

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

	int nparams = (int) extraNames.size();
	for(int iparam=0; iparam<nparams; iparam++)
	{
		param.SetName(StringUtils::ConvertMBSToWCS(extraNames[iparam]));
        param.SetValueAscii(extraValues[iparam]);
        data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
	}

	nparams = (int) paramNames.size();
	param.SetName(L"exec-guid");
	param.SetValueAscii(execId);
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

/*! Show a help message on the use of the command line parameters. */
static void show_help()
{
	const char *tab = "\t";
	cout << "The command line arguments include" << endl;
	cout << tab << "-execid <the execution id>" << endl;
	cout << tab << "-cel <the name of the parent cel file>" << endl;
	cout << tab << "-in <the name of the input data file, TSV with 1 header row>" << endl;
	cout << tab << "-out <the name of the output CHP file>" << endl;
	cout << tab << "-arrayType <the probe array type>" << endl;
	cout << tab << "-algName <the algorithm name>" << endl;
	cout << tab << "-algVersion <the algorithm version>" << endl;
	cout << tab << "-programName <the program name>" << endl;
	cout << tab << "-programVersion <the program version>" << endl;
	cout << tab << "-programCompany <the name of the company>" << endl;
	cout << tab << "-colNames <a list of column names in the input file, not including" << endl;
	cout << tab << "		   the probe set name, call and confidence>" << endl;
	cout << tab << "-colTypes <a list of column types, either \"float\" or \"int\" or \"ubyte\">" << endl;
	cout << tab << "-paramNames <a list of parameter names>" << endl;
	cout << tab << "-paramValues <a list of parameter values>" << endl;
	cout << tab << "-sumNames <a list of summary parameter names>" << endl;
	cout << tab << "-sumValues <a list of summary parameter values>" << endl;
	cout << tab << "-extraNames <a list of extra parameter names>" << endl;
	cout << tab << "-extraValues <a list of extra parameter values>" << endl;
	cout << endl;
	cout << "The input file must have the probe set name, call (AA, AB, BB, No Call)" << endl <<
			"and confidence as the first three columns. All other columns must be of" << endl <<
			"either floating point of integer type. Use the colTypes to specify the" << endl <<
			"type. The paramNames and paramValues arguments are used to specify" << endl <<
			"algorithm parameter." << endl;
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

		string execId = cmdLine.GetArgument("-execid", 0);
		string celFile = cmdLine.GetSafeArgument("-cel", 0, "");
		string inFile = cmdLine.GetArgument("-in", 0);
		string outFile = cmdLine.GetArgument("-out", 0);
		string arrayType = cmdLine.GetArgument("-arrayType", 0);
		string algName = cmdLine.GetArgument("-algName", 0);
		string algVersion = cmdLine.GetArgument("-algVersion", 0);
		string programName = cmdLine.GetArgument("-programName", 0);
		string programVersion = cmdLine.GetArgument("-programVersion", 0);
		string programCompany = cmdLine.GetArgument("-programCompany", 0);

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

		// Read the TSV input file
		int maxProbeSetNameLength = 0;

		// data collection is no longer used, so we have ReadData return the row count so that we can pass it into
		// CreateFileWithHeader call
		unsigned long rowCount = ReadData(inFile, maxProbeSetNameLength);

		// Create the CHP file.
		// Creates the header
		CreateFileWithHeader(execId, celFile, outFile, colNames, colTypes, rowCount, maxProbeSetNameLength,
			algName, algVersion, arrayType, programName, programVersion, programCompany, paramNames, paramValues,
			sumNames, sumValues, extraNames, extraValues);

		// add the body (data)
		AddFileBody(inFile,outFile, maxProbeSetNameLength, colTypes);
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
