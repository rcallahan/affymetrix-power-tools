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

\page txt2echp_description Writing expression results to a "quantification" CHP file.

This program is provided as an example of using the Fusion SDK to write expression results to an Affymetrix
"quantification" CHP file. The quantification CHP file was designed to store a set of quantification values
with their associated probe set names. The Expression Console software currently uses the "quantification" CHP
file for storing results from the RMA and PLIER algorithms.

The program is implemented as a command line application. The arguments specify items such as the TSV file
containing the expression results, algorithm parameters, summary statistics, probe array type, parent CEL file, 
and other items that are stored in the CHP file header.

 */

#include "calvin_files/converters/utils/src/CmdLine.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationFileBufferWriter.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationFileWriter.h"
#include "calvin_files/writers/src/GenericDataHeaderUpdater.h"
//
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

/*! Create a "quantification" CHP file with just the header information. The remainder of the file
 * will be created at a later time using the buffer writer technique.
 * The CHP file will contain only "quantification" results from an expression analysis.
 * @param execId The execution identifier. This identifier is used to identify the batch run that created the CHP files.
 * @param celFile The full path to the parent CEL file. The header of the CEL file is copied to the CHP file.
 * @param outFile The name of the output CHP file.
 * @param probeSetNames The probe set names.
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
	const list<string> &probeSetNames,
	const string &algName,
	const string &algVersion,
	const string &chipType,
	const string &programName,
	const string &programVersion,
	const string &programCompany,
	const vector<string>& paramNames,
	const vector<string>& paramValues,
	const vector<string>& sumNames,
	const vector<string>& sumValues
)
{
	// Determine the max probe set name.
	int numEntries = (int) probeSetNames.size();
	int maxProbeSetNameLength = 0;
	for (list<string>::const_iterator it=probeSetNames.begin(); it!=probeSetNames.end(); it++)
	{
		maxProbeSetNameLength = max(maxProbeSetNameLength, (int) it->length());
	}

	// Create the data object
	CHPQuantificationData *data = new CHPQuantificationData(outFile);
    data->SetEntryCount(numEntries, maxProbeSetNameLength);
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

	int nparams = (int) paramNames.size();
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
	CHPQuantificationFileWriter writer(*data);

	// Write the probe set names.
	writer.SeekToDataSet();
	ProbeSetQuantificationData entry;
	for (list<string>::const_iterator it=probeSetNames.begin(); it!=probeSetNames.end(); it++)
	{
        entry.name = *it;
        entry.quantification = 0.0f;
        writer.WriteEntry(entry);
    }
}

/*! Adds the expression (quantification) results to the CHP file. This function uses the buffer writer object
 * to write the data to the CHP file.
 * @param outFile The name of the output CHP file.
 * @param data The quantification (signal) results.
 */
static void UpdateQuantifications(const string &outFile, const list<float> &data)
{
	// Create the buffer writer object
	CHPQuantificationFileBufferWriter expressionQuantificationBufferWriter;
	vector<string> fileNames;
	fileNames.push_back(outFile);
	expressionQuantificationBufferWriter.Initialize(&fileNames);

	// Write the data
	for (list<float>::const_iterator it=data.begin(); it!=data.end(); it++)
	{
		expressionQuantificationBufferWriter.WriteQuantificationEntry(0, *it);
	}

	// Flush the buffer. The destructor will close the file.
	expressionQuantificationBufferWriter.FlushBuffer();
}

/*! Parse a line into a list.
 * @param line The line of data
 * @return A vector of strings from the line (separated by tabs).
 */
static vector<string> ParseLine(char *line)
{
	vector<string> tokens;
	const char *sep = "\t";
	char *token = strtok(line, sep);
	while (token != NULL)
	{
		tokens.push_back(token);
		token = strtok(NULL, sep);
	}
	return tokens;	
}

/*! Reads the expression results from a simple TSV text file.
 * @param inFile The name of the input TSV file.
 * @param names The probe set names from the input file.
 * @param quantifications The quantification values from the input file.
 */
static void ReadData(const string &inFile, list<string> &names, list<float> &quantifications)
{
	// Open the file
	ifstream instr(inFile.c_str(), ios::in);
	if (!instr)
		throw string("Unable to read the input file");

	// First line with header
	char line[1024];
	instr.getline(line, 1024);

	// Get the data
	names.clear();
	quantifications.clear();
	while (instr.getline(line, 1024))
	{
		vector<string> cols = ParseLine(line);
		names.push_back(cols[0]);
		quantifications.push_back((float)atof(cols[1].c_str()));
	}
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
	cout << tab << "-paramNames <a list of parameter names>" << endl;
	cout << tab << "-paramValues <a list of parameter values>" << endl;
	cout << tab << "-sumNames <a list of summary parameter names>" << endl;
	cout << tab << "-sumValues <a list of summary parameter values>" << endl;
	cout << endl;
	cout << "The input file must have the probe set name and quantification values" << endl <<
			"as the first two columns. The paramNames and paramValues arguments are" << endl <<
			"used to specify algorithm parameter." << endl;
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

		vector<string> paramNames;
		int n = cmdLine.GetArgumentCount("-paramNames");
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

		// Read the TSV input file
		list<string> names;
		list<float> quantifications;
		ReadData(inFile, names, quantifications);

		// Create the CHP file.
		CreateFileWithHeader(execId, celFile, outFile, names, algName, algVersion, arrayType,
			programName, programVersion, programCompany, paramNames, paramValues, sumNames, sumValues);
		UpdateQuantifications(outFile, quantifications);
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
