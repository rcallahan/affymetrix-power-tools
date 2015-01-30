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

\page Txt2DmetChp_description Writing genotype results to a "multi-data" CHP file.

This program is provided as an example of using the Fusion SDK to write genotype results to an Affymetrix
"multi-data" CHP file. The multi-data CHP file was designed to store a variety of Affymetrix results,
including genotypes and copy number. The Genotyping Console software currently uses the "multi-data" CHP
file for storing genotype results from the BRLMM-P and Birdseed algorithms as well as copy number results
from its copy number algorithms.

The program is implemented as a command line application. The arguments specify items such as the TSV file
containing the genotype results, algorithm parameters, summary statistics, probe array type, parent CEL file, 
and other items that are stored in the CHP file header.

 */

#include "calvin_files/converters/chp/src/CmdLine.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileBufferWriter.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
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

static string biAlleleIn;
static string multiAlleleIn;
static string copyNumberIn;

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

static void WriteBiAllelicData(CHPMultiDataFileWriter &writer)
{
	if(!biAlleleIn.empty())
	{
		DmetBiAllelicData entry;
		writer.SeekToDataSet(DmetBiAllelicMultiDataType);

		// Open the file
		ifstream instr(biAlleleIn.c_str(), ios::in);
		if (!instr)
		{
			throw string("Unable to read the input file");
		}

		char line[1024];
		// First line with header
		instr.getline(line, 1024);
		while (instr.getline(line, 1024))
		{
			list<string> cols = ParseLine(line);
			if(cols.size() > 0)
			{
				list<string>::const_iterator colIt = cols.begin();
				entry.name = *colIt;
				++colIt;
				entry.call = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.confidence = (float)atof(colIt->c_str());
				++colIt;
				entry.force = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.signalA = (float)atof(colIt->c_str());;
				++colIt;
				entry.signalB = (float)atof(colIt->c_str());;
				++colIt;
				entry.contextA = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.contextB = (u_int8_t)atoi(colIt->c_str());

				/*++colIt;
				entry.metrics.resize(cols.size() - 3);
				int colIdx = 0;
				for (; colIt!=cols.end(); colIt++, colIdx++)
				{
				if (extraColTypes[colIdx] == "float")
				{
				entry.metrics[colIdx].SetValueFloat((float)atof(colIt->c_str()));
				}
				else if (extraColTypes[colIdx] == "int")
				{
				entry.metrics[colIdx].SetValueInt32(atoi(colIt->c_str()));
				}
				else if (extraColTypes[colIdx] == "ubyte")
				{
				entry.metrics[colIdx].SetValueUInt8(CallToCode(colIt->c_str()));
				}
				}*/

				writer.WriteEntry(entry);
			}
			else
			{
				break;
			}
		}
	}
}

static void WriteMultiAllelicData(CHPMultiDataFileWriter &writer)
{
	if(!multiAlleleIn.empty())
	{
		DmetMultiAllelicData entry;
		writer.SeekToDataSet(DmetMultiAllelicMultiDataType);

		// Open the file
		ifstream instr(multiAlleleIn.c_str(), ios::in);
		if (!instr)
		{
			throw string("Unable to read the input file");
		}

		char line[1024];
		// First line with header
		instr.getline(line, 1024);
		while (instr.getline(line, 1024))
		{
			list<string> cols = ParseLine(line);
			if(cols.size() > 0)
			{
				list<string>::const_iterator colIt = cols.begin();
				entry.name = *colIt;
				++colIt;
				entry.call = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.confidence = (float)atof(colIt->c_str());
				++colIt;
				entry.force = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.alleleCount = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.signalA = (float)atof(colIt->c_str());;
				++colIt;
				entry.signalB = (float)atof(colIt->c_str());;
				++colIt;
				entry.signalC = (float)atof(colIt->c_str());;
				++colIt;
				entry.signalD = (float)atof(colIt->c_str());;
				++colIt;
				entry.signalE = (float)atof(colIt->c_str());;
				++colIt;
				entry.signalF = (float)atof(colIt->c_str());;
				++colIt;
				entry.contextA = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.contextB = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.contextC = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.contextD = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.contextE = (u_int8_t)atoi(colIt->c_str());
				++colIt;
				entry.contextF = (u_int8_t)atoi(colIt->c_str());

				/*++colIt;
				entry.metrics.resize(cols.size() - 3);
				int colIdx = 0;
				for (; colIt!=cols.end(); colIt++, colIdx++)
				{
				if (extraColTypes[colIdx] == "float")
				{
				entry.metrics[colIdx].SetValueFloat((float)atof(colIt->c_str()));
				}
				else if (extraColTypes[colIdx] == "int")
				{
				entry.metrics[colIdx].SetValueInt32(atoi(colIt->c_str()));
				}
				else if (extraColTypes[colIdx] == "ubyte")
				{
				entry.metrics[colIdx].SetValueUInt8(CallToCode(colIt->c_str()));
				}
				}*/

				writer.WriteEntry(entry);
			}
			else
			{
				break;
			}
		}
	}
}

static void WriteCopyNumberData(CHPMultiDataFileWriter &writer)
{
	if(!copyNumberIn.empty())
	{
		DmetCopyNumberData entry;
		writer.SeekToDataSet(DmetCopyNumberMultiDataType);
	
		// Open the file
		ifstream instr(copyNumberIn.c_str(), ios::in);
		if (!instr)
		{
			throw string("Unable to read the input file");
		}

		char line[1024];
		// First line with header
		instr.getline(line, 1024);
		while (instr.getline(line, 1024))
		{
			list<string> cols = ParseLine(line);
			if(cols.size() > 0)
			{
				list<string>::const_iterator colIt = cols.begin();

				entry.name = *colIt;
				++colIt;
				entry.call = (int16_t)atoi(colIt->c_str());
				++colIt;
				entry.confidence = (float)atof(colIt->c_str());
				++colIt;
				entry.force = (int16_t)atoi(colIt->c_str());
				++colIt;
				entry.estimate = (float)atof(colIt->c_str());
				++colIt;
				entry.lower = (float)atof(colIt->c_str());
				++colIt;
				entry.upper = (float)atof(colIt->c_str());

				/*++colIt;
				entry.metrics.resize(cols.size() - 3);
				int colIdx = 0;
				for (; colIt!=cols.end(); colIt++, colIdx++)
				{
				if (extraColTypes[colIdx] == "float")
				{
				entry.metrics[colIdx].SetValueFloat((float)atof(colIt->c_str()));
				}
				else if (extraColTypes[colIdx] == "int")
				{
				entry.metrics[colIdx].SetValueInt32(atoi(colIt->c_str()));
				}
				else if (extraColTypes[colIdx] == "ubyte")
				{
				entry.metrics[colIdx].SetValueUInt8(CallToCode(colIt->c_str()));
				}
				}*/

				writer.WriteEntry(entry);
			}
			else
			{
				break;
			}
		}
	}
}

static string InsertNum(string file, int n)
{
	size_t idx = file.find_last_of(".");
	if(idx == 0)
	{
		throw string("Improper outfile name:  " + file);
	}
	string result;
	stringstream sstrm;
	sstrm << n;
	string num = sstrm.str();
	if(idx > 0)
	{
		result = file.substr(0, idx);
		result = result.append("_");
		result = result.append(num);
		result = result.append(file.substr(idx, file.length() - idx));
	}
	else
	{
		result = file.append("_");
		result = result.append(num);
	}
	return result;
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
	{
		throw string("Unable to read the input file");
	}
	
	char line[1024];
	// First line with header
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
		if(cols.size() > 0)
		{
			maxProbeSetNameLength = max(maxProbeSetNameLength, (int) cols.begin()->size());
			++count;
		}
		else
		{
			break;
		}
	}
	return count;
}

/*! Show a help message on the use of the command line parameters. */
static void Help()
{
	const char *tab = "\t";
	cout << "The command line arguments include" << endl;
	cout << tab << "-copies <number of copies to generate - default is 1>" << endl;
	cout << tab << "-execid <the execution id>" << endl;
	cout << tab << "-cel <the name of the parent cel file>" << endl;
	cout << tab << "-biAllele <the name of the biallelic input data file, TSV with 1 header row>" << endl;
	cout << tab << "-multiAllele <the name of the multiallelic input data file, TSV with 1 header row>" << endl;
	cout << tab << "-copyNumber <the name of the copy number input data file, TSV with 1 header row>" << endl;
	cout << tab << "-out <the name of the output CHP file>" << endl;
	cout << tab << "-arrayType <the probe array type>" << endl;
	cout << tab << "-algName <the algorithm name>" << endl;
	cout << tab << "-algVersion <the algorithm version>" << endl;
	cout << tab << "-programName <the program name>" << endl;
	cout << tab << "-programVersion <the program version>" << endl;
	cout << tab << "-programCompany <the name of the company>" << endl;
	cout << tab << "-paramNames <a list of parameter names>" << endl;
	cout << tab << "-paramTypes <a list of parameter types> 1 = ascii, 2 = integer, 3 = float - default is ascii" << endl;
	cout << tab << "-paramValues <a list of parameter values>" << endl;
	cout << tab << "-sumNames <a list of summary parameter names>" << endl;
	cout << tab << "-sumTypes <a list of summary parameter types> 1 = ascii, 2 = integer, 3 = float - default is ascii" << endl;
	cout << tab << "-sumValues <a list of summary parameter values>" << endl;
	cout << tab << "-extraNames <a list of extra parameter names>" << endl;
	cout << tab << "-extraTypes <a list of extra parameter types> 1 = ascii, 2 = integer, 3 = float - default is ascii" << endl;
	cout << tab << "-extraValues <a list of extra parameter values>" << endl;
	cout << tab << "-delay <delay exit app exit in milliseconds - optional>" << endl;
	cout << endl;
	cout << "DmetMultiAllelic CHP:" << endl <<
			"The input file must have the probe set name, call, confidence, force, allele count," << endl <<
			" signal A, signal B, signal C, signal D, signal E, signal F," << endl <<
			" context A, context B, context C, context D, context E, context F as the 17 columns. " << endl << endl <<
			"DmetBiAllelic CHP:" << endl <<
			"The input file must have the probe set name, call, confidence, force, signal A, signal B," << endl <<
			" context A, context B as the 8 columns. " << endl << endl <<
			"DmetCopyNumber CHP:" << endl <<
			"The input file must have the probe set name, call, confidence, force, estimate, lower, upper as the 7 columns. ";
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
			Help();
			exit(-1);
		}

		int copies = atoi(cmdLine.GetSafeArgument("-copies", 0, "1").c_str());
		int delay = atoi(cmdLine.GetSafeArgument("-delay", 0, "0").c_str());
		string outFile = cmdLine.GetArgument("-out", 0);
		CHPMultiDataData data(outFile);
		string s = cmdLine.GetSafeArgument("-cel", 0, "");
		if (s.length() > 0 && FileUtils::Exists(s.c_str()))
		{
			FusionCELData cel;
			cel.SetFileName(s.c_str());
			cel.ReadHeader();
			GenericData *gdata = cel.GetGenericData();
			if (gdata != NULL)
			{
				data.GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr()); 
			}
			cel.Close();
		}
		biAlleleIn = cmdLine.GetSafeArgument("-biAllele", 0, "");
		if(!biAlleleIn.empty())
		{
			int maxProbeSetNameLength = 0;
			unsigned long rowCount = ReadData(biAlleleIn, maxProbeSetNameLength);
			data.SetEntryCount(DmetBiAllelicMultiDataType, rowCount, maxProbeSetNameLength);
		}
		multiAlleleIn = cmdLine.GetSafeArgument("-multiAllele", 0, "");
		if(!multiAlleleIn.empty())
		{
			int maxProbeSetNameLength = 0;
			unsigned long rowCount = ReadData(multiAlleleIn, maxProbeSetNameLength);
			data.SetEntryCount(DmetMultiAllelicMultiDataType, rowCount, maxProbeSetNameLength);
		}
		copyNumberIn = cmdLine.GetSafeArgument("-copyNumber", 0, "");
		if(!copyNumberIn.empty())
		{
			int maxProbeSetNameLength = 0;
			unsigned long rowCount = ReadData(copyNumberIn, maxProbeSetNameLength);
			data.SetEntryCount(DmetCopyNumberMultiDataType, rowCount, maxProbeSetNameLength);
		}
		s = cmdLine.GetArgument("-arrayType", 0);
		data.SetArrayType(StringUtils::ConvertMBSToWCS(s));
		s = cmdLine.GetArgument("-algName", 0);
		data.SetAlgName(StringUtils::ConvertMBSToWCS(s));
		s = cmdLine.GetArgument("-algVersion", 0);
		data.SetAlgVersion(StringUtils::ConvertMBSToWCS(s));
		s = cmdLine.GetSafeArgument("-programName", 0, "");
		if(!s.empty())
		{
			ParameterNameValueType param;
			param.SetName(L"program-name");
			param.SetValueText(StringUtils::ConvertMBSToWCS(s));
			data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
		}
		s = cmdLine.GetSafeArgument("-programVersion", 0, "");
		if(!s.empty())
		{
			ParameterNameValueType param;
			param.SetName(L"program-version");
			param.SetValueText(StringUtils::ConvertMBSToWCS(s));
			data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
		}
		s = cmdLine.GetSafeArgument("-programCompany", 0, "");
		if(!s.empty())
		{
			ParameterNameValueType param;
			param.SetName(L"program-company");
			param.SetValueText(StringUtils::ConvertMBSToWCS(s));
			data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
		}
	
		s = cmdLine.GetArgument("-execid", 0);
		ParameterNameValueTypeList params;
		ParameterNameValueType param;
		param.SetName(L"exec-guid");
		param.SetValueAscii(s);
		params.push_back(param);
		int n = cmdLine.GetArgumentCount("-paramNames");
		for(int i = 0; i < n; i++)
		{
			string paramName = cmdLine.GetArgument("-paramNames", i);
			string paramType = cmdLine.GetSafeArgument("-paramTypes", i, "1");
			string paramValue = cmdLine.GetArgument("-paramValues", i);
			param.SetName(StringUtils::ConvertMBSToWCS(paramName));
			param.SetValueAscii(paramValue);
			if(paramType == "2")
			{
				param.SetValueInt32(atoi(paramValue.c_str()));
			}
			else if(paramType == "3")
			{
				param.SetValueFloat((float)atof(paramValue.c_str()));
			}
			else
			{
				param.SetValueAscii(paramValue);
			}
			params.push_back(param);
		}
		data.AddAlgParams(params);
		params.clear();
		n = cmdLine.GetArgumentCount("-sumNames");
		for (int i = 0; i < n; i++)
		{
			string sumName = cmdLine.GetArgument("-sumNames", i);
			string sumType = cmdLine.GetSafeArgument("-sumTypes", i, "1");
			string sumValue = cmdLine.GetArgument("-sumValues", i);
			param.SetName(StringUtils::ConvertMBSToWCS(sumName));
			if(sumType == "2")
			{
				param.SetValueInt32(atoi(sumValue.c_str()));
			}
			else if(sumType == "3")
			{
				param.SetValueFloat((float)atof(sumValue.c_str()));
			}
			else
			{
				param.SetValueAscii(sumValue);
			}
			params.push_back(param);
		}
		if (!params.empty())
		{
			data.AddSummaryParams(params);
			params.clear();
		}
		n = cmdLine.GetArgumentCount("-extraNames");
		for (int i = 0; i < n; i++)
		{
			string extraName = cmdLine.GetArgument("-extraNames", i);
			string extraType = cmdLine.GetSafeArgument("-extraTypes", i, "1");
			string extraValue = cmdLine.GetArgument("-extraValues", i);
			ParameterNameValueType param;
			param.SetName(StringUtils::ConvertMBSToWCS(extraName));
			if(extraType == "2")
			{
				param.SetValueInt32(atoi(extraValue.c_str()));
			}
			else if(extraType == "3")
			{
				param.SetValueFloat((float)atof(extraValue.c_str()));
			}
			else
			{
				param.SetValueAscii(extraValue);
			}
			data.GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
		}	

		for(int i = 0; i < copies; i++)
		{
			if(i > 0)
			{
				string s = InsertNum(outFile, i);
				data.SetFilename(s);
				AffymetrixGuidType guid = AffymetrixGuid::GenerateNewGuid();
				data.GetGenericData().Header().GetGenericDataHdr()->SetFileId(guid);
			}
			CHPMultiDataFileWriter writer(data);
			WriteBiAllelicData(writer);
			WriteMultiAllelicData(writer);
			WriteCopyNumberData(writer);
		}
		if(delay > 0)
		{
			Sleep(delay);
		}
	}
	catch (string s)
	{
		cout << s << endl;
		Help();
	}
	catch (int e)
	{
		cout << "Invalid argument" << endl;
		Help();
	}
	catch (...)
	{
		cout << "Unknown error" << endl;
		Help();
	}
	return 0;
}
