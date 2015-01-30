////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#include "chipstream/EngineUtil.h"
#include "chipstream/apt-list-summary/ListSummaryEngine.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Fs.h"
#include "util/PgOptions.h"
//
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//

using namespace std;
using namespace affx;

static double median(const vector<double> & u)
{
	int len = (int) u.size();
	long half = len / 2;
	return (len % 2 == 1 ? u[half] : (u[half - 1] + u[half]) / 2.0);
}

static double PercentileCV(const vector<double> &data, int innerPercentile, double median)
{
    int startIndex = (int)(data.size() * (100 - innerPercentile) / 200.0);
    int endIndex = data.size() - startIndex;
    double stdv = 0.0;
    double diff = 0.0;
    double sumOfSquares = 0.0;
    double avg = 0.0f;
    for (int valIndex = startIndex; valIndex < endIndex; valIndex++)
    {
        avg += data[valIndex];
    }
    avg /= (endIndex - startIndex);
    for (int valIndex = startIndex; valIndex < endIndex; valIndex++)
    {
        diff = (avg - data[valIndex]);
        sumOfSquares += diff * diff;
    }
    if (sumOfSquares != 0 && (endIndex - startIndex) != 0) {
        stdv = std::sqrt(sumOfSquares / (endIndex - startIndex - 1));
    }
    return stdv / avg;
}

static string GetFileTitle(const string &file)
{
	string name;
	string::size_type start = file.rfind("\\");
	if (start == string::npos)
		start = file.rfind("/");
	if (start == string::npos)
		start = 0;
	else
		++start;
	string::size_type end = file.rfind(".");
	if (end == string::npos)
		end = file.length();
	name = file.substr(start, end-start);
    return name;
}

static vector<string> replaceTag(string replacementString, string tag, vector<string> replacementValues)
{
	vector<string> results;

	string::size_type start = replacementString.find(tag);

	// The tag must be in the replacementString if there is more than one replacementValue to ensure that the results are
	// unique.
	if (replacementValues.size() > 1 && start == string::npos)
	{
		Err::errAbort("Input parameter is missing a required tag. '" + replacementString + "' is missing: '" + tag + "'");
	}
	
	if (start != string::npos)
	{
		// Found a replacement tag.
		for (vector<string>::iterator it = replacementValues.begin(); it != replacementValues.end(); ++it)
		{
			string result = replacementString;
			result.replace(start, tag.length(), *it);
			results.push_back(result);
		}
	}
	else
	{
		// Return the string unchanged.
		results.push_back(replacementString);
	}
	return results;
}

static void trimSpaces(string& str)
{
	// Trim Both leading and trailing spaces
	size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
	size_t endpos = str.find_last_not_of(" \t");		// Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(( string::npos == startpos ) || ( string::npos == endpos))
	{
		str = "";
	}
	else
		str = str.substr( startpos, endpos-startpos+1 );
}

static vector<string> split(const string& inputString, const string& delim)
{
	vector<string> tokens;
	size_t substrBegin = 0;
	for (;;)
	{
		size_t substrEnd = inputString.find (delim, substrBegin);
		if (substrEnd == string::npos)
		{
			// No more delim - save what's left, quit.
			string subString = inputString.substr (substrBegin);
			trimSpaces(subString);
			// Avoid returning a null string from a terminating '?' or an empty inputString.
			if (! subString.empty())
				tokens.push_back (subString);
			break;
		}
		// Avoid null strings from an initial 'delim' or 'delimdelim'.
		if (substrEnd != substrBegin)
		{
			string token(inputString.substr(substrBegin, substrEnd - substrBegin));
			trimSpaces(token);
			tokens.push_back (token);
		}
		// Continue following the delim
		substrBegin = substrEnd + delim.size();
	}
	return tokens;
}

static vector<string> split(const string &inputString)
{
	return split(inputString, "?");
}

static string parseAlleleCode(string channel, string alleleCode)
{
	string delim = "//";

	vector<string> channels = split(channel, delim);
	vector<string> alleleCodes = split(alleleCode, delim);

	if (channels.size() != alleleCodes.size())
	{
		Err::errAbort("Channel-Allele Code mismatch - Order: " + channel + ", Allele Code:" + alleleCode);
	}
	
	string result;
	for (int i = 0; i < channels.size(); ++i)

	{
		int pos = atoi(channels[i].c_str());
		if (result.size() < pos)
			result.resize(pos);
		result[pos-1] = alleleCodes[i][0];
	}

	return result;
}

static int decodeChannel(const string& signalColumn)
{
	string::size_type pos = signalColumn.find_first_of("0123456789");
	if (pos == string::npos)
	{
		Err::errAbort("Unable to deduce channel number: " + signalColumn);
	}

	return atoi(signalColumn.substr(pos, 1).c_str());
}

ListSummaryEngine::Reg ListSummaryEngine::reg;

ListSummaryEngine * ListSummaryEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == ListSummaryEngine::EngineName())
		return (ListSummaryEngine *)engine;
	return NULL;
}

ListSummaryEngine::ListSummaryEngine() : annotationFileHasOneCode(true) {
    defineOptions();
}

ListSummaryEngine::~ListSummaryEngine() {
}

void ListSummaryEngine::defineOptions() {
	defineOption("", "list-files", PgOpt::STRING_OPT,
		"Text file specifying list files to process, one per line with the first line being 'list_files'.",
		"");
	defineOption("p", "algo-params", PgOpt::STRING_OPT,
		"A string specifying algo parameters. See --explain option for acceptable parameters.",
		"");
	defineOption("", "summary-out-file", PgOpt::STRING_OPT,
		"File to write result (summaries) files into.",
		"");
	defineOption("", "cv-out-file", PgOpt::STRING_OPT,
		"File to write result (cv) files into.",
		"");
	defineOption("", "count-out-file", PgOpt::STRING_OPT,
		"File to write result (counts) files into.",
		"");
	defOptMult("l", "lists", PgOpt::STRING_OPT,
		"List files to process.",
        "");
	defineOption("m", "marker-content-file", PgOpt::STRING_OPT,
		"The marker content file.",
        "");
	defineOption("", "annotation-file", PgOpt::STRING_OPT,
		"The annotation file.",
        "");
	defineOption("", "files-col", PgOpt::STRING_OPT,
		"The column name for the input files.",
        "list_files");
	defineOption("", "code-col", PgOpt::STRING_OPT,
		"The column name for the particle code.",
        "PartCode");
	defineOption("", "note-col", PgOpt::STRING_OPT,
		"The column name for the note.",
        "Note");
	defineOption("", "probe-id-col", PgOpt::STRING_OPT,
		"The column name for the probe id.",
        "Probe Id");
	defineOption("", "type-col", PgOpt::STRING_OPT,
		"The column name for the type.",
        "Type");
	defineOption("", "signal-columns", PgOpt::STRING_OPT,
		"Particle intensity column names in the LIST files.  All files must have the same names.",
        "FluorCh2");
	defineOption("", "probeset-id-col", PgOpt::STRING_OPT,
		"The column name for the probe set id.",
        "Probeset Id");
	defineOption("", "channel-col", PgOpt::STRING_OPT,
		"The column name for the channel column.",
        "Channel");
	defineOption("", "allele-code-col", PgOpt::STRING_OPT,
		"The column name for the allele code column.",
        "Allele Code");
	defineOption("", "cv-percentile", PgOpt::INT_OPT,
		"The percentile range for the CV calculation.",
        "85");
	defineOption("", "include-all-codes", PgOpt::STRING_OPT,
		"Include codes not found in the mix or annotation files.",
        "no");
}

void ListSummaryEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void ListSummaryEngine::checkOptionsImp()
{
	defineStates();
	string file;
	file = getOpt("summary-out-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an output summary file.");
	file = getOpt("cv-out-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an output cv file.");
	file = getOpt("count-out-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an output count file.");
	file = getOpt("marker-content-file");
	if (file.empty() == true)
		Err::errAbort("Must specify a marker content file.");
	file = getOpt("annotation-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an annotation file.");

	// Read in list files list from other file if specified.
	vector<string> listFiles;
	if(getOpt("list-files")!="") {
		affx::TsvFile tsv;
#ifdef WIN32
		tsv.m_optEscapeOk = false;
#endif
		std::string listFilesFile = getOpt("list-files");
		tsv.bind(0, getOpt("files-col"), &file, TSV_BIND_REQUIRED);

		try {
			if(tsv.open(listFilesFile) != TSV_OK) {
				Err::errAbort("Couldn't open list-files file: " + listFilesFile);
			}
		} catch (Except&) {
			tsv.close();	// need to close the file
			throw;
		}

		tsv.rewind();
		while(tsv.nextLevel(0) == TSV_OK) {
			listFiles.push_back(file);
		}
		tsv.close();
		Verbose::out(1, "Read " + ToStr(listFiles.size()) + " list files from: " + Fs::basename(listFilesFile));
	}
	else
	{
		listFiles = getOptVector("lists");
	}
	setOpt("lists",listFiles);
}

void ListSummaryEngine::writeHeader(std::ofstream &out, const string &signalColumn, list<pair<string, float> > *unknownCodeRates)
{
	std::vector<std::string> option_names;
	getOptionNames(option_names);
	std::vector<PgOpt::PgOptType> option_types;
	getOptionTypes(option_types);
	for (int i = 0; i < option_names.size(); i++) 
	{
		switch (option_types[i]) 
		{
		case PgOpt::BOOL_OPT:
			out << "#%" << option_names[i] << "=" << getOptBool(option_names[i]) << endl;
			break;
		case PgOpt::DOUBLE_OPT:
			out << "#%" << option_names[i] << "=" << getOptDouble(option_names[i]) << endl;
			break;
		case PgOpt::INT_OPT:
			out << "#%" << option_names[i] << "=" << getOptInt(option_names[i]) << endl;
			break;
		case PgOpt::STRING_OPT:
			if (getOpt(option_names[i]).empty() == false)
				out << "#%" << option_names[i] << "=" << getOpt(option_names[i]) << endl;
		}
	}
	if (unknownCodeRates != NULL)
	{
		for (list<pair<string, float> >::iterator it=unknownCodeRates->begin(); it!=unknownCodeRates->end(); it++)
		{
			out << "#%Unknown Code Rate." << (*it).first << "=" << (*it).second << endl;
		}
	}
	if (signalColumn.empty() == false)
	{
		out << "#%signal-column=" << signalColumn << endl;
	}
}

/**
   This is the "main()" equivalent of the engine.
*/
void ListSummaryEngine::runImp()
{
	markerOrder.clear();
	mixMap.clear();
	annotMap.clear();
	readMixFile();
	readAnnotationFile();

	// Read the input list files and compute the median signal for each code.
	vector<string> listFiles = getOptVector("lists");
	//map<string, map<string, pair<uint64_t, int> > > summaryData;
	// map <LIST file name, map<partid, pair<array of median-cv in signal col order, count> > > summaryData
	map<string, map<string, pair< vector<uint64_t>, int> > > summaryData;
	list<pair<string, float> > unknownCodeRates;
	float unknownCodeRate;
	for (vector<string>::iterator it=listFiles.begin(); it!=listFiles.end(); it++)
	{
		summaryData[*it] = computeMedianData(*it, unknownCodeRate);
		unknownCodeRates.push_back(std::make_pair<string, float>(GetFileTitle(*it), unknownCodeRate));
	}

	// Output the data to the summary and CV files.
	vector<string> signalColumns = split(getOpt("signal-columns"));
	vector<string> summaryFiles = replaceTag(getOpt("summary-out-file"), "[SIGNAL_COLUMNS]", signalColumns);
	vector<string> cvFiles = replaceTag(getOpt("cv-out-file"), "[SIGNAL_COLUMNS]", signalColumns);

	for (int isignal = 0; isignal < summaryFiles.size(); ++isignal)
	{
		int channel = isignal + 1; //decodeChannel(signalColumns[isignal]);
		ofstream outSum(summaryFiles[isignal].c_str(), ios::out);
		writeHeader(outSum, signalColumns[isignal]);
		ofstream outCV(cvFiles[isignal].c_str(), ios::out);
		writeHeader(outCV, signalColumns[isignal]);

		// Output the column headers.
		outSum << "probeset_id";
		outCV << "probeset_id";
		for (vector<string>::iterator it=listFiles.begin(); it!=listFiles.end(); it++)
		{
			const string& listFile = *it;
			string title = GetFileTitle(listFile);
			outSum << "\t" << title;
			outCV << "\t" << title;
		}
		outSum << endl;
		outCV << endl;

		// Loop over each marker/allele and output its summary data.
		for (list<string>::iterator orderIt=markerOrder.begin(); orderIt!=markerOrder.end(); orderIt++)
		{
			const string &id = *orderIt;

			if (annotMap.find(id) != annotMap.end())
			{
				const pair<string, string> &annot = annotMap[id];
				outSum << annot.first << "-" << annot.second[channel-1];
				outCV << annot.first << "-" << annot.second[channel-1];
			}
			else
			{
				outSum << id;
				outCV << id;
			}

			for (vector<string>::iterator it=listFiles.begin(); it!=listFiles.end(); it++)
			{
				const string& listFile = *it;
				outSum << "\t";
				outCV << "\t";
				if (summaryData[listFile].find(id) != summaryData[listFile].end())
				{
					uint64_t u64 = summaryData[listFile][id].first[isignal];
					uint32_t um = (uint32_t)(u64 >> 32);
					uint32_t ucv = (uint32_t)u64;
					outSum << *(float*)&um;
					outCV << *(float*)&ucv;
				}
				else
				{
					outSum << -1; //-std::numeric_limits<float>::infinity();
					outCV << -1; //-std::numeric_limits<float>::infinity();
				}
			}
			outSum << endl;
			outCV << endl;
		}

		// Close the files.
		outSum.close();
		outCV.close();
	}

	// Write the counts file.
	string countFile = getOpt("count-out-file");
	ofstream outCount(countFile.c_str(), ios::out);
	writeHeader(outCount, "", &unknownCodeRates);
	outCount << "probeset_id";
	for (vector<string>::iterator it=listFiles.begin(); it!=listFiles.end(); it++)
	{
		const string& listFile = *it;
		string title = GetFileTitle(listFile);
		outCount << "\t" << title;
	}

	outCount << endl;

	int oneChannel = -1;
	bool addAlleleCodeToCountsProbesetId = (annotationFileHasOneCode && signalColumns.size() == 1);
	if (addAlleleCodeToCountsProbesetId)
	{
		oneChannel = 1; //decodeChannel(signalColumns[0]);
	}

	// Loop over each marker/allele and output its summary data.
	for (list<string>::iterator orderIt=markerOrder.begin(); orderIt!=markerOrder.end(); orderIt++)
	{
		const string &id = *orderIt;
		if (annotMap.find(id) != annotMap.end())
		{
			if (addAlleleCodeToCountsProbesetId)
				outCount << annotMap[id].first << "-" << annotMap[id].second[oneChannel-1];
			else
				outCount << annotMap[id].first;
		}
		else
		{
			outCount << id;
		}

		for (vector<string>::iterator it=listFiles.begin(); it!=listFiles.end(); it++)
		{
			const string& listFile = *it;
			outCount << "\t";
			if (summaryData[listFile].find(id) != summaryData[listFile].end())
			{
				outCount << summaryData[listFile][id].second;
			}
			else
			{
				outCount << 0;
			}
		}
		outCount << endl;
	}
	outCount.close();
}

map<string, pair< vector<uint64_t>, int> > ListSummaryEngine::computeMedianData(const string& listFile, float &unknownCodeRate)
{
	int cvPercentile = getOptInt("cv-percentile");
	string includeAllCodesString = getOpt("include-all-codes");
	bool includeAllCodes = includeAllCodesString.length() > 0 && (includeAllCodesString[0] == 'y' || includeAllCodesString[0] == 'Y');
	map<string, vector<vector<double> > > data;
	string code;
	string id;
	string str_signal;
	double signal;
	unknownCodeRate = 0.0f;

	// Open the file and bind to the name and data columns.
	affx::TsvFile tsv;
	tsv.bind(0, getOpt("code-col"), &code, TSV_BIND_REQUIRED);
	//tsv.bind(0, getOpt("signal-col"), &str_signal, TSV_BIND_REQUIRED);

	vector<string> signalColumns = split(getOpt("signal-columns"));
	vector<string> str_signals;
	str_signals.resize(signalColumns.size());

	for (int ibind = 0; ibind < signalColumns.size(); ++ibind)
	{
		tsv.bind(0, signalColumns[ibind], &str_signals[ibind], TSV_BIND_REQUIRED);
	}

	try {
		if(tsv.open(listFile) != TSV_OK)
			Err::errAbort("Couldn't open the file: " + listFile);
	} catch (Except&) {
		tsv.close();	// need to close the file
		throw;
	}

	tsv.rewind();

	// Store the signal values for each marker.
	int totalParticles = 0;
	int totalFoundParticles = 0;
	while(tsv.nextLevel(0) == TSV_OK)
	{
		++totalParticles;
		if (mixMap.find(code) == mixMap.end())
		{
			if (includeAllCodes)
			{
				// Include all codes even if they are not in the mix or annotation files.
				id = code;
				if (std::find(markerOrder.begin(), markerOrder.end(), code) == markerOrder.end())
				{
					markerOrder.push_back(code);
				}
			}
			else
			{
				continue;
			}
		}
		else
		{
			id = mixMap[code].first;
			++totalFoundParticles;
		}

		// Get the id (from the mix file) given the code name and find it in the annotation file.
		//bool found = (annotMap.find(id) != annotMap.end());

		// If the marker is not in the annotation file and the type is > 0 then ignore this.
		// This is a way to fail markers in silico.
		//if (mixMap[code].second > 0 && found == false)
		//	continue;

		// If the type is > 0 then get the id from the annotation mapping.
		// This should be the allele name.
		//if (mixMap[code].second > 0)
		//	id = annotMap[id];

		// If the probe id from the .mix file is found in the annotation file, then
		// change the id to the annotation file code.
		//if (found == true)
		//{
		//	id = annotMap[id].first;
		//}

		for (int icol = 0; icol < signalColumns.size(); ++icol)
		{
			if (str_signals[icol] == "NaN") {
				// throw exception per Carsten 12-14-2009
				//Err::errAbort("NaN signal value found: " + listFile);
				continue;
			}

			// Add an empty data vector to the map if not found.
			if (data.find(id) == data.end())
			{
				vector<vector<double> > empty;
				empty.resize(signalColumns.size());
				data[id] = empty;
			}

			// Add the signal to the vector.
			signal = atof(str_signals[icol].c_str());
			data[id][icol].push_back(signal);
		}
	}
	tsv.close();

	// Compute the percentage of particles in the input file not found in the mix file.
	unknownCodeRate = 100.0f * (totalParticles - totalFoundParticles) / totalParticles;

	// Compute the median for each marker. If the number of particles
	// is less than a threshold then set the values to 0.
	map<string, pair< vector<uint64_t>, int> > medianData;
	for (map<string, vector<vector<double> > >::iterator it=data.begin(); it!=data.end(); it++)
	{
		vector<uint64_t> result;
		for (int icol = 0; icol < signalColumns.size(); ++icol)
		{
			std::sort(it->second[icol].begin(), it->second[icol].end());
			float m = median(it->second[icol]);
			float cv = PercentileCV(it->second[icol], cvPercentile, m);
			int n = (int) it->second[icol].size();
			uint32_t *um = (uint32_t *)&m;
			uint32_t *ucv = (uint32_t *)&cv;
			uint64_t u64 = ((uint64_t)*um << 32 | *ucv);
			result.push_back(u64);
		}

		medianData[it->first] = std::make_pair(result, (int) it->second[0].size());	// all signal columns for a particle will have the same number of values.
	}
	return medianData;
}

void ListSummaryEngine::readAnnotationFile()
{
	// Create a mapping of tag id to probe set allele codes
	string file = getOpt("annotation-file");
	TsvFile tsv;
	string probeid;
	string probesetid;
	string channel;
	string code;
	tsv.bind(0, getOpt("probe-id-col"), &probeid, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("probeset-id-col"), &probesetid, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("channel-col"), &channel, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("allele-code-col"), &code, TSV_BIND_REQUIRED);

	try {
		if(tsv.open(file) != TSV_OK)
			Err::errAbort("Couldn't open file: " + file);
	} catch (Except&) {
		tsv.close();	// need to close the file
		throw;
	}

	// If the num-channels header parameter is missing assume this is an older annotation file with only 1 channel.
	annotationFileHasOneCode = true;

	string numChannels;
	tsv.headersBegin();
	if (tsv.headersFindNext("num-channels", numChannels) != TSV_HEADER_LAST)
	{
		if (numChannels.empty() == false)
			annotationFileHasOneCode = atoi(numChannels.c_str()) == 1;
	}

	//if(tsv.nextLevel(0) == TSV_OK)
	//{
	//	// If the annotation file knows about more than one channel, the counts files should not report
	//	// the allele code as part of the probeset id.
	//	vector<string> channels = split(channel, "//");
	//	annotationFileHasOneCode = (channels.size() == 1);
	//}
	//tsv.rewind();

	while(tsv.nextLevel(0) == TSV_OK)
	{
		markerOrder.push_back(probeid);
		annotMap[probeid] = make_pair<string, string>(probesetid, parseAlleleCode(channel, code));
	}
	tsv.close();
}

void ListSummaryEngine::readMixFile()
{
	// Create a a mapping of the particle code to combination of probe id and type.
	// The probe will be stored in the first entry of the pair.
	string mixFile = getOpt("marker-content-file");
	TsvFile tsv;
	string id;
	string code;
	string note;
	int typeValue;
	tsv.bind(0, getOpt("note-col"), &note, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("code-col"), &code, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("probe-id-col"), &id, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("type-col"), &typeValue, TSV_BIND_REQUIRED);

	try {
		if(tsv.open(mixFile) != TSV_OK)
			Err::errAbort("Couldn't open file: " + mixFile);
	} catch (Except&) {
		tsv.close();	// need to close the file
		throw;
	}

	tsv.rewind();
	while(tsv.nextLevel(0) == TSV_OK)
	{
		if (id.empty() == true)
		{
			id = note;
			markerOrder.push_back(id);
		}
		mixMap[code] = make_pair<string, int>(id, typeValue);
	}
	tsv.close();
}
