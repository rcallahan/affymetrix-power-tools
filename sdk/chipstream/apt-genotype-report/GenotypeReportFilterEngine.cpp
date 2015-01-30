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

//
#include "chipstream/apt-genotype-report/GenotypeReportFilterEngine.h"
//
#include "chipstream/EngineUtil.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/PgOptions.h"
//
#include <cstring>
#include <map>
#include <string>
//
using namespace std;
using namespace affx;

GenotypeReportFilterEngine::Reg GenotypeReportFilterEngine::reg;

GenotypeReportFilterEngine * GenotypeReportFilterEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == GenotypeReportFilterEngine::EngineName())
		return (GenotypeReportFilterEngine *)engine;
	return NULL;
}

GenotypeReportFilterEngine::GenotypeReportFilterEngine() {
    defineOptions();
}

GenotypeReportFilterEngine::~GenotypeReportFilterEngine() {
}

void GenotypeReportFilterEngine::defineOptions() {
	defineOption("s", "sample-file", PgOpt::STRING_OPT,
		"Text file specifying full paths to sample files.",
		"");
	defineOption("r", "report-file", PgOpt::STRING_OPT,
		"Text file specifying the report metrics.",
		"");
	defineOption("o", "out-file", PgOpt::STRING_OPT,
		"Output text file containing the files that pass the thresholds.",
		"");
	defineOption("", "file-col", PgOpt::STRING_OPT,
		"The sample file column header.",
		"list_files");
	defineOption("", "sample-status-col", PgOpt::STRING_OPT,
		"The sample status column header.",
		"Sample Status");
	defineOption("", "sample-col", PgOpt::STRING_OPT,
		"The sample column header.",
		"Sample");
}

void GenotypeReportFilterEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void GenotypeReportFilterEngine::checkOptionsImp()
{
    defineStates();
	string str;
	str = getOpt("sample-file");
	if (str.empty() == true)
		Err::errAbort("Must specify the input sample file.");

	str = getOpt("report-file");
	if (str.empty() == true)
		Err::errAbort("Must specify the input report file.");

	str = getOpt("out-file");
	if (str.empty() == true)
		Err::errAbort("Must specify an output file.");
}

void GenotypeReportFilterEngine::writeHeader(std::ofstream &out)
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

/**
   This is the "main()" equivalent of the engine.
*/
void GenotypeReportFilterEngine::runImp()
{
	// Open the input sample file and create a map of
	// file name to full path.
	affx::TsvFile tsv;
	string file = getOpt("sample-file");
#ifdef WIN32
	tsv.m_optEscapeOk = false;
#endif
	if(tsv.open(file) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + file);
	tsv.rewind();
	string sampleFile;
	map<string, string> sampleMap;
	tsv.bind(0, getOpt("file-col"), &sampleFile, TSV_BIND_REQUIRED);
	while(tsv.nextLevel(0) == TSV_OK)
	{
		sampleMap[GetFileTitle(sampleFile)] = sampleFile;
 	}
	tsv.clear();

	// Open the input report file. Each row represents a sample with
	// call rate. These will be compared against the threshold.
	file = getOpt("report-file");
	if(tsv.open(file) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + file);
	tsv.rewind();
	string sampleStatus;
	string sampleName;
	tsv.bind(0, getOpt("sample-col"), &sampleName, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("sample-status-col"), &sampleStatus, TSV_BIND_REQUIRED);

	// Create the output file
	file = getOpt("out-file");
	ofstream out(file.c_str(), ios::out);
	writeHeader(out);
	out << getOpt("file-col") << endl;

	// Determine those files that pass the thresholds and output them
	// to the output file.
	while(tsv.nextLevel(0) == TSV_OK)
	{
		if (sampleMap.find(sampleName) == sampleMap.end())
			Err::errAbort(sampleName + " was not found in the sample file.");
		if (sampleStatus == "Pass")
			out << sampleMap[sampleName] << endl;
	}

	// Close the files.
	out.close();
	tsv.close();
}
