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
#include "chipstream/apt-sample-snp-list/SampleSNPListEngine.h"
//
#include "chipstream/EngineUtil.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/PgOptions.h"
//
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
//
using namespace std;
using namespace affx;

SampleSNPListEngine::Reg SampleSNPListEngine::reg;

SampleSNPListEngine * SampleSNPListEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == SampleSNPListEngine::EngineName())
		return (SampleSNPListEngine *)engine;
	return NULL;
}

SampleSNPListEngine::SampleSNPListEngine() : annotationFileHasOneCode(true)
{
    defineOptions();
    defineStates();
}

SampleSNPListEngine::~SampleSNPListEngine()
{
}

void SampleSNPListEngine::defineOptions()
{
	defineOption("", "count-file", PgOpt::STRING_OPT,
		"The count file.",
		"");
	defineOption("m", "marker-content-file", PgOpt::STRING_OPT,
		"The marker content file.",
        "");
	defineOption("", "annotation-file", PgOpt::STRING_OPT,
		"The annotation file.",
        "");
	defineOption("", "output-file", PgOpt::STRING_OPT,
		"The output file.",
		"");
	defineOption("", "min-particle-count", PgOpt::INT_OPT,
		"The minimum particle count.",
		"10");
	defineOption("", "note-col", PgOpt::STRING_OPT,
		"The column name for the note.",
        "Note");
	defineOption("", "probe-id-col", PgOpt::STRING_OPT,
		"The column name for the probe id.",
        "Probe Id");
	defineOption("", "type-col", PgOpt::STRING_OPT,
		"The column name for the type.",
        "Type");
	defineOption("", "probeset-id-col", PgOpt::STRING_OPT,
		"The column name for the probe set id.",
        "probeset_id");
}

void SampleSNPListEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void SampleSNPListEngine::checkOptionsImp()
{
	string val = getOpt("count-file");
	if (val.empty() == true)
		Err::errAbort("Must specify a count file.");
	
	val = getOpt("output-file");
	if (val.empty() == true)
		Err::errAbort("Must specify an output file.");

	val = getOpt("marker-content-file");
	if (val.empty() == true)
		Err::errAbort("Must specify a marker content file.");

	val = getOpt("annotation-file");
	if (val.empty() == true)
		Err::errAbort("Must specify an annotation file.");
}

void SampleSNPListEngine::writeHeader(std::ofstream &out)
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

/**
   This is the "main()" equivalent of the engine.
*/
void SampleSNPListEngine::runImp()
{
	// Check and get the options
	readMixFile();
	readAnnotationFile();
	string countFile = getOpt("count-file");
	string outFile = getOpt("output-file");
	int minCount = getOptInt("min-particle-count");

	// Open the input count file bind the data to the columns.
	Verbose::out(1, "Opening count file");
	affx::TsvFile tsv;
	if(tsv.open(countFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + countFile);
	tsv.rewind();
	string name;
	int nsamples = tsv.getColumnCount(0) - 1;
	vector<int> count(nsamples);
	tsv.bind(0, 0, &name, TSV_BIND_REQUIRED);
	for (int isample=0; isample<nsamples; isample++)
		tsv.bind(0, isample+1, &count[isample],  TSV_BIND_REQUIRED);

	// Create an output file and write the header
	ofstream out(outFile.c_str(), ios::out);
	writeHeader(out);
	out << getOpt("probeset-id-col");
	for (int isample=0; isample<nsamples; isample++)
		out << "\t" << tsv.getColumnName(0, isample+1);
	out << endl;

	// Read the counts and filter those that are below the min value.
	// Ignore the control probes.
	tsv.rewind();

	if (annotationFileHasOneCode)
	{
		// This assumes that the two allele codes are spread across two successive lines in the counts file.
		while(tsv.nextLevel(0) == TSV_OK)
		{
			if (ignoreProbes.find(name) != ignoreProbes.end())
				continue;

			vector<int> counta = count;
			tsv.nextLevel(0);
			out << name.substr(0, name.length() - 2);
			for (int isample=0; isample<nsamples; isample++)
				out << "\t" << (count[isample] < minCount || counta[isample] < minCount ? 0 : 2);
			out << endl;
		}
	}
	else
	{
		// This assumes that each line in the counts file has both allele codes.
		while(tsv.nextLevel(0) == TSV_OK)
		{
			if (ignoreProbes.find(name) != ignoreProbes.end())
				continue;

			out << name;
			for (int isample=0; isample<nsamples; isample++)
				out << "\t" << (count[isample] < minCount ? 0 : 2);
			out << endl;
		}
	}

	// Close the input and output files.
	out.close();
	tsv.close();
}

void SampleSNPListEngine::readMixFile()
{
	// Create a a mapping of the particle code to combination of probe id and type.
	// The probe will be stored in the first entry of the pair.
	string mixFile = getOpt("marker-content-file");
	TsvFile tsv;
	string id;
	string note;
	int typeValue;
	tsv.bind(0, getOpt("note-col"), &note, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("probe-id-col"), &id, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("type-col"), &typeValue, TSV_BIND_REQUIRED);
	if(tsv.open(mixFile) != TSV_OK)
		Err::errAbort("Couldn't open file: " + mixFile);
	tsv.rewind();
	while(tsv.nextLevel(0) == TSV_OK)
	{
		if (id.empty() == true)
			id = note;
		if (typeValue == -1)
			ignoreProbes[id] = true;
	}
	tsv.close();
}

void SampleSNPListEngine::readAnnotationFile()
{
	string file = getOpt("annotation-file");
	TsvFile tsv;

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
			annotationFileHasOneCode = (atoi(numChannels.c_str()) == 1);
	}

	tsv.close();
}
