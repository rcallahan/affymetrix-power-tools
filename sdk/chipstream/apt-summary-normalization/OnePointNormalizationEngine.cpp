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

#include "chipstream/apt-summary-normalization/OnePointNormalizationEngine.h"
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
#include <string>
#include <vector>
//
using namespace std;
using namespace affx;

OnePointNormalizationEngine::Reg OnePointNormalizationEngine::reg;

OnePointNormalizationEngine * OnePointNormalizationEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == OnePointNormalizationEngine::EngineName())
		return (OnePointNormalizationEngine *)engine;
	return NULL;
}

OnePointNormalizationEngine::OnePointNormalizationEngine()
{
    defineOptions();
    defineStates();
}

OnePointNormalizationEngine::~OnePointNormalizationEngine()
{
}

void OnePointNormalizationEngine::defineOptions()
{
	defineOption("", "in-summary-file", PgOpt::STRING_OPT,
		"Text file specifying the summary data.",
		"");
	defineOption("", "in-count-file", PgOpt::STRING_OPT,
		"Text file specifying the count data.",
		"");
	defineOption("", "marker-content-file", PgOpt::STRING_OPT,
		"Text file specifying the mix file.",
		"");
	defineOption("", "out-summary-file", PgOpt::STRING_OPT,
		"Text file specifying the summary data.",
		"");
	defineOption("", "min-particle-count", PgOpt::INT_OPT,
		"The minimum particle count.",
		"10");
	defineOption("", "min-value", PgOpt::DOUBLE_OPT,
		"The minimum summary value.",
		"1.0");
	defineOption("", "type-col", PgOpt::STRING_OPT,
		"The column name for the type.",
        "Type");
	defineOption("", "note-col", PgOpt::STRING_OPT,
		"The column name for the note.",
        "Note");
	defineOption("", "probe-id-col", PgOpt::STRING_OPT,
		"The column name for the probe id.",
        "Probe Id");
}

void OnePointNormalizationEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void OnePointNormalizationEngine::checkOptionsImp()
{
	string file;
	file = getOpt("in-summary-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an input summary file.");
	file = getOpt("in-count-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an input count file.");
	file = getOpt("marker-content-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an input marker file.");
	file = getOpt("out-summary-file");
	if (file.empty() == true)
		Err::errAbort("Must specify an output summary file.");
}

void OnePointNormalizationEngine::writeHeader(std::ofstream &out)
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
void OnePointNormalizationEngine::runImp()
{
	string inFile = getOpt("in-summary-file");
	string countFile = getOpt("in-count-file");
	string mixFile = getOpt("marker-content-file");
	string outFile = getOpt("out-summary-file");
	double minValue = getOptDouble("min-value");
	int minCount = getOptInt("min-particle-count");

	// Open the mix file and store those negative type entries.
	// These will be excluded from the normalized summary file.
	affx::TsvFile tsv;
	if(tsv.open(mixFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + mixFile);
	tsv.rewind();
	map<string, bool> excludeMap;
	int typeValue;
	string note;
	string id;
	tsv.bind(0, getOpt("type-col"), &typeValue, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("note-col"), &note, TSV_BIND_REQUIRED);
	tsv.bind(0, getOpt("probe-id-col"), &id, TSV_BIND_REQUIRED);
	while(tsv.nextLevel(0) == TSV_OK)
	{
		if (typeValue < 0)
		{
			if (id.empty() == true)
				id = note;
			excludeMap[id] = true;
		}
	}
	tsv.clear();

	// Open the input summary and count files
	affx::TsvFile countTsv;
	if(tsv.open(inFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + inFile);
	tsv.rewind();
	if(countTsv.open(countFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + countFile);
	countTsv.rewind();

	// Get the column names and save to the output file.
	ofstream out(outFile.c_str(), ios::out);
	writeHeader(out);
	string name;
	int n = tsv.getColumnCount(0);
	for (int i=0; i<n; i++)
	{
		name = tsv.getColumnName(0, i);
		out << name;
		if (i < n - 1)
			out << "\t";
	}
	out << endl;

	// Allocate a vector of data for the summary and count data
	// for the row and bind to the columns.
	vector<double> data(n-1);
	vector<int> count(n-1);
	tsv.bind(0, 0, &name, TSV_BIND_REQUIRED);
	for (int i=1; i<n; i++)
	{
		tsv.bind(0, i, &data[i-1],  TSV_BIND_REQUIRED);
		countTsv.bind(0, i, &count[i-1],  TSV_BIND_REQUIRED);
	}

	// Get the data from the input cound and summary files.
	vector<string> names;
	vector<vector<float> > sampleData(n-1); 
	vector<vector<int> > countData(n-1); 
	while(tsv.nextLevel(0) == TSV_OK)
	{
		countTsv.nextLevel(0);
		names.push_back(name);
		for (int i=0; i<n-1; i++)
		{
			sampleData[i].push_back(data[i]);
			countData[i].push_back(count[i]);
		}
	}
	data.erase(data.begin(), data.end());
	count.erase(count.begin(), count.end());
	tsv.clear();
	countTsv.clear();

	// Compute the min value for each sample. If no data was found
	// for the marker (all the data was below the minimum count)
	// then set the minimum marker value to the input minimum value.
	vector<float> minValues(n-1);
	minValues.assign(n-1, 9999999999.0f);
	for (int i=0; i<n-1; i++)
	{
		bool found = false;
		for (int j=0; j<(int)sampleData[i].size(); j++)
		{
			if (countData[i][j] >= minCount)
			{
				found = true;
				minValues[i] = min(minValues[i], sampleData[i][j]);
			}
		}
		if (found == false)
			minValues[i] = -1.0f;
	}

	// Output the normalized signals. The normalized values are
	// the summary value - min value for the marker + min input value.
	// If no values were found for a sample then just output a zero value.
	int nrows = (int) names.size();
	for (int j=0; j<nrows; j++)
	{
		if (excludeMap.find(names[j]) != excludeMap.end())
			continue;
		out << names[j];
		for (int i=0; i<n-1; i++)
		{
			out << "\t";
			if (countData[i][j] > 0 && minValues[i] >= 0)
				out << max(sampleData[i][j] - minValues[i] + minValue, minValue);
			else
				out << -1; //-std::numeric_limits<float>::infinity();
		}
		out << endl;
	}

	// Close the output file
	out.close();
}

