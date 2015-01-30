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

#include "chipstream/apt-summary-normalization/TwoPointNormalizationEngine.h"
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

static float average(const vector<float> &data)
{
	double sum = 0;
	int n = (int)data.size();
	for (int i=0; i<n; i++)
		sum += data[i];
	return (float) (sum / n);
}

static float percentile(const vector<float> &data, int percentile)
{
	if (data.size() == 0)
		return numeric_limits<float>::quiet_NaN();
    double dindex = percentile * ((data.size()-1) / 100.0f);
    int index = (int) (dindex + 0.5f);
	if (abs(dindex - index) <= numeric_limits<float>::epsilon())
        return data[index];
    else
        return (data[index-1] + data[index])/2;
}

TwoPointNormalizationEngine::Reg TwoPointNormalizationEngine::reg;

TwoPointNormalizationEngine * TwoPointNormalizationEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == TwoPointNormalizationEngine::EngineName())
		return (TwoPointNormalizationEngine *)engine;
	return NULL;
}

TwoPointNormalizationEngine::TwoPointNormalizationEngine()
{
    defineOptions();
    defineStates();
}

TwoPointNormalizationEngine::~TwoPointNormalizationEngine()
{
}

void TwoPointNormalizationEngine::defineOptions()
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
	defineOption("", "percentile", PgOpt::INT_OPT,
		"The percentile value to calculate.",
		"50");
	defineOption("", "neg-control-target", PgOpt::INT_OPT,
		"The negative control target.",
		"575");
	defineOption("", "target", PgOpt::INT_OPT,
		"The non-negative control target.",
		"800");
	defineOption("", "type-col", PgOpt::STRING_OPT,
		"The column name for the type.",
        "Type");
	defineOption("", "note-col", PgOpt::STRING_OPT,
		"The column name for the note.",
        "Note");
	defineOption("", "probe-id-col", PgOpt::STRING_OPT,
		"The column name for the probe id.",
        "Probe Id");
	defineOption("", "min-particle-count", PgOpt::INT_OPT,
		"The minimum particle count.",
		"10");
}

void TwoPointNormalizationEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void TwoPointNormalizationEngine::checkOptionsImp()
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

void TwoPointNormalizationEngine::writeHeader(std::ofstream &out)
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
void TwoPointNormalizationEngine::runImp()
{
	string inFile = getOpt("in-summary-file");
	string countFile = getOpt("in-count-file");
	string mixFile = getOpt("marker-content-file");
	string outFile = getOpt("out-summary-file");
	int perc = getOptInt("percentile");
	int negTarget = getOptInt("neg-control-target");
	int target = getOptInt("target");
	int minCount = getOptInt("min-particle-count");

	// Open the mix file and store those negative type entries.
	affx::TsvFile tsv;
	if(tsv.open(mixFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + mixFile);
	tsv.rewind();
	map<string, bool> controlMap;
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
			controlMap[id] = true;
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

	// Allocate a vector of data for the row and bind to the columns.
	vector<double> data(n-1);
	vector<int> count(n-1);
	tsv.bind(0, 0, &name, TSV_BIND_REQUIRED);
	for (int i=1; i<n; i++)
	{
		tsv.bind(0, i, &data[i-1],  TSV_BIND_REQUIRED);
		countTsv.bind(0, i, &count[i-1],  TSV_BIND_REQUIRED);
	}

	// Get the data. Separate the control from non-control values.
	// Only retrieve the counts for the controls.
	vector<string> names;
	vector<vector<float> > sampleData(n-1); 
	vector<vector<float> > controlData(n-1); 
	vector<vector<int> > countData(n-1);  
	vector<vector<int> > nonControlCountData(n-1); 
	while(tsv.nextLevel(0) == TSV_OK)
	{
		countTsv.nextLevel(0);
		bool isControl = (controlMap.find(name) != controlMap.end());
		if (isControl == false)
			names.push_back(name);
		for (int i=0; i<n-1; i++)
		{
			if (isControl == true)
			{
				controlData[i].push_back(data[i]);
				countData[i].push_back(count[i]);
			}
			else
			{
				sampleData[i].push_back(data[i]);
				nonControlCountData[i].push_back(count[i]);
			}
		}
	}
	data.erase(data.begin(), data.end());
	count.erase(count.begin(), count.end());
	tsv.clear();
	countTsv.clear();

	// Compute the weighted average of the controls and the
	// percentile of the non-controls.
	vector<float> wavg(n-1);
	vector<float> percentiles(n-1);
	for (int i=0; i<n-1; i++)
	{
		vector<float> temp;
		for (int j=0; j<(int)nonControlCountData[i].size(); j++)
		{
			if (nonControlCountData[i][j] >= minCount)
				temp.push_back(sampleData[i][j]);
		}
		std::sort(temp.begin(), temp.end());
		percentiles[i] = percentile(temp, perc);
		wavg[i] = 0;
		int totalCount = 0;
		for (int j=0; j<(int)controlData[i].size(); j++)
		{
			totalCount += countData[i][j];
			wavg[i] += (controlData[i][j] * countData[i][j]);
		}
		if (totalCount > 0)
			wavg[i] /= totalCount;
	}
	float avgPercentiles = target; //average(percentiles);
	float avgControl = negTarget; //average(wavg);

	// Compute the slope and intercept for each sample.
	vector<float> slope(n-1);
	vector<float> intercept(n-1);
	for (int i=0; i<n-1; i++)
	{
		slope[i] = (avgPercentiles - avgControl) / (percentiles[i] - wavg[i]);
		intercept[i] = avgControl - (slope[i] * wavg[i]);
	}

	// Output the normalized signals
	// If there is not data then output a zero.
	int nrows = (int) names.size();
	for (int j=0; j<nrows; j++)
	{
		out << names[j];
		for (int i=0; i<n-1; i++)
		{
			out << "\t";
			if (sampleData[i].empty() == false)
			{
				if (nonControlCountData[i][j] > 0)
					out << (sampleData[i][j] * slope[i]) + intercept[i];
				else
					out << -1; //-std::numeric_limits<float>::infinity();
			}
			else
				out << -1; //-std::numeric_limits<float>::infinity();
		}
		out << endl;
	}

	// Close the files
	out.close();
}

