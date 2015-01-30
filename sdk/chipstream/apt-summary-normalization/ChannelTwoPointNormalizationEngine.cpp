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

#include "chipstream/apt-summary-normalization/ChannelTwoPointNormalizationEngine.h"
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
#include <limits>

#ifdef WIN32
#ifdef max
#undef max
#endif
#endif

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

static string removeAlleleCode(const string &name)
{
	string::size_type n = name.find_last_of("-");
	if (n != string::npos)
	{
		return name.substr(0, n);
	}
	else
	{
		return name;
	}
}

static bool hasAlleleCode(const string &name)
{
	string::size_type n = name.find_last_of("-");
	return (n != string::npos);
}

// This maps between a particle name and the data for that particle.
// Used in sorting.
struct NameDataMap
{
	string name;	// particle name
	int ifile;		// channel file index
	int j;				// row index

	NameDataMap()
	{
		ifile = j = -1;
	}

	NameDataMap(string name_, int ifile_, int j_)
	{
		name = name_;
		ifile = ifile_;
		j = j_;
	}

	// comparison functor
	int operator() (const NameDataMap& lhs, const NameDataMap& rhs)
	{
		return (lhs.name < rhs.name);
	}
};

ChannelTwoPointNormalizationEngine::Reg ChannelTwoPointNormalizationEngine::reg;

ChannelTwoPointNormalizationEngine * ChannelTwoPointNormalizationEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == ChannelTwoPointNormalizationEngine::EngineName())
		return (ChannelTwoPointNormalizationEngine *)engine;
	return NULL;
}

ChannelTwoPointNormalizationEngine::ChannelTwoPointNormalizationEngine()
{
    defineOptions();
    defineStates();
}

ChannelTwoPointNormalizationEngine::~ChannelTwoPointNormalizationEngine()
{
}

void ChannelTwoPointNormalizationEngine::defineOptions()
{
	defineOption("", "in-summary-files", PgOpt::STRING_OPT,
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
	defineOption("", "compute-target", PgOpt::BOOL_OPT,
		"Compute the target by averaging medians from all channels.  target value will be ignored if compute-target is true.",
		"false");
}

void ChannelTwoPointNormalizationEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void ChannelTwoPointNormalizationEngine::checkOptionsImp()
{
	string file;
	file = getOpt("in-summary-files");
	if (file.empty() == true)
		Err::errAbort("Must specify input summary files.");
	vector<string> files = split(file);
	//if (files.size() < 2)
	//	Err::errAbort("Must specify at least two input summary files.");
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

void ChannelTwoPointNormalizationEngine::writeHeader(std::ofstream &out)
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
void ChannelTwoPointNormalizationEngine::runImp()
{
	vector<string> inFiles = split(getOpt("in-summary-files"));
	string countFile = getOpt("in-count-file");
	string mixFile = getOpt("marker-content-file");
	string outFile = getOpt("out-summary-file");
	int perc = getOptInt("percentile");
	int negTarget = getOptInt("neg-control-target");
	bool computeTarget = getOptBool("compute-target");
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

	// Open the count file and read the sample names.
	affx::TsvFile countTsv;
	if(countTsv.open(countFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + countFile);
	countTsv.rewind();

	// Get the column names from the counts file and save to the output file.
	ofstream out(outFile.c_str(), ios::out);
	writeHeader(out);
	string columnName;
	vector<string> columnNames;
	int n = countTsv.getColumnCount(0);
	for (int i=0; i<n; i++)
	{
		columnName = countTsv.getColumnName(0, i);
		columnNames.push_back(columnName);
		out << columnName;
		if (i < n - 1)
			out << "\t";
	}
	out << endl;

	// Allocate a vector of data for the row and bind to the columns.
	string name;
	vector<string> rowNames;
	vector<int> count(n-1);
	countTsv.bind(0, 0, &name, TSV_BIND_REQUIRED);
	for (int i=1; i<n; i++)
	{
		countTsv.bind(0, i, &count[i-1],  TSV_BIND_REQUIRED);
	}

	vector<bool> isControlRow;
	vector<vector<int> > countData(n-1);  
	vector<vector<int> > nonControlCountData(n-1); 
	while(countTsv.nextLevel(0) == TSV_OK)
	{
		rowNames.push_back(name);
		bool isControl = (controlMap.find(name) != controlMap.end());
		isControlRow.push_back(isControl);
		for (int i=0; i<n-1; i++)
		{
			countData[i].push_back(count[i]);
			if (isControl == false)
			{
				nonControlCountData[i].push_back(count[i]);
			}
		}
	}
	countTsv.clear();

	vector<double> data(n-1);
	vector<vector<double> > minData(inFiles.size());
	vector<vector<string> > names(inFiles.size());
	vector<NameDataMap> nameDataMap;
	vector<vector<vector<float> > > sampleData(inFiles.size()); 

	// Loop through the data files and read the data.
	for (int ifile = 0; ifile < inFiles.size(); ++ifile)
	{
		if(tsv.open(inFiles[ifile]) != TSV_OK)
			Err::errAbort("Couldn't open the file: " + inFiles[ifile]);
		tsv.rewind();

		// Check the column names
		int cols = tsv.getColumnCount(0);
		for (int i=0; i<cols; i++)
		{
			string columnName = tsv.getColumnName(0, i);
			if (columnName != columnNames[i])
			{
				// Column names in the counts file do not match those in this signals file.
				Err::errAbort("The counts and signals files do not appear to be consistent. Column names don't match. (Counts file column:'"
					+ columnNames[i] + "', " + inFiles[ifile] + " file column: '" + columnName +"')");
			}
		}

		tsv.bind(0, 0, &name, TSV_BIND_REQUIRED);
		for (int i=1; i<n; i++)
		{
			tsv.bind(0, i, &data[i-1],  TSV_BIND_REQUIRED);
		}

		minData[ifile].resize(n-1, std::numeric_limits<double>::max());

		// Get the data. Separate the control from non-control values.
		// Only retrieve the counts for the controls.
		int irow = 0;
		while(tsv.nextLevel(0) == TSV_OK)
		{
			// Check the row names.
			bool countNameHasCode = hasAlleleCode(rowNames[irow]);
			string nameNoAlleleCode = countNameHasCode ? name : removeAlleleCode(name);
			if (rowNames[irow] != nameNoAlleleCode)
			{
				// Row names in the counts file do not match those in this signals file.
				Err::errAbort("The counts and signals files do not appear to be consistent. Row names don't match. (Counts file row:'"
					+ rowNames[irow] + "', " + inFiles[ifile] + " file row (with Allele code): '" + name +"')");
			}

			sampleData[ifile].resize(n-1);

			bool isControl = (controlMap.find(name) != controlMap.end());
			if (isControl == false)
			{
				names[ifile].push_back(name);
				NameDataMap ndm(name, ifile, sampleData[ifile].size() > 0 ? sampleData[ifile][0].size() : 0);
				nameDataMap.push_back(ndm);
			}

			for (int i=0; i<n-1; i++)
			{
				if (countData[i].size() >= minCount)
					minData[ifile][i] = min(minData[ifile][i], data[i]);

				if (isControl == false)
				{
					sampleData[ifile][i].push_back(data[i]);
				}
			}

			++irow;
		}
		tsv.clear();
	}

	data.erase(data.begin(), data.end());
	count.erase(count.begin(), count.end());

	// Compute the weighted average of the controls and the
	// percentile of the non-controls.
	vector<bool> skipNorm(n-1);
	vector<float> percentiles(inFiles.size());
	vector<vector<float> > slope(n-1);
	vector<vector<float> > intercept(n-1);
	for (int i=0; i<n-1; i++)
	{
		slope[i].resize(inFiles.size());
		intercept[i].resize(inFiles.size());

		for (int ifile = 0; ifile < inFiles.size(); ++ifile)
		{
			vector<float> temp;
			for (int j=0; j<(int)nonControlCountData[i].size(); j++)
			{
				if (nonControlCountData[i][j] >= minCount)
					temp.push_back(sampleData[ifile][i][j]);
			}

			skipNorm[i] = (temp.size() < 2);	// set to the same value in each pass through the for loop. redundant but ok.

			if (skipNorm[i] == false)
			{
				std::sort(temp.begin(), temp.end());
				percentiles[ifile] = percentile(temp, perc);
			}
			else
			{
				percentiles[ifile] = 1.0f;	// just pass through the computation without blowing-up.  The value will not be used.
			}
		}

		float avgRange = target; //average(percentiles);
		float avgControl = negTarget; //average(wavg);

		if (computeTarget)
		{
			avgRange = 0;
			for (int ifile = 0; ifile < inFiles.size(); ++ifile)
			{
				avgRange += (percentiles[ifile] - minData[ifile][i]);
			}
			avgRange /= inFiles.size();
		}

		// Compute the slope and intercept for each sample.
		for (int ifile = 0; ifile < inFiles.size(); ++ifile)
		{
			slope[i][ifile] = (avgRange - avgControl) / (percentiles[ifile] - minData[ifile][i]);
			intercept[i][ifile] = avgControl - (slope[i][ifile] * minData[ifile][i]);
		}
	}

	// Sort the names
	sort(nameDataMap.begin(), nameDataMap.end(), NameDataMap());

	// Output the normalized signals
	// If there is not data then output a zero.
	for (int iname = 0; iname < nameDataMap.size(); ++iname)
	{
		int ifile = nameDataMap[iname].ifile;
		int j = nameDataMap[iname].j;

		out << nameDataMap[iname].name;
		for (int i=0; i<n-1; i++)
		{
			out << "\t";
			if (sampleData[ifile][i].empty() == false)
			{
				if (nonControlCountData[i][j] > 0)
				{
					if (skipNorm[i] == false)
					{
						out << (sampleData[ifile][i][j] * slope[i][ifile]) + intercept[i][ifile];
					}
					else
					{
						// no normalization
						out << sampleData[ifile][i][j];
					}
				}
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

