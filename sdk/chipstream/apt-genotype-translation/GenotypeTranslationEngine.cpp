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
#include "chipstream/apt-genotype-translation/GenotypeTranslationEngine.h"
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

pair<string, string> decodeReportStandAlleleCodePair(string ac, string dsa)
{
	string delim = "//";

	vector<string> alleleCodes = split(ac, delim);
	vector<string> reportAllele = split(dsa, delim);

	if (alleleCodes.size() != 2)
	{
		Err::errAbort("Report Strand Allele decoding expects two channels.");
	}

	if (alleleCodes.size() != reportAllele.size())
	{
		Err::errAbort("Allele Code-Report Strand Allele mismatch - Allele Code: " + ac + ", Report Strand Allele:" + dsa);
	}

	if (alleleCodes[0] == "A")
	{
		return make_pair<string, string>(reportAllele[0], reportAllele[1]);
	}
	else
	{
		return make_pair<string, string>(reportAllele[1], reportAllele[0]);
	}
}

GenotypeTranslationEngine::Reg GenotypeTranslationEngine::reg;

GenotypeTranslationEngine * GenotypeTranslationEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == GenotypeTranslationEngine::EngineName())
		return (GenotypeTranslationEngine *)engine;
	return NULL;
}

GenotypeTranslationEngine::GenotypeTranslationEngine()
{
    defineOptions();
}

GenotypeTranslationEngine::~GenotypeTranslationEngine()
{
}

void GenotypeTranslationEngine::defineOptions()
{
	defineOption("", "calls-file", PgOpt::STRING_OPT,
		"Text file specifying genotype calls process.",
		"");
	defineOption("", "annotation-file", PgOpt::STRING_OPT,
		"The annotation file.",
		"");
	defineOption("", "translation-out-file", PgOpt::STRING_OPT,
		"File to write allele calls results into.",
		"");
	defineOption("", "probeset-id-col", PgOpt::STRING_OPT,
		"The probe set id column header.",
		"Probeset Id");
	defineOption("", "report-strand-allele-col", PgOpt::STRING_OPT,
		"The design strand allele column header.",
		"Report Strand Allele");
	defineOption("", "allele-code-col", PgOpt::STRING_OPT,
		"The allele code column header.",
		"Allele Code");
}

void GenotypeTranslationEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void GenotypeTranslationEngine::checkOptionsImp()
{
    defineStates();
	string val = getOpt("calls-file");
	if (val.empty() == true)
		Err::errAbort("Must specify a call file.");
	
	val = getOpt("annotation-file");
	if (val.empty() == true)
		Err::errAbort("Must specify an annotation file.");

	val = getOpt("translation-out-file");
	if (val.empty() == true)
		Err::errAbort("Must specify an output file.");

	val = getOpt("probeset-id-col");
	if (val.empty() == true)
		Err::errAbort("probeset id column.");

	val = getOpt("report-strand-allele-col");
	if (val.empty() == true)
		Err::errAbort("Must specify a design strand allele column.");

	val = getOpt("allele-code-col");
	if (val.empty() == true)
		Err::errAbort("Must specify an allele code column.");
}

void GenotypeTranslationEngine::writeHeader(std::ofstream &out)
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
void GenotypeTranslationEngine::runImp()
{
	// Check and get the options
	string callFile = getOpt("calls-file");
	string annotFile = getOpt("annotation-file");
	string outFile = getOpt("translation-out-file");
	string psidCol = getOpt("probeset-id-col");
	string dsaCol = getOpt("report-strand-allele-col");
	string acCol = getOpt("allele-code-col");

	// Open the input annotation file and create a map of probeset id to design strand allele.
	Verbose::out(1, "Opening annotation file");
	affx::TsvFile tsv;
	if(tsv.open(annotFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + annotFile);
	int numChannels = 1;
	if (tsv.getHeader("num-channels", numChannels) != TSV_OK)
	{
		// If the header is missing, assume it is one channel
		numChannels = 1;
	}
	tsv.rewind();
	map<string, pair<string, string> > translationMap;
	string psid;
	string dsa;
	string ac;
	tsv.bind(0, psidCol, &psid, TSV_BIND_REQUIRED);
	tsv.bind(0, dsaCol, &dsa, TSV_BIND_REQUIRED);
	tsv.bind(0, acCol, &ac, TSV_BIND_REQUIRED);
	Verbose::out(1, "Reading annotation file");
	while(tsv.nextLevel(0) == TSV_OK)
	{
		if (translationMap.find(psid) == translationMap.end())
			translationMap[psid] = make_pair<string, string>("", "");

		if (numChannels == 1)
		{
			if (ac == "A")
				translationMap[psid].first = dsa;
			else
				translationMap[psid].second = dsa;
		}
		else
		{
			translationMap[psid] = decodeReportStandAlleleCodePair(ac, dsa);
		}
	}
	tsv.clear();

	// Open the input call file
	Verbose::out(1, "Opening call file");
	if(tsv.open(callFile) != TSV_OK)
		Err::errAbort("Couldn't open the file: " + callFile);
	tsv.rewind();

	// Create an output file and write the header
	ofstream out(outFile.c_str(), ios::out);
	writeHeader(out);
	out << "#Allele Calls" << endl;
	int n = tsv.getColumnCount(0);
	for (int i=0; i<n; i++)
	{
		out << tsv.getColumnName(0, i);
		if (i < n - 1)
			out << "\t";
	}
	out << endl;

	// Bind the calls.
	vector<int> calls(n-1);
	string marker;
	tsv.bind(0, 0, &marker, TSV_BIND_REQUIRED);
	for (int i=1; i<n; i++)
		tsv.bind(0, i, &calls[i-1], TSV_BIND_REQUIRED);

	// Translate the calls and output the results
	while(tsv.nextLevel(0) == TSV_OK)
	{
		bool found = true;
		map<string, pair<string, string> >::iterator it = translationMap.find(marker);
		if (it == translationMap.end())
		{
			out << marker;
			for (int i=0; i<n-1; i++)
			{
				out << "\t";
				switch (calls[i])
				{
				case 0:
					out << "<AA>";
					break;
				case 1:
					out << "<AB>";
					break;
				case 2:
					out << "<BB>";
					break;
				default:
					out << "<NoCall>";
					break;
				}
			}
			out << endl;
		}
		else
		{
			out << marker;
			const pair<string, string> &tt = it->second;
			for (int i=0; i<n-1; i++)
			{
				out << "\t";
				switch (calls[i])
				{
				case 0:
					out << tt.first << "/" << tt.first;
					break;
				case 1:
					out << tt.first << "/" << tt.second;
					break;
				case 2:
					out << tt.second << "/" << tt.second;
					break;
				default:
					out << "NoCall";
					break;
				}
			}
			out << endl;
		}
 	}

	// Close the input and output files.
	out.close();
	tsv.close();
}
