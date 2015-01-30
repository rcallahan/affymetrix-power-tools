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
#include "chipstream/FileListFilterEngine.h"
//
#include "chipstream/EngineUtil.h"
#include "file/TsvFile/TsvFile.h"
#include "util/AffxConv.h"
#include "util/Fs.h"
#include "util/PgOptions.h"

using namespace std;

FileListFilterEngine::Reg FileListFilterEngine::reg;

/*
 * A class to define a test.
 */
typedef struct _TestInputsType
{
	string name;	/// The name of the metric to test.
	bool has_data;	/// Flag indicating if the test should be made.
	string op;		/// The operator.
	float value;	/// The value to test.
	float thr;		/// The threshold to compare the value to.
} TestInputsType;

/*
 * Conver the pointer.
 */
FileListFilterEngine * FileListFilterEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == FileListFilterEngine::EngineName())
		return (FileListFilterEngine *)engine;
	return NULL;
}

/*
 * Construct the class by defining the options and state.
 */
FileListFilterEngine::FileListFilterEngine()
{
    defineOptions();
}

/*
 * Destruct the class
 */
FileListFilterEngine::~FileListFilterEngine()
{
}

/*
 * Defines each of the parameters for this engine.
 */
void FileListFilterEngine::defineOptions()
{
	defineOption("", "list-files-in", PgOpt::STRING_OPT,
		"Text file specifying list files to process, one per line with the first line being 'list_files'.",
		"");
	defOptMult("l", "lists", PgOpt::STRING_OPT,
		"List files to process.",
        "");
	defineOption("", "list-files-out", PgOpt::STRING_OPT,
		"Text file specifying list files to process, one per line with the first line being 'list_files'.",
		"");
	defineOption("", "files-col", PgOpt::STRING_OPT,
		"The column name for the input files.",
        "list_files");
	defOptMult("", "test", PgOpt::STRING_OPT,
		"The test to apply to filter the list files. The format is <metric>.<operator>.<threshold> where operator is eq, gt, ge, lt, le, eq, ne.",
		"Call Rate:ge:90");
	defineOption("", "test-sep", PgOpt::STRING_OPT,
		"The separator to use for delimiting tests.",
		":");
	defineOption("", "filter-lists", PgOpt::BOOL_OPT,
		"True to filter the list based on the tests.",
        "true");
}

/*
 * Check the options to make sure that each required one is properly set.
 */
void FileListFilterEngine::checkOptionsImp()
{
	string str = getOpt("list-files-out");
	if (str.empty() == true)
		Err::errAbort("Must specify an output file.");
	vector<string> tests = getOptVector("test");
	if (tests.empty() == true)
		Err::errAbort("Must specify a test.");

	vector<string> listFiles;
	if(getOpt("list-files-in").empty() == false)
	{
		affx::TsvFile tsv;
#ifdef WIN32
		tsv.m_optEscapeOk = false;
#endif
		std::string listFilesFile = getOpt("list-files-in");
		string file;
		tsv.bind(0, getOpt("files-col"), &file, affx::TSV_BIND_REQUIRED);
		try
		{
			if(tsv.open(listFilesFile) != affx::TSV_OK)
			{
				Err::errAbort("Couldn't open the input file: " + listFilesFile);
			}
		}
		catch (Except&)
		{
			tsv.close();	// need to close the file
			throw;
		}
		tsv.rewind();
		while(tsv.nextLevel(0) == affx::TSV_OK)
		{
			listFiles.push_back(file);
		}
		tsv.close();
		Verbose::out(2, "Read " + ToStr(listFiles.size()) + " list files from: " + Fs::basename(listFilesFile));
	}
	else
	{
		listFiles = getOptVector("lists");
	}
	setOpt("lists",listFiles);
}

/*
 * Performs a test on a value compared to a threshold.
 * The parameters include the operator to act on the value and its threshold.
 * The operators include:
	 ge for greater than or equal to.
	 gt for greater than.
	 le for less than or equal to.
	 lt for less than.
	 eq equal to.
	 ne not equal to.
*/
static bool TestFloatValue(float value, const string &op, float thr)
{
	if (op == "ge" && value < thr)
	{
		return false;
	}
	else if (op == "gt" && value <= thr)
	{
		return false;
	}
	else if (op == "le" && value > thr)
	{
		return false;
	}
	else if (op == "lt" && value >= thr)
	{
		return false;
	}
	else if (op == "eq" && fabs(value - thr) > std::numeric_limits<float>::epsilon())
	{
		return false;
	}
	else if (op == "ne" && fabs(value - thr) <= std::numeric_limits<float>::epsilon())
	{
		return false;
	}
	return true;
}

/*
 * Split the string into a vector based on the input separator.
 */
static vector<string> split(const string &inputString, const string &sep = ":")
{
	vector<string> tokens;
	size_t substrBegin = 0;
	for (;;)
	{
		size_t substrEnd = inputString.find (sep, substrBegin);
		if (substrEnd == string::npos)
		{
			// No more separators - save what's left, quit.
			string subString = inputString.substr (substrBegin);
			// Avoid returning a null string from a terminating separator or an empty inputString.
			if (! subString.empty())
				tokens.push_back (subString);
			break;
		}
		// Avoid null strings from an initial separator(s).
		if (substrEnd != substrBegin)
			tokens.push_back (inputString.substr (substrBegin, substrEnd - substrBegin) );
		// Continue following the separator
		substrBegin = substrEnd + 1;
	}
	return tokens;
}

/*
 * Parse the inputs into the test class.
 */
static vector<TestInputsType> ParseTests(const vector<string> &testOpt, const string &sep)
{
	vector<TestInputsType> tests;
	for (vector<string>::const_iterator it=testOpt.begin(); it!=testOpt.end(); it++)
	{
		vector<string> tokens = split(*it, sep);
		TestInputsType test;
		test.has_data = false;
		test.name = tokens[0];
		test.op = tokens[1];
		test.thr = (float)getDouble(tokens[2]);
		tests.push_back(test);
	}
	return tests;
}

/*
 * Read the QC metrics from the input file
 */
void ExtractMetricsFromFile(const string &file, vector<TestInputsType> &tests)
{
	// Initialize the has data flag
	for (vector<TestInputsType>::iterator testIt=tests.begin(); testIt!=tests.end(); testIt++)
		testIt->has_data = false;

	// If failed to open then call Err::errAbort to exit (throw exception).
	// TODO

	// Extract the metrics from the input file and match them to the tests object.
	// Set the <value> property of the TestInputType object.
}

/*
 * Run the analsis.
 */
void FileListFilterEngine::runImp()
{
    checkOptions();
	vector<string> listFiles = getOptVector("lists");
	vector<string> filteredFiles;
	bool filter = getOptBool("filter-lists");
	if (filter == true)
	{
		// Convert the string of tests to a map of parameter name to structure.
		string sep = getOpt("test-sep");
		vector<string> testOpt = getOptVector("test");
		vector<TestInputsType> tests = ParseTests(testOpt, sep);

		// Loop over each input file and test if it passes the tests. If it does
		// then add it to the output list.
		for (vector<string>::iterator listIt=listFiles.begin(); listIt!=listFiles.end(); listIt++)
		{
			bool passed = true;
			ExtractMetricsFromFile(*listIt, tests);
			for (vector<TestInputsType>::iterator testIt=tests.begin(); testIt!=tests.end() && passed == true; testIt++)
				passed = (testIt->has_data == true ? TestFloatValue(testIt->value, testIt->op, testIt->thr) : false);
			if (passed == true)
				filteredFiles.push_back(*listIt);
			Verbose::out(2, *listIt + " : " + (passed == true ? "passed" : "failed"));
		}
	}
	else
		filteredFiles = listFiles;

	// Create the output file.
	string outFile = getOpt("list-files-out");
	ofstream out(outFile.c_str(), ios::out);
	out << getOpt("files-col") << endl;
	for (vector<string>::iterator it=filteredFiles.begin(); it!=filteredFiles.end(); it++)
		out << *it << endl;
	out.close();
}
