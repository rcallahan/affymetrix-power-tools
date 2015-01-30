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

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//
#include "chipstream/SelectFilesEngine.h"
//
#include "chipstream/EngineUtil.h"
//
#include "util/PgOptions.h"
//
#include <cstring>
#include <string>
#include <stdio.h>
#include <vector>


using namespace std;
using namespace affx;

const int NUM_INPUT_FILES = 10;

SelectFilesEngine::Reg SelectFilesEngine::reg;

SelectFilesEngine * SelectFilesEngine::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == SelectFilesEngine::EngineName())
		return (SelectFilesEngine *)engine;
	return NULL;
}

SelectFilesEngine::SelectFilesEngine() {
    defineOptions();
}

SelectFilesEngine::~SelectFilesEngine() {
}

void SelectFilesEngine::defineOptions() {
	defineOption("", "input-count", PgOpt::STRING_OPT,
		"The number of input files.",
		"");

	defineOption("", "out-file", PgOpt::STRING_OPT,
		"The name of the output file.",
		"");

	char selectInput[25];
	char inputFile[25];

	for (int i = 1; i <= NUM_INPUT_FILES; ++i)
	{
		sprintf(selectInput, "select-input-%d", i);
		sprintf(inputFile, "input-file-%d", i);

		defineOption("", selectInput, PgOpt::STRING_OPT,
			"Indicates if the file should be included.",
			"");

		defineOption("", inputFile, PgOpt::STRING_OPT,
			"Name of the file.",
			"");
	}
}

void SelectFilesEngine::defineStates() { }

/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void SelectFilesEngine::checkOptionsImp()
{
	defineStates();

	// Sanity check the inputs
	string filecount = getOpt("input-count");

	int count = atoi(filecount.c_str());

	if (count == 0)
	{
			Err::errAbort("Couldn't read input-count");
	}

	string outfile = getOpt("out-file");

	if (outfile.empty())
	{
			Err::errAbort("out-file parameter value is missing");
	}

	char selectInput[25];
	char inputFile[25];

	for (int i = 1; i <= count; ++i)
	{
		sprintf(selectInput, "select-input-%d", i);
		sprintf(inputFile, "input-file-%d", i);

		string selectInputValue = getOpt(selectInput);
		string inputFileValue = getOpt(inputFile);

		if (selectInputValue.empty())
		{
			string selectInputName(selectInput);
			Err::errAbort("Missing expected parameter - " + selectInputName);
		}

		if (inputFileValue.empty())
		{
			string inputFileName(inputFile);
			Err::errAbort("Missing expected parameter - " + inputFileName);
		}
	}
}

/**
   This is the "main()" equivalent of the engine.
*/
void SelectFilesEngine::runImp()
{
	string outfilename = getOpt("out-file");
	ofstream outfile;
        Fs::aptOpen(outfilename, ios::out);
	outfile << "file" << endl;

	string filecount = getOpt("input-count");
	int count = atoi(filecount.c_str());

	char selectInput[25];
	char inputFile[25];

	for (int i = 1; i <= count; ++i)
	{
		sprintf(selectInput, "select-input-%d", i);
		string selectInputValue = getOpt(selectInput);

		if (selectInputValue[0] == 'y' || selectInputValue[0] == 'Y')
		{
			sprintf(inputFile, "input-file-%d", i);
			string inputFileValue = getOpt(inputFile);

			outfile << inputFileValue << endl;
		}
	}

	outfile.close();
}
