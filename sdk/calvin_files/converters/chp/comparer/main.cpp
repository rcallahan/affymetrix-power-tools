////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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


#include "calvin_files/Extractor/src/CmdLine.h"
#include "calvin_files/converters/chp/comparer/CHPComparer.h"
#include "calvin_files/utils/src/FileUtils.h"
//
#include <cstring>
#include <iostream>
#include <string>
//

using namespace affymetrix_calvin_utilities;
using namespace affymetrix_comparer;

/*! The input file1 tag. */
#define INPUT_FILE1_ARG "-i1"

/*! The input file2 tag. */
#define INPUT_FILE2_ARG "-i2"

/*! The help tag. */
#define HELP_ARG "-h"

/*
 * Show a help message.
 */
void ShowHelp()
{
	cout
		<< "Command line arguments:" << endl
		<< "\t" << INPUT_FILE1_ARG << " <the full path of the reference CHP file>" << endl
		<< "\t" << INPUT_FILE2_ARG << " <the full path of the comparison CHP file>" << endl << endl
		<< "Possible return codes are:" << endl
		<< "\t" << " 0 [success]" << endl
		<< "\t" << "-1 [commandline argument errors]" << endl
		<< "\t" << "-2 [files do not compare]" << endl << endl;
}

/*! The options for the comparer. */
typedef struct _ComparerOptions
{
	/*! The first file. */
	string file1;

	/*! The second file. */
	string file2;

} ComparerOptions;

/*
 * Get the file name and library path from the command line arguments.
 * Show the help if required.
 */
bool GetArguments(int argc, char **argv, ComparerOptions &options)
{
	CCmdLine cmdLine;

	// parse argc,argv 
	// no switches were given on the command line, abort
	if (cmdLine.SplitLine(argc, argv) < 1)
	{
		ShowHelp();
		return false;
	}

	// test for the 'help' case
	if (cmdLine.HasSwitch(HELP_ARG))
	{
		ShowHelp();
		return false;
	}

	// get the required argument
	try
	{
		options.file1 = cmdLine.GetArgument(INPUT_FILE1_ARG, 0);
		options.file2 = cmdLine.GetArgument(INPUT_FILE2_ARG, 0);
	}
	catch (...)
	{
		// one of the required arguments was missing, abort
		ShowHelp();
		return false;
	}
	return true;
}

/*
 * Compare the two input files
 */
int main(int argc, char* argv[])
{
	// Get the input file name.
	ComparerOptions options;
	if (GetArguments(argc, argv, options) == false)
		return -1;

	// Compare the files
	CHPComparer comp;
	int ret = 0;
	if (comp.CompareFiles(options.file1.c_str(), options.file2.c_str()) == true)
	{
		cout << "The files are identical" << endl;
	}
	else
	{
		cout << comp.Differences() << endl;
		ret = -2;
	}

	return ret;
}

