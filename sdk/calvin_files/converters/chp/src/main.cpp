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


/*! file main.cpp This file contains a command line program to convert a GCOS CHP file to Calvin format. */

#include "calvin_files/converters/chp/src/CHPFileConverter.h"
#include "calvin_files/converters/chp/src/CmdLine.h"
//
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
//

using namespace std;
using namespace affymetrix_chp_converter;

/*! The input file tag. */
#define INPUT_FILE_ARG "-i"

/*! The help tag. */
#define HELP_ARG "-h"

/*! The version tag. */
#define VERSION_ARG "-v"

/*! The library path tag. */
#define LIB_PATH_ARG "-l"

/*
 * Show a help message.
 */
void ShowHelp()
{
	cout
		<< "Command line arguments:"
		<< endl
		<< "\t" << INPUT_FILE_ARG << " <the full path of the input GCOS CHP file to convert>"
		<< endl
		<< "\t" << LIB_PATH_ARG << " <the path to where the library (CDF/PSI) files reside>"
		<< endl
		<< "\t" << VERSION_ARG << " <1 = convert to GCOS XDA, 2 = convert to Calvin>"
		<< endl;
}

/*
 * The options for the converter.
 */
typedef struct _ConverterOptions
{
	CHPFileVersionType version;
	string fileName;
	string libPath;

} ConverterOptions;

/*
 * Get the file name and library path from the command line arguments.
 * Show the help if required.
 */
bool GetArguments(int argc, char **argv, ConverterOptions &options)
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
		string str;
		options.fileName = cmdLine.GetArgument(INPUT_FILE_ARG, 0);
		options.libPath = cmdLine.GetArgument(LIB_PATH_ARG, 0);
		str = cmdLine.GetArgument(VERSION_ARG, 0);
		options.version = (CHPFileVersionType)(atoi(str.c_str()));
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
 * Extract the information from the Calvin data file and dump it to the command line.
 */
int main(int argc, char* argv[])
{
	// Get the input file name.
	ConverterOptions options;
	if (GetArguments(argc, argv, options) == false)
		return -1;

	// Test if it is a calvin file.
	CHPFileConverter converter;
	if (converter.ConvertFile(options.fileName.c_str(), options.libPath.c_str(), options.version) == false)
	{
		cout << CHPFileConverterErrorMessage(converter.ErrorCode()) << endl;
		return -1;
	}
	cout << "Conversion complete." << endl;
	return 0;
}
