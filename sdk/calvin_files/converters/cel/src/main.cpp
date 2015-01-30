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


/*! file main.cpp This file contains a command line program to convert a GCOS CEL file to Calvin format. */

#include "calvin_files/Extractor/src/CmdLine.h"
#include "calvin_files/converters/cel/src/CELFileConverter.h"
//
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
//

using namespace std;
using namespace affymetrix_cel_converter;

/*! The input file tag. */
#define INPUT_FILE_ARG "-i"

/*! The version tag. */
#define VERSION_ARG "-v"

/*! The help tag. */
#define HELP_ARG "-h"

/*
 * Show a help message.
 */
void ShowHelp()
{
	cout
		<< "Command line arguments:"
		<< endl
		<< "\t" << INPUT_FILE_ARG << " <the full path of the input GCOS CEL file to convert>"
		<< endl
		<< "\t" << VERSION_ARG << " <1 = convert to ASCII, 2 = convert to XDA, 3 = convert to Calvin>"
		<< endl;
}

/*
 * The options for the converter.
 */
typedef struct _ConverterOptions
{
	CELFileVersionType version;
	string fileName;

} ConverterOptions;

/*
 * Get the file name from the command line arguments.
 * Show the help if required.
 */
bool GetOptions(int argc, char **argv, ConverterOptions &options)
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
		str = cmdLine.GetArgument(INPUT_FILE_ARG, 0);
		options.fileName = str;

		str = cmdLine.GetArgument(VERSION_ARG, 0);
		options.version = (CELFileVersionType)(atoi(str.c_str()));
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
	if (GetOptions(argc, argv, options) == false)
		return -1;

	// Test if it is a calvin file.
	CELFileConverter converter;
	if (converter.ConvertFile(options.fileName.c_str(), options.version) == false)
	{
		cout << CELFileConverterErrorMessage(converter.ErrorCode()) << endl;
		return -1;
	}
	cout << "Conversion complete." << endl;
	return 0;
}
