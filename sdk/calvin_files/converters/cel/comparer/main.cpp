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

//#include "CmdLine.h"
#include "calvin_files/converters/cel/comparer/CELFileComparer.h"
//
#include <cstring>
#include <iostream>
#include <string>
//

#define DEFAULT_TOLERANCE 0.1f

using namespace std;
using namespace affymetrix_comparers;

/*! The input file1 tag. */
#define INPUT_FILE1_ARG "-i1"

/*! The input file2 tag. */
#define INPUT_FILE2_ARG "-i2"

/*! The help tag. */
#define HELP_ARG "-h"

/*! The possible error messages */
const char *g_Error[] = 
{
	"No errors.",
	"The input file is already a Calvin CEL file.",
	"Unable to open the input CEL file.",
	"Unable to write the output CEL file.",
	"Comparison of files failed. Check comparison file (filename + .comparison) for results"
};

/*
 * Show a help message.
 */
void ShowHelp()
{
	cout
		<< "Usage:"
		<< endl << endl
		<< "CelComparer.exe [optional arguments] "
		<< "<the full path of the input GCOS CEL file>"
		<< " "
		<< "<the full path of the input Calvin CEL file>"
		<< endl << endl;


	cout
		<< "Optional arguments:" << endl
		<< "\t" << "-f <Write detailed Diff file with all differences>" << endl
		<< "\t" << "-cs <Compare Mask values to CelSummaryReport values in comparison>" << endl
		<< "\t" << "-im <Ignore Mask section differences in comparison>" << endl
		<< "\t" << "-ih <Ignore Headers section differences in comparison>" << endl
		<< "\t" << "-t TOLERANCE <Specify tolerance for float comparisons, default is +- 0.1>" << endl << endl
		<< "Possible return codes are:" << endl
		<< "\t" << " 0 [success]" << endl
		<< "\t" << "-1 [commandline argument errors]" << endl
		<< "\t" << "-2 [files do not compare]" << endl << endl;

}

/*
 * Show the error.
 */
void ShowError(int id)
{
	cout << g_Error[id] << endl;
}

/*
 * Get the file name from the command line arguments.
 * Show the help if required.
 * 
 * Set the detailedDiffFile flag to true if the -d argument is found
 */
bool ParseInputParameters(int argc, char **argv, string& fileNameGCOS, string& fileNameCalvin, bool& detailedDiffFile, float& tolerance, bool& ignoreMask, bool& celSummary, bool& ignoreHeaders)
{
	// parse argc,argv 
	// no switches were given on the command line, abort
	if ( (argc < 3) || (argc > 6) )
	{
		ShowHelp();
		return false;
	}

	// get the required argument
	try
	{
		if (argc > 3)
		{
			for (int i=1; i < argc -2; i++)
			{
				string val = string(argv[i]);

				if ( val.compare("-f") == 0)
					detailedDiffFile = true;

				if ( val.compare("-t") == 0)
					tolerance = (float) atof(argv[i + 1]);

				if ( val.compare("-im") == 0)
					ignoreMask = true;

				if ( val.compare("-cs") == 0)
					celSummary = true;

				if ( val.compare("-ih") == 0)
					ignoreHeaders = true;
			}

			fileNameGCOS = string(argv[2]);
			fileNameCalvin = string(argv[3]);
		}

		fileNameGCOS = string(argv[argc - 2]);
		fileNameCalvin = string(argv[argc - 1]);
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
	string fileNameGCOS;
	string fileNameCalvin;
	bool detailedDiffFile = false;
	bool ignoreMask = false;
	bool celSummary = false;
	float tolerance = DEFAULT_TOLERANCE;
	bool ignoreHeaders = false;

	if (ParseInputParameters(argc, argv, fileNameGCOS, fileNameCalvin, detailedDiffFile, tolerance, ignoreMask, celSummary, ignoreHeaders) == false)
	{
		return -1;
	}

	// Test if it is a calvin file.
	CELFileComparer comparer(detailedDiffFile, tolerance, ignoreMask, celSummary, ignoreHeaders);
	comparer.SetFileNameGCOS(fileNameGCOS.c_str());
	comparer.SetFileNameCalvin(fileNameCalvin.c_str());
	int ret = 0;
	if (comparer.CompareFiles() == false)
	{
		cout << "Comparison failed." << endl;
		cout << "See comparison result file: " << comparer.GetOutputFileName() << endl;
		ret = -2;
	}
	else
	{
		cout << "Files compare: Success" << endl;
	}
	return ret;
}
