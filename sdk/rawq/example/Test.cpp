////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#include "rawq/src/RawQWorkflow.h"
//
#include <cstring>
#include <iostream>
#include <list>
#include <string.h>
#include <string>
//

std::string get_file(int argc, char *argv[], int &i)
{
	std::string file = "";

	// Get the file name.
	int j=i+1;
	while (j<argc && argv[j][0] != '-')
	{
		file += argv[j];
		if (j<argc-1 && argv[j+1][0] != '-')
			file += " ";
		++j;
	}
	i = j-1;

	return file;
}

int main(int argc, char* argv[])
{
	// Get the inputs from the command line.
	std::string libPath;
	std::list<std::string> celFiles;
	int i=1;
	while(i<argc)
	{
		if (strcmp(argv[i], "-l") == 0)
			libPath = get_file(argc, argv, i);

		else if (strcmp(argv[i], "-i") == 0)
			celFiles.push_back(get_file(argc, argv, i));

		++i;
	}

	CRawQWorkflow q;
	std::list<std::string>::iterator iter;
	for (iter = celFiles.begin(); iter != celFiles.end(); ++iter)
	{
		std::string celFile = *iter;
		CRawQWorkflow::RawQError error = q.ComputeRawQ(celFile.c_str(), libPath.c_str());
		if (error == CRawQWorkflow::NoError)
		{
			std::cout << "RawQ = " << q.GetRawQ() << std::endl;
		}
	}
	return 0;
}



