////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

//
#include <cstring>
#include <iostream>
#include <string.h>
#include <string>
//

extern void test_file_writers(std::string file, std::string cdf);
extern void test_file_readers(std::string file);

std::string get_file(int argc,const char *argv[], int &i)
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


int main(int argc, const char* argv[])
{
	if (argc == 1)
	{
		std::cout << "Synopsis: file-sdk-example <-w> <-i input_file> <-cdf cdf_file_name>" << std::endl;
		return 0 ;
	}

	int i=1;
	bool bRead = true;
	std::string file;
	std::string cdf = "C:\\Program Files\\Affymetrix\\GeneChip\\Affy_Data\\Library\\Test3.CDF";
	while(i<argc)
	{
		if (strcmp(argv[i], "-w") == 0)
			bRead = false;

		else if (strcmp(argv[i], "-i") == 0)
			file = get_file(argc, argv, i);

		else if (strcmp(argv[i], "-cdf") == 0)
			cdf = get_file(argc, argv, i);

		++i;
	}

	if (file.empty() == true)
		return 0;

	if (bRead)
		test_file_readers(file);
	else
		test_file_writers(file, cdf);

	return 0;
}
