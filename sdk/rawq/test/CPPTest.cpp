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

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <string.h>
#include <string>
//

extern std::string externalDataPath=".";

/*! Gets the file name from the command line argument.
 *
 * @param argc The number of arguments
 * @param argv The command line arguments
 * @param i The index to the command line arguments array.
 * @return The file name from the command line.
 */
std::string get_file(int argc, char *argv[], int &i)
{
	std::string file = "";
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
	int i=1;
	while(i<argc)
	{
		if (strcmp(argv[i], "-i") == 0)
			externalDataPath = get_file(argc, argv, i);
		++i;
	}

	// Get the top level suite from the registry
	CPPUNIT_NS::Test *suite = CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest();

	// Adds the test to the list of test to run
	CPPUNIT_NS::TextUi::TestRunner runner;
	runner.addTest( suite );

	// Change the default outputter to a compiler error format outputter
	runner.setOutputter( new CPPUNIT_NS::CompilerOutputter( &runner.result(), std::cerr ) );
	// Run the test.
	bool wasSucessful = runner.run();

	// Return error code 1 if the one of test failed.
	return wasSucessful ? 0 : 1;
}

