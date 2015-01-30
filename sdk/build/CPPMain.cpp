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

//
#include "util/Err.h"
#include "util/PgOptions.h"
//
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/XmlOutputter.h>
//
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

//extern std::string externalDataPath=".";
std::string externalDataPath=".";

/** 
 * Everbybodys favorite function. 
 * @param argc 
 * @param argv
 * 
 * @return 
 */
int main(int argc, char *argv[])
{
  bool wasSucessful;
  try {
    // Command line options
    PgOptions opts;
    opts.defineOption("", "xml", PgOpt::BOOL_OPT,
                        "Use XML output.",
                        "false");
    opts.defineOption("", "text", PgOpt::BOOL_OPT,
                        "Use text output.",
                        "true");
    opts.defineOption("i", "external-data", PgOpt::BOOL_OPT,
                        "Location of external data.",
                        ".");
    opts.parseArgv(argv);
    externalDataPath = opts.get("external-data");

    // this doesnt actually work, as there isnt an error handler yet...
    //   Throw exceptions rather than calling exit()
    //   Err::setThrowStatus(true);
    //   Err::setExitOnError(false);
    //   Err::setExitOnErrorValue(-1);
    // Rather we just need to create one with the correct options.
    // throw     -> true
    // verbose   -> true
    // exit      -> false
    // exitValue -> -1 (unused as we arent exiting.)
    Err::configureErrHandler(true,true,false,-1);

    // Adds the test to the list of test to run
    CppUnit::TextUi::TestRunner runner;
    CppUnit::CompilerOutputter *textOutput = NULL;
    ofstream outputFile;
    if(opts.getBool("xml")) {
        // Define the file that will store the XML output.
        //Verbose::out(0, "using xml");
        outputFile.open("CPPUnitTestResults.xml");
    
        // Specify XML output and inform the test runner of this format.
        CppUnit::XmlOutputter* outputter =
            new CppUnit::XmlOutputter(&runner.result(), outputFile);
        runner.setOutputter(outputter);
    } 
    if(opts.getBool("text")) {
        //Verbose::out(0, "using standard outputter");
        // Change the default outputter to a compiler error format outputter
        textOutput = new CppUnit::CompilerOutputter( &runner.result(), std::cerr );
    }
    
    // Get the top level suite from the registry and add tests to runner
    CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();
    runner.addTest( suite );
    
    // Run the tests.
    wasSucessful = runner.run();
    if(opts.getBool("text")) {
      textOutput->write();
      delete textOutput;
    }
        
    if(opts.getBool("xml")) {
        outputFile.close();
    }
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
    
  // Return error code 1 if the one of test failed.
  return wasSucessful ? 0 : 1;
}
