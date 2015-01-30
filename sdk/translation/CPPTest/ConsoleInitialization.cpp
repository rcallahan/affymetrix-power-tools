////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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

/**
 * @file   translation/CPPTest/ConsoleInitialization.cpp
 * @author Mybrid Spalding
 * @date   Thu Aug 28 11:41:17 PDT 2008
 * @brief  Test the set option APIs of TranslateEngine not used command line.
 */

//
#include "translation/TranslationEngine.h"
//
#include "util/Err.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

using namespace std;

/*****************************************************************************/
/**
 * translateEngineMessageHandlingCallBack
 * Synopsis:
 *
 *  There is a catch 22 where the command line message handlers need to
 * have the PGOptions information in order to set up the message handling
 * log file name, etc. Resolve this by passing in a call back to the
 * TranslateEngine constructor which initializes the PGOptions.
 *
 * @param atdOopts - filled in with all the PGOptions data.
 *
 * @return - nada, initializes message handling.
 */
/*****************************************************************************/
static void translateEngineMessageHandlingCallBack(ADTOptions & adtOpts)
{

  Err::setThrowStatus(true);

  Fs::mkdirPath(adtOpts.m_outputDir, false);

  string verboseLogName = Fs::join(adtOpts.m_outputDir,adtOpts.m_progName + ".log");

  std::ofstream *verboseOFStream = new std::ofstream;

  Fs::mustOpenToWrite(*(verboseOFStream), verboseLogName.c_str());

  ADT_VERBOSE_ENUM logVerbosity = ADT_VERBOSE_EXCEPTION;
  if (adtOpts.m_verbosity > ADT_VERBOSE_EXCEPTION) {
    logVerbosity = (ADT_VERBOSE_ENUM) adtOpts.m_verbosity;
  }

  LogStream *verboseLogStream = new
  LogStream(logVerbosity,  verboseOFStream);

  Verbose::pushMsgHandler(verboseLogStream);
  Verbose::pushProgressHandler(verboseLogStream);
  Verbose::pushWarnHandler(verboseLogStream);
  Verbose::setLevel(adtOpts.m_verbosity);

  return;
}
// end translateEngineMessageHandlingCallBack
/*****************************************************************************/
/*****************************************************************************/
/**
 * @class CopyNumberTableTest
 * @brief cppunit class for testing CallSet behavior
 */
/*****************************************************************************/
class ConsoleInitializationTest : public CppUnit::TestFixture
{

private:

  map<string, string>        m_oneToOneOptions;
  vector< vector< string > > m_sampleTable;
  vector< string >           m_experimentListVector;

  RunTimeEnvironment rte;

  CPPUNIT_TEST_SUITE(ConsoleInitializationTest);
  CPPUNIT_TEST(setUp);
  CPPUNIT_TEST(run);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void run();
};
// end ConsoleInitializationTest
/*****************************************************************************/
/*****************************************************************************/
/**
 * Initialization method called before each test routine.
 *
 */
/*****************************************************************************/
void ConsoleInitializationTest::setUp()
{

  cerr << endl;

  m_oneToOneOptions["out-dir"] = "output";
  m_oneToOneOptions["translate-file"]  = TEST_DATA_UNIT_DIR + "DMET_Plus.r3.20081008.dc_translation.txt";
  m_oneToOneOptions["annotation-file"] =   TEST_DATA_UNIT_DIR + "DMET_Plus.r3.20081008.dc_annot.csv";
  m_experimentListVector.push_back(TEST_DATA_UNIT_DIR + "NA11832__HAPMAP_01A_0807_WC_Pilot2_37.dmet.chp");

  vector < string > row;

  row.push_back("Experiment Id");
  row.push_back("Sample Type");
  row.push_back("Consented Marker List");
  row.push_back("Research Code");
  row.push_back("Facility Code");
  row.push_back("Study Code");
  row.push_back("Investigator");
  row.push_back("Patient Number");

  m_sampleTable.push_back(row);
  row.clear();

  row.push_back("f1b8c758-cadf-4744-b5b3-10e82dcd110d");
  row.push_back("sample");
  row.push_back("DMETPlus_All");
  row.push_back("ACSL");
  row.push_back("");
  row.push_back("");
  row.push_back("J.Lo");
  row.push_back("10-019");

  m_sampleTable.push_back(row);


  return;
}
// end ConsoleInitializationTest::setUp
/*****************************************************************************/
/*****************************************************************************/
/**
 * ConsoleInitializationTest::run
 *
 */
/*****************************************************************************/
void ConsoleInitializationTest::run()
{

  Util::PrintTextFunctionTitle("ConsoleInitializationTest", "run");

  map<string, string>::iterator itSS;

  TranslationEngine te(NULL, C_CMDLINE | C_CONSOLE, &translateEngineMessageHandlingCallBack);


  for (itSS  = m_oneToOneOptions.begin();
       itSS != m_oneToOneOptions.end();
       itSS++) {

    te.setOpt(itSS->first, itSS->second);
  }

  te.setOpt("experiment-list-vector", m_experimentListVector);
  te.setOpt("sample-table", m_sampleTable);

  te.run();

}
// end ConsoleInitializationTest::run
/*****************************************************************************/

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(ConsoleInitializationTest);
