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
 * @file   translation/CPPTest/ConsoleOverrideTest.cpp
 * @author Mybrid Spalding
 * @date   Thu Sep 25 09:37:22 PDT 2008
 * @brief  Test the possible reference and variant values in override.
 */

//
#include "translation/TranslationEngine.h"
//
#include "util/Err.h"
#include "util/Fs.h"
#include "util/LogStream.h"
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

  string verboseLogName = Fs::join(adtOpts.m_outputDir,adtOpts.m_progName+".log");

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
 * @class ConsoleOverrideTest
 * @brief cppunit class for testing possible reference and variant values.
 */
/*****************************************************************************/
class ConsoleOverrideTest : public CppUnit::TestFixture
{

private:

  map<string, string>        m_oneToOneOptions;
  vector< vector< string > > m_sampleTable;
  vector< string >           m_experimentListVector;
  vector< string >           m_markerList;

  RunTimeEnvironment rte;

  CPPUNIT_TEST_SUITE(ConsoleOverrideTest);
  CPPUNIT_TEST(setUp);
  CPPUNIT_TEST(run);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void run();
};
// end ConsoleOverrideTest
/*****************************************************************************/
/*****************************************************************************/
/**
 * Initialization method called before each test routine.
 *
 */
/*****************************************************************************/
void ConsoleOverrideTest::setUp()
{

  cerr << endl;

  m_oneToOneOptions["out-dir"]      = "output";
  m_oneToOneOptions["translate-file"]  = TEST_DATA_UNIT_DIR + "DMET_Plus.r3.20081008.dc_translation.txt";

  m_oneToOneOptions["annotation-file"] =   TEST_DATA_UNIT_DIR + "DMET_Plus.r3.20081008.dc_annot.csv";

  m_oneToOneOptions["genotype-override-file"] = TEST_DATA_UNIT_DIR + "NA11832__HAPMAP_01A_0807_WC_Pilot2_37.dmet_uncalled.rpt";

  m_experimentListVector.push_back(TEST_DATA_UNIT_DIR + "NA11832__HAPMAP_01A_0807_WC_Pilot2_37.dmet.chp");

  return;
}
// end ConsoleOverrideTest::setUp
/*****************************************************************************/
/*****************************************************************************/
/**
 * ConsoleOverrideTest::run
 *
 */
/*****************************************************************************/
void ConsoleOverrideTest::run()
{

  Util::PrintTextFunctionTitle("ConsoleOverrideTest", "run");

  map<string, string>::iterator itSS;

  TranslationEngine te(NULL, C_CMDLINE | C_CONSOLE, &translateEngineMessageHandlingCallBack);


  for (itSS  = m_oneToOneOptions.begin();
       itSS != m_oneToOneOptions.end();
       itSS++) {

    te.setOpt(itSS->first, itSS->second);
  }

  te.setOpt("experiment-list-vector", m_experimentListVector);

  te.run();

}
// end ConsoleOverrideTest::run
/*****************************************************************************/

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(ConsoleOverrideTest);
