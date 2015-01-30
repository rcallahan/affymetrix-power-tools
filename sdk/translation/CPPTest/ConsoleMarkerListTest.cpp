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
 * @file   translation/CPPTest/ConsoleMarkerListTest.cpp
 * @author Mybrid Spalding
 * @date   Thu Sep 25 09:38:39 PDT 2008
 * @brief  Test possible probe set id values in the marker list.
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
//

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
 * @class CopyMarkerListTest
 * @brief cppunit class for testing marker list probe sets
 */
/*****************************************************************************/
class ConsoleMarkerListTest : public CppUnit::TestFixture
{

private:

  map<string, string>        m_oneToOneOptions;
  vector< vector< string > > m_sampleTable;
  vector< string >           m_experimentListVector;
  vector< string >           m_markerList;

  RunTimeEnvironment rte;

  CPPUNIT_TEST_SUITE(ConsoleMarkerListTest);
  CPPUNIT_TEST(setUp);
  CPPUNIT_TEST(run);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void run();
};
// end ConsoleMarkerListTest
/*****************************************************************************/
/*****************************************************************************/
/**
 * Initialization method called before each test routine.
 *
 */
/*****************************************************************************/
void ConsoleMarkerListTest::setUp()
{

  cerr << endl;

  m_oneToOneOptions["out-dir"]
  = "output";

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

  row.clear();

  m_markerList.push_back("Probe Set ID");
  m_markerList.push_back("AM_14633");
  m_markerList.push_back("AM_14631");
  m_markerList.push_back("AM_14628");
  m_markerList.push_back("AM_14617");
  m_markerList.push_back("AM_14616");
  m_markerList.push_back("AM_14612");
  m_markerList.push_back("AM_14610");
  m_markerList.push_back("AM_14609");
  m_markerList.push_back("AM_14605");
  m_markerList.push_back("AM_14592");
  m_markerList.push_back("AM_14592");
  m_markerList.push_back("AM_14582");
  m_markerList.push_back("AM_14581");
  m_markerList.push_back("AM_14577");
  m_markerList.push_back("AM_14575");
  m_markerList.push_back("AM_14630");
  m_markerList.push_back("AM_14627");
  m_markerList.push_back("AM_14625");
  m_markerList.push_back("AM_14624");
  m_markerList.push_back("AM_14622");
  m_markerList.push_back("AM_14621");
  m_markerList.push_back("AM_14620");
  m_markerList.push_back("AM_14619");
  m_markerList.push_back("AM_14607");
  m_markerList.push_back("AM_14606");
  m_markerList.push_back("AM_14604");
  m_markerList.push_back("AM_14603");
  m_markerList.push_back("AM_14602");
  m_markerList.push_back("AM_14599");
  m_markerList.push_back("AM_14598");
  m_markerList.push_back("AM_14596");
  m_markerList.push_back("AM_14595");
  m_markerList.push_back("AM_14594");
  m_markerList.push_back("AM_14593");
  m_markerList.push_back("AM_14591");
  m_markerList.push_back("AM_14590");
  m_markerList.push_back("AM_14589");
  m_markerList.push_back("AM_14588");
  m_markerList.push_back("AM_14587");
  m_markerList.push_back("AM_14579");
  m_markerList.push_back("AM_14578");
  m_markerList.push_back("AM_10143");
  m_markerList.push_back("AM_10152");
  m_markerList.push_back("AM_10154");
  m_markerList.push_back("AM_10162");
  m_markerList.push_back("AM_10167");
  m_markerList.push_back("AM_10175");
  m_markerList.push_back("AM_10177");
  m_markerList.push_back("AM_10178");
  m_markerList.push_back("AM_10179");
  m_markerList.push_back("AM_10183");
  m_markerList.push_back("AM_10184");
  m_markerList.push_back("AM_10185");
  m_markerList.push_back("AM_10188");
  m_markerList.push_back("AM_10190");
  m_markerList.push_back("AM_10191");
  m_markerList.push_back("AM_10193");
  m_markerList.push_back("AM_10144");
  m_markerList.push_back("AM_10145");
  m_markerList.push_back("AM_10146");
  m_markerList.push_back("AM_10147");
  m_markerList.push_back("AM_10149");
  m_markerList.push_back("AM_10150");
  m_markerList.push_back("AM_10151");
  m_markerList.push_back("AM_10153");
  m_markerList.push_back("AM_10153");
  m_markerList.push_back("AM_10155");
  m_markerList.push_back("AM_10156");
  m_markerList.push_back("AM_10158");
  m_markerList.push_back("AM_10159");
  m_markerList.push_back("AM_10160");
  m_markerList.push_back("AM_10161");
  m_markerList.push_back("AM_10163");
  m_markerList.push_back("AM_10164");
  m_markerList.push_back("AM_10166");
  m_markerList.push_back("AM_10168");
  m_markerList.push_back("AM_10169");
  m_markerList.push_back("AM_10172");
  m_markerList.push_back("AM_10174");
  m_markerList.push_back("AM_10176");
  m_markerList.push_back("AM_10180");
  m_markerList.push_back("AM_10181");
  m_markerList.push_back("AM_10189");
  m_markerList.push_back("AM_10189");
  m_markerList.push_back("AM_10192");
  m_markerList.push_back("AM_11498");
  m_markerList.push_back("AM_11499");
  m_markerList.push_back("AM_11500");
  m_markerList.push_back("AM_11501");
  m_markerList.push_back("AM_11506");
  m_markerList.push_back("AM_11520");
  m_markerList.push_back("AM_10774");
  m_markerList.push_back("AM_10771");
  m_markerList.push_back("AM_10770");
  m_markerList.push_back("AM_10769");
  m_markerList.push_back("AM_10768");
  m_markerList.push_back("AM_10766");
  m_markerList.push_back("AM_10766");
  m_markerList.push_back("AM_10765");
  m_markerList.push_back("AM_10762");
  m_markerList.push_back("AM_10778");
  m_markerList.push_back("AM_10776");
  m_markerList.push_back("AM_10775");
  m_markerList.push_back("AM_10772");
  m_markerList.push_back("AM_10767");




  return;
}
// end ConsoleMarkerListTest::setUp
/*****************************************************************************/
/*****************************************************************************/
/**
 * ConsoleMarkerListTest::run
 *
 */
/*****************************************************************************/
void ConsoleMarkerListTest::run()
{

  Util::PrintTextFunctionTitle("ConsoleMarkerListTest", "run");

  map<string, string>::iterator itSS;

  TranslationEngine te(NULL, C_CMDLINE | C_CONSOLE, &translateEngineMessageHandlingCallBack);


  for (itSS  = m_oneToOneOptions.begin();
       itSS != m_oneToOneOptions.end();
       itSS++) {

    te.setOpt(itSS->first, itSS->second);
  }

  te.setOpt("experiment-list-vector", m_experimentListVector);
  te.setOpt("marker-list-vector", m_markerList);

  te.run();

}
// end ConsoleMarkerListTest::run
/*****************************************************************************/

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(ConsoleMarkerListTest);
