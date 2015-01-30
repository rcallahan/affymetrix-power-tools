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
 * @file   CallElementTest.cpp
 * @author Mybrid Spalding
 * @date   Wed Apr  9 10:37:52 PDT 2008
 *
 * @brief  Testing the CallSet object.
 *
 */

#include "translation/CallSet.h"
#include "translation/RunTimeEnvironment.h"
#include "translation/TranslationTableModel.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h" // includes Verbose.h
#include "util/Util.h"
#include "util/CPPTest/Setup.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <map>
//


using namespace std;

/*****************************************************************************/
/**
 * @class CallSetTest
 * @brief cppunit class for testing CallSet behavior
 */
/*****************************************************************************/
class CallSetTest : public CppUnit::TestFixture
{

private:

  map<string, CallSet*> csTestCases;


  string m_defaultTranslationTableFileName;
  string m_programName;
  string m_outputDir;

  TranslationTableModel *m_ttm;
  RunTimeEnvironment m_rte;

  CPPUNIT_TEST_SUITE(CallSetTest);
  CPPUNIT_TEST(describeVerbose);
  CPPUNIT_TEST(operator_equal);
  CPPUNIT_TEST(addCallElement);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void operator_equal();
  void describeVerbose();
  void addCallElement();

};
// end class CallSetTest
/*****************************************************************************/
/*****************************************************************************/
/**
 * Initialization method called before each test routine.
 */
/*****************************************************************************/
void CallSetTest::setUp()
{


  cerr << endl;


  m_programName = "CallSetTest_setUp";
  m_defaultTranslationTableFileName =
    TEST_DATA_UNIT_DIR  + "TTable_v20080110_EarlyAccess.txt";

  m_outputDir = "output";

  m_rte.m_adtOpts.m_progName = m_programName;
  m_rte.m_adtOpts.m_outputDir = m_outputDir;

  m_rte.m_adtOpts.m_verbosity = ADT_VERBOSE_INPUT_FILES;

  m_rte.initializeRunTimeEnvironment();

  m_rte.m_adtOpts.m_inputTTableType = ADT_TRANSLATION_TABLE_TYPE_DMET2;

  m_ttm = new TranslationTableModel(m_rte,
                                    m_defaultTranslationTableFileName,
                                    ADT_TRANSLATION_TABLE_TYPE_DMET2);


  // Add MARKER Allele call sets as test cases.
  for (int row = 1; row < 24; row ++) {

    stringstream setNameSStr;
    //cerr << ttm->rows[row][0] << ttm->rows[row][1] << endl;

    setNameSStr << "marker" << row;
    csTestCases[setNameSStr.str()] = new CallSet(setNameSStr.str());
    csTestCases[setNameSStr.str()]
    ->addCallElement(*(m_ttm), row, m_ttm->getColumnIndex(ADT_DMET3_TT_REFERENCE)
                     , true, false);
  }

  // Add a HAPLOTYPE Allele call set as test cases.

  csTestCases["haplotype1"] = new CallSet("haplotype1");
  csTestCases["haplotype2"] = new CallSet("haplotype2");


  for (int row = 73; row < 80; row++) {

    //cerr << m_ttm->m_rows[row][0] << m_ttm->m_rows[row][1] << endl;
    vector<int> alleleColumns = m_ttm->getHaplotypeAlleleColumns(row);


    if (alleleColumns.size() > 0) {

      for (int column = 0; column < alleleColumns.size(); column++) {

        csTestCases["haplotype1"]
        ->addCallElement(*(m_ttm), row, alleleColumns[column], true, true);
        csTestCases["haplotype2"]
        ->addCallElement(*(m_ttm), row, alleleColumns[column], true, true);
      }
    }
  }

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::operator== test cases.
 *
 */
/*****************************************************************************/

void CallSetTest::operator_equal()
{
  cout << endl;
  Util::PrintTextFunctionTitle("CallSetTest", "operator_equal");
  Err::setThrowStatus(true);


  // isDescriptive, always false if either is not descriptive.
  cerr << "operator==, isDescriptive: a != a, !a.m_isDescriptive...";
  CPPUNIT_ASSERT(*csTestCases["marker1"] != *csTestCases["marker1"]);
  cerr << "ok" << endl;

  cerr << "operator==, isDescriptive: a == a, a.m_isDescriptive...";
  csTestCases["marker1"]->m_isDescriptive = true;
  CPPUNIT_ASSERT(*csTestCases["marker1"] == *csTestCases["marker1"]);
  cerr << "ok" << endl;


  csTestCases["marker3"]->m_isDescriptive = true;
  //csTestCases["marker3"]->describeVerbose(rte );
  //csTestCases["marker4"]->describeVerbose(rte );

  cerr << "operator==, isDescriptive: a != b, a.m_isDescriptive and !b.m_isDescriptive...";
  CPPUNIT_ASSERT(*csTestCases["marker3"] != *csTestCases["marker4"]);
  cerr << "ok" << endl;


  cerr << "operator==, isDescriptive: a == b, a.m_isDescriptive and b.m_isDescriptive...";
  csTestCases["marker4"]->m_isDescriptive = true;
  CPPUNIT_ASSERT(*csTestCases["marker3"] == *csTestCases["marker4"]);
  cerr << "ok" << endl;


  cerr << "operator==, set size not equal: a != b, a.size() != b.size()...";
  csTestCases["haplotype1"]->m_isDescriptive = true;
  CPPUNIT_ASSERT(*csTestCases["marker1"] != *csTestCases["haplotype1"]);
  cerr << "ok" << endl;

  cerr << "operator==, set size = 1: a == b...";
  CPPUNIT_ASSERT(*csTestCases["marker3"] == *csTestCases["marker4"]);
  cerr << "ok" << endl;

  cerr << "operator==, set size > 1: a == b...";
  csTestCases["haplotype2"]->m_isDescriptive = true;
  CPPUNIT_ASSERT(*csTestCases["haplotype1"] == *csTestCases["haplotype2"]);
  cerr << "ok" << endl;

  return;
}
// end CallSetTest::operator_equal
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::describeVerboe test cases.
 *
 */
/*****************************************************************************/
void CallSetTest::describeVerbose()
{
  cout << endl;
  Util::PrintTextFunctionTitle("CallSetTest", "describeVerbose");
  Err::setThrowStatus(true);

  m_rte.m_currentVerbosity = ADT_VERBOSE_INPUT_FILES ;

  Verbose::setLevel(ADT_VERBOSE_INPUT_FILES);

  cerr << "describeVerbose: marker1->describeVerbose(rte)..." << endl;
  csTestCases["marker1"]->describeVerbose(m_rte);
  cerr << "ok" << endl;

  cerr << "describeVerbose: haplotype1->describeVerbose(rte)..." << endl;
  csTestCases["haplotype1"]->describeVerbose(m_rte);
  cerr << "ok" << endl;

  return;
}
// end CallSetTest::describeVerbose
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::addCallElement test cases.
 *
 */
/*****************************************************************************/
void CallSetTest::addCallElement()
{
  cout << endl;
  Util::PrintTextFunctionTitle("CallSetTest", "addCallElement");
  Err::setThrowStatus(true);

  CallSet callset1("callset1");

  cerr << "addCallElement, duplicate CallElement, allowMultiAllelic = false: !callset1.addCallElement(duplicate)...";
  CPPUNIT_ASSERT(callset1.addCallElement(*(m_ttm), 1, m_ttm->getColumnIndex(ADT_DMET3_TT_REFERENCE)));
  CPPUNIT_ASSERT(!callset1.addCallElement(*(m_ttm), 1, m_ttm->getColumnIndex(ADT_DMET3_TT_REFERENCE)));
  cerr << "ok" << endl;

  cerr << "addCallElement, duplicate CallElement, allowMultiAllelic = true: callset1.addCallElement(duplicate)...";
  CPPUNIT_ASSERT(callset1.addCallElement(*(m_ttm), 1, m_ttm->getColumnIndex(ADT_DMET3_TT_REFERENCE), true));
  cerr << "ok" << endl;

  cerr << "addCallElement, invalid column dbSNP: callset1.addCallElement(DBSNPS_INDEX)...";
  NEGATIVE_TEST(callset1.addCallElement(*(m_ttm), 1, m_ttm->getColumnIndex(ADT_DMET3_TT_DBSNP)), Except);
  cerr << "ok" << endl;


  return;
}
// end CallSetTest::addCallElement
/*****************************************************************************/


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(CallSetTest);

////////////////////////////////////////////////////////////////
