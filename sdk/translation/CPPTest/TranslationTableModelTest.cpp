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
 * @file   translation/CPPTest/TranslationTableModelTest.cpp
 * @author Mybrid Spalding
 * @date   Thu Apr 10 12:44:56 PDT 2008
 *
 * @brief  Testing the TranslationTableModel.cpp class
 *
 */

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
 * @class TranslationTableModelTest
 * @brief cppunit class for testing CallSet behavior
 */
/*****************************************************************************/
class TranslationTableModelTest : public CppUnit::TestFixture
{

private:


  string m_defaultTranslationTableFileName;
  string m_missingHeaderLineFileName;
  string m_invaliddbSNPColumnHeaderFileName;
  string m_invalidGeneDataFileName;
  string m_invalidProbeSetDataFileName;
  string m_invaliddbSNPDataFileName;
  string m_invalidDefiningDataFileName;
  string m_invalidValidatedDataFileName;
  string m_invalidHaplotypeDataFileName;
  string m_invalidReferenceDataFileName;
  string m_invalidVariantDataFileName;

  string m_emptyGeneDataFileName;

  string m_programName;
  string m_outputDir;

  TranslationTableModel *m_ttm;
  RunTimeEnvironment m_rte;

  CPPUNIT_TEST_SUITE(TranslationTableModelTest);
  CPPUNIT_TEST(constructor);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void constructor();
};
// end class TranslationTableModelTest
/*****************************************************************************/
/*****************************************************************************/
/**
 * Initialization method called before each test routine.
 */
/*****************************************************************************/
void TranslationTableModelTest::setUp()
{

  cerr << endl;

  m_programName = "TranslationTableModeTest_setUp";

  m_defaultTranslationTableFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_EarlyAccess.txt";
  m_missingHeaderLineFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_missing_header.txt";

  m_invaliddbSNPColumnHeaderFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_dbSNP_column.txt";

  m_invalidGeneDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_gene_data.txt";
  m_invalidProbeSetDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_probeSet_data.txt";
  m_invaliddbSNPDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_dnsnp_data.txt";
  m_invalidDefiningDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_defining_data.txt";
  m_invalidValidatedDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_validated_data.txt";
  m_invalidHaplotypeDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_haplotype_data.txt";
  m_invalidReferenceDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_reference_data.txt";
  m_invalidVariantDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_invalid_variant_data.txt";

  m_emptyGeneDataFileName =
    TEST_DATA_UNIT_DIR + "TTable_v20080110_empty_gene_data.txt";

  m_outputDir = "output";

  m_rte.m_adtOpts.m_progName = m_programName;
  m_rte.m_adtOpts.m_outputDir = m_outputDir;

  m_rte.m_adtOpts.m_verbosity = ADT_VERBOSE_INPUT_FILES;

  m_rte.initializeRunTimeEnvironment();
  m_rte.m_adtOpts.m_inputTTableType = ADT_TRANSLATION_TABLE_TYPE_DMET2;

  return;
}
// end TranslationTableModelTest::setUp
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::TranslationTableModel constructor test cases.
 *
 */
/*****************************************************************************/

void TranslationTableModelTest::constructor()
{
  cout << endl;
  Util::PrintTextFunctionTitle("TranslationTableModelTest", "constructor");
  Err::setThrowStatus(true);


  // Missing TsvFile HEADER line.

  cerr << "constructor: missing TsvFile header line..." << endl;;
  //CPPUNIT_ASSERT( *csTestCases["marker1"] != *csTestCases["marker1"] );
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_missingHeaderLineFileName,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);

  cerr << "ok" << endl;

  cerr << "constructor: invalid dbSNP column header..." << endl;;
  //CPPUNIT_ASSERT( *csTestCases["marker1"] != *csTestCases["marker1"] );
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invaliddbSNPColumnHeaderFileName,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);

  cerr << "ok" << endl;

  cerr << "constructor: invalid Gene data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invalidGeneDataFileName,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid ProbeSet data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invalidProbeSetDataFileName,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid dbSNP data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invaliddbSNPDataFileName,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid Defining data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invalidDefiningDataFileName,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid Validated data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invalidValidatedDataFileName,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid Haplotype data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invalidHaplotypeDataFileName ,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid Reference data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invalidReferenceDataFileName ,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid Variant data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_invalidVariantDataFileName,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  cerr << "constructor: empty Gene data value..." << endl;
  NEGATIVE_TEST(TranslationTableModel(m_rte,
                       m_emptyGeneDataFileName ,
                       ADT_TRANSLATION_TABLE_TYPE_DMET2), Except);
  cerr << "ok" << endl;

  return;
}
// end TranslationTableModelTest::constructor
/*****************************************************************************/


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(TranslationTableModelTest);

////////////////////////////////////////////////////////////////
