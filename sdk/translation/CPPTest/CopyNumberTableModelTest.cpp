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
 * @file   translation/CPPTest/CopyNumberTableModelTest.cpp
 * @author Mybrid Spalding
 * @date   Fri Apr 11 11:21:33 PDT 2008
 *
 * @brief  Testing the CopyNumberTableModel.cpp class
 *
 */


#include "translation/CopyNumberTableModel.h"
#include "translation/RunTimeEnvironment.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h" // includes Verbose.h
#include "util/Util.h"
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
 * @class CopyNumberTableTest
 * @brief cppunit class for testing CallSet behavior
 */
/*****************************************************************************/
class CopyNumberTableModelTest : public CppUnit::TestFixture
{

private:


  string m_defaultCopyNumberTableFileName;
  string m_missingHeaderLineFileName;
  string m_emptySampleNameDataFileName;
  string m_invalidAllele1DataFileName;
  string m_invalidAllele2DataFileName;
  string m_invalidProbeSetDataFileName;
  string m_invalidGeneDataFileName;

  string m_programName;
  string m_outputDir;

  CopyNumberTableModel *m_ttm;
  RunTimeEnvironment m_rte;

  CPPUNIT_TEST_SUITE(CopyNumberTableModelTest);
  CPPUNIT_TEST(constructor);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void constructor();
};
// end class CopyNumberTableModelTest
/*****************************************************************************/
/*****************************************************************************/
/**
 * Initialization method called before each test routine.
 *
 */
/*****************************************************************************/
void CopyNumberTableModelTest::setUp()
{

  cerr << endl;

  m_programName = "CopyNumberTableModeTest_setUp";

  /*
  m_defaultCopyNumberTableFileName =
    TEST_DATA_UNIT_DIR + "20080109_31set_DMET3_Genotypes_Short_qc_16various.txt";
  m_missingHeaderLineFileName =
    TEST_DATA_UNIT_DIR + "20080109_31set_DMET3_Genotypes_Short_qc_missing_header.txt";
  m_emptySampleNameDataFileName =
    TEST_DATA_UNIT_DIR + "20080109_31set_DMET3_Genotypes_Short_qc_empty_sample_data.txt";

  m_invalidAllele1DataFileName =
    TEST_DATA_UNIT_DIR + "20080109_31set_DMET3_Genotypes_Short_qc_invalid_allele1_data.txt";
  m_invalidAllele2DataFileName =
    TEST_DATA_UNIT_DIR + "20080109_31set_DMET3_Genotypes_Short_qc_invalid_allele2_data.txt";
  m_invalidProbeSetDataFileName =
    TEST_DATA_UNIT_DIR + "20080109_31set_DMET3_Genotypes_Short_qc_invalid_probeSet_data.txt";
  m_invalidGeneDataFileName =
    TEST_DATA_UNIT_DIR + "20080109_31set_DMET3_Genotypes_Short_qc_invalid_gene_data.txt";
  */

  m_outputDir = "output";

  m_rte.m_adtOpts.m_progName = m_programName;
  m_rte.m_adtOpts.m_outputDir = m_outputDir;

  m_rte.m_adtOpts.m_verbosity = ADT_VERBOSE_INPUT_FILES;

  m_rte.initializeRunTimeEnvironment();

  return;
}
// end CopyNumberTableModelTest::setUp
/*****************************************************************************/
/*****************************************************************************/
/**
 * CopyNumberTableModel::CopyNumberTableModel constructor test cases.
 *
 */
/*****************************************************************************/

void CopyNumberTableModelTest::constructor()
{
  cout << endl;
  Util::PrintTextFunctionTitle("CopyNumberTableModelTest", "constructor");
  Err::setThrowStatus(true);

  return;

  // Missing TsvFile HEADER line.

  cerr << "constructor: missing TsvFile header line..." << endl;;
  cerr << "FATAL ERROR expected" << endl;
  CPPUNIT_ASSERT_THROW(CopyNumberTableModel(m_rte,
                       m_missingHeaderLineFileName), Except);

  cerr << "ok" << endl;

  cerr << "constructor: invalid Gene data value..." << endl;
  cerr << "FATAL ERROR expected" << endl;
  CPPUNIT_ASSERT_THROW(CopyNumberTableModel(m_rte,
                       m_invalidGeneDataFileName), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid ProbeSet data value..." << endl;
  cerr << "FATAL ERROR expected" << endl;
  CPPUNIT_ASSERT_THROW(CopyNumberTableModel(m_rte,
                       m_invalidProbeSetDataFileName), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid Allele1 data value..." << endl;
  cerr << "FATAL ERROR expected" << endl;
  CPPUNIT_ASSERT_THROW(CopyNumberTableModel(m_rte,
                       m_invalidAllele1DataFileName), Except);
  cerr << "ok" << endl;

  cerr << "constructor: invalid Allele2 data value..." << endl;
  cerr << "FATAL ERROR expected" << endl;
  CPPUNIT_ASSERT_THROW(CopyNumberTableModel(m_rte,
                       m_invalidAllele2DataFileName), Except);
  cerr << "ok" << endl;

  cerr << "constructor: empty Sample Name data value..." << endl;
  cerr << "FATAL ERROR expected" << endl;
  CPPUNIT_ASSERT_THROW(CopyNumberTableModel(m_rte,
                       m_emptySampleNameDataFileName), Except);
  cerr << "ok" << endl;

  return;
}
// end CopyNumberTableModelTest::constructor
/*****************************************************************************/


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(CopyNumberTableModelTest);

////////////////////////////////////////////////////////////////
