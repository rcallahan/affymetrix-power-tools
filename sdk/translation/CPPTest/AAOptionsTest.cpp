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
 * @file   translation/CPPTest/AAOptionsTest.cpp
 * @author Mybrid Spalding
 * @date   Fri Aug 29 12:44:43 PDT 2008
 * @brief  Test invalide console options.
 */

//

#include "translation/TranslationEngine.h"
//
#include "util/Err.h"        // includes Verbose.h
#include "util/CPPTest/Setup.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//

using namespace std;

/*****************************************************************************/
/**
 * @class AAOptionsTest
 * @brief cppunit class for testing invalid console options.
 */
/*****************************************************************************/
class AAOptionsTest : public CppUnit::TestFixture
{

private:

  map<string, string>        m_oneToOneOptions;
  vector< vector< string > > m_sampleTable;
  vector< string >           m_experimentListVector;

  RunTimeEnvironment rte;

  CPPUNIT_TEST_SUITE(AAOptionsTest);
  CPPUNIT_TEST(setUp);
  CPPUNIT_TEST(noOptions);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void noOptions();
};
// end AAOptionsTest
/*****************************************************************************/
/*****************************************************************************/
/**
 * Initialization method called before each test routine.
 */
/*****************************************************************************/
void AAOptionsTest::setUp()
{

  return;
}
// end AAOptionsTest::setUp
/*****************************************************************************/
/*****************************************************************************/
/**
 * OptionsTest::run
 *
 */
/*****************************************************************************/
void AAOptionsTest::noOptions()
{

  Util::PrintTextFunctionTitle("AAOptionsTest", "noOptions");

  Err::setThrowStatus(true);

  TranslationEngine te;

  te.m_controllerMask |= C_CMDLINE;

  NEGATIVE_TEST(te.run(), Except);
  cerr << "ok" << endl;

  CPPUNIT_ASSERT(te.m_rte.m_invalidOptionsMask & ADT_INVALID_OPTION_MASK_MISSING_TRANSLATION_TABLE);

  CPPUNIT_ASSERT(te.m_rte.m_invalidOptionsMask & ADT_INVALID_OPTION_MASK_MISSING_EXPERIMENT_CHP);

}
// end AAOptionsTest::noOptions
/*****************************************************************************/

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(AAOptionsTest);
