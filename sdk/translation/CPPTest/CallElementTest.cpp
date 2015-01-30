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
 * @date   Mon Mar 17 12:46:31 PDT 2008
 *
 * @brief  Testing the CallElement object.
 *
 */


//
#include "translation/CallElement.h"
#include "translation/RunTimeEnvironment.h"
//
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

/**
 * @class CallElementTest
 * @brief cppunit class for testing call element behavior
 */
class CallElementTest : public CppUnit::TestFixture

{

private:

  map<string, CallElement*> ceTestCases;

  CPPUNIT_TEST_SUITE(CallElementTest);
  CPPUNIT_TEST(operator_equal);
  CPPUNIT_TEST(isValidXXXX);
  CPPUNIT_TEST(describeVerbose);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void operator_equal();
  void isValidXXXX();
  void describeVerbose();

};

void CallElementTest::setUp()
{


  ceTestCases.insert(pair<string, CallElement*>("equal1", new CallElement("366793", "G")));
  ceTestCases.insert(pair<string, CallElement*>("equal2", new CallElement("366793", "G")));
  ceTestCases.insert(pair<string, CallElement*>("equal3", new CallElement("366793", "T")));
  ceTestCases.insert(pair<string, CallElement*>("equal4", new CallElement("466793", "G")));
  ceTestCases.insert(pair<string, CallElement*>("wildcardNC", new CallElement("466793", "NC")));
  ceTestCases.insert(pair<string, CallElement*>("wildcardPRA", new CallElement("466793", "PRA")));

}

void CallElementTest::operator_equal()
{
  cout << endl;
  Util::PrintTextFunctionTitle("CallElementTest", "operator_equal");
  Err::setThrowStatus(true);
  CallElement *equal1 = ceTestCases["equal1"];
  CallElement *equal2 = ceTestCases["equal2"];
  CallElement *equal3 = ceTestCases["equal3"];
  CallElement *equal4 = ceTestCases["equal4"];
  CallElement *wildcardNC = ceTestCases["wildcardNC"];
  CallElement *wildcardPRA = ceTestCases["wildcardPRA"];

  // operator==, operator!=
  cerr << "a == b: a.m_probeSet == b.m_probeSet, a.m_base == b.m_base...";
  CPPUNIT_ASSERT((*equal1) == (*equal2));
  cerr << "ok" << endl;

  cerr << "a != b: a.m_probeSet == b.m_probeSet, a.m_base != b.m_base...";
  CPPUNIT_ASSERT((*equal1) != (*equal3));
  cerr << "ok" << endl;

  cerr << "a != b: a.m_probeSet != b.m_probeSet, a.m_base == b.m_base...";
  CPPUNIT_ASSERT((*equal1) != (*equal4));
  cerr << "ok" << endl;

  cerr << "a == b, b.m_base = NC: a.m_probeSet == b.m_probeSet, b.m_base = NC...";
  CPPUNIT_ASSERT((*equal4) == (*wildcardNC));
  cerr << "ok" << endl;

  cerr << "a == b, b.m_base = PRA: a.m_probeSet == b.m_probeSet, b.m_base = PRA...";
  CPPUNIT_ASSERT((*equal4) == (*wildcardPRA));
  cerr << "ok" << endl;

  cerr << "a != b, b.m_base = PRA: a.m_probeSet != b.m_probeSet, b.m_base = PRA...";
  CPPUNIT_ASSERT((*equal1) != (*wildcardPRA));
  cerr << "ok" << endl;


  cerr << endl;
  return;
}


void CallElementTest::isValidXXXX()
{

  cout << endl;
  Util::PrintTextFunctionTitle("CallElementTest", "isValidXXXX");
  Err::setThrowStatus(true);
  CallElement *valid1 = ceTestCases["equal1"];

  // isProbeSet
  cerr << "isValidProbeSet: valid1->isValidProbeSetd(valid1->m_probeSet)...";
  CPPUNIT_ASSERT(valid1->isValidProbeSet(valid1->m_probeSet));
  cerr << "ok" << endl;

  cerr << "isValidProbeSet: !valid1->isValidProbeSetd(\"badProbeSetId\")...";
  CPPUNIT_ASSERT(!valid1->isValidProbeSet("badProbeSetd"));
  cerr << "ok" << endl;

  cerr << "isValidProbeSet: !valid1->isValidProbeSetd(\"03939p\")...";
  CPPUNIT_ASSERT(!valid1->isValidProbeSet("03939p"));
  cerr << "ok" << endl;

  // isValidBaseId
  cerr << "isValidBaseId: valid1->isValidBase(valid1->m_base)...";
  CPPUNIT_ASSERT(valid1->isValidBase(valid1->m_bases[0]));
  cerr << "ok" << endl;

  cerr << "isValidBaseId: valid1->isValidBase(\"G\")...";
  CPPUNIT_ASSERT(valid1->isValidBase("G"));
  cerr << "ok" << endl;

  cerr << "isValidBaseId: !valid1->isValidBase(\"g\")...";
  CPPUNIT_ASSERT(!valid1->isValidBase("g"));
  cerr << "ok" << endl;

  cerr << "isValidBaseId: valid1->isValidBase(\"INS\")...";
  CPPUNIT_ASSERT(valid1->isValidBase("INS"));
  cerr << "ok" << endl;

  cerr << "isValidBaseId: valid1->isValidBase(\"DEL\")...";
  CPPUNIT_ASSERT(valid1->isValidBase("DEL"));
  cerr << "ok" << endl;

  cerr << "isValidBaseId: valid1->isValidBase(\"NC\")...";
  CPPUNIT_ASSERT(valid1->isValidBase("NC"));
  cerr << "ok" << endl;

  cerr << "isValidBaseId: valid1->isValidBase(\"PRA\")...";
  CPPUNIT_ASSERT(valid1->isValidBase("PRA"));
  cerr << "ok" << endl;

  cerr << endl;
  return;
}

void CallElementTest::describeVerbose()
{
  cout << endl;
  Util::PrintTextFunctionTitle("CallElementTest", "describeVerbose");
  Err::setThrowStatus(true);

  RunTimeEnvironment rte;

  CallElement *describe1 = ceTestCases["equal1"];

  rte.m_currentVerbosity = ADT_VERBOSE_INPUT_FILES ;

  Verbose::setLevel(ADT_VERBOSE_INPUT_FILES);

  cerr << "describeVerbose: describe1->describeVerbose(rte)..." << endl;
  CPPUNIT_ASSERT_NO_THROW(describe1->describeVerbose(rte));
  cerr << "ok" << endl;

  cerr << endl;
  return;
}

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(CallElementTest);

////////////////////////////////////////////////////////////////
