////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

/**
 * @file   PgOptionsTest.cpp
 * @author Chuck Sugnet
 * @date   Tue May  3 13:14:45 2005
 *
 * @brief  Testing program for the PgOptions module.
 *
 */

//
#include "util/AffxString.h"
#include "util/Convert.h"
#include "util/PgOptions.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <string>
//
#include "util/CPPTest/Setup.h"

using namespace std;

void define_pgoptionstest_options(PgOptions* opts)
{
  opts->setUsage("This is a test set of options.");

  PgOpt* o;

  opts->defineOption("i", "integer-opt",PgOpt::INT_OPT,
                     "Tests to see what an integer value looks like. "
                     "The description of this option should wrap.",
                     "10");

  opts->defineOption("d", "double-opt", PgOpt::DOUBLE_OPT,
                     "Tests to see what an double value looks like.",
                     "10.5");

  opts->defineOption("s1", "string-opt-1", PgOpt::STRING_OPT,
                     "Tests to see what an string value looks like.",
                     "stringOption-1");

  o=opts->defineOption("s2", "string-opt-2", PgOpt::STRING_OPT,
                       "Tests to see what an string value looks like.",
                       "stringOption-2");
  o->allowMutipleValues(1);

  opts->defineOption("b", "bool-opt", PgOpt::BOOL_OPT,
                     "Tests to see what an bool value looks like.",
                     "1");
};



/**
 * @class PgOptionsTest
 * @brief Cppunit class for testing PgOptions.
 */
class PgOptionsTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( PgOptionsTest );
  //
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testConstructorCopy );
  CPPUNIT_TEST( testDefineOption );
  CPPUNIT_TEST( testMultiples );
  CPPUNIT_TEST( testParsing );
  CPPUNIT_TEST( testNegativeParsing );
  CPPUNIT_TEST( testResetToDefaults );
  CPPUNIT_TEST( testFlag1 );
  CPPUNIT_TEST( testClear );
  CPPUNIT_TEST( testClearOption );
  CPPUNIT_TEST( testAppendOptions );
  CPPUNIT_TEST( testOperatorAssign );
  CPPUNIT_TEST( testAddPgOpt );
  CPPUNIT_TEST( testSet );
  CPPUNIT_TEST( testPush );
  CPPUNIT_TEST( testIsOptDefined );
  //*********vliber**********
  CPPUNIT_TEST( testFunctions1);
  CPPUNIT_TEST( testFunctions2);
  CPPUNIT_TEST( testConstructor2 );
  CPPUNIT_TEST( testFunctions3);
  CPPUNIT_TEST( testmatchOneArg);
  //
  CPPUNIT_TEST_SUITE_END();

public:
  PgOptions *opts;

  /// @brief Setup before each unit test.
  void setUp();

  /// @brief dispose of the setup.
  void tearDown();

  /// @brief Test the construction of arguments.
  void testConstructor();

  /// @brief Test the copy constructor (rsatin 04may09)
  void testConstructorCopy();

  /// @brief Test the parsing of arguments.
  void testParsing();

  /// @brief Test negative cases for invalid arguments (rsatin 05may09)
  void testNegativeParsing();

  /// @brief Test reset to defaults.
  void testResetToDefaults();

  /// @brief Test defineOption()
  void testDefineOption();
  
  /// @brief Test multiple args working
  void testMultiples();

  /// @brief Test flags being set correctly
  void testFlag1();

  /// @brief test whole-object clear method (rsatin 04may09)
  void testClear();

  /// @brief test single option clear method (rsatin 04may09)
  void testClearOption();

  /// @brief test appendOptions method (rsatin 05may09)
  void testAppendOptions();

  /// @brief test operator= method (rsatin 05may09)
  void testOperatorAssign();

  /// @brief test addPgOpt method (rsatin 04may09)
  void testAddPgOpt();

  /// @brief test set method (rsatin 05may09)
  void testSet();

  /// @brief test push method (rsatin 05may09)
  void testPush();

  /// @brief test isOptDefined method (rsatin 05may09)
  void testIsOptDefined();

  //vliber
  void testFunctions1();
  void testFunctions2();
  void testConstructor2();
  void testFunctions3();
  void testmatchOneArg();
};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( PgOptionsTest );

/// @brief     set opts to a known set of options.
void PgOptionsTest::setUp() {
  //printf("PgOptionsTest::setup()\n");
  opts = new PgOptions;
  define_pgoptionstest_options(opts);
}

/// @brief dispose of the fixture.
void PgOptionsTest::tearDown() {
  //printf("PgOptionsTest::teardown()\n");
  delete opts;
  opts=NULL;
}

/// @brief Make sure that our options were initialized correctly.
void PgOptionsTest::testConstructor() {
  Verbose::out(1, "");
  Verbose::out(1, "***PgOptionsTest testcases***");
  Verbose::out(1, "PgOptionsTest::testConstructor");
 //printf("testConstructor\n");
  PgOpt *o;
  
  //
  //CPPUNIT_ASSERT(o==NULL);
  o=opts->findOpt("integer-opt");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="10");
  CPPUNIT_ASSERT(o->m_longName==string("integer-opt"));//vliber
  POSITIVE_TEST(o->checkParseIsOk(o->getDefaultValue()));//vliber
  
  o=opts->findOpt("double-opt");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="10.5");
  CPPUNIT_ASSERT(o->getValue()=="10.5");//vliber
  POSITIVE_TEST(o->checkParseIsOk(o->getDefaultValue()));//vliber
    
  o=opts->findOpt("s1");
  CPPUNIT_ASSERT(o->m_type==PgOpt::STRING_OPT);//vliber
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="stringOption-1");
  POSITIVE_TEST(o->checkParseIsOk(o->getDefaultValue()));//vliber
  //Goal: FATAL ERROR: Option 'string-opt-1' already defined.
  NEGATIVE_TEST(opts->defineOption("s1", "string-opt-1", PgOpt::STRING_OPT, "string value.", "stringOption"),Except);//vliber

  o=opts->findOpt("b");
  CPPUNIT_ASSERT(o->m_type==PgOpt::BOOL_OPT);//vliber
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="1");
  CPPUNIT_ASSERT(o->getValue()=="1");//vliber
  // should find the same option under both names.
  CPPUNIT_ASSERT(opts->findOpt("b")==opts->findOpt("bool-opt")); 
  POSITIVE_TEST(o->checkParseIsOk(o->getDefaultValue()));//vliber
  //opts->usage();//todo vliber  test usage
}

/// @brief test copy constructure PgOptions(const PgOptions &options) (rsatin 04may09)
void PgOptionsTest::testConstructorCopy() {
  Verbose::out(1, "PgOptionsTest::testConstructorCopy");
  PgOptions *opts2 = new PgOptions(*opts);
  PgOpt *o;
  //
  o=opts2->findOpt("integer-opt");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="10");
  o=opts2->findOpt("double-opt");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="10.5");
  o=opts2->findOpt("s1");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="stringOption-1");
  o=opts2->findOpt("b");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="1");
  CPPUNIT_ASSERT(o->getValue()=="1");//vliber
  // should find the same option under both names.
  CPPUNIT_ASSERT(opts2->findOpt("b")==opts2->findOpt("bool-opt"));
}


/// @brief are the command line args parsed correctly
void PgOptionsTest::testParsing() {
  Verbose::out(1, "PgOptionsTest::testParsingy");
  //printf("testParsing\n");

  const char* argv[] = {"myProg",
                  "-i", "10",
                  "-b",
                  "hello", "world",
                  "--double-opt", "-10.5",
                  "-s1", "somestring",
                  NULL};
 
  opts->resetToDefaults();
  opts->parseArgv(argv);
  //opts->dump();
  //
  CPPUNIT_ASSERT( opts->getBool("b") == true );
  CPPUNIT_ASSERT( opts->getInt("integer-opt") == 10 );
  CPPUNIT_ASSERT( Convert::doubleCloseEnough(opts->getDouble("d"), -10.5));
  CPPUNIT_ASSERT( opts->get("string-opt-1")=="somestring");
  CPPUNIT_ASSERT( opts->get("s2") == "stringOption-2");
  //
  CPPUNIT_ASSERT( opts->getArgCount() == 2 );
  CPPUNIT_ASSERT( opts->getArg(0)=="hello");
  CPPUNIT_ASSERT( opts->getArg(1)=="world");
}

/// @brief negative test cases for invalid command line args (rsatin 05may09)
void PgOptionsTest::testNegativeParsing() {
 Verbose::out(1, "PgOptionsTest::testNegativeParsing");
  //printf("testNegativeParsing\n");

  const char* argv1[] = {"myProg",
                   "--=10",
                   "hello", "world",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv1), Except );             // argument "=..." 

  const char* argv2[] = {"myProg",
                   "--no-integer-opt",
                   "hello", "world",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv2), Except );             // '--no-' with non-boolean option

  const char* argv3[] = {"myProg",
                   "--alkjskfjkd-opt",
                   "hello", "world",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv3), Except );             // unrecognized option

  
/*BUG rsatin 05may09 ==> type-specific value format errors are not detected (validation missing or not functional)
  const char* argv4[] = {"myProg",
                   "--bool-opt=kjhjahd",
                   "hello", "world",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv4), Except );             // boolean option with invalid value 

  const char* argv5[] = {"myProg",
                   "--integer-opt=kjhjahd",
                   "hello", "world",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv5), Except );             // integer option with invalid value 

  const char* argv6[] = {"myProg",
                   "--double-opt=kjhjahd",
                   "hello", "world",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv6), Except );             // double option with invalid value 
*/
  PgOpt* o;

  const char* argv7[] = {"myProg",
                   "--bool-opt=true",
                   "--bool-opt=false",
                   "hello", "world",
                   NULL};

  CPPUNIT_ASSERT( o = opts->findOpt("bool-opt") );
  o->allowMutipleValues(false);
  POSITIVE_TEST( opts->parseArgv(argv7) );                     // multiple values for single valued option
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="false" );

  const char* argv8[] = {"myProg",
                   "--integer-opt=13",
                   "--integer-opt=29",
                   "hello", "world",
                   NULL};

  CPPUNIT_ASSERT( o = opts->findOpt("integer-opt") );
  o->allowMutipleValues(false);
  POSITIVE_TEST( opts->parseArgv(argv8) );                     // multiple values for single valued option
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="29" );

  const char* argv9[] = {"myProg",
                   "--string-opt-1=first",
                   "--string-opt-1=second",
                   "hello", "world",
                   NULL};

  CPPUNIT_ASSERT( o = opts->findOpt("string-opt-1") );
  o->allowMutipleValues(false);
  POSITIVE_TEST( opts->parseArgv(argv9) );                     // multiple values for single valued option
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="second" );

  opts->resetToDefaults();
  opts->defineOption("", "xml-file",PgOpt::STRING_OPT,
                     "Path to XML file passed to PgOptions::setOptionsFromXMLFile.",
                     "");
  opts->defineOption("", "xml-file-append-only",PgOpt::STRING_OPT,
                     "Path to append XML file passed to PgOptions::setOptionsFromXMLFile.",
                     "");
  const char* argv10[] = {"myProg",
                   "--xml-file=input/valid_sample.xml",
                   "--xml-file-append-only=input/valid_sample.xml",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv10), Except );             // two xml files. 

  const char* argv11[] = {"myProg",
                   "--xml-file=input/valid_sample.xml",
                   "--xml-file=input/valid_sample.xml",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv11), Except );             // two xml files. 

  const char* argv12[] = {"myProg",
                   "--xml-file-append-only=input/valid_sample.xml",
                   "--xml-file-append-only=input/valid_sample.xml",
                   NULL};
  
  NEGATIVE_TEST( opts->parseArgv(argv12), Except );             // two xml files. 
  
}

/**
 * Check resetToDefaults does so.
 */
void PgOptionsTest::testResetToDefaults() {
  Verbose::out(1, "PgOptionsTest::testResetToDefaults");
  //printf("testResetToDefaults\n");

  PgOpt* o;

  // set them to non-default values.
  o = opts->findOpt("integer-opt");
  o->setValue("123456");
  o = opts->findOpt("string-opt-1");
  o->setValue("not the default value");

  // boink!
  opts->resetToDefaults();

  // check
  o = opts->findOpt("integer-opt");
  CPPUNIT_ASSERT(o->getValue(0)==o->getDefaultValue());
  o = opts->findOpt("string-opt-1");
  CPPUNIT_ASSERT(o->getValue(0)==o->getDefaultValue());
}

/**
 * Check that defineOption works.
 * This uses the "defineOption" method to dynamicly allocate options.
 */
void PgOptionsTest::testDefineOption() {
  cout<<endl;
  Verbose::out(1, "PgOptionsTest::testDefineOption");
  //printf("testDefineOption\n");

  PgOptions* my_opts=new PgOptions;
  my_opts->setUsage("Usage info.");

  my_opts->defineOption("-h","--help",PgOpt::BOOL_OPT,
                        "A help message.",
                        "false");
  my_opts->defineOption("-d","--double",PgOpt::DOUBLE_OPT,
                        "A float.",
                        "0.0");
  my_opts->defineOption("-s","--string",PgOpt::STRING_OPT,
                        "A string.",
                        "the quick brown fox.");
  my_opts->defineOption("-int","--integer",PgOpt::INT_OPT,
                        "A integer.",
                        "0");
  // negative cases: duplicate names (rsatin 04may09)
  NEGATIVE_TEST( my_opts->defineOption("-h","--help",PgOpt::INT_OPT,"h help","70"), Except );
  NEGATIVE_TEST( my_opts->defineOption("dX","--double",PgOpt::DOUBLE_OPT,"d help","3.1415"), Except );
//BUG rsatin 05may09 ==> does not detect option with duplicate short name
//NEGATIVE_TEST( my_opts->defineOption("-d","--doubleX",PgOpt::DOUBLE_OPT,"d help","2.7184"), Except );
  // negative cases: invalid value for numeric value types (rsatin 04may09)
  NEGATIVE_TEST( my_opts->defineOption("i2","i-name2",PgOpt::INT_OPT,"i2 help","j3u3h"), Except );
  NEGATIVE_TEST( my_opts->defineOption("d2","d-name2",PgOpt::DOUBLE_OPT,"d2 help","2ihg4"), Except );
//BUG rsatin 05may09 ==> PgOpt::checkParseIsOk does not call Convert::toBoolCheck to validate BOOL_OPT values
//NEGATIVE_TEST( my_opts->defineOption("b2","b-name2",PgOpt::BOOL_OPT,"b2 help","2kje4"), Except );
  //
  delete my_opts;
}

/// @brief     Are mutiple args working?
void PgOptionsTest::testMultiples() {
  Verbose::out(1, "PgOptionsTest::testMultiples");
  PgOptions* my_opts;
  PgOpt* opt;

  my_opts=new PgOptions;
  my_opts->setUsage("Usage info.");
  opt=my_opts->defineOption("m","mult",PgOpt::STRING_OPT,
                            "Multiples",
                            "");
  opt->allowMutipleValues(1);
  // my_opts->dump();

  const char* argv[]={
    "progname",
    "-m","a0",
    "-m","b1",
    "-m","c2",
    "--mult","d3",
    "--mult","e4",
    "--mult","f5",
    "-m=equal", // 6
    "-m=", // 7
    "-m==", // 8
    "-m=9",
    "hello","world",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv);
  //
  //my_opts->dump();
  //
  CPPUNIT_ASSERT(opt->getValueCount()==10);
  CPPUNIT_ASSERT(opt->getValue()=="a0");
  CPPUNIT_ASSERT(opt->getValue(0)=="a0");
  CPPUNIT_ASSERT(opt->getValue(1)=="b1");
  CPPUNIT_ASSERT(opt->getValue(5)=="f5");
  CPPUNIT_ASSERT(opt->getValue(6)=="equal");
  CPPUNIT_ASSERT(opt->getValue(7)=="");
  CPPUNIT_ASSERT(opt->getValue(8)=="=");
  //
  delete my_opts;
}

/// @brief  Are the flags being set correctly?
void PgOptionsTest::testFlag1() {
  Verbose::out(1, "PgOptionsTest::testFlag1");
  PgOptions* my_opts;
  PgOpt* opt;

  my_opts=new PgOptions;
  my_opts->setUsage("Usage info.");
  opt=my_opts->defineOption("","flag",PgOpt::BOOL_OPT,
                            "","false");
  opt=my_opts->defineOption("","print",PgOpt::BOOL_OPT,
                            "","true");

  const char* argv1[]={
    "progname",
    "-flag",
    "last",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv1);
  //my_opts->dump();
  //
  CPPUNIT_ASSERT(my_opts->getBool("flag")==true);
  CPPUNIT_ASSERT(my_opts->getArgCount()==1);
  CPPUNIT_ASSERT(my_opts->getArg(0)=="last");

  const char* argv2[]={
    "progname",
    "-flag","1",
    "last",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv2);
  //my_opts->dump();
  //
  CPPUNIT_ASSERT(my_opts->getBool("flag")==true);
  CPPUNIT_ASSERT(my_opts->getArgCount()==1);
  CPPUNIT_ASSERT(my_opts->getArg(0)=="last");

  const char* argv3[]={
    "progname",
    "-flag","false",
    "last",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv3);
  //my_opts->dump();
  //
  CPPUNIT_ASSERT(my_opts->getBool("flag")==false);
  CPPUNIT_ASSERT(my_opts->getArgCount()==1);
  CPPUNIT_ASSERT(my_opts->getArg(0)=="last");

  const char* argv4[]={
    "progname",
    "-flag=0",
    "last",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv4);
  //my_opts->dump();
  //
  CPPUNIT_ASSERT(my_opts->getBool("flag")==false);
  CPPUNIT_ASSERT(my_opts->getArgCount()==1);
  CPPUNIT_ASSERT(my_opts->getArg(0)=="last");

  const char* argv5[]={
    "progname",
    "last",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv5);
  // my_opts->dump();
  //
  CPPUNIT_ASSERT(my_opts->getBool("print")==true);
  CPPUNIT_ASSERT(my_opts->getArgCount()==1);
  CPPUNIT_ASSERT(my_opts->getArg(0)=="last");

  const char* argv6[]={
    "progname",
    "-no-print",
    "last",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv6);
  // my_opts->dump();
  //
  CPPUNIT_ASSERT(my_opts->getBool("print")==false);
  CPPUNIT_ASSERT(my_opts->getArgCount()==1);
  CPPUNIT_ASSERT(my_opts->getArg(0)=="last");

  const char* argv7[]={
    "progname",
    "--",
    "--last",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv7);
  // my_opts->dump();
  //
  CPPUNIT_ASSERT(my_opts->getBool("flag")==false);
  CPPUNIT_ASSERT(my_opts->getBool("print")==true);
  CPPUNIT_ASSERT(my_opts->getArgCount()==0);

  const char* argv8[]={
    "progname",
    "--no-flag",
    "true",
    NULL};
  my_opts->resetToDefaults();
  my_opts->parseArgv(argv8);
  // my_opts->dump();
  //
  CPPUNIT_ASSERT(my_opts->getBool("flag")==false);
  CPPUNIT_ASSERT(my_opts->getBool("print")==true);
  CPPUNIT_ASSERT(my_opts->getArgCount()==1);
  CPPUNIT_ASSERT(my_opts->getArg(0)=="true");

  //
  delete my_opts;
}

/// @brief test whole-object clear method (rsatin 04may09)
void PgOptionsTest::testClear() {
  Verbose::out(1, "PgOptionsTest::testClear");
  POSITIVE_TEST( opts->clear() );
  CPPUNIT_ASSERT( opts->m_progName == "" );
  CPPUNIT_ASSERT( opts->m_usageMsg == "" );
  CPPUNIT_ASSERT( opts->m_strXMLParameterFileName == "" );
  CPPUNIT_ASSERT( opts->m_strXMLParameterFileGuid == "" );
  CPPUNIT_ASSERT( opts->m_argv.empty() );
  CPPUNIT_ASSERT( opts->m_args.empty() );
  CPPUNIT_ASSERT( opts->m_option_vec.empty() );
  CPPUNIT_ASSERT( opts->m_option_section.empty() );
  CPPUNIT_ASSERT( opts->m_option_map.empty() );
}

/// @brief test single option clear method (rsatin 04may09)
void PgOptionsTest::testClearOption() {
  Verbose::out(1, "PgOptionsTest::testClearOption");
  PgOpt *o;
  // setup option with multiple values
  o = opts->findOpt("integer-opt");
  CPPUNIT_ASSERT( o != NULL );
  o->allowMutipleValues(true);
  o->pushValue("25");
  o->pushValue("50");
  CPPUNIT_ASSERT( o->getValueCount() >= 2 );
  CPPUNIT_ASSERT( o->getDefaultValue()=="10" );
  o = NULL;
  // positive case: clear values
  POSITIVE_TEST( opts->clear("integer-opt") );                    // clear option with non-empty values list
  o = opts->findOpt("integer-opt");
  CPPUNIT_ASSERT( o != NULL );
  CPPUNIT_ASSERT( o->getValueCount() == 0 );
  CPPUNIT_ASSERT( o->getDefaultValue()=="10" );
  o = NULL;
  // positive case: clear values
  POSITIVE_TEST( opts->clear("integer-opt") );                    // clear option with empty values list
  o = opts->findOpt("integer-opt");          
  CPPUNIT_ASSERT( o != NULL );
  CPPUNIT_ASSERT( o->getValueCount() == 0 );
  CPPUNIT_ASSERT( o->getDefaultValue()=="10" );
  o = NULL;
  // negative case: clear values for non-existent option
  NEGATIVE_TEST( opts->clear("kjhsdfkjsh"), Except );             // clear non-existent option
}

/// @brief test appendOptions method (rsatin 05may09)
void PgOptionsTest::testAppendOptions() {
  Verbose::out(1, "PgOptionsTest::testAppendOptions");
  PgOpt *o;
  // initialize non-empty initial collection of options
  PgOptions *opts2 = new PgOptions;
  opts2->defineOption("d0", "double-opt-0", PgOpt::DOUBLE_OPT,
                      "Initial double option.",
                      "50.1");
  opts2->defineOption("s0", "string-opt-0", PgOpt::STRING_OPT,
                      "Initial string option.",
                      "stringOption-0");
  // append more options to non-empty collection of options
  size_t opts_end_appended = opts->m_option_vec.size();
  size_t opts2_end_initial = opts2->m_option_vec.size();
  POSITIVE_TEST( opts2->appendOptions(*opts) );                   // test appendOptions method
  size_t opts2_end_final = opts2->m_option_vec.size();
  // validate result is union of initial and appended contents
  CPPUNIT_ASSERT( opts2_end_final == opts2_end_initial+opts_end_appended );
  o=opts2->findOpt("double-opt-0");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="50.1");
  o=opts2->findOpt("s0");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="stringOption-0");
  o=opts2->findOpt("integer-opt");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="10");
  o=opts2->findOpt("double-opt");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="10.5");
  o=opts2->findOpt("s1");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="stringOption-1");
  o=opts2->findOpt("b");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="1");
  // validate appending empty options list has no effect
  PgOptions *opts_empty = new PgOptions;
  POSITIVE_TEST( opts2->appendOptions(*opts_empty) );             // test appendOptions method
  opts2_end_final = opts2->m_option_vec.size();
  CPPUNIT_ASSERT( opts2_end_final == opts2_end_initial+opts_end_appended );
  delete opts_empty;
  delete opts2;
}

/// @brief test operator= method (rsatin 05may09)
void PgOptionsTest::testOperatorAssign() {
   Verbose::out(1, "PgOptionsTest::testOperatorAssign");
  PgOpt *o;
  // initialize non-empty initial collection of options
  PgOptions *opts2 = new PgOptions;
  opts2->defineOption("d0", "double-opt-0", PgOpt::DOUBLE_OPT,
                      "Initial double option.",
                      "50.1");
  opts2->defineOption("s0", "string-opt-0", PgOpt::STRING_OPT,
                      "Initial string option.",
                      "string value 0");
  // append more options to non-empty collection of options
  size_t opts_end_assigned = opts->m_option_vec.size();
  POSITIVE_TEST(*opts2 = *opts );                             // test operator= method
  size_t opts2_end_final = opts2->m_option_vec.size();
  // validate contains assigned options only 
  CPPUNIT_ASSERT( opts2_end_final == opts_end_assigned );
  CPPUNIT_ASSERT( opts2->findOpt("d0") == NULL );
  CPPUNIT_ASSERT( opts2->findOpt("double-opt-0") == NULL );
  CPPUNIT_ASSERT( opts2->findOpt("s0") == NULL );
  CPPUNIT_ASSERT( opts2->findOpt("string-opt-0") == NULL );
  o=opts2->findOpt("integer-opt");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="10");
  o=opts2->findOpt("double-opt");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="10.5");
  o=opts2->findOpt("s1");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="stringOption-1");
  o=opts2->findOpt("b");
  CPPUNIT_ASSERT(o!=NULL);
  CPPUNIT_ASSERT(o->getDefaultValue()=="1");
  // validate assigned empty options list clears collection of options
  PgOptions *opts_empty = new PgOptions;
  POSITIVE_TEST( *opts2 = *opts_empty );                      // test operator= method
  opts2_end_final = opts2->m_option_vec.size();
  CPPUNIT_ASSERT( opts2_end_final == 0 );
  delete opts_empty;
  delete opts2;
}

/// @brief test addPgOpt method (rsatin 05may09)
void PgOptionsTest::testAddPgOpt() {
   Verbose::out(1, "PgOptionsTest::testAddPgOpt");
  // initialize new option with multiple values
  PgOpt *o_orig = new PgOpt;
  o_orig->m_shortName = "ia";
  o_orig->m_longName = "integer-added"; 
  o_orig->m_help = "ia help text"; 
  o_orig->m_type = PgOpt::INT_OPT;
  o_orig->m_allowMultiple = true;
  o_orig->m_defaultValue = "35"; 
  o_orig->pushValue("10");
  o_orig->pushValue("20");
  // positive case: add new option  
  PgOpt *o_copy = NULL;
  POSITIVE_TEST( o_copy = opts->addPgOpt(o_orig) );             // test addPgOpt method
  CPPUNIT_ASSERT( o_orig != o_copy );
  delete o_orig;                                                // insure not referencing original object
  // validate newly added option content
  o_copy = opts->findOpt("integer-added");
  CPPUNIT_ASSERT( o_copy != NULL );
  CPPUNIT_ASSERT( o_copy->m_shortName == "ia" );
  CPPUNIT_ASSERT( o_copy->m_longName == "integer-added" );
  CPPUNIT_ASSERT( o_copy->m_help == "ia help text" );
  CPPUNIT_ASSERT( o_copy->m_type == PgOpt::INT_OPT );
  CPPUNIT_ASSERT( o_copy->m_allowMultiple );
  CPPUNIT_ASSERT( o_copy->getDefaultValue() == "35" );
  CPPUNIT_ASSERT( o_copy->getValueCount() == 2 );
  CPPUNIT_ASSERT( o_copy->getValue(0) == "10" );
  CPPUNIT_ASSERT( o_copy->getValue(1) == "20" );
}


/// @brief test set method (rsatin 05may09)
void PgOptionsTest::testSet() {
  Verbose::out(1, "PgOptionsTest::testSet");

  PgOpt* o;

  // setup for test case
  opts->defineOption("x", "xml-file",PgOpt::STRING_OPT,
                     "Path to XML file passed to PgOptions::setOptionsFromXMLFile.",
                     "10");

  // positive cases: set options to valid values 
  POSITIVE_TEST( opts->set("bool-opt", "true" ) );                  // set boolean single-valued option
  CPPUNIT_ASSERT( o = opts->findOpt("bool-opt") );
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="true" );
  POSITIVE_TEST( opts->set("integer-opt", "23" ) );                 // set integer single-valued option
  CPPUNIT_ASSERT( o = opts->findOpt("integer-opt") );
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="23" );
  POSITIVE_TEST( opts->set("double-opt", "76.3" ) );                // set double single-valued option
  CPPUNIT_ASSERT( o = opts->findOpt("double-opt") );
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="76.3" );
  POSITIVE_TEST( opts->set("string-opt-1", "string1") );            // set string single-valued option
  CPPUNIT_ASSERT( o = opts->findOpt("string-opt-1") );
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="string1" );

  // positive cases: set multi-valued option to new single value
  o = opts->findOpt("string-opt-2");                                // initialize option with multiple values
  o->allowMutipleValues(true);
  o->clearValues();
  o->pushValue("strA");
  o->pushValue("strB");
  opts->set("string-opt-2", "str-set-1" );                          // set option with multiple values to single-value
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="str-set-1" );
  opts->set("string-opt-2", "str-set-2" );                          // set again, insure option stays single-value 
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="str-set-2" );

  // negative cases: invalid option name value pairs 
  NEGATIVE_TEST( opts->set("dkjhdskjs", "value" ), Except ); 
//NEGATIVE_TEST( opts->set("bool-opt", "dysj3q" ), Except );        //BUG rsatin 06may09 boolean value not validated
//NEGATIVE_TEST( opts->set("integer-opt", "asdfdsf" ), Except );    //BUG rsatin 06may09 boolean value not validated
//NEGATIVE_TEST( opts->set("double-opt", "dafdsa" ), Except );      //BUG rsatin 06may09 boolean value not validated

  // positive special case: set path for setOptionsFromXMLFile (first attempt)     // TO DO: need XML options file to load
//POSITIVE_TEST( opts->set("xml-file", "c:\\temp\\XMLFile.xml") );                 // set xml-file option first time

  // negative special case: set path for setOptionsFromXMLFile (repeat attempts)
//NEGATIVE_TEST( opts->set("xml-file", "c:\\temp\\filename.xml"), Except );        // set xml-file option repeat attempt
}


/// @brief test push method (rsatin 05may09)
void PgOptionsTest::testPush() {
  Verbose::out(1, "PgOptionsTest::testPush");

  PgOpt* o;

  // positive cases: push single valid values for each option type
  POSITIVE_TEST( opts->push("bool-opt", "true" ) );                    // push boolean single-valued option value
  CPPUNIT_ASSERT( o = opts->findOpt("bool-opt") );
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="true" );
  POSITIVE_TEST( opts->push("integer-opt", "23" ) );                   // push integer single-valued option value
  CPPUNIT_ASSERT( o = opts->findOpt("integer-opt") );
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="23" );
  POSITIVE_TEST( opts->push("double-opt", "76.3" ) );                  // push double single-valued option value
  CPPUNIT_ASSERT( o = opts->findOpt("double-opt") );
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="76.3" );
  POSITIVE_TEST( opts->push("string-opt-1", "string1") );              // push string single-valued option value
  CPPUNIT_ASSERT( o = opts->findOpt("string-opt-1") );
  CPPUNIT_ASSERT( o->getValueCount()==1 );
  CPPUNIT_ASSERT( o->getValue()=="string1" );

  // positive cases: push another valid values for single-value only options
  POSITIVE_TEST( opts->push("bool-opt", "false" ) );                   // push boolean single-valued option value (again)
  CPPUNIT_ASSERT( o = opts->findOpt("bool-opt") );
//CPPUNIT_ASSERT( o->getValueCount()==1 );                             //BUG rsatin 06may09 second value pushed onto list
//CPPUNIT_ASSERT( o->getValue()=="false" );
  POSITIVE_TEST( opts->push("integer-opt", "73" ) );                   // push boolean single-valued option value (again)
  CPPUNIT_ASSERT( o = opts->findOpt("integer-opt") );
//CPPUNIT_ASSERT( o->getValueCount()==1 );                             //BUG rsatin 06may09 second value pushed onto list
//CPPUNIT_ASSERT( o->getValue()=="73" );
  POSITIVE_TEST( opts->push("double-opt", "3.1415" ) );                // push boolean single-valued option value (again)
  CPPUNIT_ASSERT( o = opts->findOpt("double-opt") );
//CPPUNIT_ASSERT( o->getValueCount()==1 );                             //BUG rsatin 06may09 second value pushed onto list
//CPPUNIT_ASSERT( o->getValue()=="3.1415" );
  POSITIVE_TEST( opts->push("string-opt-1", "string2" ) );             // push boolean single-valued option value (again)
  CPPUNIT_ASSERT( o = opts->findOpt("string-opt-1") );
//CPPUNIT_ASSERT( o->getValueCount()==1 );                             //BUG rsatin 06may09 second value pushed onto list
//CPPUNIT_ASSERT( o->getValue()=="string2" );

  // positive cases: push multiple valid values for multiple-value enabled option
  POSITIVE_TEST( opts->push("string-opt-2", "stringA") );   
  POSITIVE_TEST( opts->push("string-opt-2", "stringB") );   
  POSITIVE_TEST( opts->push("string-opt-2", "stringC") );   
  CPPUNIT_ASSERT( o = opts->findOpt("string-opt-2") );
  CPPUNIT_ASSERT( o->getValueCount()==3 );
  CPPUNIT_ASSERT( o->getValue(0)=="stringA" );
  CPPUNIT_ASSERT( o->getValue(1)=="stringB" );
  CPPUNIT_ASSERT( o->getValue(2)=="stringC" );

  // negative cases: invalid option name value pairs 
  NEGATIVE_TEST( opts->push("dkjhdskjs", "value"), Except ); 
//NEGATIVE_TEST( opts->push("bool-opt", "dysj3q" ), Except );          //BUG rsatin 06may09 boolean value not validated
//NEGATIVE_TEST( opts->push("integer-opt", "asdfdsf" ), Except );      //BUG rsatin 06may09 boolean value not validated
//NEGATIVE_TEST( opts->push("double-opt", "dafdsa" ), Except );        //BUG rsatin 06may09 boolean value not validated
}

/// @brief test isOptDefined method (rsatin 05may09)
void PgOptionsTest::testIsOptDefined() {
   Verbose::out(1, "PgOptionsTest::testIsOptDefined");
  // positive cases
  CPPUNIT_ASSERT( opts->isOptDefined("b") );
  CPPUNIT_ASSERT( opts->isOptDefined("bool-opt") );
  CPPUNIT_ASSERT( opts->isOptDefined("i") );
  CPPUNIT_ASSERT( opts->isOptDefined("integer-opt") );
  CPPUNIT_ASSERT( opts->isOptDefined("d") );
  CPPUNIT_ASSERT( opts->isOptDefined("double-opt") );
  CPPUNIT_ASSERT( opts->isOptDefined("s1") );
  // negative cases
  CPPUNIT_ASSERT( opts->isOptDefined("string-opt-1") );
  CPPUNIT_ASSERT( !opts->isOptDefined("didshiuhef") );
  CPPUNIT_ASSERT( !opts->isOptDefined("/n/r/t") );
  CPPUNIT_ASSERT( !opts->isOptDefined("-b") );
  CPPUNIT_ASSERT( !opts->isOptDefined("*") );
  CPPUNIT_ASSERT( !opts->isOptDefined("") );
}

// test function vliber
void PgOptionsTest::testFunctions1() {

   Verbose::out(1, "PgOptionsTest::testFunctions1");
  PgOptions* opts;
  PgOpt* o;

  opts=new PgOptions;
  CPPUNIT_ASSERT(opts->getUsage()== string(""));
  opts->setUsage("Usage info.");
  CPPUNIT_ASSERT(opts->getUsage()== string("Usage info."));
  o=opts->defineOption("i", "integer-opt",PgOpt::INT_OPT,
                     "Tests to see what an integer value looks like. "
                     "The description of this option should wrap.",
                     "10");
  CPPUNIT_ASSERT(o->getDefaultValue()=="10");
  CPPUNIT_ASSERT(o->getValue()=="10");
  POSITIVE_TEST(o->getValue());
  CPPUNIT_ASSERT(o->m_shortName==string("i"));
  CPPUNIT_ASSERT(o->m_help==string("Tests to see what an integer value looks like. The description of this option should wrap."));
  CPPUNIT_ASSERT(o->m_type==PgOpt::INT_OPT);
  CPPUNIT_ASSERT(o->m_allowMultiple==0);
  CPPUNIT_ASSERT_NO_THROW(o->checkParseIsOk(o->getDefaultValue()));
  //Goal: FATAL ERROR: Don't recognize option with name: 'integer'.
  NEGATIVE_TEST(opts->mustFindOpt("integer"),Except);
  POSITIVE_TEST(opts->mustFindOpt("integer-opt"));
  POSITIVE_TEST(opts->bind(string(""),o));
  CPPUNIT_ASSERT(opts->mustFindOpt("integer-opt")->getDefaultValue()==opts->mustFindOpt("i")->getDefaultValue());

  CPPUNIT_ASSERT(o->isSet()==false);
  o->allowMutipleValues(1);
  CPPUNIT_ASSERT(o->m_allowMultiple==1);
  o->pushValue(string("one"));
  CPPUNIT_ASSERT(o->getValueCount()==1);
  o->pushValue(string("two"));
  CPPUNIT_ASSERT(o->getValueCount()==2);
  o->allowMutipleValues(0);
  CPPUNIT_ASSERT(o->m_allowMultiple==0);
  // ???? todo vliber Allow to push for option with multiple=0 Why? Is it correct?
  //May be it is OK, because app. verifyes multiple in the matchOneArg function before add to the vector.
  //look in the test case testmatchOneArg
  o->pushValue(string("three"));
  CPPUNIT_ASSERT(o->getValueCount()==3);
  CPPUNIT_ASSERT(o->isSet()==true);

  o->clearValues();
  CPPUNIT_ASSERT(o->getValueCount()==0);
  o->allowMutipleValues(1);
  CPPUNIT_ASSERT(o->m_allowMultiple==1);
  o->pushValue(string("one"));
  CPPUNIT_ASSERT(o->getValueCount()==1);
  o->pushValue(string("two"));
  CPPUNIT_ASSERT(o->getValueCount()==2);

  o->resetToDefault();
  CPPUNIT_ASSERT(o->getValueCount()==0);
  //o->dump();
  o->pushValue(string("one"));
  CPPUNIT_ASSERT(o->getValue()=="one");
  CPPUNIT_ASSERT(o->getValueCount()==1);
  o->pushValue(string("two"));
  CPPUNIT_ASSERT(o->getValueCount()==2);
  CPPUNIT_ASSERT(o->getValue(1)=="two");
  //o->dump();
  o->setValue("four");
  CPPUNIT_ASSERT(o->getValueCount()==1);
  CPPUNIT_ASSERT(o->getValue()=="four");
  //Goal is to receive a message: FATAL ERROR: Out of bounds. (idx>size) -OK
  NEGATIVE_TEST(o->getValue(-1),std::exception); 

  o->resetToDefault();
  CPPUNIT_ASSERT(o->getValue()==o->getDefaultValue());
  //o->dump();

  o->pushValue(string("one"));
  o->pushValue(string("two"));
  o->pushValue(string("three"));
  o->pushValue(string("four"));
  CPPUNIT_ASSERT(o->getValueCount()==4);
  //o->dump();
  vector<string> out;
  o->push_user_values_into(out);
  CPPUNIT_ASSERT(out.size()==4);
  CPPUNIT_ASSERT(out[0]=="one");
  CPPUNIT_ASSERT(out[1]=="two");
  CPPUNIT_ASSERT(out[2]=="three");
  CPPUNIT_ASSERT(out[3]=="four");
  delete opts;
}

PgOptions *getNewOpt() {
  PgOptions *my_opts=new PgOptions;
  my_opts->setUsage("Usage info.");
  PgOpt* opt;
  opt = my_opts->defineOption("m","mult",PgOpt::STRING_OPT,
                            "Multiples",
                            "");
  opt->allowMutipleValues(1);
  //cout<<endl;
  //my_opts->dump();

  return my_opts;
}

//test function vliber
void PgOptionsTest::testFunctions2() {

  Verbose::out(1, "PgOptionsTest::testFunctions2");
  PgOptions* my_opts;

  const char* argv[]={
    "progname",
    "-m","a0",
    "-m","b1",
    "-m","c2",
    "--mult","d3",
    "--mult","e4",
    "--mult","f5",
    "-m=equal", // 6
    "-m=", // 7
    "-m==", // 8
    "-m=9",
    "hello","world",
    NULL};

  my_opts = getNewOpt();

  CPPUNIT_ASSERT((my_opts->argc()==0)==true);
  POSITIVE_TEST(my_opts->argc());
  //my_opts->dump();

  my_opts->parseArgv(argv);
  CPPUNIT_ASSERT((my_opts->argc()==19)==true);
  CPPUNIT_ASSERT(my_opts->commandLine()=="progname -m a0 -m b1 -m c2 --mult d3 --mult e4 --mult f5 -m=equal -m= -m== -m=9 hello world");
  CPPUNIT_ASSERT(my_opts->argv(0)=="progname");
  CPPUNIT_ASSERT(my_opts->argv(1)=="-m");
  CPPUNIT_ASSERT(my_opts->argv(2)=="a0");
  CPPUNIT_ASSERT(my_opts->argv(3)=="-m");
  CPPUNIT_ASSERT(my_opts->argv(4)=="b1");
  CPPUNIT_ASSERT(my_opts->argv(16)=="-m=9");
  //my_opts->dump();

  my_opts->resetToDefaults();
  CPPUNIT_ASSERT((my_opts->argc()==0)==true);
  POSITIVE_TEST(my_opts->argc());
  //my_opts->dump();

  delete my_opts;
  my_opts = getNewOpt();

  my_opts->parseArgv(argv);
  CPPUNIT_ASSERT((my_opts->argc()==19)==true);
  CPPUNIT_ASSERT(my_opts->commandLine()=="progname -m a0 -m b1 -m c2 --mult d3 --mult e4 --mult f5 -m=equal -m= -m== -m=9 hello world");
  CPPUNIT_ASSERT(my_opts->argv(0)=="progname");
  CPPUNIT_ASSERT(my_opts->argv(1)=="-m");
  CPPUNIT_ASSERT(my_opts->argv(2)=="a0");
  CPPUNIT_ASSERT(my_opts->argv(3)=="-m");
  CPPUNIT_ASSERT(my_opts->argv(4)=="b1");
  CPPUNIT_ASSERT(my_opts->argv(16)=="-m=9");
  //my_opts->dump();

  my_opts->clearValues();
  CPPUNIT_ASSERT((my_opts->argc()==19)==true);
  CPPUNIT_ASSERT(my_opts->argv(1)=="-m");
  CPPUNIT_ASSERT(my_opts->argv(2)=="a0");
  CPPUNIT_ASSERT(my_opts->argv(3)=="-m");
  CPPUNIT_ASSERT(my_opts->argv(4)=="b1");
  CPPUNIT_ASSERT(my_opts->argv(16)=="-m=9");
  CPPUNIT_ASSERT(my_opts->getArg(0)==string("hello"));
  CPPUNIT_ASSERT(my_opts->getArg(1)=="world");
  //my_opts->dump();

  delete my_opts;
  my_opts = getNewOpt();

  my_opts->parseArgv(argv);
  CPPUNIT_ASSERT(my_opts->getProgName()=="progname");
  CPPUNIT_ASSERT((my_opts->getArgCount()==2)==true);
  CPPUNIT_ASSERT(my_opts->getArg(0)==string("hello"));
  CPPUNIT_ASSERT(my_opts->getArg(1)=="world");
  //my_opts->dump();

  CPPUNIT_ASSERT(my_opts->get("m")=="a0");
  CPPUNIT_ASSERT(my_opts->get("mult")=="a0");
  
  delete my_opts;
}

//test functions vliber
void PgOptionsTest::testFunctions3()
{
	Verbose::out(1, "PgOptionsTest::testFunctions3");
   PgOptions* my_opts=new PgOptions;
   my_opts->setUsage("test.");
   PgOpt* opt1=new PgOpt ;
   opt1->m_shortName="i";
   opt1->m_longName="integer-opt";
   opt1->m_help="Test1.";
   opt1->m_type=PgOpt::INT_OPT;
   opt1->m_defaultValue="10";
   PgOpt* opt1Out;
   opt1Out=my_opts->addPgOpt_nocopy(opt1);
   CPPUNIT_ASSERT(opt1Out->m_longName=="integer-opt");
   CPPUNIT_ASSERT(opt1Out->getValue()=="10");
   CPPUNIT_ASSERT(my_opts->findOpt("integer-opt")!=NULL);                  
   CPPUNIT_ASSERT(my_opts->findOpt("integer-opt-test")==NULL);

   std::map<std::string, PgOpt *>::iterator ii=my_opts->m_option_map.begin();
   //cout << (*ii).first << ": " << (*ii).second << endl;
   CPPUNIT_ASSERT((*ii).first=="i");
   CPPUNIT_ASSERT((*ii).second==opt1Out);
   ++ii;
   //cout << (*ii).first << ": " << (*ii).second << endl;
   CPPUNIT_ASSERT(((*ii).first=="integer-opt")==true);
   CPPUNIT_ASSERT((*ii).second==opt1Out);
   

   PgOpt* opt2=new PgOpt ;
   opt2->m_shortName="i1";
   opt2->m_longName="integer-opt1";
   opt2->m_help="Test2.";
   opt2->m_type=PgOpt::INT_OPT;
   opt2->m_defaultValue="11";
   PgOpt* opt2Out;
   opt2Out=my_opts->addPgOpt_nocopy(opt2);
   CPPUNIT_ASSERT(opt2Out->m_longName=="integer-opt1");
   CPPUNIT_ASSERT(opt2Out->getValue()=="11");
   CPPUNIT_ASSERT(my_opts->findOpt("i1")!=NULL);                  
   CPPUNIT_ASSERT(my_opts->findOpt("integer-opt-test1")==NULL);

   std::map<std::string, PgOpt *>::iterator i=my_opts->m_option_map.begin();
   //cout << (*i).first << ": " << (*i).second << endl;
   CPPUNIT_ASSERT((*i).first=="i");
   CPPUNIT_ASSERT((*ii).second==opt1Out);
   ++i;
   //cout << (*i).first << ": " << (*i).second << endl;
   CPPUNIT_ASSERT(((*i).first=="i1")==true);
   CPPUNIT_ASSERT((*i).second==opt2Out);
   ++i;
   //cout << (*i).first << ": " << (*i).second << endl;
   CPPUNIT_ASSERT(((*i).first=="integer-opt")==true);
   CPPUNIT_ASSERT((*i).second==opt1Out);
   ++i;
   //cout << (*i).first << ": " << (*i).second << endl;
   CPPUNIT_ASSERT(((*i).first=="integer-opt1")==true);
   CPPUNIT_ASSERT((*i).second==opt2Out);
   

   //my_opts->usage();

  delete my_opts;
  
}

//test function vliber

//test constructor vliber
void PgOptionsTest::testConstructor2() {

  Verbose::out(1, "PgOptionsTest::testConstructor2");
  PgOpt* o1;
  PgOpt* o2;
  PgOpt* o3;
  PgOptions* my_opts1=new PgOptions;
  PgOptions* my_opts2=new PgOptions;
  PgOptions* my_opts3=new PgOptions;

  my_opts1->setUsage("Constructors test 2.");

  o1=my_opts1->defineOption("-h","--help",PgOpt::BOOL_OPT,"A help message.","false");
  CPPUNIT_ASSERT(o1->m_allowMultiple==0);
  //my_opts1->dump();
  o2=my_opts2->defOpt("-h","--help",PgOpt::BOOL_OPT,"A help message.","false");
  CPPUNIT_ASSERT(o2->m_allowMultiple==0);
  //my_opts2->dump();
  CPPUNIT_ASSERT(o1->m_shortName==o2->m_shortName);
  CPPUNIT_ASSERT(o1->m_longName==o2->m_longName);
  CPPUNIT_ASSERT(o1->m_help==o2->m_help);
  CPPUNIT_ASSERT(o1->m_type==o2->m_type);
  CPPUNIT_ASSERT(o1->m_defaultValue==o2->m_defaultValue);
  CPPUNIT_ASSERT(my_opts1->findOpt("--help")->m_type==my_opts2->findOpt("-h")->m_type);

  o3=my_opts3->defOptMult("-h","--help",PgOpt::BOOL_OPT,"A help message.","false");
  CPPUNIT_ASSERT(o1->m_shortName==o3->m_shortName);
  CPPUNIT_ASSERT(o1->m_longName==o3->m_longName);
  CPPUNIT_ASSERT(o1->m_help==o3->m_help);
  CPPUNIT_ASSERT(o1->m_type==o3->m_type);
  CPPUNIT_ASSERT(o1->m_defaultValue==o3->m_defaultValue);
  CPPUNIT_ASSERT(o3->m_allowMultiple==1);
  //my_opts3->dump();
  //
  delete my_opts1;
  delete my_opts2;
  delete my_opts3;
}




//test function vliber
void PgOptionsTest::testmatchOneArg()
{
  Verbose::out(1, "PgOptionsTest::testmatchOneArg");
  PgOptions* my_opts1;
  
  my_opts1=new PgOptions;
  my_opts1->setUsage("Usage info.");
  PgOpt* opt1=my_opts1->defineOption("m","mult",PgOpt::STRING_OPT,"Multiples","");
    
  const char* argv1[]={
    "progname",
    "--b","9", //should be m or mult
    "hello","world",
    NULL};
  //Goal: FATAL ERROR: Don't recognize option: 'b'
  NEGATIVE_TEST(my_opts1->parseArgv(argv1),std::exception);

  const char* argv2[]={
    "progname",
    "--=9","b",//should be m or mult
    "hello","world",
    NULL};
  //Goal: FATAL ERROR: Shouldnt have a blank argument.
  NEGATIVE_TEST(my_opts1->parseArgv(argv2),std::exception);

  const char* argv3[]={
    "progname",
    "--no-m=1",//should not be no- for not boolean type (current type is PgOpt::STRING_OPT)
    "hello","world",
    NULL};
  //Goal: FATAL ERROR: Cant use '--no-' with 'm': Not a boolean option.
  NEGATIVE_TEST(my_opts1->parseArgv(argv3),std::exception);

  //multiple values and  allowMultiple==1 - positive result add a0, b1, and c2 to args
  const char* argv4[]={
    "progname",
    "-m","a0",
    "-m","b1",
	"-m","c2", NULL};
 opt1->allowMutipleValues(1); 
 CPPUNIT_ASSERT(opt1->m_allowMultiple==1);
 my_opts1->parseArgv(argv4);
 POSITIVE_TEST(my_opts1->parseArgv(argv4));
 //my_opts1->dump();

 //multiple values and  allowMultiple==0 - positive result add only last value (c2) to args
 //todo vliber I am not sure about that. Is it what app. suppose to do.
 //if argv4 array does not have NULL at the end, then exception will be thrown and app. does not catch it.
 opt1->allowMutipleValues(0);
 CPPUNIT_ASSERT(opt1->m_allowMultiple==0);
 my_opts1->parseArgv(argv4);
 POSITIVE_TEST(my_opts1->parseArgv(argv4));
 //my_opts1->dump();

 //one value and  allowMultiple==1 - positive result add only first value (a0) to args
 const char* argv5[]={
    "progname",
    "-m","a0",
    NULL};
 opt1->allowMutipleValues(1);
 CPPUNIT_ASSERT(opt1->m_allowMultiple==1);
 my_opts1->parseArgv(argv5);
 POSITIVE_TEST(my_opts1->parseArgv(argv5));
 //my_opts1->dump();

 //one value and  allowMultiple==0 - positive result add only first value (a0) to args
 opt1->allowMutipleValues(0);
 CPPUNIT_ASSERT(opt1->m_allowMultiple==0);
 my_opts1->parseArgv(argv5);
 POSITIVE_TEST(my_opts1->parseArgv(argv5));
 //my_opts1->dump();
  
}


