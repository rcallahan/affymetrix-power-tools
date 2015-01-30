////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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


//
#include "util/AffxString.h"
#include "util/TextFileCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
//
#include "util/CPPTest/Setup.h"

using namespace std;


/**
 * @class TextFileCheckTest
 * @brief cppunit class for testing functions from Util/TextFileCheck class.
 */


class TextFileCheckTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( TextFileCheckTest );
  CPPUNIT_TEST( test );
  
  CPPUNIT_TEST_SUITE_END();

public:
  
  void test();
  
};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( TextFileCheckTest );




void TextFileCheckTest::test()
{
  //positive
  cout<<endl;
  Verbose::out(1, "***TextFileCheck testcases***");
  Verbose::out(1, "TextFileCheckTest::test");
  TextFileCheck textCheck1 (INPUT+"/test1.txt","./input/test1Gold.txt", 0);
  std::string errorMsg ("");
  CPPUNIT_ASSERT(textCheck1.check(errorMsg));
  POSITIVE_TEST(textCheck1.check(errorMsg));

  //negative
  // goal is to receive a message: Unable to open generated file ./input1/test1.txt -OK
  TextFileCheck textCheck2 ("./input1/test1.txt","./input/test1Gold.txt", 0);
  CPPUNIT_ASSERT(!textCheck2.check (errorMsg));
  if (! textCheck2.check (errorMsg))
    cout << errorMsg << endl;
  
  // goal is to receive a message: Unable to open gold file input1\test1Gold.txt - OK
  TextFileCheck textCheck3 (INPUT+"/test1.txt","input1/test1Gold.txt", 0);
  CPPUNIT_ASSERT(!textCheck3.check(errorMsg)); 
  if (! textCheck3.check (errorMsg))
    cout << errorMsg << endl;
  
  // goal is to receive a message:e The generated file, ./input/test5Bad.txt, has fewer lines than the gold file, ./input/test1Gold.txt -OK
  TextFileCheck textCheck4 (INPUT+"/test5Bad.txt","./input/test1Gold.txt", 0);
  CPPUNIT_ASSERT(!textCheck4.check (errorMsg));
  if (! textCheck4.check (errorMsg))
    cout << errorMsg << endl;

  //goal is to receive a message: The generated file, ./input/test1.txt, has more lines than the gold file, ./input/test1Gold_Bad1.txt -OK
  TextFileCheck textCheck5 (INPUT+"/test1.txt","./input/test1Gold_Bad1.txt", 0);
  CPPUNIT_ASSERT(!textCheck5.check (errorMsg));
  if (! textCheck5.check (errorMsg))
    cout << errorMsg << endl;

  //goal is to receive a message: Mismatch reading generated file input\test1Bad.txt:
  //gold line: '3 0 3.000123 4.234 5.3456789 6.123456789 3'
  //generated line: '3 0 1.000123 4.234 5.3456789 6.123456789 3' -OK
  TextFileCheck textCheck6 (INPUT+"/test1Bad.txt","input/test1Gold.txt", 0);
  CPPUNIT_ASSERT(!textCheck6.check (errorMsg));
  if (! textCheck6.check (errorMsg))
    cout << errorMsg << endl;

  //goal is to receive a message: Mismatch reading generated file input\test1Bad.txt:
  //gold line: '6 0 6.000123 7.234 8.3456789 9.123456789 6'
  //generated line: '6 0 6.000123 7.234 8.3456789 3.123456789 6' -OK
  TextFileCheck textCheck7 (INPUT+"/test1Bad.txt","input/test1Gold.txt", 6);
  CPPUNIT_ASSERT(!textCheck7.check (errorMsg));
  if (! textCheck7.check (errorMsg))
    cout << errorMsg << endl;

  //goal is to receive a message: Unable to open generated file -OK 
  TextFileCheck textCheck8 ("","", 0);
  CPPUNIT_ASSERT(!textCheck8.check(errorMsg));
  if (! textCheck8.check (errorMsg)) 
    cout << errorMsg << endl;

  //additional empty white space at the end of the generated file
  //goal is to receive a message:
  //Mismatch reading generated file input/test1_space.txt:
  //gold line: '0 0 0.000123 1.234 2.3456789 3.123456789 1E-10'
  //generated line: '0 0 0.000123 1.234 2.3456789 3.123456789 1E-10  '
  TextFileCheck textCheck9 (INPUT+"/test1_space.txt","./input/test1Gold.txt", 0);
  CPPUNIT_ASSERT(!textCheck9.check(errorMsg));
  if (! textCheck9.check (errorMsg))
    cout << errorMsg << endl;

}
