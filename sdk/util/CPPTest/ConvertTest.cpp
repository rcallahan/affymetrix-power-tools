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
 * @file   ConvertTest.cpp
 * @author Chuck Sugnet
 * @date   Wed May  4 15:45:33 2005
 * 
 * @brief  Testing the convert functions.
 * last change by vliber on 01/22/08
 */

//
#include "util/Convert.h"
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <climits>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>

//
#include "util/CPPTest/Setup.h"

using namespace std;

#define NaN numeric_limits<double>::quiet_NaN()
#define Inf numeric_limits<double>::infinity()

/**
 * @class ConvertTest
 * @brief cppunit class for testing conversion functions.
 */
class ConvertTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ConvertTest );
  CPPUNIT_TEST( testInteger );
  CPPUNIT_TEST( testUnsignedInteger );
  CPPUNIT_TEST( testFloat );
  CPPUNIT_TEST( testBool );
  CPPUNIT_TEST( testToString );
  CPPUNIT_TEST( testToStr );
  CPPUNIT_TEST( testStrToIntVec );
  CPPUNIT_TEST_SUITE_END();

public:
  /** Test converting integers. */
  void testInteger();
  /** Test converting unsigned integers. */
  void testUnsignedInteger();
  /** Test converting floats. */
  void testFloat();
  /** Test converting booleans. */
  void testBool();
  /** Test for converting things to strings. */
  void testToString();
  /** Test for converting strings template. */
  void testToStr();
  /** */
  void testStrToIntVec();
};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ConvertTest );

/**
 * Make sure that toString functions work.
 */
void ConvertTest::testToString() {
  Verbose::out(1, "ConvertTest::testToString");
  CPPUNIT_ASSERT( Convert::toString(10) == "10" );
  CPPUNIT_ASSERT( Convert::toString(-10) == "-10" );
  CPPUNIT_ASSERT( Convert::toString(+10) == "10" );

  CPPUNIT_ASSERT( Convert::toString(10.5) == "10.5" );
  CPPUNIT_ASSERT( Convert::toString(-10.5) == "-10.5" );
  CPPUNIT_ASSERT( Convert::toString(+10.5) == "10.5" );
  //************vliber***************  
  CPPUNIT_ASSERT( (Convert::toString(0) == "0" )==true);
  CPPUNIT_ASSERT(( Convert::toString(-2147483647) == "-2147483647")==true );
  CPPUNIT_ASSERT(( Convert::toString(2147483647) == "2147483647")==true );
  CPPUNIT_ASSERT( (Convert::toString(0.0) == "0")==true ); 
  //negative
  // input parameter is not a double
  //NEGATIVE_TEST(Convert::toString('a'),std::exception);// todo vliber function does not throw any exceptions
  long i=-2147483648;
  POSITIVE_TEST(Convert::toString((double)i));
  //NEGATIVE_TEST(Convert::toString((int)i),std::exception);//todo vliber function does not throw any exceptions
}
  
/** 
 * Make sure that the integer conversion works.
 */
void ConvertTest::testInteger() {
	Verbose::out(1, "");
	Verbose::out(1, "ConvertTest::testInteger");
  int i = 0;
  bool success = 0;
  i = Convert::toIntCheck("1", &success);
  CPPUNIT_ASSERT( i == 1 && success == true);

  i = Convert::toIntCheck("Q", &success);
  CPPUNIT_ASSERT( i == 0 && success == false);

  i = Convert::toIntCheck("Q", &success);
  CPPUNIT_ASSERT( i == 0 && success == false);

  i = Convert::toIntCheck("-12", &success);
  CPPUNIT_ASSERT( i == -12 && success == true);

  i = Convert::toIntCheck("+12", &success);
  CPPUNIT_ASSERT( i == 12 && success == true);
  //************vliber***************
  i = Convert::toIntCheck("2147483647", &success);
  CPPUNIT_ASSERT( i == 2147483647 && success == true);
  i = Convert::toIntCheck("-2147483648", &success);
  CPPUNIT_ASSERT( i == -2147483648 && success == true);

  i = Convert::toIntCheck("q", &success);
  CPPUNIT_ASSERT( i == 0 && success == false);
  POSITIVE_TEST(Convert::toIntCheck("q", &success));

  i = Convert::toIntCheck("0", &success);
  CPPUNIT_ASSERT( i == 0 && success == true);

  i = Convert::toIntCheck("#", &success);
  CPPUNIT_ASSERT( i == 0 && success == false);
  POSITIVE_TEST(Convert::toIntCheck("#", &success));

  i = Convert::toIntCheck("2147483648", &success);
  CPPUNIT_ASSERT( i == 0 && success == false);
  POSITIVE_TEST(Convert::toIntCheck("2147483648", &success));
  
  CPPUNIT_ASSERT_EQUAL(0,Convert::toInt("0"));
  CPPUNIT_ASSERT_EQUAL(-2147483647,Convert::toInt("-2147483647"));
  CPPUNIT_ASSERT_EQUAL(2147483647,Convert::toInt("2147483647"));
}

/** 
 * Make sure that the integer conversion works.
 */
void ConvertTest::testUnsignedInteger() {
	Verbose::out(1, "ConvertTest::testUnsignedInteger");
  unsigned int i = 0;
  bool success = 0;
  string s = ToStr(UINT_MAX);
  i = Convert::toUnsignedIntCheck("1", &success);
  CPPUNIT_ASSERT( i == 1 && success == true);

  i = Convert::toUnsignedIntCheck("Q", &success);
  CPPUNIT_ASSERT( i == 0 && success == false);

  i = Convert::toUnsignedIntCheck(s.c_str(), &success);
  CPPUNIT_ASSERT( i == 0 && success == false);

  i = Convert::toUnsignedIntCheck("+12", &success);
  CPPUNIT_ASSERT( i == 12 && success == true);
  //************vliber***************
  i = Convert::toUnsignedIntCheck("0", &success);
  CPPUNIT_ASSERT( i == 0 && success == true);
  i = Convert::toUnsignedIntCheck("-12", &success);
  //CPPUNIT_ASSERT( i == -12 && success == true); // todo vliber: fails on Linux
  i = Convert::toUnsignedIntCheck("g", &success);
  CPPUNIT_ASSERT( i == 0 && success == false);
}


/**
 * Make sure that float conversions work when they should.
 * Currently not testing for exact equality would need to 
 * write a "close enough" comparitor.
 */
void ConvertTest::testFloat() {
  Verbose::out(1, "ConvertTest::testFloat");
  float f = 0;
  bool success = true;
  f = Convert::toFloatCheck("1.77864e-322", &success);
  CPPUNIT_ASSERT( Convert::doubleCloseEnough(f, 0.0) && success == true);

  // bool sameVal = false;
  // f = Convert::toFloatCheck("-1.77864e3222", &success);
  // int infX =  ( isinf(f) ) ;
  // sameVal = (infX != 0 && f < 0);
  // CPPUNIT_ASSERT( sameVal && success == true);

  // sameVal = false;
  // f = Convert::toFloatCheck("1.77864e3222", &success);
  // infX =  ( isinf(f) ) ;
  // sameVal = infX != 0 && f > 0;
  // CPPUNIT_ASSERT( sameVal && success == true);

  f = Convert::toFloatCheck("1.0", &success);
  CPPUNIT_ASSERT( Convert::doubleCloseEnough(f, 1.0) && success == true);


  f = Convert::toFloatCheck("-1.0", &success);
  CPPUNIT_ASSERT( Convert::doubleCloseEnough(f, -1.0) && success == true);

  f = Convert::toFloatCheck("-1.G", &success);
  CPPUNIT_ASSERT( Convert::doubleCloseEnough(f, 0) && success == false);
  //************vliber***************
  f = Convert::toFloatCheck("0.0", &success);
  CPPUNIT_ASSERT( Convert::doubleCloseEnough(f, 0.0) && success == true);
}

/** 
 * Make sure that we understand bools and fail on junk.
 */
void ConvertTest::testBool() {
	Verbose::out(1, "ConvertTest::testBool");
  bool success = true, value = true;

  // T!
  value = Convert::toBoolCheck("true", &success);
  CPPUNIT_ASSERT( value && success);
  value = Convert::toBoolCheck("TRUE", &success);
  CPPUNIT_ASSERT( value && success);
  value = Convert::toBoolCheck("1", &success);
  CPPUNIT_ASSERT( value && success);

  // F!
  value = Convert::toBoolCheck("false", &success);
  CPPUNIT_ASSERT( !value && success);
  value = Convert::toBoolCheck("FALSE", &success);
  CPPUNIT_ASSERT( !value && success);
  value = Convert::toBoolCheck("0", &success);
  CPPUNIT_ASSERT( !value && success);

  // bad!
  value = Convert::toBoolCheck("junk", &success);
  CPPUNIT_ASSERT( !value && !success);
  value = Convert::toBoolCheck("", &success);
  CPPUNIT_ASSERT( !value && !success);
}

/** Test for converting strings template. */
void ConvertTest::testToStr() {
	Verbose::out(1, "ConvertTest::testToStr");
  CPPUNIT_ASSERT( string("-5.295") == ToStr(-5.295) );
  CPPUNIT_ASSERT( string("100") == ToStr(100) );
  CPPUNIT_ASSERT( string("Some stupid string") == ToStr("Some stupid string") );
  CPPUNIT_ASSERT( string("Cast to string") == ToStr(string("Cast to string")) );
   //************vliber***************
   CPPUNIT_ASSERT( string("4294967295") == ToStr(4294967295) );
   CPPUNIT_ASSERT( string("0") == ToStr(0) );
   CPPUNIT_ASSERT( string("-2147483647") == ToStr(-2147483647) );
   CPPUNIT_ASSERT( string("true") == ToStr(true) );
   //************rsatin***************
   CPPUNIT_ASSERT( string("-inf") == ToStr(-Inf) );
   CPPUNIT_ASSERT( string("inf") == ToStr(Inf) );
   //CPPUNIT_ASSERT( string("nan") == ToStr(NaN) );  // TO DO rsatin: returns "1.#QNAN" in Windows

}

void ConvertTest::testStrToIntVec() {
	Verbose::out(1, "ConvertTest::testStrToIntVec");
  std::vector<int> vec;

  Convert::strToIntVec("",',',vec);
  CPPUNIT_ASSERT(vec.size()==0);

  Convert::strToIntVec("1",',',vec);
  CPPUNIT_ASSERT(vec.size()==1);
  CPPUNIT_ASSERT(vec[0]==1);

  Convert::strToIntVec("0,1,2,3",',',vec);
  CPPUNIT_ASSERT(vec.size()==4);
  CPPUNIT_ASSERT(vec[0]==0);
  CPPUNIT_ASSERT(vec[1]==1);
  CPPUNIT_ASSERT(vec[2]==2);
  CPPUNIT_ASSERT(vec[3]==3);
}

