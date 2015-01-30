////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
 * @file   BioTypesTest.cpp
 * @author csugne
 * @date   Mon Nov  7 10:04:34 PST 2005
 * 
 * @brief  Testing the PmAdjuster functions.
 * 
 */

//
#include "chipstream/BioTypes.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;
/**
 * @class BioTypesTest
 * @brief cppunit class for testing conversion functions.
 */
class BioTypesTest : public CppUnit::TestFixture {
public:
  BioTypesTest() {}

  CPPUNIT_TEST_SUITE( BioTypesTest );
  CPPUNIT_TEST( gtype_test );
  CPPUNIT_TEST( gender_test );
  CPPUNIT_TEST_SUITE_END();

  ///
  void gtype_test();
  ///
  void gender_test();
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( BioTypesTest );

void
BioTypesTest::gtype_test()
{
  printf("GType Test: ");

  // If we use more than one byte per GType, we will use
  // too much memory.
  CPPUNIT_ASSERT(sizeof(affx::GType)==1);

  // Heaven help us if they are equal!
  CPPUNIT_ASSERT(affx::AA!=affx::AB);
  CPPUNIT_ASSERT(affx::AA!=affx::BB);
  CPPUNIT_ASSERT(affx::AA!=affx::NN);
  //
  CPPUNIT_ASSERT(affx::AB!=affx::BB);
  CPPUNIT_ASSERT(affx::AB!=affx::NN);
  //
  CPPUNIT_ASSERT(affx::BB!=affx::NN);

  // should be inverses
  CPPUNIT_ASSERT(affx::GType_from_int(affx::GType_to_int(affx::AA))==affx::AA);
  CPPUNIT_ASSERT(affx::GType_from_int(affx::GType_to_int(affx::AB))==affx::AB);
  CPPUNIT_ASSERT(affx::GType_from_int(affx::GType_to_int(affx::BB))==affx::BB);
  CPPUNIT_ASSERT(affx::GType_from_int(affx::GType_to_int(affx::NN))==affx::NN);

  //
  CPPUNIT_ASSERT(!affx::GType_called(affx::NN));
  CPPUNIT_ASSERT(affx::GType_called(affx::AA));
  CPPUNIT_ASSERT(affx::GType_called(affx::AB));
  CPPUNIT_ASSERT(affx::GType_called(affx::BB));
  CPPUNIT_ASSERT(!affx::GType_called(affx::INVALID));
  CPPUNIT_ASSERT(!affx::GType_called(100));
  //
  CPPUNIT_ASSERT(affx::GType_valid(affx::NN));
  CPPUNIT_ASSERT(affx::GType_valid(affx::AA));
  CPPUNIT_ASSERT(affx::GType_valid(affx::AB));
  CPPUNIT_ASSERT(affx::GType_valid(affx::BB));
  CPPUNIT_ASSERT(!affx::GType_valid(affx::INVALID));
  CPPUNIT_ASSERT(!affx::GType_valid(100));

  printf("ok.\n");
}

void
BioTypesTest::gender_test()
{
  printf("Gender Test: ");

  // shouldnt be equal
  CPPUNIT_ASSERT(affx::Female!=affx::Male);
  CPPUNIT_ASSERT(affx::Female!=affx::UnknownGender);
  //
  CPPUNIT_ASSERT(affx::Male!=affx::UnknownGender);

  CPPUNIT_ASSERT(strlen(getGenderString(affx::Female))<MAX_GENDER_STRING_LENGTH);
  CPPUNIT_ASSERT(strlen(getGenderString(affx::Male  ))<MAX_GENDER_STRING_LENGTH);
  CPPUNIT_ASSERT(strlen(getGenderString(affx::UnknownGender))<MAX_GENDER_STRING_LENGTH);
  //
  printf("ok.\n");
}
