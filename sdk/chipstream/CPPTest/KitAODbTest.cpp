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
 * @file   KitAODbTest.cpp
 * @brief  Testing the PmAdjuster functions.
 */

//
#include "chipstream/KitAODb.h"
//
#include "util/Fs.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <stdio.h>

class KitAODbTest : public CppUnit::TestFixture {
public:
  KitAODbTest() {}

  CPPUNIT_TEST_SUITE( KitAODbTest );
  CPPUNIT_TEST( test_read );
  CPPUNIT_TEST_SUITE_END();

  ///
  void test_clear();
  void test_read();
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( KitAODbTest );

void
KitAODbTest::test_clear()
{
  printf("KitAODbTest: test_clear()\n");

  KitAODb kitdb;

  kitdb.clear();
}

void
KitAODbTest::test_read()
{
  printf("KitAODbTest: test_read()\n");


  KitAODb kitdb;

  kitdb.readTsv("./input/kit-ao/kit-ao-probeset-pc1.txt");

  printf("kitdb.kitdb.getClassiferAValue()==%f\n",kitdb.getClassiferAValue());  
  CPPUNIT_ASSERT(kitdb.getClassiferAValue()==-2.0);
  printf("kitdb.kitdb.getClassiferOValue()==%f\n",kitdb.getClassiferOValue());  
  CPPUNIT_ASSERT(kitdb.getClassiferOValue()==1.0);

  // check the counts.
  printf("kitdb counts: a=%d e=%d n=%d\n",
         kitdb.addedCount(),kitdb.entryCount(),kitdb.nameCount());

  // counts are for this file.
  CPPUNIT_ASSERT(kitdb.addedCount()==300);
  CPPUNIT_ASSERT(kitdb.entryCount()==300);

  // try a dump and reload.
  Fs::ensureWriteableDirPath("./output", false);
  std::string tmp_db="./output/kit-ao-probeset-pc1.out";
  kitdb.writeTsv(tmp_db);
  kitdb.clear();

  // make a fresh db
  KitAODb kitdb_tmp;
  kitdb_tmp.readTsv(tmp_db);
  CPPUNIT_ASSERT(kitdb_tmp.getClassiferAValue()==-2.0);
}
