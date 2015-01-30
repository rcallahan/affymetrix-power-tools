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
 * @file   SparseMartTest.cpp
 * @author csugne
 * @date   Mon Nov  7 20:16:12 PST 2005
 * 
 * @brief  Testing the SparseMart functions.
 * 
 */
#ifndef SPARSEMARTTEST_H
#define SPARSEMARTTEST_H

#include "SparseMart.h"
#include "TableReader.h"
#include "util/TableFile.h"
#include "util/Convert.h"
#include "util/Err.h"
#include <string>
#include <vector>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

using namespace std;
/**
 * @class SparseMartTest
 * @brief cppunit class for testing conversion functions.
 */
class SparseMartTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SparseMartTest );
  CPPUNIT_TEST( testReadEntireTable );
  CPPUNIT_TEST( testReadHalfTable );
  CPPUNIT_TEST_SUITE_END();

public:
  // blank test.
  void testReadEntireTable();
  void testReadHalfTable();
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( SparseMartTest );

void SparseMartTest::testReadEntireTable() {
  /* Read in data for answer */
  TableFile toNorm('\t','#',false,false);
  toNorm.open("input/norm-data.txt");

  vector<bool> probes(100);
  for(int i = 0; i < probes.size(); i++) 
    probes[i] = true;
  SparseMart iMart(probes);

  vector<string> fileNames;
  TableReader reader(100);  
  reader.registerIntensityMart(&iMart);
  fileNames.push_back("input/norm-data.txt");
  reader.setFiles(fileNames);
  reader.readFiles();

  CPPUNIT_ASSERT( iMart.getCelFileCount() == 3);
  CPPUNIT_ASSERT( iMart.getProbeCount() == 100);

  const vector<string> names = iMart.getRowNames();
  for(int i = 0; i < names.size(); i++) {
    CPPUNIT_ASSERT( names[i] == ToStr(i) );
  }
  for(int probeIx = 0; probeIx < iMart.getProbeCount(); probeIx++) {
    for(int chipIx = 0; chipIx < iMart.getCelDataSetCount(); chipIx++) {
      float value = iMart.getProbeIntensity(probeIx, chipIx);
      float expected = Convert::toFloat(toNorm.getData(chipIx, probeIx).c_str());
      CPPUNIT_ASSERT( value == expected );
    }
  }
}

void SparseMartTest::testReadHalfTable() {
  /* Read in data for answer */
  TableFile toNorm('\t','#',false,false);
  toNorm.open("input/norm-data.txt");

  vector<bool> probes(100);
  for(int i = 0; i < probes.size(); i+=2) 
    probes[i] = true;
  SparseMart iMart(probes);

  vector<string> fileNames;
  TableReader reader(100);  
  reader.registerIntensityMart(&iMart);
  fileNames.push_back("input/norm-data.txt");
  reader.setFiles(fileNames);
  reader.readFiles();

  CPPUNIT_ASSERT( iMart.getCelFileCount() == 3);
  CPPUNIT_ASSERT( iMart.getProbeCount() == 100);

  const vector<string> names = iMart.getRowNames();
  for(int i = 0; i < names.size(); i++) {
    CPPUNIT_ASSERT( names[i] == ToStr(i) );
  }
  for(int probeIx = 0; probeIx < iMart.getProbeCount(); probeIx++) {
    for(int chipIx = 0; chipIx < iMart.getCelDataSetCount(); chipIx++) {
      if(probeIx % 2 == 0) {
        float value = iMart.getProbeIntensity(probeIx, chipIx);
        float expected = Convert::toFloat(toNorm.getData(chipIx, probeIx).c_str());
        CPPUNIT_ASSERT( value == expected );
      }
      else
        CPPUNIT_ASSERT( !iMart.isProbeAvailable(probeIx));
    }
  }
}

#endif 
