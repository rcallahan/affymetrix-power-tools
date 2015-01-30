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
 * @file   DataStoreTest.cpp
 * @author csugne
 * @date   Thu Oct 22 11:06:29 PDT 2009
 * 
 * @brief  Testing the DataStore functions.
 * 
 */
#ifndef DATASTORETEST_H
#define DATASTORETEST_H

#include "chipstream/DataStore.h"
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
 * @class DataStoreTest
 * @brief cppunit class for testing conversion functions.
 */
class DataStoreTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( DataStoreTest );
  CPPUNIT_TEST( testReadEntireTable );
  CPPUNIT_TEST_SUITE_END();

public:
  // blank test.
  void testReadEntireTable();
};

// Registers the fixture into the registry
CPPUNIT_TEST_SUITE_REGISTRATION( DataStoreTest );

void DataStoreTest::testReadEntireTable() {
  /* Read in data for answer */
  TableFile toNorm('\t','#',false,false);
  toNorm.open("input/norm-data.txt");
  vector<int> probeOrder(toNorm.numRows());
  for(int i = 0; i < probeOrder.size(); i++) 
    probeOrder[i] = probeOrder.size() - i - 1;
  DataStore dsMart("data-store.f5");
  dsMart.setProbeOrder(probeOrder);
  dsMart.initIntensity(toNorm.numRows(), 1);

  vector<float> data(toNorm.numRows());
  vector<char> gc(toNorm.numRows());
  for(int i = 0; i < toNorm.numRows(); i++) {
    data[i] = Convert::toFloat(toNorm.getData(i,0).c_str());
    gc[i] = (char) (i % 10);
  }
  dsMart.writeColumn(0, "col1", data);
  dsMart.setProbeGc(gc);
  vector<float> dataCheck(toNorm.numRows());
  vector<char> gcCheck(toNorm.numRows());
  dsMart.getProbeGc(gcCheck);
  dsMart.fillInCelData(0, dataCheck);
  for(int i = 0; i < toNorm.numRows(); i++) {
    bool dataOk = true, gcOk = true;
    dataOk &= dataCheck[i] == data[i];
    dataOk &= data[i] == dsMart.getProbeIntensity(i, 0, 0);
    gcOk &= gcCheck[i] == gc[i];
    gcOk &= gc[i] == dsMart.getProbeGc(i);
    CPPUNIT_ASSERT( dataOk );
    CPPUNIT_ASSERT( gcOk );
  }
}


#endif 
