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



#include "util/BaseEngine.h"
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
//
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif
//
#include "util/CPPTest/Setup.h"

using namespace std;

class BaseEngineTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( BaseEngineTest );
  CPPUNIT_TEST( getMetaDataDescriptionTest );
  CPPUNIT_TEST_SUITE_END();

public:
  
  void getMetaDataDescriptionTest();
  
};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( BaseEngineTest );


void BaseEngineTest::getMetaDataDescriptionTest() {
  Err::setThrowStatus(true);
  vector<pair<string, string> > data;

  const char* argv[] = {
    "dummy",
    "--meta-data-info",
    "key1=val1",
    "--meta-data-info",
    "key2=val2",
    "--meta-data-info",
    "key3=val3",
    NULL};
  BaseEngine base;
  
  base.parseArgv(argv);;
  data = base.getMetaDataDescription();
  for (int i = 0; i < data.size(); i++) {
    pair<string,string> p = data[i];
    if (p.first == "key1")
      CPPUNIT_ASSERT(p.second == "val1");
    if (p.first == "key2")
      CPPUNIT_ASSERT(p.second == "val2");
    if (p.first == "key3")
      CPPUNIT_ASSERT(p.second == "val3");
  }


  const char *argvFail1[] = {
    "dummy",
    "--meta-data-info",
    "=val1",
    "--meta-data-info",
    "key2=val2",
    "--meta-data-info",
    "key3=val3",
    NULL};
    
  base.parseArgv(argvFail1);
  NEGATIVE_TEST(base.getMetaDataDescription(), Except);


  const char *argv2[] = {
    "dummy",
    "--meta-data-info",
    "key1=val=1",
    "--meta-data-info",
    "key2=val2",
    "--meta-data-info",
    "key3=val3",
    NULL};

  BaseEngine base2;
  base2.parseArgv(argv2);
  data = base2.getMetaDataDescription();
  // actually we dont throw anymore.
  // CPPUNIT_ASSERT(false); // should throw
  for (int i = 0; i < data.size(); i++) {
    pair<string,string> p = data[i];
    if (p.first == "key1")
      CPPUNIT_ASSERT(p.second == "val=1");
    if (p.first == "key2")
      CPPUNIT_ASSERT(p.second == "val2");
    if (p.first == "key3")
      CPPUNIT_ASSERT(p.second == "val3");
  }

  const char *argvFail2[] = {
    "dummy",
    "--meta-data-info",
    "key1=",
    "--meta-data-info",
    "key2=val2",
    "--meta-data-info",
    "key3=val3",
    NULL};
  
  base2.parseArgv(argvFail2);
  NEGATIVE_TEST(base2.getMetaDataDescription(), Except);
}

