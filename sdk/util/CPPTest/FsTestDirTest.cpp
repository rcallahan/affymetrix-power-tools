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
 * @file   FsTestDirTest.cpp
 * @author Mybrid Spalding
 * @date   2011/10/28
 *
 * @brief  Class for implementing APT_TEST_DIR environment variable. 
 */

//
//
#include "util/FsTestDir.h"
#include "util/CPPTest/Setup.h"
//
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cstdlib>
#include <iostream>
#include <limits>
#undef max

class FsTestDirTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( FsTestDirTest );
  //
  CPPUNIT_TEST( test_setTestDir );
  //
  CPPUNIT_TEST_SUITE_END();

  void test_setTestDir();

};

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( FsTestDirTest );

void FsTestDirTest::test_setTestDir()
{
  std::string testpath("test_setTestDir");

  // Not testing the underlying file IO capabilities
  // as all file IO is done by Fs.cpp
  const char * savetestdir = getenv("APT_TEST_DIR");
  std::string testdir = "test_setTestDir";

#ifdef _WIN32
  _putenv_s("APT_TEST_DIR", testdir.c_str());
#else
  setenv("APT_TEST_DIR", testdir.c_str(), 1);
#endif
  
  Verbose::out(1,"");
  bool clean = false;
  bool done = false;
  while ( !done ) {
    std::string apt_test_dir;
    const char *env = getenv("APT_TEST_DIR");
    if ( env != NULL ) {
      apt_test_dir = env;
    }
    
    Verbose::out(1, "APT_TEST_DIR environment variable is: \"" + apt_test_dir + ToStr("\""));
    NEGATIVE_TEST(FsTestDir().setTestDir(testpath, clean = false, std::numeric_limits<int64_t>::max()  ), std::exception);
    NEGATIVE_TEST(FsTestDir().setTestDir(""), std::exception);
    NEGATIVE_TEST(FsTestDir().setTestDir("."), std::exception);
    NEGATIVE_TEST(FsTestDir().setTestDir("aa"), std::exception);
    NEGATIVE_TEST(FsTestDir().setTestDir("/var/tmp"), std::exception);
    NEGATIVE_TEST(FsTestDir().setTestDir("C:\temp"), std::exception);
    if ( getenv("APT_TEST_DIR") == NULL ) {
      done = true;
    }
    else {
#ifdef _WIN32      
      _putenv_s("APT_TEST_DIR", ""); // same as unsetenv
#else
      unsetenv("APT_TEST_DIR");
#endif      
    }
  }

  if ( savetestdir != NULL ) {
#ifdef _WIN32
    _putenv_s("APT_TEST_DIR", savetestdir);
#else    
    setenv("APT_TEST_DIR", savetestdir, 1);
#endif    
  }
}
