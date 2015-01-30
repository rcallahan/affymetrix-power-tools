////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
// sdk/util/FsTestDir.cpp ---
//

//
#include "util/AffxString.h"
#include "util/AptErrno.h"
#include "util/Fs.h"
#include "util/FsTestDir.h"

//////////

AptErr_t FsTestDir::setTestDir( const std::string & testcasepath, const bool clean, const int64_t bytescheck ) {

  AffxString requestTest(testcasepath);
  const std::string api_name = "FsTestDir::setTestDir";
  requestTest.strip();
  
  if( requestTest.isBlankOrEmpty() ) {
    return m_path.setErr(APT_ERR, api_name + ": testcasepath cannot be an empty string.");
  }
  FsPath test(requestTest);
  if ( test.isAbsolute() || (requestTest[0] == '.'))  {
    return m_path.setErr(APT_ERR, api_name + ": testcasepath cannont be an absolute path or start with \".\".");
  }
  if ( requestTest.length() < 3 ) {
    return m_path.setErr(APT_ERR, api_name + ": testcasepath must be 3 or more characters.");
  }
  
  std::string topdir;

  if (getenv("APT_TEST_DIR") != NULL)  {
    topdir = getenv("APT_TEST_DIR");
    if ( ! Fs::dirExists( topdir ) ){
      return m_path.setErr(APT_ERR, topdir + "..." + api_name + ", APT_TEST_DIR directory not found.");
    }
  }

  if ( bytescheck > 0 ) {
    int64_t bytesfree = Fs::getFreeDiskSpace(topdir.empty()  ? "." : topdir);
    if ( bytescheck > bytesfree ) {
      std::string errmsg = api_name + ": byte check failed with bytes free \"" +
        ToStr(bytesfree) + "\" < \"" + ToStr(bytescheck) + "\" bytes requested.";
      return m_path.setErr(APT_ERR, errmsg);
    }
  }

  std::string completeDir;
  completeDir = Fs::join( topdir, APT_TEST_ROOT_DIR, requestTest );


  if ( Fs::dirExists( completeDir ) ) {
    if ( clean ) {
      Verbose::out(1, completeDir + ", " + api_name + " clean.");
      Fs::rm_rf(completeDir);
      Fs::mkdir(completeDir );
    }
    /// TODO: warn users directory exists?
  }
  else {
    Fs::mkdirPath(completeDir);
    Verbose::out(1, completeDir + ", " + api_name + ": directory created.");
  }

  
  AptErr_t aptErr = m_path.setDirPath(completeDir);

  if (aptErr == APT_OK) {
    std::string rootdir = Fs::join(topdir, APT_TEST_ROOT_DIR);
    std::string cmd=  "echo " + rootdir + "> " + APT_TEST_ROOT_DIR + ".txt";
    system(cmd.c_str());
#ifndef _WIN32
    if ( ! Fs::exists(APT_TEST_ROOT_DIR) ) {
      symlink(rootdir.c_str(),  APT_TEST_ROOT_DIR.c_str());
    }
#endif    
  }

  return aptErr;

}
