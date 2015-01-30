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
// sdk/util/FsTestDir.h
//

#ifndef _UTIL_FSTESDIR_H_
#define _UTIL_FSTESTDIR_H_

#include "util/FsPath.h"

//
#include <string>

/// @file   util/FsTestDir.h
/// @brief  A class for implementing APT_TEST_DIR enviroment variable. 

const std::string APT_TEST_ROOT_DIR = "affy-test-tmp";

/**
 * @class FsTestDir.h
 * @brief   A class for implementing APT_TEST_DIR enviroment variable. 
 * Full regression tests are now taking upwards of 50GB of test data space
 * that needs to be offloaded to temp space that is not backed up. 
 * In other words developer source code lives on a backed up network
 * drive whereas test data can be local for performance and cost.
 * An environment variable APT_TEST_DIR was chosen for ease of
 * portability. Any run time location of APT_TEST_DIR should
 * be favored or configuration of APT_TEST_DIR via "./configure". 
 *
 * When APT_TEST_DIR is set then the test output directory is:
 * ${APT_TEST_DIR}/test-generated/<testname>
 *
 * When APT_TEST_DIR is not set then the test output directory
 * is relative the to run time path:
 * test-generated/<testname>
 *
 * The usage here is commenserate with historical 'test-generated'
 * output convention.
 * 
 * Finally, deletion of the test directory is part of the API
 * and was left out of the existing Make process.
 * Using "make clean" will work when APT_TEST_DIR is not set
 * and "test-generated" is relative as it has been
 * historically. However "make clean" will not delete
 * the APT_TEST_DIR. Therefore developers should favor
 * using the clean option of this API. 
 *
 */
class FsTestDir {
  public:

  FsPath m_path;

  std::string asString() const {
    return m_path.asString();
  }
  
  AptErr_t setTestDir( const std::string & testcasepath, const bool clean = false, const int64_t bytescheck = 1073741824);


private:
  
};

#endif // _UTIL_FSTESTDIR_H_
