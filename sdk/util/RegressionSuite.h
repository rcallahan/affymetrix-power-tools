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
 * @file   RegressionSuite.h
 * @author Chuck Sugnet
 * @date   Wed Mar  3 13:13:29 2010
 * 
 * @brief  Base class for holding basic options for regression data.
 */

#ifndef REGRESSION_SUITE_H
#define REGRESSION_SUITE_H

#include "util/PgOptions.h"
#include <string>

class RegressionSuite {

public:

  bool doValgrind() const;

  std::string getDatabase() const;

  std::string getSuiteName() const { return m_Name; }

  void parseArgv(const char * const * const argv);

  int leftoverArgs() {
    return m_Opts.getArgCount();
  }

private:

  void defineOptions();
  std::string m_Name;
  mutable PgOptions m_Opts;
  
};

#endif /* REGRESSION_SUITE_H */
