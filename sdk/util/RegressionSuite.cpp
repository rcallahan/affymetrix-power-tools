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

//
#include "util/RegressionSuite.h"
//
#include "util/Fs.h"
#include "util/Util.h"

bool RegressionSuite::doValgrind() const {
  return m_Opts.getBool("valgrind");
}

std::string RegressionSuite::getDatabase() const {
  return m_Opts.get("database");
}

void RegressionSuite::defineOptions() {
    m_Opts.defineOption("", "database", PgOpt::STRING_OPT,
                     "Path to sqlite database to add results to.",
                     "");
    m_Opts.defineOption("", "valgrind", PgOpt::BOOL_OPT,
                     "Should the tests be run under valgrind.",
                     "false");
}

void RegressionSuite::parseArgv(const char * const * const argv) {
  defineOptions();
  m_Name = Fs::basename(argv[0]);
  m_Opts.parseArgv(argv);
}

