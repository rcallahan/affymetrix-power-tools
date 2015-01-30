////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
 * @file   DMET2Check.h
 * @author Mybrid Spalding
 * @date   Wed Jun  4 10:18:22 PDT 2008
 * @brief  Allele translation table regression class for DMET2 files.
 */

#ifndef TRANSLATION_REGRESSION_DMET2_CHECK_H
#define TRANSLATION_REGRESSION_DMET2_CHECK_H

//
#include "translation/RunTimeEnvironment.h"
//
#include "util/RegressionCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cassert>
#include <cstring>
#include <string>
//

//

using namespace std;


class DMET2Check : public RegressionCheck
{
public:


  string m_dmet2File;
  string m_dmet3File;
  string m_genoFile;
  RunTimeEnvironment *m_rte;

  DMET2Check(RunTimeEnvironment & rte,
             const std::string &dmet2File,
             const std::string &dmet3File,
             const std::string &genoFile);

  DMET2Check() {
    m_rte = NULL;
  }

};

class DMET2MarkerCheck : public DMET2Check
{
public:

  bool check(string & msg);

  DMET2MarkerCheck(RunTimeEnvironment & rte,
                   const std::string &dmet2MarkerFile,
                   const std::string &dmet3MarkerFile,
                   const std::string &genoFile):
      DMET2Check(rte, dmet2MarkerFile, dmet3MarkerFile, genoFile) {};

  DMET2MarkerCheck() {};
  ~DMET2MarkerCheck() {};

};


class DMET2HaplotypeCheck : public DMET2Check
{
public:

  bool check(string & msg);

  DMET2HaplotypeCheck(RunTimeEnvironment & rte,
                      const std::string &dmet2HaplotypeFile,
                      const std::string &dmet3HaplotypeFile,
                      const std::string &genoFile):
      DMET2Check(rte, dmet2HaplotypeFile, dmet3HaplotypeFile, genoFile) {};

  DMET2HaplotypeCheck() {};
  ~DMET2HaplotypeCheck() {};

};


#endif /* TRANSLATION_REGRESSION_DMET2_CHECK_H */
