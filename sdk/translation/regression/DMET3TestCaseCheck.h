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
 * @file   DMET3TestCaseCheck.h
 * @author Mybrid Spalding
 * @date   Thu Nov 13 16:54:29 PST 2008
 * @brief  DMET3 test case report check class. 
 */

#ifndef TRANSLATION_REGRESSION_DMET3_TEST_CASE_CHECK_H
#define TRANSLATION_REGRESSION_DMET3_TEST_CASE_CHECK_H

//
#include "translation/RunTimeEnvironment.h"
//
#include "util/RegressionCheck.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cassert>
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

class DMET3ReportTestCheckFiles {
public:
  std::string m_chpFiles;
  std::string m_translationFile;
  std::string m_annotationFile;
  std::string m_markerFile;
  std::string m_overrideFile;
  std::map< std::string, std::string> m_resultFiles;
 

};


class DMET3TestCaseCheck : public RegressionCheck
{
public:

  bool check(std::string & msg);

  std::string m_testCase;
  std::vector< std::string > m_testCaseRow;
  RunTimeEnvironment m_rte;
  DMET3ReportTestCheckFiles m_tcf;
    
  DMET3TestCaseCheck(RunTimeEnvironment & rte,
                     std::string testCaseRow,
                     std::vector< std::string > &testCase,
                     DMET3ReportTestCheckFiles &tcf );

private:
  bool checkAllColumns(std::string & msg);
  

};




#endif /* TRANSLATION_REGRESSION_DMET3_REPORT_CHECK_H */
