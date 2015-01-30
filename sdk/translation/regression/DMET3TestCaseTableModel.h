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
 * @file   DMET3TestCaseTableModel.h
 * @author Mybrid Spalding
 * @date   Thu Nov 13 15:34:22 PST 2008
 * @brief  Class thet represents a DMET3 test case. 
 */

#ifndef TRANSLATION_REGRESSION_DMET3_TEST_CASE_TABLE_MODEL_H
#define TRANSLATION_REGRESSION_DMET3_TEST_CASE_TABLE_MODEL_H


//
#include "translation/TranslationInputTsvTableModel.h"
#include "translation/regression/ATDRegression.h"
//
#include <cstring>
#include <set>
#include <string>
//


class DMET3TestCaseTableModel : public TranslationInputTsvTableModel
{
public:

  enum {
    CHP_FILES,
    TRANSLATION_FILE,
    ANNOTATION_FILE,
    MARKER_FILE,
    OVERRIDE_FILE,
    BASE_NAME, 
    EXPECTED_RESULT_REPORTS,
    TEST_CODES,
    LOG_MESSAGE,
    FAILED_MESSAGE,
  };


  DMET3TestCaseTableModel() {}
  DMET3TestCaseTableModel(RunTimeEnvironment &rte,
                          const std::string testCaseFile);


};



#endif /* TRANSLATION_REGRESSION_DMET3_TEST_CASE_TABLE_MODEL_H */

