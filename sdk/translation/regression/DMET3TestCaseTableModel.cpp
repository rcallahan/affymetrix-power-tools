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
 * @brief  Class representing a DMET3 Test Case
 */

#include "translation/regression/DMET3TestCaseTableModel.h"
//
#include <cassert>
#include <iostream>
#include <sstream>
//

///////////////////////////////////////////////////////////////////////////////
// DMET3TestCaseTableModel.cpp
///////////////////////////////////////////////////////////////////////////////

//Experiment Gene ProbeSet A1 A2 Ref Var Call
const TittmColumnDefinition DMET3_TEST_CASE_TABLE_MODEL_COLUMN_DEFINITIONS[] = { // column_name, empty_ok, index, valid_regex
  { std::string("CHP Files"),      0,  0,
    std::string("(?:.*\\.)(?:(?:[xX][mM][lL])|(?:[cC][Hh][pP])|(?:[tT][xX][tT]))$"), NULL },
  { std::string("Translation File"),
                                   0,  1,
    std::string(".*(?:translation)|(?:\\.csv)|(?:\\.txt)$"), NULL},
  { std::string("Annotation File"),0,  2,
    std::string("^.*(?:\\.csv)$"), NULL },
  { std::string("Marker File"),    1,  3,
    std::string("^.*(?:\\.[Tt][Xx][Tt])$"), NULL },
  { std::string("Override File"),  1,  4,
    std::string("^.*(?:\\.rpt)$"), NULL },
  { std::string("Base Name"),      0,  5, std::string(".*"), NULL },
  { std::string("Expected Result Reports"),
                                   0,  6,
    std::string("^(?:(?:comprehensive|summary|uncalled),?)+$"), NULL },
  { std::string("Test Codes"),     1,  7, std::string(".*"), NULL },
  { std::string("Log Message"),    0,  8, std::string(".*"), NULL },
  { std::string("Failed Message"), 0,  9, std::string(".*"), NULL},
};


/*****************************************************************************/
/**
 * DMET3TestCaseTableModel::DMET3TestCaseTableModel
 * Synopsis: Default constructor.
 *
 * @param rte - the runtime environment
 * @param testCaseFile - the TsvFile with the test case data.
 * 
 */
/*****************************************************************************/
DMET3TestCaseTableModel::DMET3TestCaseTableModel(RunTimeEnvironment & rte,
                                   const std::string testCaseFile) :
    TranslationInputTsvTableModel(rte,
                                  testCaseFile,
                                  DMET3_TEST_CASE_TABLE_MODEL_COLUMN_DEFINITIONS,
                                  (size_t) sizeof(DMET3_TEST_CASE_TABLE_MODEL_COLUMN_DEFINITIONS) / sizeof(DMET3_TEST_CASE_TABLE_MODEL_COLUMN_DEFINITIONS[0]),
                                  true)
{



}
// end DMET3TestCaseTableModel::DMET3TestCaseTableModel
/*****************************************************************************/
