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
 * @file   CopyNumberTableModel.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  DEPRECATED DMET2 ONLY: Class for wrapping the DMET2 copy number report TSV file.
 */


#include "translation/CopyNumberTableModel.h"
//
#include <sstream>

using namespace affx;
using namespace std;

const std::string DMET2_COPY_NUMBER_GENES[] = { "CYP2D6" };

const int NUM_DMET2_COPY_NUMBER_GENES = 1;

const TittmColumnDefinition CNTM_COLUMN_DEFINITIONS[] = {
// experiment sample prediction notes
// 16x6_NA18507_Rep_1 NA18507 1 CN > 0: that is, CN=1,2,3,...
  // column name,         empty?, column, valid regular expersion
  { std::string("experiment"), 0,     0,  std::string(".*"), NULL },
  { std::string("sample"),     0,     1,  std::string(".*"), NULL },
  { std::string("prediction"), 0,     2,  std::string("^-?[0-9]$"), NULL },
  { std::string("notes"),      0,     3,  std::string(".*"), NULL },
};

/*****************************************************************************/
/**
 * isGeneCopyNumberSensitive
 * Synopsis:
 *
 * DEPRECATED: DMET2 only
 *
 * HACK: This should be fixed. We currently check a hard coded C list
 * of gene names to see if they are copy number sensitive. This should
 * changed to come from some input file.
 *
 * @param gene - the gene name to check
 *
 * @return true - if copy number data applies to this gene.
 */
/*****************************************************************************/
bool isGeneCopyNumberSensitive(std::string gene)
{


  for (int i = 0; i < NUM_DMET2_COPY_NUMBER_GENES; i++) {
    if (gene == DMET2_COPY_NUMBER_GENES[i]) {
      return true;
    }
  }

  return false;

}
// end isGeneCopyNumberSensitive
/*****************************************************************************/
/*****************************************************************************/
/**
 * CopyNumberTableModel::CopyNumberTableModel
 * Synopsis:
 *
 * DEPRECATED: DMET2 only
 *
 * This is the TransferObject to the TsvFile DataAccessObject design
 * pattern. There is only ever one of these in existence, i.e. this
 * is a single instance class.
 *
 * copy number file -> copy number table model
 *
 *
 * Main constructor to read the copy number report file and stash the data
 * into the model.
 *
 * @param rte - Single instance RunTimeEnvironment
 * @param cntmFileName - the input file to read.
 */
/*****************************************************************************/
CopyNumberTableModel::CopyNumberTableModel(const RunTimeEnvironment &rte,
    const std::string & cntmFileName) :
    TranslationInputTsvTableModel(rte, cntmFileName, CNTM_COLUMN_DEFINITIONS,
                                  (sizeof(CNTM_COLUMN_DEFINITIONS) / sizeof(CNTM_COLUMN_DEFINITIONS[0])), false)
{


  for (int i = 0; i < m_rows.size(); i++) {

    if (m_rows[i][PREDICTION_INDEX] == "0") {

      std::string call;
      pcrecpp::StringPiece notesConsume(m_rows[i][NOTES_INDEX]);

      if (pcrecpp::RE("(\\S+\\/\\S+)").FindAndConsume(&notesConsume, &call)) {
        m_experimentCopyNumberCall[m_rows[i][EXPERIMENT_INDEX]] = call;
      }
    }
  }

}
// end CopyNumberTableModel::CopyNumberTableModel
/*****************************************************************************/



