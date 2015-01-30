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
 * @file   ExperimentFileTableModel.cpp
 * @author Mybrid Spalding
 * @date   Mon Jun  9 11:46:42 PDT 2008
 * @brief  Class for wrapping the probe set/ marker list file.
 */


#include "calvin_files/utils/src/FileUtils.h"
#include "translation/ExperimentFileTableModel.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include "pcrecpp.h"
//
#include <cstring>
#include <string>
#include <vector>
//

using namespace std;

// This file only has one column, the file name of a CHP file.
const TittmColumnDefinition EXPERIMENT_FILE_TABLE_COLUMN_DEFINITIONS[] = {
  // column name,         empty?, column, valid regular expersion
  { std::string("chp_files"),  0,  0,  std::string(".*\\.[cC][hH][pP]$"), NULL },
};


/*****************************************************************************/
/**
 * ExperimentFileTableModel::ExperimentFileTableModel
 * Synopsis:
 *
 * An input list of CHP files can be given at the command line and if so
 * is the file is  read via the TsvFile object, in keeping with Affy practice.
 *
 *
 * Main constructor to read the file list and stash the data
 * into the standard TranslationInputTsvTableModel.
 *
 * @param rte - RunTimeEnvironment which contains options
 * @param probeSetFilterFileName - the input file to read.
 */
/*****************************************************************************/
ExperimentFileTableModel::ExperimentFileTableModel(const RunTimeEnvironment &rte,
    const std::string & experimentListFileName) :
    TranslationInputTsvTableModel(rte, experimentListFileName,
                                  EXPERIMENT_FILE_TABLE_COLUMN_DEFINITIONS,
                                  (sizeof(EXPERIMENT_FILE_TABLE_COLUMN_DEFINITIONS) / sizeof(EXPERIMENT_FILE_TABLE_COLUMN_DEFINITIONS[0])), false, NULL)
{

  pcrecpp::RE reCommentedOut("^\\s*#");

  for (int i = 0; i < m_rows.size(); i++) {
    if (! reCommentedOut.PartialMatch(m_rows[i][0])) {
      m_experimentFiles.push_back(m_rows[i][0]);
    }
  }

}
// end ExperimentFileTableModel::ExperimentFileTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * ExperimentFileTableMode::validateFiles:
 * Synopsis:
 *
 *  Check if  CHP files exist and if not output a message and return false.
 *
 *
 * @return - true if all files exist.
 */
/*****************************************************************************/
bool ExperimentFileTableModel::validateFiles()
{

  bool okFiles = true;

  for (int i = 0; i < m_experimentFiles.size(); i++) {

    if (!  Fs::fileExists(m_experimentFiles[i].c_str())) {
      okFiles = false;
      Verbose::out(ADT_VERBOSE_NORMAL, m_experimentFiles[i] + ": file not found.");
    }

  }


  return okFiles;

}
// end ExperimentFileTableModel::validateFiles
/*****************************************************************************/
