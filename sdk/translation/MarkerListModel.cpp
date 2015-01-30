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
 * @file   MarkerListModel.cpp
 * @author Mybrid Spalding
 * @date   Mon Jun  9 11:46:42 PDT 2008
 * @brief  Class for wrapping the probe set/marker list/file.
 */


#include "translation/MarkerListModel.h"
//
#include <cstring>
#include <sstream>
#include <string>
#include <vector>
//

using namespace affx;
using namespace std;

// This file only has one column, the probe set to filter on.
const TittmColumnDefinition MARKER_LIST_COLUMN_DEFINITIONS[] = {
  // column name,         empty?, column, valid regular expersion
  { std::string("probe-set"),      0,     0,  std::string("^[\\w\\d]+$"), NULL },
};

const TittmColumnDefinition MARKER_LIST_FIRST_COLUMN = { std::string("probe-set"), 1, 0, std::string("^[\\w\\d]+$"), NULL };

/*****************************************************************************/
/**
 * MarkerListModel::MarkerListModel
 *
 * FILE onstructor to read the probe set filter as a command line, TSV file.
 *
 * @param rte - RunTimeEnvironment which contains options
 * @param markerListFileName - the input file to read.
 */
/*****************************************************************************/
MarkerListModel::MarkerListModel(const RunTimeEnvironment &rte,
                                 const std::string & markerListFileName) :
    TranslationInputTsvTableModel(rte, markerListFileName,
                                  MARKER_LIST_COLUMN_DEFINITIONS,
                                  (sizeof(MARKER_LIST_COLUMN_DEFINITIONS) / sizeof(MARKER_LIST_COLUMN_DEFINITIONS[0])), false, NULL)
{

  for (int i = 0; i < m_rows.size(); i++) {
    m_probeSetList.push_back(m_rows[i][0]);
  }


}
// end MarkerListModel::MarkerListModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * MarkerListModel::MarkerListModel
 *
 * VECTOR constructor for the console std::vector data and convert the data
 * into the model.
 *
 * @param rte - RunTimeEnvironment which contains options
 *
 */
/*****************************************************************************/
MarkerListModel::MarkerListModel(RunTimeEnvironment &rte)
{


  if (rte.m_adtOpts.m_inputProbeSetVector.size() == 0) {
    return;
  }

  // HEADER

  TittmColumnDefinition probeSetTCD = MARKER_LIST_FIRST_COLUMN;

  m_columnDefinition.push_back(probeSetTCD);

  // DATA

  rte.m_adtOpts.m_inputMarkerListFile = rte.m_adtOpts.m_inputProbeSetVector[0];

  for (int i = 1; i < rte.m_adtOpts.m_inputProbeSetVector.size(); i++) {
    std::vector < std::string > newRow;
    m_probeSetList.push_back(rte.m_adtOpts.m_inputProbeSetVector[i]);
    newRow.push_back(rte.m_adtOpts.m_inputProbeSetVector[i]);
    m_rows.push_back(newRow);
  }

  return;

}
// end MarkerListModel::MarkerListModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * MarkerListModel::describeVerbose:
 * Synopsis:
 *
 * Run time debug routine that can be invoked with various verbose levels.
 * ADT_VERBOSE_INPUT_FILES is the default level.
 *
 * @param rte       - run time environment
 * @param override  - ADT_VERBOSE_INPUT_FILES by default,
 *                    otherwise the ADT_VERBOSE_ENUM level to output
 *
 *
 * @return -  nada
 */
/*****************************************************************************/
void  MarkerListModel::describeVerbose(const RunTimeEnvironment & rte,
                                       ADT_VERBOSE_ENUM overrideLevel)
{

  ADT_VERBOSE_ENUM level = overrideLevel == ADT_VERBOSE_NULL ? ADT_VERBOSE_INPUT_FILES : overrideLevel;

  if (rte.m_currentVerbosity < level) {
    return;
  }

  std::stringstream msgSStr;

  msgSStr << "MarkerListModel: count of specified probe sets: " << m_probeSetList.size() << endl;


  msgSStr << "MarkerListModel: specified probe sets: ";

  for (int i = 0; i < m_probeSetList.size(); i++) {
    msgSStr << m_probeSetList[i] << " ";
  }

  msgSStr << endl;

  Verbose::out(level, msgSStr.str(), false);

  return;
}
// end MarkerListModel::describeVerbose
/*****************************************************************************/



