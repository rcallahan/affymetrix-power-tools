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
 * @file   DMET2Model.cpp
 * @author Mybrid Spalding
 * @date   Wed Jun  4 11:23:13 PDT 2008
 * @brief  Allele translation table regression class for all DMET2 input files.
 */

#include "translation/regression/DMET2Model.h"
//
#include <cassert>
#include <iostream>
#include <sstream>
//

///////////////////////////////////////////////////////////////////////////////
// DMET2MarkerModel
///////////////////////////////////////////////////////////////////////////////

//Experiment Gene ProbeSet A1 A2 Ref Var Call
const TittmColumnDefinition DMET2_MARKER_MODEL_COLUMN_DEFINITIONS[] = { // column_name, empty_ok, index, valid_regex
  { std::string("Experiment"), 0,  0, std::string(".*"), NULL },
  { std::string("Gene"),       0,  1, std::string("^[A-Z\\d]+$"), NULL},
  { std::string("AssayId"),    0,  2, std::string("^(?:\\d+)$"), NULL },
  { std::string("A1"),         0,  3, std::string(".*"), NULL },
  { std::string("A2"),         1,  4, std::string(".*"), NULL },
  { std::string("Ref"),        0,  5, std::string(".*"), NULL },
  { std::string("Var"),        0,  6, std::string(".*"), NULL },
  { std::string("Call"),       1,  7,
    std::string("^(?:Ref|Var)/(?:Ref|Var)"), NULL },
};

const TittmColumnDefinition DMET3_MARKER_MODEL_COLUMN_DEFINITIONS[] = { // column_name, empty_ok, index, valid_regex
  { std::string("Experiment"), 0,  0, std::string(".*"), NULL },
  { std::string("Gene"),       0,  1, std::string("^[A-Z\\d]+$"), NULL},
  { std::string("ProbeSet"),   0,  2, std::string("^(?:\\d+)$"), NULL },
  { std::string("A1"),         0,  3, std::string(".*"), NULL },
  { std::string("A2"),         1,  4, std::string(".*"), NULL },
  { std::string("Ref"),        0,  5, std::string(".*"), NULL },
  { std::string("Var"),        0,  6, std::string(".*"), NULL },
  { std::string("Call"),       1,  7,
    std::string("^(?:Ref|Var)/(?:Ref|Var)"), NULL },
};


/*****************************************************************************/
/**
 * DMET2MarkerModel::DMET2MarkerModel
 * Synopsis: Default constructor.
 *
 * @param rte - the runtime environment
 * @param markerFile - the TsvFile with the data.
 * @return - description
 */
/*****************************************************************************/
DMET2MarkerModel::DMET2MarkerModel(RunTimeEnvironment *rte,
                                   const string markerFile,
                                   bool isDMET2) :
    TranslationInputTsvTableModel(*rte,
                                  markerFile,
                                  isDMET2 ? DMET2_MARKER_MODEL_COLUMN_DEFINITIONS : DMET3_MARKER_MODEL_COLUMN_DEFINITIONS,
                                  (size_t) sizeof(DMET2_MARKER_MODEL_COLUMN_DEFINITIONS) / sizeof(DMET2_MARKER_MODEL_COLUMN_DEFINITIONS[0]),
                                  true)
{

  // Create a hash of the experiment, gene probeSet for comparison with the
  // DMET3 output.

  for (int i = 0; i < m_rows.size(); i++) {

    string experiment = m_rows[i][0];
    string gene       = m_rows[i][1];
    string probeSet    = m_rows[i][2];

    string ega        = experiment + "|" + gene + "|" + probeSet;

    if (m_egaIndex.count(ega) > 0) {
      APT_ERR_ABORT(markerFile  + ": duplicate experiment, gene, probeSet detected: " + experiment + ", " + gene + ", " + probeSet);

    }
    m_egaIndex[ega] = i;

    if (m_experiments.find(experiment) == m_experiments.end()) {
      m_experiments.insert(experiment);
    }

  }


}
// end DMET2MarkerModel::DMET2MarkerModel
/*****************************************************************************/
///////////////////////////////////////////////////////////////////////////////
// DMET2HaplotypeModel
///////////////////////////////////////////////////////////////////////////////

// Experiment Gene Call Call_Count Known_Count
const TittmColumnDefinition DMET2_HAPLOTYPE_MODEL_COLUMN_DEFINITIONS[] = { // column_name, empty_ok, index, valid_regex
  { std::string("Experiment"),  0,  0, std::string(".*"), NULL },
  { std::string("Gene"),        0,  1, std::string("^[A-Z\\d]+$"), NULL},
  { std::string("Call"),        0,  2, std::string(".*"), NULL },
  { std::string("Call_Count"),  0,  3, std::string("^(?:\\d+)$"), NULL },
  { std::string("Known_Count"), 0,  4, std::string("^(?:\\d+)$"), NULL },
};

/*****************************************************************************/
/**
 * DMET2HaplotypeModel::DMET2HaplotypeModel
 * Synopsis: Default constructor.
 *
 * @param rte - the runtime environment
 * @param haplotypeFile - the TsvFile with the data.
 * @return - description
 */
/*****************************************************************************/
DMET2HaplotypeModel::DMET2HaplotypeModel(RunTimeEnvironment *rte,
    const string haplotypeFile) :
    TranslationInputTsvTableModel(*rte,
                                  haplotypeFile,
                                  DMET2_HAPLOTYPE_MODEL_COLUMN_DEFINITIONS,
                                  (size_t) sizeof(DMET2_HAPLOTYPE_MODEL_COLUMN_DEFINITIONS) / sizeof(DMET2_HAPLOTYPE_MODEL_COLUMN_DEFINITIONS[0]),
                                  true)
{

  // Create a hash of the experiment, gene call for comparison with the
  // DMET3 output.

  for (int i = 0; i < m_rows.size(); i++) {

    string experiment = m_rows[i][0];
    string gene       = m_rows[i][1];
    string call       = m_rows[i][2];
    string egc        = experiment + "|" + gene + "|" + call;

    if (m_egcIndex.count(egc) > 0) {
      APT_ERR_ABORT(haplotypeFile  + ": duplicate experiment gene call detected: " + experiment + ", " + gene + ", " + call);
    }

    m_egcIndex[egc] = i;

    if (m_experiments.find(experiment) == m_experiments.end()) {
      m_experiments.insert(experiment);
    }
  }


}
// end DMET2HaplotypeModel::DMET2HaplotypeModel
/*****************************************************************************/

