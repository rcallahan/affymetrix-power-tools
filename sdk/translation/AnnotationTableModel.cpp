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
 * @file   AnnotationTableModel.cpp
 * @author Mybrid Spalding
 * @date   Thu Oct 30 15:39:56 2008
 *
 * @brief  Audit input tranlation and annoatation files for errors.
 *
 *
 */


#include "translation/AnnotationTableModel.h"
//
#include <cstring>
#include <sstream>
#include <string>
#include <vector>
//


const TittmColumnDefinition ANNOTATION_TABLE_COLUMN_DEFINITIONS[] = {
  // m_columnName,        m_emptyOk, m_index, m_validRE
  { std::string("Probe Set ID"),    0, 0,  std::string("^[\\w\\d]+$"), NULL },
  { std::string("Associated Gene"), 0, 1, std::string("^(?:[\\/\\s_A-Z0-9]+)|(?:N\\/A)$"), NULL },
  { std::string("dbSNP RS ID"),  1, 2,
    std::string("^(?:rs\\d+(?:[|;\\s])*)+|(?:N/?A)$"), NULL },
  { std::string("Common Name"), 1, 3, std::string(".*") , NULL },
  { std::string("Chromosome"), 0, 4, std::string("^(?:\\d+|[Xx]|[yY])$"), NULL },
  { std::string("Physical Position"), 0, 5, std::string("^\\d+$"), NULL },
  { std::string("Genome Context"), 0, 6, std::string("^[\\-\\<\\>a-zA-Z\\[\\]\\/]+$"), NULL },
  { std::string("Design Strand"), 0, 7, std::string("^-?1$"), NULL },
  { std::string("Alleles Design Strand"), 1, 8, std::string("^[\\s\\/ACTGN\\-]+$"), NULL },
  { std::string("Allele Code"),  0, 9, std::string("^[\\s\\/A-Z\\-]+$"), NULL },
  { std::string("Gene Strand"), 1, 10, std::string("^(?:\\+)|(?:-)|(?:N\\/?A)$"), NULL },
  { std::string("Switch Design Strand to Report"), 0, 11, std::string("^[01]|NA$"), NULL },
  { std::string("Alleles Reported Strand"), 0, 12,  std::string("^[\\s\\/A-Z\\-]+$"), NULL },
  { std::string("Reported Strand"), 1, 13, std::string("^\\+|-|N\\/?A$"), NULL },
  { std::string("Type"), 0, 14, std::string("^snp|in-del|tiling|copy_number$"), NULL },
  { std::string("Validated"), 0, 15, std::string("^N|Y|N\\/?A$"), NULL },
  { std::string("DbSNP annot"), 1, 16, std::string(".*"), NULL },
  { std::string("PharmGKB Code"), 1, 17, std::string("^N\\/A|PA\\d+$"), NULL },
  { std::string("Alleles-Alias Reported Strand"), 1, 18, std::string(".*"), NULL },
  { std::string("Context Code"), 0, 19,  std::string("^(?:N\\/A)|(?:[\\s\\/\\d]+)$"), NULL },
  { std::string("Context Sequence"), 0, 20, std::string("^[\\-\\s\\[\\]\\<\\>ACGNTacgt\\/]+$"), NULL },
};




/*****************************************************************************/
/**
 * AnnotationTableModel::AnnotationTableModel
 * Synopsis:
 * Main constructor to read the probe set filter file and stash the data
 * into the model.
 *
 * @param rte - RunTimeEnvironment which contains options
 * @param probeSetFilterFileName - the input file to read.
 */
/*****************************************************************************/
AnnotationTableModel::AnnotationTableModel(const RunTimeEnvironment &rte,
    const std::string & annotationFileName) :
    TranslationInputTsvTableModel(rte, annotationFileName,
                                  ANNOTATION_TABLE_COLUMN_DEFINITIONS,
                                  (sizeof(ANNOTATION_TABLE_COLUMN_DEFINITIONS) / sizeof(ANNOTATION_TABLE_COLUMN_DEFINITIONS[0])), false, NULL)
{

  for (int row = 0; row < size(); row++) {

    m_probeSetIndex[m_rows[row][PROBE_SET_ID]].push_back(row);

  }


}
// end AnnotationTableModel::AnnotationTableModel
/*****************************************************************************/
