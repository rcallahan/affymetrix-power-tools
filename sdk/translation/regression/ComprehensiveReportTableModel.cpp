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
 * @file   ComprehensiveReportTableModel.cpp
 * @author Mybrid Spalding
 * @date   Fri Nov 14 14:16:02 PST 2008
 * @brief  Comprehensive report wrapper for testing purposes.
 */


#include "translation/regression/ComprehensiveReportTableModel.h"
//
#include "translation/CallElement.h"
#include "translation/RunTimeEnvironment.h"
#include "translation/TranslationTableModel.h"
//
#include "util/Err.h"
//
#include "pcrecpp.h"
//
#include <sstream>
//

using namespace affx;
using namespace std;

const TittmColumnDefinition CRTM_COLUMN_DEFINITIONS[] = {
  {std::string("Index"),              1, 0, std::string(".*"), NULL },
  {std::string("CHP Filename"),       1, 1, std::string(".*"), NULL},
  {std::string("Gene"),               1, 2, std::string(".*"), NULL },
  {std::string("Known Call"),         1, 3, std::string(".*"), NULL },
  {std::string("Unknown Call"),       1, 4, std::string(".*"), NULL},
  {std::string("Interpretation Code"),1, 5, std::string(".*"), NULL},
  {std::string("Summary Flag"),       1, 6, std::string(".*"), NULL },
  {std::string("Relevant Alleles"),   1, 7 , std::string(".*"), NULL},
  {std::string("Common Name"),        1, 8, std::string(".*"), NULL },
  {std::string("Probe Set ID"),       1, 9, std::string(".*"), NULL }, 
  {std::string("Basecall"),           1, 10, std::string(".*"), NULL },
  {std::string("Reference Base"),     1, 11, std::string(".*"), NULL },
  {std::string("Variant Base"),       1, 12, std::string(".*"), NULL },
  {std::string("Call"),               1, 13, std::string(".*"), NULL },
  {std::string("Haplotype Marker"),   1, 14, std::string(".*"), NULL },
  {std::string("Change for Variant"), 1, 15, std::string(".*"), NULL },
  {std::string("cDNA Change"),        1, 16, std::string(".*"), NULL },
  {std::string("Genome Position"),    1, 17, std::string(".*"), NULL },
  {std::string("dbSNP RS ID"),        1, 18, std::string(".*"), NULL },
  {std::string("Original Basecall"),  1, 19, std::string(".*"), NULL },
  {std::string("Override Comment"),   1, 20, std::string(".*"), NULL },
};

/*****************************************************************************/
/**
 * ComprehensiveReportTableModel::ComprehensiveReportTableModel
 * Synopsis:
 *
 *
 * @param rte - RunTimeEnvironment which contains options
 * @param reportFileName - the input file to read in.
 * @param ttm - The translation table instance used for validation.
 */
/*****************************************************************************/
ComprehensiveReportTableModel::ComprehensiveReportTableModel(const RunTimeEnvironment &rte,  const string &reportFileName) :
    TranslationInputTsvTableModel(rte, reportFileName, CRTM_COLUMN_DEFINITIONS, (sizeof(CRTM_COLUMN_DEFINITIONS) / sizeof(CRTM_COLUMN_DEFINITIONS[0])), true, NULL, NULL, false)
{




}
// end ComprehensiveReportTableModel::ComprehensiveReportTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * ComprehensiveReportTableModel::getColumnName
 * Synopsis:
 *
 *
 */
/*****************************************************************************/
std::string ComprehensiveReportTableModel::getColumnName(int col)  {

  assert( col < m_rows[0].size() );
  
  return (CRTM_COLUMN_DEFINITIONS[col].m_columnName);
  
}
// end ComprehensiveReportTableModel::getColumnName
/*****************************************************************************/
/*****************************************************************************/
/**
 * ComprehensiveReportTableModel::getColumnIndex
 * Synopsis:
 *
 *
 */
/*****************************************************************************/
int ComprehensiveReportTableModel::getColumnIndex(std::string columnName)  {

  if ( m_columnIndex.count(columnName) != 0 ) {
    return m_columnIndex[columnName];
  }

  for ( int col = 0; col < m_rows[0].size(); col++ ) {
    if ( CRTM_COLUMN_DEFINITIONS[col].m_columnName == columnName ) {
      m_columnIndex[columnName] = col;
      return col;
    }
  }

  APT_ERR_ABORT("Programmer Error: invalid column name given for comprehensive report.");
  
  
  return -1;
  
}
// end ComprehensiveReportTableModel::getColumnIndex
/*****************************************************************************/
