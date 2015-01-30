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
 * @file   SampleInfoTableModel.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  A single instance class that wraps sample information either passed in by the console as a std::vector or when read from a file. 
 */


#include "translation/SampleInfoTableModel.h"
//
#include <sstream>

using namespace affx;
using namespace std;

const TittmColumnDefinition SITM_COLUMN_DEFINITIONS[] = {
// experiment sample prediction notes
// 16x6_NA18507_Rep_1 NA18507 1 CN > 0: that is, CN=1,2,3,...
  // column name,         empty?, column, valid regular expersion
  { std::string("Experiment Id"), 0,     0,  std::string("^[\\w\\d-]+$"), NULL },
};

const TittmColumnDefinition SAMPLE_INFO_COLUMN_TEMPLATE = { std::string("NULL"), 1, 0, std::string(".*"), NULL };

/*****************************************************************************/
/**
 * convertTCDArrayToVector
 * Synopsis:
 * 
 * Callback function required for TranslationInputTsvFile constructor
 * to generate the dynamic column definitions for SITM_COLUMN_DEFINTIONS.
 * There is only one fixed column, the first field, GUID. This field
 * is used to join with the CHP guid to marry up the correct sample
 * info record with the corresponding CHP file.
 *
 * @param rte                - the single instance RunTimeEnvironment
 * @param sitmTCD            - the model used to represent column meta data
 * @param sampleInfoFileName - the TSV file name
 *
 * @return - true if the file is ok.
 *
 */
/*****************************************************************************/
static bool convertTCDArrayToVector(const class RunTimeEnvironment & rte, std::vector< TittmColumnDefinition > & sitmTCD, const std::string & sampleInfoFileName)
{

  sitmTCD.clear();

  std::stringstream msgSStr;

  size_t       numStaticDefinitions
  = sizeof(SITM_COLUMN_DEFINITIONS) / sizeof(SITM_COLUMN_DEFINITIONS[0]);

  bool         okHeaderLine = true;

  // STATIC column definitions
  for (int i = 0; i < numStaticDefinitions; i++) {

    APT_ERR_ASSERT(SITM_COLUMN_DEFINITIONS[i].m_index == i, "");
    sitmTCD.push_back(SITM_COLUMN_DEFINITIONS[i]);
    sitmTCD.back().m_ignoreRE = new pcrecpp::RE("^\\#");

  }

  // The sample info columns are dynamic and not fixed.
  // Scan the TsvFile headers and pick up all the columns.

  TsvFile tsv;

  //
  tsv.m_optAutoTrim   = true; // remove '"'s
  tsv.m_optQuoteChar1 = 0;    // ignore "'"s

  if (tsv.open(sampleInfoFileName) != TSV_OK) {
    APT_ERR_ABORT(sampleInfoFileName + ": failed opening input Tsv file.");
  }

  // DYNAMIC column definitions.
  for (int i = numStaticDefinitions; i < tsv.getColumnCount(0) ; i++) {

    std::string columnName;

    if (tsv.cidx2cname(0, i, columnName) == TSV_OK) {

      TittmColumnDefinition sampleTCD = SAMPLE_INFO_COLUMN_TEMPLATE;
      sampleTCD.m_index = i;
      sampleTCD.m_columnName = columnName;
      sitmTCD.push_back(sampleTCD);

    } // if the tsv file parsed ok
    else {
      APT_ERR_ABORT("tsv.cidx2cname returned unexpected error");
    } // if there was a TSV_ERROR


  } // for each column, append it to the column definitions.

  if (sitmTCD.size() == numStaticDefinitions) {
    okHeaderLine = false;
    msgSStr << sampleInfoFileName << ": invalid TsvFile header first-line, ";
    msgSStr << "no columns found." << endl;

  }

  tsv.close();

  if (!okHeaderLine) {
    if (sitmTCD.size() > 5) {
      Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str(), false);
    } else {
      Verbose::out(ADT_VERBOSE_NORMAL, sampleInfoFileName + ": invalid TsvFile is missing the header line.", false);
    }
  }

  return okHeaderLine;

}
// end convertTCDArrayToVector
/*****************************************************************************/
/*****************************************************************************/
/**
 * SampleInfoTableModel::SampleInfoTableModel
 * Synopsis:
 *
 * Command line constructor to read the sample info report file and
 * stash the data into the TranslationInputTsvTableModel.
 *
 * @param rte                - single instance RunTimeEnvironment 
 * @param sampleInfoFileName - the input file to read.
 */
/*****************************************************************************/
SampleInfoTableModel::SampleInfoTableModel(const class RunTimeEnvironment &rte,
    const std::string & sampleInfoFileName) :
    TranslationInputTsvTableModel(rte, sampleInfoFileName, SITM_COLUMN_DEFINITIONS,
                                  (sizeof(SITM_COLUMN_DEFINITIONS) / sizeof(SITM_COLUMN_DEFINITIONS[0])), false,  &convertTCDArrayToVector),    m_sampleInfoFileName(sampleInfoFileName)
{



}
// end SampleInfoTableModel::SampleInfoTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * SampleInfoTableModel::SampleInfoTableModel
 * Synopsis:
 *
 * Console constructor to coerce the sample info from a PgOpts
 * std::vector into the TSV model TranslationInputTsvTableModel.
 * 
 *
 * @param rte - the single instance RunTimeEnvironment which contains the
 * std::vector in rte.m_adtOpts read from PgOpts:
 * std::vector< std::vector < std::string > > m_sampleTable
 *
 */
/*****************************************************************************/
SampleInfoTableModel::SampleInfoTableModel(const class RunTimeEnvironment &rte)
{


  int numDynamicRows = rte.m_adtOpts.m_sampleTable.size();

  if (numDynamicRows == 0) {
    return;
  }

  // HEADER

  for (int i = 0; i < rte.m_adtOpts.m_sampleTable[0].size() ; i++) {

    TittmColumnDefinition sampleTCD = SAMPLE_INFO_COLUMN_TEMPLATE;
    sampleTCD.m_index = i;
    sampleTCD.m_columnName = rte.m_adtOpts.m_sampleTable[0][i];
    m_columnDefinition.push_back(sampleTCD);
  }


  // DATA
  for (int i = 1; i < numDynamicRows; i++) {
    m_rows.push_back(rte.m_adtOpts.m_sampleTable[i]);
  }


  return;

}
// end SampleInfoTableModel::SampleInfoTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * SampleInfoTableModel::getSampleInfo
 * Synopsis:
 *
 *  Given the the experiment id (name) for
 *  geno data (CHP or Genotype Short Report),
 *  returns the corresponding sample row.
 *
 *
 * @param experimentId - experiment name, CHP file name or as report in
 * the Genotype Short Report. 
 *
 * @return - a std::vector< std::string > with the experiment's sample info
 */
/*****************************************************************************/
std::vector< std::string > SampleInfoTableModel::getSampleInfo(const std::string & experimentId)
{

  std::vector< std::string > results;

  if (experimentId.empty()) {
    return results;
  }

  //cerr << "getSampleInfo: " << experimentId << endl;

  for (int row = 0; (row < m_rows.size()) && (results.size() == 0); row++) {
    if (m_rows[row][0] == experimentId) {
      for (int i = 1; i < m_rows[row].size(); i++) {
        results.push_back(m_rows[row][i]);
      }
    }
  }

  return results;
}
// end SampleInfoTableModel::getSampleInfo
/*****************************************************************************/

