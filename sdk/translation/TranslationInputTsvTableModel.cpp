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
 * @file   TranslationInputTsvTableModel.cpp
 * @author Mybrid Spalding
 * @date   Wed Apr 23 10:09:02 PDT 2008
 * @brief  Base class for the multitude of single instance models of input TsvFile or Console PgOptions abstracted as a 'vector of std::vector of std::strings' table in memory. 
 */

//
#include "translation/TranslationInputTsvTableModel.h"
//
#include "util/Err.h" // includes "util/Verbose.h"
#include "util/Util.h"
//
#include <iostream>
#include <sstream>

using namespace affx;
using namespace std;



///////////////////////////////////////////////////////////////////////////////
// BEGIN TranslationInputTsvTableModel
///////////////////////////////////////////////////////////////////////////////

/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::TranslationInputTsvTableModel
 * Synopsis:
 * Main constructor that slurps in the TsvFile and creates the in memory
 * version of same said file. Do not use this object if the file is
 * too large to fit into main memory, instead see
 * TranslationInputStreamTableModel.
 *
 * There are two callback functions passed in:
 * 1.) createDynamicColumnDefinitions, this function takes in the set of
 * static columns, or static column schema, already initialized and then
 * modifies the schema based on columns in the input file itself.
 * Dynamic columns need to be appended to the schem in order for
 * enums defined in the super class to work.
 *
 * 2.) filterRow, this function takes in the row and returns true or
 * false. If true is  returned then the row is skipped and
 * not appended to the 'm_rows' of data. 
 *
 * @param rte - the single instance run time environment (obiously not a TranslationInputTsvTableModel) 
 * @param inputTsvFileName - the TsvFile object, tab delimited. 
 * @param columnDefinitions - the schema for the table
 * @param cdSize - the number of columns in the schema (columnDefinitions)
 * @param strictTsvColumnOrder  - whether to enforce strict column count and order
 * @param createDynamicColumnDefinitions - call back function to append dynamic columns to the static schema where the dynamic columns are taken from input file itself. Primarily implemented for the A1-ANN fields in the translation table file. 
 * @param filterRow - callback function applied to each row during the slurp. 
 * @param validateColumnRegex - boolean to indicate if validation of the columns should occur. 
 *
 */
/*****************************************************************************/
TranslationInputTsvTableModel::TranslationInputTsvTableModel(
  const RunTimeEnvironment& rte,
  const std::string& inputTsvFileName,
  const TittmColumnDefinition columnDefinitions[],
  size_t cdSize,
  bool strictTsvColumnOrder,
  bool (*createDynamicColumnDefinitions)(const RunTimeEnvironment & rte, std::vector<TittmColumnDefinition> & columnDefs, const std::string & inputTsvFile), bool (*filterRow)(std::vector< std::string > & row), bool validateColumnRegex) :  m_strictTsvColumnOrder(strictTsvColumnOrder), m_validateColumnRegex(validateColumnRegex)
{

  std::vector< TittmColumnDefinition > columnDefinitions_vec;

  if (createDynamicColumnDefinitions) {
    if (!(*createDynamicColumnDefinitions)(rte, columnDefinitions_vec, inputTsvFileName)) {
      APT_ERR_ABORT(inputTsvFileName + ": errors detected.");
    }
    APT_ERR_ASSERT(columnDefinitions_vec.size() > 0, "");
  } else {
    for (int i = 0; i < cdSize; i++) {
      columnDefinitions_vec.push_back(columnDefinitions[i]);
    }
  }

  _initializeColumnDefinitions(rte, columnDefinitions_vec);

  clearData();

  readTsvFile(rte, inputTsvFileName, filterRow);

}
// end TranslationInputTsvTableModel::TranslationInputTsvTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::readTsvFile
 * Synopsis:
 * 
 * The part of the main constructor which reads in the TsvFile
 * into main memory. Public in case multiple reads are desired.
 * The data is read into the 'm_rows" table of std::strings.
 *
 * When a data error is encournted the file is fully parsed and all
 * data errors are reported on before an exception is thrown.
 *
 * Regular expression validation is supported via the regular expression
 * attached to the schema. More involved validation can be accomplished
 * via the filterRow callback function. 
 *
 * @param rte - RunTimeEnvironment
 * @param inputTsvFileName - the file in question.
 * @param filterRow = callback function, return of true ignores row, see the constructor for details.
 */
/*****************************************************************************/
void TranslationInputTsvTableModel::readTsvFile(const RunTimeEnvironment & rte,
    const std::string& inputTsvFileName,
    bool (*filterRow)(std::vector<std::string> & row))
{

  // This should never happen unless a new constructor is created
  // that doesn't initialize columns.
  APT_ERR_ASSERT(m_columnDefinition.size(), "");

  clearData();

  TsvFile tsv;

  Verbose::out(ADT_VERBOSE_INPUT_FILES, ToStr("### TranslationInputTsvTableModel::readFile('") + inputTsvFileName + ToStr(")\n"));


  //
  tsv.m_optAutoTrim   = true; // remove '"'s
  tsv.m_optQuoteChar1 = 0;    // ignore "'"s
  tsv.m_optCheckFormatOnOpen=false; // cause we accept files which have 8bit chars.

  if (tsv.open(inputTsvFileName) != TSV_OK) {
    APT_ERR_ABORT(inputTsvFileName + ": failed opening input Tsv file.");
  }

  bool okColumnDefinitions = true;

  if (m_strictTsvColumnOrder) {
    okColumnDefinitions =
      _matchTsvColumnsWithColumnDefinitions(inputTsvFileName, tsv);
  }

  if (okColumnDefinitions) {
    okColumnDefinitions = _initializeTsvColumns(inputTsvFileName, tsv);
  }

  int columnCount = tsv.getColumnCount(0);
  bool okColumnValues = true;
  int row = 1;
  std::string whiteSpace = " ";

  while (okColumnDefinitions && (tsv.nextLevel(0) == TSV_OK)) {
    //
    std::vector<std::string> newRow;
    std::string columnValue;
    row++;

    newRow.resize(columnCount, "");
    bool okRow = true;
    for (int tsvColumn = 0; tsvColumn < columnCount; tsvColumn++) {

      std::stringstream lineSStr;
      std::string       columnName = getColumnName(m_tsvIndexToColumnDefinitionIndex[tsvColumn]);
      // get the data

      if (tsv.get(0, tsvColumn, columnValue) != TSV_OK) {
        APT_ERR_ABORT("tsv.get != TSV_OK");
      }

      Util::trimString(columnValue, whiteSpace.c_str());

      // empty and allowed?
      if (columnValue.empty()) {
        if (m_columnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_emptyOk == false) {
          lineSStr << inputTsvFileName << ": ROW " << row  << " : ";
          lineSStr << columnName << " : ERROR missing required column.";
          Verbose::out(rte.m_currentVerbosity, lineSStr.str());
          okColumnValues = false;
        }
      }
      // validate it against regexp
      else if (m_validateColumnRegex && !m_columnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_validRE.FullMatch(columnValue)) {
        lineSStr << inputTsvFileName << ": COLUMN [";
        lineSStr << tsvColumn << "] " << columnName << " : ERROR invalid value: \"";
        lineSStr << columnValue << "\"";
        lineSStr << " (expecting " << m_columnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_validRE.pattern();
        lineSStr << ")";
        lineSStr << endl;

        for (int j = 0 ; j < columnCount;  j++) {
          tsv.get(0, j, columnValue);
          lineSStr << columnValue << ",";
        }
        Verbose::out(rte.m_currentVerbosity, lineSStr.str());
        okColumnValues = false;
      }
      // Ingore row based on data
      else if ((m_columnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_ignoreRE != NULL) &&
               (m_columnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_ignoreRE->PartialMatch(columnValue))) {
        okRow = false;

      }
      newRow[m_tsvIndexToColumnDefinitionIndex[tsvColumn]] = columnValue;
    }

    // Ingore row based on filter
    if ((filterRow != NULL) && (*filterRow)(newRow)) {
      okRow = false;
    }

    // save our row to the data.
    if (okRow) {
      m_fileRowIndex[m_rows.size()] = row;
      m_rows.push_back(newRow);
    }
  }
  //

  tsv.close();

  if (! okColumnDefinitions || ! okColumnValues) {
    APT_ERR_ABORT(inputTsvFileName + ": invalid input file or data detected.");
  }

}
// end TranslationInputTsvTableModel::readTsvFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::getColumnName
 * Synopsis:
 *
 * Get the column name from the Transltation Table Tsv File. Used
 * for dynamic columns not defined in the static schema.
 *
 * @param column - the name of the 0 based column index column in question.
 */
/*****************************************************************************/
std::string TranslationInputTsvTableModel::getColumnName(int column)
{

  // This should always be true unless some constructor is created
  // that doesn't initialize this.
  APT_ERR_ASSERT(m_columnDefinition.size(), "");

  APT_ERR_ASSERT((column >= 0) && (column < columnCount()), "");

  return m_columnDefinition[column].m_columnName;

}
// end TranslationInputTsvTableModel::getColumnName
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::getHeaderAsString
 * Synopsis:
 *
 * Debug routine that returns the table schema, header line in the
 * TSV file, as a comma (delimiter) separated list.
 *
 * @param delimiter - comma by default
 *
 * @return header - the list of column names
 */
/*****************************************************************************/
std::string TranslationInputTsvTableModel::getHeaderAsString(const std::string delimiter)  const
{

  std::stringstream headerSStr;

  for (int j = 0 ; j < columnCount(); j++) {
    headerSStr << m_columnDefinition[j].m_columnName;

    if (j < (columnCount() - 1)) {
      headerSStr << delimiter;
    }

  } // for each column

  return headerSStr.str();

}
// end TranslationInputTsvTableModel::getHeaderAsString
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::getHeaderAsVector
 * Synopsis:
 * 
 * Debug conveinance routine that returns the schema as a std::vector of std::strings. 
 *
 * @return std::vector<std::string> - the column names as a std::vector
 */
/*****************************************************************************/
std::vector< std::string >  TranslationInputTsvTableModel::getHeaderAsVector() const
{

  std::vector< std::string > header;

  for (int j = 0 ; j < columnCount(); j++) {
    header.push_back(m_columnDefinition[j].m_columnName);
  } // for each column

  return header;

}
// end TranslationInputTsvTableModel::getHeaderAsVector
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::getRowAsString
 * Synopsis:
 * 
 * Debug convienance routine that std::stringifies the row. Also used in
 * error messages. All columns are "separated" by a single comma (or passed
 * in delimiter).
 * If the data within a column contains a column there is no attempt made
 * quote the columns. Empty columns will have comma sequences as ",,,". 
 * column 1 2   3 
 * data   A B,C D => "A,B,C,D"
 *
 * The delimiter option can be used if the data contains commas but another
 * delimiter would not, so as to  insure the count of delimiters
 * is the column count. 
 *
 * @param row       - the index into 'm_rows'
 * @param delimiter - comma by default.
 *
 * @return std::string - the comma (delimiter) separated std::string
 */
/*****************************************************************************/
std::string TranslationInputTsvTableModel::getRowAsString(const int row, const std::string delimiter)  const
{

  APT_ERR_ASSERT((row >= 0) && (row < size()), "");

  std::stringstream rowSStr;

  for (int j = 0 ; j < columnCount(); j++) {
    rowSStr << m_rows[row][j];

    if (j < (columnCount() - 1)) {
      rowSStr << delimiter;
    }

  } // for each column

  return rowSStr.str();

}
// end TranslationInputTsvTableModel::getRowAsString
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::columnCount
 * Synopsis:
 * 
 * Interogate an already instantiated table for the number of columns.
 * This method assumes ALL rows have the same column width and just
 * returns the width of the column defintion structure. .
 *
 * @return - the number of columns in a row. 
 */
/*****************************************************************************/
int TranslationInputTsvTableModel::columnCount() const
{

  if (m_columnDefinition.size()  > 0) {
    return m_columnDefinition.size();
  }

  return 0;

}
// end TranslationInputTsvTableModel::columnCount
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::_initializeColumnDefinitions
 * Synopsis:
 * Helper function to the constructor,
 * TranslationInputTsvTableModel::TranslationInputTsvTableModel
 *
 * Initialize the definition of columns data structure. This is compared
 * against the input TSV file when it is opened for integrity. 
 *
 * @param rte - the single instance RunTimeEnvironment
 * @param columnDefinitions - the table schema 
 *
 */
/*****************************************************************************/
void TranslationInputTsvTableModel::_initializeColumnDefinitions(const RunTimeEnvironment & rte, const std::vector< TittmColumnDefinition> columnDefinitions)
{

  APT_ERR_ASSERT(columnDefinitions.size() > 0, "");

  m_columnDefinition.clear();

  // Iterate so we can doublecheck the m_index value.
  for (int i = 0; i < columnDefinitions.size(); i++) {

    // Ok, columnDefintions are created at compile time
    // and are NOT dynamic at run-time, therefore we assert to the following
    // is true.

    APT_ERR_ASSERT(columnDefinitions[i].m_index == i, "");

    m_columnNameIndex[columnDefinitions[i].m_columnName] = i;
    m_columnDefinition.push_back(columnDefinitions[i]);
  }

  return;
}
// end TranslationInputTsvTableModel::_initializeColumnDefinitions
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::_initializeTsvColumns
 * Synopsis:
 *
 * Helper function to TranslationInputTsvTableModel::readTsvFile
 *
 * Initialize the TsvFile object columns from the column definitions
 * and validate against the table schema. 
 *
 * @param inputTsvFileName - The name of the TsvFile in question.
 * @param tsv              - The TsvFile object in question.
 *
 * @return false - on error.
 *
 */
/*****************************************************************************/
bool TranslationInputTsvTableModel::_initializeTsvColumns(const std::string & inputTsvFileName, affx::TsvFile & tsv)
{

  bool         okTsvColumns = true;
  std::stringstream msgSStr;

  m_tsvIndexToColumnDefinitionIndex.resize(tsv.getColumnCount(0), -1);


  for (int column = 0; column < m_columnDefinition.size(); column++) {

    // find the tsv column index for this
    int tsvColumn
    = tsv.cname2cidx(0, m_columnDefinition[column].m_columnName);

    if (tsvColumn < 0) {
      okTsvColumns = false;
      msgSStr << inputTsvFileName << ": invalid TsvFile header first-line, ";
      msgSStr << "column [" << (column + 1) << "], \"";
      msgSStr << m_columnDefinition[column].m_columnName << "\"  not found." << endl;
      continue;
    }

    m_tsvIndexToColumnDefinitionIndex[tsvColumn]
    = m_columnDefinition[column].m_index;

  }

  if (!okTsvColumns) {
    Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str(), false);
  }

  return okTsvColumns;

}
// end TranslationInputTsvTableModel::_initializeTsvColumns
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputTsvTableModel::_matchTsvColumnsWithColumnDefinitions
 * Synopsis:
 *
 * Helper function to TranslationInputStreamTableModel::_openTsvFile
 *
 * It is not normal Affymetrix practice to require to have only the
 * set of columns required nor is it normal practice to have those
 * columns in any order. The TsvFile object can access columns by
 * name independent of order. However, there exists requirements
 * for this application to enforce column order as well as column order
 * enforce only columns diseried.
 * 
 * This helper is used to enforce strict column ordering when required.
 * Inspect the first line of the TsvFile and insure the appropriate header
 * column names exist one-to-one and in the same order
 * with all the column definitions. Note this is only required for certain
 * files such as the human generated TranslationTable.xls spreadsheet.
 *
 * @param inputTsvFileName - the Tsv file in question
 * @param tsv - the TsvFile object that is opened and corresponds to
 *              inputTsvFileName
 *
 * @return false - the header line is invalid.
 */
/*****************************************************************************/
bool TranslationInputTsvTableModel::_matchTsvColumnsWithColumnDefinitions(const std::string & inputTsvFileName, affx::TsvFile & tsv)
{

  bool okHeaderLine = true;

  int invalidHeaderCount = 0;

  std::stringstream msgSStr;

  for (int i = 0; (i < tsv.getColumnCount(0)) &&
       (i < m_columnDefinition.size()) ; i++) {

    std::string columnName;

    if (tsv.cidx2cname(0, i, columnName) == TSV_OK) {

      pcrecpp::RE columnRE(m_columnDefinition[i].m_columnName);

      if (! columnRE.FullMatch(columnName)) {
        okHeaderLine = false;
        msgSStr << inputTsvFileName << ": invalid TsvFile header first-line, ";
        msgSStr << "column [" << (i + 1) << "]: found \"" << columnName << "\" but was expecting \"";
        msgSStr << m_columnDefinition[i].m_columnName << "\"." << endl;
        invalidHeaderCount++;
      }
    } else {

      APT_ERR_ABORT("tsv.cidx2cname returned unexpected error");

    } // if there was a TSV_ERROR

  } // for each column, inspect it

  if (invalidHeaderCount > 0) {
    if (invalidHeaderCount < 4 || true) {
      Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str(), false);
    } else {
      Verbose::out(ADT_VERBOSE_NORMAL, inputTsvFileName + ": invalid TsvFile is missing the header line.", false);
    }
  }

  return okHeaderLine;

}
// end TranslationInputTsvTableModel::_matchTsvColumnsWithColumnDefinitions
/*****************************************************************************/
