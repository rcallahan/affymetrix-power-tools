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
 * @file   TranslationTableModel.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  Single instance data model class that wraps the translation file in the TranslationInputTsvModel.
 */


#include "translation/TranslationTableModel.h"
//
#include "translation/CallSet.h"
//
#include "calvin_files/utils/src/GenoCallCoder.h"
//
#include <iostream>
#include <list>
#include <set>
#include <sstream>
//

using namespace std;
using namespace affx;


//Gene Reference Link ProbeSet Switch Design Strand dbSNP Defining cDNA Genomic Change External ID Validated Haplotype Reference Variant A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16 A17 A18 A19 A20 A21 A22 A23 A24 A25 A26 A27 A28 A29
//ABCB1 Reference Link ProbeSet refSNP ID Defining cDNA Nucleotide Position Genomic Nucleotide Position Change External ID Validated Haplotype Reference Variant WT *14 A893S A893T  N44S I829V S400N A599T V801M

pcrecpp::RE RE_IGNORE_ROW("^\\s*[#]");

const TittmColumnDefinition TTM_DMET2_COLUMN_DEFINTIONS[] = { // m_columnName,        m_emptyOk, m_index, m_validRE
  { std::string("Gene"),           0,  0, std::string("^[A-Z\\d]+$"), NULL },
  { std::string("Reference Link"), 1,  1, std::string(".*"), NULL },
  { std::string("AssayID"),        0,  2,
    std::string("^(?:AssayID)|(?:[\\w\\d]+)|(?:\\#.*)$"), NULL },
  { std::string("dbSNP"),          1,  3,
    std::string("^(?:(?:rs\\d+[|;\\s]*)+|(?:N/?A)|(?:refSNP\\sID)|(?:dbSNP))$"),
    NULL },
  { std::string("Defining"),       1,  4,
    std::string("^(?:(?:N)|(?:Defining)|(?:\\*?[>A-Za-z\\d\\-\\s\\+]+))$"), NULL },
  { std::string("cDNA"),           1,  5, std::string(".*"), NULL },
  { std::string("Genomic"),        1,  6, std::string(".*"), NULL },
  { std::string("Change"),         1,  7, std::string(".*"), NULL },
  { std::string("External ID"),    1,  8, std::string(".*"), NULL },
  { std::string("Validated"),      1,  9,
    std::string("^(?:Validated)|(?:[YN])$"), NULL },
  { std::string("Haplotype"),      1, 10,  std::string("(?:Haplotype)|(?:[YN])"), NULL },
  { std::string("Reference"),      1, 11,
    std::string("^(?:Reference)|(?:Minor Allele)|(?:Major Allele)|(?:[ACGT0]+)|(?:INS)|(?:DEL)|(?:-)$"), NULL },
  { std::string("Variant"),        1, 12,
    std::string("^(?:Variant)|(?:Minor Allele)|(?:Major Allele)|(?:[ACGT0]+)|(?:INS)|(?:DEL)|(?:-)$"), NULL },
};

const TittmColumnDefinition TTM_DMET3_COLUMN_DEFINTIONS[] = { // m_columnName,        m_emptyOk, m_index, m_validRE
  { std::string("Gene"),           0,  0, std::string("^[A-Z\\d]+$"), NULL },
  { std::string("Reference Link"), 1,  1, std::string(".*"), NULL },
  { std::string("Probe Set ID"),        1,  2,
    std::string("^(?:Probe\\sSet\\sID)|(?:[_\\-\\w\\d]*\\d+)|(?:\\#.*)$"), NULL },
  { std::string("Switch Design Strand to Report"), 1,  3, std::string("(?:Switch.Design.Strand.to.Report)|(?:[YN])|(?:\\#.*)$"), NULL },
  { std::string("dbSNP RS ID"),          1,  4,
    std::string("^(?:(?:rs\\d+[|;\\s]*)+|(?:N/?A)|(?:refSNP\\sID)|(?:dbSNP.RS.ID))$"),
    NULL },
  { std::string("Defining"),       1,  5,
    std::string("^(?:(?:N)|(?:Defining)|(?:\\*?[>A-Za-z\\d\\-\\s\\+]+))$"), NULL },
  { std::string("cDNA Nucleotide Position"), 1,  6, std::string(".*"), NULL },
  { std::string("Genome Position"),        1,  7, std::string(".*"), NULL },
  { std::string("Change"),         1,  8, std::string(".*"), NULL },
  { std::string("Common Name"),    1,  9, std::string(".*"), NULL },
  { std::string("Haplotype"),      1, 10,  std::string("(?:Haplotype)|(?:[YN])"), NULL },
  { std::string("Reference"),      1, 11,
    std::string("^(?:Reference)|(?:Minor Allele)|(?:Major Allele)|(?:[ACGT0]+)|(?:INS)|(?:DEL)|(?:-)$"), NULL },
  { std::string("Variant"),        1, 12,
    std::string("^(?:Variant)|(?:Minor Allele)|(?:Major Allele)|(?:[ACGT0]+)|(?:INS)|(?:DEL)|(?:-)$"), NULL },
};


const TittmColumnDefinition ALLELE_COLUMN_TEMPLATE = { std::string("ANULL"), 1, 0, std::string(".*"), NULL };


class Dmet3ToDmet2ColumnMap
{
public:
  ADT_DMET3_TT_COLUMNS_ENUM m_dmet3Index;
  ADT_DMET2_TT_COLUMNS_ENUM m_dmet2Index;
};

const Dmet3ToDmet2ColumnMap DMET3_TO_DMET2_COLUMN_MAP[] = {
  { ADT_DMET3_TT_GENE,  ADT_DMET2_TT_GENE },
  { ADT_DMET3_TT_REFERENCE_LINK,  ADT_DMET2_TT_REFERENCE_LINK },
  { ADT_DMET3_TT_PROBE_SET_ID,  ADT_DMET2_TT_ASSAYID },
  { ADT_DMET3_TT_DBSNP,  ADT_DMET2_TT_DBSNP },
  { ADT_DMET3_TT_DEFINING,  ADT_DMET2_TT_DEFINING },
  { ADT_DMET3_TT_CDNA, ADT_DMET2_TT_CDNA },
  { ADT_DMET3_TT_GENOME_POSITION,  ADT_DMET2_TT_GENOMIC },
  { ADT_DMET3_TT_CHANGE,  ADT_DMET2_TT_CHANGE },
  { ADT_DMET3_TT_COMMON_NAME,  ADT_DMET2_TT_EXTERNAL_ID },
  { ADT_DMET3_TT_HAPLOTYPE, ADT_DMET2_TT_HAPLOTYPE },
  { ADT_DMET3_TT_REFERENCE, ADT_DMET2_TT_REFERENCE },
  { ADT_DMET3_TT_VARIANT, ADT_DMET2_TT_VARIANT },
  { ADT_DMET3_TT_ALLELE_START, ADT_DMET2_TT_ALLELE_START },
};

/*****************************************************************************/
/**
 * convertTCDArrayToVector
 * Synopsis:
 
 * Callback function for TranslationInputTsvFile constructor
 * to generate the dynamic allele column definitions for TTM_COLUMN_DEFINTIONS.
 * The reason being is that the
 * number of Allele columns is dynamic and determined by the input file
 * at run time and not before.
 * 
 * @param rte - the single instance run time environment
 * @param ttmTCD - the already filled in std::vector to build on the static plus dynamic stuff.
 * @param ttableFileName - the TsvFile for the Translation Table
 *
 * @return true - if ok, parent should throw an exception on false.
 *
 */
/*****************************************************************************/
static bool convertTCDArrayToVector(const RunTimeEnvironment & rte, std::vector< TittmColumnDefinition > & ttmTCD, const std::string & ttableFileName)
{


  APT_ERR_ASSERT(rte.m_adtOpts.m_inputTTableType != ADT_TRANSLATION_TABLE_TYPE_INVALID, "");

  bool isDmet2 = (rte.m_adtOpts.m_inputTTableType == ADT_TRANSLATION_TABLE_TYPE_DMET2);

  ttmTCD.clear();

  std::stringstream msgSStr;
  size_t       numStaticDefinitions
  = isDmet2 ? (size_t) sizeof(TTM_DMET2_COLUMN_DEFINTIONS) / sizeof(TTM_DMET2_COLUMN_DEFINTIONS[0]) : (size_t) sizeof(TTM_DMET3_COLUMN_DEFINTIONS) / sizeof(TTM_DMET3_COLUMN_DEFINTIONS[0]);
  bool         okHeaderLine = true;

  if (isDmet2) {
    Verbose::out(ADT_VERBOSE_INPUT_FILES, "DMET2 translation file format detected.");
  } else {
    Verbose::out(ADT_VERBOSE_INPUT_FILES, "DMET3 translation file format detected.");
  }

  // STATIC column definitions
  for (int i = 0; i < numStaticDefinitions; i++) {

    if (isDmet2) {
      APT_ERR_ASSERT(TTM_DMET2_COLUMN_DEFINTIONS[i].m_index == i, "");
      ttmTCD.push_back(TTM_DMET2_COLUMN_DEFINTIONS[i]);
    } else {
      APT_ERR_ASSERT(TTM_DMET3_COLUMN_DEFINTIONS[i].m_index == i, "");
      ttmTCD.push_back(TTM_DMET3_COLUMN_DEFINTIONS[i]);
    }
    ttmTCD.back().m_ignoreRE = new pcrecpp::RE("^\\#");
  }

  // The allele columns are dynamic and not fixed.
  // Scan the TsvFile headers and pick up the columns that match A\d+.

  TsvFile tsv;

  //
  tsv.m_optAutoTrim   = true; // remove '"'s
  tsv.m_optQuoteChar1 = 0;    // ignore "'"s

  if (tsv.open(ttableFileName) != TSV_OK) {
    APT_ERR_ABORT(ttableFileName + ": failed opening input Tsv file.");
  }

  pcrecpp::RE re("^A[\\d]+$");

  int numAlleleColumns = 0;

  // DYNAMIC Allele column definitions.
  for (int i = numStaticDefinitions; i < tsv.getColumnCount(0) ; i++) {

    std::string columnName;
    bool okAlleleColumn = true;

    if (tsv.cidx2cname(0, i, columnName) == TSV_OK) {

      if (re.FullMatch(columnName)) {

        int alleleNum = atoi(columnName.substr(1).c_str());

        // Columns must be in consecutive order A1, A2, ... AN
        if (alleleNum == (i - numStaticDefinitions + 1)) {
          numAlleleColumns++;
          TittmColumnDefinition alleleTCD = ALLELE_COLUMN_TEMPLATE;
          alleleTCD.m_index = i;
          alleleTCD.m_columnName = columnName;
          ttmTCD.push_back(alleleTCD);

        } else {
          okAlleleColumn = false;
        } // if the column count matches.
      } else {
        okAlleleColumn = false;
      } // if the column name matches
    } // if the tsv file parsed ok
    else {
      APT_ERR_ABORT("tsv.cidx2cname returned unexpected error");
    } // if there was a TSV_ERROR

    if (!okAlleleColumn) {
      okHeaderLine = false;
      msgSStr << ttableFileName << ": invalid TsvFile header first-line, ";
      msgSStr << "column [" << (i + 1) << "]: found \"" << columnName << "\" but was expecting a consecutive Allele column name with format of A[0-9]+.";
      msgSStr << endl;
    }

  } // for each column, inspect it for an Allele name.

  if (numAlleleColumns == 0) {
    okHeaderLine = false;
    msgSStr << ttableFileName << ": invalid TsvFile header first-line, ";
    msgSStr << "no Allele Columns found." << endl;

  }

  tsv.close();

  if (!okHeaderLine) {
    if (numAlleleColumns > 4) {
      Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str(), false);
    } else {
      Verbose::out(ADT_VERBOSE_NORMAL, ttableFileName + ": invalid TsvFile is missing the header line.", false);
    }
  }

  return okHeaderLine;

}
// end convertTCDArrayToVector
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::getTranslationTableFileType (STATIC)
 * Synopsis:
 *
 * DMET2 or DMET3 indicator method.
 * 
 * The translation table can be either DMET2 or DMET3 format, return
 * the format. Determine the format by inspection. A single DMET common
 * column "Gene" is used to validate the file is a translation table.
 * A specific DMET3 column is used to deferentiate between DMET2 and DMET3.
 *
 * @param ttableFileName - the file to inspect.
 *
 * @return - DMET2, DMET3 or invalid.
 */
/*****************************************************************************/
ADT_TRANSLATION_TABLE_TYPE_ENUM TranslationTableModel::getTranslationTableFileType(const std::string & ttableFileName)
{

  if (ttableFileName.empty()) {
    return ADT_TRANSLATION_TABLE_TYPE_INVALID;
  }

  TsvFile tsv;

  //
  tsv.m_optAutoTrim   = true; // remove '"'s
  tsv.m_optQuoteChar1 = 0;    // ignore "'"s

  if (tsv.open(ttableFileName) != TSV_OK) {
    APT_ERR_ABORT(ttableFileName + ": failed opening input Tsv file.");
  }

  const std::string commonName = "Gene";
  bool commonNameFound = false;

  for (int i = 0; i < tsv.getColumnCount(0) ; i++) {

    std::string columnName;

    if (tsv.cidx2cname(0, i, columnName) == TSV_OK) {

      if (columnName == commonName) {
        commonNameFound = true;
      } else if (commonNameFound && (pcrecpp::RE("Switch.Design.Strand").PartialMatch(columnName))) {
        return ADT_TRANSLATION_TABLE_TYPE_DMET3;
      }

    }

  }


  if (commonNameFound) {

    return ADT_TRANSLATION_TABLE_TYPE_DMET2;
  }

  return ADT_TRANSLATION_TABLE_TYPE_INVALID;

}
// end TranslationTableModel::getTranslationTableFileType
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::TranslationTableModel
 * Synopsis:
 *
 * Main constructor that slurps in the input allele
 * translation table file into main memory a table.
 *
 * @param rte - the run single instance time envrionemnt
 * @param ttableFileName - the translation table to slurp into memory
 * @param ttableFileType - DMET2 or DMET3
 *
 */
/*****************************************************************************/
/*****************************************************************************/
TranslationTableModel::TranslationTableModel(const RunTimeEnvironment &rte,
    const std::string & ttableFileName,
    ADT_TRANSLATION_TABLE_TYPE_ENUM ttableFileType)
    :
    TranslationInputTsvTableModel(rte,
                                  ttableFileName,
                                  (ttableFileType == ADT_TRANSLATION_TABLE_TYPE_DMET2) ?  TTM_DMET2_COLUMN_DEFINTIONS : TTM_DMET2_COLUMN_DEFINTIONS,
                                  (ttableFileType == ADT_TRANSLATION_TABLE_TYPE_DMET2) ? (size_t) sizeof(TTM_DMET2_COLUMN_DEFINTIONS) / sizeof(TTM_DMET2_COLUMN_DEFINTIONS[0]) : (size_t) sizeof(TTM_DMET3_COLUMN_DEFINTIONS) / sizeof(TTM_DMET3_COLUMN_DEFINTIONS[0]),
                                  false,
                                  &convertTCDArrayToVector), m_ttableFileType(ttableFileType)
{

  std::string prevGene;
  std::string gene;

  // DMET2 translation table converstion from DMET3 column

  int numDmet2Columns = (int) sizeof(DMET3_TO_DMET2_COLUMN_MAP) / sizeof(DMET3_TO_DMET2_COLUMN_MAP[0]);

  for (int i = 0; i < numDmet2Columns; i++) {
    m_dmet3ToDmet2ColumnIndex[DMET3_TO_DMET2_COLUMN_MAP[i].m_dmet3Index]
    = DMET3_TO_DMET2_COLUMN_MAP[i].m_dmet2Index;
  }


  for (int i = 0, headerRow = 0; i < m_rows.size(); i++, prevGene = gene) {

    gene = m_rows[i][getColumnIndex(ADT_DMET3_TT_GENE)];

    if (prevGene != gene) {
      m_geneCopyNumberIndicator[gene] = false;
      headerRow = i;
    } else {
      // Ok, for tri-allelic markers an probeSet will repeat and not be Unique.
      // However, an probeSet can never be in both a Hapoltype group and
      // non-Haplotype so repeats are guaranteed to always be the same value.
      // Setting multiple times is no harm, no foul.
      std::string probeSet = m_rows[i][getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];
      m_probeSetToIsHaplotypeMarker[probeSet]
      = (m_rows[i][getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)] == "Y");

      m_probeSetRowIndex[probeSet].push_back(i);
      // Used in probe set filtering.
      m_probeSetHeaderRowIndex[probeSet] = headerRow;
    }
    m_rowHeaderRowIndex[i]             = headerRow;

  }

}
// end TranslationTableModel::TranslationTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::probeSetFilter:
 * Synopsis:
 * Filter the translation table with a set of probe set ids passed in
 * using the marker list option (-m, --marker-list).
 *
 * Suppress ProbeSet Ids not in the probe set list within the translation
 * table by ignoring rows in the  "m_rows" data structure. T
 * that would be it pretty much, easy peasy,
 * except we need to collapse overlapping allele names.
 *
 * @param rte - the single instance run time environment
 * @param mlm - The single instance marker list table model (1 column) of probes to use.
 *
 * @return true - if it's ok to continue runnning.
 */
/*****************************************************************************/
bool TranslationTableModel::probeSetFilter(const RunTimeEnvironment & rte, MarkerListModel & mlm)
{

  bool okToContinue = true;


  // deletion in the m_rows is tricky. Deleting the row will
  // move the relative positions.

  std::list< std::string > deleteProbeSetList;
  std::list< int > deleteRows;
  std::list< int > deletedHaplotypeMarkers;

  std::map< std::string, std::vector< int >  >::const_iterator itSVS;

  // Put all probe sets into the delete list so they can then be deleted.  
  for (itSVS = m_probeSetRowIndex.begin(); itSVS != m_probeSetRowIndex.end(); itSVS++) {
    deleteProbeSetList.push_back(itSVS->first);
  }

  // Delete the indicated markers in the delete list.
  for (int i = 0; okToContinue && (i < mlm.m_probeSetList.size()); i++) {

    if (m_probeSetRowIndex.count(mlm.m_probeSetList[i])) {

      deleteProbeSetList.remove(mlm.m_probeSetList[i]);
    }
  }

  if (! okToContinue) {
    return okToContinue;
  }

  std::list< std::string >::const_iterator itLS;

  for (itLS = deleteProbeSetList.begin();
       (itLS != deleteProbeSetList.end()) && okToContinue; itLS++) {

    std::string probeSet = *itLS;

    // Multi-allelic probe sets have multiple rows. 
    for (int j = 0; j < m_probeSetRowIndex[probeSet].size(); j++) {

      // Do NOT delete the row erstwhile we have indexes that would need
      // to be rebuilt. Instead, just mark the row as one to be ignored.
      ignoreRow(m_probeSetRowIndex[probeSet][j], true);

      if (isHaplotypeMarker(probeSet)) {
        deletedHaplotypeMarkers.push_back(m_probeSetRowIndex[probeSet][j]);
      }
    }
    m_probeSetRowIndex.erase(probeSet);
  }

  if (okToContinue && !deletedHaplotypeMarkers.empty()) {
    okToContinue = _reconcileDeletedHaplotypeRows(rte, deletedHaplotypeMarkers);
  }

  if (m_probeSetRowIndex.size() == 0) {
    APT_ERR_ABORT("The translation file has no markers in the supplied marker list.");
  }

  Verbose::out(ADT_VERBOSE_NORMAL, "The translation file has " + ToStr(m_probeSetRowIndex.size()) + " markers in the supplied marker list.");

  // Fix up any indexes that were based upon table row position.

  return okToContinue;
}
// end TranslationTableModel::probeSetFilter
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::getColumnIndex:
 * Synopsis:
 * DMET2 & DMET3 columns can be retrieved via this API without knowing
 * the table type if the column is shared. For columns distinct to the
 * type, use the ENUM provided for the table. 
 *

 * @param dmet3Column - enum in the header file.
 *
 * @return - the index
 */
/*****************************************************************************/
int TranslationTableModel::getColumnIndex(ADT_DMET3_TT_COLUMNS_ENUM dmet3Column) const
{

  if (getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3) {
    return dmet3Column;
  }


  if (m_dmet3ToDmet2ColumnIndex.count(dmet3Column) == 0) {
    APT_ERR_ABORT(ToStr(dmet3Column) + ":invalid column index passed to getColumnIndex.");
  }
  APT_ERR_ASSERT(m_dmet3ToDmet2ColumnIndex.count(dmet3Column) > 0, "");

  std::map< ADT_DMET3_TT_COLUMNS_ENUM, ADT_DMET2_TT_COLUMNS_ENUM >::const_iterator it;

  it = m_dmet3ToDmet2ColumnIndex.find(dmet3Column);

  APT_ERR_ASSERT(it != m_dmet3ToDmet2ColumnIndex.end(), "");

  return it->second;

}
// end getColumnIndex
/*****************************************************************************/
/*****************************************************************************/
/**
 * getColumnRegex:
 * Synopsis:
 *
 * Returns the validation regular expression std::string for the column.
 *
 * @param index - column index
 *
 * @return - the regex as pcrecpp::RE
 */
/*****************************************************************************/
pcrecpp::RE TranslationTableModel::getColumnRegex(int index)
{

  if (index < m_columnDefinition.size()) {

    return m_columnDefinition[index].m_validRE;
  }

  // convert int to std::string
  std::stringstream ss;
  ss << index;
  std::string indexAsString = ss.str();

  APT_ERR_ABORT(indexAsString + std::string(": programing Error: bad column index passed to TranslationTableModel::getColumnIndex"));

  return pcrecpp::RE("");

}
// end TranslationTableModel::getColumnRegex
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::getCopyNumberColumn:
 * Synopsis:
 *
 * The allele columns A1 - AN do not normally parsed for the non-haploypte
 * markers. However, the copy number column is the
 * exception. This API returns the ANN column of the copy number designation
 * or -1 if no copy number column is found.
 *
 *
 * @param row - A data row in the translation table model
 *
 * @return - the column index if found or -1 else.
 */
/*****************************************************************************/
int TranslationTableModel::getCopyNumberColumn(int row)
{

  if (m_rowCopyNumberColumn.count(row) > 0) {
    return m_rowCopyNumberColumn[row];
  }

  m_rowCopyNumberColumn[row] = -1;

  for (int i = getColumnIndex(ADT_DMET3_TT_ALLELE_START); i < columnCount(); i++) {

    if (m_rows[row][i] == "0") {
      m_rowCopyNumberColumn[row] = i;
      break;
    }
  }

  return  m_rowCopyNumberColumn[row] ;

}
// end TranslationTableModel::getCopyNumberColumn
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::getHaplotypeAlleleColumns
 * Synopsis:
 * 
 * If a  Translation Table row is Haplotype then there will be
 * one or more allele name columns (A1 - A29) with an entry.
 *
 * If this method returns a size() of 0 then this means the row
 * was non-descriptive, technically an illegal condition for a Haplotype.
 * The exception is  not thrown  here but
 * rather the responsiblity of the callee.
 *
 * @param row - integer of the row, 0 based
 * @param sartColumn - typically "Reference".
 * 
 * @return alleleCols- the std::vector<int> of columns that have stuff.
 */
/*****************************************************************************/
std::vector<int> TranslationTableModel::getHaplotypeAlleleColumns(const int row ,
    const int startColumn)
{

  std::vector<int> alleleCols;

  int startIndex = getColumnIndex(ADT_DMET3_TT_REFERENCE);

  if (startColumn >= 0) {
    startIndex = startColumn;
  }

  pcrecpp::RE re("^\\s*#");
  // Skip commented out column names.

  for (int j =  startIndex;j < columnCount(); j++) {

    std::string alleleName = m_rows[row][j];
    if (! alleleName.empty() && !re.PartialMatch(alleleName)) {
      alleleCols.push_back(j);
    }
  }

  return alleleCols;

}
// end TranslationTableModel::getHaplotypeAlleleColumns
/*****************************************************************************/
/**
 * TranslationTableModel::getHeaderRowgetHeaderRow:
 * Synopsis:
 *
 * Given some arbitrary row in the table this routine looks
 *  back (up) for the first lower index row header corresponsing to the
 *  later child row. An invalid row will case an assert error.
 *
 * @param childRow - the child row in question.
 *
 * @return row - an integer in the set [0-childRow) corresponding to the gene heeader for that marker. 
 */
/*****************************************************************************/
int TranslationTableModel::getHeaderRow(const int childRow)
{

  std::map< int, int>::const_iterator itSI
  = m_rowHeaderRowIndex.find(childRow);

  if (itSI == m_rowHeaderRowIndex.end()) {
    cerr << "getHeaderRow: programming error with invalid childRow: ";
    cerr << childRow << endl;
    APT_ERR_ASSERT(itSI != m_rowHeaderRowIndex.end(), "");
  }

  return itSI->second;

}
// end TranslationTableModel::getHeaderRow
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::getHeaderRowColumnName
 * Synopsis:
 * A header row repeats the descriptive headers for common columns and
 * has custom headings for allele columns. Sometimes we can't just use the
 * headings provided in the Translation Table for business purposes. For
 * those cases we std::map the column to the business requirement.
 *
 *
 * @param headerRow - integer of the row, 0 based which contains a header
 * @param column   - integer of the column
 *
 * @returns columnName - may be different than the data for business purposes.
 */
/*****************************************************************************/
std::string TranslationTableModel::getHeaderRowColumnName(int headerRow,
    int column)
{


  validateRowIsHeader(headerRow);

  std::string columnName;

  if (column == getColumnIndex(ADT_DMET3_TT_REFERENCE)) {
    columnName = "Ref";
  } else if (column == getColumnIndex(ADT_DMET3_TT_VARIANT)) {
    columnName = "Var";
  } else {
    columnName = m_rows[headerRow][column];
  }

  return columnName;

}
// end TranslationTableModel::getHeaderRowColumnName
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::getProbeSetRowIndex
 * Synopsis:
 *
 *  Gets a probe set (assay id) row index. For multi-allelic probes with
 *  multiple rows in the translation table then a relative row (0,1,..count)
 *  needs to be supplied.
 *
 * @param probeSet - the assay id
 * @param multiAllelicRow - the relative row based upon the count, typically
 *                          just 0 or 1.
 *
 * @return - the row in the translation table, -1 of the probeSet is not found.
 */
/*****************************************************************************/
int TranslationTableModel::getProbeSetRowIndex(const std::string & probeSet,
    int multiAllelicRow) const
{

  std::map<std::string, std::vector<int> >::const_iterator itSI;

  itSI = m_probeSetRowIndex.find(probeSet);

  if (itSI == m_probeSetRowIndex.end()) {
    return -1;
  }

  APT_ERR_ASSERT(multiAllelicRow < itSI->second.size(), "");

  return (itSI->second[multiAllelicRow]);

}
// end TranslationTableModel::getProbeSetRowIndex
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::getProbeSetRowIndexSize
 * Synopsis:
 *
 * Counts the number rows a probe set (assay id) appears in the translation
 * table. Multi-allelic will appear more than once. 
 *
 * @param probeSet - the assay id
 * @param multiAllelicRow - the relative row based upon the count, typically
 *                          just 0 or 1.
 *
 * @return - the size (number of rows, i.e. 1 for bi-allelic), -1 if the
 *           probeSet is not found.
 */
/*****************************************************************************/
int TranslationTableModel::getProbeSetRowIndexSize(const  std::string & probeSet) const
{

  std::map<std::string, std::vector<int> >::const_iterator itSI;

  itSI = m_probeSetRowIndex.find(probeSet);

  if (itSI == m_probeSetRowIndex.end()) {
    return -1;
  }

  return (itSI->second.size());

}
// end TranslationTableModel::getProbeSetRowIndexSize
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::getReportAllele
 * Synopsis:
 * 
 * A lookup into the previously set "m_reportAlleles" index. 
 *
 * @param probeSet - half the key
 * @param allele -  the other haf of the key 
 * 
 *
 * @return reportAllele - the report allele if it exists. 
 */
/*****************************************************************************/
std::string TranslationTableModel::getReportAllele( const std::string probeSet, const std::string allele ) {

  std::string key = probeSet + ":" + allele;

  if ( m_reportAlleles.count( key ) > 0 ) {
    return m_reportAlleles[key];
  }

  return allele;

}
// end TranslationTableModel::getReportAllele
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::ignoreRow
 * Synopsis:
 * A setter and getter for ignoring rows in the translation table
 * via a commented out probe set id column.
 * 
 * Business logic surrounding when a translation table row is to be
 * ignored. Currently this is only when the ProbeSet column begins with
 * '#'. If this every changes add the business logic here.
 *
 * @param row - the row in question, 0 based.
 * @param bool - assignIgnore, set rather than get which results in '#' being prepended to the probe set id.
 *
 * @return true - indicates to ignore this row.
 */
/*****************************************************************************/
bool TranslationTableModel::ignoreRow(const int row, bool assignIgnore)
{


  bool okRow = true;

  APT_ERR_ASSERT((row >= 0) && (row < m_rows.size()), "");

  if (assignIgnore) {
    m_rows[row][getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)] = "#" + m_rows[row][getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];
  } else {

    okRow = RE_IGNORE_ROW.PartialMatch(m_rows[row][getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)]);

  }

  return okRow;

}
// end TranslationTableModel::ignoreRow
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::isHaplotypeMarker
 * Synopsis:
 * Is a row a Haplotype row? Yes if the 'Haplotype' column contains a 'Y'.
 *  
 * Vagaries of the STL std::map. If an access to a std::map key doesn't exist it
 * gets created. We don't want that. Therefore we break out this
 * function to do the iterator look up.
 * If an unknown probeSet is passed then false is returned.
 *
 *
 * @param probeSet - probeSet within the translation table.
 * 
 * @return true - the probeSet exists and Haplotype column is marked 'Y'.
 */
/*****************************************************************************/
bool TranslationTableModel::isHaplotypeMarker(const std::string & probeSet) const
{

  if (m_probeSetToIsHaplotypeMarker.count(probeSet) == 0) {
    return false;
  }

  std::map< std::string, bool >::const_iterator it =  m_probeSetToIsHaplotypeMarker.find(probeSet);

  if ((it != m_probeSetToIsHaplotypeMarker.end()) && (it->second  == true)) {
    return true;
  }

  return false;

}
// end TranslationTableModel::isHaplotypeMarker
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::isRowHaplotype
 * Synopsis:
 * 
 * Used only to validate that the input file for the Translation Table
 * has a header which names the gene Allele possible translations
 * before reading any records for a gene.
 *
 *
 * @param ttm - TranslationTableModel.
 * @param row - the row in question, 0 based.
 *
 * @return true - if the row is indeed a header row.
 */
/*****************************************************************************/
bool TranslationTableModel::isRowHaplotype(const int row)
{

  APT_ERR_ASSERT((row >= 0) && (row < m_rows.size()), "");

  pcrecpp::RE re("[Y]");

  return re.FullMatch(m_rows[row][getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)]);

}
// end TranslationTableModel::isRowHaplotype
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::setReportAlleles
 * Synopsis: 
 *
 * Hash a probeset:TranslationTableAllele - ReportAllele.
 * Note that this method cannot be put in the constructer because
 * the requisite GenoCallCoder must be instiated with headers from the
 * first CHP file. Parsing the first CHP file comes far later down stream
 * then when the constructor is called. 
 *
 *
 * @param gcc - The genoCallCoder initialized with the first CHP
 *              headers. 
 */
/*****************************************************************************/
void TranslationTableModel::setReportAlleles( GenoCallCoder & gcc ) {

  int markerCount = 0;
  int reportCount = 0;
  int switchStrandCount = 0;
  int differenceCount = 0;

  std::stringstream reportSStr;
  
  for ( int row = 0; row < size() ; row++ ) {

    if (ignoreRow(row)) {
      continue;
    }

    // Ignore header rows.
    std::string switchStrand = m_rows[row][getColumnIndex(ADT_DMET3_TT_SWITCH_STRAND)];
    if ( !((switchStrand == "N") || (switchStrand == "Y")) ) {
      continue;
    }


    markerCount++;
    
    std::string probeSet    = m_rows[row][getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];
    std::string reference   = m_rows[row][getColumnIndex(ADT_DMET3_TT_REFERENCE)];
    std::string variant     = m_rows[row][getColumnIndex(ADT_DMET3_TT_VARIANT)];
    bool useComplement = switchStrand == "Y";

    std::string ttReference = reference;
    std::string ttVariant = variant;

    if ( useComplement ) {
      switchStrandCount++;
      reference = CallElement::reverseComplement(reference);
      variant  = CallElement::reverseComplement(variant);
    }
    
    std::string key = probeSet + ":" + ttReference;
    std::string reportAllele;

    
    if ( m_reportAlleles.count(key) == 0 ) {

      try {
        reportAllele = gcc.referenceAlleleToReportAllele( probeSet, reference );
        if ( !reportAllele.empty() && reportAllele != reference ) {
          reportCount++;
          m_reportAlleles[key] = reportAllele;
        }
        if ( ttReference != reportAllele ) {
          differenceCount++;
        }
        reportSStr << markerCount << ":" << reportCount << ":" << differenceCount << "{" << switchStrand << ":" << switchStrandCount << "} Ref " << key << "(" << ttReference << ") => " << reportAllele << endl;
      }
      catch (...) {
      }
    }
    key = probeSet + ":" + ttVariant;
    if ( m_reportAlleles.count(key) == 0 ) {
      try {
        reportAllele = gcc.referenceAlleleToReportAllele( probeSet, variant );
        if ( !reportAllele.empty() && reportAllele != variant ) {
          reportCount++;
          m_reportAlleles[key] = reportAllele;
        }
        if ( ttVariant != reportAllele ) {
          differenceCount++;
        }
        reportSStr << markerCount << ":" << reportCount << ":" << differenceCount << "{" << switchStrand << ":" << switchStrandCount << "} Var " << key << "(" << ttVariant << ") => " << reportAllele << endl;
          
      }
      catch (...) {
      }
    }
  }

  reportSStr << "Total Markers: " << markerCount << endl;
  reportSStr << "Total Report Makers different from translation table: " << differenceCount << endl;
  reportSStr << "Total Switch Strand = 'Y' count: " << switchStrandCount << endl;
  Verbose::out(ADT_VERBOSE_INPUT_FILES, reportSStr.str());
  
  
}
// end TranslationTableModel::setReportAlleles
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::validateRowIsHeader
 * Synopsis:
 * 
 * Used to validate that the input file for the Translation Table has a header
 * which names the gene Allele possible translations before reading any
 * records for a gene.
 *
 * An member attribute caches the result after the first assessment so
 * use with impunity. 
 *
 *
 * SAMPLE:
 * ABCB1 Reference Link ProbeSet refSNP ID Defining cDNA Nucleotide Position Genomic Nucleotide Position Change External ID Validated Haplotype Reference Variant WT *14 A893S A893T  N44S I829V S400N A599T V801M
 *
 *
 * @param row - the row in question, 0 based.
 * @param giveWarning - boolean indicator for messaging
 *
 * @return true - if the row is indeed a header row.
 */
/*****************************************************************************/
bool TranslationTableModel::validateRowIsHeader(const int row, bool giveWarning)
{

  APT_ERR_ASSERT((row >= 0) && (row < m_rows.size()), "");

  bool okToContinue = true;

  if (m_validatedHeaderRows.count(row)) {
    return okToContinue;
  }

  for (int j = 0; (j < columnCount()) && okToContinue; j++) {

    pcrecpp::RE regex = getColumnRegex(getColumnIndex(ADT_DMET3_TT_ALLELE_START));

    if (j < getColumnIndex(ADT_DMET3_TT_ALLELE_START)) {
      regex = getColumnRegex(j);
    }

    if (!(regex.FullMatch(m_rows[row][j]))) {
      if (giveWarning) {
        std::stringstream msg;

        msg << "Invalid row(" << m_fileRowIndex[row] << "), column(";
        msg << j << " ) entry of \"" << m_rows[row][j] << "." << endl;
        msg << "[ ROW " << m_fileRowIndex[row] << "] ";
        msg << getRowAsString(row) << endl;
        Verbose::warn(ADT_VERBOSE_EXCEPTION,  msg.str());
      }
      okToContinue = false;
    }
  }

  if (okToContinue) {
    m_validatedHeaderRows.insert(row);
  }

  return okToContinue;

}
// end TranslationTableModel::validateRowIsHeader
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::_reconcileDeletedHaplotypeRows
 * Synopsis:
 *
 * When haplotype rows are deleted, like with the probe set filter,
 * this routine fixes up the haplotype allele names, condensing duplicate
 * calls to one call. Two calls are equivalent after a marker has been
 * deleted if two  new groups are equivalent. In that case, one allele
 * call column is designated as the combined column and any matching
 * columns are commented out (i.e. the allele name). The allele name
 * for commented out columns is then appended to the designated call
 * column.
 *
 * @param rte - the single instance run time environment. a
 * @param deletedRows - those haplotype group markers that were deleted.
 *
 * @return - true if program execution is ok to continue.
 */
/*****************************************************************************/
bool TranslationTableModel::_reconcileDeletedHaplotypeRows(const RunTimeEnvironment & rte, std::list< int > deletedRows)
{

  bool okToContinue = true;

  // Multiple rows may participate in the same gene.
  // Therefore first we isolate the distinct headers.
  set< int > reconcileHeaderRows;



  pcrecpp::RE reSkip("^\\s*#");
  std::list< int >::iterator itLI;

  for (itLI = deletedRows.begin(); itLI !=  deletedRows.end(); itLI++) {
    int headerRow = getHeaderRow(*itLI);
    if (reconcileHeaderRows.count(headerRow) == 0) {
      reconcileHeaderRows.insert(headerRow);
    }
  }

  // Now we have a list of just the genes we care about. Create CallSet objects
  // for each group. We capitalize on the operator== in the CallSet to
  // compare groups.

  set< int >::iterator itSI;
  for (itSI = reconcileHeaderRows.begin();
       itSI !=  reconcileHeaderRows.end(); itSI++) {

    int headerRow = *itSI;

    std::vector< int > alleleCols = getHaplotypeAlleleColumns(headerRow, getColumnIndex(ADT_DMET3_TT_ALLELE_START));
    std::map< int,  CallSet > haplotypeGroup;


    for (int row = headerRow + 1;
         (row < m_rows.size()) && (getHeaderRow(row) == headerRow)
         ; row++) {

      if (ignoreRow(row)) {
        continue;
      }

      for (int j = 0; j < alleleCols.size(); j++) {

        int col = alleleCols[j];

        // Wierd side effect of std::maps. Accessing an element that doesn't
        // exist creates it. In this case we want to create the
        // empty CallSet. This is a side effect of the translation
        // table as well where the "reference" allele has all empty
        // columns. We need to create the empty CallSet in that case.
        // Also, the common case will be that a group with a filtered
        // probe set will have an empty CallSet for the allele that
        // corresponded to that probe set.

        if (haplotypeGroup.count(col) == 0) {
          haplotypeGroup[col].m_name = m_rows[headerRow][col];
          haplotypeGroup[col].m_isDescriptive = true;
        }
        if (! m_rows[row][col].empty()) {
          haplotypeGroup[col].addCallElement(*this, row, col, true, true);
        }

      } // for each gene column

    } // for each gene row

    std::map< int, CallSet >::iterator itMIC;

    for (itMIC = haplotypeGroup.begin(); itMIC != haplotypeGroup.end(); itMIC++) {


      std::map< int, CallSet >::iterator jtMIC;

      for (jtMIC = itMIC, jtMIC++; jtMIC != haplotypeGroup.end(); jtMIC++) {

        // Skip this column if it has already been commented out,
        // i.e. combined with another group.
        if (reSkip.PartialMatch(m_rows[headerRow][jtMIC->first]) ||
            (!jtMIC->second.isEmpty() && ( jtMIC->second.getCallType() != ADT_CALL_TYPE_HAPLOTYPE_GROUP ) )) {
          continue;
        }

        // Compare the CallSet groups.

        if (itMIC->second == jtMIC->second) {

          // Rolling up the names was added back in as an option
          // on 10/15/2008.

          if (rte.m_adtOpts.m_useFirstDupAlleleDef) {

            //Comment out the combined allele name so it will be skipped
            // for any processing down stream.
            m_rows[headerRow][jtMIC->first] = "#" + m_rows[headerRow][jtMIC->first];
          } else {
            // Rolling up the names was canceled
            // by Elaine and Carsten on 10/01/2008
            // Code left in place in case they change their mind.
            std::stringstream alleleCombinedNameSStr;
            alleleCombinedNameSStr << m_rows[headerRow][itMIC->first] << "_or_";
            alleleCombinedNameSStr << m_rows[headerRow][jtMIC->first];
            m_rows[headerRow][itMIC->first] = alleleCombinedNameSStr.str();
            Verbose::out(ADT_VERBOSE_INPUT_FILES, "Combined allele " +  m_rows[headerRow][jtMIC->first] + " with " +  m_rows[headerRow][itMIC->first]);
          }

        }
      }
    }

  } // foreach gene row header that needs reconciled



  return okToContinue;

}

// end TranslationTableModel::_reconcileDeletedHaplotypeRows
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationTableModel::generateGeneMarkerList
 * Synopsis:
 *
 * Copy Number Zero operations require the set of all markers or probe
 * sets for a gene during validation and reporting. This API builds
 * this list. Broken out as a separate API and not part of the
 * constructor because the MarkerList filter needs to be applied
 * before this API is called. The MarkerList filter happens outside
 * the constructor and therefore so does this operation. 
 *
 *
 */
/*****************************************************************************/
void TranslationTableModel::generateGeneMarkerList()  {

  std::map< std::string, std::vector< int >  >::const_iterator itSVI;

  // The m_probeSetRowIndex is filtered with only the filtered ids. 
  for (itSVI = m_probeSetRowIndex.begin(); itSVI != m_probeSetRowIndex.end(); itSVI++) {

    m_geneMarkers[m_rows[itSVI->second[0]][getColumnIndex(ADT_DMET3_TT_GENE)]].push_back( itSVI->first ) ;

  }

  return;

}
// end TranslationTableModel::generateGeneMarkerList
/*****************************************************************************/
