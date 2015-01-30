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
 * @file   TranslationInputStreamTableModel.cpp
 * @author Mybrid Spalding
 * @date   Fri Jun 13 13:09:30 PDT 2008
 * @brief  Started as base class data model for either CHP or TSV files, but ended as an exclusive class to GenotypeTableModel and contains GenotypeTableModel data. 
 */

//
#include "translation/TranslationInputStreamTableModel.h"
//
#include "translation/CallElement.h"
#include "translation/TranslationTableModel.h"
//
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/utils/src/AffymetrixGuid.h"
#include "calvin_files/utils/src/GenoCallCoder.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "file/CHPFileData.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <iostream>
#include <sstream>


using namespace affx;
using namespace std;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;

///////////////////////////////////////////////////////////////////////////////
// BEGIN TranslationInputStreamTableModel
///////////////////////////////////////////////////////////////////////////////

/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::TranslationInputStreamTableModel
 * Synopsis:
 *
 * CHP constructor.
 *
 * Note that most CHP data is completely ignored excepting the alleles
 * and the probe set id.
 * Be only using three pieces of  CHP data from the CHP file
 * (probe set id, allele1 and allele2) allows the design to be coallece CHP
 * files into the TSV data model for the DMET2, Genotype Short Report.
 * This is signficant because research
 * at Affy after DMET3 is release will able to use the Genotype Short Report
 * input for DMET3. Creaing text files of sample data can be done in Excel
 *  and this avoids creating CHP files soley for the purpose of
 *  running this translation program. 
 * 
 * @param rte - the single instance run time environment.
 * @param tsvColumnDefinitions  - CHP data is coerced into TSV fields. 
 * @param tcdSize - size of tsvColumnDefinitions
 * @param gotm  - the genotype override table model.
 *
 */
/*****************************************************************************/
TranslationInputStreamTableModel::TranslationInputStreamTableModel(
  const RunTimeEnvironment& rte,
  const TittmColumnDefinition tsvColumnDefinitions[],
  size_t tcdSize,
  GenotypeOverrideTableModel *gotm) : m_gotm(gotm)
{

  m_type                 = ADT_EXPERIMENT_STREAM_TYPE_CHP;
  m_chpData              = NULL;
  m_genoCallCoder        = NULL;

  _initializeColumnDefinitions(rte, tsvColumnDefinitions, tcdSize);


  // Ok, we shouldn't be trapping this exception here, but
  // during the input validation. This is a programming error
  // so just assert.
  APT_ERR_ASSERT(rte.m_adtOpts.m_inputExperimentFiles.size() != 0, "");

  for (size_t i = 0; i < rte.m_adtOpts.m_inputExperimentFiles.size(); i++) {
    m_chpExperimentFileQueue.push_back(
      rte.m_adtOpts.m_inputExperimentFiles[i]);
  }


  return;

}
// end TranslationInputStreamTableModel::TranslationInputStreamTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::TranslationInputStreamTableModel
 * Synopsis:
 *
 * Tsv constructor.
 *
 * A DMET2 constructor for the Genotype Short Report input of genotype
 * data. The schema for this object is the schema for the Genotype Short
 * Report even for CHP files. 
 *
 * @param rte - the single instance run time environment. 
 * @param inputTsvFileName - the file with the Genotype Short Report data. 
 * @param tsvColumnDefinitions - the Genotype Short Report columns
 * @param tcdSize - the tsvColumnDefinitions size
 * @param strictTsvColumnOrder - when true the columns must be in order. 
 * @param tsvExperimentColumn - if column order is not strict, which column contains the experiment. 
 * @param tsvGeneColumn - if column order is not strict, which column contains the gene. 
 * @param tsvProbeSetColumn - if column order is not strict, which column contains the probe set id.. 
 * @param gotm - the single instance Genotype Override Table Model.
 *
 */
/*****************************************************************************/
TranslationInputStreamTableModel::TranslationInputStreamTableModel(
  const RunTimeEnvironment& rte,
  const std::string& inputTsvFileName,
  const TittmColumnDefinition tsvColumnDefinitions[],
  size_t tcdSize,
  bool strictTsvColumnOrder,
  int tsvExperimentColumn,
  int tsvGeneColumn,
  int tsvProbeSetColumn,
  GenotypeOverrideTableModel *gotm) : m_gotm(gotm)
{

  m_type                 = ADT_EXPERIMENT_STREAM_TYPE_TSV;
  m_strictTsvColumnOrder = strictTsvColumnOrder;
  m_tsvExperimentColumn  = tsvExperimentColumn;
  m_tsvGeneColumn        = tsvGeneColumn;
  m_tsvFileName          = inputTsvFileName;
  m_tsvProbeSetColumn    = tsvProbeSetColumn;
  m_chpData              = NULL;
  m_genoCallCoder        = NULL;

  _initializeColumnDefinitions(rte, tsvColumnDefinitions, tcdSize);


  _openTsvFile(rte);


}
// end TranslationInputStreamTableModel::TranslationInputStreamTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::columnCount
 * Synopsis:
 * 
 * Interogate an already instantiated table for the number of columns.
 * This method assumes ALL rows have the same column width and just
 * returns the width of the column defintion structure. .
 *
 */
/*****************************************************************************/
int TranslationInputStreamTableModel::columnCount() const
{

  if (m_tsvColumnDefinition.size()  > 0) {
    return m_tsvColumnDefinition.size();
  }

  return 0;

}
// end TranslationInputStreamTableModel::columnCount
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::clearData
 * Synopsis
 * 
 * It is required that memory usage should constant as files and
 * experiments are streamed. This method resets the object to
 * no data other than the meta data required to open the next stream
 * (file). 
 *
 * Calls to this API should clear all memory in between experiments.
 *
 */
/*****************************************************************************/
void  TranslationInputStreamTableModel::clearData()
{

  m_rows.clear() ;

  if (m_chpData != NULL) {
    delete m_chpData;
    m_chpData = NULL;
  }


  return;

}
// end TranslationInputStreamTableModel::clearData
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::getColumnName
 * Synopsis:
 *
 * Get the column name for a column index, i.e. integer -> std::string.
 * It's an assertion error to pass in an out of bound index. 
 *
 * @param column - the 0 based column index for the column in question.
 * 
 * @return m_columnName - the name corresponding to column. 
 */
/*****************************************************************************/
std::string TranslationInputStreamTableModel::getColumnName(int column)
{

  // This should always be true unless some constructor is created
  // that doesn't initialize this.
  APT_ERR_ASSERT(m_tsvColumnDefinition.size(), "");

  APT_ERR_ASSERT((column >= 0) && (column < columnCount()), "");

  return m_tsvColumnDefinition[column].m_columnName;

}
// end TranslationInputStreamTableModel::getColumnName
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::getRowAsString
 * Synopsis:
 * 
 * Debug convienance routine that std::stringifies and entire row. Also
 * used in error messages. All columns are "separated" by a single comma.
 * If the data within a column contains a column there is no attempt made
 * quote the columns. Empty columns will have comma sequences as ",,,". 
 * column 1 2   3 
 * data   A B,C D => "A,B,C,D"
 * 
 * @param row - integer of the row, 0 based
 * 
 * @return std::string - the comma separated std::string
 */
/*****************************************************************************/
std::string TranslationInputStreamTableModel::getRowAsString(const int row)  const
{

  APT_ERR_ASSERT((row >= 0) && (row < size()), "");

  std::stringstream rowSStr;

  for (int j = 0 ; j < columnCount(); j++) {
    rowSStr << m_rows[row][j];

    if (j < (columnCount() - 1)) {
      rowSStr << ",";
    }

  } // for each column

  return rowSStr.str();

}
// end TranslationInputStreamTableModel::getRowAsString
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::readNextExperiment:
 * Synopsis:
 *
 *  A wrapper that wraps the appropriate TSV for CHP read next experiment
 * method. The TranslationEngine calls this method.       
 *
 * 
 * 
 * @param rte - the single instance run time environment. 
 * @param ttm - the single instance translation table model. 
 * @param geneExperimentCopyNumberCall
 *            - The copy number 0 designation per gene std::map. 
 * 
 * @return - the number of records read or 0 if no more records are left.
 */
/*****************************************************************************/
int TranslationInputStreamTableModel::readNextExperiment(const RunTimeEnvironment & rte, TranslationTableModel & ttm, const std::map< std::string, std::string> & geneExperimentCopyNumberCall, std::map< std::string, int> & geneCopyNumber)
{


  if (m_type == ADT_EXPERIMENT_STREAM_TYPE_TSV) {
    return _readNextExperimentTsv(rte, ttm, geneExperimentCopyNumberCall);
  } else if (m_type == ADT_EXPERIMENT_STREAM_TYPE_CHP) {
    return _readNextExperimentCHP(rte, ttm, geneExperimentCopyNumberCall, geneCopyNumber );
  }

  return 0;

}
// end TranslationInputStreamTableModel::readNextExperiment
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::sortTsvFileByExperiment:
 * Synopsis:
 *
 * Genotype Short Report files are NOT sorted by experiment, but need
 * to be in order to be streamed. If this program detects an unsorted
 * Genotype Short Report then it aborts. However, before it aborts
 * it sorts the TSV file and creates a new copy of the file
 * with ".sorted" as an extension. In this way the user can
 * simply call the program a second time using the ".sorted" file
 * and then continue on. 
 * 
 * This is utility routine used independent of translating calls. Use this
 * routine with its own option switch as a convienance to sort
 * an input genotype short report file. Sorting is a requirement
 * for streaming. 
 *
 *
 * @rte - the single instance RunTimeEnvironment which contains the name of the file to be sorted. 
 *
 * @return - void, creates a new file of sorted records. Throws exceptions
 *           on IO errors.
 */
/*****************************************************************************/
void TranslationInputStreamTableModel::sortTsvFileByExperiment(const RunTimeEnvironment & rte)
{

  // This should never happen unless a new constructor is created
  // that doesn't initialize columns.
  APT_ERR_ASSERT(m_tsvColumnDefinition.size(), "");

  TsvFile sortTsv;

  std::string sortName = m_tsvFileName + ".sort.txt";

  Verbose::out(ADT_VERBOSE_NORMAL, "sortTsvFileByExperiment: " + m_tsvFileName
               +  " -> " + sortName);

  for (size_t i = 0; i < m_tsvColumnDefinition.size(); i++) {

    sortTsv.defineColumn(0, i, m_tsvColumnDefinition[i].m_columnName);
  }

  //
  m_tsv.m_optAutoTrim   = true; // remove '"'s
  m_tsv.m_optQuoteChar1 = 0;    // ignore "'"s


  if (m_tsv.open(m_tsvFileName) != TSV_OK) {
    APT_ERR_ABORT(m_tsvFileName + ": failed opening input Tsv file.");
  }

  m_tsv.headersBegin();

  std::string key;
  std::string value;
  sortTsv.addHeaderComment(" File generated by sortTsvFileByExperiment.");
  while (m_tsv.headersNext(key, value) == TSV_OK) {
    sortTsv.addHeader(key, value);
  }

  sortTsv.writeTsv_v1(sortName);

  int row = 0;
  while (m_tsv.nextLevel(0) == TSV_OK) {

    m_rows.resize(row + 1);
    m_rows[row].resize(m_tsvColumnCount, "");

    // Read in a row, one column at a time.
    for (int tsvColumn = 0; tsvColumn < m_tsvColumnCount; tsvColumn++) {
      std::string columnValue;
      if (m_tsv.get(0, tsvColumn, columnValue) != TSV_OK) {
        APT_ERR_ABORT("tsv.get != TSV_OK");
      }
      if (!columnValue.empty()) {
        m_rows[row][m_tsvIndexToColumnDefinitionIndex[tsvColumn]] = columnValue;
      }
    }

    m_experimentRowIndex[m_rows[row][m_tsvExperimentColumn]].push_back(row);
    row++;
  }

  std::map<std::string, std::vector< int > >::iterator itSV;

  row = 0;
  for (itSV = m_experimentRowIndex.begin(); itSV != m_experimentRowIndex.end(); itSV++) {

    for (size_t i = 0; i < itSV->second.size() ; i++) {
      for (int tsvColumn = 0; tsvColumn < m_tsvColumnCount; tsvColumn++) {
        sortTsv.set(0, tsvColumn,  m_rows[itSV->second[i]][tsvColumn]);
      }
      sortTsv.writeLevel(0);
      row++;
    }
  }

  Verbose::out(ADT_VERBOSE_NORMAL, ToStr(row) + " sorted Genotype records written to " + sortName);

  sortTsv.close();

  m_tsv.close();

  clearData();

  return;

}
// end TranslationInputStreamTableModel::sortTsvFileByExperiment
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_closeTsvFile:
 * Synopsis:
 *
 * Helper function to TranslationInputStreamTableModel::_close().
 *
 * Call tsv.close() on the m_tsv. Made a wrapper in case messaging
 * or further clean up work is ever  needed. 
 *
 *
 * @return - void
 */
/*****************************************************************************/
void TranslationInputStreamTableModel::_closeTsvFile()
{

  m_tsv.close();

}
// end TranslationInputStreamTableModel::_closeTsvFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_getAllelesFromProbeSet:
 * Synopsis:
 *
 * Helper function to TranslationInputStreamTableModel::_readNextExperimentCHP
 *
 * Translates a probe set id's textual 'allele1' and 'allele2'
 * from the CHP file using the APIs found in the GenoCallCoder. 
 *
 *  The CHP calls are encoded as "A/A" or "A/B" or "B/B", etc. These
 *  need to be translated into the actual alleles in the annotation file.
 *  Do the look up in the GenoCallCoder and return the bases.
 * 
 * Note: the report alleles are not retrieved here but instead are
 * imported at start up time using all the probe sets and alleles in the
 * translation table.
 *
 * Also note: the CHP will contain many more probe sets then are defined in
 * the translation table. When this API returns false then this is
 * the deciding factor when to ignore a probe set marker in the
 * CHP file and the data is dropped on the floor. 
 *
 * @param rte - the single instance RunTimeEnvironment.
 * @param ttm - the single instance translation table.
 * @param dbad - the DmetBiAllelicData with probe set and call data.
 *               Note: BiAllelic is a misnomer, multi allelic is represented as
 *               well.
 * @param genoCallDecoder - the object that does translation from "A/A" -> "TGT/TTT".
 * @param allele1 - returned base for the ALLELE1 column.
 * @param allele2 - returned base for the ALLELE2 column.
 *
 * @return - true if the probeSet exists  in the translation table,
 *           the allele1 and allele2 std::strings
 */
/*****************************************************************************/
bool TranslationInputStreamTableModel::_getAllelesFromProbeSet(const class RunTimeEnvironment & rte, TranslationTableModel & ttm,  affymetrix_calvin_data::DmetBiAllelicData & dbad,  std::string & allele1, std::string & allele2)
{



  APT_ERR_ASSERT(!dbad.name.empty(), "");

  std::string probeSet    = dbad.name;

  if (ttm.m_probeSetRowIndex.count(probeSet) == 0) {
    return false;
  }


  std::string decodedAbstractCall;
  std::string decodedCall;
  std::string reportCall;

  allele1.clear();
  allele2.clear();

  bool probeSetSwitchStrand = ((ttm.getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3) && ttm.m_rows[ttm.m_probeSetRowIndex[probeSet][0]][ttm.getColumnIndex(ADT_DMET3_TT_SWITCH_STRAND)] == "Y");


  decodedCall = m_genoCallCoder->genotypeCallNumToReferenceAllele(probeSet, dbad.call);

  if (decodedCall.empty()) {
    APT_ERR_ABORT("Empty reference call returned by abstractAlleleToReferenceAllele.");
  }

  std::string first, second, copyNumber;


  if (!pcrecpp::RE("^([\\w\\d\\-]+)$").FullMatch(decodedCall, &first)) {

    if (! pcrecpp::RE("^(\\w+|\\-)/(\\w+|\\-)$").FullMatch(decodedCall, &first, &second)) {
      APT_ERR_ABORT(probeSet + ": unknown call? " + decodedCall);
    }
  }

  // Copy Number 0 translation

  if ((first == "17") || (first == "ZeroCopyNumber")) {
    first = "0";
  }

  allele1 = first;

  if (first == "NotAvailable") {
    // Returning false causes the marker to be ignored, exactly the behavior
    // required.
    //    return false;
  }

  if (!second.empty()) {
    allele2 = second;
  }


  if (probeSetSwitchStrand) {

    allele1 = CallElement::reverseComplement(allele1);
    allele2 = CallElement::reverseComplement(allele2);

  }

  return true;

}
// end TranslationInputStreamTableModel::_getAllelesFromProbeSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_getGeneFromProbeSet:
 * Synopsis:
 * Helper function to TranslationInputStreamTableModel::_readNextExperimentCHP
 *
 * The CHP file doesn't contain the gene. However, the assay id (probe set)
 * is distinct per translation table so it can be looked up in the
 * translation table.
 *
 * @param probeSet - the marker for which to get the gene
 * @param ttm - the single instance translation table containing the gene
 *
 * @return - gene as std::string, empty if the probeSet doesn't exist
 */
/*****************************************************************************/
std::string TranslationInputStreamTableModel::_getGeneFromProbeSet(const std::string probeSet,   const TranslationTableModel & ttm)
{

  std::string gene;

  if (ttm.getProbeSetRowIndex(probeSet) < 0) {
    return gene;
  }

  gene = ttm.m_rows[ttm.getProbeSetRowIndex(probeSet)][ttm.getColumnIndex(ADT_DMET3_TT_GENE)];

  return gene;
}
// end TranslationInputStreamTableModel::_getGeneFromProbeSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * searchParams
 * Synopsis:
 *
 * Helper function to TranslationInputStreamTableModel::_getCHPGuid and
 * TranslationInputStreamTableModel::_initializeNextExperimentCHPFile.
 *
 * The calvin_files directory contains APIs to search for headers  or
 * parameters in the CHP file's various inheritance tree. However, these
 * APIs are specific the the class implemented for. This search
 * searchs the entire inheritance graph in one convienant API.
 * This function also convienantly handles all the conversion from
 * WCS to std::string. 
 *
 * Note: This is not the most efficient API but the two functions
 * being helped are
 * initialization functions only used at start up and hence performance
 * is not a factor. 
 *
 * @param gdh - GenericDataHeader from m_chpData
 * @params searchName - what to search for.
 * @params searchValue - the return value. 
 *
 * @return - searchValue, if found
 */
/*****************************************************************************/
static void searchParams(affymetrix_calvin_io::GenericDataHeader & gdh, const std::string searchName ,  std::string &searchValue )
{

  if (!searchValue.empty()) {
    return;
  }

  ParameterNameValueTypeIt itPNVT;
  ParameterNameValueTypeIt jtPNVT;

  gdh.GetNameValIterators(itPNVT, jtPNVT);

  for (; itPNVT != jtPNVT; itPNVT++) {

    if (StringUtils::ConvertWCSToMBS(itPNVT->GetName()) == searchName ) {
      searchValue = StringUtils::ConvertWCSToMBS(itPNVT->ToString());
      return;
    }
  }

  std::vector<GenericDataHeader>::iterator itGDH;
  std::vector<GenericDataHeader>::iterator jtGDH;

  gdh.GetParentIterators(itGDH, jtGDH);

  for (; itGDH != jtGDH ; itGDH++) {
    searchParams(*itGDH, searchName, searchValue);
  }

  return;

}
// end searchParams
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_getCHPGuid:
 * Synopsis:
 *
 * Helper function to _initializeNextExperimentCHPFile
 * 
 *   Returns the CHP file assay id. Broken out due to confusion about
 * which GUID is the right guid. Also, some prototype CHP files stuff
 * the GUID in a different place than that required by the 
 * canonical gd.ArrayIdentifier(). If gd.ArrayIdentifier fails then
 * searchParams is called. Only prototype CHPS files have exhibited this
 * behavior. 
 *
 * 
 * 
 * @param rte - the single instance run time environment. 
 * @param gd - GenericData from m_chpData
 * 
 * @return - guid, The GUID std::string.
 */
/*****************************************************************************/
/*****************************************************************************/
std::string TranslationInputStreamTableModel::_getCHPGuid(const RunTimeEnvironment & rte, affymetrix_calvin_io::GenericData & gd)
{


  //m_chpGuid = m_chpData->GetFileHeader()->GetGenericDataHdr()->GetFileId();
  //m_chpGuid = gd.ArrayIdentifier();

  AffymetrixGuidType guid;

  guid = gd.ArrayIdentifier();


  if (guid.empty()) {
    GenericDataHeader* gdh = gd.Header().GetGenericDataHdr();

    if (gdh  == NULL)  {
      APT_ERR_ABORT("Can't find header: " + ToStr(ARRAY_TYPE_IDENTIFIER));
    }
    searchParams(*gdh, StringUtils::ConvertWCSToMBS(ARRAY_ID_PARAM_NAME) ,guid);
    if (guid.empty()) {
      APT_ERR_ABORT(m_chpExperimentFile + ": invalid CHP file, missing Array ID GUID (" + ToStr(ARRAY_TYPE_IDENTIFIER) + ").");
    }
  }

  return guid;

}
// end TranslationInputStreamTableModel::_getCHPGuid
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_getExternalIdFromProbeSet:
 * Synopsis:
 * Helper function to TranslationInputStreamTableModel::_readNextExperimentCHP
 *
 * The CHP file doesn't contain the external id. However, the probe set id
 * is distinct per translation table so it can be looked up in the
 * translation table.
 *
 * @param probeSet - the marker
 * @param ttm      - the single instance translation table
 *
 * @return - external id (common name)  as std::string, empty if the probeSet doesn't exist
 */
/*****************************************************************************/
std::string TranslationInputStreamTableModel::_getExternalIdFromProbeSet(const std::string probeSet,   const TranslationTableModel & ttm)
{

  std::string externalId;

  if (ttm.getProbeSetRowIndex(probeSet) < 0) {
    return externalId;
  }

  externalId = ttm.m_rows[ttm.getProbeSetRowIndex(probeSet)][ttm.getColumnIndex(ADT_DMET3_TT_COMMON_NAME)];

  return externalId;
}
// end TranslationInputStreamTableModel::_getExternalIdFromProbeSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_initializeColumnDefinitions
 * Synopsis: 
 *
 * Helper function to the constructers
 * TranslationInputStreamTableModel::TranslationInputStreamTableModel
 * 
 * Initialize the definition of columns data structure. This is compared
 * against the input TSV file when it is opened for integrity. 
 * 
 * @param rte - the single instance run time environment. 
 * @param columnDefinitions - the typically hard coded columns. 
 * @param tcdSize - the size of columnDefinitions. 
 *
 */
/*****************************************************************************/
void TranslationInputStreamTableModel::_initializeColumnDefinitions(const RunTimeEnvironment & rte, const TittmColumnDefinition columnDefinitions[], int tcdSize)
{

  m_tsvColumnDefinition.clear();

  // Iterate so we can doublecheck the m_index value.
  for (int i = 0; i < tcdSize; i++) {

    APT_ERR_ASSERT(columnDefinitions[i].m_index == i, "");

    m_tsvColumnDefinition.push_back(columnDefinitions[i]);
  }

  return;
}
// end TranslationInputStreamTableModel::_initializeColumnDefinitions
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_initializeNextExperimentCHPFile:
 * Synopsis:
 * Helper function to TranslationInputStreamTableModel::_readNextExperimentCHP
 * 
 * Open the next CHP file, read the headers and postion things to read the
 * first data record.
 *
 * Initialize the global report alleles.
 * The report alleles global to all experiments are 
 * initialized here. There are two reasons for this:
 * 1.) The GenoCallCoder requires parameters taken from a CHP file
 * specific to this set of CHP files.
 * 2.) The first CHP file is read for the first time here. The parameters
 * are take from that file here. 
 * 
 *
 * @param rte - the RunTimeEnvironment with the command line list
 *              of all the experiment CHP files.
 * @param ttm - the single instance translation table model updated with the
 * report alleles on the first call. 
 *
 * @return - true if ok, false otherwise. Throws Exception on IO error. Updates the ttm with the report alleles on first call. 
 */
/*****************************************************************************/
bool TranslationInputStreamTableModel::_initializeNextExperimentCHPFile(const RunTimeEnvironment & rte, TranslationTableModel & ttm)
{

  bool okToContinue = true;

  if (m_chpExperimentFileQueue.size() == 0) {
    return false;
  }

  m_chpExperimentFile = m_chpExperimentFileQueue.front();
  m_chpExperimentFileQueue.pop_front();


  // Open the file and read the headers.
  // I traced this "Read" call back to GenericFileReader.cpp
  // and "Read" only reads the headers.


  m_chpData = new affymetrix_calvin_io::CHPMultiDataData;

  affymetrix_calvin_io::CHPMultiDataFileReader reader;
  reader.SetFilename(m_chpExperimentFile);

  reader.Read(*m_chpData);

  FileHeader *chpFileHeader = m_chpData->GetFileHeader();

  DataGroupHdrIt itDG, jtDG;

  chpFileHeader->GetDataGroupIts(itDG, jtDG);

  for (; itDG != jtDG; itDG++) {

    DataSetHdrIt itDSH, jtDSH;

    itDG->GetDataSetIterators(itDSH, jtDSH);

  }

  ParameterNameValueTypeList summaryList = m_chpData->GetSummaryParams();


  ParameterNameValueTypeList::iterator  itPNV;

  for (itPNV = summaryList.begin(); itPNV != summaryList.end(); itPNV++) {
    std::string name = StringUtils::ConvertWCSToMBS(itPNV->GetName());
  }

  ParameterNameValueTypeList algorithmList = m_chpData->GetAlgParams();

  GenericData & gd = m_chpData->GetGenericData();

  WStringVector names;
  gd.DataGroupNames(names);

  for (int i = 0; i < 3; i++) {

    std::string dataTypeName = StringUtils::ConvertWCSToMBS(MultiDataDataSetNames[ADT_CHP_DATA_TYPES[i]]);

    m_chpDataTypeSize[ADT_CHP_DATA_TYPES[i]] = m_chpData->GetEntryCount(ADT_CHP_DATA_TYPES[i]);

    if (m_chpDataTypeSize[ADT_CHP_DATA_TYPES[i]] == 0) {
      Verbose::warn(ADT_VERBOSE_EXCEPTION, dataTypeName  + ": group has no records (0 size)\n");
    }

  }

  // EXPERIMENT NAME, set to the file name, minus directories
  m_experimentName = Fs::basename(m_chpExperimentFile);
  bool okMatch = pcrecpp::RE("\\.[cC][hH][Pp]$").PartialMatch(m_experimentName);

  APT_ERR_ASSERT(okMatch, "");

  if (m_experimentName.empty()) {
    APT_ERR_ABORT(m_chpExperimentFile + ": invalid experiment file name.");
  }

  m_chpGuid = _getCHPGuid(rte, gd);

  //call-encoding-version=1.0 (ASCII)
  //call-encoding-max-alleles=6 (Int32)
  //call-encoding-data-size=UCHAR (ASCII)
  std::string versionHeader = "call-encoding-version";
  std::string maxAllelesHeader = "call-encoding-max-alleles";
  std::string dataSizeHeader = "call-encoding-data-size";
  std::string version, maxAlleles, dataSize;

  
  GenericDataHeader* gdh = gd.Header().GetGenericDataHdr();
  
  
  searchParams(  *gdh, versionHeader, version );
  searchParams(  *gdh, maxAllelesHeader, maxAlleles );
  searchParams(  *gdh, dataSizeHeader, dataSize );

  Verbose::out(ADT_VERBOSE_INPUT_FILES, Fs::basename(m_chpExperimentFile) + " " + versionHeader +  ": {" + version + "}" );
  Verbose::out(ADT_VERBOSE_INPUT_FILES, Fs::basename(m_chpExperimentFile) + " " +  maxAllelesHeader + ": {" + maxAlleles +  "}" );
  Verbose::out(ADT_VERBOSE_INPUT_FILES, Fs::basename(m_chpExperimentFile) + " " + dataSizeHeader + ": {" + dataSize+  "}" );

  if ( m_genoCallCoder == NULL ) {
    if ( rte.m_adtOpts.m_prototypeCHPFiles ) {
      m_genoCallCoder = new GenoCallCoder(6, "UCHAR", "1.0", '/', rte.m_adtOpts.m_inputAnnotationFile);
    }
    else {
      m_genoCallCoder = new GenoCallCoder(atoi(maxAlleles.c_str()), dataSize, version, '/', rte.m_adtOpts.m_inputAnnotationFile);
    }

    if ((ttm.getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3) &&  !rte.m_adtOpts.m_ignoreReportAllele) {

      ttm.setReportAlleles( *m_genoCallCoder);
    }

    m_chpVersion = version;
    m_chpDataEncodingSize = dataSize;
    m_chpMaxAlleles = maxAlleles;
  }
  else if ( ! rte.m_adtOpts.m_prototypeCHPFiles && ((version != m_chpVersion) || (dataSize != m_chpDataEncodingSize) ||  (maxAlleles != m_chpMaxAlleles ) ) ){
      Verbose::warn(ADT_VERBOSE_NORMAL, std::string("Conflicting CHP type: ") + Fs::basename(m_chpExperimentFile) + " " + versionHeader +  ": {" + version + "},  " +  maxAllelesHeader + ": {" + maxAlleles +  "}, " + dataSizeHeader + ": {" + dataSize +  "}" );
      APT_ERR_ABORT(std::string("CHP files with different versions detected, was expecting: ") + versionHeader + ": {" + m_chpVersion + "}, " + maxAllelesHeader + ": {" + m_chpMaxAlleles + "}, " + dataSizeHeader + ": {" + m_chpDataEncodingSize + "}");

  }



  return okToContinue;


}
// end TranslationInputStreamTableModel::_initializeNextExperimentCHPFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_initializeTsvColumns
 * Synopsis:
 * Helper function to TranslationInputStreamTableModel::_openTsvFile.
 * 
 * Initialize the TsvFile object columns from the column definitions.
 * The TSV file does not have to be in column order, but this
 * method always puts them in column order internally so that
 * indexing by column index integer can take place. 
 *
 *
 * @param m_tsvFileName - The name of the TsvFile in question used for messaging.
 * @param tsv           - The TsvFile object in question corresponding m_tsvFileName. 
 *
 * @returns false - on error.
 *
 */
/*****************************************************************************/
bool TranslationInputStreamTableModel::_initializeTsvColumns(const std::string & inputTsvFileName, affx::TsvFile & tsv)
{

  bool         okTsvColumns = true;
  std::stringstream msgSStr;

  m_tsvIndexToColumnDefinitionIndex.resize(tsv.getColumnCount(0), -1);


  for (size_t column = 0; column < m_tsvColumnDefinition.size(); column++) {

    // find the tsv column index for this
    int tsvColumn
    = tsv.cname2cidx(0, m_tsvColumnDefinition[column].m_columnName);

    if (tsvColumn < 0) {
      okTsvColumns = false;
      msgSStr << inputTsvFileName << ": invalid TsvFile header first-line, ";
      msgSStr << "column [" << (column + 1) << "], \"";
      msgSStr << m_tsvColumnDefinition[column].m_columnName << "\"  not found." << endl;
      continue;
    }

    m_tsvIndexToColumnDefinitionIndex[tsvColumn]
    = m_tsvColumnDefinition[column].m_index;

  }

  if (!okTsvColumns) {
    Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str(), false);
  }

  return okTsvColumns;

}
// end TranslationInputStreamTableModel::_initializeTsvColumns
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_matchTsvColumnsWithColumnDefinitions
 * Synopsis:
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
 * @return false - the TSV columns are not complete or in column order. 
 */
/*****************************************************************************/
bool TranslationInputStreamTableModel::_matchTsvColumnsWithColumnDefinitions(const std::string & inputTsvFileName, affx::TsvFile & tsv)
{

  bool okHeaderLine = true;

  int invalidHeaderCount = 0;

  std::stringstream msgSStr;

  for (size_t i = 0; (i < tsv.getColumnCount(0)) &&
       (i < m_tsvColumnDefinition.size()) ; i++) {

    std::string columnName;

    if (tsv.cidx2cname(0, i, columnName) == TSV_OK) {

      if (columnName != m_tsvColumnDefinition[i].m_columnName) {
        okHeaderLine = false;
        msgSStr << inputTsvFileName << ": invalid TsvFile header first-line, ";
        msgSStr << "column [" << (i + 1) << "]: found \"" << columnName << "\" but was expecting \"";
        msgSStr << m_tsvColumnDefinition[i].m_columnName << "\"." << endl;
        invalidHeaderCount++;
      }
    } else {

      APT_ERR_ABORT("tsv.cidx2cname returned unexpected error");

    } // if there was a TSV_ERROR

  } // for each column, inspect it

  if (invalidHeaderCount > 0) {
    if (invalidHeaderCount < 4) {
      Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str(), false);
    } else {
      Verbose::out(ADT_VERBOSE_NORMAL, inputTsvFileName + ": invalid TsvFile is missing the header line.", false);
    }
  }

  return okHeaderLine;

}
// end TranslationInputStreamTableModel::_matchTsvColumnsWithColumnDefinitions
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_openTsvFile:
 * Synopsis:
 * Helper function to the TSV file constructor,
 * TranslationInputStreamTableModel::TranslationInputStreamTableModel
 *
 * Opens the Genotype Short Report TSV file, ready to read the first row
 * using the streaming readNextExperimentTSV API. Clears any data
 * from a previous open. 
 *
 *
 * @param rte - the single instance RunTimeEnviroment.
 *
 * @return - void, throws an exception on open error.
 */
/*****************************************************************************/
void TranslationInputStreamTableModel::_openTsvFile(const RunTimeEnvironment & rte)
{

  // This should never happen unless a new constructor is created
  // that doesn't initialize columns.

  APT_ERR_ASSERT(m_tsvColumnDefinition.size(), "");

  clearData();

  if (Verbose::getParam().m_Verbosity > ADT_VERBOSE_INPUT_FILES) {
    Verbose::out(ADT_VERBOSE_INPUT_FILES, ToStr("### TranslationInputStreamTableModel::readFile('") + m_tsvFileName + "')\n");
  }

  //
  m_tsv.m_optAutoTrim   = true; // remove '"'s
  m_tsv.m_optQuoteChar1 = 0;    // ignore "'"s

  if (m_tsv.open(m_tsvFileName) != TSV_OK) {
    APT_ERR_ABORT(m_tsvFileName + ": failed opening input Tsv file.");
  }

  bool okColumnDefinitions = true;

  if (m_strictTsvColumnOrder) {
    okColumnDefinitions =
      _matchTsvColumnsWithColumnDefinitions(m_tsvFileName, m_tsv);
  }

  if (okColumnDefinitions) {
    okColumnDefinitions = _initializeTsvColumns(m_tsvFileName, m_tsv);
  }

  m_tsvColumnCount = m_tsv.getColumnCount(0);
  m_tsvRow = 0;

  if (! okColumnDefinitions) {
    APT_ERR_ABORT(m_tsvFileName + ": invalid input file detected.");
  }



  return;

}
// TranslationInputStreamTableModel::_openTsvFile
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_readNextExperimentCHP:
 * Synopsis:
 * Helper function to TranslationInputStreamTableModel::readNextExperiment.
 *
 * Streaming CHP files, where each CHP file is one experiment.
 * Read in the entire files worth of experiment gene data.
 * Fills in the 'm_rows' table with experiment data. 
 *
 * This is the place where override data from the override file supplants each
 * record of the CHP file as it is read in. This is how translation
 * is accomplished via the override file. The override data is
 * swapped for the CHP data here. 
 *
 * 
 * @param rte - the single instance run time environment 
 * @param ttm - the single instance translation table model
 * @param geneExperimentCopyNumberCall
 *            - the DMET2 std::map per gene copy number alleles
 * @param geneCopyNumber
 *            - the DMET3 std::map of gene copy number
 * 
 * @return - The number for records (markers) read for the experiment.
 */
/*****************************************************************************/
int TranslationInputStreamTableModel::_readNextExperimentCHP(const RunTimeEnvironment & rte, TranslationTableModel & ttm, const std::map<std::string, std::string> &geneExperimentCopyNumberCall, std::map<std::string, int> & geneCopyNumber)
{

  int chpRow       = 0;
  int recordsRead  = 0;

  clearData();

  if (!_initializeNextExperimentCHPFile(rte, ttm)) {
    return 0;
  }

  std::vector< affymetrix_calvin_io::MultiDataType > dataTypes;
  dataTypes.push_back(ADT_CHP_DATA_TYPES[BI_TYPE]);
  dataTypes.push_back(ADT_CHP_DATA_TYPES[MULTI_TYPE]);

  for (size_t j = 0; j <  dataTypes.size(); j++) {

    for (int i = 0; i < m_chpDataTypeSize[dataTypes[j]] ; i++) {

      affymetrix_calvin_data::DmetBiAllelicData dbad;

      m_chpData->GetEntry(dataTypes[j], i , dbad);

      // Not all probe sets are going to exist in the translation table,
      // if not then just drop the probeSet data on the floor.
      if (ttm.getProbeSetRowIndex(dbad.name) < 0) {
        continue;
      }


      std::string allele1, allele2, originalAllele1, originalAllele2;

      originalAllele1.clear();
      originalAllele2.clear();

      if (! _getAllelesFromProbeSet(rte, ttm, dbad, allele1, allele2)) {
        continue;
      }

      // Genotype Override
      if (m_gotm != NULL) {
        std::vector< std::string > alleleOverride = m_gotm->getOverrideAlleles(m_experimentName, dbad.name);
        if (alleleOverride.size() > 0) {
          if ( CallElement::isWildCardBase(allele1)  ) {
            m_gotm->m_chpFileProbeSetOriginalBasecall[m_experimentName + "::" + dbad.name] = allele1;
          } else {
            m_gotm->m_chpFileProbeSetOriginalBasecall[m_experimentName + "::" + dbad.name] = allele1 + "/" + allele2;
          }
          
          originalAllele1 = allele1;
          originalAllele2 = allele2;
          allele1 = alleleOverride[0];
          if ( alleleOverride.size() == 2 ) {
            allele2 = alleleOverride[1];
          }
          else {
            allele2 = "";
          }
        }
      }

      std::vector<std::string> newRow;
      chpRow++;

      newRow.resize(GT_END_FIELD_INDEX, "");

      newRow[GT_SAMPLE_INDEX]      = "N/A";
      newRow[GT_EXPERIMENT_INDEX]  = m_experimentName;
      newRow[GT_GENES_INDEX]       = _getGeneFromProbeSet(dbad.name, ttm);
      newRow[GT_EXTERNAL_ID_INDEX] = _getExternalIdFromProbeSet(dbad.name, ttm);
      newRow[GT_PROBE_SET_INDEX]   = dbad.name;
      newRow[GT_ALLELE1_INDEX]     = allele1;
      newRow[GT_ALLELE2_INDEX]     = allele2;
      newRow[GT_ORIGINAL_ALLELE1_INDEX] = originalAllele1;

      newRow[GT_ORIGINAL_ALLELE2_INDEX] = originalAllele2;

      m_fileRowIndex[m_rows.size()] = i;

      recordsRead++;
      m_experimentRowIndex[m_experimentName].push_back(chpRow - 1);

      m_probeSetRowIndex[dbad.name] = m_rows.size();
      m_rows.push_back(newRow);

      if ( geneCopyNumber.count(newRow[GT_GENES_INDEX]) == 0 ) {
        geneCopyNumber[newRow[GT_GENES_INDEX]] = 2;
      }
      
      if ( allele1  == "0" ) {
        geneCopyNumber[newRow[GT_GENES_INDEX]]  = 0;
      }
      else if ( !CallElement::isWildCardBase(allele1) &&
                (geneCopyNumber[newRow[GT_GENES_INDEX]] == 2)
                && ( allele2.empty() ) ){
        geneCopyNumber[newRow[GT_GENES_INDEX]] = 1;
      }
    }

  }


  return recordsRead;
}
// end  TranslationInputStreamTableModel::_readNextExperimentCHP
/*****************************************************************************/
/*****************************************************************************/
/**
 * TranslationInputStreamTableModel::_readNextExperimentTsv
 * Synopsis:
 * Helper function to TranslationInputStreamTableModel::readNextExperiment.
 *
 * Stream experiments from the Genotype Short Report file one experiment
 * at a time. Fills in the 'm_rows' table with experiment data. 
 *
 * @param rte - the single instance run time environment. 
 * @param ttm  - the single instance translation table model
 * @param geneExperimentCopyNumberCall 
 *            - the DMET2 std::map of per gene copy number alleles
 * @param geneCopyNumber
 *            - the DMET3 std::map of gene copy number
 * 
 * @return - The number for records (markers) read for the experiment.
 */
/*****************************************************************************/
int TranslationInputStreamTableModel::_readNextExperimentTsv(const RunTimeEnvironment & rte, const TranslationTableModel & ttm, const std::map<std::string, std::string> &geneExperimentCopyNumberCall)
{

  bool   okExperimentGene   = true;
  int    recordsRead        = 0;
  std::string currentExperiment;

  clearData();

  if (m_tsvNextRow.size() > 0) {
    currentExperiment = m_tsvNextRow[m_tsvExperimentColumn];
    if (m_experimentRowIndex.count(currentExperiment) > 0) {
      Verbose::out(ADT_VERBOSE_NORMAL, "\n\nDMET2 genotype short report input files need to be sorted by experiment before they can be input.");
      sortTsvFileByExperiment(rte);
      APT_ERR_ABORT(currentExperiment + ": duplicate experiment detected.");
    }

    m_fileRowIndex[m_rows.size()] = m_tsvRow;
    m_rows.push_back(m_tsvNextRow);
    m_tsvNextRow.clear();
    recordsRead++;
  }

  std::string whiteSpace = " ";
  bool okColumnValues = true;


  for (int row = 0; okColumnValues && okExperimentGene &&
       (m_tsv.nextLevel(0) == TSV_OK); row++) {
    //
    std::vector<std::string> newRow;
    m_tsvRow++;

    newRow.resize(m_tsvColumnCount + 4, "");

    bool okRow = true;

    // Read in a row, one column at a time.
    for (int tsvColumn = 0; tsvColumn < m_tsvColumnCount; tsvColumn++) {

      std::stringstream lineSStr;
      std::string       columnName = getColumnName(m_tsvIndexToColumnDefinitionIndex[tsvColumn]);
      // get the data
      std::string columnValue;

      if (m_tsv.get(0, tsvColumn, columnValue) != TSV_OK) {
        APT_ERR_ABORT("tsv.get != TSV_OK");
      }


      Util::trimString(columnValue, whiteSpace.c_str());

      // empty and allowed?
      if (columnValue.empty()) {
        if (m_tsvColumnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_emptyOk == false) {
          lineSStr << m_tsvFileName << ": ROW " << m_tsvRow  << " : ";
          lineSStr << columnName << " : ERROR missing required column.";
          Verbose::out(rte.m_currentVerbosity, lineSStr.str());
          okColumnValues = false;
        }
      }
      // validate it against regexp
      else if (!m_tsvColumnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_validRE.FullMatch(columnValue)) {
        lineSStr << m_tsvFileName << ": ROW " << m_tsvRow  << " : COLUMN [";
        lineSStr << tsvColumn << "] " << columnName << " : ERROR invalid value: \"";
        lineSStr << columnValue << "\"";
        lineSStr << " (expecting " << m_tsvColumnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_validRE.pattern();
        lineSStr << ")";
        lineSStr << endl;

        for (int j = 0 ; j < m_tsvColumnCount;  j++) {
          m_tsv.get(0, j, columnValue);
          lineSStr << columnValue << ",";
        }
        Verbose::out(rte.m_currentVerbosity, lineSStr.str());
        okColumnValues = false;
      }
      // Ingore row
      else if ((m_tsvColumnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_ignoreRE != NULL) &&
               (m_tsvColumnDefinition[m_tsvIndexToColumnDefinitionIndex[tsvColumn]].m_ignoreRE->PartialMatch(columnValue))) {
        okRow = false;

      }
      // put the column data in our row.
      newRow[m_tsvIndexToColumnDefinitionIndex[tsvColumn]] = columnValue;
    }
    // save our row to the data.

    // Genotype Override
    if (m_gotm != NULL) {
      std::vector< std::string > override = m_gotm->getOverrideAlleles(newRow[m_tsvExperimentColumn],  newRow[m_tsvProbeSetColumn]);
      if (override.size() ) {
        cerr << newRow[m_tsvExperimentColumn] << " " << newRow[m_tsvProbeSetColumn] << " has override " << override[0] << endl;
        m_gotm->m_chpFileProbeSetOriginalBasecall[newRow[m_tsvExperimentColumn] + "::" + newRow[m_tsvProbeSetColumn]] = newRow[GT_ALLELE1_INDEX];
        newRow[GT_ALLELE1_INDEX] = override[0];
        newRow[GT_ALLELE2_INDEX] = "";
        if ( override.size() == 2 ) {
          newRow[GT_ALLELE2_INDEX] = override[1];
        }
      }
    }
    // DMET2 Copy Number Zero
    if ( ttm.getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3 ) {
      std::string key = newRow[GT_GENES_INDEX] + ":" + newRow[GT_EXPERIMENT_INDEX];
      if ( geneExperimentCopyNumberCall.count(key) > 0 ) {
        newRow[GT_ALLELE1_INDEX] = "0";
        newRow[GT_ALLELE2_INDEX] = "";
      }
        
    }
    
    if (okRow) {

      if (currentExperiment.empty()) {
        currentExperiment = newRow[m_tsvExperimentColumn];
      }
      if (currentExperiment == newRow[m_tsvExperimentColumn]) {
        m_fileRowIndex[m_rows.size()] = m_tsvRow;
        recordsRead++;
        m_experimentRowIndex[currentExperiment].push_back(m_tsvRow - 1);

        m_rows.push_back(newRow);
      } else {
        m_tsvNextRow = newRow;
        okExperimentGene = false;
      }
    }
  }


  if (! okColumnValues) {
    APT_ERR_ABORT(m_tsvFileName + ": invalid input file data detected.");
  }

  if (m_rows.size() != recordsRead) {
    m_rows.resize(recordsRead);
  }

  return recordsRead;

}
// end TranslationInputStreamTableModel::_readExperimentNextGeneTsv
/*****************************************************************************/
