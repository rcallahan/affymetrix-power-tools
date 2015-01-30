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
 * @file   GenotypeOverrideTableModel.cpp
 * @author Mybrid Spalding
 * @date   Tue Jul 22 08:33:25 PDT 2008
 * @brief  Class to supplant Genotype data (CHP) with user supplied data.
 */

//
#include "translation/GenotypeOverrideTableModel.h"
//
#include "translation/CallElement.h"
#include "translation/RunTimeEnvironment.h"
#include "translation/TranslationTableModel.h"
//
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include "pcrecpp.h"
//
#include <sstream>

using namespace affx;
using namespace std;

const TittmColumnDefinition GOTM_COLUMN_DEFINITIONS[] = { // column_name, empty_ok, index, valid_regex
  { std::string("CHP Filename"),      0, 0, std::string(".*"), NULL },
  { std::string("Gene"),              0, 1,
    std::string("^[A-Z\\d]+$"), NULL},
  { std::string("Common Name"),       0, 2, std::string(".*"), NULL  },
  { std::string("Probe Set ID"),      0, 3,
    std::string("^(?:[\\w\\d]+)$"), NULL },
  { std::string("Basecall"),          1, 4,
    std::string("^(?:0|ZeroCopyNumber|PossibleRareAllele|NoCall|NotAvailable|PRA|NC|-)|(?:(?:[ACGT0]+|INS|Ins|DEL|Del|-)/(?:[ACGT]+|INS|Ins|DEL|Del|-))$"), NULL },
  { std::string("Override Comment"),  1, 5, std::string(".*"), NULL },
  { std::string("Reference Allele"),  0, 6,
    std::string("^(?:[ACGT0]+)|(?:INS)|(?:DEL)|(?:-)$"), NULL },
  { std::string("Variant Allele"),    0, 7,
    std::string("^(?:[ACGT]+|INS|Ins|DEL|Del|-)(?:[,](?:[ACGT]+|INS|Ins|DEL|Del|-))*$"), NULL },
};

/*****************************************************************************/
/**
 * filterRowCallback:
 * Synopsis:
 *
 *  A callback function passed to TranslationInputTsvTableModel.
 *
 *  DEPRECATED: As of 10/08/2008 filtering is deprecated. If someone
 *  wants to override a specific call with a NoCall they can now do so.
 *  This function filters out NoCall values such that they would never
 *  be used from the override file.
 *
 *
 *
 * @param name1 - description
 * @param name1 - description
 * @return - description
 */
/*****************************************************************************/
/*
static bool filterRowCallback(std::vector< std::string > & row)
{

  pcrecpp::RE filter("No[cC]all|PossibleRareAllele|NC|PRA");

  APT_ERR_ASSERT(row.size() == GenotypeOverrideTableModel::END_GOTM_FIELD_INDEX, "");

  if (row[GenotypeOverrideTableModel::BASECALL_INDEX].empty() ||
      filter.PartialMatch(row[GenotypeOverrideTableModel::BASECALL_INDEX])) {
    return true;
  }


  return false;

}
*/
// end filterRowCallback
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeOverrideTableModel::GenotypeOverrideTableModel
 * Synopsis:
 *
 * This is the TransferObject to the TsvFile DataAccessObject design
 * pattern. There is only ever one of these in existence, i.e. this
 * is a single instance class.
 *
 * genotype override file -> genotype override table model
 *
 * Default TSV constructor for the GenotypeOverrideTableModel.
 *
 * This constructor is beefy for one reason: to validate
 * rows and then to swap the row information with error messaging
 * in the event that the row is invalid.
 *
 * Unlike most other input files for the translation algorithm,
 * the override file is permitted to have bad data and translation
 * still runs. When bad rows are detected, an error message is stuffed
 * into the override comment field and the override row is otherwise
 * ignored.
 *
 * Validation occurs per allele or probe set with the exception of
 * Zero Copy Number. With Zero Copy Number all probe sets for a gene
 * must be present or the override is rejected.
 * 
 * @param rte - RunTimeEnvironment which contains options
 * @param genoOverrideFileName - the input file to read in.
 * @param ttm - The translation table instance used for validation.
 */
/*****************************************************************************/
GenotypeOverrideTableModel::GenotypeOverrideTableModel(const class RunTimeEnvironment &rte,  const std::string &genoOverrideFileName, class TranslationTableModel & ttm) :
    TranslationInputTsvTableModel(rte, genoOverrideFileName, GOTM_COLUMN_DEFINITIONS, (sizeof(GOTM_COLUMN_DEFINITIONS) / sizeof(GOTM_COLUMN_DEFINITIONS[0])), true, NULL, NULL, false)
{

  std::vector< std::string > alleles;
  std::set< std::string >  copyNumberZeroGenes;
  std::map< std::string, std::set< std::string > > geneMarkers;
  std::stringstream        invalidSStr;
  pcrecpp::RE              zeroCopyRE("ZeroCopyNumber|0");
  pcrecpp::RE              zeroCopyValidSetRE("ZeroCopyNumber|0|NotAvailable");
  
  for (size_t row = 0; row < m_rows.size(); row++) {

    // VALIDATION
    // We requested that validation not be done by the constructor so
    // that special action can be taken here with respect to validation.
    // Specifically, invalid rows are appended to a reporting structure
    // and does not cause the translation to stop.

    bool okRow         = true;
    bool okRowProbeSet = true;
    std::string probeSet    = m_rows[row][PROBE_SET_INDEX];
    std::string gene        = m_rows[row][GENE_INDEX];
    std::string chpFile     = m_rows[row][CHP_FILE_INDEX] ;
    std::string key         = chpFile + "::" + probeSet;
    std::string geneKey     = chpFile + "::" + gene;
    alleles.clear();
    
    invalidSStr.str("");
      
    invalidSStr << "Invalid basecall override ignored";

    if ( (geneMarkers.count( geneKey )  == 0 ) || ( !geneMarkers[geneKey].count(probeSet)) )  {
      geneMarkers[geneKey].insert(probeSet);
    }
    else {
      okRow = false;
      invalidSStr << " | " << chpFile << " " << gene << " " << " duplicate  probe set " << probeSet << " detected." <<  endl;
    }

    for (int col = 0; (col < columnCount()) && okRow; col++) {

      std::string columnValue = m_rows[row][col];
        
      if (! GOTM_COLUMN_DEFINITIONS[col].m_validRE.FullMatch(columnValue)) {
        okRow = false;

        invalidSStr << " | Unrecognized " << GOTM_COLUMN_DEFINITIONS[col].m_columnName << " (" << columnValue << ") " ;

        if (col == PROBE_SET_INDEX) {
          okRowProbeSet = false;
        }

      } // if column is valid

    } // for each column


    // BASECALL VALIDATION
    // Make sure the base call is comprised of values found in the
    // reference and variant columns.

    if (okRow) {

      alleles = getOverrideAlleles(m_rows[row][CHP_FILE_INDEX], m_rows[row][PROBE_SET_INDEX], row);

      if ((alleles.size() != 2) && !zeroCopyValidSetRE.FullMatch(alleles[0] ) ) {
        okRow = false;
        invalidSStr << " | basecall format must be a copy number two designation  allele/allele " << m_rows[row][BASECALL_INDEX];
      }
 
      for (size_t i = 0; okRow && (i < alleles.size()); i++) {

        if (!CallElement::isWildCardBase(alleles[i]) &&
            !zeroCopyRE.FullMatch(alleles[i]) &&
            !(alleles[i] == m_rows[row][REFERENCE_INDEX]) &&
            !(pcrecpp::RE(alleles[i]).PartialMatch("((?:^)|(?:,))" + m_rows[row][VARIANT_INDEX] + "(?:(?:,)|(?:$))"))) {
          okRow = false;
          invalidSStr << " | basecall " << m_rows[row][BASECALL_INDEX] << " has alleles not found in either the reference {" << m_rows[row][REFERENCE_INDEX] <<  "} or variant list {" <<  m_rows[row][VARIANT_INDEX] << "} ";
        }
      }
    }


    // TRANSLATION TABLE VALIDATION
    if (okRow) {
      bool geneValid     = false;
      bool variantValid  = false;
      bool referenceValid = false;
      std::string rowGene;
      std::string rowReference;
      std::string rowVariants;

      if (ttm.m_probeSetRowIndex.count(probeSet) == 0) {
        okRow = false;
        okRowProbeSet = false;
        invalidSStr << " | "  << rte.m_adtOpts.m_inputTTableFile << " does not contain probe set(" << probeSet << ") ";
      } else {
        for (size_t i = 0 ; i < ttm.m_probeSetRowIndex[probeSet].size(); i++) {

          int ttRow = ttm.m_probeSetRowIndex[probeSet][i];
          if (m_rows[row][GENE_INDEX] == ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_GENE)]) {
            geneValid = true;
          }
          rowGene = ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_GENE)];
          if (m_rows[row][REFERENCE_INDEX] == ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE)]) {
            referenceValid = true;
          }
          rowReference = ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE)];
          if (pcrecpp::RE(ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)]).PartialMatch("((?:^)|(?:,))" + m_rows[row][VARIANT_INDEX] + "(?:(?:,)|(?:$))")) {
            variantValid = true;
          }
          if (rowVariants.empty()) {
            rowVariants = ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)];
          } else {
            rowVariants = rowVariants + "," + ttm.m_rows[ttRow][ttm.getColumnIndex(ADT_DMET3_TT_VARIANT)];
          }
        } // for each allele of a possible multiallelic probe set.

        if (!geneValid) {
          okRow = false;
          invalidSStr << " | "  << rte.m_adtOpts.m_inputTTableFile << " probe set " << probeSet << " gene " << rowGene << " does not match (" << m_rows[row][GENE_INDEX] << ") ";

        }
        if (!referenceValid) {
          okRow = false;
          invalidSStr << " | "  << rte.m_adtOpts.m_inputTTableFile << " probe set " << probeSet << " reference " << rowReference << " does not match (" << m_rows[row][REFERENCE_INDEX] << ") ";

        }
        if (!variantValid) {
          okRow = false;
          invalidSStr << " | "  << rte.m_adtOpts.m_inputTTableFile << " probe set " << probeSet << " variant set {" << rowVariants << "} does not contain (" << m_rows[row][REFERENCE_INDEX] << ") ";

        }

        if ( okRow && alleles.size() && zeroCopyRE.FullMatch( alleles[0] ) &&
             !(ttm.getCopyNumberColumn(ttm.m_probeSetRowIndex[probeSet][0]) > 0) ) {
            okRow = false;
            invalidSStr << " | "  << Fs::basename(chpFile)
                        << " gene probe set " << gene << ", " << probeSet 
                        << " given ZeroCopyNumber override but has no such copy number zero designation in the translation table.";
        }
      }
    }

    if (okRow) {
      m_chpFileProbeSetToRowIndex[ key ] = row;
      if ( zeroCopyRE.FullMatch( alleles[0] ) && !copyNumberZeroGenes.count(geneKey) ) {
        copyNumberZeroGenes.insert(geneKey);
      }
    } else {
      if (okRowProbeSet) {

        m_chpFileProbeSetInvalidBasecallComment[key].push_back(m_rows[row][BASECALL_INDEX]);
        m_chpFileProbeSetInvalidBasecallComment[key].push_back(invalidSStr.str());
        
      }

      invalidSStr << " ROW {" << getRowAsString(row) << "}";

      Verbose::out(ADT_VERBOSE_NORMAL, invalidSStr.str());
    }

  } // for each row

  if ( copyNumberZeroGenes.size() == 0 ) {
    return;
  }

  // ZERO COPY NUMBER VALIDATION
  
  // Validate that all copy number zero genes have the same set of markers
  // as the translation table.
  // Note there is an important order of occurence here.
  // The marker list which filters the translation table needs to be applied
  // to the translation table before this constructor is called.
  // The m_probeSetIndex is maintained to only have the filtered probe sets
  // and not the complete set of probe sets. The translation table contains
  // all probe sets where the filtered probe sets have been ignored or
  // commented out.

  std::map< std::string, std::string > validCopyNumberZeroCheck;
  bool invalidCopyNumberZeroGenesExist = false;

  for (size_t row = 0; row < m_rows.size(); row++) {

    std::string probeSet    = m_rows[row][PROBE_SET_INDEX];
    std::string gene        = m_rows[row][GENE_INDEX];
    std::string chpFile     = m_rows[row][CHP_FILE_INDEX] ;
    std::string key         = chpFile + "::" + probeSet;
    std::string geneKey     = chpFile + "::" + gene;

    if ( ! copyNumberZeroGenes.count(geneKey) || ! m_chpFileProbeSetToRowIndex.count( key ) || validCopyNumberZeroCheck.count(geneKey) ) {
      continue;
    }

    alleles = getOverrideAlleles(m_rows[row][CHP_FILE_INDEX], m_rows[row][PROBE_SET_INDEX], row);

    if ( !zeroCopyValidSetRE.FullMatch(alleles[0]) ) {
        invalidSStr.str("");
        invalidSStr << "Invalid basecall override ignored | ";
        invalidSStr << Fs::basename(chpFile) << " gene " << gene << " probe set " << probeSet << " is not ZeroCopyNumber, but " << m_rows[row][BASECALL_INDEX] << ". ZeroCopyNumber override requires all probe sets to be present and ZeroCopyNumber basecall." ;
        validCopyNumberZeroCheck[geneKey] = invalidSStr.str();
        invalidCopyNumberZeroGenesExist = true;
        break;

    }

    if ( geneMarkers[geneKey].size() < ttm.m_geneMarkers[gene].size() ) {
      invalidSStr.str("");
      invalidSStr << "Invalid basecall override ignored | ";
      invalidSStr << " Gene " << gene << " ZeroCopyNumber requires " << ttm.m_geneMarkers[gene].size() << " probe set markers but only " << geneMarkers[geneKey].size() << " probe set markers were provided.";
      validCopyNumberZeroCheck[geneKey] = invalidSStr.str();
      invalidCopyNumberZeroGenesExist = true;
      continue;
    }      

    
    if ( ! validCopyNumberZeroCheck.count(geneKey) ) {
      for (size_t i = 0; i < ttm.m_geneMarkers[gene].size(); i++ ) {
        if ( ! geneMarkers[geneKey].count( ttm.m_geneMarkers[gene][i] ) ) {
          invalidSStr.str("");
          invalidSStr << "Invalid basecall override ignored | ";
          invalidSStr << Fs::basename(chpFile) << " gene " << gene << " ZeroCopyNumber markers provided differ from the translation table where " << ttm.m_geneMarkers[gene][i] << " was not found.";
          validCopyNumberZeroCheck[geneKey] = invalidSStr.str();
          invalidCopyNumberZeroGenesExist = true;
        }
      }
    }

  }

  // Remove the overrides from the valid list and set up the invalid
  // message.
  
  if ( invalidCopyNumberZeroGenesExist ) {

    for (size_t row = 0; row < m_rows.size(); row++) {

      std::string probeSet    = m_rows[row][PROBE_SET_INDEX];
      std::string gene        = m_rows[row][GENE_INDEX];
      std::string chpFile     = m_rows[row][CHP_FILE_INDEX] ;
      std::string key         = chpFile + "::" + probeSet;
      std::string geneKey     = chpFile + "::" + gene;

      if ( validCopyNumberZeroCheck.count(geneKey) &&
           (validCopyNumberZeroCheck[geneKey] != std::string("VALID")) ) {

        if ( m_chpFileProbeSetToRowIndex.count(key)) {
          m_chpFileProbeSetToRowIndex.erase( key );
          m_chpFileProbeSetInvalidBasecallComment[key].push_back(m_rows[row][BASECALL_INDEX]);
          m_chpFileProbeSetInvalidBasecallComment[key].push_back(validCopyNumberZeroCheck[geneKey]);
        }
      }
    }

  }

  return;

}
// end GenotypeOverrideTableModel::GenotypeOverrideTableModel
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeOverrideTableModel::describeVerbose:
 * Synopsis:
 *
 * Run time debug routine that can be invoked with various verbose levels.
 * ADT_VERBOSE_INPUT_FILES is the default level.
 *
 * @param rte       - run time environment
 * @param override  - ADT_VERBOSE_INPUT_FILES by default,
 *                    otherwise the ADT_VERBOSE_ENUM level to output
 *
 */

/*****************************************************************************/
void GenotypeOverrideTableModel::describeVerbose(ADT_VERBOSE_ENUM level)
{

  std::stringstream msgSStr;

  msgSStr << "GenotypeOverrideTableModel::describeVerbose" << endl;

  msgSStr << "m_chpFileProbeSetToRowIndex:" << endl;

  std::map< std::string, int >::iterator itSI = m_chpFileProbeSetToRowIndex.begin();

  for (; itSI != m_chpFileProbeSetToRowIndex.end(); itSI++) {

    msgSStr << itSI->first << " row [" << itSI->second << "]" << endl;

  }

  std::map< std::string, std::string >::iterator itSS = m_chpFileProbeSetOriginalBasecall.begin();

  msgSStr << "m_chpFileProbeSetOriginalBasecall:" << endl;
  for (; itSS != m_chpFileProbeSetOriginalBasecall.end(); itSS++) {

    msgSStr << itSS->first << " : " << itSS->second << endl;

  }


  msgSStr << "m_chpFileProbeSetInvalidBasecallComment:" << endl;

  std::map< std::string,  std::vector< std::string > >::iterator itSVS
  = m_chpFileProbeSetInvalidBasecallComment.begin();

  for (; itSVS != m_chpFileProbeSetInvalidBasecallComment.end(); itSVS++) {

    msgSStr << itSVS->first << " : ";
    for (size_t i = 0; i < itSVS->second.size(); i++) {
      msgSStr << itSVS->second[i] << ", ";
    }

    msgSStr << endl;
  }

  Verbose::out(level, msgSStr.str());

}
// end GenotypeOverrideTableModel::describeVerbose
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeOverrideTableModel::getOverrideAlleles:
 * Synopsis:
 *
 *    If an override exists, returns a non-empty std::vector pair of std::strings
 * corresponding to the override basecall alleles.
 *
 * NOTE: Override files only have "NoCall", "PossibleRareAllele" or "
 * "NotAvailable" and not the format of "NoCall/NoCall" as
 * would be returned by the GenoCallCoder. Therefore the return result
 * for wildcard overrides are reformated to the GenoCallCoder format
 * ("NoCall/NoCall").
 *
 *
 *
 * @param chpFileName - experiment chp file name (minus .CHP)
 * @param probeSetId  - probeset id in question.
 * @param row         - the row in "m_rows".
 *
 * @return - alleles as a std::string pair inside a std::vector.
 */
/*****************************************************************************/
std::vector< std::string > GenotypeOverrideTableModel::getOverrideAlleles(const std::string & chpFileName, const std::string & probeSetId, int row)
{

  std::vector< std::string > results;

  if (chpFileName.empty() || probeSetId.empty()) {
    return results;
  }
  bool printResults = false;
  
  if (row < 0) {

    printResults = true;
    std::string key = chpFileName + "::" + probeSetId;

    std::map< std::string , int >::const_iterator itSI;
    itSI = m_chpFileProbeSetToRowIndex.find(key) ;

    if (itSI == m_chpFileProbeSetToRowIndex.end()) {
      return results;
    }
    row = itSI->second;
  }

  std::string basecall = m_rows[row][BASECALL_INDEX];

  pcrecpp::RE split("([^/]+)/([^/]+)");

  std::string allele1;
  std::string allele2;

  if (split.FullMatch(basecall, &allele1, &allele2)) {
    results.push_back(allele1);
    results.push_back(allele2);
  }
  
  if (pcrecpp::RE("0|ZeroCopyNumber|PossibleRareAllele|NoCall|NotAvailable|NC|PCRA|NA").FullMatch(basecall)) {
    if ( basecall  == "ZeroCopyNumber" ) {
      basecall = "0";
    }
    results.push_back(basecall);
    if ( basecall != "0" ) {
      results.push_back(basecall);
    }
  }


  return results;

}
// end GenotypeOverrideTableModel::getOverrideAlleles
/*****************************************************************************/
/*****************************************************************************/
/**
 * GenotypeOverrideTableModel::getOverrideInfo:
 * Synopsis:
 *
 *    If an override exists, returns a non-empty std::vector PAIR of std::strings
 * corresponding to the override basecall original basecall and the
 * override comment (override basecall, override comment).
 *
 *
 *
 * @param chpFileName - experiment chp file name (minus .CHP)
 * @param probeSetId  - assay id in question.
 * @param includeInvalidOriginalBaseCall - boolean that returns that includes the original base call, even for invalid calls. 
 *
 * @return - (original basecall, commment) as a std::vector pair of std::strings.
 */
/*****************************************************************************/
std::vector< std::string > GenotypeOverrideTableModel::getOverrideInfo(const std::string & chpFileName, const std::string & probeSetId, bool includeInvalidOriginalBaseCall )
{

  std::vector< std::string > results;

  if (chpFileName.empty() || probeSetId.empty()) {
    return results;
  }

  std::string key = chpFileName + "::" + probeSetId;

  std::map< std::string , int >::const_iterator itSI;

  itSI = m_chpFileProbeSetToRowIndex.find(key) ;


  if (itSI == m_chpFileProbeSetToRowIndex.end()) {

    // INVALID uncalled file row check.
    std::map< std::string , std::vector< std::string >  >::const_iterator jtSVS;
    jtSVS =  m_chpFileProbeSetInvalidBasecallComment.find(key);

    if (jtSVS == m_chpFileProbeSetInvalidBasecallComment.end()) {
      // COMMON CASE, most probe set rows will not be in the uncalled file.
      return results;
    }

    // ROW was INVALID
    APT_ERR_ASSERT(jtSVS->second.size() == 2, "");

    if ( includeInvalidOriginalBaseCall ) {
      results.push_back(jtSVS->second[0]);
    }
    else {
      results.push_back("");
    }

    results.push_back(jtSVS->second[1]);

    return results;
  }

  std::string comment = m_rows[itSI->second][OVERRIDE_COMMENT_INDEX];

  std::map< std::string, std::string >::const_iterator itSS;

  itSS = m_chpFileProbeSetOriginalBasecall.find(key);

  if (itSS == m_chpFileProbeSetOriginalBasecall.end()) {
    return results;
  }

  std::string basecall = itSS->second;

  results.push_back(basecall);
  results.push_back(comment);

  return results;

}
// end GenotypeOverrideTableModel::getOverrideInfo
/*****************************************************************************/
