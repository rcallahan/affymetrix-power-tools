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
 * @file   GeneCall.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:04:13 PDT 2008
 * @brief  Business object that is the heart of the matter, class with the algorithm that translates experiment data.
 */

#include "translation/GeneCall.h"
//
#include "translation/CallResults.h"
#include "translation/GenotypeTableModel.h"
#include "translation/TranslationTableModel.h"
//
#include "util/Err.h" // includes "util/Verbose.h"
//
#include <sstream>
//


using namespace std;


/*****************************************************************************/
/**
 * GeneCall::GeneCall
 * Synopsis:
 *
 * Translation Table constructer. The translation table row order is
 * per marker. The column order starting with Reference, Variant
 * and then A1 - ANN is what is being defined here as a call.
 * The translation table is rows of markes with columns of calls.
 *
 * Each gene has a number of calls designated. The GeneCall object
 * represents the entire set of calls for a gene. The GeneCall is
 * a set of sets, where the sets are CallSets. The CallSets are sets
 * of CallElements where the CallElements represent the probset id
 * with reference and variant alleles.
 *
 * A GeneCall then represents the set of columns A1 - A50 (ANN) of the
 * translation table for a particular gene.
 *
 * The columns are chopped up at the Gene level. Within
 * a gene, the column is subdivided by haplotype or marker.
 * A column, like A1, then will have one GeneCall per marker where
 * the "Haplotype" column is "N". For those markers of the gene
 * where "Haplotype" column is "Y" then exactly one GeneCall is created.
 *
 * NOTE: While in general the "Haplotype" column is "Y" records are grouped
 * together for a particular gene in the translation table, this is not
 * always the case. Probe sets with "Haplotype" column of "Y" can be
 * interpersed with the "Haplotype" column "N" probesets.
 *
 * For a particular gene there is ever only one haplotype group.
 * That group is designated by the set of markers for the gene with
 * "Haplotype" column indicator of "Y".
 * A gene haplotype call (GeneCall) then corresponds to a column for
 * all those rows of a gene in the translation table within a haplotype group.
 *
 * This constructor takes in the header row for a particular gene
 * and then instantiates all the CallSets for all the columns.
 *
 * The CallSets for each call is subsequently filled with the probe set ids,
 * reference and variant alleles when the translation table is scanned.
 *
 * In other words, multiple data structures are being created when the
 * translation table is being scanned. To prevent scanning the table
 * multiple times this constructor simply instantiates the empty CallSet
 * per call and the callee (TranslationTable) is responsible for filling
 * in the CallSets per call as the subsequent rows are being scanned.
 *
 *
 * @param rte - RunTimeEnvironment which contains options
 * @param ttm - TranslationTableModel with the meta data about the call.
 * @param headerRow - the row in the ttm that has the meta data.
 * @param isHaplotypeCall - if the data row has been predetermined to
 *   have "Y" in the "Haplotype" column.
 *
 */
/*****************************************************************************/
GeneCall::GeneCall(const RunTimeEnvironment & rte,
                   TranslationTableModel & ttm,
                   const int headerRow,
                   bool isHaplotypeCall)

{

  APT_ERR_ASSERT((headerRow >= 0) && (headerRow < ttm.size()), "");

  m_gene = ttm.m_rows[headerRow][ttm.getColumnIndex(ADT_DMET3_TT_GENE)];

  Verbose::out(ADT_VERBOSE_INPUT_FILES, "GeneCall::GeneCall initializing gene: " + m_gene);

  // Initialize the m_alleleSet for each possible call.

  std::stringstream alleleNamesSStr;

  pcrecpp::RE reCommentedOut("^\\s*#");
  // Skip commented out column names.
  std::map<int, int> columnToAlleleIndex;
  std::map<std::string, int> numHaplotypeMarkers;


  for (int i = ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE); i < ttm.columnCount(); i++) {


    // For markers, only add the allele if the header is copy number 0.
    // We can only tell if the header is copy number zero if the data
    // in the following row equals "0";
    if (!isHaplotypeCall && (i >= ttm.getColumnIndex(ADT_DMET3_TT_ALLELE_START)) && (((headerRow + 1) >= ttm.size()) || (ttm.m_rows[headerRow +1][i] != "0"))) {
      continue;
    }

    std::string alleleName = ttm.getHeaderRowColumnName(headerRow, i);

    if (existsAlleleCallSet(alleleName)) {
      // Copy Number 0 applies to all markers, ignore it.
      if (((headerRow + 1) < ttm.size()) && (ttm.m_rows[headerRow +1][i] == "0")) {
        continue;
      }

      Verbose::warn(ADT_VERBOSE_NORMAL, "Duplicate allele name found: " + alleleName + " for translation file row: " + ttm.getRowAsString(headerRow));
    }

    // All allele name columns will have a name until we get to the last one.
    if (!alleleName.empty() && !reCommentedOut.PartialMatch(alleleName)) {
      columnToAlleleIndex[i] = m_alleleSet.size();
      m_alleleSet.push_back(CallSet(alleleName, i));
      alleleNamesSStr << alleleName << " ";
      numHaplotypeMarkers[alleleName] = 0;
    }
  }

  if (isHaplotypeCall) {
    for (int i = headerRow; (i < ttm.size()) && (ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_GENE)] == m_gene); i++) {
      for (int j = ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE);
           j < ttm.columnCount(); j++) {

        if (!ttm.m_rows[i][j].empty() &&
            (ttm.m_rows[i][ttm.getColumnIndex(ADT_DMET3_TT_HAPLOTYPE)] == "Y")) {
          numHaplotypeMarkers[m_alleleSet[columnToAlleleIndex[j]].m_name]++;
        }
      }
    }
    // If an Allele doesn't belong to a haplotype group, then just remove it.
    for (int i = m_alleleSet.size() - 1; i >= 0 ; i--) {
      if ((m_alleleSet[i].m_columnInTranslationTable >
           ttm.getColumnIndex(ADT_DMET3_TT_ALLELE_START)) &&
          (numHaplotypeMarkers[m_alleleSet[i].m_name] == 0)) {
        m_alleleSet.erase(m_alleleSet.begin() + i);
      }
    }
  }

  Verbose::out(ADT_VERBOSE_INPUT_FILES, "GeneCall::GeneCall gene: " + m_gene + ": initializing allele names: " + alleleNamesSStr.str());
  //describeVerbose(rte, ADT_VERBOSE_NORMAL );

  return;

}
// end GeneCall::GeneCall
/*****************************************************************************/
/*****************************************************************************/
/**
 * GeneCall::addCallElementToAlleleCallSet
 * Synopsis:
 *
 * Initialize a CallElement from the translation table (and not
 * an experiment data in the Genotype data (CHP)).
 * To do this we need to identify the ProbeSet and
 * the reference and variant base. The variant base will be either in the
 * column that resides between A1-A29 or in the "Variant" column itself
 * if no such animal exists.
 *
 *
 * @param ttm - TranslationTableModel with the meta data about the call.
 * @param headerRow - the row in the ttm that has the meta data.
 * @param row - the row in the TranslationTableModel in question.
 * @param column - the column in the TranslationTableModel in question.
 * @param allowMultiAllelic - probeSet may appear multiple times.
 * @param rowIsHaplotype - "Haplotype" column has "Y"
 *
 * @return - true upon successfull add.
 */
/*****************************************************************************/
bool GeneCall::addCallElementToAlleleCallSet(TranslationTableModel & ttm,
    const int headerRow,
    const int row,
    const int column,
    const bool allowMultiAllelic,
    const bool rowIsHaplotype)

{

  APT_ERR_ASSERT((headerRow >= 0) && (headerRow < ttm.size()) &&
                 (row >= 0) && (row < ttm.size()) &&
                 (headerRow < row), "");
  APT_ERR_ASSERT((column >= ttm.getColumnIndex(ADT_DMET3_TT_REFERENCE)) &&
                 (column < ttm.columnCount()), "");


  std::string alleleName = ttm.getHeaderRowColumnName(headerRow, column);

  std::vector< class CallSet >::iterator itCS;

  for (itCS = m_alleleSet.begin();
       itCS != m_alleleSet.end() ; ++itCS) {

    if (itCS->m_name == alleleName) {
      itCS->m_type = rowIsHaplotype ? ADT_CALL_TYPE_HAPLOTYPE_GROUP :
                     ADT_CALL_TYPE_MARKER;
      itCS->m_isDescriptive = true;
      return itCS->addCallElement(ttm, row, column, allowMultiAllelic, rowIsHaplotype);
    }
  }

  std::stringstream msgSStr;

  msgSStr << "addCallElementToAlleleCallSet: can't find allele name: " <<  alleleName <<  " from the following possibilities: ";
  for (itCS = m_alleleSet.begin();
       itCS != m_alleleSet.end() ; ++itCS) {
    msgSStr << itCS->m_name << " ";
  }

  Verbose::out(ADT_VERBOSE_NORMAL, msgSStr.str());

  for (itCS = m_alleleSet.begin();
       itCS != m_alleleSet.end() ; ++itCS) {
    //itCS->describeVerbose(ADT_VERBOSE_NORMAL );
  }

  return false;

}
// end GeneCall::addCallElementToAlleleCallSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * _getNumDescriptive
 * Synopsis:
 *
 * Helper function for describeVerbose only.
 * Return the number of CallSets in the m_alleleSet which are descriptive.
 *
 *
 * @return int - the number of descriptive
 */
/*****************************************************************************/
static int _getNumDescriptive(const std::vector<CallSet> & alleleSet)
{

  int i = 0;
  int numDescriptive = 0;

  for (i = 0, numDescriptive = 0 ; i < alleleSet.size(); i++) {
    if (alleleSet[i].m_isDescriptive)
      numDescriptive++;
  }

  return numDescriptive;

}
// end _getNumDescriptive
/*****************************************************************************/
/*****************************************************************************/
/**
 * _getNumNonDescriptive
 * Synopsis:
 *
 * Helper function for describeVerbose only.
 * Return the number of CallSets in the m_alleleSet which are not descriptive.
 *
 *
 * @return int - the number of non descriptive
 *
 */
/*****************************************************************************/
static int _getNumNonDescriptive(const std::vector<CallSet> & alleleSet)
{
  return (alleleSet.size() - _getNumDescriptive(alleleSet));
}
// end _getNumNonDescriptive
/*****************************************************************************/
/*****************************************************************************/
/**
 * GeneCall::describeVerbose
 * Synopsis:
 *
 * Run time debug routine that can be invoked with various verbose levels.
 * ADT_VERBOSE_INPUT_FILES is the default level.
 *
 * @param rte       - run time environment
 * @param level     - ADT_VERBOSE_INPUT_FILES by default,
 *                    otherwise the ADT_VERBOSE_ENUM level to output
 *
 */
/*****************************************************************************/
void GeneCall::describeVerbose(const RunTimeEnvironment & rte, ADT_VERBOSE_ENUM level)
{

  if (!level && (rte.m_currentVerbosity < ADT_VERBOSE_INPUT_FILES))
    return;

  ADT_VERBOSE_ENUM outLevel = level ? level : ADT_VERBOSE_INPUT_FILES;

  std::stringstream describe1SStr;

  describe1SStr << "GeneCall::describeVerbose gene: " << m_gene << ":";
  describe1SStr << " alleleSet total count : ";
  describe1SStr << m_alleleSet.size() << endl;
  describe1SStr << "GeneCall::describeVerbose gene: " << m_gene << ":";
  describe1SStr << " alleleSet descriptive count : ";
  describe1SStr << _getNumDescriptive(m_alleleSet) << endl;
  describe1SStr << "GeneCall::describeVerbose gene: " << m_gene << ":";
  describe1SStr << " alleleSet non-descriptive count : ";
  describe1SStr << _getNumNonDescriptive(m_alleleSet);
  Verbose::out(outLevel, describe1SStr.str());

  std::stringstream describe4SStr;
  describe4SStr << "GeneCall::describeVerbose gene: " << m_gene << ": ";
  describe4SStr << "vector<CallSet> alleleSet follows: " ;

  Verbose::out(outLevel, describe4SStr.str());

  int i;
  for (i = 0; i < m_alleleSet.size(); i++) {
    m_alleleSet[i].describeVerbose(rte, level);
  }

  return;

}
// end GeneCall::describeVerbose
/*****************************************************************************/
/*****************************************************************************/
/**
 * GeneCall::existsAlleleCallSet
 * Synopsis:
 *
 * Returns true if a CallSet in m_alleleSet has the passed in allele name.
 *
 * @params alleleName - the call name to check (i.e. "Ref", "*4")
 *
 */
/*****************************************************************************/
bool GeneCall::existsAlleleCallSet(const std::string & alleleName)
{

  std::vector<CallSet>::iterator itCS;

  for (itCS = m_alleleSet.begin(); itCS != m_alleleSet.end(); itCS++) {
    if (itCS->m_name == alleleName) {
      return true;
    }
  }

  return false;
}
// end GeneCall::existsAlleleCallSet
/*****************************************************************************/
/**
 * GeneCall::getAlleleCallSet
 * Synopsis:
 *
 * Convienance function to return either the Reference or Variant allele set.
 *
 * @param index - between 0 and  num m_alleleSet.size()
 *
 * @return CallSet reference corresponding to the index.
 */
/*****************************************************************************/
CallSet &  GeneCall::getAlleleCallSet(unsigned int index)
{

  APT_ERR_ASSERT(m_alleleSet.size() > 0, "");

  return m_alleleSet[index];

}
// end GeneCall::getAlleleCallSet
/*****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//                      TRANSLATION ALGORITHM                               //
//////////////////////////////////////////////////////////////////////////////

/*****************************************************************************/
/**
 * _getProbeSetBaseCount
 * Synopsis:
 *
 * Helper function to be called only once by _calculateTHPR, broken out
 * for clarity.
 *
 * Any ProbeSet Id in the translation table will have two possible bases
 * for the first row, and one additional for every row thereafter.
 * This is represented in the completeSet passed in as one CallElement
 * per ProbeSet Id and a list of bases associated therein.
 *
 *
 *
 * @params chromatid1 - CallSet generated from ALLELE1 column from Genotype data.
 * @params completeSet - CallSet generated from the translation table
 *                       with possible multi-allelic bases.
 * @params ch1Index - CallElement index to compare. Note that CallElements
 *                   are keyed by ProbeSet. The i'th element refers to the
 *                   i'th element in chromatid1 and the same probe set in
 *                   completeSet that corresponds to the i'th element in
 *                   chromatid1.
 *
 * @returns count - the number of bases in completeSet.
 */
/*****************************************************************************/
static int _getProbeSetBaseCount(const CallSet & chromatid1,
                                 const CallSet & completeSet,
                                 const int ch1Index)
{

  APT_ERR_ASSERT(ch1Index < chromatid1.size(), "");

  std::map< std::string, CallElement>::const_iterator iCEit
  = chromatid1.m_ceSet.begin();

  for (int i = 0; i < ch1Index; i++, iCEit++) {};

  const CallElement & ceCompleteSet
  = completeSet.getProbeSetCallElement(iCEit->second.m_probeSet);

  APT_ERR_ASSERT(ceCompleteSet.getCompleteSetSize() > 1, "");

  return (ceCompleteSet.getCompleteSetSize());

}
// end GeneCall::_getProbeSetBaseCount
/*****************************************************************************/
/*****************************************************************************/
/**
 * _isHeteroZygousPair
 * Synopsis:
 *
 * Helper function to be called only once by _calculateTHPR, broken out
 * for clarity.
 *
 * If the CallElements between the two CallSets for the i'th element have
 * the same base then they are consider homozygous, else heterozygous.
 *
 *
 * @params chromatid1 - CallSet generated from ALLELE1 column from Genotype data.
 * @params chromatid2 - CallSet generated from ALLELE2 column from Genotype data.
 * @params ch1Index - CallElement index to compare. Note that CallElements
 *                   are keyed by ProbeSet. The i'th element refers to the
 *                   i'th element in chromatid1 and the same probe set in
 *                   chromatid2 that corresponds to the i'th element in
 *                   chromatid1.
 *
 * @returns true - if the bases for the CallElements compared are different.
 */
/*****************************************************************************/
static bool _isHeteroZygousPair(const CallSet & chromatid1,
                                const CallSet & chromatid2,
                                const int ch1Index)

{
  APT_ERR_ASSERT(ch1Index < chromatid1.size(), "");
  std::map< std::string, CallElement>::const_iterator iCEit =
    chromatid1.m_ceSet.begin();

  for (int i = 0; i < ch1Index; i++, iCEit++) {}

  const CallElement & ceChromatid1 = chromatid1.getProbeSetCallElement(iCEit->second.m_probeSet);
  const CallElement & ceChromatid2 = chromatid2.getProbeSetCallElement(iCEit->second.m_probeSet);

  // Ok, even though there can be multiple bases in a CallElement, that
  // only applies to CallElements constructed from the Translation Table.
  // For CallElements constructed from the Genotype data there is always
  // ever exactly one base.

  std::string base1 = ceChromatid1.m_bases[0];
  std::string base2 = ceChromatid2.m_bases[0];

  APT_ERR_ASSERT(!base1.empty() && !base2.empty(), "");


  return (base1 != base2);
}
// end GeneCall::_isHeteroZygousPair
/*****************************************************************************/
/*****************************************************************************/
/**
 * _calculateTHPR
 * Synopsis:
 *
 * Helper function to be called only once by GeneCall::translateExperimentCall,
 * broken out for clarity.
 *
 * Calculate the requred number of haplotype matches given tri-allelic and
 * wild-cards are in the mix. Note that MARKERs are considered HAPLOTYPE GROUP
 * of size 1 for the purposes of this API.
 *
 *
 * Key:
 *
 * 0 : false, no heterozygous base pairs have been detected.
 * 1 : true, at least one heterozygous base pair has been detected.
 *
 * THP0 : Total Haplotype Marker Pairs, no heterozygous base pairs have been detected.
 * THP1 : Total Haplotype Marker Pairs, at least one heterozygous base pair has been detected.
 * i : iteration count of the number of haplotype marker pairs analyzed so far.
 * Ni : for a Wildcard, iteration translation table possible allele bases
 *      corresponding to the
 *      probe set Haplotype marker pairs with the wildcard.
 *      For bi-allelic this is two bases,
 *      the reference and the variant. For tri-allelic this number is 3, for
 *      the reference and two variants.
 *
 *  * : wildcard pair of Haplotype Markers.
 *
 * Formuala:
 *
 *
 * THP0(0) = 1
 * THP1(0) = 0
 *
 *             { i+1: homozygous   => THP0(i)
 * THP0(i+1) = | i+1: heterozygous => 0
 *             { i+1: *            => THP0(i) * Ni
 *
 *             { i+1: homozygous   => THP1(i)
 * THP1(i+1) = | i+1: heterozygous => 2*THP1(i) + THP0(i)
 *             { i+1: *            => Ni * THP1(i) + ( Ni(Ni-1)/2 ) *
 *                                    2*THP1(i) + ( Ni(Ni - 1) )/2) * THP0(i)
 *
 * THPR = THP1(total) + THP0(total)
 *
 * @param chromatid1 - ALLELE1 column in the genotype file.
 * @param chromatid2 - ALLELE2 column in the genotype file.
 * @param completeSet - from the translation table, the complete set
 *                       of possible probe sets, possible multi-allelic.
 * @param HPPA      - passed in for return.
 *
 * @return HPPA - the number of haplotype pairs possible per allele.
 * @return THPR - the total number of required haplotype paris.
 *                 If this total is not found then UNK should be used
 *                 to generate the required total.
 */
/*****************************************************************************/
static int _calculateTHPR(const CallSet & chromatid1,
                          const CallSet & chromatid2,
                          const CallSet & completeSet,
                          int & HPPA)
{



  // NON-HAPLOTYPE GROUP MARKERS and HAPLOTYPE GROUPS of size 1.
  if (chromatid1.m_copyNumber == 0) {
    return 1;
  }

  // HAPLOTYPE GROUP MARKERS
  int THPR = 0;

  int THP0 = 1;
  int THP1 = 0;

  HPPA = 1;
  std::map< std::string, CallElement>::const_iterator iCEit = chromatid1.m_ceSet.begin();

  for (int i = 0; i < chromatid1.size(); i++, iCEit++) {

    bool isHet = _isHeteroZygousPair(chromatid1, chromatid2, i);
    int  Ni    = _getProbeSetBaseCount(chromatid1, completeSet, i);


    if (iCEit->second.hasWildCards()) {
      // HPPA
      HPPA *= Ni;
      // THP1
      THP1 = Ni * THP1 + (Ni) * (Ni - 1) * THP1 + ((Ni) * (Ni - 1) / 2) * THP0;
      // THP0
      THP0 = THP0 * Ni;
    } else if (isHet) {
      // THP1
      THP1 = 2 * THP1 + THP0;
      // THP0
      THP0 = 0;
    } else { // is homozygous
      // THP1 = THP1;
      // THP0 = THP0;
    }

  }

  THPR = THP0 + THP1;

  return THPR;

}
// end GeneCall::_calculateTHPR
/*****************************************************************************/
/*****************************************************************************/
/**
 * _translateExperimentCallCheck
 * Synopsis:
 *
 * Helper function to be called only once by GeneCall::translateExperimentCall,
 * broken out for clarity.
 *
 * Does all the upfront validation and exception throwing.
 * Note the multiple catches. This means all exceptions will be logged
 * if multiple exist and not just the first exception.
 *
 * @param rte - RunTimeEnvironment which contains options
 * @param completeSet - the known universe of accepted probeSets/bases.
 * @param chromatid1 - ALLELE1 from the Genotype Short Report.
 * @param chromatid2 - ALLELE2 from the Genotype Short Report.
 * @param alleleSet  - the CallSets of calls.
 *
 *
 * @returns true - if no exceptions were logged and its ok to continue
 *                 call translation.
 */
/*****************************************************************************/
static bool _translateExperimentCallCheck(RunTimeEnvironment & rte,
    const CallSet & completeSet,
    const std::string & geneName,
    const std::string & experimentName,
    CallSet & chromatid1,
    CallSet & chromatid2,
    const std::vector<CallSet> & alleleSet)
{

  bool okToCall = true;

  std::vector<class CallSet>::const_iterator testCallSetIt = alleleSet.begin();

  // All experiments should have a minimum of 1 CallSet, the reference set.
  APT_ERR_ASSERT(testCallSetIt != alleleSet.end(), "");

  if (!completeSet.contains(chromatid1)) {
    std::stringstream msgSStr;
    completeSet.describeVerbose(rte, ADT_VERBOSE_NORMAL);
    chromatid1.describeVerbose(rte, ADT_VERBOSE_NORMAL);
    msgSStr << experimentName << " " << geneName << " " << chromatid1.m_name <<  ": invalid probe set / base detected.";
    if (rte.m_adtOpts.m_ignoreUnknownAlleles) {
      okToCall = false;
    } else {
      APT_ERR_ABORT(msgSStr.str());
    }

  }

  if (!completeSet.contains(chromatid2)) {
    std::stringstream msgSStr;
    msgSStr << experimentName << " " << geneName << " " << chromatid2.m_name <<  ": invalid probe set / base detected.";
    if (rte.m_adtOpts.m_ignoreUnknownAlleles) {
      okToCall = false;
    } else {
      APT_ERR_ABORT(msgSStr.str());
    }
  }

  std::stringstream completeSStr;
  completeSStr <<  experimentName << ", " <<  geneName <<  ", " << testCallSetIt->m_name << ": incomplete probe set for haplotype group detected.";

  if (chromatid1.size() != testCallSetIt->size()) {

    if (rte.m_adtOpts.m_enforceCompleteHaplotypeGroup) {
      APT_ERR_ABORT(completeSStr.str());
    }
    Verbose::out(ADT_VERBOSE_EXCEPTION, completeSStr.str());
    okToCall = false;
  }

  if (chromatid2.m_ceSet.size() != testCallSetIt->m_ceSet.size()) {

    if (rte.m_adtOpts.m_enforceCompleteHaplotypeGroup) {
      APT_ERR_ABORT(completeSStr.str());
    }
    okToCall = false;
    Verbose::out(ADT_VERBOSE_EXCEPTION, completeSStr.str());
  }

  // Ok, this is a HAPLOTYPE GROUP or MARKER and all items are NC
  // then just return an empty call.
  if (chromatid1.hasWildCards() &&
      (chromatid1.numWildCards() == chromatid1.size())) {

    okToCall = false;
  }

  return okToCall;
}
// end _translateExperimentCallCheck
/*****************************************************************************/
/*****************************************************************************/
/**
 * _translateExperimentCallCopyNumberOne
 * Synopsis:
 *
 *  Helper function for translateExperimentCall broken out for brevity.
 *  Specifically translates to a copy number 1 call result.
 *
 * @param rte         - single instance of the run time environment
 * @param results     - for return
 * @param chromatid1  - from the experiment data for the single chromatid
 * @param geneCopyNumberZeroCallSet - for the specific gene
 *
 * @return results - updated with copy number "data" for pairing with
 * the chromatid1 passed in.
 */
/*****************************************************************************/
static void _translateExperimentCallCopyNumberOne(
  const RunTimeEnvironment & rte,
  CallResults & results,
  const CallSet & chromatid1,
  const CallSet & geneCopyNumberZeroCallSet)
{

  // COPY NUMBER 1
  if (chromatid1.m_copyNumber != 1) {
    return;
  }
  std::string copy1Name;

  if (geneCopyNumberZeroCallSet.getCallType() != ADT_CALL_TYPE_NULL) {
    copy1Name = geneCopyNumberZeroCallSet.m_name;
  } else {
    APT_ERR_ABORT(results.m_experimentName + " " + results.m_geneName + ": copy number 1 values detected with no copy number 0 designation in the translation file.");
  }

  if (results.size() == 0) {
    results.appendAlleleCall(copy1Name , 1, geneCopyNumberZeroCallSet);
    return;
  }

  std::map< std::string, bool> callCount;
  std::vector<CallResults::AlleleCall> newResults;


  for (int i = results.m_alleleCalls.size() - 1; i >= 0 ; i--) {

    std::string alleleCall = copy1Name +  "/" +  results.m_alleleCalls[i].m_allele;
    //cerr << experimentName << " " <<  geneName << " : copy number 0 call"  << alleleCall << endl;
    alleleCall = alleleCall.substr(0, alleleCall.rfind('/'));

    if (callCount[alleleCall]) {
      results.m_alleleCalls.erase(results.m_alleleCalls.begin() + 1);
      results.m_markerCallCount--;
    } else {
      callCount[alleleCall] = true;
      results.m_alleleCalls[i].m_allele = alleleCall;
    }
  }

  return;

}
// end _translateExperimentCallCopyNumberOne
/*****************************************************************************/
/*****************************************************************************/
/**
 * Translating an experiment call means checking for allele names that
 * correspond in the DMET3 translation table. For reporting purposes all
 * input experiment records from the Genotype CHIP or Short Report files
 * are to have a corresponding record somewhere in the reporting, i.e.
 * be included in the result set. There only exception are for unrecognized
 * genes. An unrecognized gene is simply skipped over. Records for unrecognized
 * genes will not be in the experiment results.
 *
 *
 * @param rte - RunTimeEnvironment which contains options
 * @param ttm - TranslationTableModel with the meta data about the call.
 * @param headerRow - the row in the ttm that has the meta data.
 *
 */
/*****************************************************************************/
CallResults GeneCall::translateExperimentCall(RunTimeEnvironment & rte,
    const CallSet & completeSet,
    GenotypeTableModel & gtm,
    CallSet & chromatid1,
    CallSet & chromatid2,
    const std::string & dmet2CopyNumberCall,
    const CallSet & geneCopyNumberZeroCallSet)
{


  bool isHaplotypeGroup
  = (chromatid1.getCallType() == ADT_CALL_TYPE_HAPLOTYPE_GROUP);

  unsigned int geneCopyNumber = 2;

  CallResults results(gtm.m_geneName, gtm.m_experimentName, chromatid1, chromatid2);


  // The callee passes in possibly invalid CallSets for the
  // 2 chromosomes derived from the
  // experiment records in the Genotype Short Report file.
  if (!  _translateExperimentCallCheck(rte, completeSet, gtm.m_geneName,
                                       gtm.m_experimentName, chromatid1,
                                       chromatid2, m_alleleSet)) {
    return results;
  }

  //DMET2 COPY NUMBER 0
  if (rte.m_adtOpts.m_dmet2Calling) {
    if (((rte.m_adtOpts.m_inputTTableType == ADT_TRANSLATION_TABLE_TYPE_DMET3) || isHaplotypeGroup) && !dmet2CopyNumberCall.empty()) {
      results.appendAlleleCall(dmet2CopyNumberCall, 0);

      Verbose::out(ADT_VERBOSE_CALL, gtm.m_geneName + ": " + gtm.m_experimentName + ": copy number call: " + dmet2CopyNumberCall);
      return results;
    }
  } else if (chromatid1.m_copyNumber == 0) {
    if (geneCopyNumberZeroCallSet.getCallType() != ADT_CALL_TYPE_NULL) {
      std::string alleleCall = geneCopyNumberZeroCallSet.m_name +  "/" +  geneCopyNumberZeroCallSet.m_name;
      results.appendAlleleCall(alleleCall, 0, geneCopyNumberZeroCallSet, geneCopyNumberZeroCallSet);
      return results;
    }

    else {
      geneCopyNumberZeroCallSet.describeVerbose(rte, ADT_VERBOSE_NORMAL);
      APT_ERR_ABORT(gtm.m_experimentName + " " + gtm.m_geneName + ": copy number 0 values detected with no copy number 0 designation in the translation file.");
    }
  }

  geneCopyNumber = chromatid1.m_copyNumber;

  // THPS = HPPA * AHPS - (AHPS * ( AHPS-1) )/ 2;
  int                HPPA  = 1; // Haplotype Pairs per allele
  int                THPS  = 0; // Total Haplotype pairs seen
  std::map< std::string, int > AHPS;      // Allele haplotype pairs seen.

  int                THPR  = 0; // Total Haplotype pairs required.
  // THPR is a complex calculation broken out into another routine.
  THPR = _calculateTHPR(chromatid1, chromatid2, completeSet, HPPA);

  for (std::vector<CallSet>::const_iterator iCSit = m_alleleSet.begin();
       iCSit != m_alleleSet.end(); iCSit++) {
    AHPS[iCSit->m_name] = 0;
  }


  if (rte.m_currentVerbosity == ADT_VERBOSE_CALL) {
    std::stringstream msgSStr;
    msgSStr << "Total Haplotype Pairs Possible/Haplotype Pairs Per Allele: ";
    msgSStr << THPR << "/" << HPPA << endl;
    msgSStr << "chromatid1: " << (*chromatid1.m_ceSet.begin()).first << ", ";
    msgSStr << (*chromatid1.m_ceSet.begin()).second.basesToString() << endl;
    msgSStr << "chromatid2: " << (*chromatid2.m_ceSet.begin()).first << ", ";
    msgSStr << (*chromatid2.m_ceSet.begin()).second.basesToString() << endl;
    Verbose::out(ADT_VERBOSE_CALL, msgSStr.str());
  }

  // The marker reference and variant sets are always the first two
  // sets, skip over them for haplotype groups;
  int offset = isHaplotypeGroup ? 2 : 0;

  for (std::vector<CallSet>::const_iterator iCSit = m_alleleSet.begin() + offset;
       iCSit != m_alleleSet.end() && (THPS < THPR);  ++iCSit) {

    // The following condition is possible if the allele name was
    // added as the second half of the call. If the number possible pairs
    // has already been seen for this allele then just optimize
    // by continuing.
    if (!(AHPS[iCSit->m_name] < HPPA))  continue;

    const CallSet & firstSet = *iCSit;

    Verbose::out(ADT_VERBOSE_CALL, "First CallSet: " + firstSet.m_name, false);

    // Check if any combination of the two chromsomes match each
    // element. If so we have a match and the second match
    // will be a result of the remaining bases not used to make the first
    // match.

    if (! firstSet.match(chromatid1, chromatid2)) {
      Verbose::out(ADT_VERBOSE_CALL, " not found!");
      continue;
    }

    Verbose::out(ADT_VERBOSE_CALL, " match!");

    // Build the second set from the difference between the
    // the test set and the inputs.

    CallSet secondSet = firstSet.buildSecondSet(rte, chromatid1, chromatid2);

    for (std::vector<CallSet>::const_iterator jCSit = iCSit;
         (jCSit != m_alleleSet.end()) && (results.size() < THPR) &&
         (AHPS[firstSet.m_name] < HPPA);   ++jCSit) {

      Verbose::out(ADT_VERBOSE_CALL, "Second CallSet: " + jCSit->m_name, false);

      if (secondSet == *jCSit) {

        if (jCSit->m_name != firstSet.m_name) {
          AHPS[jCSit->m_name]++;
        }
        AHPS[firstSet.m_name]++;

        std::string alleleCall = firstSet.m_name +  "/" +  jCSit->m_name;
        results.appendAlleleCall(alleleCall, geneCopyNumber, firstSet, *jCSit);

        Verbose::out(ADT_VERBOSE_CALL, " match (" + alleleCall + ")!");
      } else {
        Verbose::out(ADT_VERBOSE_CALL, " not found!");
      }

    } // for each allele in the translation table get a second match.

    if (AHPS[firstSet.m_name] < HPPA) {
      std::string alleleCall = firstSet.m_name + "/UNK" ;
      results.appendAlleleCall(alleleCall, geneCopyNumber, firstSet);
    }
    THPS += HPPA;

  } // For each allele in the translation table, check for a match.

  // If the number of allele pairs seen is less than the total required
  // then we add "UNK/UNK" to represent the difference.
  if (THPS < THPR) {
    results.appendAlleleCall("UNK/UNK", geneCopyNumber);
  }

  _translateExperimentCallCopyNumberOne(rte, results, chromatid1, geneCopyNumberZeroCallSet);

  results.orderCalls();

  return results;

}
// end GeneCall::translateExperimentCall
/*****************************************************************************/
