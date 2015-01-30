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
 * @file   CallSet.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:00:56 PDT 2008
 * @brief  A little more functionaliy than using the STL 'set' container for CallElement objects.
 */

#include "translation/CallSet.h"
//
#include "translation/GenotypeTableModel.h"
#include "translation/TranslationTableModel.h"
//
#include "calvin_files/array/src/ArrayId.h"
#include "util/Err.h" // includes Verbose.h
//

using namespace std;

/*****************************************************************************/
/** 
 * CallSet::CallSet
 * Synopsis
 *
 * 1.) Default constructor for STL containers.
 * 2.) Translation table constructor used in "TranslationTable.cpp".
 * 
 * @param name - the allele or marker name "*1". 
 * @param columnInTranslationTable - A1-ANN index, "Ref" or "Var" index
 * for markers. 
 * 
 */
/*****************************************************************************/
CallSet::CallSet(const std::string & name, const int & columnInTranslationTable ) :  m_columnInTranslationTable(columnInTranslationTable), m_name(name)
{
  m_isDescriptive        = false;
  m_isMultiAllelic       = false;
  m_hasWildCards         = false;
  m_type                 = ADT_CALL_TYPE_NULL;
  m_copyNumber           = 2;
  m_hasZeroCopyNumberNotAvailable = false;
}
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::CallSet
 * Synopsis:
 * 
 * Constructor for CHP experiment data to be compared with translation table
 * call sets instantiated with the other constructor. 
 * A call set represents one of two possible chromatids.
 *
 * The GenoTypeTableModel experiment data from the CHP file is
 * passed in preconfigured to a projection of sorted markers for
 * one particular gene. Further the reference set represents
 * all the probesets (markers) to be included for this gene call set.
 *   For a SNP, all the gene markers are scanned but only one is used.
 * If a haplotype group for, say, "*7", requires six markers then
 * all gene markers are scanned and only the six are used. 
 *
 * This constructor could be improved by building a probe set index (std::map)
 * in the GenotypeTableModel. 
 *
 * 
 *
 * @param rte          - The single instance RunTimeEnvironment
 * @param gtm          - The single instance experiment model 
 * @param referenceSet - The reference set of accepted probeSets for this
 *                      experiment gene (from the Translation Table).
 *
 * @param chromoTid - 1 or 2.
 */
/*****************************************************************************/
CallSet::CallSet(const class RunTimeEnvironment & rte,
                 class GenotypeTableModel & gtm,
                 const CallSet & referenceSet, const int chromaTid)
{

  APT_ERR_ASSERT((chromaTid == 1) || (chromaTid == 2), "");

  m_isDescriptive        = true;
  m_isMultiAllelic       = false;
  m_hasWildCards         = false;
  m_type                 = referenceSet.getCallType();
  m_copyNumber           = 2;
  m_hasZeroCopyNumberNotAvailable = false;

  std::stringstream nameSStr;

  nameSStr << "chromatid" << chromaTid;

  m_name = nameSStr.str();


  for (int i  = 0; i < gtm.getGeneRowSize(gtm.m_geneName); i++) {

    int row = gtm.getGeneRow(gtm.m_geneName, i);

    // IGNORE previously invalidated experiment rows due to bad probe sets or bases.

    std::string exprProbeSet = gtm.m_rows[row][GT_PROBE_SET_INDEX];


    // FILTER, only add the element if the GeneCall this CallSet is for
    // includes the CallElement in its CompleteSet.

    if (referenceSet.hasProbeSet(exprProbeSet)) {

      std::string allele1 = gtm.m_rows[row][GT_ALLELE1_INDEX] ;

      bool isWildCardBase = CallElement::isWildCardBase(allele1);

      if (!isWildCardBase) {
        _setCopyNumber(rte, gtm, row);
      }

      if ( (m_copyNumber == 0 ) &&
           pcrecpp::RE("NA|NotAvailable").FullMatch(allele1) ) {
        m_hasZeroCopyNumberNotAvailable = true;
      }

      std::string exprBase;

      if (isWildCardBase || (chromaTid == 1) || (m_copyNumber < 2)) {   // ALLELE1 column
        exprBase = gtm.m_rows[row][GT_ALLELE1_INDEX];

        APT_ERR_ASSERT(!exprBase.empty(), "");

        // Set the base
      } else {
        // ALLELE2 column
        exprBase          = gtm.m_rows[row][GT_ALLELE2_INDEX];
      }

      m_ceSet[exprProbeSet]             = CallElement(exprProbeSet, exprBase);

      if (m_ceSet[exprProbeSet].hasWildCards()) {
        m_hasWildCards = true;
      }

    }
  }


}
// end CallSet::CallSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::addCallElement
 * Synopsis:
 * 
 * Given a row from the Translation Table, adds a CallElement
 * ProbeSet/base pair to the CallSet. The element added
 * will correspond to the passed in name (take from columns A1-ANN)
 * unless reference is specified.
 *
 * Returns false if the set element being added is a duplicate.
 *
 * @param ttm - TranslationTableModel raw data that mirors the input TsvFile.
 * @param row - the row in question.
 * @param base_column - the Allele Name column or Reference column.
 * @param allowMultiAllelic - if the element already exists treat the
 *                            duplicate ProbeSet as multiAllelic.
 *
 * @returns true/false - true if successful
 */
/*****************************************************************************/
bool CallSet::addCallElement(TranslationTableModel & ttm,
                             const int row,
                             const int base_column,
                             const bool allowMultiAllelic,
                             const bool isHaplotype)
{


  APT_ERR_ASSERT((row >= 0) && (row < ttm.size()), "");
  APT_ERR_ASSERT((base_column >= 0) && (base_column < ttm.columnCount()), "");

  bool okToAdd = true;

  std::string probeSet = ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_PROBE_SET_ID)];

  std::string base = ttm.m_rows[row][base_column];

  bool hasReverseComplement = ((ttm.getType() == ADT_TRANSLATION_TABLE_TYPE_DMET3) && (ttm.m_rows[row][ttm.getColumnIndex(ADT_DMET3_TT_SWITCH_STRAND)] == "Y"));

  if (isHaplotype) {
    m_type = ADT_CALL_TYPE_HAPLOTYPE_GROUP;
  }

  std::map<std::string, CallElement>::const_iterator itFind;

  CallElement testCE(probeSet, base, 2, hasReverseComplement);

  if ((itFind = m_ceSet.find(probeSet)) == m_ceSet.end()) {

    m_ceSet[probeSet] = testCE;

    if (CallElement::isWildCardBase(base)) {
      m_hasWildCards = true;
    }

  } else if (allowMultiAllelic) {
    if (!contains(testCE)) {
      okToAdd = m_ceSet[probeSet].addBase(base);
    }
  } else {
    okToAdd = false;
  }

  if (okToAdd) {
    _setCopyNumber(base);
  }

  return okToAdd;

}
// end CallSet::addCallElement
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::buildSecondSet
 * Synopsis:
 *
 * A method used for translation. The idea is that during translation
 * alleles in the two chromatid CallSets can be mixed and matched.
 * One allele from chromatid 1, the rest from chromatid2 for example.
 * This combination is instantiated as a CallSet. This method is
 * invoked on the combination CallSet and then the difference
 * between the elments in this set and the initial chromatid1 and
 * chromatid2 sets is return.
 * 
 * In other words, given a matched allele set from any combination
 * of probeSet/bases from the two chromatid experiment sets for a gene, then
 * constructs difference "second" set using the components of both
 * not found in the matched set.
 *
 *
 * @param rte - the single instance of the run time enviroment
 * @param chromatid1 - from the Genotype experiment data set.
 * @param chromatid2 - from the Genotype experiment data set.
 */
/*****************************************************************************/
CallSet CallSet::buildSecondSet(const RunTimeEnvironment & rte,
                                const CallSet & chromatid1,
                                const CallSet & chromatid2) const
{

  CallSet secondSet("Difference");

  // Simple case, no wild cards.

  secondSet.m_isDescriptive = true;

  std::map<std::string, CallElement>::const_iterator iCEit;
  std::map<std::string, CallElement>::const_iterator jCEit;
  std::map<std::string, CallElement>::const_iterator kCEit;

  for (iCEit = m_ceSet.begin(); (iCEit != m_ceSet.end()); iCEit++) {

    jCEit = chromatid1.m_ceSet.find(iCEit->first);
    kCEit = chromatid2.m_ceSet.find(iCEit->first);
    // Ok, all sets are the same size and contain the same ProbeSets in
    // precisely the same sorted order. In other words, each set should
    // have exactly the same ProbeSets. We assert this.

    APT_ERR_ASSERT((jCEit != chromatid1.m_ceSet.end()) &&
                    (kCEit != chromatid2.m_ceSet.end()), "");


    // Ok, if one is a wild card then both must be a wild card.
    // In other words, wild card calls are always homozygous.
    if (jCEit->second.hasWildCards()) {
      APT_ERR_ASSERT(kCEit->second.hasWildCards(), "");
      secondSet.m_ceSet[iCEit->first] = jCEit->second;
    } else if (iCEit->second == jCEit->second) {
      CallElement newCE = kCEit->second;
      newCE.m_chromatid = 2;
      secondSet.m_ceSet[iCEit->first] = newCE;
    } else if (iCEit->second == kCEit->second) {
      CallElement newCE = jCEit->second;
      newCE.m_chromatid = 1;
      secondSet.m_ceSet[kCEit->first] = newCE;
    } else {
      APT_ERR_ASSERT(false, "");   // Progammer error, never reached condition.
    }

  }
  APT_ERR_ASSERT(secondSet.size() == this->size(), "");

  return secondSet;

}
// CallSet::buildSecondSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::clear
 * 
 * Re-iniatialize the CallSet to default conditions.
 *
 */
/*****************************************************************************/
void CallSet::clear()
{

  m_isDescriptive        = false;
  m_isMultiAllelic       = false;
  m_hasWildCards         = false;
  m_ceSet.clear();


}
// end CallSet::clear
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::contains( CallSet )
 * Synopsis
 *
 * A set operator. Detemine if all the CallElements of a CallSet are
 * all contained
 * within the one being interogated. Calls the "contains" operator
 * of each CallElement. See CallSet::contains( CallElement). 
 *
 *
 * @param c - the presumed subset.
 *
 * @returns true - if c is contained in this CallSet.
 */
/*****************************************************************************/
bool CallSet::contains(const CallSet & c) const
{

  bool contained = true;

  std::map<std::string, CallElement>::const_iterator iCEit;

  for (iCEit = c.m_ceSet.begin();
       contained && (iCEit != c.m_ceSet.end()); iCEit++) {

    contained = contains(iCEit->second);

  }

  return contained;

}
// end CallSet::contains
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::contains( CallElement )
 * Synopsis:
 *
 * If a CallElement has the same name and the CallElement::operator==
 * matches one in this CallSet then the CallElement is said to be
 * contained. 
 *
 * Primarily used in CallSet::contains( CallSet )
 *
 * @param ce - the CallElement to compare with this CallSet's own list. 
 *
 * @returns true - if contained
 */
/*****************************************************************************/
bool CallSet::contains(const CallElement & ce) const
{

  bool contained = false;

  std::map<std::string, CallElement>::const_iterator iCEit;

  for (iCEit = m_ceSet.begin();
       !contained && (iCEit != m_ceSet.end()); iCEit++) {

    if (iCEit->second == ce) {
      contained = true;
    }

  }

  return contained;

}
// end CallSet::contains
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::describeVerbose
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
void CallSet::describeVerbose(const RunTimeEnvironment & rte,
                              ADT_VERBOSE_ENUM override) const
{
  if (!override && (rte.m_currentVerbosity < ADT_VERBOSE_INPUT_FILES))
    return;

  if (!override) {
    describeVerbose(rte.m_currentVerbosity) ;
  } else {
    describeVerbose(override);
  }
}
// end CallSet::describeVerbose
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::describeVerbose
 *     
 * Run time debug routine that can be invoked with various verbose levels.
 * ADT_VERBOSE_INPUT_FILES is the default level.
 *
 * @param level    - ADT_VERBOSE_INPUT_FILES by default
 *                    otherwise the ADT_VERBOSE_ENUM level to output
 *
 */
/*****************************************************************************/
/*****************************************************************************/
void CallSet::describeVerbose(ADT_VERBOSE_ENUM level) const
{


  std::stringstream msgSStr;

  msgSStr << "CallSet::describeVerbose name: " << m_name << ": ";

  if (m_ceSet.size() == 0) {
    msgSStr << " EMPTY ";
  } else {
    msgSStr << "isDescriptive : " << m_isDescriptive << endl;

    msgSStr << "CallSet::describeVerbose name: " << m_name << ": ";
    msgSStr << "CallElement Set size: " << m_ceSet.size() << endl;
    msgSStr << "CallSet::describeVerbose name: " << m_name << ": ";
    msgSStr << "type: ";
    std::string type = "ADT_CALL_TYPE_NULL";
    if (m_type == ADT_CALL_TYPE_MARKER) {
      type = "ADT_CALL_TYPE_MARKER" ;
    } else if (m_type == ADT_CALL_TYPE_HAPLOTYPE_GROUP) {
      type = "ADT_CALL_TYPE_HAPLOTYPE_GROUP";
    }
    msgSStr << type;
  }

  Verbose::out(level, msgSStr.str());

  std::map<std::string, CallElement>::const_iterator itCE;

  for (itCE = m_ceSet.begin();
       itCE != m_ceSet.end();
       ++itCE) {

    itCE->second.describeVerbose(level);

  }


  return;

}
// end CallSet::describeVerbose
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::getProbeSetCallElement
 * Synopsis:
 *
 * A selector to detemine if the set of CallElements contains a
 * particular ProbeSet and if so then return a reference to the
 * matching CallElement. 
 *
 * @param probeSet - the probeSet (marker) in question
 *
 * @return true - if the ProbeSet std::string passed in is in the m_ceSet.
 */
/*****************************************************************************/
const CallElement & CallSet::getProbeSetCallElement(const std::string & probeSet) const
{


  APT_ERR_ASSERT(hasProbeSet(probeSet), "");

  std::map<std::string, CallElement >::const_iterator iCEit = m_ceSet.begin();

  for (; iCEit != m_ceSet.end(); iCEit++) {

    if (iCEit->second.m_probeSet == probeSet) {
      break;
    }
  }

  if (iCEit == m_ceSet.end()) {
    APT_ERR_ASSERT(false, "");
  }

  return iCEit->second;

}
// end CallSet::getProbeSetCallElement
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::getMarkerProbeSet
 * Synopsis
 *
 * A convienance method that takes advantage of the fact that a
 * marker will have only one item in the set and hence just
 * returns the first CallElement's name (probeSet id).
 * 
 * It is an exception to call this API on a non-ADT_CALL_TYPE_MARKER
 * CallSet. 
 * 
 *
 * @returns probeSet - the probeSet  name of the first element in the ceSet.
 */
/*****************************************************************************/
std::string CallSet::getMarkerProbeSet() const
{

  APT_ERR_ASSERT(getCallType() == ADT_CALL_TYPE_MARKER, "");

  std::map<std::string, CallElement>::const_iterator iit;

  iit = m_ceSet.begin();

  return iit->first;

}
// end CallSet::getMarkerProbeSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::hasProbeSet
 * Synopsis:
 *
 *  Detemine if the set of CallElements contains a particular ProbeSet.
 *
 * @param probeSet - the probeSet in question
 * 
 * @returns true - if the ProbeSet std::string passed in is in the m_ceSet.
 */
/*****************************************************************************/
bool CallSet::hasProbeSet(const std::string & probeSet) const
{


  bool hasProbeSet = false;

  std::map<std::string, CallElement>::const_iterator iCEit = m_ceSet.find(probeSet);

  hasProbeSet = (iCEit != m_ceSet.end());

  return hasProbeSet;

}
// end CallSet::hasProbeSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::match
 * Synopsis:
 *
 * The heart of the translation process deep in the guts of
 * GeneCall::translate. Detemine if the CallSet can match any combination
 * of CallElements given a pair of chromatid CallSets.
 *
 * If a translation table CallSet matches some combination of chromatid1
 * and chromatid2 experiment data CallSet then a match is made. Then
 * a call to CallSet::buildSecondSet is made and the diffenerce
 * is matched to any of the remaining translation table CallSets for
 * for this gene. If not than an "UNK" determination is made. 
 *
 * @param a - chromatid1 from the CHP 
 * @param b - chromatid2 from the CHP
 *
 * @returns true - if a combination can be made.
 */
/*****************************************************************************/
bool CallSet::match(const CallSet & a, const CallSet & b) const
{

  bool match = true;

  std::map<std::string, CallElement>::const_iterator iCEit;
  std::map<std::string, CallElement>::const_iterator c1CEit;
  std::map<std::string, CallElement>::const_iterator c2CEit;

  for (iCEit = m_ceSet.begin(); (iCEit != m_ceSet.end()) && match; iCEit++) {

    c1CEit = a.m_ceSet.find(iCEit->first);
    c2CEit = b.m_ceSet.find(iCEit->first);

    APT_ERR_ASSERT((c1CEit != a.m_ceSet.end()) && (c2CEit != b.m_ceSet.end()), "");

    if ((iCEit->second != c1CEit->second) &&
        (iCEit->second != c2CEit->second)) {
      match = false;
    }
  }

  return match;

}
// end CallSet::match
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::numWildCars
 * Synopsis:
 *
 * Count the number of wildcard CallElements within the CallSet.
 *
 * @return numWildCards
 *
 */
/*****************************************************************************/
int CallSet::numWildCards() const
{

  int numWildCards = 0;

  std::map<std::string, CallElement>::const_iterator iCEit;

  for (iCEit = m_ceSet.begin(); (iCEit != m_ceSet.end()) ; iCEit++) {

    if (iCEit->second.hasWildCards()) {
      numWildCards++;
    }

  }

  return numWildCards;

}
// end CallSet::numWildCards
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::operator==
 *
 * Essentially just walk the list of CallElements and comparing all for
 * equality, with some added business rules for Copy Number.
 *
 *
 * @param c - The CallSet to compare.
 * @returns true - if the CallSets are "equal".
 *
 */
/*****************************************************************************/
bool CallSet::operator==(const CallSet & c) const
{

  // CallSets can only be equal if they are valid.
  if ((! m_isDescriptive) || (! c.m_isDescriptive)) {
    return false;
  }

  if (m_ceSet.size() != c.m_ceSet.size()) {
    return false;
  }
  // Use the operator== for each CellElement to
  // detect if there is a difference.

  std::map<std::string, CallElement>::const_iterator it1;
  std::map<std::string, CallElement>::const_iterator it2;

  for (it1 = m_ceSet.begin();  it1 != m_ceSet.end();   it1++, it2++) {

    it2 = c.m_ceSet.find(it1->first);

    if (it2 == c.m_ceSet.end()) {
      return false;
    }

    if ((it1->second != it2->second)) {
      return false;
    }
  }

  return true;
}
// CallSet::operator==
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallSet::_setCopyNumber:
 * Synopsis:
 *
 *   Sets the call number to 0, 1 or 2 for the entire gene based on
 * experiment data. One probeset's copy number sets the entire
 * CallSet because copy number is
 * the same for the the entire gene set of markers.
 * 
 * Experiment CHP data should never contain varying copy number
 * types. However, hand crafted data may contain varying copy number
 * states. This routine will set the copy state to the lowest
 * number possible (0) when called repeatedly. Once copy number 0
 * is established then that's it. All markers will then be treated
 * as copy nuber 0 regardless of actual data values. 
 *
 *
 * @param rte - The RunTimeEnvironment
 * @param gtm - the experiment records for one gene.
 * @param row - the row in the gtm
 *
 * @return - void, sets member value "m_copyNumber".
 */
/*****************************************************************************/
void CallSet::_setCopyNumber(const RunTimeEnvironment & rte,
                             GenotypeTableModel & gtm,
                             const int row)
{


  // If any previous determination of copy number 0 was reached
  // then it applies to an entire haplotpye group.


  if ( gtm.m_geneCopyNumber.count(gtm.m_geneName) > 0 ) {
    m_copyNumber = gtm.m_geneCopyNumber[gtm.m_geneName];
    return;
  }
  
  if (m_copyNumber == 0) {
    return;
  }

  std::string allele1Base = gtm.m_rows[row][GT_ALLELE1_INDEX];
  std::string allele2Base = gtm.m_rows[row][GT_ALLELE2_INDEX];

  if (allele1Base.empty() || CallElement::isWildCardBase(allele1Base)) {
    return;
  }

  if (allele1Base == "0") {
    m_copyNumber = 0;
  } else if (allele2Base.empty()) {
    m_copyNumber = 1;
  }

  return;
}
// end CallSet::_setCopyNumber
/*****************************************************************************/
/*****************************************************************************/
/**
 * _setCopyNumber:
 * Synopsis:
 *
 *   Sets the call number to 0 if base is "0";
 *
 *
 * @param base - the allele "base", could be "0".
 *
 * @return - void, sets member value "m_copyNumber".
 */
/*****************************************************************************/
void CallSet::_setCopyNumber(const std::string & base)
{


  if (base == "0") {
    m_copyNumber = 0;
  }

  return;
}
// end CallSet::_setCopyNumber
/*****************************************************************************/
