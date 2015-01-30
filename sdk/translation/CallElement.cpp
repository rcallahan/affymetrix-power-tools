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
 * @file   CallElement.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:02:01 PDT 2008
 * @brief  Atomic class element for translation.
 *
 */

#include "translation/CallElement.h"
//
#include "util/Err.h"        // includes Verbose.h
//
#include "pcrecpp.h"
//

using namespace std;

/*****************************************************************************/
/**
 * CallElement::CallElement()
 * Synopsis:
 *
 * Empty constructor for STL storage containers or for
 * translation file call elements.
 *
 * NOTE: base is used interchangabley to mean allele. In DMET 2
 * alleles where always represented as just one base.
 *
 *
 */
/*****************************************************************************/
CallElement::CallElement() :  m_chromatid(0), m_probeSet("")
{
  m_hasWildCards         = false;
  m_hasReverseComplement = false;
  m_completeSetSize      = -1;
  
}
// end CallElement::CallElement
/*****************************************************************************/

/*****************************************************************************/
/**
 * CallElement::CallElement()
 * Synopsis:
 * GenotypeTableModel (CHP) version of the CallElement constructor.
 *
 * This is not obvious but when a CallElement is created for
 * the GenotypeTableModel, it will always have exactly one base.
 * Multiple bases only apply to mult-allelic defintions in the
 * Translation table.
 *
 * NOTE: base is used interchangabley to mean allele. In DMET two
 * alleles where always represented as just one base.
 *
 * @param probeSet    - aka "marker", ie.e AM_14633, etc.
 * @param base        - aka "alleles", A, or ATATGGA. For multiallelic
 *                      markers intialize with one base and then
 *                      append markers using CallElement::addBase().
 * @param copyNumber  - 2 by default. 0 and 1 are the only other legal values.
 * @param hasRevComp  - reverse complement, if the "Switch" field in translation
 *                      file or annotation file is "Y/1" then the reverse
 *                      compliment of the base will be used durring translation.
 *                      See "CallElement::reverseComplement".
 *
 */
/*****************************************************************************/
CallElement::CallElement(const std::string & probeSet,  const std::string & base, const unsigned int copyNumber, const bool hasRevComp):    m_chromatid(0), m_copyNumber( copyNumber),  m_hasReverseComplement(hasRevComp),   m_probeSet(probeSet)
{

  m_completeSetSize      = -1;

  APT_ERR_ASSERT(m_copyNumber <= 2, "Invalid parameters to CallElement constructor");

  if (! isValidProbeSet(probeSet)) {
    APT_ERR_ABORT(probeSet + ": Invalid probeSet not recognized by translation file");
  }

  if (! isValidBase(base)) {
    APT_ERR_ABORT(probeSet + ": Base ("  + base + ") not recognized by translation file.");
  }

  addBase(base);

  m_hasWildCards = isWildCardBase(base);

}
// end CallElement::CallElement
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::addBase
 * Synopsis:
 *
 * Append a base (allele) to the std::vector of bases (alleles). Insure the
 * added base is not a duplicate as well that the base is valid.
 * Also adds the passed in base to its reverse complement list if the
 * element the 'hasReverseComplement' attribute is true;
 *
 * @param base - the base to add.
 * @return true - if base was added successfully.
 *
 */
/*****************************************************************************/
bool CallElement::addBase(const std::string & base)
{

  if (!isValidBase(base)) {
    return false;
  }

  for (int i = 0 ; i < m_bases.size(); i++) {
    if (m_bases[i] == base) {
      return false;
    }
  }

  m_bases.push_back(base);

  if (m_hasReverseComplement) {
    m_rbases.push_back(reverseComplement(base));
  }

  return true;

}
// end CallElement::addBase
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::basesToString
 * Synopsis:
 *
 * Return a list of the bases (alleles) for this CallElement as a
 * comma separated list. For most bi-allelic markers in the translation
 * table the list size will be one.
 *
 * @param reverseComplement - return the list of reverse complements and not
 * the original bases (alleles).
 *
 * @returns basesString - i.e. "A,T,G".
 *
 */
/*****************************************************************************/
std::string CallElement::basesToString(bool reverseComplement) const
{

  std::string basesString;

  if (reverseComplement) {
    for (int i = 0; i < m_rbases.size(); i++) {
      if (i == 0) {
        basesString = m_rbases[0];
      } else {
        basesString = basesString + "," + m_rbases[i];
      }
    }

  } else {
    for (int i = 0; i < m_bases.size(); i++) {
      if (i == 0) {
        basesString = m_bases[0];
      } else {
        basesString = basesString + "," + m_bases[i];
      }
    }
  }

  return basesString;

}
// end basesToString
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::setCompleteSetSize
 * Synopsis:
 *
 * A CallElement typically represents one marker, but in one special
 * case the complete call element set represents all alleles (bases)
 * possible for a marker. This API is used in that one special
 * case to get the complete set size. The complete set size
 * is used in the translation algorithm.
 * 
 * Complete set size, n: 1 The count of all allels (bases) minus any copy
 * number bases. 
 *
 * return size - the count
 *
 */
/*****************************************************************************/
int CallElement::setCompleteSetSize() {


  m_completeSetSize = 0;
  
  for ( int i = 0; i < m_bases.size(); i++ ) {
    if (m_bases[i] != "0") {
      m_completeSetSize++;
    }
  }

  return m_completeSetSize;
  
}
// end CallElement::setCompleteSetSize()
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::describeVerbose
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
void CallElement::describeVerbose(const RunTimeEnvironment & rte,
                                  ADT_VERBOSE_ENUM override) const
{

  if (!override && (rte.m_currentVerbosity < ADT_VERBOSE_INPUT_FILES))
    return;

  if (override) {
    describeVerbose(override);
  } else {
    describeVerbose(ADT_VERBOSE_INPUT_FILES);
  }

  return;

}
/*****************************************************************************/
void CallElement::describeVerbose(ADT_VERBOSE_ENUM level) const
{


  std::stringstream msgSStr;

  msgSStr << "CallElement::describeVerbose (probeSet : bases): (";
  msgSStr << m_probeSet << ":" << basesToString();

  if (m_hasReverseComplement) {
    msgSStr << " [ reverse complement bases : ";
    msgSStr << basesToString(true) << " ]";
  }

  msgSStr << ")";

  Verbose::out(level, msgSStr.str());

  return;

}
// end CallElement::describeVerbose
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::isValidProbeSet
 * Synopsis:
 *
 * Regular expression test to examine if a probe set (marker) std::string is
 * valid.
 *
 *
 * @param testProbeSet - the std::string to validate
 *
 * @return true - if valid
 */
/*****************************************************************************/
bool CallElement::isValidProbeSet(const std::string & testProbeSet)
{

  pcrecpp::RE re("^[_\\-\\w\\d]*\\d+$");

  return re.FullMatch(testProbeSet);

}
// end CallElement::isValidProbeSet
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::isValidBase
 * Synopsis:
 *
 * Regular expression test to examine if a base (allele) is valid.
 * Wildcard designations such as "PossibleRareAllele" and "NoCall"
 * are handled correctly and will return true.
 *
 *
 * @param testBase - the base (allele) std::string to validate.
 *
 * @return - true if valid
 *
 */
/*****************************************************************************/
bool CallElement::isValidBase(const std::string & testBase)
{


  pcrecpp::RE re("[ACTG0]+|-|INS|DEL|Ins|Del");

  return (isWildCardBase(testBase) || re.FullMatch(testBase));

}
// end CallElement::isValidBase
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::isWildCardBase
 * Synopsis:
 *
 * Regular exxpression test to examine if a base (allele) is a wild card.
 *
 * Wildcard list:
 * 1.) NoCall (NC)
 * 2.) PossibleRareAllele (PRA)
 * 3.) NotAvailable
 *
 *
 * @param base - the base (allele) to test.
 *
 * @return true - if the base (allele) is a wild card.
 *
 */
/*****************************************************************************/
bool CallElement::isWildCardBase(std::string base)
{


  if (pcrecpp::RE("NC|NoCall|PRA|PossibleRareAllele|NotAvailable").FullMatch(base)) {
    return(true);
  }

  return false;

}
// end CallElement::isWildCardBase
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::operator==
 * Synopsis:
 *
 * The purpose of the CallElement class is to provide this one operator.
 *
 * This operator treats the probe set (marker) id and a corresponding
 * list of bases (alleles) as an atomic unit.
 *
 * Comparing call elements as equivalent means the atomic consituents
 * must also be equivalent according to a set of rules. 
 *
 * Business rules treat no call and  possible rare
 * allele as wildcards that match any base. It also compares
 * multi-allelic CallElements and if they share a common base then
 * they are equal. Note that for CallElements derived from
 * the Genotype Short Report data they can never be multi-allelic.
 * This means that the only multi-allelic CallElement being
 * compared will be the Translation Table file CallElement with the
 * GenotypeTableModel (CHP) CallElement.
 *
 * @param ce1 - the CallElement to compare
 *
 * @return true - if equal.
 */
/*****************************************************************************/
bool CallElement::operator==(CallElement const & ce1) const
{


  if (m_probeSet != ce1.m_probeSet)  return false;

  for (int i = 0; i < m_bases.size(); i++) {

    if (isWildCardBase(m_bases[i]))     return true;

    for (int j = 0; j < ce1.m_bases.size(); j++) {

      if (isWildCardBase(ce1.m_bases[j])) return true;

      if (m_bases[i] == ce1.m_bases[j])  {
        return true;
      }

    }

  }

  return false;

}
// end CallElement::operator==
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallElement::reverseComplement
 *
 * Synopsis:
 *
 *   Switch to the other strand by taking the reverse complement.
 *
 * @param seq - the base (allele) sequence to take the complement of.
 *
 * @return rev - the reverse complement
 */
/*****************************************************************************/
std::string CallElement::reverseComplement(const std::string seq)
{

  std::string rev;

  rev.clear();

  if (CallElement::isWildCardBase(seq)) {
    return seq;
  }

  for (int i = seq.size() - 1; i >= 0; i--) {
    switch (seq[i]) {
    case '-':
      rev.append("-");
      break;
    case 'a':
      rev.append("t");
      break;
    case 'c':
      rev.append("g");
      break;
    case 'g':
      rev.append("c");
      break;
    case 't':
      rev.append("a");
      break;
    case 'A':
      rev.append("T");
      break;
    case 'C':
      rev.append("G");
      break;
    case 'G':
      rev.append("C");
      break;
    case 'T':
      rev.append("A");
      break;
    case '0':   // Copy Number 0
      rev.append("0");
      break;
    default:
      std::string temp(seq[i], 1);
      rev.append(temp);
    }
  }
  return rev;

}
// end CallElement::reverseComplement
/*****************************************************************************/
