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
 * @file   CallResults.cpp
 * @author Mybrid Spalding
 * @date   Wed Apr 30 07:52:49 PDT 2008
 * @brief  Binary representation of allele translation call.
 */

#include "translation/CallResults.h"
//
#include "util/Err.h" // includes "util/Verbose.h"
//
#include "pcrecpp.h"
//
#include <cassert>
#include <cstring>
#include <string>
//

using namespace std;

/*****************************************************************************/
/**
 * CallResults::CallResults
 * Synopsis:
 *
 * The process is to instantiate the object with the experiment data
 * prior to translation, pass this object to translation which then updates
 * the object with the translation results.
 *
 *
 *
 * @param geneName - the gene being translated
 * @param experimentName  - typically the CHP file name minus ".chp"
 * @param experimentChromatid1  - from the CHP file
 * @param experimentChromatid2  - from the CHP file

 *
 */
/*****************************************************************************/
CallResults::CallResults(const std::string & geneName,
                         const std::string & experimentName,
                         const CallSet & experimentChromatid1,
                         const CallSet & experimentChromatid2) :
    m_experimentChromatid1(experimentChromatid1),
    m_experimentChromatid2(experimentChromatid2),
    m_experimentName(experimentName),
    m_geneName(geneName),
    m_unknownCall("UNK/UNK")
{

  m_markerCallCount = 0;
  return;

}
// end CallResults::CallResults
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallResults::appendAlleleCall
 * Synopsis:
 *
 * Append an allele call to the std::vector of allele CallSets
 * and initialize the pair of chromatid CallSets to "Unknown".
 *
 * Typically this method is only called when the basecall is "UNK/UNK".
 *
 * @param call                 - The base call
 * @param geneCopyNumber - 0, 1 or 2
 *
 *
 */
/*****************************************************************************/
void CallResults::appendAlleleCall(const std::string & call,
                                   const unsigned int geneCopyNumber)
{

  m_alleleCalls.push_back(AlleleCall(call,
                                     m_unknownCall,
                                     m_unknownCall));


  m_experimentChromatid1.m_copyNumber = geneCopyNumber;


  return;

}
// end CallResults::appendAlleleCall
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallResults::appendAlleleCall
 * Synopsis:
 *
 * Append an allele call to the std::vector of allele CallSets
 * and intialize the second chromatid CallSet to "Unknown".
 *
 * Typically call this method only when the basecall is "?/UNK".
 *
 * The passed in CallSets represent some combination of alleles
 * from the translation table that matched the experiment data.
 * The result chromatids will never contain wildcards but represents
 * the exact allele values in the translation table
 * for a haplotype group.
 * In other words, the CallSet is simply the
 * values in the translation table corresponding to the base call
 * value, i.e "*1" is "T,T,A" for three probesets "1,2,3".
 * The reason for this "copy" is to simply enable
 * easy reporting where theses values are immediately at hand and do not need
 * to be looked up in the translation table. The translation table data
 * is shipped with each base call.
 *
 *
 * @param call                 - The base call, i.e. "*1/UNK"
 * @param geneCopyNumber       - 0, 1 or 2
 * @param resultChromatid1     - CallSet from the translation table.
 *
 */
/*****************************************************************************/
void CallResults::appendAlleleCall(const std::string & call,
                                   const unsigned int geneCopyNumber,
                                   const CallSet & resultChromatid1)
{

  m_markerCallCount++;
  m_alleleCalls.push_back(AlleleCall(call,
                                     resultChromatid1,
                                     m_unknownCall));

  m_experimentChromatid1.m_copyNumber = geneCopyNumber;

  return;
}
// end CallResults::appendAlleleCall
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallResults::appendAlleleCall
 * Synopsis:
 *
 * Append an allele call to the std::vector of allele CallSets.
 *
 * The passed in CallSets represent some combination of alleles
 * from the translation table that matched the experiment data.
 * The result chromatids will never contain wildcards but represents
 * the exact allele values in the translation table
 * for a haplotype group.
 * In other words, the CallSet is simply the
 * values in the translation table corresponding to the base call
 * value, i.e "*1" is "T,T,A" for three probesets "1,2,3".
 * The reason for this "copy" is to simply enable
 * easy reporting where theses values are immediately at hand and do not need
 * to be looked up in the translation table. The translation table data
 * is shipped with each base call.
 *
 * @param call                 - The base call
 * @param resultChromatid1     - from the translation file.
 * @param resultChromatid2     - from the translation file.
 *
 */
/*****************************************************************************/
void CallResults::appendAlleleCall(const std::string & call,
                                   const unsigned int geneCopyNumber,
                                   const CallSet & resultChromatid1,
                                   const CallSet & resultChromatid2)
{

  m_markerCallCount++;
  m_alleleCalls.push_back(AlleleCall(call,
                                     resultChromatid1,
                                     resultChromatid2));


  m_experimentChromatid1.m_copyNumber = geneCopyNumber;

  return;

}
// end CallResults::appendAlleleCall
/*****************************************************************************/
/*****************************************************************************/
/**
 * CallResults::orderCalls
 * Synopsis:
 *
 * The algorithm modifies the basecalls (i.e. "*3/ *1") to be in sorted order
 * (i.e. "*1/ *3").
 *
 * This method enforces a numerical order followed by alphabetical order
 * when applicable as "UNK" will always be guaranteed to be second.
 *
 */
/*****************************************************************************/
void CallResults::orderCalls()
{


  // Copy Number 1, cant't reorder a value like "*1"
  if (m_experimentChromatid1.m_copyNumber == 1) {
    return;
  }

  pcrecpp::RE reTrap("([^/]+)");

  for (int i = 0; i < size(); i++) {

    AlleleCall & a = m_alleleCalls[i];


    std::string call = a.m_allele;
    std::string first;
    std::string second;

    pcrecpp::StringPiece callConsume(call);

    bool okFind = reTrap.FindAndConsume(&callConsume, &first);
    APT_ERR_ASSERT(okFind, "");
    okFind = reTrap.FindAndConsume(&callConsume, &second);
    APT_ERR_ASSERT(okFind, "");

    APT_ERR_ASSERT(!first.empty() && ! second.empty(), "");

    std::string sortFirst;
    int sortFirstInt = 0;
    pcrecpp::StringPiece firstConsume(first);
    if (pcrecpp::RE("(\\d+)").FindAndConsume(&firstConsume, &sortFirst)) {
      sortFirstInt = atoi(sortFirst.c_str());
    } else {
      sortFirst = first;
    }

    std::string sortSecond;
    int sortSecondInt = 0;
    pcrecpp::StringPiece secondConsume(second) ;
    if (pcrecpp::RE("(\\d+)").FindAndConsume(&secondConsume, &sortSecond)) {
      sortSecondInt = atoi(sortSecond.c_str());
    } else {
      sortSecond = second;
    }

    if (first == "UNK") {
      call = second + "/" + first;
    } else if (second == "UNK") {
      call = first + "/" + second;
    } else if (sortFirstInt == sortSecondInt) {
      if (first.compare(second) <= 0) {
        call = first + "/" + second;
      } else {
        call = second + "/" + first;
      }
    } else if (sortFirstInt &&
               ((sortSecondInt == 0) || (sortFirstInt < sortSecondInt))) {
      call = first + "/" + second;
    } else {
      call = second + "/" + first;
    }
    a.m_allele = call;
  }

  return;

}
// end CallResults::orderCall
/*****************************************************************************/

