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
 * @file   ExperimentGeneResults.cpp
 * @author Mybrid Spalding
 * @date   Tue Mar 18 14:42:35 PDT 2008
 * @brief  CallResults wrapper class for all gene translated experiment results to be consumed by various report objects.
 */

#include "translation/ExperimentGeneResults.h"
//
#include "translation/GeneCall.h"
//

using namespace std;

/*****************************************************************************/
/**
 * ExperimentGeneResults::ExperimentGeneResults
 * Synopsis:
 *
 * Main constructor that will simply contain CallResults in a list
 * built by using "appendCallResults".
 *
 *
 * @param geneName        - the gene name
 * @param experimentName  - the corresponding experiment
 * @param zeroCopyHaplotype - a gene can only ever have one "allele" name for the zero copy number designation, and this is it.
 *
 */
/*****************************************************************************/
ExperimentGeneResults::ExperimentGeneResults(const std::string & geneName, const std::string & experimentName, const std::string & zeroCopyHaplotype):
    m_experimentName(experimentName), m_geneName(geneName), m_geneZeroCopyHaplotype(zeroCopyHaplotype)
{
  m_markerCallCount = 0;

}
// end ExperimentGeneResults::ExperimentGeneResults
/*****************************************************************************/
/*****************************************************************************/
/**
 * ExperimentGeneResults::appendCallResults
 * Synopsis:
 *
 * Append the CallResults for a pair of chromatids to the std::vector of
 * CallResults. Does some book keeping as well by building indexes
 * in place to be used later at report time.
 *
 *
 * @param newCallResults - The results to append.
 * @param geneCall      - The variantSet and reference CallSet are needed
 *                         for multi-allelic reporting on variant bases.
 *
 */
/*****************************************************************************/
void ExperimentGeneResults::appendCallResults(const CallResults & newResults,
    GeneCall & geneCall)
{

  m_callResults.push_back(newResults);
  m_markerCallCount += newResults.m_markerCallCount;

  std::map<std::string, CallElement>::iterator iCEit;
  std::map<std::string, std::vector<std::string> >::iterator jit;
  std::map<std::string, CallElement>::iterator kCEit;

  if (newResults.getCallType() == ADT_CALL_TYPE_MARKER) {
    m_probeSetToMarkerCallResults[ newResults.getMarkerProbeSet()]
    = m_callResults.size()  - 1;
  }

  for (iCEit = geneCall.getVariantCallSet().m_ceSet.begin();
       iCEit != geneCall.getVariantCallSet().m_ceSet.end();
       iCEit++) {

    jit = m_probeSetVariants.find(iCEit->first);

    if (jit ==  m_probeSetVariants.end()) {
      for (int i = 0; i < iCEit->second.m_bases.size(); i++) {
        m_probeSetVariants[iCEit->first].push_back(iCEit->second.m_bases[i]);
      }
    } else {
      for (int i = 0; i < iCEit->second.m_bases.size(); i++) {
        std::vector<std::string>::iterator kit;
        bool okToAdd = true;

        for (kit = m_probeSetVariants[iCEit->first].begin();
             (kit != m_probeSetVariants[iCEit->first].end()) &&
             okToAdd;
             kit++) {
          if (*kit == iCEit->second.m_bases[i]) {
            okToAdd = false;
          }
        }
        if (okToAdd) {
          m_probeSetVariants[iCEit->first].push_back(iCEit->second.m_bases[i]);
        }
      }
    }
  }

  for (kCEit = geneCall.getReferenceCallSet().m_ceSet.begin();
       kCEit != geneCall.getReferenceCallSet().m_ceSet.end();
       kCEit++) {

    jit = m_probeSetReference.find(kCEit->first);

    if (jit ==  m_probeSetReference.end()) {
      for (int i = 0; i < kCEit->second.m_bases.size(); i++) {
        m_probeSetReference[kCEit->first].push_back(kCEit->second.m_bases[i]);
      }
    } else {
      for (int i = 0; i < kCEit->second.m_bases.size(); i++) {
        std::vector<std::string>::iterator kit;
        bool okToAdd = true;

        for (kit = m_probeSetReference[kCEit->first].begin();
             (kit != m_probeSetReference[kCEit->first].end()) &&
             okToAdd;
             kit++) {
          if (*kit == kCEit->second.m_bases[i]) {
            okToAdd = false;
          }
        }
        if (okToAdd) {
          m_probeSetReference[kCEit->first].push_back(kCEit->second.m_bases[i]);
        }
      }
    }
  }


  return;
}

// end ExperimentGeneResults::appendCallResults
/*****************************************************************************/

