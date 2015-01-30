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
 * @file   CallResults.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 10:06:39 PDT 2008
 * @brief  Binary representation of allele translation call.
 */

#ifndef TRANSLATION_CALLRESULTS_H
#define TRANSLATION_CALLRESULTS_H

#include "translation/CallSet.h"


class CallResults
{
public:

  class AlleleCall
  {
  public:
    std::string  m_allele;           // The translation result.
    CallSet      m_resultChromatid1; // From the Translation Table.
    CallSet      m_resultChromatid2; // From the Translation Table.

    AlleleCall(const std::string & allele, const CallSet & resultChromatid1,
             const CallSet & resultChromatid2) :
      m_allele(allele), m_resultChromatid1(resultChromatid1),
      m_resultChromatid2(resultChromatid2) {}

  };


  // With wild cards there can be many calls (m_alleleCalls).
  // Without wild cards there will ever only be one.
  std::vector <AlleleCall> m_alleleCalls;
  CallSet                  m_experimentChromatid1; // From the Genotype data
  CallSet                  m_experimentChromatid2; // From the Genotype data
  std::string              m_experimentName;
  std::string              m_geneName;
  int                      m_markerCallCount;

  CallResults(const std::string & geneName,
              const std::string & experimentName,
              const CallSet & experimentChromatid1,
              const CallSet & experimentChromatid2);

  void         appendAlleleCall(const std::string & allele,
                        const unsigned int geneCopyNumber);
  void         appendAlleleCall(const std::string & allele,
                        const unsigned int geneCopyNumber,
                        const CallSet & resultChromatid1);
  void         appendAlleleCall(const std::string & allele,
                        const unsigned int geneCopyNumber,
                        const CallSet & resultChromatid1,
                        const CallSet & resultChromatid2);
  ADT_CALL_TYPE_ENUM
              getCallType() const {
    return m_experimentChromatid1.getCallType();
  }
  std::string getMarkerProbeSet() const {
    return m_experimentChromatid1.getMarkerProbeSet();
  }
  void        orderCalls();

  int         size() const {  return m_alleleCalls.size(); }

private:

  CallSet m_unknownCall;

};

#endif /* TRANSLATION_CALLRESULTS_H */
