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
 * @file   ExperimentGeneResults.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 10:06:39 PDT 2008
 * @brief  CallResults wrapper class for all gene translated experiment results to be consumed by various report objects.
 */

#ifndef TRANSLATION_EXPERIMENTGENERESULTS_H
#define TRANSLATION_EXPERIMENTGENERESULTS_H

#include "translation/CallResults.h"
#include "translation/ExperimentReport.h"
#include "translation/GenotypeTableModel.h"
//


class ExperimentGeneResults
{
public:


  // Common Data Elements
  std::vector<CallResults>       m_callResults;
  std::string                    m_experimentName;
  ExperimentReportTypeEnum       m_experimentType;
  std::string                    m_geneName;
  std::string                    m_geneZeroCopyHaplotype;
  int                            m_markerCallCount;
  std::map < std::string,
  std::vector< std::string > > m_probeSetReference;
  std::map < std::string,
  std::vector< std::string > > m_probeSetVariants;
  std::map< std::string, int>    m_probeSetToMarkerCallResults;


  ExperimentGeneResults() {
    m_markerCallCount = 0;
  }

  ExperimentGeneResults(const std::string & geneName,
                        const std::string & experimentName,
                        const std::string & zeroCopyHaplotype);


  void   appendCallResults(const CallResults & newResults,
                           class GeneCall & geneCall) ;

  size_t size() const {
    return m_callResults.size();
  }


};


#endif /* TRANSLATION_EXPERIMENTGENERESULTS_H */
