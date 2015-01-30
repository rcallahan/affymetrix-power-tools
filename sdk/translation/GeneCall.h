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
 * @file   GeneCall.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:04:13 PDT 2008
 * @brief  Business object that is the heart of the matter, class with the algorithm that translates experiment data.
 */

#ifndef TRANSLATION_GENECALL_H
#define TRANSLATION_GENECALL_H


#include "translation/CallSet.h"
#include "translation/RunTimeEnvironment.h"

class GeneCall
{
public:



  std::vector<CallSet> m_alleleSet;
  std::string          m_gene;

  GeneCall(const RunTimeEnvironment & rte,
           class TranslationTableModel & ttm,
           const int headerRow,
           bool isHaplotypeCall);


  bool               addCallElementToAlleleCallSet(
    class TranslationTableModel & ttm,
    const int headerRow,
    const int row,
    const int column,
    bool allowMultiAllelic,
    bool rowIsHaplotype);
  void               describeVerbose(const RunTimeEnvironment & rte, ADT_VERBOSE_ENUM level = ADT_VERBOSE_NULL);
  bool               existsAlleleCallSet(const std::string & alleleName);
  CallSet &          getAlleleCallSet(unsigned int index);
  CallSet &          getReferenceCallSet() {
    return getAlleleCallSet(0);
  }
  CallSet &          getVariantCallSet()   {
    return getAlleleCallSet(1);
  }
  class CallResults  translateExperimentCall(RunTimeEnvironment & rte,
            const CallSet & completeSet,
            GenotypeTableModel & gtm,
            CallSet & chromatid1,
            CallSet & chromatid2,
            const std::string & copyNumberCall,
            const CallSet & copyNumberZeroCallSet);


private:


};

#endif /* TRANSLATION_GENECALL_H */
