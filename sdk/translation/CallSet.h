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
 * @file   CallSet.h
 * @author Mybrid Spalding
 * @date   Tue Mar 18 09:00:56 PDT 2008
 * @brief  A little more functionaliy than using the STL 'set' container for CallElement objects.
 */


#ifndef TRANSLATION_CALLSET_H
#define TRANSLATION_CALLSET_H

#include "translation/CallElement.h"
//
#include <set>


//! The call type as dictated by the translation file.
enum ADT_CALL_TYPE_ENUM {
  ADT_CALL_TYPE_NULL,
  ADT_CALL_TYPE_MARKER,
  ADT_CALL_TYPE_HAPLOTYPE_GROUP,
  ADT_CALL_TYPE_COPY_NUMBER,
};


class CallSet
{
public:

  //! \cond IGNORE_OBVIOUS
  std::map<std::string, CallElement >   m_ceSet; 
  int                                   m_columnInTranslationTable;
  unsigned int                          m_copyNumber;
  bool                                  m_hasZeroCopyNumberNotAvailable;
  bool                                  m_hasWildCards;
  bool                                  m_isDescriptive;
  bool                                  m_isMultiAllelic;
  std::string                           m_name;
  ADT_CALL_TYPE_ENUM                    m_type;
  //! \endcond

  CallSet( const std::string & name = "", const int & columnInTranslationTable = -1);
  
  CallSet(const class RunTimeEnvironment & rte,
          class GenotypeTableModel & gtm,
          const CallSet & referenceSet,
          const int chromaTid);


  bool               addCallElement(class TranslationTableModel & ttm,
                                    const int row,
                                    const int base_column,
                                    const bool allowMultiAllelic = false,
                                    const bool isHaplotype = false);


  CallSet            buildSecondSet(const RunTimeEnvironment & rte,
                                   const CallSet & chromatid1,
                                   const CallSet & chromatid2) const;
  void               clear();
  bool               contains(const CallSet & c) const;
  bool               contains(const CallElement & ce) const;
  void               describeVerbose(ADT_VERBOSE_ENUM level) const;
  void               describeVerbose(const RunTimeEnvironment & rte,
                                     ADT_VERBOSE_ENUM override =
                                       ADT_VERBOSE_NULL) const;
  //! inline
  ADT_CALL_TYPE_ENUM getCallType() const { return m_type; }; 
  std::string        getMarkerProbeSet() const;
  const CallElement & getProbeSetCallElement(const std::string & probeSet) const;
  bool               hasProbeSet(const std::string & probeSet) const;
  //! inline
  bool               hasWildCards() const { return m_hasWildCards; } 
  //! inline
  bool               isEmpty() const {  return (m_ceSet.size() == 0); } 
  bool               match(const CallSet & a, const CallSet & b) const;
  int                numWildCards() const;
  bool               operator==(const CallSet & cs) const;
  //! inline
  bool               operator!=(const CallSet & cs) const { return !(*this == cs); } 
  //! inline
  int                size() const { return m_ceSet.size(); } 

private:
  void               _setCopyNumber(const RunTimeEnvironment & rte,
                                    class GenotypeTableModel & gtm,
                                    const int row);
  void               _setCopyNumber(const std::string & base);

};

#endif /* TRANSLATION_CALLSET_H */
