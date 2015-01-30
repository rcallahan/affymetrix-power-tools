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
  * @file   CallElement.h
  * @author Mybrid Spalding
  * @date   Thu Oct 30 15:43:35 2008
  *
  * @brief  Atomic class element for translation.
  *
  *
  */

#ifndef TRANSLATION_CALL_ELEMENT_H
#define TRANSLATION_CALL_ELEMENT_H

#include "translation/RunTimeEnvironment.h"
//
#include <cstring>
#include <string>
#include <vector>
//

class CallElement
{

public:


  std::vector<std::string>  m_bases;
  unsigned int              m_chromatid;
  unsigned int              m_copyNumber;
  bool                      m_hasReverseComplement;
  std::string               m_probeSet;
  std::vector<std::string>  m_rbases;

  CallElement();

  CallElement(const std::string & probeSet,  const std::string & base,
              const unsigned int copyNumber = 2, const bool hasRevComp = false);

  ~CallElement() {};

  bool               addBase(const std::string & base);
  std::string        basesToString(bool reverseComplement = false) const;
  void               describeVerbose(const RunTimeEnvironment & rte,
                                     ADT_VERBOSE_ENUM override = ADT_VERBOSE_NULL) const;
  void               describeVerbose(ADT_VERBOSE_ENUM level) const;
  int                getCompleteSetSize()  const { return m_completeSetSize; };
  bool               isValidBase(const std::string &);
  bool               isValidProbeSet(const std::string &);
  int                setCompleteSetSize();
  bool               hasWildCards() const {
    return m_hasWildCards;
  }
  static bool        isWildCardBase(std::string base);
  bool               operator==(CallElement const &) const;
  bool               operator!=(CallElement const & b) const {
    return !(*this == b);
  }
  static std::string reverseComplement(const std::string seq);


private:
  bool               m_hasWildCards;
  int                m_completeSetSize;

};


#endif /* TRANSLATION_CALLELEMENT_H */
