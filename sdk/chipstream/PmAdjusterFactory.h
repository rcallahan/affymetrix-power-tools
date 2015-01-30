////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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
 * @file   PmAdjusterFactory.h
 * @author Chuck Sugnet
 * @date   Tue Oct 25 11:56:37 2005
 * 
 * @brief Factory class for making chip streams based on a string
 * representation.
 */

#ifndef _PMADJUSTERFACTORY_H_
#define _PMADJUSTERFACTORY_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/PmAdjuster.h"
#include "chipstream/PmOnlyAdjust.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * Factory class for making chip streams based on a string representation. 
 */
class PmAdjusterFactory {

public:

  /** 
   * @brief Constructor. Registers the objects we know how to create and how
   * they describe themselves.
   */  
  PmAdjusterFactory();

  /** 
   * @brief Create a pointer to a new PmAdjuster object as described
   * in the string specification.
   * @param spec - Specification string i.e. pm-gcbg
   * 
   * @return Pointer to new PmAdjuster objects, must be deleted when
   * finished.
   */
  PmAdjuster *adjusterForString(const std::string &spec);

  /** 
   * @brief Factory for creating PmAdjuster from string representation.
   * 
   * @param spec - Name of PmAdjuster.
   * @param layout - Probe and probe set info.
   * @return PmAdjuster requested.
   */
  PmAdjuster *pmAdjusterForString(const std::string &spec, ChipLayout &layout);

  /** 
   * @brief Factory for creating PmAdjuster from string representation.
   * 
   * @param spec - Name of PmAdjuster.
   * @param board - Blackboard with various state info
   * 
   * @return PmAdjuster requested.
   */
  PmAdjuster *pmAdjusterForString(const std::string &spec, PsBoard &board);

  /** 
   * Set the probes to be used as a background distribution.
   * @param controlProbes - probes to use as controls.
   */
  void setControlProbes(const std::vector<Probe *> controlProbes) {
    m_GcControlProbes = controlProbes;
  }

  /** 
   * Get the vector of documentation objects that this factory is aware of.
   * @return - vector of documentation objects.
   */
  std::vector<SelfDoc> getDocs() {
    return m_Docs;
  }

  /** 
   * Can this factory be used to build from this specification?
   * @param spec - Specification in form "chipstream.key=value.key=value"
   * @return true if this factory knows how to build from this specification, false otherwise
   */
  bool canBuild(const std::string &spec) {
    return SelfCreate::canMake(spec, m_Docs);
  }
  
protected:

  /// Background probes to use. Note that memory is owned elsewhere.
  std::vector<Probe *> m_GcControlProbes;
  /// Self documentation
  std::vector<SelfDoc> m_Docs;
  /// Self creation
  std::vector<SelfCreate::selfCreator> m_Creators;

};

#endif /* _PMADJUSTERFACTORY_H_ */
  
