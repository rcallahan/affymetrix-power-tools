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
 * @file   PmAdjusterFactory.cpp
 * @author Chuck Sugnet
 * @date   Wed Jan 11 14:42:45 2006
 * 
 * @brief Factory class for making chip streams based on a string
 * representation.
 */

//
#include "chipstream/PmAdjusterFactory.h"
//
#include "chipstream/GcAdjust.h"
#include "chipstream/MmAdjust.h"
#include "chipstream/PmSumAdjust.h"
//
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Util.h"
//

using namespace std;

/** 
 * @brief Constructor. Registers the objects we know how to create and how
 * they describe themselves.
 */  
PmAdjusterFactory::PmAdjusterFactory() {
  // No adjustment
  m_Docs.push_back(PmOnlyAdjust::explainSelf());
  m_Creators.push_back(&PmOnlyAdjust::newObject);
  
  // MM adjustment
  m_Docs.push_back(MmAdjust::explainSelf());
  m_Creators.push_back(&MmAdjust::newObject);

  // GCBG adjustment
  m_Docs.push_back(GcAdjust::explainSelf());
  m_Creators.push_back(&GcAdjust::newObject);

  // PM Allele sum adjustment
  m_Docs.push_back(PmSumAdjust::explainSelf());
  m_Creators.push_back(&PmSumAdjust::newObject);
}

/** 
 * @brief Create a pointer to a new PmAdjuster object as described
 * in the string specification.
 * @param spec - Specification string i.e. pm-gcbg
 * 
 * @return Pointer to new PmAdjuster objects, must be deleted when
 * finished.
 */
PmAdjuster *PmAdjusterFactory::adjusterForString(const std::string &spec) {
  PmAdjuster *stream = NULL;
  SelfCreate *create = NULL;
  create = SelfCreate::selfCreateFromString(spec, m_Docs, m_Creators, "PmAdjuster");
  /* Check class type. */
  if(InstanceOf(create, PmAdjuster)) {
    stream = static_cast<PmAdjuster *>(create);
  }
  else {
    Err::errAbort("Class doesn't appear to be of type PmAdjuster.");
  }
  return stream;
}

/** 
 * @brief Factory for creating PmAdjuster from string representation.
 * 
 * @param spec - Name of PmAdjuster.
 * @param layout - Probe and probe set info.
 * 
 * @return PmAdjuster requested.
 */
PmAdjuster *PmAdjusterFactory::pmAdjusterForString(const std::string &spec, ChipLayout &layout) {
  PmAdjusterFactory factory;
  PmAdjuster *adjuster = factory.adjusterForString(spec);
  assert(adjuster);
  /* Some post processing for specific objects. */
  if(InstanceOf(adjuster, GcAdjust)) {
    if(m_GcControlProbes.empty()) 
      Err::errAbort("Must specify a .bgp file when using GcAdjust");
    static_cast<GcAdjust *>(adjuster)->setLayout(layout, m_GcControlProbes);
  }
  if(InstanceOf(adjuster, PmSumAdjust)) {
    static_cast<PmSumAdjust *>(adjuster)->setLayout(layout);
  }
  else if(InstanceOf(adjuster, MmAdjust)) {
    static_cast<MmAdjust *>(adjuster)->setLayout(layout);
  }
  return adjuster;
}

/** 
 * @brief Factory for creating PmAdjuster from string representation.
 * 
 * @param spec - Name of PmAdjuster.
 * @param board - Blackboard with various state info
 * 
 * @return PmAdjuster requested.
 */
PmAdjuster *PmAdjusterFactory::pmAdjusterForString(const std::string &spec, PsBoard &board) { // ChipLayout &layout) {
  PmAdjusterFactory factory;
  PmAdjuster *adjuster = factory.adjusterForString(spec);
  assert(adjuster);
  /* Some post processing for specific objects. */
  adjuster->setParameters(board);
  return adjuster;
}

