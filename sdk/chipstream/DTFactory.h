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
 * @file   DTFactory.h
 * @author Chuck Sugnet
 * @date   Mon Dec 14 14:12:32 2009
 * 
 * @brief  Wrapper class for all of our DataTransform factories
 * 
 */

#ifndef _DTFACTORY_H_
#define _DTFACTORY_H_

#include <string>

#include "chipstream/DataTransform.h"
#include "chipstream/ChipStreamFactory.h"
#include "chipstream/QuantMethodFactory.h"
#include "chipstream/PmAdjusterFactory.h"

/**
 * Class for creating DataTransform objects by string description. Usual
 * chipstream, pmadjuster and quantmethod string descriptions should work.
 * Classes must implement the setParameter() function to be initialized
 * properly.
 */
class DTFactory {

public:

  /** 
   * @brief Create a pointer to a new DataTransform object as described
   * in the string specification.
   * @param spec - Specification string i.e. quant-norm.sketch=1000000,
   * or med-polish.expon=true, etc.
   * 
   * @return Pointer to new DataTransform objects, must be deleted when
   * finished.
   */
   DataTransform *dataTransformForString(const std::string &spec, PsBoard &board);

private: 

  /// When building a ChipStreamDataTransform use this factory for chipstream
  ChipStreamFactory m_CSFactory;
  /// Legacy factory for building quantification Method factories
  QuantMethodFactory m_QMFactory;
  /// Legacy factory for build
  PmAdjusterFactory m_PmAdjFactory;

};

#endif /* _DTFACTORY_H_ */
