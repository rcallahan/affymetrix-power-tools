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
 * @file   ChipStreamDataTransform.h
 * @author Chuck Sugnet
 * @date   Wed Oct  7 16:45:33 PDT 2009
 *
 * @brief  Wrapper for ChipStream objects that makes them into DataTransform objects
 */
#ifndef _CHIPSTREAMDATATRANSFORM_H_
#define _CHIPSTREAMDATATRANSFORM_H_

//
#include "chipstream/DataStore.h"
#include "chipstream/ChipStream.h"
#include "chipstream/DataTransform.h"
#include "chipstream/SketchQuantNormTran.h"
#include "util/Util.h"
//

/**
 * Wrapper for ChipStream objects that makes them into DataTransform objects
 */
class ChipStreamDataTransform : public DataTransform {

public:

  /** 
   * Constructor that takes ChipStream to use for algorithmic processing.
   * 
   * @param cStream ChipStream object we are wrapping.
   */
  ChipStreamDataTransform(ChipStream *cStream) {
    m_ChipStream = cStream;
  }
  
  /**
   * Destructor - Delete our chipstream object.
   */
  ~ChipStreamDataTransform() {
    Freez(m_ChipStream);
  }

  /** 
   * Core API. Take in a blackboard, a data store and output an updated blackboard and
   * possibly modified data store.
   * 
   * @param board Blackboard state input, to be updated as dictated by trasnform
   * @param in    Intensity level data to be transformed
   * @param out   Any intensity data post transform
   */
  virtual bool transformData(PsBoard &board, const DataStore &in, DataStore &out);

  /** 
   * Core API. Take in a blackboard, a data store and output an updated blackboard and
   * possibly modified data store. This class will use the ChipStream object it wraps
   * to do the actual numerics and handles the blackboard as well as marshalling data
   * in and out of the data store.
   * 
   * @param board Blackboard state input, to be updated as dictated by transform
   * @param in    Intensity level data to be transformed
   * @param out   Any intensity data post transform
   */  
  virtual std::string getName() { return m_ChipStream->getDocName(); }
  
private:

  /// ChipSteam object that we are wrapping to get algorithmic implementation.
  ChipStream *m_ChipStream;

}; 

#endif /* CHIPSTREAMDATATRANSFORM_H */
