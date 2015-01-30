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
 * @file   DataTransform.h
 * @author Chuck Sugnet
 * @date   Wed Oct  7 15:44:24 PDT 2009
 *
 * @brief Base class for objects that will manipulate data and produce a new
 *        set. For example statistics might be calculated, background might be
 *        subtracted, normalization might be performed, or data might be reordered.
 */
#ifndef _DATATRANSFORM_H_
#define _DATATRANSFORM_H_

//
#include "chipstream/DataStore.h"
#include "chipstream/PsBoard.h"
//


/**
 * Base class for objects that will manipulate data and produce a new set. For
 * example statistics might be calculated, background might be subtracted,
 * normalization might be performed, or data might be reordered.
 */
class DataTransform {

public:

  /** Virtual destructor for a virtual class. */
  virtual ~DataTransform() {}

  /** 
   * Core API. Take in a blackboard, a data store and output an updated blackboard and
   * possibly modified data store.
   * 
   * @param board Blackboard state input, to be updated as dictated by transform
   * @param in    Intensity level data to be transformed
   * @param out   Any intensity data post transform
   */
  virtual bool transformData(PsBoard &board, const DataStore &in, DataStore &out) = 0;

  /** 
   * Get a unique identifier for this data transform.
   * 
   * @return string that identifies the type of transform
   */
  virtual std::string getName() = 0;
}; 

#endif /* DATATRANSFORM_H */
