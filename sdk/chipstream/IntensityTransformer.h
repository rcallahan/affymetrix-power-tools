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
 * @file   IntensityTransformer.h
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:19:29 2005
 * 
 * @brief  Abstract base class for objects that transform intensity data.
 */
#ifndef _INTENSITYTRANSFORMER_H_
#define _INTENSITYTRANSFORMER_H_

//
#include "chipstream/AptTypes.h"
//
#include <cstring>
#include <string>
//

/**
 *  IntensityTransformer - Abstract base class for objects that transform
 * intensity data.  
 */
class IntensityTransformer {

public:

  /** 
   * @brief Virtual destructor for a virtual class.
   */
  virtual ~IntensityTransformer() {}

  /** 
   * @brief transform the intensity point supplied coming from a particular
   * probe in a particular microarray.
   * 
   * @param probeIx - Probe index from the cel file.
   * @param chipIxs - Set of chip indexes from same sample.
   * @param intensity - Original intensities.
   * @param return - Transformed intensities.
   */
  virtual float transform(int probeIx, int chipIx, float intensity) = 0;
  
  /** 
   * @brief What is the name of the transformer?
   * @return name of transformer.
   */
  virtual std::string getType() { return m_Type; }

protected:
  /// String description of our type.
  std::string m_Type;
  
};

#endif /* _INTENSITYTRANSFORMER_H_ */
