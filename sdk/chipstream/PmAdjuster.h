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
 * @file   PmAdjuster.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 09:16:57 2005
 * 
 * @brief Interface for determining a change based on intensity of perfect match
 * probe.
 */
#ifndef _PMADJUSTER_H_
#define _PMADJUSTER_H_

//
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/IntensityTransformer.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
//
#include <vector>

/**
 *  Interface for determining a change based on intensity of perfect match
 */
class PmAdjuster : public SelfDoc, public SelfCreate {

public:

  /** 
   * @brief Virtual destructor for a virtual class.
   */
  virtual ~PmAdjuster() {};

  /** 
   * @brief Function to determine how much to adjust a perfect match intensity.
   * @param probeIx - Index of probe in cel file data.
   * @param chipIx - Microarray or chip index.
   * @param iMart - IntensityMart which holds raw data.
   * @param iTrans - Vector of transformations that should be performed on raw data.
   * @param pmIintensity - Intensity of perfect match probe to be adjusted, may
   * be modified from original value depending on adjuster
   * @param bgrdAdjust - Background adjustment, if any, recommended (i.e. MM intensity)
   */
  virtual void pmAdjustment(int probeIx, int chipIx, 
                            const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                            float &pmIntensity, float &bgrdAdjust) = 0;

  /** 
   * Subset of probes that should are needed for estimating parameters.
   * @param probes - bitmask of probes to be set.
   */
  virtual void setProbes(std::vector<bool> &probes) {};

  /** 
   * @brief What is the name of the adjuster?
   * @return name of adjuster.
   */  
  std::string getType() { return m_Type;}

  /**
   * Custom configuration for this PmAdjuster
   */
  virtual void setParameters(PsBoard &board) {
    Err::errAbort("setParameters() is not supported in class: " + getDocName());
  }

protected:
  /// String description of adjuster
  std::string m_Type;
};

#endif /* _PMADJUSTER_H_ */
