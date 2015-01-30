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
 * @file   PmSumAdjust.h
 * @author Chuck Sugnet
 * @date   Mon May 14 11:45:08 2007
 * 
 * @brief Class for doing an additive adjustment based on a probe that
 * may be hybridizing to another allele
 * 
 */

#ifndef _PMSUMADJUST_H_
#define _PMSUMADJUST_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/PmAdjuster.h"
#include "chipstream/QuantMethod.h"
//
#include "stats/stats.h"
#include "util/Verbose.h"
//
#include <map>
#include <vector>
//

/// String describing gc adjust algorithm/module
#define PMSUMADJUSTSTR "pm-sum"

/**
 *  Class for doing an adjustment based on the median intensity of probes
 * with similar GC content.
 */
class PmSumAdjust : public PmAdjuster {

public:
  
  /** Constructor. */
  PmSumAdjust(bool noMatchOk);

  /** 
   * @brief Constructor that takes the list of probe sets to remember
   * the pm/pm pairs.
   * 
   * @param layout - Annotation of probes on microarray.
   */
  void setLayout(ChipLayout &layout);

  /**
   * @brief Constructor that takes the list of probe sets to remember
   * the pm/pm pairs.
   *
   * @param board - blackboard with various state to access
   */
  void setParams(PsBoard &board);
  
  /** 
   * Subset of probes that should be loaded into memory representation.
   * @param probes - bitmask of probes to be used.
   */
  void setProbes(std::vector<bool> &probes);

  /** 
   * @brief Given a PM probe intensity add the matching PM intensity
   * for the other allele to the pmIntensity.
   *
   * @param probeIx - Index of probe in cel file data.
   * @param chipIx - Microarray or chip index.
   * @param iMart - IntensityMart which holds raw data.
   * @param iTrans - Vector of transformations that should be performed on raw data.
   * @param pmIintensity - Intensity of perfect match probe to be adjusted, may
   * be modified from original value depending on adjuster
   * @param bgrdAdjust - Background adjustment, if any, recommended (i.e. MM intensity)
   */
  void pmAdjustment(int probeIx, int chipIx, 
                    const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                    float &pmIntensity, float &bgrdAdjust);

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc);

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf();

  /** 
   * @brief This static function should be overridden by child classes
   * to return an object of the correct type initialized correctly
   * with the parameters in the string, string map. All objects
   * created this way should be deleted when finished using.
   * 
   * @param param - Map of key/value pairs to initialize the object.
   * 
   * @return Pointer toCreate object, this should be sub casted as necessary.
   */
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

private:
  /// If don't find matching allele probe to add ok to just add zero?
  bool m_NoMatchOk;
  /// Map of PM index on chip to PM probes from other allele index on chip.
  std::vector<int> m_Vec;

};

#endif /* _PMSUMADJUST_H_ */
