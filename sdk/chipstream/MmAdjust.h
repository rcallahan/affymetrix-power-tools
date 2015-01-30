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
 * @file   MmAdjust.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 09:30:28 2005
 * 
 * @brief  Class for supplying the mismatch intensity as pm adjustment.
 */
#ifndef _MMADJUST_H_
#define _MMADJUST_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/PmAdjuster.h"

/// String describing mm adjust algorithm/module
#define MMADJUSTSTR "pm-mm"

/**
 *  Class for supplying the mismatch intensity as pm adjustment.
 */
class MmAdjust : public PmAdjuster {

public:

  MmAdjust(){
    setupSelfDoc(*this);
    m_Type = getDocName();
  }

  /** 
   * @brief Constructor that takes the list of probe sets to remember
   * the pm/mm pairs.
   * 
   * @param layout - Annotation of probes on microarray.
   */
  void setLayout(ChipLayout &layout);

  /** 
   * @brief Setup ourselves based on data in the blackboard, specifically the
   * the pm/mm pairs.
   * 
   * @param board - blackboard with various state
   */
  void setParameters(PsBoard &board);

  /** 
   * Subset of probes that should be loaded into memory representation.
   * @param probes - bitmask of probes to be used.
   */
  void setProbes(std::vector<bool> &probes);

  /** 
   * @brief Given a PM probe intensity, supply the MM probe intensity.
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
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName(MMADJUSTSTR);
    doc.setDocDescription("Use mismatch probe as adjustment for perfect match. Has strength of being unbiased, but often the mismatch probe binds the match target.");
    // No parameters...
  }

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf() { 
    SelfDoc doc;
    setupSelfDoc(doc);
    return doc;
  }

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
  static SelfCreate *newObject(std::map<std::string,std::string> &param) {
    if(param.size() != 0) { Err::errAbort(ToStr("No parameters for ") + ToStr(MMADJUSTSTR)); }
    return new MmAdjust();
  }

private:
  /// Map of PM index on chip to MM index on chip.
  //  std::map<int, int> m_Map;
  std::vector<int> m_Vec;

};

#endif /* _MMADJUST_H_ */
