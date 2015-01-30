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
 * @file   PmOnlyAdjust.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 09:28:27 2005
 * 
 * @brief  Pm adjustment class that does no modification.
 */
#ifndef _PMONLYADJUST_H_
#define _PMONLYADJUST_H_

//
#include "chipstream/PmAdjuster.h"
//
#include "util/Convert.h"
//

/// String describing pm only adjust algorithm/module
#define PMONLYSTR "pm-only"

/**
 *  Pm adjustment class that does no modification.
 * 
 */
class PmOnlyAdjust : public PmAdjuster {

public:

  /** 
   * @brief Constructor.
   */
  PmOnlyAdjust() {
    setupSelfDoc(*this);
    m_Type = getDocName();
  }

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
  void pmAdjustment(int probeIx, int chipIx, 
                    const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                    float &pmIntensity, float &bgrdAdjust) {
    bgrdAdjust = 0;
  }

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName(PMONLYSTR);
    doc.setDocDescription("No adjustment. Just uses unmodified PM intensity values.");
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
    if(param.size() != 0) { Err::errAbort(ToStr("No parameters for ") + ToStr(PMONLYSTR)); }
    return new PmOnlyAdjust();
  }

  /**
   * Custom configuration for this PmAdjuster
   */
  virtual void setParameters(PsBoard &board) {
    // deliberate no-op here to avoid calling parent classes version
  }

};

#endif /* _PMONLYADJUST_H_ */
