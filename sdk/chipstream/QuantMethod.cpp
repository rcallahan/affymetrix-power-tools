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
 * @file   QuantMethod.cpp
 * @author Chuck Sugnet
 * @date   Tue Jan 17 21:32:26 2006
 *
 * @brief  Abstract interface for a quantification method.
 */

#include "chipstream/QuantMethod.h"
//

/** 
 * Convert a Scale enumeration to a text description. 
 * @param s - Enumerated value.
 * @return - String version of enumeration value.
 */
std::string QuantMethod::scaleToTxt(enum QuantMethod::Scale s) {
  switch (s) {
  case Linear : 
    return "linear";
  case Log2 :
    return "log2";
  case Pvalue :
    return "p-value";
  case NegLog10 :
    return "neg-log10";
  default:
    Err::errAbort("QuantMethod::scaleToTxt() - Unknown type: '" + ToStr(s) + "'");
  }
  return ""; // avoid compiler warnings
}


/** 
 * Convert a QuantType enumeration to a human readable string representation.
 * @param q - QuantType to convert.
 * @return - Text version of QuantType enumeration.
 */
std::string QuantMethod::quantTypeToTxt(enum QuantMethod::QuantType q) {
  switch (q) {
  case Summarize : 
    return "signal";
  case Detection :
    return "detection";
  case GenoCall :
    return "genotype";
  default :
    Err::errAbort("QuantMethod::quantTypeToTxt() - Unknown type: '" + ToStr(q) +"'");
  }
  return ""; // avoid compiler warnings
}

/** 
 * @brief Virtual destructor for a virtual class.
 */
QuantMethod::~QuantMethod() {
};

/** 
 * @brief Set up the quantification method given all the data about the probe
 * set, chip layout and data. Avoids issues of transformations and
 * adjustments for a more streamlined approach.
 * 
 * @param psGroup - Probes to be used for final estimate.
 * @param iMart - Raw data from chips.
 * @return True if setup sucessful, false otherwise.
 */
bool QuantMethod::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart) {
  std::vector<ChipStream *> iTrans;
  PmOnlyAdjust pmAdjust;
  return setUp(psGroup, iMart, iTrans, pmAdjust);
}

/** 
 * @brief Transform raw data for probe index on chip and chip index using the
 * chip stream intensity transformers.
 * 
 * @param probeIx - Probe index on chip.
 * @param chipIx - Chip index in experiment.
 * @param iMart - Raw data. 
 * @param iTrans - Transformations to be applied to data.
 * 
 * @return 
 */
float QuantMethod::transformPrimaryData(probeid_t probeIx, 
                                        chipid_t chipIx, 
                                        const IntensityMart &iMart, 
                                        std::vector<ChipStream *> &iTrans,
                                        unsigned int channelIx) {
  float intensity;
  if (iTrans.empty()) {
    intensity = iMart.getProbeIntensity(probeIx, chipIx, channelIx);
  }
  else {
    intensity = iTrans[iTrans.size() - 1]->getTransformedIntensity(probeIx, chipIx, channelIx);
  }
  return intensity;
}

/** 
 * @brief Returns PM minus attenuated mismatch value so such that the return value is 
 * greater than or equal to zero. In short, we return (x+sqrt(x^2+H))/2 where 
 * H = 4*pM*MM*L and X = PM-MM. This is in effect generalized log without the log
 * transformation.
 * 
 * @param pmI - perfect match intensity
 * @param mmI - mismatch intensity
 * @param l - optional tunable parameter, default is 0.005 must be between 0 and 1 inclusive
 * @param h - optional parameter, allows one to override computation of H from
 * data. If less than zero computed from data as 4*PM*MM*l
 * 
 * @return 
 */
float QuantMethod::attenuateIntensity(float pmI, float mmI, float l, float h) {
  h = h < 0 ? 4*pmI*mmI*l : h; // if h is less than zero, calculate default value: H = 4*PM*MM*L
  // X = PM-MM
  // glog(X,H) = log( (X+sqrt(X^2+H))/2)
  // we return the non log value here
  return ( (pmI-mmI)+sqrt(((pmI-mmI)*(pmI-mmI)) + h) )/2;

}

/** 
 * What scale are the target effects on? Are they linear (like
 * plier) or log2 (like rma)
 * @return - What scale are the results provided in?
 */
enum QuantMethod::Scale QuantMethod::getScale() {
  return Linear;
}

/** 
 * What type of quantification method is this? Summarization like
 * plier detection like dabg.
 * @return - Type of quantification method.
 */
enum QuantMethod::QuantType QuantMethod::getQuantType() {
  return Summarize;
}

/** 
 * Get set up for a run of reporting probesets. Often used to open file
 * streams and print headers to files etc.
 * 
 * @param layout - Where the probesets, probes, etc are on the chip.
 * @param iMart  - intensity holder
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethod::prepare(const IntensityMart &iMart) {
  return true;
}

/** 
 * No more probesets will be processed, this is a chance to finish outputting
 * results and clean up.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethod::finish() {
  return true;
}
