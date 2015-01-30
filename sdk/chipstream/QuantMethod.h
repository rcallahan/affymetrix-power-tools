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
 * @file   QuantMethod.h
 * @author Chuck Sugnet
 * @date   Tue Jan 17 21:32:26 2006
 *
 * @brief  Abstract interface for a quantification method.
 */

#ifndef _QUANTMETHOD_H_
#define _QUANTMETHOD_H_

//
//#include <ctime>
#include "chipstream/AnalysisInfo.h"
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/PmAdjuster.h"
#include "chipstream/PmOnlyAdjust.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/PsBoard.h"
//
#include "util/Verbose.h"

//
#include <cstring>
#include <string>
#include <vector>


/**
 *  QuantMethod - Interface for computing quantification summaries from PM intensities
 * grouped into probe set groups.
 * 
 */
class QuantMethod : public SelfDoc, public SelfCreate {

public:

  /** Different scales that the summary results can be in. */
  enum Scale {
    Linear,   ///< Same scale as cel file, like plier results.
    Log2,     ///< Log base 2 scale, like rma results.
    Pvalue,   ///< Scale is p-value (between 0 and 1 inclusive)
    NegLog10  ///< Scale is -1 * log_10(value), useful for uncompressing p-values
  };

  /** 
   * Convert a Scale enumeration to a text description. 
   * @param s - Enumerated value.
   * @return - String version of enumeration value.
   */
  static std::string scaleToTxt(enum Scale s);

  /** Different types of expression quantification methods. */
  enum QuantType {
    Summarize, ///< Produces a summarized intensity estimate like rma or plier.
    Detection, ///< Produces a absent/present call and confidence
    GenoCall   ///< Produces a genotyping call.
  };

  /** 
   * Convert a QuantType enumeration to a human readable string representation.
   * @param q - QuantType to convert.
   * @return - Text version of QuantType enumeration.
   */
  static std::string quantTypeToTxt(enum QuantType q);

 /** 
   * @brief Virtual destructor for a virtual class.
   */
  virtual ~QuantMethod();

  /** 
   * @brief Set up the quantification method given all the data about the probe
   * set, chip layout and data.
   * 
   * @param psGroup - Probes to be used for final estimate.
   * @param iMart - Raw data from chips.
   * @param iTrans - Transformations to be applied to data before use.
   * @param pmAdjust - How to estimate background, or MM probe.
   * @return True if setup sucessful, false otherwise.
   */
  virtual bool setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart,
                     std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust) = 0;

  /** 
   * @brief Set up the quantification method given all the data about the probe
   * set, chip layout and data. Avoids issues of transformations and
   * adjustments for a more streamlined approach.
   * 
   * @param psGroup - Probes to be used for final estimate.
   * @param iMart - Raw data from chips.
   * @return True if setup sucessful, false otherwise.
   */
  virtual bool setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart);

  /** 
   * @brief Do the heavy lifting of estimation.
   */
  virtual void computeEstimate() = 0;

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
  static float transformPrimaryData(probeid_t probeIx, 
                                    chipid_t chipIx, 
                                    const IntensityMart &iMart, 
                                    std::vector<ChipStream *> &iTrans,
				    unsigned int channelIx = 0);

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
  static float attenuateIntensity(float pmI, float mmI, float l=0.005, float h=-1);

  /** 
   * @brief What is the name of the quantification method?
   * @return name of adjuster.
   */  
  virtual std::string getType() = 0;
  
  /** 
   * What version of the expression algorithm are we implementing?
   * @return version string. something like "1.0" but different versions
   * must be implemented in derive algorithm classes
   */
  virtual std::string getVersion() = 0;

  /** 
   * What scale are the target effects on? Are they linear (like
   * plier) or log2 (like rma)
   * @return - What scale are the results provided in?
   */
  virtual enum Scale getScale();

  /** 
   * What type of quantification method is this? Summarization like
   * plier detection like dabg.
   * @return - Type of quantification method.
   */
  virtual enum QuantType getQuantType();

  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param layout - Where the probesets, probes, etc are on the chip.
   * @param iMart  - intensity holder
   * 
   * @return true if success, false otherwise.
   */
  virtual bool prepare(const IntensityMart &iMart);

  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool finish();

  /**
   * Custom configuration for this QuantMethod
   */
  virtual void setParameters(PsBoard &board) {
    Err::errAbort("setParameters() is not supported in class: " + getDocName());
  }

  void setAnalysisInfo(AnalysisInfo &info) {
    m_Info = info;
  }

protected:
  AnalysisInfo m_Info;

};

#endif /* _QUANTMETHOD_H_ */
