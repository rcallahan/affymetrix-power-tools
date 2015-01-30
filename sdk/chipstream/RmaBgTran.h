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
 * @file   RmaBgTran.h
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:32:52 2005
 * 
 * @brief  Class for doing RMA style background subtraction.
 */
#ifndef _RMABGTRAN_H_
#define _RMABGTRAN_H_

//
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityTransformer.h"
//
#include "portability/affy-base-types.h"
#include "rma/RMA.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
//

/// String describing rma background subtraction method
#define RMABGSTR "rma-bg"
/**
 *  RmaBgTran - Class for doing RMA style background subtraction.
 */
class RmaBgTran : public ChipStream {

public:

  /** 
   * @brief Constructor.
   * @param numBins - Number of bins to use in density estimation 16384 is magic
   *                  number from biocoductor, must be a power of two.
   */
  RmaBgTran(int numBins=16384);

  /** 
   * @brief Set which probes to use for the background estimation.
   * 
   * @param pmProbes Bitmask where PM probes are set to true.
   */
  void setPmProbes(const std::vector<bool> &pmProbes);

  /**
   * Setup the background probes and background gc indexes by reading
   * them from the blackboard
   */
  void setParameters(PsBoard &board);

  /** 
   * @brief transform the intensity point supplied coming from a particular
   * probe in a particular microarray.
   * 
   * @param probeIx - Probe index from the cel file.
   * @param chipIx - Set of chips from same sample.
   * @param intensity - Original intensities.
   * @param return - transformed intensities.
   */
  float transform(int probeIx, int chipIx, float intensity);
  void transform(int chipIx, std::vector<float>& intensity);

  void newChip(std::vector<float> &data);

  /** 
   * @brief Method for being passed a new cel file worth of data.
   * @param data - vector of vectors of cel file data from same sample.
   */
  void newDataSet(IntensityMart* iMart);

/** 
 * @brief Method for adding cel file data to chipstream normalization method.
 * @param data - cel file vector of data.
 */
  void initializeData(const std::vector<float>& data);

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName(RMABGSTR);
    doc.setDocDescription("Performs an RMA style background adjustment as described in Irizarry et al 2003.");
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
    if(param.size() != 0) { Err::errAbort(ToStr("No parameters for ") + ToStr(RMABGSTR)); }
    return new RmaBgTran();
  }

private: 

  /// Parameters for an indivdual chip.
  struct ChipParam {
    double mu;     /// mode of "background" signal
    double sigma;  /// stdev of "background" signal
    double alpha;  /// mean of exponential "true" signal
  };

  /// Where parameters for each chip are stored.
  std::vector<struct ChipParam> m_Params;
  /// Bit mask where each pm probe is marked. 
  std::vector<bool> m_PmProbes;
  /// How many pm probes are there?
  int m_PmProbeCount;
  /// How many bins to use for RMA density estimation, must be power of two.
  int m_NumBins;
  
};

#endif /* _RMABGTRAN_H_ */
