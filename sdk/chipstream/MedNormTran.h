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
 * @file   MedNormTran.h
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:40:14 2005
 * 
 * @brief  Class for doing median, or average, normalization
 */
#ifndef _MEDNORMTRAN_H_
#define _MEDNORMTRAN_H_

//
#include "chipstream/ChipStream.h"
//
#include "util/Convert.h"
//
#include <vector>
//

/// String describing median normalization algorithm
#define MEDNORMSTR "med-norm"

/**
 *  MedNormTran -  Class for doing median, or average, normalization
 */
class MedNormTran : public ChipStream {

public:

  /** 
   * Constructor
   * 
   * @param target - Value to set median or average value to.
   * @param doAverage - Do average instead of median.
   * @param calcTarget - Calculate target as average or median of all data.
   * @param lowPrecision - Should we pretend like this got wrote to a cel file and 
   *                       read back again (meaning truncation)
   */
  MedNormTran(float target, bool doAverage, bool calcTarget, bool lowPrecision) ;

  /** 
   * @brief Only use a subset of all probes for normalization. For
   * example RMA only uses PM probes or might want to just use control
   * genes.
   * @param subsetProbes - bitmask indicating which probes to use for
   * normalization.
   */
  void setSubProbes(const std::vector<bool> &subsetProbes);

  /** 
   * @brief transform the intensity point supplied coming from a particular
   * probe in a particular microarray.
   * 
   * @param probeIx - Probe index from the cel file.
   * @param chipIx - Set of chip indexes from same sample.
   * @param intensity - Original intensities. 
   * @param return - transformed intensities.
   */
  float transform(int probeIx, int chipIx, float intensity);

  void transform(int chipIx, std::vector<float>& intensity);

  void newChip(std::vector<float> &data);

  /** 
   * @brief Method for being passed a new cel file worth of data.
   * @param data - Vector of vectors of cel file data from same sample.
   */  
  void newDataSet(IntensityMart* iMart);

  /** 
   * @brief Method for being passed a new cel file worth of data. Calculates
   * and stores the median (or average) of data supplied.
   * @param data - cel file vector of data.
   */  
  void initializeData(const std::vector<float> &data);

  /** 
   * @brief Signal that no more data is coming (i.e. newChip() will
   * not be called anymore. Currently the class builds up all the
   * summaries and figures out the normalization factores when this
   * function is called.  This makes it difficult to pass through data
   * to downstream.
   */
  void setTarget();

  void finishedChips() {
    setTarget();
  }

  /**
   * @brief Perform any necessary cleanup
   */
  void endDataSet();

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions();

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

protected:
  /// Summary (median or average) for each chip seen.
  std::vector<double> m_Summary;
  /// Scale factore to multiply each chip by to get target
  std::vector<double> m_Scale;
  /// Bitmask indicating which probes to use for normalization
  std::vector<bool> m_SubProbes;
  /// Number of probes that will be used for normalization
  int m_SubProbeCount;
  /// What we want the summary for each chip to be.
  double m_TargetNorm;
  /// Do the average rather than the median.
  bool m_DoAverage;
  /// Has the target been set yet?
  bool m_TargetUnset;
  /// Should we mimic the truncation seen in normalized cel files?
  bool m_LowPrecision;
};

#endif /* _MEDNORMTRAN_H_ */
