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
 * @file   QuantIterPlier.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 11:44:24 2005
 * 
 * @brief Class for doing probe set quantification estimate by iteratively
 * calling PLIER with the probes that best correlate with signal estimate.
 */
#ifndef _QUANTITERPLIER_H_
#define _QUANTITERPLIER_H_

//
#include "chipstream/QuantPlierBase.h"
//
#include "plier/IterPlier.h"
//

/// String describing iterplier summary method
#define QUANTITERPLIERSTR "iter-plier"

/**
 *  Class for doing probe set quantification estimate by iteratively
 * calling PLIER with the probes that best correlate with signal estimate.
 * 
 */
class QuantIterPlier : public QuantPlierBase {

public:

  /** 
   * What version of the algorithm are we implementing?
   * @return version string.
   */
  std::string getVersion() {
    return "2.0";
  }

  /** 
   * Constructor.  
   */
  QuantIterPlier(QuantPlierParams &plierParams);
  /**
   * Destructor.
   */
  ~QuantIterPlier();

  /** 
   * @brief How many times and with how many probes are we running with.
   * 
   * @param iterations - Vector of number of probes to
   * use at each iteration, i.e. 22,11. Code implicitly uses all
   * probes for first iteration. So supplying {22,11} actually does
   * all probes, then 22, then 11.
   */
  void setIterations(vector<int> &iterations) {
    m_Iterations = iterations;
  }

  /**
   * Custom configuration for this QuantMethod
   */
  void setParameters(PsBoard &board) {
    vector<int> iters;
    iters.push_back(22);
    iters.push_back(11);
    setIterations(iters);
  }


  /** 
   * @brief Do the heavy lifting of plier.
   */
  void computeEstimate();

  /** 
   * @brief Return number of features (probes).
   * @return feature count.
   */
  inline unsigned int getNumFeatures()  { return m_ProbeCount; }

  /** 
   * @brief Was this probe used for the estimate?
   * @param probeIx - index of probe in probe set.
   * @return true if feature (probe) was used, false otherwise.
   */
  inline bool featureUsed(unsigned int probeIx) {
    assert(probeIx < m_ProbeCount);
    return m_ProbesUsed[probeIx] == 1;
  }

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName(QUANTITERPLIERSTR);
    doc.setDocDescription("Do probe set quantification estimate by iteratively calling PLIER with the probes that best correlate with signal estimate. The version of PLIER used by IterPLIER differs from the previous version by the addition of a SafteyZero, NumericalTolerance, and FixPrecomputed. These options are intended to improve the stability of PLIER results when using precomputed feature reponse values. To get the older PLIER behavior set SafetyZero to 0.0, NumericalTolerance to 0.0, and FixPrecomputed to false.");
    /* Set common parameters and documentation used by both QuantPlier and QuantIterPlier. */
    doc.setDocOptions(getDefaultDocOptions());
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
    QuantPlierParams plierParams;
    SelfDoc doc = explainSelf();
    fillPlierParams(param, plierParams, doc);
    QuantIterPlier *plier = new QuantIterPlier(plierParams);
    plier->setDocOptions(doc.getDocOptions());
    return plier;
  }

private :
  void initMemory();
  /** 
   * @brief Allocate enough memory for at least maxChips and
   * maxProbes.
   * 
   * @param maxChips - Number of chips that can be computed.
   * @param maxProbes - Number of probes or features that can be computed.
   */
  void allocMemory(unsigned int maxChips, unsigned int maxProbes);

  /** 
   * @brief Free up all the memory that has bee allocated.
   */
  void freeMemory();

  /** 
   * @brief Specify precomputed feature effects to be used as map from
   * probe id index to the feature effect.  This is not supported
   * for iter-plier.
   * 
   * @param featureEffects - mapping from probe idex to feature effects.
   */
  void setFeaturePriorEffects(double* featureEffects, int iSize) {
    Err::errAbort("setFeaturePriorEffects can not be called in iter-plier");
  }

  /// How many iterations to do at each level.
  vector<int> m_Iterations;
  /// Array with 1 if probes used, 0 otherwise.
  int *m_ProbesUsed;
  /// Instance of Iterplier to be used.
  IterPlier m_IterPlier;
};

#endif /* _QUANTITERPLIER_H_ */
