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
 * @file   QuantPlier.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 10:45:35 2005
 * 
 * @brief Class for doing probe set quantification using the PLIER (Probe
 * Logarithmic Error Intensity Estimate) method.
 */
#ifndef _QUANTPLIER_H_
#define _QUANTPLIER_H_

//
#include "chipstream/QuantPlierBase.h"

/// String description of quant plier chipstream.
#define QUANTPLIERSTR "plier"

/**
 *  Class for doing probe set quantification using the PLIER (Probe
 * Logarithmic Error Intensity Estimate) method.
 */
class QuantPlier : public QuantPlierBase {

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
  QuantPlier(QuantPlierParams &plierParams);

  /**
   * Destructor.
   */
  virtual ~QuantPlier();

  /** 
   * @brief Do the heavy lifting of plier.
   */
  void computeEstimate();

  /** 
   * @brief Return number of features (probes).
   * @return feature count.
   */
  unsigned int getNumFeatures();

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc);


  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.  The default
   * parameters are from the inputs.txt file in the sdk example.
   * @return SelfDoc
   */
  static SelfDoc explainSelf();

  /** 
   * @brief This static function should be overridden by child classes
   * to return an object of the correct type initialized correctly
   * with the parameters in the string, string map. All objects
   * created this way should be deleted when finished using.
   * The default values used for the plier parameters are from the
   * inputs.txt file in the sdk example.
   * 
   * @param param - Map of key/value pairs to initialize the object.
   * 
   * @return Pointer toCreate object, this should be sub casted as necessary.
   */
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

  /** 
   * @brief Specify precomputed feature effects to be used as map from
   * probe id index to the feature effect.  This make plier much
   * faster to do as doesn't have to iterate through anything, but
   * have to make sure that the feature effects being used match the
   * analysis that is being done. For example the feature effects
   * might be different using pm-only, pm-gcbg and pm-mm adjustments.
   * 
   * @param featureEffects - mapping from probe idex to feature effects.
   */
  void setFeaturePriorEffects(double* featureEffects, int iSize);

private:
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
   * @brief Set memory variables to zero/null
   */
  void zeroMemoryVars();
 

};

#endif /* _QUANTPLIER_H_ */
