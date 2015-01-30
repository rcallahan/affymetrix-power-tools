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
 * @file   QuantExprMethod.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 10:11:43 2005
 * 
 * @brief Interface for computing quantification summaries from PM intensities
 * grouped into probe set groups.
 * 
 */
#ifndef _QUANTEXPRMETHOD_H_
#define _QUANTEXPRMETHOD_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/PmAdjuster.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantMethod.h"
//
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 *  QuantExprMethod - Interface for computing quantification summaries from PM intensities
 * grouped into probe set groups.
 * 
 */
class QuantExprMethod : public QuantMethod {

public:
  
  /** 
   * @brief Virtual destructor for a virtual class.
   */
  virtual ~QuantExprMethod() {};

  /** 
   * @brief Set the PM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   * @param data - Data to be set.
   */
  virtual void setPMDataAt(unsigned int probeIx, unsigned int chipIx, double data) = 0;

  /** 
   * @brief Get the PM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   */
  virtual double getPMDataAt(unsigned int probeIx, unsigned int chipIx) = 0;

  /** 
   * @brief Set the MM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   * @param data - Data to be set.
   */
  virtual void setMMDataAt(unsigned int probeIx, unsigned int chipIx, double data) = 0;

  /** 
   * @brief Clear out all the data.
   */
  virtual void clear()  = 0;

  /** 
   * Set the number of probes and number of chips to be used for estimation.
   * @param numProbes - Probe count.
   * @param numChips - Microarray count.
   */
  virtual void setBounds(unsigned int numProbes, unsigned int numChips) = 0;

  /** 
   * @brief Return number of features (probes).
   * @return feature count.
   */
  virtual unsigned int getNumFeatures() = 0;

  /** 
   * @brief Return number of targets (experiments or chips).
   * @return target count
   */
  virtual unsigned int getNumTargets() = 0;

  /** 
   * @brief Does this quantification method estimate feature effects?
   * @return true or false if getFeatureEffects() supplies estimates.
   */
  virtual bool haveFeatureEffects() { return false; }

  /** 
   * @brief Get the feature effect of a particular probe.
   * @param probeIx - Index of probe in probe sets.
   * @return feature (probe) effects.
   */
  virtual double getFeatureEffect(unsigned int probeIx) = 0;

  /** 
   * @brief Get the estimated intensity for a chip. This is usually the data of
   * primary interest.
   * @param chipIx - Index of chip in experiment.
   * @return estimated intensity.
   */
  virtual double getTargetEffect(unsigned int chipIx) = 0; 

  /** 
   * @brief Does this quantification method supply residuals?
   * @return Return true if getResidual() works, false otherwise.
   */
  virtual bool haveResiduals() { return false; }

  /** 
   * @brief Get the residual (the intensity not explained by model). This is
   * often supplied as log() transformed.
   *
   * @param probeIx - Index of probe in probe sets.
   * @param chipIx - Chip index in experiment.
   * @return intensity not explained by model.
   */
  virtual double getResidual(unsigned int probeIx, unsigned int chipIx) = 0;

  /** 
   * @brief Get the estimated intensity for a chip. This is usually the data of
   * primary interest. Alias for getTargetEffect()
   * @param chipIx - Index of chip in experiment.
   * @return estimated intensity.
   */
  virtual double getSignalEstimate(unsigned int chipIx) = 0;

  /** 
   * @brief What is the name of the quantification method?
   * @return name of adjuster.
   */  
  virtual std::string getType() { return m_Type; }

  /** 
   * @brief Was this probe used for the estimate?
   * @param probeIx - index of probe in probe set.
   * @return true if feature (probe) was used, false otherwise.
   */
  virtual bool featureUsed(unsigned int probeIx) { return true; }
  
  /** 
   * @brief Get the Probe for a particular index in the probe set.
   * @param probeIx - index of probe in probe set.
   * @return - Probe pointer of interest.
   */
  virtual const Probe *getFeature(unsigned int probeIx) { return m_Probes[probeIx]; }

  /**
   * @brief What should the filename for summaries be called?
   */
  virtual const std::string &getSummarySuffix() {
    static std::string suffix = ".summary.txt";
    return suffix;
  }
  
  /**
   * @brief What should the feature response file be called?
   */
  virtual const std::string &getFeatureResponseSuffix() {
    static std::string suffix = ".feature-response.txt";
    return suffix;
  }

  /**
   * @brief What should the residual file be called?
   */
  virtual const std::string &getResidualSuffix() {
    static std::string suffix = ".residuals.txt";
    return suffix;
  }

   /** 
   * What version of the algorithm are we implementing?
   * @return version string. abstract method which specific
   * algorithm must implement. 
   */
  virtual inline std::string getVersion()=0;

protected :

  /// Type of method
  std::string m_Type; 
  /// Probes used, memory not owned here. 
  std::vector<Probe *> m_Probes;
};

#endif /* _QUANTMETHOD_H_ */
