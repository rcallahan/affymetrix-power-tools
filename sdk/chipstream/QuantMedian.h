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
 * @file   QuantMedian.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 10:30:02 2005
 * 
 * @brief Class for doing probe set quantification using median
 * as does MEDIAN.
 */
#ifndef _QUANTMEDIAN_H_
#define _QUANTMEDIAN_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/PmAdjuster.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantExprMethod.h"
//
#include <cmath>
#include <vector>
//

/// String describing iterplier summary method
#define QUANTMEDIANSTR "median"

#define QUANTMEDIAN_VERSION "1.0"

/**
 *  Class for doing probe set quantification using median
 * as does MEDIAN.
 */
class QuantMedian : public QuantExprMethod {

public:

  /** 
   * What version of the algorithm are we implementing?
   * @return version string.
   */
  std::string getVersion();

  /**
   * Basic constructor.
   */
  QuantMedian(bool doMean=false, bool attenuate=true, float l=0.005, float h=-1);

  /** 
   * @brief Does this quantification method estimate feature effects?
   * @return true or false if getFeatureEffects() supplies estimates.
   */
  bool haveFeatureEffects();

  /** 
   * @brief Does this quantification method supply residuals?
   * @return Return true if getResidual() works, false otherwise.
   */
  bool haveResiduals();

  /** 
   * @brief Get the feature effect of a particular probe.
   * @param probeIx - Index of probe in probe sets.
   * @return feature (probe) effects.
   */
  double getFeatureEffect(unsigned int probeIx);
  
  /** 
   * @brief Get the estimated intensity for a chip. This is usually the data of
   * primary interest.
   * @param chipIx - Index of chip in experiment.
   * @return estimated intensity.
   */
  double getTargetEffect(unsigned int chipIx);

  /** 
   * @brief Get the residual (the intensity not explained by
   * model). Supplied as log_2() transformed.
   * @param probeIx - Index of probe in probe sets.
   * @param chipIx - Chip index in experiment.
   * @return intensity not explained by model.
   */
  double getResidual(unsigned int probeIx, unsigned int chipIx);

  /** 
   * @brief Get the estimated intensity for a chip. This is usually the data of
   * primary interest. Alias for getTargetEffect()
   * @param chipIx - Index of chip in experiment.
   * @return estimated intensity.
   */
  double getSignalEstimate(unsigned int chipIx);

  /** 
   * @brief Clear out all the data.
   */
  void clear();

  /** 
   * @brief Set the PM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   * @param data - Data to be set.
   */
  void setPMDataAt(unsigned int probeIx, unsigned int chipIx, double data);

/** 
   * @brief Get the PM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   */
  double getPMDataAt(unsigned int probeIx, unsigned int chipIx);

  /** 
   * @brief Set the MM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   * @param data - Data to be set.
   */
  void setMMDataAt(unsigned int probeIx, unsigned int chipIx, double data);

  /** 
   * Set the number of probes and number of chips to be used for estimation.
   * @param numProbes - Probe count.
   * @param numChips - Microarray count.
   */
  void setBounds(unsigned int numProbes, unsigned int numChips);

  /** 
   * What scale are the target effects on? Are they linear (like
   * plier) or log2 (like median)
   * @return - What scale are the results provided in?
   */
  virtual Scale getScale();

  /** 
   * @brief Compute the median. If PM-MM is zero
   * threshold at a small positive number to avoid taking logs of
   * numbers that aren't positive.
   */
  void computeEstimate();

  /** 
   * @brief Set up the quantification method given all the data about the probe
   * set, chip layout and data.
   * 
   * @param psGroup - Probes to be used for final estimate.
   * @param layout - Chip layout annotation.
   * @param iMart - Raw data from chips.
   * @param iTrans - Transfomediantions to be applied to data before use.
   * @param pmAdjust - How to estimate background, or MM probe.
   * @return True if setup sucessful, false otherwise.
   */
  bool setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart,
             std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust);

  /** 
   * @brief Return number of features (probes).
   * @return feature count.
   */
  unsigned int getNumFeatures();

  /** 
   * @brief Return number of targets (experiments or chips).
   * @return target count
   */
  unsigned int getNumTargets();

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions();

  /** 
   * Fill in the infomediantion for Self documentation.
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

  /**
   * Custom configuration for this QuantMethod
   */
  virtual void setParameters(PsBoard &board) {
    // no-op here to override the error from parent class
  }

private:
  /// Number of chips in experiment
  unsigned int m_ChipCount;
  /// Number of probes (features) in probe sets
  unsigned int m_ProbeCount;
  /// Data matrix for PM probes.
  std::vector<std::vector<float> > m_PM;
  /// Data matrix for MM estimates.
  std::vector<std::vector<float> > m_MM;
  /// Chip effects (medians)
  std::vector<float> m_Summaries;
  /// Should we use the mean rather than the median
  bool m_DoMean;
  /// Should we use an attenuated mismatch, or simply cap the lower bounds when using non-PM-only adjuster
  bool m_Attenuate;
  /// Tunable parameter for handling background value in non-PM-only adjusters
  float m_L;
  /// Override use of computed H from data with fixed constant in background attenuation
  float m_H;
};

#endif /* _QUANTMEDIAN_H_ */
