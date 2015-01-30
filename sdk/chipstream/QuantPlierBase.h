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
 * @file   QuantPlierBase.h
 * @author Pete Klosterman
 * @date   Wed Feb 22 10:55:22 2006
 * 
 * @brief  Base class for doing probe set quantification using either
 * the PLIER (Probe Logarithmic Error Intensity Estimate) method or
 * the ITER-PLIER method, which iteratively calls PLIER with the
 * probes that best correlate with signal estimate.
 */
#ifndef _QUANTPLIERBASE_H_
#define _QUANTPLIERBASE_H_

//
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantPlierParams.h"
//
#include "plier/affyplier.h"
#include "util/Convert.h"
#include "util/Err.h"
//
#include <cstring>
#include <map>
#include <string>
//

/**
 *  Base class for doing probe set quantification using either
 * the PLIER (Probe Logarithmic Error Intensity Estimate) method or
 * the ITER-PLIER method, which iteratively calls PLIER with the
 * probes that best correlate with signal estimate.
 */
class QuantPlierBase : public QuantExprMethod {

public:

  /**
   * Constructor. 
   */
  QuantPlierBase(QuantPlierParams &plierParams);

  /**
   * Destructor.
   */
  virtual ~QuantPlierBase();

  /** 
   * @brief Set up the quantification method given all the data about the probe
   * set, chip layout and data.
   * 
   * @param psGroup - Probes to be used for final estimate.
   * @param layout - Chip layout annotation.
   * @param iMart - Raw data from chips.
   * @param iTrans - Transformations to be applied to data before use.
   * @param pmAdjust - How to estimate background, or MM probe.
   * @return bool - True if setup sucessful, false otherwise.
   */
  bool setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
             std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust);

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
   * @brief Do the heavy lifting of plier or iter-plier.
   */
  virtual void computeEstimate() = 0;

  /** 
   * @brief Return number of features (probes).
   * @return feature count.
   */
  virtual unsigned int getNumFeatures() = 0;

  /** 
   * @brief Return number of targets (experiments or chips).
   * @return target count
   */
  unsigned int getNumTargets();

  /** 
   * @brief Specify precomputed feature effects to be used as map from
   * probe id index to the feature effect.  This make plier much
   * faster to do as doesn't have to iterate through anything, but
   * have to make sure that the feature effects being used match the
   * analysis that is being done. For example the feature effects
   * might be different using pm-only, pm-gcbg and pm-mm adjustments.
   * 
   * Note: this is only used by plier, not by iter-plier.
   * 
   * @param featureEffects - mapping from probe idex to feature effects.
   */
  virtual void setFeaturePriorEffects(double* featureEffects, int iSize) = 0;

  /** 
   * @brief Fill the plier parameters object based on the user-defined
   * key/value pairs.
   * The default values used for the plier parameters are from the
   * inputs.txt file in the sdk example.
   * 
   * @param param - Map of key/value pairs to initialize the object.
   * @param plierParams - Plier parameters object to fill.
   */
  static void fillPlierParams(std::map<std::string,std::string> &param, 
                              QuantPlierParams &plierParams,
                              SelfDoc &doc);

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions();

protected :

  /** 
   * @brief Allocate enough memory for at least maxChips and
   * maxProbes.
   * 
   * @param maxChips - Number of chips that can be computed.
   * @param maxProbes - Number of probes or features that can be computed.
   */
  virtual void allocMemory(unsigned int maxChips, unsigned int maxProbes) = 0;

  /** 
   * @brief Free up all the memory that has been allocated.
   */
  virtual void freeMemory() = 0;

  /** 
   * @brief Use this to toggle whether or not QuantPlier should be using
   * precomputed feature effects (aka feature response, sometimes know incorrectly as 'affinity').
   * Much faster to do as plier doesn't have to iterate through anything, but
   * have to make sure that the feature effects being used match the analysis that
   * is being done. For example the feature effects might be different using
   * pm-only, pm-gcbg and pm-mm adjustments.
   * 
   * @param value - true if using precomputed feature effects false otherwise.
   */
  void setUsePrecompFeatureEffects(bool value);

  /**
   * Custom configuration for this QuantMethod
   */
  virtual void setParameters(PsBoard &board);


  /// Instance of plier to be used.
  caffyplier m_Plier;
  /// PM data matrix (double** as  that is what plier uses.)
  double **m_PM;
  /// MM data matrix (double** as  that is what plier uses.)
  double **m_MM;
  /// Residual matrix (double** as  that is what plier uses.) log_2 transformed.
  double **m_Residuals;
  /// Effects attributed to probe (feature) and chips by model
  double *m_ProbeEffects, *m_ChipEffects;
  /// Number of chips and counts in this probe set analysis.
  unsigned int m_ChipCount, m_ProbeCount;
  /// Map of probe ids to feature effects.
  double* m_FeaturePriorEffect;
  /// Maximum number of chips and probes that can be done without
  /// reallocating more memory.
  unsigned int m_MaxChips, m_MaxProbes;

private :

  /** 
   * @brief Plier has a lot of parameters. Lets set as many as possible in
   * one place... Got default parameters from the inputs.txt file
   * in the sdk example.
   * @param plier - plier class to set parameters for.
   */
  void setPlierParams(iaffyplier &plier, QuantPlierParams &plierParams);

  /// Should we use the precomputed probe effects in m_ProbePriorEffect?
  bool m_UsePrecompEffects;
};

#endif /* _QUANTPLIERBASE_H_ */
