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
 * @file   QuantRma.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 10:30:02 2005
 * 
 * @brief Class for doing probe set quantification using median polish
 * as does RMA.
 */
#ifndef _QUANTRMA_H_
#define _QUANTRMA_H_

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

/// Min value that RMA uses fo median polish. Must be positive as log() will be
/// taken.
#define RMA_MIN_VALUE 0.000000000001
/// String describing iterplier summary method
#define QUANTRMASTR "med-polish"
/**
 *  Class for doing probe set quantification using median polish
 * as does RMA.
 */
class QuantRma : public QuantExprMethod {

public:

  /** 
   * What version of the algorithm are we implementing?
   * @return version string.
   */
  std::string getVersion() {
    return "1.0";
  }

  /**
   * Basic constructor.
   */
  QuantRma(bool fixFeatureEffect=false, bool useInputModel=false, bool fitFeatureResponse=true, bool exponentiate=false, bool attenuate=true, float l=0.005, float h=-1) { 
    setupSelfDoc(*this);
    m_Type = std::string(QUANTRMASTR);
	m_FixFeatureEffect = fixFeatureEffect;
    m_UseInputModel=useInputModel;
    m_FitFeatureResponse=fitFeatureResponse;
    m_Exponentiate = exponentiate;
    m_Attenuate = attenuate;
    m_L = l;
    m_H = h;
    setOptValue("FixFeatureEffect", m_FixFeatureEffect);
    setOptValue("UseInputModel", m_UseInputModel);
    setOptValue("FitFeatureResponse", m_FitFeatureResponse);
    setOptValue("expon", m_Exponentiate);
    setOptValue("attenuate", m_Attenuate);
    setOptValue("l", ToStr(m_L));
    setOptValue("h", ToStr(m_H));
    m_UsePrecompEffects = false;
    m_UseInputModel = false;
    m_FitFeatureResponse = true;
    m_FeaturePriorEffect = NULL;
    m_ChipCount = 0;
    m_ProbeCount = 0;
  }

  virtual ~QuantRma()
  {
	if (m_FeaturePriorEffect != NULL) {delete[] m_FeaturePriorEffect; m_FeaturePriorEffect = NULL;}
  }


  /** 
   * @brief Does this quantification method estimate feature effects?
   * @return true or false if getFeatureEffects() supplies estimates.
   */
  inline bool haveFeatureEffects() { return true;}

  /** 
   * @brief Does this quantification method supply residuals?
   * @return Return true if getResidual() works, false otherwise.
   */
  inline bool haveResiduals() { return true;}

  /** 
   * @brief Get the feature effect of a particular probe.
   * @param probeIx - Index of probe in probe sets.
   * @return feature (probe) effects.
   */
  inline double getFeatureEffect(unsigned int probeIx) { 
    assert(probeIx < m_ProbeCount);
    return m_ProbeEffects[probeIx];
  }
 
  /**
   * @brief Specify precomputed feature effects.
   *
   * @param featureEffects - an array of feature effects.
   */
  void setFeaturePriorEffects(double *featureEffects, int iSize){
    Verbose::out(3, "QuantRma::setFeaturePriorEffects()");
	if (m_FeaturePriorEffect != NULL) {delete[] m_FeaturePriorEffect;}
	m_FeaturePriorEffect = new double[iSize];
	for (int iIndex = 0; (iIndex < iSize); iIndex++)
	{
		m_FeaturePriorEffect[iIndex] = featureEffects[iIndex];
	}
    setUsePrecompFeatureEffects(true);
  }

  /**
   * @brief Use this to toggle whether or not QuantRma should be using
   * precomputed feature effects (aka feature response, sometimes know incorrectly as 'affinity').
   * To mimic what occurs in pliers we do not dynamically fit the model using the precomputed feature
   * effects as input.  Essentially we do just "one pass" of the med-polish algorithm using the 
   * precomputed feature effects and output the resulting chip effects. When this note was written
   * the above procedure is the default behaviour when a precomputed feature effects file is read. 
   *
   * @param value - true if using precomputed feature effects false otherwise.
   */
  void setUsePrecompFeatureEffects(bool value) {
      m_UsePrecompEffects = value;
      if(m_UsePrecompEffects) { 
        if(m_FeaturePriorEffect == NULL) {
          Err::errAbort("QuantRma::setUsePrecompFeatureEffects() - Can't use precomputed feature effects without supplying them first.");
        }
        setFitFeatureResponse(false);
        setUseInputModel(true);  
      }
      else {
        setFitFeatureResponse(true);
        setUseInputModel(false); 
      }
  }
  /** 
   * @brief Set the boolean value which determines whether or not we dynaimically fit 
   * the model to the feature effect values. 
   */
  void setFitFeatureResponse(bool value){
    m_FitFeatureResponse = value;
    setOptValue("FitFeatureResponse", value);
  }
   
  /** 
   * @brief Get the boolean value which determines whether or not we dynaimically fit 
   * the model to the feature effect values. 
   */
  bool getFitFeatureResponse(){
   return  m_FitFeatureResponse;
  }
   
  /** 
   * @brief Set the boolean value which determines whether or not we use precomputed 
   * feature effect values. 
   */
  void setUseInputModel(bool value){
    m_UseInputModel=value;
    setOptValue("UseInputModel", value);
  }
   
  /** 
   * @brief Get the boolean value which determines whether or not we use  precomputed
   * feature effect values. 
   */
  bool getUseInputModel (){
   return  m_UseInputModel;
  }


  /** 
   * @brief Get the estimated intensity for a chip. This is usually the data of
   * primary interest.
   * @param chipIx - Index of chip in experiment.
   * @return estimated intensity.
   */
  inline double getTargetEffect(unsigned int chipIx) {
    assert(chipIx < m_ChipCount);
    if(m_Exponentiate)
        return pow((double)2.0,(double)m_ChipEffects[chipIx]);
    return m_ChipEffects[chipIx];
  }

  /** 
   * @brief Get the residual (the intensity not explained by
   * model). Supplied as log_2() transformed.
   * @param probeIx - Index of probe in probe sets.
   * @param chipIx - Chip index in experiment.
   * @return intensity not explained by model.
   */
  inline double getResidual(unsigned int probeIx, unsigned int chipIx) {
    assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
    return m_PM[chipIx][probeIx];
  }

  /** 
   * @brief Get the estimated intensity for a chip. This is usually the data of
   * primary interest. Alias for getTargetEffect()
   * @param chipIx - Index of chip in experiment.
   * @return estimated intensity.
   */
  inline double getSignalEstimate(unsigned int chipIx) {
    return getTargetEffect(chipIx);
  }

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
  inline void setPMDataAt(unsigned int probeIx, unsigned int chipIx, double data) {
    assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
    m_PM[chipIx][probeIx] = (float)data;
  }

/** 
   * @brief Get the PM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   */
  inline double getPMDataAt(unsigned int probeIx, unsigned int chipIx) {
    assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
	return (double) m_PM[chipIx][probeIx];
  }

  /** 
   * @brief Set the MM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   * @param data - Data to be set.
   */
  inline void setMMDataAt(unsigned int probeIx, unsigned int chipIx, double data) {
    assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
    m_MM[chipIx][probeIx] = (float)data;
  }

  /** 
   * Set the number of probes and number of chips to be used for estimation.
   * @param numProbes - Probe count.
   * @param numChips - Microarray count.
   */
  void setBounds(unsigned int numProbes, unsigned int numChips);

  /** 
   * What scale are the target effects on? Are they linear (like
   * plier) or log2 (like rma)
   * @return - What scale are the results provided in?
   */
  virtual Scale getScale() { return m_Exponentiate ? Linear : Log2; }

  /** 
   * @brief Do the heavy lifting of median polish. If PM-MM is zero
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
   * @param iTrans - Transformations to be applied to data before use.
   * @param pmAdjust - How to estimate background, or MM probe.
   * @return True if setup sucessful, false otherwise.
   */
  bool setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
             std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust);

  /** 
   * @brief Return number of features (probes).
   * @return feature count.
   */
  inline unsigned int getNumFeatures()  { return m_ProbeCount; }

  /** 
   * @brief Return number of targets (experiments or chips).
   * @return target count
   */
  inline unsigned int getNumTargets() { return m_ChipCount; }

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> opts;
    SelfDoc::Opt fixFeatureEffect = {"FixFeatureEffect", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                          "Force the calculation of target estimate using calculated feature effects."};
    opts.push_back(fixFeatureEffect);
    SelfDoc::Opt useInputModel = {"UseInputModel", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", 
                                       "Use provided values as the initial model of Feature Responses."};
    opts.push_back(useInputModel);     
    SelfDoc::Opt fitFeatureResponse = {"FitFeatureResponse", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",  
                                            "Fit Feature Response dynamically or don't update from initial values."};
    opts.push_back(fitFeatureResponse);
    SelfDoc::Opt expon = {"expon", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                          "Convert back from log space by exponentiating the target estimate (i.e. 2^x)."};
    opts.push_back(expon);
    SelfDoc::Opt attenuate = {"attenuate", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                              "Indicate whether or not to attenuate mismatch value when using non-PM-only adjuster."};
    opts.push_back(attenuate);
    SelfDoc::Opt l = {"l", SelfDoc::Opt::Double, "0.005", "0.005", "0", "1",
                      "Tunable parameter for attenuating mismatch value for non-PM-only adjusters."};
    opts.push_back(l);
    SelfDoc::Opt h = {"h", SelfDoc::Opt::Double, "-1", "-1", "-1", "NA",
                      "Used fixed constant to attenuate mismatch. Is set to 4*PM*MM*L when default is supplied."};
    opts.push_back(h);
    
    return opts;
  }

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName(QUANTRMASTR);
    doc.setDocDescription("Performs a median polish to estimate target and probe effects. Resulting summaries are in log2 space by default. Used in summary step of RMA as described in Irizarry et al 2003.");
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
    SelfDoc doc = explainSelf();
    bool fixFeatureEffect = false;
    bool useInputModel = false;
    bool fitFeatureResponse = true;
    bool expon = false;
    bool attenuate = true;
    float l = 0.005f;
    float h = -1.0f;
    fillInValue(fixFeatureEffect, "FixFeatureEffect", param, doc);
    fillInValue(useInputModel, "UseInputModel", param, doc);
    fillInValue(fitFeatureResponse, "FitFeatureResponse", param, doc);
    fillInValue(expon, "expon", param, doc);
    fillInValue(attenuate, "attenuate", param, doc);
    fillInValue(l, "l", param, doc);
    fillInValue(h, "h", param, doc);
    QuantRma *rma = new QuantRma(fixFeatureEffect,useInputModel,fitFeatureResponse,expon,attenuate,l,h);
    return rma;
  }

  /**
   * Custom configuration for this QuantMethod
   */
  virtual void setParameters(PsBoard &board);

private:
  /// Number of chips in experiment
  unsigned int m_ChipCount;
  /// Number of probes (features) in probe sets
  unsigned int m_ProbeCount;
  /// Data matrix for PM probes.
  std::vector<std::vector<float> > m_PM;
  /// Data matrix for MM estimates.
  std::vector<std::vector<float> > m_MM;
  /// Effects of columns (probes)
  std::vector<float> m_ProbeEffects;
  /// Effects of rows (chips)
  std::vector<float> m_ChipEffects;
  /// Should we force the calculation of the target effects using calculated feature effects?
  bool m_FixFeatureEffect;
  /// Should we exponentiate the target effects? i.e. 2^X to compensate for taking log_2
  bool m_Exponentiate;
  /// Should we use an attenuated mismatch, or simply cap the lower bounds when using non-PM-only adjuster
  bool m_Attenuate;
  /// Tunable parameter for handling background value in non-PM-only adjusters
  float m_L;
  /// Override use of computed H from data with fixed constant in background attenuation
  float m_H;
  /// Precomputed feature effects. 
  double* m_FeaturePriorEffect;
  /// Should we use the precomputed probe/feature effects in m_FeaturePriorEffect?
  bool m_UsePrecompEffects;
  // Do we dynamically fit the model to the feature effects (either default values or values provided  by user).   
  bool m_FitFeatureResponse;
  // Do we use values provided by the user    
  bool m_UseInputModel;

};

#endif /* _QUANTRMA_H_ */
