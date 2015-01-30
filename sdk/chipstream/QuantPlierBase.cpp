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
 * @file   QuantPlierBase.cpp
 * @author Pete Klosterman
 * @date   Wed Feb 22 10:55:22 2006
 * 
 * @brief  Base class for doing probe set quantification using either
 * the PLIER (Probe Logarithmic Error Intensity Estimate) method or
 * the ITER-PLIER method, which iteratively calls PLIER with the
 * probes that best correlate with signal estimate.
 */

//
#include "chipstream/QuantPlierBase.h"
//
#include "chipstream/QuantMethodFactory.h"

using namespace std;

/** 
 * @brief Does this quantification method estimate feature effects?
 * @return true or false if getFeatureEffects() supplies estimates.
 */
bool QuantPlierBase::haveFeatureEffects() {
  return true;
}

/** 
 * @brief Does this quantification method supply residuals?
 * @return Return true if getResidual() works, false otherwise.
 */
bool QuantPlierBase::haveResiduals() {
  return true;
}

/** 
 * @brief Get the feature effect of a particular probe.
 * @param probeIx - Index of probe in probe sets.
 * @return feature (probe) effects.
 */
double QuantPlierBase::getFeatureEffect(unsigned int probeIx) { 
  assert(probeIx < m_ProbeCount);
  return m_ProbeEffects[probeIx];
}
  
/** 
 * @brief Get the estimated intensity for a chip. This is usually the data of
 * primary interest.
 * @param chipIx - Index of chip in experiment.
 * @return estimated intensity.
 */
double QuantPlierBase::getTargetEffect(unsigned int chipIx) {
  assert(chipIx < m_ChipCount);
  return m_ChipEffects[chipIx];
}

/** 
 * @brief Get the residual (the intensity not explained by
 * model). Supplied as log_2() transformed.
 * @param probeIx - Index of probe in probe sets.
 * @param chipIx - Chip index in experiment.
 * @return intensity not explained by model.
 */
double QuantPlierBase::getResidual(unsigned int probeIx, unsigned int chipIx) {
  assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
  return m_Residuals[chipIx][probeIx];
}

/** 
 * @brief Get the estimated intensity for a chip. This is usually the data of
 * primary interest. Alias for getTargetEffect()
 * @param chipIx - Index of chip in experiment.
 * @return estimated intensity.
 */
double QuantPlierBase::getSignalEstimate(unsigned int chipIx) {
  return getTargetEffect(chipIx);
}


/** 
 * @brief Set the PM data at probe set probe index and chip index.
 * 
 * @param probeIx - Index of probe in probe set.
 * @param chipIx - Index of chip.
 * @param data - Data to be set.
 */
void QuantPlierBase::setPMDataAt(unsigned int probeIx, unsigned int chipIx, double data) {
  assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
  m_PM[chipIx][probeIx] = data;
}

/** 
 * @brief Get the PM data at probe set probe index and chip index.
 * 
 * @param probeIx - Index of probe in probe set.
 * @param chipIx - Index of chip.
 */
double QuantPlierBase::getPMDataAt(unsigned int probeIx, unsigned int chipIx) {
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
void QuantPlierBase::setMMDataAt(unsigned int probeIx, unsigned int chipIx, double data) {
  assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
  m_MM[chipIx][probeIx] = data;
}

/** 
 * @brief Return number of targets (experiments or chips).
 * @return target count
 */
unsigned int QuantPlierBase::getNumTargets() {
  return m_ChipCount;
}

/** 
 * @brief Fill the plier parameters object based on the user-defined
 * key/value pairs.
 * The default values used for the plier parameters are from the
 * inputs.txt file in the sdk example.
 * 
 * @param param - Map of key/value pairs to initialize the object.
 * @param plierParams - Plier parameters object to fill.
 */
void QuantPlierBase::fillPlierParams(std::map<std::string,std::string> &param, 
                                     QuantPlierParams &plierParams,
                                     SelfDoc &doc) {
  // unused
  // std::map<std::string,std::string>::const_iterator iter;

  float attenuation = 0;
  fillInValue(attenuation, "atten", param, doc);
  plierParams.setSeaAttenuation(attenuation);

  double initAugmentation = 0;
  fillInValue(initAugmentation, "InitAugmentation", param, doc);
  plierParams.setInitAugmentation(initAugmentation);

  double initDefaultFeatureResponse = 0;
  fillInValue(initDefaultFeatureResponse, "InitDefaultFeatureResponse", param, doc);
  plierParams.setInitDefaultFeatureResponse(initDefaultFeatureResponse);

  double initDefaultTargetResponse = 0;
  fillInValue(initDefaultTargetResponse, "InitDefaultTargetResponse", param, doc);
  plierParams.setInitDefaultTargetResponse(initDefaultTargetResponse);
    
  double seaOptConvergence = 0;
  fillInValue(seaOptConvergence, "SeaOptConvergence", param, doc);
  plierParams.setSeaOptConvergence(seaOptConvergence);

  int seaOptIteration = 0;
  fillInValue(seaOptIteration, "SeaOptIteration", param, doc);
  plierParams.setSeaOptIteration((long)seaOptIteration);

  float gmCutoff = 0;
  fillInValue(gmCutoff, "PlierGmCutoff", param, doc);
  plierParams.setPlierGmCutoff(gmCutoff);

  float differentialFeaturePenalty = 0;
  fillInValue(differentialFeaturePenalty, "PlierDifferentialFeaturePenalty", param, doc);
  plierParams.setPlierDifferentialFeaturePenalty(differentialFeaturePenalty);

  float differentialTargetPenalty = 0;
  fillInValue(differentialTargetPenalty, "PlierDifferentialTargetPenalty", param, doc);
  plierParams.setPlierDifferentialTargetPenalty(differentialTargetPenalty);

  bool useMMLikelihood = false;
  fillInValue(useMMLikelihood, "PlierUseMMLikelihood", param, doc);
  plierParams.setPlierUseMMLikelihood(useMMLikelihood);

  bool useInputModel = false;
  fillInValue(useInputModel, "PlierUseInputModel", param, doc);
  plierParams.setPlierUseInputModel(useInputModel);

  bool fitFeatureResponse = false;
  fillInValue(fitFeatureResponse, "PlierFitFeatureResponse", param, doc);
  plierParams.setPlierFitFeatureResponse(fitFeatureResponse);

  double optConvergence = 0;
  fillInValue(optConvergence, "PlierOptConvergence", param, doc);
  plierParams.setPlierOptConvergence(optConvergence);

  int optIteration = 0;
  fillInValue(optIteration, "PlierOptIteration", param, doc);
  plierParams.setPlierOptIteration((long)optIteration);

  double optDropMax = 0;
  fillInValue(optDropMax, "PlierOptDropMax", param, doc);
  plierParams.setPlierOptDropMax(optDropMax);

  double optLambdaLimit = 0;
  fillInValue(optLambdaLimit, "PlierOptLambdaLimit", param, doc);
  plierParams.setPlierOptLambdaLimit(optLambdaLimit);

  int optimization = 0;
  fillInValue(optimization, "optmethod", param, doc);
  plierParams.setPlierOptOptimizationMethod((long)optimization);
    
  int optBalanceMethod = 0;
  fillInValue(optBalanceMethod, "PlierOptBalanceMethod", param, doc);
  plierParams.setPlierOptBalanceMethod((long)optBalanceMethod);

  bool optFixPrecomputed = false;
  fillInValue(optFixPrecomputed, "FixPrecomputed", param, doc);
  plierParams.setFixPrecomputed(optFixPrecomputed);

  double optNumericalTolerance = 0;
  fillInValue(optNumericalTolerance, "NumericalTolerance", param, doc);
  plierParams.setNumericalTolerance(optNumericalTolerance);

  double optSafetyZero = 0;
  fillInValue(optSafetyZero, "SafetyZero", param, doc);
  plierParams.setSafetyZero(optSafetyZero);

  bool optFixFeatureEffect = false;
  fillInValue(optFixFeatureEffect, "FixFeatureEffect", param, doc);
  plierParams.setFixFeatureEffect(optFixFeatureEffect);
}

/** 
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> QuantPlierBase::getDefaultDocOptions() { 
  std::vector<SelfDoc::Opt> opts;
  SelfDoc::Opt optmethod = {"optmethod", SelfDoc::Opt::Integer, "0", "0", "0", "1",
                            "Optimization method to use for plier 1 for SEA (Simplified Expression Analysis), 0 for full Plier optimization."};
  opts.push_back(optmethod);
  SelfDoc::Opt atten = {"atten", SelfDoc::Opt::Double, "0.005", "0.005", "0", "1",
                        "Attenuation to use for background with SEA."};
  opts.push_back(atten);
  SelfDoc::Opt initAugmentation = {"InitAugmentation", SelfDoc::Opt::Double, "0.1", "0.1", "0", "NA", 
                                   "Positive number added to all values to void zero values in input data."};
  opts.push_back(initAugmentation);
  SelfDoc::Opt initDefaultFeatureResponse = {"InitDefaultFeatureResponse", SelfDoc::Opt::Double, "1.0", "1.0", "0", "NA",
                                             "Default FeatureResponse if not supplied."};
  opts.push_back(initDefaultFeatureResponse);
  SelfDoc::Opt initDefaultTargetResponse = {"InitDefaultTargetResponse", SelfDoc::Opt::Double, "1.0", "1.0", "0", "NA",
                                            "Default TargetResponse if not supplied."};
  opts.push_back(initDefaultTargetResponse);
  SelfDoc::Opt seaOptConvergence = {"SeaOptConvergence", SelfDoc::Opt::Double, "0.000001", "0.000001", "0", "NA",
                                    "Change in log-value at which to stop."};
  opts.push_back(seaOptConvergence);
  SelfDoc::Opt seaOptIteration = {"SeaOptIteration", SelfDoc::Opt::Integer, "2000", "2000", "1", "NA",
                                  "Max number of SEA iteration to avoid infinite loops in SEA."};
  opts.push_back(seaOptIteration);
  SelfDoc::Opt plierGmCutoff = {"PlierGmCutoff", SelfDoc::Opt::Double, "0.15", "0.15", "0", "NA",
                                "Controls discounting outliers, larger values indicate that fewer outliers are expected."};
  opts.push_back(plierGmCutoff);
  SelfDoc::Opt plierDifferentialFeaturePenalty = {"PlierDifferentialFeaturePenalty", SelfDoc::Opt::Double, "0.001", "0.001", "0", "NA",
                                                  "Bayes penalty for peculiar features."};
  opts.push_back(plierDifferentialFeaturePenalty);
  SelfDoc::Opt plierDifferentialTargetPenalty = {"PlierDifferentialTargetPenalty", SelfDoc::Opt::Double, "0.000001", "0.000001", "0", "NA",
                                                 "Bayes penalty for really peculiar TargetResponses."};
  opts.push_back(plierDifferentialTargetPenalty);
  SelfDoc::Opt plierUseMMLikelihood = {"PlierUseMMLikelihood", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                                       "Use mm or background based likelihood."};
  opts.push_back(plierUseMMLikelihood);
  SelfDoc::Opt plierUseInputModel = {"PlierUseInputModel", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", 
                                     "Use provided values as the initial model of Feature Responses."};
  opts.push_back(plierUseInputModel);     
  SelfDoc::Opt plierFitFeatureResponse = {"PlierFitFeatureResponse", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",  
                                          "Fit Feature Response dynamically or don't update from initial values."};
  opts.push_back(plierFitFeatureResponse);
  SelfDoc::Opt plierOptConvergence = {"PlierOptConvergence", SelfDoc::Opt::Double, "0.000001", "0.000001", "0", "NA", 
                                      "Value of PLIER convergence, change in likelihood."};
  opts.push_back(plierOptConvergence);
  SelfDoc::Opt plierOptIteration = {"PlierOptIteration", SelfDoc::Opt::Integer, "3000", "3000", "0", "NA",
                                    "Max number of PLIER iteration to avoid infinite loops in PLIER."};
  opts.push_back(plierOptIteration);
  SelfDoc::Opt plierOptDropMax = {"PlierOptDropMax", SelfDoc::Opt::Double, "3.0", "3.0", "1.0", "NA", 
                                  "Used during descent to avoid negative or zero values."};
  opts.push_back(plierOptDropMax);
  SelfDoc::Opt plierOptLambdaLimit = {"PlierOptLambdaLimit", SelfDoc::Opt::Double, "0.01", "0.01", "0", "1", 
                                      "Minimum step multiplier in method."};
  opts.push_back(plierOptLambdaLimit);
  SelfDoc::Opt plierOptOptimizationMethod = {"PlierOptOptimizationMethod", SelfDoc::Opt::Integer, "0", "0", "0", "1",
                                             "Optimization method to use for plier 1 for SEA (Simplified Expression Analysis), 0 for full Plier optimization."}; 
  opts.push_back(plierOptOptimizationMethod); 
  SelfDoc::Opt plierOptBalanceMethod = {"PlierOptBalanceMethod", SelfDoc::Opt::Integer, "0", "0", "0", "5",
                                        "Identifiability method: 0 (sum), 1 (median), and 2 (SEA signal)."};
  opts.push_back(plierOptBalanceMethod);
  SelfDoc::Opt fixPrecomputed = {"FixPrecomputed", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                                 "Recompute signal after computing feature response. Set to true to get more consistent results with those generated from precomputed feature responses."};
  opts.push_back(fixPrecomputed);
  SelfDoc::Opt numericalTolerance = {"NumericalTolerance", SelfDoc::Opt::Double, "0.1", "0.1", "NA", "NA",
                                     "Set how strict the results should agree when using FixPrecomputed."};
  opts.push_back(numericalTolerance);
  SelfDoc::Opt safetyZero = {"SafetyZero", SelfDoc::Opt::Double, "0.000001", "0.000001", "NA", "NA",
                             "Set a small positive number to use in place of zero."};
  opts.push_back(safetyZero);
  SelfDoc::Opt FixFeatureEffect = {"FixFeatureEffect", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                                   "Fix the feature effect calculation to be more exect. Break when limit reached in plieralg::FitAdditiveModel()."};
  opts.push_back(FixFeatureEffect);
  return opts;
}


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
void QuantPlierBase::setUsePrecompFeatureEffects(bool value) {
  m_UsePrecompEffects = value;
  if(m_UsePrecompEffects) {
    if(m_FeaturePriorEffect == NULL) {
      Err::errAbort("QuantPlier::setUsePrecompFeatureEffects() - Can't use precomputed feature effects without supplying them first.");
    }
    m_Plier.setPlierFitFeatureResponse(false);
    m_Plier.setPlierUseInputModel(true);
    setOptValue("PlierFitFeatureResponse", false);
    setOptValue("PlierUseInputModel", true);
  }
  else {
    m_Plier.setPlierFitFeatureResponse(true);
    m_Plier.setPlierUseInputModel(false);
    setOptValue("PlierFitFeatureResponse", true);
    setOptValue("PlierUseInputModel", false);
  }
}

/** 
 * Plier has a lot of parameters. Lets set as many as possible in
 * one place.  The plierParams object is initialized in the
 * QuantPlier newObject() method.
 * @param plier - plier class to set parameters for.
 * @param plierParams - plier parameters object.
 */
void QuantPlierBase::setPlierParams(iaffyplier &plier, QuantPlierParams &plierParams) {
  plier.setInitAugmentation(plierParams.getInitAugmentation());
  plier.setInitDefaultFeatureResponse(plierParams.getInitDefaultFeatureResponse());
  plier.setInitDefaultTargetResponse(plierParams.getInitDefaultTargetResponse());
  
  plier.setSeaAttenuation(plierParams.getSeaAttenuation());
  plier.setSeaOptConvergence(plierParams.getSeaOptConvergence());
  plier.setSeaOptIteration(plierParams.getSeaOptIteration());

  plier.setPlierGmCutoff(plierParams.getPlierGmCutoff());
  plier.setPlierDifferentialFeaturePenalty(plierParams.getPlierDifferentialFeaturePenalty());

  plier.setPlierUseMMLikelihood(plierParams.getPlierUseMMLikelihood());
  plier.setPlierUseInputModel(plierParams.getPlierUseInputModel());
  plier.setPlierFitFeatureResponse(plierParams.getPlierFitFeatureResponse());

  plier.setPlierOptConvergence(plierParams.getPlierOptConvergence());
  plier.setPlierOptIteration(plierParams.getPlierOptIteration());
  plier.setPlierOptDropMax(plierParams.getPlierOptDropMax());
  plier.setPlierOptLambdaLimit(plierParams.getPlierOptLambdaLimit());
  plier.setPlierOptOptimizationMethod(plierParams.getPlierOptOptimizationMethod());
  plier.setPlierOptBalanceMethod(plierParams.getPlierOptBalanceMethod());

  plier.setNumericalTolerance(plierParams.getNumericalTolerance());
  plier.setSafetyZero(plierParams.getSafetyZero());
  plier.setFixPrecomputed(plierParams.getFixPrecomputed());
  plier.setFixFeatureEffect(plierParams.getFixFeatureEffect());
}

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
bool QuantPlierBase::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
                           std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust) {
  unsigned int chipIx = 0, probeIx = 0, psIx, atomIx = 0;
  unsigned int chipCount = 0, atomPmCount = 0;
  float intensity = 0.0;
  assert(psGroup.probeSets.size() > 0);

  chipCount = iMart.getCelFileCount();
  atomPmCount = psGroup.countPmProbes();

  if(atomPmCount == 0)
    return false;

  /* Prepare the quantification method for this many atoms and chips. */
  setBounds(atomPmCount,chipCount);

  atomPmCount = 0;
  m_Probes.clear();

  /* Loop through and fill in the data as coming from the intensity mart. */
  for(psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
    const ProbeSet *ps = psGroup.probeSets[psIx];
    if(ps == NULL) {
      Verbose::out(1, "No probeset for index: " + ToStr(psIx));
      continue;
    }
    for(atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom &atom = *(ps->atoms[atomIx]);
      unsigned int channelIx = atom.getChannelCode();
      for(probeIx = 0; probeIx < atom.probes.size(); probeIx++) {
        Probe *p = atom.probes[probeIx];
        if(p->type == Probe::PMST || p->type == Probe::PMAT) {
          unsigned int probeIndex = p->id;
          /* If we are using precomputed probe effects fill in each probe here. */
          if(m_UsePrecompEffects) {
	    if (m_FeaturePriorEffect == NULL) {
	      Err::errAbort("No feature effects available.");
	    }
	    else {
	      /// Expected feature effect is 1.0 +/- some minor adjustment. 
	      try {
		m_ProbeEffects[atomPmCount] = m_FeaturePriorEffect[p->getApid()];
	      } 
	      catch(...) {
		Err::errAbort("Feature effect is out of bounds for probe "+ ToStr(p->id + 1));
	      }
	      if (m_ProbeEffects[atomPmCount] <= 0) {
		Err::errAbort("The feature effect values must be positive. probe_id may be missing from the input feature effects.  probeset: "+ ToStr(ps->name) +" probe_id: " + ToStr(p->id + 1) + ", FeatureEffect = " + ToStr(m_ProbeEffects[atomPmCount]));
	      } // Safety against log zero.
			  }
          }
          /* Add probe to our probes and fill in intensity data for each chip. */
          m_Probes.push_back(p);
          for(chipIx = 0; chipIx < chipCount; chipIx++) {
            /* PM transformation. */
            intensity = transformPrimaryData(probeIndex, chipIx, iMart, iTrans, channelIx);
            /* MM transformation. */
            float mmIntensity = 0;
            pmAdjust.pmAdjustment(probeIndex, chipIx, iMart, iTrans, intensity, mmIntensity);
            setPMDataAt(atomPmCount, chipIx, intensity);
            setMMDataAt(atomPmCount, chipIx, mmIntensity);
          }
          atomPmCount++;
        }
      }
    } /* atoms. */
  } /* probe sets. */
  return true;
}

/** Constructor. */
QuantPlierBase::QuantPlierBase(QuantPlierParams &plierParams) {
	m_FeaturePriorEffect = NULL;
  /* Note: precomputed effects are not to used by iter-plier. */
  m_UsePrecompEffects = false;
  //aw: this needs to be handled by the child as self doc is not setup yet
  //setUsePrecompFeatureEffects(m_UsePrecompEffects);
  setPlierParams(m_Plier, plierParams);
  m_MaxProbes = 0;
  m_MaxChips = 0;
  m_ChipCount = 0;
  m_ProbeCount = 0;
  m_PM = NULL;
  m_MM = NULL;
  m_Residuals = NULL;
  m_ProbeEffects = NULL;
  m_ChipEffects = NULL;
}

/** 
 * Destructor
 */
QuantPlierBase::~QuantPlierBase() 
{
	//Verbose::out(1, "~QuantPlierBase() is being called.");
	if (m_FeaturePriorEffect != NULL)
	{
		//Verbose::out(1, "m_FeaturePriorEffect is being destroyed.");
		delete[] m_FeaturePriorEffect;
		m_FeaturePriorEffect = NULL;
	}
}

/** 
 * Set the number of probes and number of chips to be used for estimation.
 * @param numProbes - Probe count.
 * @param numChips - Microarray count.
 */
void QuantPlierBase::setBounds(unsigned int numProbes, unsigned int numChips) {
  if ((numProbes > m_MaxProbes) || (numChips > m_MaxChips)) {
    allocMemory(numChips, numProbes);
  }
  m_ChipCount = numChips;
  m_ProbeCount = numProbes;
}

/** 
 * @brief Clear out all the data.
 */
void QuantPlierBase::clear() {
  unsigned int probeIx = 0, chipIx = 0;
  for(chipIx = 0; chipIx < m_MaxChips; chipIx++) {
    for(probeIx = 0; probeIx < m_MaxProbes; probeIx++) {
      m_PM[chipIx][probeIx] = 0;
        m_MM[chipIx][probeIx] = 0;
    }
  }
  
  for(chipIx = 0; chipIx < m_MaxChips; chipIx++) 
    m_ChipEffects[chipIx] = 0;
  
  for(probeIx = 0; probeIx < m_MaxProbes; probeIx++)
    m_ProbeEffects[probeIx] = 0;
}

void QuantPlierBase::setParameters(PsBoard &board) {
  string featEff = board.getOptions()->getOpt("use-feat-eff");
  if (!featEff.empty()) {
      Verbose::out(1, "Plier loading feature effects.");
        ChipLayout *layout = NULL;
        board.get("chiplayout", &layout);
      int size = layout->m_PlFactory.getApidMax();
      double *effects = new double[size];
      fill_n(effects, size, 0);
      int probeCount = QuantMethodFactory::openPrecompFeatureEffects(featEff, effects, *layout);
      Verbose::out(1, "Loaded " + ToStr(probeCount) + " feature effects.");
      setFeaturePriorEffects(effects, size);
      FreezArray(effects);
    }
}
