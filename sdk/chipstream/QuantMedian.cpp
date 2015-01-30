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
 * @file   QuantMedian.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 10:35:48 2005
 * 
 * @brief Class for doing probe set quantification using median
 */

//
#include "chipstream/QuantMedian.h"
//
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/ProbeSet.h"
//
#include "stats/stats.h"
#include "util/Verbose.h"
//


using namespace std;

/** 
 * What version of the algorithm are we implementing?
 * @return version string.
 */
std::string QuantMedian::getVersion() {
  return QUANTMEDIAN_VERSION;
}

/**
 * Basic constructor.
 */
QuantMedian::QuantMedian(bool doMean, bool attenuate, float l, float h) { 
  setupSelfDoc(*this);
  m_Type = std::string(QUANTMEDIANSTR);
  m_DoMean = doMean;
  m_Attenuate = attenuate;
  m_L = l;
  m_H = h;
  setOptValue("mean", m_DoMean);
  setOptValue("attenuate", m_Attenuate);
  setOptValue("l", ToStr(m_L));
  setOptValue("h", ToStr(m_H));
  m_ChipCount = 0;
  m_ProbeCount = 0;
}

/** 
 * @brief Does this quantification method estimate feature effects?
 * @return true or false if getFeatureEffects() supplies estimates.
 */
bool QuantMedian::haveFeatureEffects() { return false;}

/** 
 * @brief Does this quantification method supply residuals?
 * @return Return true if getResidual() works, false otherwise.
 */
bool QuantMedian::haveResiduals() { return false;}

/** 
 * @brief Get the feature effect of a particular probe.
 * @param probeIx - Index of probe in probe sets.
 * @return feature (probe) effects.
 */
double QuantMedian::getFeatureEffect(unsigned int probeIx) { 
  Err::errAbort("QuantMedian::getFeatureEffect() - Don't have feature effects.");
  return 0.0;
}
  
/** 
 * @brief Get the estimated intensity for a chip. This is usually the data of
 * primary interest.
 * @param chipIx - Index of chip in experiment.
 * @return estimated intensity.
 */
double QuantMedian::getTargetEffect(unsigned int chipIx) {
  assert(chipIx < m_ChipCount);
  return m_Summaries[chipIx];
}

/** 
 * @brief Get the residual (the intensity not explained by
 * model). Supplied as log_2() transformed.
 * @param probeIx - Index of probe in probe sets.
 * @param chipIx - Chip index in experiment.
 * @return intensity not explained by model.
 */
double QuantMedian::getResidual(unsigned int probeIx, unsigned int chipIx) {
  Err::errAbort("QuantMedian::getResidual() - Don't have residuals.");
  return 0.0;
}

/** 
 * @brief Get the estimated intensity for a chip. This is usually the data of
 * primary interest. Alias for getTargetEffect()
 * @param chipIx - Index of chip in experiment.
 * @return estimated intensity.
 */
double QuantMedian::getSignalEstimate(unsigned int chipIx) {
  return getTargetEffect(chipIx);
}


/** 
 * @brief Set the PM data at probe set probe index and chip index.
 * 
 * @param probeIx - Index of probe in probe set.
 * @param chipIx - Index of chip.
 * @param data - Data to be set.
 */
void QuantMedian::setPMDataAt(unsigned int probeIx, unsigned int chipIx, double data) {
  assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
  m_PM[chipIx][probeIx] = (float)data;
}

/** 
 * @brief Get the PM data at probe set probe index and chip index.
 * 
 * @param probeIx - Index of probe in probe set.
 * @param chipIx - Index of chip.
 */
double QuantMedian::getPMDataAt(unsigned int probeIx, unsigned int chipIx) {
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
void QuantMedian::setMMDataAt(unsigned int probeIx, unsigned int chipIx, double data) {
  assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
  m_MM[chipIx][probeIx] = (float)data;
}

/** 
 * What scale are the target effects on? Are they linear (like
 * plier) or log2 (like median)
 * @return - What scale are the results provided in?
 */
QuantMethod::Scale QuantMedian::getScale() {
  return  Linear;
}

  /** 
   * @brief Return number of features (probes).
   * @return feature count.
   */
unsigned int QuantMedian::getNumFeatures()  { return m_ProbeCount; }

/** 
 * @brief Return number of targets (experiments or chips).
 * @return target count
 */
unsigned int QuantMedian::getNumTargets() { return m_ChipCount; }

/** 
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> QuantMedian::getDefaultDocOptions() { 
  std::vector<SelfDoc::Opt> opts;
  SelfDoc::Opt expon = {"mean", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                        "Use the mean rather than the median."};
  opts.push_back(expon);
  SelfDoc::Opt attenuate = {"attenuate", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                            "Indicate whether or not to attenuate mismatch value when using non-PM-only adjuster."};
  opts.push_back(attenuate);
  SelfDoc::Opt l = {"l", SelfDoc::Opt::Double, "0.005", "0.005", "0", "1",
                    "Tunable parameter for attenuating mismatch value for non-PM-only adjusters."};
  opts.push_back(l);
  SelfDoc::Opt h = {"h", SelfDoc::Opt::Double, "-1", "-1", "-1", "NA",
                    "Used fixed constant to attenuate mismatch."};
  opts.push_back(h);
    
  return opts;
}

/** 
 * Fill in the infomediantion for Self documentation.
 * @param doc - Self documenter to be filled in.
 */
void QuantMedian::setupSelfDoc(SelfDoc &doc) {
  doc.setDocName(QUANTMEDIANSTR);
  doc.setDocDescription("Use the median of probes for a particular chip as the summary.");
  doc.setDocOptions(getDefaultDocOptions());
}

/** 
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc QuantMedian::explainSelf() { 
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
SelfCreate *QuantMedian::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  bool mean = false;
  bool attenuate = true;
  float l = 0.005;
  float h = -1;
  fillInValue(mean, "mean", param, doc);
  fillInValue(attenuate, "attenuate", param, doc);
  fillInValue(l, "l", param, doc);
  fillInValue(h, "h", param, doc);
  QuantMedian *median = new QuantMedian(mean,attenuate,l,h);
  return median;
}

/** 
 * @brief Clear out all the data.
 */
void QuantMedian::clear() {
  unsigned int probeIx = 0, chipIx = 0;
  m_ChipCount = 0;
  m_ProbeCount = 0;
  // Set the world to 0
  for(chipIx = 0; chipIx < m_PM.size(); chipIx++) {
    for(probeIx = 0; probeIx < m_PM[0].size(); probeIx++) {
      m_PM[chipIx][probeIx] = 0;
      m_MM[chipIx][probeIx] = 0;
    }
  }
  for(chipIx = 0; chipIx < m_Summaries.size(); chipIx++) 
    m_Summaries[chipIx] = 0;
}

/** 
 * Set the number of probes and number of chips to be used for estimation.
 * @param numProbes - Probe count.
 * @param numChips - Microarray count.
 */
void QuantMedian::setBounds(unsigned int numProbes, unsigned int numChips) {
  unsigned int i = 0;
  m_ProbeCount = numProbes;
  m_ChipCount = numChips;
  /* Check empty case. */
  if(m_PM.empty()) {
    m_Summaries.resize(numChips);
    while(m_PM.size() < numChips) {
      m_PM.push_back(vector<float>(numProbes));
      m_MM.push_back(vector<float>(numProbes));
    }
  }
  /* Check number of columns (probes). */
  if(numProbes > m_PM[0].size()) {
    for(i = 0; i < m_PM.size(); i++) {
      m_PM[i].resize(numProbes);
      m_MM[i].resize(numProbes);
    }
  }
  /* Check number of rows (chips). */
  if(numChips > m_PM.size()) {
    m_Summaries.resize(numChips);
    while(m_PM.size() < numChips) {
      m_PM.push_back(vector<float>(m_PM[0].size()));
      m_MM.push_back(vector<float>(m_MM[0].size()));
    }
  }
}

/** 
 * @brief Compute median. If PM-MM is zero
 * threshold at a small positive number to avoid taking logs of
 * numbers that aren't positive.
 */
void QuantMedian::computeEstimate() {
  float intensity = 0;
  unsigned int colIx = 0, rowIx = 0;
  
  for(rowIx = 0; rowIx < m_ChipCount; rowIx++) {
    for(colIx = 0; colIx < m_ProbeCount; colIx++) {
      if(m_Attenuate) {
        intensity = attenuateIntensity(m_PM[rowIx][colIx],m_MM[rowIx][colIx],m_L,m_H);
      }
      else {
        intensity = m_PM[rowIx][colIx] - m_MM[rowIx][colIx];
      }
      m_PM[rowIx][colIx] = intensity;
    }
  }
  // calculate the medians
  for(rowIx = 0; rowIx < m_ChipCount; rowIx++) {
    if(m_DoMean) {
      m_Summaries[rowIx] = average(m_PM[rowIx].begin(), m_PM[rowIx].begin() + m_ProbeCount);
    }
    else {
      m_Summaries[rowIx] = median(m_PM[rowIx].begin(), m_PM[rowIx].begin() + m_ProbeCount);
    }
  }
}

/** 
 * @brief Set up the quantification method given all the data about the probe
 * set, chip layout and data.
 * 
 * @param psGroup - Probes to be used for final estimate.
 * @param layout - Chip layout annotation.
 * @param iMart - Raw data from chips.
 * @param iTrans - Transfomediantions to be applied to data before use.
 * @param pmAdjust - How to estimate background, or MM probe.
 * @return bool - True if setup sucessful, false otherwise.
 */
bool QuantMedian::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
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
  setBounds(atomPmCount, chipCount);
  
  atomPmCount = 0;
  m_Probes.clear();
  /* Loop through and fill in the data as coming from the intensity mart. */
  for(psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
    const ProbeSet *ps = psGroup.probeSets[psIx];
    if(ps == NULL)
      continue;
    for(atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom &atom = *(ps->atoms[atomIx]);
      unsigned int channelIx = atom.getChannelCode();
      for(probeIx = 0; probeIx < atom.probes.size(); probeIx++) {
        Probe *p = atom.probes[probeIx];
        if(p->type == Probe::PMST || p->type == Probe::PMAT) {
          unsigned int probeIndex = p->id;
          m_Probes.push_back(p);
          for(chipIx = 0; chipIx < chipCount; chipIx++) {
            /* PM transformation. */
            intensity = transformPrimaryData(probeIndex, chipIx, iMart, iTrans, channelIx);
            float mmIntensity = 0;
            /* MM transformation. */
	    ///@todo use channel to fetch/set correct values
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
