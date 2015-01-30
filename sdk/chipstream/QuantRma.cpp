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
 * @file   QuantRma.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 10:35:48 2005
 * 
 * @brief Class for doing probe set quantification using median polish
 * as does RMA.
 */

//
#include "chipstream/QuantRma.h"
//
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantMethodFactory.h"
//
#include "rma/RMA.h"
#include "util/Verbose.h"
//

using namespace std;
/** 
 * @brief Clear out all the data.
 */
void QuantRma::clear() {
  unsigned int probeIx = 0, chipIx = 0;
 
  // Set the world to 0
  for(chipIx = 0; chipIx < m_PM.size(); chipIx++) {
    for(probeIx = 0; probeIx < m_PM[0].size(); probeIx++) {
      m_PM[chipIx][probeIx] = 0;
      m_MM[chipIx][probeIx] = 0;
    }
  }
  for(chipIx = 0; chipIx < m_ChipEffects.size(); chipIx++) 
    m_ChipEffects[chipIx] = 0;
  
  for(probeIx = 0; probeIx < m_ProbeEffects.size(); probeIx++)
    m_ProbeEffects[probeIx] = 0;
}

/** 
 * Set the number of probes and number of chips to be used for estimation.
 * @param numProbes - Probe count.
 * @param numChips - Microarray count.
 */
void QuantRma::setBounds(unsigned int numProbes, unsigned int numChips) {
  unsigned int i = 0;
  m_ProbeCount = numProbes;
  m_ChipCount = numChips;
  /* Check empty case. */
  if(m_PM.empty()) {
    m_ChipEffects.resize(numChips);
    m_ProbeEffects.resize(numProbes);
    while(m_PM.size() < numChips) {
      m_PM.push_back(vector<float>(numProbes));
      m_MM.push_back(vector<float>(numProbes));
    }
  }
  /* Check number of columns (probes). */
  if(numProbes > m_PM[0].size()) {
    m_ProbeEffects.resize(numProbes);
    for(i = 0; i < m_PM.size(); i++) {
      m_PM[i].resize(numProbes);
      m_MM[i].resize(numProbes);
    }
  }
  /* Check number of rows (chips). */
  if(numChips > m_PM.size()) {
    m_ChipEffects.resize(numChips);
    while(m_PM.size() < numChips) {
      m_PM.push_back(vector<float>(m_PM[0].size()));
      m_MM.push_back(vector<float>(m_MM[0].size()));
    }
  }
}

/** 
 * @brief Do the heavy lifting of median polish. If PM-MM is zero
 * threshold at a small positive number to avoid taking logs of
 * numbers that aren't positive.
 */
void QuantRma::computeEstimate() {
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
      /* Don't want negative values, threshold at RMA_MIN_VALUE */
      m_PM[rowIx][colIx] = max(intensity, (float)RMA_MIN_VALUE);
    }
  }

  if(m_FitFeatureResponse){
    if(m_UseInputModel){
	  if (m_FixFeatureEffect) {
	    std::vector<std::vector<float> > PM = m_PM;
        RMA::medianPolishWithPrecomputedEffectsUsedAsSeedValues(m_PM, m_ChipCount, m_ProbeCount, m_ProbeEffects, m_ChipEffects);
        RMA::medianPolishWithPrecomputedEffectsOnePass(PM, m_ChipCount, m_ProbeCount, m_ProbeEffects, m_ChipEffects);
	  } else {
        RMA::medianPolishWithPrecomputedEffectsUsedAsSeedValues(m_PM, m_ChipCount, m_ProbeCount, m_ProbeEffects, m_ChipEffects);
	  }
    }else{
	  if (m_FixFeatureEffect) {
	    std::vector<std::vector<float> > PM = m_PM;
        RMA::medianPolishPsetFromMatrix(m_PM, m_ChipCount, m_ProbeCount, m_ProbeEffects, m_ChipEffects);
        RMA::medianPolishWithPrecomputedEffectsOnePass(PM, m_ChipCount, m_ProbeCount, m_ProbeEffects, m_ChipEffects);
	  } else {
        RMA::medianPolishPsetFromMatrix(m_PM, m_ChipCount, m_ProbeCount, m_ProbeEffects, m_ChipEffects);
	  }
    }
  } else{
    if(m_UseInputModel){
      RMA::medianPolishWithPrecomputedEffectsOnePass(m_PM, m_ChipCount, m_ProbeCount, m_ProbeEffects, m_ChipEffects);
    } else{
      Err::errAbort("Invalid combination of Fit Feature Response and Use Input Model, both cannot be false.");
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
 * @param iTrans - Transformations to be applied to data before use.
 * @param pmAdjust - How to estimate background, or MM probe.
 * @return bool - True if setup sucessful, false otherwise.
 */
bool QuantRma::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
                     std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust) {
  unsigned int chipIx = 0, probeIx = 0, psIx, atomIx = 0;
  unsigned int chipCount = 0, atomPmCount = 0;
  float intensity = 0.0;
  assert(psGroup.probeSets.size() > 0);
  chipCount = iMart.getCelFileCount();
  int numberOfProbes = psGroup.countPmProbes();
  if(numberOfProbes == 0)
    return false;
  /* Prepare the quantification method for this many atoms and chips. */
  setBounds(numberOfProbes, chipCount);
  
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
          /* If we are using precomputed probe effects fill in each probe here. */
          if(m_UsePrecompEffects) {    
	    m_ProbeEffects.resize(numberOfProbes);
	    if (m_FeaturePriorEffect == NULL) {
			Verbose::out(1, "trying to use feature effects but null.");
	      Err::errAbort("No feature effects available.");
	    }
	    else {
	      /// Expected feature effect is 1.0 +/- some minor adjustment.
	      try {
		m_ProbeEffects[atomPmCount] = m_FeaturePriorEffect[p->getApid()];
	      } 
	      catch(...) {
			  Verbose::out(1, "Feature effect for probe "+ ToStr(p->id + 1) +" is out of bounds.");
		Err::errAbort("Feature effect for probe "+ ToStr(p->id + 1) +" is out of bounds.");
	      }

	      if (m_ProbeEffects[atomPmCount] <= 0) {
		Err::errAbort("The feature effect values must be positive.  Feature effect for probe "+ ToStr(p->id + 1) +" = "+ ToStr(m_ProbeEffects[atomPmCount]));
	      } // Safety against log zero.
	    }
          }

          m_Probes.push_back(p);
          for(chipIx = 0; chipIx < chipCount; chipIx++) {
            /* PM transformation. */
            intensity = transformPrimaryData(probeIndex, chipIx, iMart, iTrans, channelIx);
            /* MM transformation. */
            float mmIntensity = 0;

            pmAdjust.pmAdjustment(p->id, chipIx, iMart, iTrans, intensity, mmIntensity);
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

/**
 * Custom configuration for this QuantMethod
 */
void QuantRma::setParameters(PsBoard &board) {
  string featEff = board.getOptions()->getOpt("use-feat-eff");
    if (!featEff.empty()) {
      ChipLayout *layout = NULL;
      board.get("chiplayout", &layout);
      int size = layout->m_PlFactory.getApidMax();
      double *effects = new double[size];
      fill(effects, effects+size, 0);
      int probeCount = QuantMethodFactory::openPrecompFeatureEffects(featEff, effects, *layout);
      Verbose::out(1, "Loaded " + ToStr(probeCount) + " feature effects.");
      setFeaturePriorEffects(effects, size);
      FreezArray(effects);
    }
  }
