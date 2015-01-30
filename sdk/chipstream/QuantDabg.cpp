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
 * @file   QuantDabg.cpp
 * @author Chuck Sugnet
 * @date   Thu Dec 29 17:22:18 2005
 * 
 * @brief  Class for doing Detected Above BackGround (DABG) analysis. Idea is to
 * use the intensity of probes in a probeset to calculate a p-value for the
 * entire probeset. The individual intensities are compared to a set of
 * background probes.
 * 
 */

//
#include "chipstream/QuantDabg.h"
//
#include "portability/affy-base-types.h"
#include "util/Util.h"
#include "util/Verbose.h"
#include "util/md5sum.h"
//
#include <cfloat>
//

using namespace std;
using namespace affx;

/** 
 * Set the number of probes and number of chips to be used for estimation.
 * @param numProbes - Probe count.
 * @param numChips - Microarray count.
 */
void QuantDabg::setBounds(unsigned int numProbes, unsigned int numChips) {
  unsigned int i = 0;
  m_ProbeCount = numProbes;
  m_ChipCount = numChips;
  /* Check empty case. */
  if (m_PM.empty()) {
    m_SetPvals.resize(numChips);
    m_GcProbeCounts.resize(numProbes);
    while (m_PM.size() < numChips) {
      m_PM.push_back(vector<double>(numProbes));
      m_ProbePvals.push_back(vector<double>(numProbes));
    }
  }
  /* Check number of columns (probes). */
  if (numProbes > m_PM[0].size()) {
    m_GcProbeCounts.resize(numProbes);
    for (i = 0; i < m_PM.size(); i++) {
      m_PM[i].resize(numProbes);
      m_ProbePvals[i].resize(numProbes);
    }
  }
  /* Check number of rows (chips). */
  if (numChips > m_PM.size()) {
    m_SetPvals.resize(numChips);

    while (m_PM.size() < numChips) {
      m_PM.push_back(vector<double>(m_PM[0].size()));
      m_ProbePvals.push_back(vector<double>(m_ProbePvals[0].size()));
    }
  }
}

/** 
 * Setup the dabg objects with intensity data from background probes. This should
 * be called only once.
 * 
 * @param iMart - Contains data for background probes for each chip.
 * @param iTrans - Chipstream transformations for raw data in iMart.
 * @param pmAdjust - Adjuster for perfect match probes.
 */
void QuantDabg::initializeDabg(const IntensityMart &iMart, 
                               std::vector<ChipStream *> &iTrans, 
                               PmAdjuster &pmAdjust) {
  m_Dabgs.resize(m_ChipCount);
  if (m_BgProbes.empty()) {
    Err::errAbort("QuantDabg::initializeDabg() - Must specify at least some probes for background distribution.");
  }

  /* load up the dabg objects for each chip with background data. */
  for (unsigned int chipIx = 0; chipIx < m_ChipCount; chipIx++) {
    Dabg &dabg = m_Dabgs[chipIx];
    dabg.intensity_reserve(m_BgProbes.size());
    for (unsigned int probeIx = 0; probeIx < m_BgProbes.size(); probeIx++) {
      int id = m_BgProbes[probeIx];
      float intensity = transformPrimaryData(id, chipIx, iMart, iTrans);
      assert(id < m_ProbeGcVec.size());
      if (m_ProbeGcVec[id] == (char)NULLPROBEGC) {
        Err::errAbort("QuantDabg - Need probe index: " + ToStr(id) + " GC count for background, but wasn't loaded.");
      }
      dabg.add_intensity(m_ProbeGcVec[id], intensity);
    }
    dabg.prepare();
  }
  m_DabgInit = true;
}

/** 
 * @brief Set up the quantification method given all the data about the probe
 * set, chip and data.
 * 
 * @param psGroup - Probes to be used for final estimate.
 * @param iMart - Raw data from chips.
 * @param iTrans - Transformations to be applied to data before use.
 * @param pmAdjust - How to estimate background, or MM probe.
 * @return bool - True if setup sucessful, false otherwise.
 */
bool QuantDabg::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
                      std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust) {
  unsigned int chipIx = 0, probeIx = 0, psIx, atomIx = 0;
  unsigned int chipCount = 0, atomPmCount = 0;
  double intensity = 0.0;
  assert(psGroup.probeSets.size() > 0);

  chipCount = iMart.getCelFileCount();
  atomPmCount = psGroup.countPmProbes();
  if (atomPmCount == 0)
    return false;
  /* Prepare the quantification method for this many atoms and chips. */
  setBounds(atomPmCount, chipCount);
  if (!m_DabgInit) {
    initializeDabg(iMart, iTrans, pmAdjust);
  }
  atomPmCount = 0;
  m_Probes.clear();
  m_GcProbeCounts.clear();
  /* Loop through and fill in the data as coming from the intensity mart. */
  for (psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
    const ProbeSet *ps = psGroup.probeSets[psIx];
    if (ps == NULL)
      continue;
    for (atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom &atom = *(ps->atoms[atomIx]);
      unsigned int channelIx = atom.getChannelCode();
      for (probeIx = 0; probeIx < atom.probes.size(); probeIx++) {
        Probe *p = atom.probes[probeIx];
        if (p->type == Probe::PMST || p->type == Probe::PMAT) {
          unsigned int probeIndex = p->id;
          int gcCount = m_ProbeGcVec[probeIndex];
          if (gcCount == NULLPROBEGC) {
            Err::errAbort("Unable to figure out GC count for probe with id: " + ToStr(probeIndex + 1));
          }
          //if (m_Dabgs[0].intensity_count(gcCount)==0) {
          //  std::string errorMessage = "No background probes found for GC bin " + ToStr(gcCount) + " for probe with id: " + ToStr(probeIndex + 1);
          //  if (m_BgpFile != "") {
          //    errorMessage.append("\nPerhaps there is no probe with GC count " + ToStr(gcCount) + " in bgp-file: " + m_BgpFile);
          //  }       
          //  Err::errAbort(errorMessage);
          //}
          m_GcProbeCounts.push_back(gcCount);
          m_Probes.push_back(p);
          for (chipIx = 0; chipIx < chipCount; chipIx++) {
            /* PM transformation. */
            intensity = transformPrimaryData(probeIndex, chipIx, iMart, iTrans, channelIx);
            ///@todo use channel to set correct pmdata
            setPMDataAt(atomPmCount, chipIx, intensity);
          }
          atomPmCount++;
        }
      }
    } /* atoms. */
  } /* probe sets. */
  return true;
}


/** 
 * Which probes should be used as the background distribution.
 * @param bgProbes - vector of probes to be used, memory 
 * owned elsewhere.
 */
void QuantDabg::setBgProbes(const std::vector<int> &bgProbes) {
  m_BgProbes = bgProbes;
  vector<int>::iterator probeIx;
  md5sum md5;
  for (probeIx = m_BgProbes.begin(); probeIx != m_BgProbes.end(); ++probeIx) {
    md5.update_nbo((*probeIx));
  }
  md5.final(m_ProbeMd5Sum);
  setOptValue("subsetmd5", m_ProbeMd5Sum);
}

/** 
 * Which probes should be used as the background distribution.
 * @param bgProbes - vector of probes to be used, memory 
 * owned elsewhere.
 */
void QuantDabg::setBgProbes(const std::vector<Probe *> &bgProbes) {
  m_BgProbes.reserve(bgProbes.size());
  int probeIx = 0;
  md5sum md5;
  for (probeIx = 0; probeIx < bgProbes.size(); ++probeIx) {
    md5.update_nbo(bgProbes[probeIx]->id);
    m_BgProbes.push_back(bgProbes[probeIx]->id);
  }
  md5.final(m_ProbeMd5Sum);
  setOptValue("subsetmd5", m_ProbeMd5Sum);
}

/** 
 * Return the -1 * log_10(val). Val is capped at DBL_MIN to avoid taking lots of values <= 0. 
 * @param val - Regular p-value.
 * @return -1 * log_10(val)
 */
double QuantDabg::negLog10Prob(double val) {
  val = Max(val, DBL_MIN);
  return -1 * log(val)/log(10.0);
}

void QuantDabg::computeEstimate() {
  for (unsigned int chipIx = 0; chipIx < m_ChipCount; chipIx++) {
    double setPvalue = 0;
    Dabg &dabg = m_Dabgs[chipIx];
    if (m_DoFisher) {
      dabg.compute_target_fisher(m_ProbeCount, &m_GcProbeCounts[0], 
                                 &m_PM[chipIx][0], m_Adjustment, &setPvalue,
                                 &m_ProbePvals[chipIx][0]);
    }
    else if (m_DoPercentile) {
      dabg.compute_target_percentile(m_ProbeCount, &m_GcProbeCounts[0], 
                                 &m_PM[chipIx][0], m_PercentileEst,
                                 &setPvalue, &m_ProbePvals[chipIx][0]);
    }
    else {
      Err::errAbort("QuantDabg::computeEstimate() - Must set m_DoPercentile or m_DoFisher.");
    }
    if (m_ReportLog) {
      setPvalue = negLog10Prob(setPvalue);
      for (unsigned int probeIx = 0; probeIx < m_ProbeCount; probeIx++) {
        m_ProbePvals[chipIx][probeIx] = negLog10Prob(m_ProbePvals[chipIx][probeIx]);
      }
    }
    m_SetPvals[chipIx] = setPvalue;
  }
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
SelfCreate *QuantDabg::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  bool doFisher = true, doPercentile = false, reportLog = false;
  int adjustment= 0;
  double percentile = .75;
  fillInValue(doFisher, "chisq", param, doc);
  fillInValue(percentile, "percentile", param, doc);
  fillInValue(doPercentile, "usepercentile", param, doc);
  fillInValue(reportLog, "neglog10", param, doc);
  fillInValue(adjustment, "adjust-dof", param, doc);
  if (doPercentile) 
    doFisher = false;
  else if (doFisher) 
    doPercentile = false;
  QuantDabg *dabg = new QuantDabg(doFisher, doPercentile, percentile, reportLog, adjustment);
  return dabg;
}
  
