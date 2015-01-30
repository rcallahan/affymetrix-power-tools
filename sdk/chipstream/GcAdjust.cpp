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
 * @file   GcAdjust.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 09:47:42 2005
 * 
 * @brief Class for doing an adjustment based on the median intensity of probes
 * with similar GC content.
 */

//
#include "chipstream/GcAdjust.h"
//
#include "chipstream/ChipLayout.h"
#include "chipstream/QuantMethod.h"
//
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Verbose.h"
#include "util/md5sum.h"
//

using namespace std;
using namespace affx;

/** 
 * @brief Constructor taking information about both the chip and
 * probes to be used for estimating parameters.
 * 
 * @param layout - Annotation of probes on microarray.
 * @param controlProbes - Probes to be used for estimating parameters.
 */
void GcAdjust::setLayout(ChipLayout &layout, vector<Probe *> &controlProbes) {
  m_Probes.resize(controlProbes.size());
  /* Calculate md5 sum of probe ids. */
  md5sum md5;
  for(int bgIx = 0; bgIx != controlProbes.size(); ++bgIx) {
    //md5.update(&(*bgIx)->id, sizeof((*bgIx)->id));
    md5.update_nbo(controlProbes[bgIx]->id);
    m_Probes[bgIx] = controlProbes[bgIx]->id;
  }
  md5.final(m_ProbeMd5Sum);
  setOptValue("subsetmd5", m_ProbeMd5Sum);

  /* Remember that haven't filled the bins yet. */
  m_BinsFilled = false;
  /* Get map of probe to probe id to gc content. */
  m_ProbeGcVec = layout.getGcProbeVec();
}

/** 
 * @brief Constructor taking information about both the chip and
 * probes to be used for estimating parameters.
 * 
 * @param board - Blackboard containing lots of state.
 */
void GcAdjust::setParameters(PsBoard &board) {
  DataStore *info = board.getProbeInfo();
  info->getGcControlProbes(m_Probes);
  info->getProbeGc(m_ProbeGcVec);

  /* Calculate md5 sum of probe ids. */
  md5sum md5;
  for(int bgIx = 0; bgIx != m_Probes.size(); ++bgIx) {
    md5.update_nbo(m_Probes[bgIx]);
  }
  md5.final(m_ProbeMd5Sum);
  setOptValue("subsetmd5", m_ProbeMd5Sum);

  /* Remember that haven't filled the bins yet. */
  m_BinsFilled = false;
}

/** 
 * @brief Set the probes necessary for estimating GC background to
 * true so they will be loaded.
 * 
 * @param probes - bitmask for probe to be loaded.
 */
void GcAdjust::setProbes(vector<bool> &probes) {
  md5sum md5;
  unsigned int i = 0;
  for(i = 0; i < m_Probes.size(); i++) {
    md5.update_nbo(m_Probes[i]);
    probes[m_Probes[i]] = true;
  }
  md5.final(m_ProbeMd5Sum);
}

/** 
 * @brief Fill in the GC bins using the raw data in intensity mart
 * and processed by iTrans ChipStream.
 * 
 * @param chipIx - Index of chip to be filled in.
 * @param iMart - Raw data to be used.
 * @param iTrans - Chipstreams that will be used to transform the raw data.
   */
void GcAdjust::fillChipBins(int chipIx, const IntensityMart &iMart, std::vector<ChipStream *> &iTrans) {
  assert(chipIx >= 0 && chipIx < iMart.getCelDataSetCount());
  unsigned int probeIx = 0, binIx = 0;
  int pIx = -1;
  unsigned int binCount = m_Bins[chipIx].size();
  vector<vector<float> > gcBins;

  for(binIx = 0; binIx < binCount; binIx++) {
    gcBins.push_back(vector<float>());
  }

  // Fill in the bins with raw data processed by chip streams.
  for(probeIx = 0; probeIx < m_Probes.size(); probeIx++) {
    float intensity = 0;
    unsigned int bin = 0;
    pIx = m_Probes[probeIx];
    if(m_ProbeGcVec[pIx] == (char)NULLPROBEGC)
      Err::errAbort("GcAdjust - Need probe index: " + ToStr(pIx) + " GC count for background, but wasn't loaded.");
    bin = m_ProbeGcVec[pIx];
    ///@todo iterate over each channel
    intensity = QuantMethod::transformPrimaryData(pIx, chipIx, iMart, iTrans);
    gcBins[bin].push_back(intensity);
  }
  
  // Estimate the medians for each bin.
  for(binIx = 0; binIx < gcBins.size(); binIx++) {
    float med =  -1.0;
    if(gcBins[binIx].empty()) {
      if(chipIx == 0)
        Verbose::out(2, "Warning: GC Bin " + ToStr(binIx) + " has no data.");
    }
    else
      med = median_in_place(gcBins[binIx].begin(), gcBins[binIx].end());
    m_Bins[chipIx][binIx] = med;
  }
}

/** 
 * @brief Calculate the parameters for each GC count.
 * 
 * @param iMart - Raw data from chips.
 * @param iTrans - Transformations to apply to data.
 */
void GcAdjust::calcParams(const IntensityMart &iMart, std::vector<ChipStream *> &iTrans) {
  int chipIx = 0;
  int binCount = m_MaxGc;
  int chipCount = iMart.getCelDataSetCount();
  assert(chipCount > 0);
  assert(binCount > 0);
  if(m_Probes.empty()) {
    Err::errAbort("GcAdjust::calcParams() - Doesn't look like setLayout() has been called, no control probes.");
  }

  /* Make the correct size of bin matrix. */
  m_Bins.resize(chipCount);
  for(chipIx = 0; chipIx < chipCount; chipIx++) 
    m_Bins[chipIx].resize(binCount);
  
  for(chipIx = 0; chipIx < chipCount; chipIx++) 
    fillChipBins(chipIx, iMart, iTrans);
}


/** 
 * @brief Function to determine how much to adjust a perfect match intensity.
 * @param probeIx - Index of probe in cel file data.
 * @param chipIx - Microarray or chip index.
 * @param iMart - IntensityMart which holds raw data.
 * @param iTrans - Vector of transformations that should be performed on raw data.
 * @param pmIintensity - Intensity of perfect match probe to be adjusted, may
 * be modified from original value depending on adjuster
 * @param bgrdAdjust - Background adjustment, if any, recommended (i.e. MM intensity)
 */
void GcAdjust::pmAdjustment(int probeIx, 
                            int chipIx, 
                            const IntensityMart &iMart, 
                            std::vector<ChipStream *> &iTrans, 
                            float &pmIntensity, 
                            float &bgrdAdjust)
{
  int gcBin = 0;

  // If we haven't calculated parameters yet, do it now.
  if(!m_BinsFilled) {
    calcParams(iMart, iTrans);
    m_BinsFilled = true;
  }

  // Look up probe gc count.
  if(m_ProbeGcVec[probeIx] == (char)NULLPROBEGC) 
    Err::errAbort("Unable to figure out GC count for probe with id: " + ToStr(probeIx + 1));
  gcBin = m_ProbeGcVec[probeIx];
  if(gcBin < 0 || gcBin >= m_MaxGc)
    Err::errAbort("GC count out of accepted range for probe with id: " + ToStr(probeIx + 1));
  if(m_Bins[chipIx][gcBin] < 0.0) 
    //Err::errAbort("No background probe correction available for GC count " + ToStr(gcBin) + ". Failed on GC correction of probe with id: " + ToStr(probeIx + 1));
	bgrdAdjust = 0.0f;

  else // Return the pre calculated median for GC count.
    bgrdAdjust = m_Bins[chipIx][gcBin];
}
