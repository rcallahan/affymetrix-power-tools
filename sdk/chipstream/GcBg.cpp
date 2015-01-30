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
 * @file   GcBg.cpp
 * @author Chuck Sugnet
 * @date   Mon Sep 25 16:27:31 PDT 2006
 * 
 * @brief Class for doing an background adjustment on all probes based
 * on the median intensity of probes with similar GC content.
 */


//
#include "chipstream/GcBg.h" 
//
#include "chipstream/DataStore.h"
#include "chipstream/QuantMethod.h"
//
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Verbose.h"
#include "util/md5sum.h"
//
#include <vector>
//


using namespace std;
using namespace affx;

void GcBg::learnParameters(const std::vector<float> &data, int binCount, std::vector<vector<float> > &chipBins,
                           std::vector<int> &probes, std::vector<char> &probeGcVec, bool warnings) {
  unsigned int probeIx = 0, binIx = 0;
  vector<vector<float> > gcBins(binCount);
  if(probeGcVec.empty()) {
    Err::errAbort("GcBg::learnParameters() - Must have some gc control probes set to calculate background.");
  }
  // Fill in the bins with raw data processed by chip streams.
  for(probeIx = 0; probeIx < probes.size(); probeIx++) {
    float intensity = 0;
    unsigned int bin = 0;
    int probeIndex = probes[probeIx];
    if(probeGcVec[probeIndex] == (char)NULLPROBEGC)
      Err::errAbort("GcBg - Need probe index: " + ToStr(probeIndex) + " GC count for background, but wasn't loaded.");
    bin = probeGcVec[probeIndex];
    if(bin >= binCount || bin < 0) {
      Err::errAbort("Probe: " + ToStr(probeIndex) + " has gc content of: " + ToStr(probeGcVec[probeIndex]) + 
                    " expecting between 0 and " + ToStr(binCount));
    }
    intensity = data[probeIndex];
    gcBins[bin].push_back(intensity);
  }
  vector<float> binMeds(binCount, 0);
  // Estimate the medians for each bin.
  for(binIx = 0; binIx < gcBins.size(); binIx++) {
    // Unlike GcAdjust, we will default to 0 rather than error. Otherwise
    // one would have to have the GC count for all probes on the array
    float med =  0;
    if(gcBins[binIx].empty()) {
      if(warnings)
        Verbose::out(2, "Warning: GC Bin " + ToStr(binIx) + " has no data.");
    }
    else
      med = median_in_place(gcBins[binIx].begin(), gcBins[binIx].end());
    binMeds[binIx] = med;
  }
  chipBins.push_back(binMeds);
}

float GcBg::transform(int probeIx, int chipIx, float intensity, int gcBin) {
  float result = 0;
  if(gcBin == (char)NULLPROBEGC) {
    Verbose::out(4,"Unable to figure out GC count for probe with id: " + ToStr(probeIx + 1) + ". Using 0 for GC correction.");
  } 
  else {
    if(gcBin < 0 || gcBin >= m_MaxGc) {
      Err::errAbort("GC count out of accepted range for probe with id: " + ToStr(probeIx + 1) +
                    " Cel file: " + m_TransformedIMart->getCelFileNames()[chipIx]);
    }
    // Return the pre calculated median for GC count.
    result = m_Bins[chipIx][gcBin];
  }
  if(m_Attenuate) {
    result = QuantMethod::attenuateIntensity(intensity, result, m_L, m_H);
  }
  else { 
    result = intensity - result;
  }
  return result;
}

float GcBg::transform(int probeIx, int chipIx, float intensity) {
  // Look up probe gc count.
  char gc = m_BgProbeGc[probeIx];
  float result = transform(probeIx, chipIx, intensity, gc);
  return result;
}

void GcBg::transform(int chipIx, std::vector<float>& intensity) {
  for (uint32_t probeIx = 0; probeIx < intensity.size(); probeIx++) {
    intensity[probeIx] = transform(probeIx, chipIx, intensity[probeIx]);
  }
}

void GcBg::newChip(std::vector<float> &data) {
   learnParameters(data, m_MaxGc, m_Bins, m_BgProbes, m_BgProbeGc, m_ChipCount == 0);
}

void GcBg::newDataSet(IntensityMart* iMart) {
  int dataSetCount = iMart->getCelDataSetCount();
  m_TransformedIMart = iMart;
  // if there aren't any chipstream nodes after this, then don't store
  // all of the intensities
  if (m_Streams.empty()) {
    m_TransformedIMart->setStoreAllCelIntensities(false);
  }
  std::vector<float> data;
  for (int d = 0; d < dataSetCount; d++) {
    data = iMart->getCelData(d);
    newChip(data);
    m_ChipCount++;
    transform(d, data);
    m_TransformedIMart->setProbeIntensity(d, data);
  }

  chipStreamPassNewChip(m_TransformedIMart);
}

void GcBg::setControlProbes(std::vector<int> &vec) {
  md5sum md5;
  for(int i = 0; i < vec.size(); i++) {
    md5.update_nbo(vec[i]);
  }
  md5.final(m_ProbeMd5Sum);
  setOptValue("subsetmd5", m_ProbeMd5Sum);
  m_BgProbes = vec;
}

void GcBg::setControlProbes(std::vector<Probe *> &vec) {
  m_BgProbes.resize(vec.size());

  for(int i = 0; i < vec.size(); i++) {
    m_BgProbes[i] = vec[i]->id;
  }

  setControlProbes(m_BgProbes);
}

/** 
 * Setting the GC counts for all the probes on the array.
 * @param vec - GC count for every probe on array indexed by id.
 */
void GcBg::setProbeGcVec(const std::vector<char> &vec) {
  if(m_BgProbes.empty()) {
    APT_ERR_ABORT("GcBG::setProbeGcVec() - Must specify background probe ids first.");
  }
  m_BgProbeGc.resize(vec.size());
  for(int i = 0; i < m_BgProbeGc.size(); i++) {
    m_BgProbeGc[i] = (char) vec[i];
  }
}

/** 
 * Setting the GC counts for all the probes on the array.
 * @param vec - GC count for every probe on array indexed by id.
 */
void GcBg::setProbeGcVec(const std::vector<unsigned char> &vec) {
  if(m_BgProbes.empty()) {
    APT_ERR_ABORT("GcBG::setProbeGcVec() - Must specify background probe ids first.");
  }
  m_BgProbeGc.resize(vec.size());
  for(int i = 0; i < m_BgProbeGc.size(); i++) {
    m_BgProbeGc[i] = (char) vec[i];
  }
}

/**
 * Setup the background probes and background gc indexes by reading
 * them from the blackboard
 */
void GcBg::setParameters(PsBoard &board) {
  std::vector<int> vec;
  board.getProbeInfo()->getGcControlProbes(vec);
  setControlProbes(vec);
  std::vector<char> gcVec;
  board.getProbeInfo()->getProbeGc(gcVec);
  setProbeGcVec(gcVec);
}
