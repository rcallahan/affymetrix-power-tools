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
 * @file   RmaBgTran.cpp
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:34:04 2005
 * 
 * @brief  Class for doing RMA style background subtraction.
 */

//
#include "chipstream/RmaBgTran.h"
//
#include "chipstream/DataStore.h"
#include "chipstream/DiskIntensityMart.h"
//
#include "rma/RMA.h"
#include "chipstream/PsBoard.h"
#include "util/Err.h"
#include "util/Verbose.h"
//

using namespace std;

/** 
 * @brief Constructor.
 * @param numBins - Number of bins to use in density estimation 16384 is magic
 *                  number from biocoductor, must be a power of two.
 */
RmaBgTran::RmaBgTran(int numBins) {
  int pow2, root2;
  setupSelfDoc(*this);
  m_Type = RMABGSTR;
  /* Check to make sure that we have a power of two number of bins. */
  root2 = Util::round(log((double)numBins) / log(2.0));
  pow2 = Util::round(pow(2.0, root2));
  if(pow2 == numBins)
    m_NumBins = numBins;
  else 
    Err::errAbort("RmaBgTran() - Number of bins must be a power of two.");
  m_PmProbeCount = 0;
}

/** 
 * @brief Set which probes to use for the background estimation.
 * 
 * @param pmProbes Bitmask where PM probes are set to true.
 */
void RmaBgTran::setPmProbes(const vector<bool> &pmProbes) {
  unsigned int i = 0;
  m_PmProbes = pmProbes;
  m_PmProbeCount = 0;
  for(i = 0; i < m_PmProbes.size(); i++) {
    if(m_PmProbes[i] == true) 
      m_PmProbeCount++;
  }
  Verbose::out(2, string("Using: ") + ToStr(m_PmProbeCount) + string(" PM probes for analysis."));
}

/**
 * Setup the background probes and background gc indexes by reading
 * them from the blackboard
 */
void RmaBgTran::setParameters(PsBoard &board) {
  std::vector<bool> pm;
  board.getProbeInfo()->getProbePm(pm);
  setPmProbes(pm);
}

/** 
 * @brief transform the intensity point supplied coming from a particular
 * probe in a particular microarray.
 * 
 * @param probeIx - Probe index from the cel file.
 * @param chipIx - Set of chips from same sample.
 * @param intensity - Original intensities.
 * @param return - transformed intensities.
 */
float RmaBgTran::transform(int probeIx, int chipIx, float intensity) {
  bool error = false;
  intensity = RMA::bgSubData(intensity, 
                             m_Params[chipIx].mu, 
                             m_Params[chipIx].sigma,
                             m_Params[chipIx].alpha, 
                             error);
  if (error) {
    Err::errAbort("RmaBgTran::transform() - Error in background subtraction. "
                  "This is normally a result of an intensity distribution which is outside what RMA can process. "
                  "Remove this CEL file and retry. "
                  "ChipIdx: " + ToStr(chipIx) + "(file: " + m_TransformedIMart->getCelFileNames()[chipIx]+") "
                  "Probe: " + ToStr(probeIx));
  }
  return intensity;
}


void RmaBgTran::transform(int chipIx, std::vector<float>& intensity) {
  for (int probeIx = 0; probeIx < intensity.size(); probeIx++) {
    intensity[probeIx] = transform(probeIx, chipIx, intensity[probeIx]);
  }
}


void RmaBgTran::newChip(std::vector<float> &data) {
  initializeData(data);
}

/** 
 * @brief Method for being passed a new cel file worth of data.
 * @param data - cel file vector of data.
 */
void RmaBgTran::newDataSet(IntensityMart* iMart) {
  int dataset_count = iMart->getCelDataSetCount();
  m_TransformedIMart = iMart;
  // if there aren't any chipstream nodes after this, then don't store
  // all of the intensities
  if (m_Streams.empty()) {
    m_TransformedIMart->setStoreAllCelIntensities(false);
  }

  std::vector<float> data;
  for (int d = 0; d < dataset_count; d++) {
    data = iMart->getCelData(d);
    assert(data.size() > 0);
    newChip(data);
    transform(d, data);
    m_TransformedIMart->setProbeIntensity(d, data);
  }

  chipStreamPassNewChip(m_TransformedIMart);
}


/** 
 * @brief Method for adding cel file data to chipstream normalization method.
 * @param data - cel file vector of data.
 */
void RmaBgTran::initializeData(const std::vector<float>& data) {
  struct ChipParam param;
  vector<float> copy(m_PmProbeCount);
  unsigned int currentPm = 0;
  unsigned int i = 0;
  
  /* Sanity checks. */  
  if(m_PmProbes.empty()) 
    Err::errAbort("RmaBgTran::initializeData() - Doesn't appear that PM probes hve been set.");

  if(data.size() != m_PmProbes.size()) 
    Err::errAbort("RmaBgTran::initializeData() - Chip Data size (" + ToStr(data.size()) + ") different than Pm Probe Vector (" + ToStr(m_PmProbes.size()) + ")");
  
  /* First filter out PM probes. */
  for(i = 0; i < m_PmProbes.size(); i++) {
    if(m_PmProbes[i] == true) {      
      copy[currentPm++] = data[i];
    }
  }

  /* Estimate parameters for current chip from pm only probes. */
  RMA::estimateBgParam(copy, &param.mu, &param.sigma, &param.alpha, m_NumBins);
  Verbose::out(2, "Chip " + ToStr(m_Params.size()) + " using RMA params:"+
               " alpha: " + ToStr(param.alpha) +
               " mu: " + ToStr(param.mu) + 
               " sigma: " + ToStr(param.sigma));
  // point out a possble source of error with RMA...
  // this might help the end user remove many bad files in one go.
  if (10.0<param.alpha) {
    Verbose::out(2,"The alpha is rather high.  If the run crashes, this could be the problem.");
  }

  m_Params.push_back(param);  
}

