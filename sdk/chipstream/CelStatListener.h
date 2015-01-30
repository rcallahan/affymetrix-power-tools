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

#ifndef _CELSTATLISTENER_H_
#define _CELSTATLISTENER_H_

#include "chipstream/CelListener.h"
#include "chipstream/CumulativeStats.h"
//

/** 
 * Class for objects that need to actually access real cel file data.
 */
class CelStatListener : public CelListener, public ChipSummary {

public:

  /** Constructor. */
  CelStatListener(std::map<std::string,std::vector<bool> > &masks) {

      std::map<std::string, std::vector<bool> >::iterator iter;

      for(iter=masks.begin(); iter != masks.end(); iter++) {
        addProbeMask(iter->first, iter->second);
      }

  }

  /** Virtual class, virtual destructor */
  virtual ~CelStatListener() {}
  
  /** 
   * Cel file of data.
   *
   * @param cel - Handle to an open cel phone.
   */
  virtual void newChip(affymetrix_fusion_io::FusionCELData *cel) {
    std::vector<std::vector<bool> >::iterator maskIx;
    std::vector<CumulativeStats<double> > chipStats(m_ProbeMasks.size());
    calcStats(m_ProbeMasks, cel, chipStats);

    std::vector<ChipSummary::Metric> metrics;
    for(size_t i = 0; i<chipStats.size(); i++) {
        metrics.push_back(ChipSummary::Metric(m_MaskNames[i] + "_mean",chipStats[i].getMean()));
    }
    m_SummaryStats.push_back(metrics);
    setValid(true);
  }

private:

  /** 
   * Add a new probe subset via a mask associated with identifier name.
   *
   * @param name - string identifier associated with this mask (i.e. "all", "pm")
   * @param mask - bit mask with probes to calculate stats for set to true.
   */
  void addProbeMask(const std::string &name, const std::vector<bool> &mask) {
    m_MaskNames.push_back(name);
    m_ProbeMasks.push_back(mask);
    declareMetric(name + "_mean",ChipSummary::Metric::Double);
  }

  /** 
   * Loop though an open cel file and fill in the statistics for each
   * probe subset.
   * 
   * @param masks - Different probe subsets to calculate stats for.
   * @param cel - Cel file to get data from.
   * @param chipStats - Vector to fill in with stats for each mask.
   */
  void calcStats(const std::vector<std::vector<bool> > &masks, 
                 affymetrix_fusion_io::FusionCELData *cel,
                 std::vector<CumulativeStats<double> > &chipStats) {
    assert(chipStats.size() == masks.size());
    int numCells = cel->GetNumCells();
    /* Get each entry in a cel file and then add it to different masks
       as appropriate. */
    for(int probeIx = 0; probeIx < numCells; probeIx++) {
      std::vector<std::vector<bool> >::const_iterator maskIx;
      std::vector<CumulativeStats<double> >::iterator statIx;
      float intensity = cel->GetIntensity(probeIx);
      /* Check all the masks to see if they include this data point. */
      for(maskIx = masks.begin(), statIx = chipStats.begin(); 
          maskIx != masks.end() && statIx != chipStats.end(); 
          maskIx++, statIx++) {
        if(*(maskIx->begin() + probeIx) == true) {
          statIx->addData(intensity);
        }
      }
    }
  }
  
private:
  /// The individual identifiers for particular subsets of probes.
  std::vector<std::string>  m_MaskNames;
  /// Subsets of probes in the form of bit masks.
  std::vector<std::vector<bool> > m_ProbeMasks;
};

#endif /* _CELSTATLISTENER_H_ */
