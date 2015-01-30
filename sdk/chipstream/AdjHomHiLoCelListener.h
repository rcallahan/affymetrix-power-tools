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

#ifndef ADJHOMHILOCELLISTENER_H
#define ADJHOMHILOCELLISTENER_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/CelListener.h"
#include "chipstream/ChipSummary.h"
#include "chipstream/HomHiLoCelListener.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeSet.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @brief Class for computing hilo qc metric from a cel set of cel files.
 */
class AdjHomHiLoCelListener : public ChipSummary, public CelListener {

public:

  /** 
   * Constructor.
   * @param probeSets - Vector of probesets to compute metric on
   */
  AdjHomHiLoCelListener(std::vector<ProbeListPacked> &randProbesets, 
                        std::vector<ProbeListPacked> &set1Probesets,
                        std::vector<ProbeListPacked> &set2Probesets,
                        std::string randLabel,
                        std::string set1Label,
                        std::string set2Label,
                        double k = 2.0f, 
                        double emThresh = 0.05f, 
                        double binSize = 0.02f, 
                        std::string label = "adj-minhilo",
                        double hiLoDiffCutOff = 2.0f,
                        double hiLoFailedValue = 0.0) :
      m_RandProbesets(randProbesets), 
      m_Set1Probesets(set1Probesets), 
      m_Set2Probesets(set2Probesets),
      m_LabelRand(randLabel),
      m_LabelSet1(set1Label),
      m_LabelSet2(set2Label),
      m_K(k), 
      m_EmThresh(emThresh), 
      m_BinSize(binSize),
      m_Label(label),
      m_HiLoDiffCutoff(hiLoDiffCutOff),
      m_HiLoFailedValue(hiLoFailedValue) { 

        m_RandHiLo = new HomHiLoCelListener(m_RandProbesets, randLabel);
        m_Set1HiLo = new HomHiLoCelListener(m_Set1Probesets, set1Label);
        m_Set2HiLo = new HomHiLoCelListener(m_Set2Probesets, set2Label); 

        declareMetric(m_Label,ChipSummary::Metric::Double);
        
        declareMetrics(m_RandHiLo->getMetricDefs());
        declareMetrics(m_Set1HiLo->getMetricDefs());
        declareMetrics(m_Set2HiLo->getMetricDefs());
  }
  
  /**
   * Virtual destructor.
   */
  virtual ~AdjHomHiLoCelListener(){
      delete m_RandHiLo;
      delete m_Set1HiLo;
      delete m_Set2HiLo;
  }

  /** 
   * Process another cel files worth of data.
   * @param cel - Filehandle to open cel file.
   */
  void newChip(affymetrix_fusion_io::FusionCELData *cel);

  /** 
   * Get the names for the cel files that have been seen.
   * @return - vector of cel file names, one for each cel file in order seen.
   */
  std::vector<std::string> getCelNames() {
    return m_CelNames;
  }

protected:

  /// Vector of probesets to compute HiLo on
  /// Pointer memory is owned elsewhere...
  std::vector<ProbeListPacked> m_RandProbesets;
  std::vector<ProbeListPacked> m_Set1Probesets;
  std::vector<ProbeListPacked> m_Set2Probesets;
  /// HomHiLoCelListeners for underlying metrics
  HomHiLoCelListener *m_RandHiLo;
  HomHiLoCelListener *m_Set1HiLo;
  HomHiLoCelListener *m_Set2HiLo;
  /// Labels for the underlying HiLo listeners
  std::string m_LabelRand;
  std::string m_LabelSet1;
  std::string m_LabelSet2;
  /// Name of the cel files that have been called.
  std::vector<std::string> m_CelNames;
  /// Our K parameter for the contrast centers transformation.
  double m_K; 
  /// EM Threshold to use
  double m_EmThresh;
  /// Size of contrast bins
  double m_BinSize;
  /// Label for the chip summary output
  std::string m_Label;
  /// Cut-off for when to use rand vs set to m_HiLoFailedValue
  double m_HiLoDiffCutoff;
  /// Value to report when the cutoff is not meet
  double m_HiLoFailedValue;
};

#endif /* ADJHOMHILOCELLISTENER_H */
