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

#ifndef VARIABILITYSCORECELLISTENER_H
#define VARIABILITYSCORECELLISTENER_H

#include "chipstream/AnalysisStreamFactory.h"
#include "chipstream/BioTypes.h"
#include "chipstream/CelReader.h"
#include "chipstream/ChipSummary.h"
#include "chipstream/GenoUtility.h"
#include "chipstream/MultiQuantMethodListener.h"
#include "chipstream/PmOnlyAdjust.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QCProbesetOptions.h"
#include "chipstream/QuantLabelZ.h"
#include "chipstream/QuantMethodFactory.h"
#include "chipstream/QuantMethodGTypeReport.h"
#include "chipstream/QuantMethodReportListener.h"
#include "chipstream/SparseMart.h"
//
#include "algorithm/em/PrimeEM.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @brief Class for computing hilo qc metric from a cel set of cel files.
 */
class VariabilityScoreCelListener : public ChipSummary, public CelListener{

public:

  /**
   * Constructor.
   * @param probeSets - Vector of probesets to compute metric on
   */
  VariabilityScoreCelListener(std::vector<ProbeListPacked>& probeSets,
                              QCProbesetOptions& psOpts,
                              ChipLayout* layout,
                              const std::string& label = "D_score");

  /**
       * Virtual destructor.
   */
  virtual ~VariabilityScoreCelListener(){}

  /**
   * Process another cel files worth of data.
   * @param cel - Filehandle to open cel file.
   * @param probeSets - list of probelist
   * @param chipStats - cumulative stats object to compute statistics
   * @param channel - purely for logging purposes report the channel
   */

  void newMultiChannelChip(int celIx); //index of metaCelFileNames

  /**
   * calculate the summary stats for a given probeset
   * @param cel - Filehandle to open cel file
   * @param probeSets - list of probelist
   * @param chipStats - cumulative stats object to compute statistics
   * @param channel - purely for logging purposes report the channel
   */
  void calcStats(vector<ProbeListPacked>& probeSets,
                 affymetrix_fusion_io::FusionCELData* cel,
                 CumulativeStats<double>& chipStats,
                 int channel);

  /**
   * calculate the summary stats for a given probeset
   * @param cel - Filehandle to open cel file
   * @param probeSets - list of probelist
   * @param chipStats1 - cumulative stats object to compute statistics - v_score
   * @param chipStats2 - cumulative stats object to compute statistics2 - cv_score
   * @param channel - purely for logging purposes report the channel
   */
  void calcStats(vector<ProbeListPacked>& probeSets,
                 affymetrix_fusion_io::FusionCELData* cel,
                 CumulativeStats<double>& chipStats1,
                 CumulativeStats<double>& chipStats2,
                 int channel);


  /**
   * Get the names for the pair of cel files that have been seen.
   * @return - vector of cel file names, one for each cel file in order seen.
   */
  std::vector<std::string> getCelNames() {
    return m_CelNames;
  }

  void newChip(affymetrix_fusion_io::FusionCELData *cel);
  /**
   * Set Label
   */
  void setLabel(const std::string &label) {
      m_Label = label;
  }

 protected:
   /** fill in the probe set with the allele values */
  bool fillInAlleleProbeSets(const ProbeSet &gtPs, ProbeSet &aAllele, ProbeSet &bAllele);

  /// Vector of probesets to compute HiLo on
  /// Pointer memory is owned elsewhere...
  std::vector<ProbeListPacked> m_ProbeSets;
  /// Name of the cel files that have been called.
  QCProbesetOptions m_psOpts;
  ChipLayout* m_layout;
  std::vector<std::string> m_CelNames;
  std::string m_Label;
};

#endif /* VARIABILITYSCORECELLISTENER*/
