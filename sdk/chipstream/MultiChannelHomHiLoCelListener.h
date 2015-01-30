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

#ifndef MULTICHANNELHOMHILOCELLISTENER_H
#define MULTICHANNELHOMHILOCELLISTENER_H

//
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
class MultiChannelHomHiLoCelListener : public ChipSummary, public CelListener {

public:

  // shouldnt have a basic constructor.
  // MultiChannelHomHiLoCelListener();

  /**
   * Constructor.
   * @param probeSets - Vector of probesets to compute metric on
   */
  MultiChannelHomHiLoCelListener(std::vector<ProbeListPacked>& probeSets,
                                 QCProbesetOptions& psOpts,
                                 ChipLayout* layout,
                                 double k = 2.0f,
                                 double emThresh = 0.05f,
                                 double binSize = 0.02f,
                                 std::string label = "minhilo");

  /**
   * Virtual destructor.
   */
  virtual ~MultiChannelHomHiLoCelListener();

  void declareMetrics();

  /**
   * Get the names for the pair of cel files that have been seen.
   * @return - vector of cel file names, one for each cel file in order seen.
   */
  std::vector<std::string> getCelNames() {
    return m_CelNames;
  };

  void newChip(affymetrix_fusion_io::FusionCELData* cel);
  /**
   * Set Label
   */
  void setLabel(const std::string &label) {
    m_Label = label;
  };

  /**
   * Loop through the probesets provided and calculate a contrast value
   * for each one using the median of PM probes for A allele and B
   * allele. (CES constrast space)
   * @param cel - Cel File to get data from.
   * @param chrXProbeSets - Probesets to process (should be chrX non-pseudo autosomal).
   * @param contrastValues - Contrast values vector to be filled in.
   */
  static void fillInContrastValues(int celIx,
                                   std::vector<ProbeList> &probeSets,
                                   std::vector<double> &contrastValues,
                                   double k = 2);
  /**
   * Calculate contrast values
   */
//  static double CalculateEMContrast(affymetrix_fusion_io::FusionCELData* celA,
//                                    affymetrix_fusion_io::FusionCELData* celB,
//                                    const ProbeSet* ps,
//                                    const double k = 2);

  /**
   * Compute the cluster peaks using EM
   */
  static void computeClusterPeaks(const vector<double>& contrastValues,
                                  double &peak1,
                                  double &peak2,
                                  double &peak3,
                                  double emThresh = 0.05);

  /**
   * Compute the Valley Between Two Points
   */
  static double computeClusterValley(const std::vector<double>& contrast,
                                     double x1,
                                     double x2,
                                     double bin = 0.02);

  /**
   * Compute the density of values for a given bin
   */
  static double getContrastDensity(const std::vector<double>& contrast,
                                   double x,
                                   double bin = 0.02);

  /* stole from QuantGTypeMethod.cpp, mainly to extract out summarized allele signals */
  bool summarizeAllele(ProbeSet* pSet,
                       vector<double>& summaryValues,
                       ChipLayout& layout,
                       SparseMart& iMart,
                       std::vector<ChipStream*>* iTrans,
                       PmAdjuster& pmAdjust,
                       QuantMethod* quantMethod,
                       bool lowPrecision);

  void GimmeContrast (std::vector<ProbeListPacked>& probeSets,
                      std::vector<double>& contrastValues,
                      double k);
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
  /// Our K parameter for the contrast centers transformation.
  double m_K;
  /// EM Threshold to use
  double m_EmThresh;
  /// Size of contrast bins
  double m_BinSize;
  /// Label for the chip summary output
  std::string m_Label;
  ///
  affymetrix_fusion_io::FusionCELData* m_cel;

  /** clear a probe set*/
  void clearProbeSet(ProbeSet &ps);
};

#endif /* HOMHILOCELLISTENER_H */
