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

#ifndef HOMHILOCELLISTENER_H
#define HOMHILOCELLISTENER_H
//
#include "chipstream/BioTypes.h"
#include "chipstream/CelListener.h"
#include "chipstream/ChipSummary.h"
#include "chipstream/GenoUtility.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeSet.h"
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
class HomHiLoCelListener : public ChipSummary, public CelListener {

public:

  /** 
   * Constructor.
   * @param probeSets - Vector of probesets to compute metric on
   */
  HomHiLoCelListener(std::vector<ProbeListPacked> &probeSets, 
                     std::string label = "minhilo",
                     double k = 2.0f, 
                     double emThresh = 0.05f,
                     double binSize = 0.02f) : 
      m_ProbeSets(probeSets), 
      m_K(k), 
      m_EmThresh(emThresh), 
      m_BinSize(binSize),
      m_Label(label) 
  {
    declareMetric(m_Label,ChipSummary::Metric::Double);
  }

  
  /**
   * Virtual destructor.
   */
  virtual ~HomHiLoCelListener(){}
  
  /** 
   * Process another std::vector<double> worth of data.
   * @param vData - The std::vector<double> of data to process.
   * @param dEmThreshold - The EM Threshold.
   * @param dBinSize - The Bin Size.
   * @return - The statistic.
   */
  static double computeStatistic(std::vector<double>& vData, double dEmThreshold = 0.05f, double dBinSize = 0.04f);

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

  /** 
   * Loop through the probesets provided and calculate a contrast value
   * for each one using the median of PM probes for A allele and B
   * allele. (CES constrast space)
   * @param cel - Cel File to get data from.
   * @param chrXProbeSets - Probesets to process (should be chrX non-pseudo autosomal).
   * @param contrastValues - Contrast values vector to be filled in.
   */
  static void fillInContrastValues(affymetrix_fusion_io::FusionCELData *cel, 
                                   std::vector<ProbeListPacked> &chrXProbeSets, 
                                   std::vector<double> &contrastValues, 
                                   double k = 2);
  /**
   * Calculate contrast values
   */
  static double CalculateEMContrast(affymetrix_fusion_io::FusionCELData* cel, 
                                    const ProbeSet* ps, 
                                    const double k = 2);

  /**
   * Compute the cluster peaks using EM
   */
  static void computeClusterPeaks(vector<double> &contrastValues, 
                                    double &peak1, 
                                    double &peak2, 
                                    double &peak3, 
                                    double emThresh = 0.05,
									double dMinMu = 0.25);

  /**
   * Compute the Valley Between Two Points
   */ 
  static double computeClusterValleyFixedEndPoint(const std::vector<double>& contrast, 
                                     double x1, 
                                     double x2, 
                                     double bin = 0.02);

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
  /**
   * Fit the model N(-delta,sigma), N(0,K*sigma), N(delta,sigma)
   * using EM to calculate parameters sigma, delta, P, K
   */
  static void computeMixtureModel(const std::vector<double>& xdata,
                                  double& sigma,
                                  double& delta,
                                  double& P,
                                  double& K,
                                  bool symmetrize   = true,
                                  double sigmaInit  = 0.2,
                                  double deltaInit  = 1.0,
                                  double PInit      = 0.6,
                                  double KInit      = 0.2);

   /**
    * Calculate snpqc metric
    */
   static float computeSNPQC(const vector<double>& xdata, bool symmetrize = true);

private:

  /**
   * Helper functions for computeMixtureModel()
   */
  static double logLikelihood(const std::vector<double>& xdata,
                              double sigma,
                              double delta,
                              double P,
                              double K);

  static double mixSum(double x, double sigma, double delta, double P, double K);
  static double dnorm(double x, double sigma, double delta);

protected:

  /// Vector of probesets to compute HiLo on
  /// Pointer memory is owned elsewhere...
  std::vector<ProbeListPacked> m_ProbeSets;
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
};

#endif /* HOMHILOCELLISTENER_H */
