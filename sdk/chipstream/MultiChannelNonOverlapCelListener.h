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

#ifndef MULTICHANNELNONOVERLAPCELLISTENER_H
#define MULTICHANNELNONOVERLAPCELLISTENER_H

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
class MultiChannelNonOverlapCelListener : public ChipSummary, public CelListener {

public:

  // this doenst exist!
  //MultiChannelNonOverlapCelListener();

  // only this one.
  MultiChannelNonOverlapCelListener(std::vector<ProbeListPacked>& probeSets, 
                                    QCProbesetOptions& psOpts,
                                    ChipLayout* layout, 
                                    std::string label);
  //
  void declareMetrics();
  
  virtual ~MultiChannelNonOverlapCelListener();

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

  void addProbeMask(const std::string &name, const std::vector<bool> &mask);
  
  /**
   * Loop through the probesets provided and calculate a contrast value
   * for each one using the median of PM probes for A allele and B
   * allele. (CES constrast space)
   * @param cel - Cel File to get data from.
   * @param chrXProbeSets - Probesets to process (should be chrX non-pseudo autosomal).
   * @param contrastValues - Contrast values vector to be filled in.
   */
  static void fillInContrastValues(int celIx,
                                   std::vector<ProbeListPacked> &probeSets,
                                   std::vector<double> &contrastValues,
                                   double k = 2);
  
  /**
   * Go thru GC contrast values and find the mean and stdev
   * go thru AT contrast values and count how many are within 2stdev of GC mean
   * calculate the fraction of non-overlap
   * @param ATcontrast Values
   * @param GCcontrast Values
   * @param GCstats vector to fill in with stats
   */
  double calcNonOverlap( const std::vector<double> &ATcontrastValues,
                         const std::vector<double> &GCcontrastValues,
                         CumulativeStats<double>  &GCstats);
  

  /**
   * Go thru one channel's contrast values and calculate the fraction 
   * of probesets above or under given threshold.
   * @param contrastValues for a given channel
   * @param threshold to decide boundary
   * @para greater count probeset greater than threshold or not
   */
  double calcCrossThreshold(const std::vector<double> &contrastValues,
                            double threshold, 
                            bool greater = true );
  
  
  /**
   * calculate average_log_diff and log_diff_score (average_log_diff/log_diff_stdev)
   * @para logdiffvalues for one group
   * @para avg average_log_diff
   * @para score log_diff_score
   * @para IQRS  inter quantile range
   **/
  void calcLogDiff (std::vector<double> &logdiffValues, 
                    double &avg, 
                    double &score,
                    double &IQRS);
  
  /**
   * calculate average_log_diff and log_diff_score (average_log_diff/log_diff_stdev)
   * @para logdiffvalues for one group
   * @para score log_diff_score
   **/
  void calcLogDiff (std::vector<double> &logdiffValues, 
                    double &score);
  
  /**
   * Go thru GC contrast values and find the median and IQR
   * go thru AT contrast values and find the median and IQR
   * calculate the abs(GC_median - AT_median) / average_IQR 
   * @param ATcontrast Values
   * @param GCcontrast Values
   */
  double calcMedianDiff(const std::vector<double> &ATcontrastValues,
                        const std::vector<double> &GCcontrastValues,
                        CumulativeStats<double> &ATstats,
                        CumulativeStats<double> &GCstats);
  
  bool summarizeAllele(ProbeSet* pSet,
                       vector<double> &summaryValues,
                       SparseMart& dMart,
                       std::vector<ChipStream *>* iTrans,
                       PmAdjuster& pmAdjust,
                       QuantMethod* quantMethod,
                       bool lowPrecision);
  
  void GimmeContrast(std::vector<ProbeListPacked>& probeSets,
                     std::vector<double>& ATcontrastValues, 
                     std::vector<double>& GCcontrastValues, 
                     std::vector<double> &ATLogdiffValues, 
                     std::vector<double> &GCLogDiffValues, 
					 std::vector<double> &ATSaturationValues,
					 std::vector<double> &GCSaturationValues,
                     double k);

  void readSampleAlleleMap(const std::string& filename,
                           std::map<std::string,std::vector<int> >& sampleAlleleMap);

protected:
  // fill in the probe set with the allele values
  bool fillInAlleleProbeSets(const ProbeSet &gtPs, ProbeSet &aAllele, ProbeSet &bAllele); 
  
  affymetrix_fusion_io::FusionCELData* m_cel;

  /// Vector of probesets to compute HiLo on
  /// Pointer memory is owned elsewhere...
  std::vector<ProbeListPacked> m_ProbeSets;
  /// Name of the cel files that have been called.
  QCProbesetOptions m_psOpts;
  ChipLayout* m_layout;
  std::vector<std::string> m_CelNames;
  /// Our K parameter for the contrast centers transformation.
  double m_K;
  /// Label for the chip summary output
  std::string m_Label;
  /// probesetOptions
  
  /// The individual identifiers for particular subsets of probes.
  std::vector<std::string>  m_MaskNames;
  /// Subsets of probes in the form of bit masks.
  std::vector<std::vector<bool> > m_ProbesetMasks;
  
  ///** clear a probe set*/
  //void clearProbeSet(ProbeSet& ps); 
};

#endif /* MULTICHANNELNONOVERLAPCELLISTENER_H */
