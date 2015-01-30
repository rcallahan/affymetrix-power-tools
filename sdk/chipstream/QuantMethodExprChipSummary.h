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
 * @file   QuantMethodExprChipSummary.h
 * @author Alan Williams
 * 
 * @brief Summarize basic statistics on the performance of
 * a particular summarization or detection method across a number of
 * cel files.
 */
#ifndef _QUANTMETHODEXPRCHIPSUMMARY_H_
#define _QUANTMETHODEXPRCHIPSUMMARY_H_

//
#include "chipstream/ChipSummary.h"
#include "chipstream/CumulativeStats.h"
#include "chipstream/MetaProbeset.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantMethodExprReport.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Util.h"
//
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <set>
//

/**
 * Summarize basic statistics on the performance of a
 * particular summarization or detection method across a number of cel
 * files. Detection based methods output the percentage above a certain
 * score. Summarization methods output the average and stdev for
 * summarized values, RLE (relative log expression), and MAD (median
 * absolute residual) when method provides residuals.
 * 
 */
class QuantMethodExprChipSummary : public QuantMethodExprReport, public ChipSummary {

public :

  /** Constructor. */
  QuantMethodExprChipSummary(const std::vector<std::string> &chipNames,
                             bool doThreshold, double minThreshold, double maxThreshold,
                             std::string qccFile,
                             std::vector<MetaProbeset *> &metaSets);

  /** 
   * Get set up for a run of reporting probesets. 
   * @param qMethod - Quantification method to be used.
   * @return true if success, false otherwise.
   */
  bool prepare(QuantMethod &qMethod, const IntensityMart &iMart);

  /** 
   * After every probeset computation this function is called an is an opportunity
   * to query the quantification method for results, residuals, etc.
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * @return true if success, false otherwise.
   */  
  bool report(ProbeSetGroup &psGroup, QuantMethod &qMethod,
              const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
              PmAdjuster &pmAdjust);

  /** 
   * Finish outputting results and clean up.
   * @param qMethod - Quantification method that was used.
   * @return true if success, false otherwise.
   */
  bool finish(QuantMethod &qMethod);

private :
  /** 
   * Load up the groups from a simple tab delimited text file. Need to
   * have a column called 'group_name' and then at least one column
   * named 'probeset_id' or 'probeset_name'. If groups called
   * 'pos_control' and 'neg_control' are both supplied then the AUC
   * from the ROC is estimated and reported as well. Looks something
   * like:
   *
   * \verbatim
     group_name   probeset_name
     pos_control  AFFX-123
     pos_control  AFFX-456
     neg_control  AFFX-789
     \endverbatim
   * @param fileName - path to text file to read from.
   * @param doThreshold - Should we be doing thresholding (detection based)?
   * @param minThreshold - Minimum value to pass threshold.
   * @param maxThreshold - Maximum value to pass threshold.
   */
  void readGroupsFile(const std::string &fileName, bool doThreshold, 
                      double minThreshold, double maxThreshold,
                      std::set<std::string> &psToRun);

  /**
   * @brief Class for function predicate, to sort pairs by the first coordinate,
   * or the second coordinate if the first coordinates are equal.
   */
  class CoordSort {
  public:
    bool operator() (const std::pair<float, float>* p1, const std::pair<float, float>* p2) const;
  };

  /**
   * Keep track of intensity, RLE and MAD statistics for probesets that
   * are provided iteratively via the report() function.
   */
  class ProbeSetGroupStats {
    
  public:
    
    /** 
     * Constructor. If the groupNames is empty then
     * all probesets are included in statistics.
     */
    ProbeSetGroupStats(const std::string &name, uint32_t numChips, 
                       bool cacheData, bool doThreshold,
                       double minThreshold, 
                       double maxThreshold,
                       const std::set<std::string> &groupNames);
    /** 
     * Check to see if we are including a particular ProbeSetGroup in
     * our statistics.
     * @param psGroup - ProbeSetGroup to check for inclusion.
     * @return bool - true if should be included, false otherwise.
     */
    bool inGroup(ProbeSetGroup &psGroup);
    
    /** 
     * Given a probeset group, extract and keep track of summary data
     * statistics.
     * @param qMethod - The quantification method to extract statistics from.
     */
    void report(ProbeSetGroup &psGroup, QuantExprMethod &qMethod);

    std::string m_Name; ///< Symbolic name for this group of probe sets.
    bool m_DoThreshold;  ///< Should we just do counts of things that pass threshold? 
    double m_MinThreshold; ///< What is the minimum threshold? Anything >= to this is accepted.
    double m_MaxThreshold; ///< What is the maximum threshold? Anything <= to this is accepted.
    uint32_t m_PsCount; ///< How many probeset groups have we seen?
    uint32_t m_AtomCount; ///< How many pm probes have been used?
    bool m_CacheData; ///< Are we keeping around the data for later examination?
    std::set<std::string> m_ValidNames; ///< What names should be included in this group.
    std::vector<CumulativeStats<double> > m_Summaries; ///< Summary for each chip.
    std::vector<CumulativeStats<double> > m_RLEs; ///< Relative log expression for each chip to median.
    std::vector<CumulativeStats<double> > m_MADs; ///< Median absolute residual for each chip.
    std::vector<std::vector<float> > m_Data; ///< Data cache for each chip, if requested.
  }; /* end of ProbeSetGroupStats class */

  /** 
   * Add a group that we wish to keep statistics for. If a probeset's name is in
   * the name set the results for that probeset will be added to the summary
   * stats. If both the name set is empty it acts as a wildcard and
   * statistics are kept for every probeset.
   * 
   * @param name - Reference name for this group.
   * @param cacheData - Should we keep the data in RAM as we go? (i.e. for ROC curve).
   * @param names - Set containing names of probesets to be kept.
   * @param doThreshold - Should we be doing thresholding (detection based)?
   * @param minThreshold - Minimum value to pass threshold.
   * @param maxThreshold - Maximum value to pass threshold.
   * @param isPosControl - Should this group of stats be considered a positive control
   * @param isNegControl - Should this group of stats be considered a negative control
   */
  void addGroupStat(const std::string &name, 
                    bool cacheData,
                    const std::set<std::string> &names,
                    bool doThreshold, double minThresh, double maxThresh,
                    bool isPosControl, bool isNegControl);
  
  /** 
   * Calculate the area under the ROC curve generated by the scores
   * for the positive and negative vectors. This code comes directly from
   * Peter's translation of the original exact perl code
   * @param pos - Vector of scores for positive examples.
   * @param neg - Vector of scores for negative examples.
   * @return - Estimate of area under the ROC curve.
   */
  double calcAuc(std::vector<float> &pos, std::vector<float> &neg, bool increasing = true);

private :
  /// Name for each of our chips.
  std::vector<std::string> m_ChipNames;
  /// Index of positive and negative controls in stats if they exist.
  int32_t m_posControlIx, m_negControlIx;
  /// Are we going to process pos and negative controls?
  bool m_HaveNegControls, m_HavePosControls;
  /// Our sets of probesets that we are keeping statistics for.
  std::vector<ProbeSetGroupStats> m_Stats;
  /// Is this run report using thresholds (i.e. detection)?
  bool m_DoThreshold;
  /// Map probeset IDs to probeset names for those probesets which go into summary stats
  std::map<std::string, std::string > m_QuantInHeader;
  /// Map probeset IDs/names to what group for those going into header
  std::map<std::string, std::string > m_QuantInHeaderGroup;
  /// Order of qc probesets in header
  std::vector<std::string> m_QcPsNameVec;
  /// Signal estimates for QC probesets to report in header
  std::vector<std::vector<double> > m_QcPsResults;
  /// Which QC probesets for the header did we actually process
  std::vector<std::string> m_QcPsSeen;
};

#endif /* _QUANTMETHODEXPRCHIPSUMMARY_H_ */
