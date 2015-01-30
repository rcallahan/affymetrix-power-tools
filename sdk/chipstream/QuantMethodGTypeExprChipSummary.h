////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
 * @file   QuantMethodGTypeExprChipSummary.h
 * @author Alan Williams
 * 
 * @brief  Class for accumulating some run statistics from a 
 *         genotyping run from the expression quant method
 * 
 */

#ifndef QUANTMETHODGTYPEEXPRCHIPSUMMARY_H
#define QUANTMETHODGTYPEEXPRCHIPSUMMARY_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/ChipSummary.h"
#include "chipstream/GenderCalls.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethodReport.h"
#include "chipstream/SummaryStats.h"
//
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <cfloat>
#include <iostream>
#include <vector>
//

class QuantMethodGTypeExprChipSummary : public QuantMethodReport, public ChipSummary {

public: 

  QuantMethodGTypeExprChipSummary(const std::string& analysisName);

  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param qMethod - Quantification method to be used.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool prepare(QuantMethod &qMethod, const IntensityMart &iMart) ;

  /** 
   * After every probeset computation this function is called an is an opportunity
   * to query the quantification method for results, residuals, etc.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool report(ProbeSetGroup &psGroup, QuantMethod &qMethod,
              const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
              PmAdjuster &pmAdjust);

  /** 
   * If a probeset fails to compute for whatever reason then this method is
   * called rather than the normal report call above. By default does nothing.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool reportFailure(ProbeSetGroup &psGroup, QuantMethod &qMethod,
                     const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                     PmAdjuster &pmAdjust);

  /** 
   * Add the summary statistics for the raw intensity of probes in this probeset group.
   *  
   * @param psGroup - Probes to be added.
   * @param iMart - Data to get intensity from.
   */
  virtual void addIntensityValues(ProbeSetGroup &psGroup, const IntensityMart &iMart);
  
  virtual void addAlleleSummary(QuantExprMethod &qMethod);

  virtual void addMadResiduals(QuantExprMethod &qMethod);

  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * @param qMethod - Quantification method that was used.
   * @return true if success, false otherwise.
   */
  virtual bool finish(QuantMethod &qMethod);

private:
  /// The name of our analysis.
  std::string m_AnalysisName;
  // These objects are for keeping track of various summarization based statistics.
  /// Raw chip intensity.
  SummaryStats m_IntensitySummary;
  /// Allele summarization
  SummaryStats m_AlleleSummary;
  /// Residuals summaries. median absolute deviation
  SummaryStats m_MadResidualSummary;
  /// Deviation of alleles from mean. (euclidean distance)
  SummaryStats m_AbsDeviation; 

};

#endif /* QUANTMETHODGTYPEEXPRCHIPSUMMARY_H */
