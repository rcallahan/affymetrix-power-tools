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
 * @file   QuantMethodRunReport.h
 * @author Chuck Sugnet
 * @date   Fri Jul 28 12:13:55 2006
 * 
 * @brief Summarize and output basic statistics on the performance of
 * a particular summarization or detection method across a number of
 * cel files.
 */
#ifndef _QUANTMETHODRUNREPORT_H_
#define _QUANTMETHODRUNREPORT_H_


//
#include "chipstream/ChipSummary.h"
#include "chipstream/CumulativeStats.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodReport.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Util.h"
//
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <set>


/**
 * Summarize and output basic statistics on the performance of a
 * particular summarization or detection method across a number of cel
 * files. All reports have statistics for raw probe intensity across a
 * chip. Detection based methods output the percentage above a certain
 * score. Summarization methods output the average and stdev for
 * summarized values, RLE (relative log expression), and MAD (median
 * absolute residual) when method provides residuals.
 * 
 */
class QuantMethodRunReport : public QuantMethodReport {

public :

  /** Constructor. */
  QuantMethodRunReport(const std::vector<std::string> &chipNames);

  /** Override the default file name **/
//  void setFileName(const std::string &filename) {
//      m_FileName = filename;
//  }

  /** 
   * Get set up for a run of reporting probesets. 
   * @param qMethod - Quantification method to be used.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * @return true if success, false otherwise.
   */
  bool prepare(QuantMethod &qMethod, const IntensityMart &iMart);

  /** 
   * After every probeset computation this function is called an is an opportunity
   * to query the quantification method for results, residuals, etc.
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * @param layout - Where the probesets, probes, etc are on the chip.
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
  bool finish(QuantMethod &qMethod) {
      return finish();
  }

  /** 
   * Finish outputting results and clean up.
   * @return true if success, false otherwise.
   */
  bool finish();

  /**
   * Register a Chip Summary interface instance to pull metrics from
   */
  void registerChipSummary(ChipSummary *summary);

  /**
   * Override base class to catch meta info and comments from engine
   */
  // this is sooo bougus
  // void out(const std::string &s, bool comment);

private :
  /// File name for printing report.
  //std::string m_FileName;
  /// Name for each of our chips.
  std::vector<std::string> m_ChipNames;
  /// Various Chip Summary Data Sources
  std::vector<ChipSummary *> m_ChipSummaries;
  /// Keys in order for header. -- META INFO
  //std::vector<std::string> m_HeaderKeys;
  /// Values matching the above keys in same order as above. -- META INFO
  //std::vector<std::string> m_HeaderVals;
  /// Other descriptive string and comments to put in header. -- META INFO
  //std::vector<std::string> m_OtherLines;
};

#endif /* _QUANTMETHODRUNREPORT_H_ */
