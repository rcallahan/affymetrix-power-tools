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
 * @file   QuantMethodReport.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 12:09:42 2005
 * 
 * @brief  Class for reporting results of quantification methods.
 */
#ifndef _QUANTMETHODREPORT_H_
#define _QUANTMETHODREPORT_H_

//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/TsvReport.h"
//
#include "util/Util.h"
//
#include <cstdlib>
#include <iostream>
//

/**
 *   Class for reporting results of quantification methods. Current idea is that
 *   the quantification methods will have at least one matching reporter to form
 *   computation/report pairs.
 */
class QuantMethodReport : public affx::TsvReport {

public:

  /**
   * Virtual destructor for a virtual class.
   */
  virtual ~QuantMethodReport();

  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param qMethod - Quantification method to be used.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool prepare(QuantMethod &qMethod, const IntensityMart &iMart) = 0;

  /** 
   * After every probeset computation this function is called an is an opportunity
   * to query the quantification method for results, residuals, etc.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool report(ProbeSetGroup &psGroup, 
                      QuantMethod &qMethod, 
                      const IntensityMart &iMart, 
                      std::vector<ChipStream *> &iTrans, 
                      PmAdjuster &pmAdjust) = 0;

  /** 
   * If a probeset fails to compute for whatever reason then this method is
   * called rather than the normal report call above. By default does nothing.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool reportFailure(ProbeSetGroup &psGroup, 
                             QuantMethod &qMethod,
                             const IntensityMart &iMart, 
                             std::vector<ChipStream *> &iTrans, 
                             PmAdjuster &pmAdjust);


  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * 
   * @param qMethod - Quantification method that was used.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool finish(QuantMethod &qMethod) = 0;
    
  /** 
   * @brief Print a message out to the stream.
   * @param s - Message to be printed.
   * @param comment - prefix this with a comment character?
   */
  // virtual void out(const std::string &s, bool comment=false) {}

  /** 
   * @brief Print out a comment delimiter for separation.
   */
  // virtual void commentDelim() {}

  /** 
   * @brief Print out the header for reports.
   * @param s - List of column names (usually file names.)
   * @param timeStr - Time string.
   */
  // virtual void printHeader(const std::vector<std::string> &s, const std::string& timeStr) {}

  /** 
   * @brief Finish text header.
   */
  //virtual void finishTextHeader() {}

  virtual void addStdHeaders(QuantMethodReport *qReport,
                             const std::string& execGuid, 
                             const std::string& reportGuid,
                             const std::string& timeStr,
                             const std::string& commandLine,
                             const std::string& execVersion,
                             const AnalysisInfo& info);
};

#endif /* _QUANTMETHODREPORT_H_ */
