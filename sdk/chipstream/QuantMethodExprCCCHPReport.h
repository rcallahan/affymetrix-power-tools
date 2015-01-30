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
 * @file   QuantMethodExprCCCHPReport.h
 * @author David Le
 * @date   Mon May 15 12:09:42 2006
 * 
 * @brief  Class for reporting results of quantification methods.
 */
#ifndef _QUANTMETHODEXPRCCCHPREPORT_H_
#define _QUANTMETHODEXPRCCCHPREPORT_H_

//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/ChipSummary.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantMethodReport.h"
//
#include "calvin_files/writers/src/CalvinCHPQuantificationFileBufferWriter.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationFileWriter.h"
#include "util/Util.h"
//
#include <cstdlib>
#include <iostream>
//

/**
 *   Class for reporting results of quantification methods.
 */
class QuantMethodExprCCCHPReport : public QuantMethodReport {

public:
  
  /** Constructor. */
  QuantMethodExprCCCHPReport(AnalysisInfo& chpInfo,
                             const std::string& prefix,
                             const std::string algName);
  

  /** Destructor. */
  ~QuantMethodExprCCCHPReport();

  /**
   * Set the CHP filenames to use for output. This overrides any prefix and allows
   * the caller to place the chp files anyplace they are desired.
   * @param fileNames - Vector of fileNames with full path to output
   * file, one for each cel file.
   */
  void setChpFileNames(std::vector<std::string> &fileNames);

  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param qMethod - Quantification method to be used.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  bool prepare(QuantMethod &qMethod, const IntensityMart &iMart);

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
  bool report(ProbeSetGroup &psGroup, QuantMethod &qMethod,
              const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
              PmAdjuster &pmAdjust);

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
   bool reportFailure(ProbeSetGroup &psGroup, 
                      QuantMethod &qMethod, 
                      const IntensityMart &iMart, 
                      std::vector<ChipStream *> &iTrans, 
                      PmAdjuster &pmAdjust);

  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * @param qMethod - Quantification method that was used.
   * @return true if success, false otherwise.
   */
  bool finish(QuantMethod &qMethod);
  
  /**
   * Register a Chip Summary interface instance to pull metrics from
   */
  void registerChipSummary(ChipSummary *summary);

private:
  void removeAllChps();
  void removeTmpChps();
  /** Swap out the .cel for .chp */
  void setupFileNames(const IntensityMart &iMart);

  /** Do a check that the original id we wrote out and the one that
      actually came with the probeset group match. */
  void checkCurrentId(ProbeSetGroup &psGroup);

  /// Our prefix to filenames, often a path.
  std::string m_Prefix;

  /// Our algorithm name. Will be appended before .chp, i.e. mycel.cel-> mycel.algName.chp
  std::string m_AlgName;

  /// How many probesets have we seen this far?
  int m_CurrentProbeSetCount;
  /// Use for writing signals to a buffer. When the buffer is full, write it to CHP files.
  affymetrix_calvin_io::CHPQuantificationFileBufferWriter m_ExpressionQuantificationBufferWriter;
  /// Names of all of our cel files
  std::vector<std::string> m_CELFileNames;
  /// Names of all of our chp files, one for each cel file.
  std::vector<std::string> m_CHPFileNames;
  /// Names of all of our chp files used by buffer writer
  std::vector<std::string> m_TmpChpFiles;
  /// All the infomation about our algorithm, etc.
  AnalysisInfo m_Info;
  /// Various Chip Summary Data Sources
  std::vector<ChipSummary *> m_ChipSummaries;
  /// Vector of probeset results which did not fail
  /// this is an index into m_Info.m_ProbesetNames
  std::vector<int> m_GoodProbesets;
  /// Vector of GUIDs for the respective CEL files
  std::vector<std::string> m_celGuids;
};

#endif /* _QUANTMETHODEXPRCCCHPREPORT_H_ */
