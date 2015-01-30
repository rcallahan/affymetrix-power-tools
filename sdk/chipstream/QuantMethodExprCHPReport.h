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
 * @file   QuantMethodExprCHPReport.h
 * @author David Le
 * @date   Mon May 15 12:09:42 2006
 * 
 * @brief  Class for reporting results of quantification methods.
 */
#ifndef _QUANTMETHODEXPRCHPREPORT_H_
#define _QUANTMETHODEXPRCHPREPORT_H_

//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantMethodReport.h"
//
#include "file/CHPFileBufferWriter.h"
#include "file/CHPFileWriter.h"
#include "util/Util.h"
//
#include <cstdlib>
#include <iostream>
//

/**
 *   Class for reporting results of quantification methods.
 */
class QuantMethodExprCHPReport : public QuantMethodReport {

public:
  
  /** Constructor. */
  QuantMethodExprCHPReport(AnalysisInfo& chpInfo, 
                           const std::string& prefix, 
                           const std::string& algName);

  /**
   * Set the CHP filenames to use for output. This overrides any prefix and allows
   * the caller to place the chp files anyplace they are desired.
   * @param fileNames - Vector of fileNames with full path to output
   * file, one for each cel file.
   */
  void setChpFileNames(std::vector<std::string> &fileNames);

  /** Destructor. */
  ~QuantMethodExprCHPReport();

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
  
private:
  /** Swap out the .cel for .chp */
  void setupFileNames(const IntensityMart &iMart);

  /** Fill in a chp header file info. */
  void setupChpFile(affxchpwriter::CCHPFileWriter &chp, AnalysisInfo &info);

  /// Our prefix to filenames, often a path.
  std::string m_Prefix;

  /// Our algorithm name. Will be appended before .chp, i.e. mycel.cel-> mycel.algName.chp
  std::string m_AlgName;

  /// How many probesets have we seen this far?
  int m_CurrentProbeSetCount;
  /// Use for writing signals to a buffer. When the buffer is full, write it to CHP files.
  affxchpwriter::CCHPFileBufferWriter m_ExpressionEntryBufferWriter;
  /// Names of all of our cel files
  std::vector<std::string> m_CELFileNames;
  /// Names of all of our chp files, one for each cel file.
  std::vector<std::string> m_CHPFileNames;
  std::vector<std::string> m_FilesForWriter;
  /// All the infomation about our algorithm, etc.
  AnalysisInfo m_Info;
};

#endif /* _QUANTMETHODEXPRCHPREPORT_H_ */
