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
 * @file   QuantMethodGTypeCHPReport.h
 * @author Chuck Sugnet
 * @date   Wed Mar 15 12:46:53 2006
 * 
 * @brief Reporter for genotyping probe sets that outputs chp files. There is one
 * GCOS XDA CHP file produced for each cel file and the order of the GCOS XDA CHP
 * file is guaranteed to be the same as the cdf file. It is very important that
 * the number and order of the probeset results is the exact same as the CDF
 * file. There have been a number of high impact bugs on this front and it is
 * important to be very careful.
 */

#ifndef _QUANTMETHODGTYPECHPREPORT_H_
#define _QUANTMETHODGTYPECHPREPORT_H_

//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethodReport.h"
//
#include "file/CHPFileBufferWriter.h"
#include "file/CHPFileWriter.h"
#include "util/Util.h"
//
#include <cfloat>
#include <iostream>
//

class QuantMethodGTypeCHPReport : public QuantMethodReport {
public: 
  /** Constructor. */
  QuantMethodGTypeCHPReport(AnalysisInfo &chpInfo, 
                            const std::string& prefix, 
                            const std::string& algName) 
  {
      m_Info = chpInfo;
      m_AlgName = algName;
      m_Prefix = prefix;
      m_CurrentProbeSetCount = 0;
      Err::check(chpInfo.m_ProbesetNames.size() == chpInfo.m_NumProbeSets, 
          "Error: QuantMethodGTypeCHPReport::QuantMethodGTypeCHPReport() - m_NumProbeSets != m_ProbesetNames.size()");
  }

  /** Destructor. */
  ~QuantMethodGTypeCHPReport() { }

  /**
   * Set the CHP filenames to use for output. This overrides any prefix and allows
   * the caller to place the chp files anyplace they are desired.
   * @param fileNames - Vector of fileNames with full path to output
   * file, one for each cel file.
   */
  void setChpFileNames(std::vector<std::string> &fileNames) {
    m_CHPFileNames = fileNames;
  }

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
  bool report(ProbeSetGroup &psGroup, 
              QuantMethod &qMethod,
              const IntensityMart &iMart, 
              std::vector<ChipStream *> &iTrans, 
              PmAdjuster &pmAdjust);

  
  /** 
   * For genotyping CHP files GTYPE likes to have the median of the raw intensity
   * probes stuffed into the confidence field of the CHP files. I'm not sure why
   * but Richard C asked for it so here it is.
   * 
   * @param psGroup - Groupt of probes.
   * @param iMart - Raw intensity data for chips.
   * @param chipIx - Which chip we are taking the median for.
   * 
   * @return median of all raw intensities for the PM probes in psGroup.
   */
  float medianOfPmProbes(ProbeSetGroup &psGroup, const IntensityMart &iMart, int chipIx);

  /** 
   * If a probeset fails to compute for whatever reason then this method is
   * called rather than the normal report call above. By default does nothing.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
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
  affxchpwriter::CCHPFileBufferWriter m_GenotypeEntryBufferWriter;
  /// Names of all of our cel files
  std::vector<std::string> m_CELFileNames;
  /// Names of all of our chp files, one for each cel file.
  std::vector<std::string> m_CHPFileNames;
  std::vector<std::string> m_FilesForWriter;
  /// All the infomation about our algorithm, etc.
  AnalysisInfo m_Info;
};

#endif /* _QUANTMETHODGTYPECHPREPORT_H_ */
