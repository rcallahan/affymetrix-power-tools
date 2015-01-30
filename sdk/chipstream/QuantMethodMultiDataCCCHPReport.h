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
 * @file   QuantMethodMultiDataCCCHPReport.h
 * @author David Le
 * @date   Wed Mar 15 12:46:53 2006
 * 
 * @brief Reporter for genotyping probe sets that outputs chp files. There is
 * one CC CHP file produced for each cel file and the order of the CC CHP file is
 * guarenteed to be the same as the cdf file. It is very important that the
 * number and order of the probeset results is the exact same as the CDF
 * file. There have been a number of high impact bugs on this front and it is
 * important to be very careful.
 */

#ifndef QUANTMETHODGTYPECCCHPREPORT_H
#define QUANTMETHODGTYPECCCHPREPORT_H

//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/ChipSummary.h"
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethodReport.h"
//
#include "calvin_files/writers/src/CalvinCHPMultiDataFileBufferWriter.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <cfloat>
#include <iostream>
//

class QuantMethodMultiDataCCCHPReport : public QuantMethodReport {
public: 
  /** Constructor. */
QuantMethodMultiDataCCCHPReport(const std::string& prefix)
  {
      m_Prefix = prefix;
      m_ProcessedProbeSetCount = 0;
      m_ProcessedExprProbeSetCount = 0;
      m_ProcessedGTypeProbeSetCount = 0;
      nGenotype = 0;
      nExpression = 0;
      nReporting = 0;
      genoTypeOnly = true;
      m_output_chip_view_md5sum=true;
      m_ProbeSetsToReport = NULL;
  }

  /** Destructor. */
  ~QuantMethodMultiDataCCCHPReport() { }

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
  
  /** The list of types to load */
  void SetGenoTypeOnly(bool genoOnly) { genoTypeOnly = genoOnly; }

  /**
   * Register a Chip Summary interface instance to pull metrics from
   */
  void registerChipSummary(ChipSummary *summary) {
      m_ChipSummaries.push_back(summary);
  }

  virtual void addStdHeaders(QuantMethodReport *qReport,
                             const std::string& execGuid, 
                             const std::string& reportGuid,
                             const std::string& timeStr,
                             const std::string& commandLine,
                             const std::string& execVersion,
                             const AnalysisInfo& info);

  void registerProbeSetsToReport(std::set<const char *, Util::ltstr> *probeSetsToReport) {
      m_ProbeSetsToReport = probeSetsToReport;
  }


private:
  void removeAllChps();
  void removeTmpChps();
  /** Do a check that the original id we wrote out and the one that
      actually came with the probeset group match. */
  void checkCurrentId(ProbeSetGroup &psGroup) {
      if(ToStr(psGroup.probeSets[0]->name) != ToStr(m_Info.m_ProbesetNames[m_ProcessedProbeSetCount])) {
          Err::errAbort(ToStr("QuantMethodMultiDataCCCHPReport::report() - Expecting probeset: ") +
                  m_Info.m_ProbesetNames[m_ProcessedProbeSetCount] + ToStr(" got: ") +
                  ToStr(psGroup.probeSets[0]->name) + " at index: " + ToStr(m_ProcessedProbeSetCount));
      }
    }

  /** The number of expression probe sets */
  int nExpression;
  
  /** The number of genotyping probe sets */
  int nGenotype;

  /** The number of probe sets being reported*/
  int nReporting;

  /** The list of types to load */
  bool genoTypeOnly;

  /** Swap out the .cel for .chp */
  void setupFileNames(const IntensityMart &iMart);

  /// Our prefix to filenames, often a path.
  std::string m_Prefix;

  /// Our algorithm name. Will be appended before .chp, i.e. mycel.cel-> mycel.algName.chp
  std::string m_AlgName;

  /// How many probesets have we seen this far?
  int m_ProcessedProbeSetCount;
  int m_ProcessedExprProbeSetCount;
  int m_ProcessedGTypeProbeSetCount;
  /// Names of all of our cel files
  std::vector<std::string> m_CELFileNames;
  /// Names of all of our chp files, one for each cel file.
  std::vector<std::string> m_CHPFileNames;
  /// Names of all of our chp files used by buffer writer
  std::vector<std::string> m_TmpChpFiles;
  /// All the infomation about our algorithm, etc.
  AnalysisInfo m_Info;
  /// Use for writing signals to a buffer. When the buffer is full, write it to CHP files.
  affymetrix_calvin_io::CHPMultiDataFileBufferWriter m_GenotypeEntryBufferWriter;
  /// Various Chip Summary Data Sources
  std::vector<ChipSummary *> m_ChipSummaries;
  /// Vector of probeset results which did not fail
  /// this is an index into m_Info.m_ProbesetNames
  std::vector<int> m_GoodProbesetsExpr;
  std::vector<int> m_GoodProbesetsGType;
  /// These vectors parallel the ones above, but are indexes
  /// into the specific data type section in chp file
  std::vector<int> m_GoodProbesetsExprIndex;
  std::vector<int> m_GoodProbesetsGTypeIndex;
  /// Vector of GUIDs for the respective CEL files
  std::vector<std::string> m_celGuids;

  /// put the "chip-view-md5sum" in the header?
  /// this is because we need to supress it for legacy sw. *sigh*
public:
  bool m_output_chip_view_md5sum;
  std::set<const char *, Util::ltstr> *m_ProbeSetsToReport;

};

#endif /* QUANTMETHODGTYPECCCHPREPORT_H */
