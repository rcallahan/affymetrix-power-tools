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
 * @file   QuantMethodGTypeReport.h
 * @author Chuck Sugnet
 * @date   Fri Feb 24 15:15:14 2006
 * 
 * @brief  Reporter for outputting genotype results.
 * 
 */
#ifndef _QUANTMETHODGTYPEREPORT_H_
#define _QUANTMETHODGTYPEREPORT_H_

//
#include "chipstream/BioTypes.h"
#include "chipstream/QuantMethodReport.h"
#include "chipstream/TsvReport.h"
//
#include "util/Err.h"
#include "util/Util.h"
//
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
//

class QuantMethodGTypeReport : public QuantMethodReport {
public:

  /** 
   * Constructor takes the prefix that will be used for output files.
   * @param maxNameLength      How long is the longest probeset name for determining field size in hdf5
   * @param outputForcedCalls  Should the forced calls file be written
   * @param outputContext      Should the context of the snp be output i.e. near another snp
   */
  QuantMethodGTypeReport(bool outputForcedCalls=false, bool outputContext=false, bool outputProbabilities=false, int probFileSampleCount=500);

  /** Destructor. */
  ~QuantMethodGTypeReport();

  /** 
   * Use a compact version of the file5 format with structure of:
   * /calls - with genotypes as single byte
   * /confidences - with confidences as a 4 byte float
   * @param maxNameLength Maximum probeset name length
   */
  void setCompactFile5Format(int maxNameLength=TSVREPORT_PROBESET_STRLEN);

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
  bool report(ProbeSetGroup& psGroup,
              QuantMethod& qMethod, 
              const IntensityMart& iMart, 
              std::vector<ChipStream *>& iTrans, 
              PmAdjuster& pmAdjust);
    
  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * 
   * @param qMethod - Quantification method that was used.
   * 
   * @return true if success, false otherwise.
   */
  bool finish(QuantMethod &qMethod);
    
  void registerProbeSetsToReport(std::set<const char *, Util::ltstr> *probeSetsToReport) {
      m_ProbeSetsToReport = probeSetsToReport;
  }

private:
  /// Output stream for calls and confidences.
  affx::TsvReport m_CallsOutTsv;
  affx::TsvReport m_ConfsOutTsv;
  affx::TsvReport m_ContextOutTsv;
  affx::TsvReport m_ForcedCallsOutTsv;
  std::vector<affx::TsvReport> m_ProbabilitiesOutTsvVec;

  bool m_outputForcedCalls;
  bool m_outputContext;
  bool m_outputProbabilities;
  int m_probFileSampleCount;

  std::set<const char *, Util::ltstr> *m_ProbeSetsToReport;
  affx::File5_File *m_File5;
  affx::File5_Group* m_F5Group;
  int m_MaxNameLength;
  /// Should we output confidences as a float rather than a double and calls as a char
  bool m_DoCompact;
};

#endif /* _QUANTMETHODGTYPEREPORT_H_ */
