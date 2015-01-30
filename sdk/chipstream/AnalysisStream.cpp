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
 * @file   AnalysisStream.cpp
 * @author Chuck Sugnet
 * @date   Wed Jan 11 14:49:39 2006
 * 
 * @brief  
 */

//
#include "chipstream/AnalysisStream.h"
//
#include "chipstream/QuantMethodExprReport.h"

using namespace std;

bool AnalysisStream:: doAnalysis(               ProbeSetGroup &psGroup,
                          			IntensityMart &iMart,
                          			bool doReport,
                                                bool alleleSummaryOnly) {
  bool success = true;
  if (m_QMethod->setUp(psGroup, iMart, m_CStreams, *m_PmAdjust)) {
    if (!alleleSummaryOnly) {
      m_QMethod->computeEstimate();
      if (doReport) {
        for (unsigned int i = 0; i < m_Reporters.size(); i++) {
          m_Reporters[i]->report(psGroup, *m_QMethod, iMart, m_CStreams, *m_PmAdjust);
        }
      }
    }
  }
  else {
    Verbose::out(5, "Warning setup failed for name: " + ToStr(psGroup.name));
    success = false;
  }
  if (!success && doReport) {
    for (unsigned int i = 0; i < m_Reporters.size(); i++) {
      m_Reporters[i]->reportFailure(psGroup, *m_QMethod, iMart, m_CStreams, *m_PmAdjust);
    }
  }
  return success;
}

/** 
 * Destructor.
 */
AnalysisStream::~AnalysisStream() {
  for (unsigned int i = 0; i < m_CStreams.size(); i++) {
    Freez(m_CStreams[i]);
  }
  for (unsigned int i = 0; i < m_CSStages.size(); i++) {
    cleanupStage(i);
  }
  delete m_PmAdjust;
  delete m_QMethod;
  for (unsigned int i = 0; i < m_Reporters.size(); i++) {
    delete m_Reporters[i];
  }
}

/** 
 * Set the reporter for an analysis stream.
 *
 * @param outDir - directory to output things to.
 * @param precision - number of digits after decimal place.
 * @param doResiduals - output residual errors.
 * @param doFeatureEffects - output feature responses.
 */
QuantMethodExprReport*
AnalysisStream::addTextReporter(const std::string& outDir,
                                     int precision,
                                     bool doResiduals,
                                     bool doFeatureEffects) {

  return addTextReporter3(outDir,
                          precision,precision,precision,
                          doResiduals,
                          doFeatureEffects);
}

QuantMethodExprReport*
AnalysisStream::addTextReporter3(const std::string& outDir,
                                 int precision_sum,
                                 int precision_feff,
                                 int precision_res,
                                 bool doResiduals,
                                 bool doFeatureEffects) {

  QuantMethodExprReport* qReport=new QuantMethodExprReport(m_Info.getNumCols());
  qReport->m_qMethod=m_QMethod;
  //
  qReport->setDirPath(outDir);
  qReport->setFilename("AddTextReporter-DEBUG"); // debugging not used
  qReport->setFileprefix(m_Name);
  qReport->setIsHeaderBuffer(1);
  //
  qReport->m_DoSummary=true;
  qReport->m_DoFeatureEffects=doFeatureEffects;
  qReport->m_DoResiduals=doResiduals;
  // overall precision
  qReport->setPrecision(precision_sum);
  // set the precisions for each of the reports
  qReport->m_summary_tsv_precision=precision_sum;
  qReport->m_feffects_tsv_precision=precision_feff;
  qReport->m_residuals_tsv_precision=precision_res;
  //
  qReport->setFormat(affx::TsvReport::FMT_TSV);
  //
  addReporter(qReport);
  //
  return qReport;
}

/**
 * Write out the initial text header for a quantification method.
 * 
 * @param execGuid - Job unique id.
 * @param timeStr - When job was run.
 * @param commandLine - Command line specified to program.
 * @param colNames - Names of columns (usually cel file names).
 * @param execVersion - Version, cvs id string.
 */
void AnalysisStream::addStdHeaders(const std::string& execGuid, 
                                   const std::string& timeStr,
                                   const std::string& commandLine,
                                   // const std::vector<std::string> &colNames,
                                   const std::string& execVersion,
                                   const AnalysisInfo& info) {
  // printf("### AnalysisStream::addStdHeaders()\n"); // debug
  //
  string qMethodGuid = affxutil::Guid::GenerateNewGuid();
  for (unsigned int i = 0; i < m_Reporters.size(); i++) {
    QuantMethodReport *qReport = m_Reporters[i];
    qReport->addStdHeaders(qReport, 
                           execGuid,
                           qMethodGuid, 
                           timeStr,
                           commandLine,
                           // colNames,
                           execVersion, 
                           info);
  }
}

/**
 * Write out the initial text header for a quantification method.
 * 
 * @param layout - Chip information.
 * @param execGuid - Job unique id.
 * @param timeStr - When job was run.
 * @param commandLine - Command line specified to program.
 * @param colNames - Names of columns (usually cel file names).
 * @param execVersion - Version, cvs id string.
 */
void AnalysisStream::addStdHeaders(ChipLayout &layout, 
                                   const std::string& execGuid, 
                                   const std::string& timeStr,
                                   const std::string& commandLine,
                                   // const std::vector<std::string> &colNames,
                                   const std::string& execVersion,
                                   const AnalysisInfo &info) {
  // printf("### AnalysisStream::addStdHeaders()\n"); // debug
  //
  string qMethodGuid = affxutil::Guid::GenerateNewGuid();
  for (unsigned int i = 0; i < m_Reporters.size(); i++) {
    QuantMethodReport *qReport = m_Reporters[i];
    qReport->addStdHeaders(qReport, 
                           execGuid,
                           qMethodGuid, 
                           timeStr,
                           commandLine,
                           // colNames,
                           execVersion, info);
  }
}
