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
 * @file   AnalysisStream.h
 * @author Chuck Sugnet
 * @date   Tue Jan 10 16:48:31 2006
 * 
 * @brief Object for managing chipstreams, perfect match adjusters and
 * summarization methods.
 */

#ifndef _ANALYSISSTREAM_H_
#define _ANALYSISSTREAM_H_

//
#include "chipstream/AnalysisInfo.h"
#include "chipstream/ChipStream.h"
#include "chipstream/ChipStreamDataTransform.h"
#include "chipstream/IntensityReader.h"
#include "chipstream/QuantMethod.h"
#include "chipstream/QuantMethodExprReport.h"
#include "chipstream/QuantMethodReport.h"
//
#include "util/Guid.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

/**
\page analysisstream Architecture: AnalysisStream Tutorial

The AnalysisStream class is meant to enable routine analysis workflows. The idea
is that a ChipStream list of data transformers are connected to a perfect match
adjuster (PmAdjuster) which then feeds into a summarization method (QuantMethod)
for each ProbeSet, or group of ProbeSets. These analyses can be specified by a
text description and then used relatively easily. Besides the AnalysisStream a
 object to specify the physical layout of the probes on the chip and
groups of probes is necessary. Also needed is an IntensityMart to supply the raw
microarray data from all of the chips in the analysis. A schematic of the overall AnalysisStream
process can be seen in Figure 1 below.

@image html analysisStreamDataProcessSnps.png "Figure 1: Example of an AnalysisStream for genotyping probesets. Raw data from cel files is normalized and then processed with the BRLMM algorithm in the QuantBrlmm class. Actual genotype calls and summary statistics are supplied by reporters that listen for results from Brlmm algorithm."

For expression GeneChips, the BRLMM quantification method in Figure 1 would be
replaced with another quantification method that summarizes the individual
probes to better estimate overall intensity attributed to the abundance of a
particular transcript in the cell. Figure 2 gives some examples for the
different types of ChipStream data modifications, PmAdjuster methods, and
expression QuantMethods for summarizing results.

@image html analysisstreamExamples.png "Figure 2: Examples for the different types of ChipStream data modifications, PmAdjuster methods, and expression QuantMethods for summarizing results."

An example of an AnalyisStream 
process for doing RMA analysis can be found in the
sdk/chipstream/example/simpleRmaExample.cpp file.

*/

/** Utility class for centralizing analysis pathway */
class AnalysisStream : public SelfDoc, public SelfCreate {

public:
  
  /** Constructor. */
  AnalysisStream() { 
    m_PmAdjust = NULL;
    m_QMethod = NULL;
    analysisGuid = affxutil::Guid::GenerateNewGuid();
  }

  /** 
   * Destructor.
   */
  virtual ~AnalysisStream();

  /**
   * Delete a stage as it is done processing
   */
  void cleanupStage(int i) {
    Freez(m_CSStages[i]);
  }


  /** 
   * Add a reporter that prints a QC report.
   * 
   * @param detectedP - P value threshold for detection.
   * @param startTime - Time of the start of the run.
   * @param QCOut - ofstream to write report to.
   */
  void addQCReporter(const double detectedP, const time_t startTime, std::ofstream &QCOut);

  /**
   * Add a reporter that prints to simple tab delimited file for an
   * analysis stream.
   *
   * @param outDir - directory to output things to.
   * @param precision - number of digits after decimal place.
   * @param doResiduals - output residual errors.
   * @param doFeatureEffects - output feature responses.
   */
  virtual QuantMethodExprReport* addTextReporter(const std::string& outDir,
                                int precision,
                                bool doResiduals,
                                bool doFeatureEffects);
  virtual QuantMethodExprReport* addTextReporter3(const std::string& outDir,
                                int precision_sum,
                                int precision_feff,
                                int precision_res,
                                bool doResiduals,
                                bool doFeatureEffects);

  /** 
   * Write out the initial text header for a quantification method.
   * 
   * @param execGuid - Job unique id.
   * @param timeStr - When job was run.
   * @param commandLine - Command line specified to program.
   * @param colNames - Names of columns (usually cel file names).
   * @param execVersion - Version, cvs id string.
   */
  virtual void addStdHeaders(const std::string& execGuid, 
                             const std::string& timeStr,
                             const std::string& commandLine,
                             const std::string& execVersion,
                             const AnalysisInfo& info);

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
  virtual void addStdHeaders(ChipLayout& layout,
                             const std::string& execGuid, 
                             const std::string& timeStr,
                             const std::string& commandLine,
                             const std::string& execVersion,
                             const AnalysisInfo& info);

  /** 
   * Connect a reporter to an analysis stream.
   * 
   * @param reporter - object to use for output.
   */
  virtual void addReporter(QuantMethodReport* reporter) {
    m_Reporters.push_back(reporter);
  }

  /** 
   * Set the object for adusting perfect match intensities.
   * @param pmAdjuster - object for specifying adustments.
   */
  virtual void setPmAdjuster(PmAdjuster *pmAdjuster) {
    m_PmAdjust = pmAdjuster;
  }

  /** 
   * Set the object for summarizing an individual probe set.
   * @param qMethod - Quantification object.
   */
  virtual void setQuantMethod(QuantMethod *qMethod) {
    m_QMethod = qMethod;
  }

  /** 
   * Add a new chipstream object to the end of the current
   * list of chipstream transformation objects.
   * 
   * @param cStream - New chipstream transformation object.
   */
  virtual void addChipStream(ChipStream *cStream) {
    /* Register as next in chipstream linked list for receiving data. */
    if(!m_CStreams.empty()) {
      m_CStreams[m_CStreams.size() - 1]->registerStream(cStream);
      cStream->registerParent(m_CStreams[m_CStreams.size() - 1]);
    }
    m_CStreams.push_back(cStream);
  }

  /** 
   * What should this particular analysis stream be called?
   * @param s - New name for analysis stream.
   */
  void setName(const std::string &s) {
    m_Name = s;
  }

  /** 
   * What is the name of this particular analysis stream?
   */
  std::string getName() {
    return m_Name;
  }

  /** 
   * Get the beginning of the chipstream list.
   * @return - First in the list of chipstream objects.
   */
  ChipStream *getChipStreamHead() {
    if(m_CStreams.empty())
      return NULL;
    return m_CStreams[0];
  }

  /**
   * Tell a CelReader what chipstream objects it needs to know about.
   * @param reader - 
   */
  virtual void registerChipStreamObjs(IntensityReader &reader) {
    if(getChipStreamHead() != NULL)
      reader.registerStream(getChipStreamHead());
  }

  /** 
   * Do the analysis for a particular group of probe sets.
   * 
   * @param psGroup - Collection of probe sets to get probes from.
   * @param iMart - Object containing raw data values for all chips.
   * @param doReport - Should the quantification report object be called?
   * @param alleleSummaryOnly - A hack to get apt-probeset-genotype to output Axiom .summary files without genotyping.  "People" decided that it was easier to do this than to get apt-probeset-summarize to read dual-channel CEL files.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool doAnalysis(ProbeSetGroup &psGroup, 
                          IntensityMart &iMart, 
                          bool doReport,
                          bool alleleSummaryOnly = false);
  /** 
   * Get the quantification method used. Useful for getting results after
   * doAnalysis() has been called.
   * @return - Quantification method for this stream.
   */
  virtual QuantMethod *getQuantMethod() {
    return m_QMethod;
  }

  virtual void prepare(const IntensityMart &iMart) {
    m_QMethod->prepare(iMart);
    for(unsigned int i = 0; i < m_Reporters.size(); i++) {
      m_Reporters[i]->prepare(*m_QMethod, iMart);
    }
  }
  
  virtual void finish() {
    m_QMethod->finish();
    for(unsigned int i = 0; i < m_Reporters.size(); i++) {
      m_Reporters[i]->finish(*m_QMethod);
    }
  }

  virtual std::vector<ChipStream *> *getChipStream() {
    return &m_CStreams;
  }

  virtual PmAdjuster *getPmAdjuster() {
    return m_PmAdjust;
  }

  std::string getGuid() {
    return analysisGuid;
  }
  
  /** 
   * Get a formal description of each object in the analysis stream
   * including chipstream objects, pm adjuster and quantification
   * method. The string description is of the form that SelfCreate
   * objects would use to create themselves.
   * @return - text description of objects and their state.
   */
  std::string getAnalysisSpec() { 
    std::string description;
    std::vector<ChipStream *>::iterator csIx;
    for(csIx = m_CStreams.begin(); csIx != m_CStreams.end(); ++csIx) {
      description += (*csIx)->getState() + ",";
    }
    description += m_PmAdjust->getState() + ",";
    description += m_QMethod->getState() + ",";
    description += getState();
    return description;
  }

  void setInfo(AnalysisInfo &info) {
    m_Info = info;
    Verbose::out(1,"Setting analysis info.");
    m_QMethod->setAnalysisInfo(info);
  }

  AnalysisInfo getInfo() const {
    return m_Info;
  }

  std::string getOutPrefix() const { 
    return m_OutPrefix; 
  }

  void setOutPrefix(const std::string &prefix) { 
    m_OutPrefix = prefix;
  }

  /**
   * Get the number of registered reporters.
   * @return The number of reporters.
   */
  int getReporterSize() { return (int) m_Reporters.size(); }

  /**
   * Get the reporter, given by an index.
   * @return The reporter.
   */
  QuantMethodReport * getReporter(int index) { return m_Reporters[index]; }

  /// Data transform stages
  std::vector<ChipStreamDataTransform *> m_CSStages;

protected: 

  /// Name of this analysis stream.
  std::string m_Name;
  /// Chip stream transformers for this analsyis (i.e bg sub)
  std::vector<ChipStream *> m_CStreams;

  /// Adjustment to use for this analysis.
  PmAdjuster *m_PmAdjust;
  /// Quantification method (i.e. plier)
  QuantMethod *m_QMethod;
  /// Report object for outputting results
  std::vector<QuantMethodReport *> m_Reporters;
  /// Pseudorandom identifier for this analysis 
 std::string analysisGuid;
  /// A collection of information about program and analyisis state.
  AnalysisInfo m_Info;
  /// Output directory prefix where sometimes diagnostic files go.
  std::string m_OutPrefix;
};

#endif /* _ANALYSISSTREAM_H_ */
