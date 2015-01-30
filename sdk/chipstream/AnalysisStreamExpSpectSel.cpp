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
 * @file   AnalysisStreamExpSpectSel.cpp
 * @author Chuck Sugnet
 * @date   Sun Dec 24 22:52:10 2006
 * 
 * @brief  
 * 
 * 
 */

//
#include "chipstream/AnalysisStreamExpSpectSel.h"
//
#include "algorithm/spectclust/SpectClust.h"
#include "util/Fs.h"
#include "util/Err.h"
//
#include <cfloat>

bool AnalysisStreamExpSpectSel::doFeatureSelection(std::set<probeid_t> &goodIds, std::vector<double> &confVals,
                                                   Matrix &PM, std::vector<probeid_t> &probeIds, const char *psName) {
  Matrix aicData;
  aicData = PM;
  if(m_DoRatio) {
    SpectClust::rowMedianDivide(PM);
  }
  PM = PM.t();
  if(m_DoLog) {
    AnalysisStreamExpPcaSel::log2Matrix(PM);
  }
  DistanceMetric *metric = NULL;
  if(m_Metric == AnalysisStreamExpSpectSel::Angle) {
    metric = new AngleMetric(PM);
  }
  else if(m_Metric == AnalysisStreamExpSpectSel::Corr) {
    metric = new CorrelationMetric();
  }
  else if(m_Metric == AnalysisStreamExpSpectSel::GaussianRadial) {
    metric = new GuassianRadial(m_Sigma);
  }
  else {
    Err::errAbort("Don't recognized metric of type: " + ToStr(m_Metric));
  }
  SymmetricMatrix D;
  Matrix EVec;

  std::vector<double> eVals;
  std::vector<int> clusters;
  SpectClust::fillInDistance(D, PM, *metric, false);
  delete metric;
  if(m_NormDist) {
    SpectClust::normalizeSum(D);
  }
  bool converged = false;
  if(m_FullEigen) {
    converged = SpectClust::findNLargestSymEvals(D, 2, eVals,EVec);
  }
  else {
    converged = SpectClust::findNLargestEvals(D, 2, eVals, EVec, m_MaxEigIterations);
  }
  if(m_ProbeWeights.is_open()) {
    for(int i = 0; i < probeIds.size(); i++) {
      m_ProbeWeights << psName << "\t" << probeIds[i] << "\t" << EVec.element(i,1) << endl;
    }
  }
  SpectClust::partitionClusters(D, EVec, eVals, 2, clusters, m_HardMinimum, m_Partition, m_Margin);
  SpectClust::orderIntraClusterCorr(PM, clusters, 2, confVals);

  // if we're doing filtering on information criteria then calculate the confidence values.
  if(m_InfoFilter != AnalysisStreamExpPcaSel::NoFilter) {
    vector<double> aicVals, bicVals;
    SpectClust::AICcRmaMedianCluster(aicData, 2, clusters, aicVals, bicVals, true);
    if(m_InfoFilter == AnalysisStreamExpPcaSel::BIC)
      confVals = bicVals;
    else if(m_InfoFilter == AnalysisStreamExpPcaSel::AIC)
      confVals = aicVals;
    else
      Err::errAbort("AnalysisStreamExpSpectSel::doFeatureSelection() - Don't recognize info filter: '" + ToStr(m_InfoFilter) + "'");
  }
  //  SpectClust::AICc(aicData, 2, clusters, aicVals, confVals);
  for(int i = 0; i < probeIds.size(); i++) {
    if(clusters[i] == 1) {
      goodIds.insert(probeIds[i]);
    }
  }
  return converged;
}
                                                 

bool AnalysisStreamExpSpectSel::fillInSelectProbes(ProbeSetGroup &selectGroup, std::vector<ChipStream *> &cStream, 
                                                 IntensityMart &iMart, ProbeSetGroup &psGroup, std::vector<double> &confVals) {
  int chipCount = iMart.getCelFileCount();
  int atomPmCount = psGroup.countPmProbes();
  int newAtomPmCount = atomPmCount;
  bool success = true;
  if(atomPmCount <= m_HardMinimum) {
    success = false;
  }
  else {
    vector<bool> probesToUse(atomPmCount, false);
    vector<probeid_t> probeIds(atomPmCount, -1); // Probe ids in order as they are in psGroup
    set<probeid_t> goodIds;
    Matrix PM(atomPmCount, chipCount); // rows are probes, columns are chips
    ///@todo need to use correct channel for multi channel expression probesets
    AnalysisStreamExpPcaSel::fillInPmData(PM, probeIds, psGroup, cStream, iMart, false);

    string tempName;
    tempName = psGroup.name;
    doFeatureSelection(goodIds, confVals, PM, probeIds, tempName.c_str());
    int minPossibleProbes = max(m_HardMinimum, (int)(atomPmCount * m_HardProportion));
    int maxPossibleProbes = (int)(atomPmCount * m_HardMaxProportion);
    bool infoGood = true;
    if(confVals[0] <= confVals[1] && m_InfoFilter != AnalysisStreamExpPcaSel::NoFilter) {
      infoGood = false;
      success = false;
    }
    if(minPossibleProbes >= goodIds.size() || goodIds.size() >= maxPossibleProbes) {
      fill(confVals.begin(), confVals.end(), 0);
      success = false;
    }
    else {
      AnalysisStreamExpPcaSel::makeSelectProbesetGroup(selectGroup, goodIds, psGroup);
      newAtomPmCount = selectGroup.countPmProbes();
    }
  }
  return success;
}

/** 
 * Do the analysis for a particular group of probe sets.
 * 
 * @param psGroup - Collection of probe sets to get probes from.
 * @param layout - How probes/probesets are laid out on chip.
 * @param iMart - Object containing raw data values for all chips.
 * @param doReport - Should the quantification report object be called?
 * @param alleleSummaryOnly - this is a parameter that makes sense in the base class AnalysisStream::doAnalysis, but not here.  Included here only to make the inheritance work.  Feel free to ignore.
 * 
 * @return true if success, false otherwise.
 */
bool AnalysisStreamExpSpectSel::doAnalysis(ProbeSetGroup &psGroup,
                                           IntensityMart &iMart, 
                                           bool doReport,
                                           bool alleleSummaryOnly) {
  if(m_Debug && !m_AllProbes.is_open()) {
    /// @todo use TsvFile
    Fs::mustOpenToWrite(m_AllProbes, Fs::join(m_OutPrefix,getName() + ".spect-select.data.txt"));
    Fs::mustOpenToWrite(m_UsedProbes,Fs::join(m_OutPrefix,getName() + ".spect-select.useddata.txt"));
    m_AllProbes << "probeset\tprobe";
    m_UsedProbes << "probeset\tprobe";

    std::vector<std::string> celFiles = iMart.getCelFileNames();
    for(int i = 0; i < celFiles.size(); i++) {
      m_AllProbes << "\t" << celFiles[i];
      m_UsedProbes << "\t" << celFiles[i];
    }
    m_AllProbes << endl;
    m_UsedProbes << endl;

    Fs::mustOpenToWrite(m_ProbeWeights,Fs::join(m_OutPrefix,getName() + ".spect-select.weights.txt"));
    m_ProbeWeights << "probeset\tprobe\tweight" << endl;
  }
  if(!m_Log.is_open()) {
    Fs::mustOpenToWrite(m_Log,Fs::join(m_OutPrefix,getName() + ".spect-select.report.txt"));
    /// @todo Use TsvFile
    m_Log << "probeset_id\ttotal\tused\tclust1_info\tclust2_info\torig_probes\tused_probes" << endl;
  }
  bool success = true;
  ProbeSetGroup psSelectGroup, *toUse = NULL;
  vector<double> confVals;
  confVals.push_back(0);
  confVals.push_back(0);
  if(fillInSelectProbes(psSelectGroup, m_CStreams, iMart, psGroup, confVals)) {
    toUse = &psSelectGroup;
  }
  else {
    toUse = &psGroup;
  }
  if(m_Log.is_open()) {
    AnalysisStreamExpPcaSel::reportProbesUsed(m_Log, *toUse, psGroup, confVals);
  }
  if(m_AllProbes.is_open()) {
    AnalysisStreamExpPcaSel::reportProbeLevelData(m_AllProbes, psGroup, iMart, m_CStreams);
    AnalysisStreamExpPcaSel::reportProbeLevelData(m_UsedProbes, *toUse, iMart, m_CStreams);
  }
  if(m_QMethod->setUp(*toUse, iMart, m_CStreams, *m_PmAdjust)) {
    m_QMethod->computeEstimate();
    if(doReport) {
      for(unsigned int i = 0; i < m_Reporters.size(); i++) {
        m_Reporters[i]->report(*toUse, *m_QMethod, iMart, m_CStreams, *m_PmAdjust);
      }
    }
  }
  else {
    Verbose::out(5, "Warning setup failed for name: " + ToStr(toUse->name));
    success = false;
  }
  if(!success && doReport) {
    for(unsigned int i = 0; i < m_Reporters.size(); i++) {
      m_Reporters[i]->reportFailure(*toUse, *m_QMethod, iMart, m_CStreams, *m_PmAdjust);
    }
  }
  return success;
}
