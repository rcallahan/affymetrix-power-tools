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
 * @file   AnalysisStreamExpSpectSel.h
 * @author Chuck Sugnet
 * @date   Sun Dec 24 14:17:23 2006
 * 
 * @brief  AnalysisStream that selects features to cluster based on spectral clustering.
 */

#ifndef _ANALYSISSTREAMEXPSPECTSEL_H_
#define _ANALYSISSTREAMEXPSPECTSEL_H_

//
#include "chipstream/AnalysisStreamExpPcaSel.h"
//
#include "algorithm/spectclust/SpectClust.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <fstream>
#include <iostream>
#include <vector>
//

class AnalysisStreamExpSpectSel : public AnalysisStreamExpression {

public: 
  enum DistMetric {
    Corr,
    Angle,
    GaussianRadial
  };
  
  enum DistMetric stringToDist(const char *dist) {
    if(Util::sameString("corr", dist)) {
      return Corr;
    }
    else if(Util::sameString("angle", dist)) {
      return Angle;
    }
    else if(Util::sameString("gauss-radial", dist)) {
      return GaussianRadial;
    }
    Err::errAbort("Don't recognize DistMetric of type: '" + ToStr(dist) + "'");
    return Corr; // for compiler...
  }

  ~AnalysisStreamExpSpectSel() {
    Fs::carefulClose(m_Log);
  }

  AnalysisStreamExpSpectSel(bool doDebug, bool doFullEigen, int maxEigIter, int hardMin, 
                            double minProportion, const std::string &cutVal, bool logData,
                            double sigma, const std::string &metric, bool normDist, bool doRatio,
                            double margin, const std::string &infoFilter) {
    setupSelfDoc(*this);
    m_DoLog = logData;
    m_Debug = doDebug;
    m_HardMinimum = hardMin;
    m_HardProportion = minProportion;
    m_HardMaxProportion = 1;
    m_MaxEigIterations = maxEigIter;
    m_FullEigen = doFullEigen;
    m_Sigma = sigma;
    m_NormDist = normDist;
    m_Metric = stringToDist(metric.c_str());
    m_InfoFilter = AnalysisStreamExpPcaSel::stringToFilter(infoFilter.c_str());
    m_DoRatio = doRatio;
    m_Margin = margin;
    if(cutVal == "ncut") 
      m_Partition = DBL_MAX;
    else if(cutVal == "zero")
      m_Partition = 0;
    else 
      Err::errAbort("'" + cutVal + "' is not a valid parameter for 'cut-val' try 'ncut' or 'zero'");
    setOptValue("log2", m_DoLog);
    setOptValue("cut-val", cutVal);
    setOptValue("min-percent", ToStr(m_HardProportion));
    setOptValue("hard-min", ToStr(m_HardMinimum));
    setOptValue("debug", m_Debug);
    setOptValue("max-eig-iter", ToStr(m_MaxEigIterations));
    setOptValue("full-eigen", m_FullEigen);
    setOptValue("normdist", m_NormDist);
    setOptValue("metric", metric);
    setOptValue("ratio", m_DoRatio);
    setOptValue("margin", ToStr(m_Margin));
    setOptValue("info-criterion", infoFilter);
    if (m_FullEigen == true && m_DoLog == false && m_NormDist == true && metric == "gauss-radial") {
      Err::errAbort("The combination of options: full-eigen=true, log2=false, normdist=true, and metric=gauss-radial results in a process that does not converge.  Changing at least one of these parameters is advised.");
    }
  }

  bool doFeatureSelection(std::set<probeid_t> &goodIds, std::vector<double> &confVals,
                          Matrix &PM, std::vector<probeid_t> &probeIds, const char *psName);

  bool fillInSelectProbes(ProbeSetGroup &selectGroup, std::vector<ChipStream *> &cStream, 
                          IntensityMart &iMart, ProbeSetGroup &psGroup, std::vector<double> &confVals);

  /** 
   * Do the analysis for a particular group of probe sets. First
   * selecting probes that are closest to the principal component of
   * the data.
   * 
   * @param psGroup - Collection of probe sets to get probes from.
   * @param layout - How probes/probesets are laid out on chip.
   * @param iMart - Object containing raw data values for all chips.
   * @param doReport - Should the quantification report object be called?
 * @param alleleSummaryOnly - this is a parameter that makes sense in the base class AnalysisStream::doAnalysis, but not here.  Included here only to make the inheritance work.  Feel free to ignore.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool doAnalysis(ProbeSetGroup &psGroup,  
                          IntensityMart &iMart, 
                          bool doReport, 
                          bool alleleSummaryOnly = false);

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName("spect-select");
    doc.setDocDescription("Picks probes that are similar to each other based on spectral cluster and normalized cut.");
    doc.setDocOptions(getDefaultDocOptions());
  }

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf() { 
    SelfDoc doc;
    setupSelfDoc(doc);
    return doc;
  }
  
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> opts;
    SelfDoc::Opt doDebug = {"debug", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                            "Print out debugging data files (can be very large)."};
    opts.push_back(doDebug);
    SelfDoc::Opt doFullEigen = {"full-eigen", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                            "Explicitly calculate all eigen vectors rather than use power method to quickly get top N"};
    opts.push_back(doFullEigen);

    SelfDoc::Opt maxEigenIter = {"max-eig-iter", SelfDoc::Opt::Integer, "200", "200", "NA", "NA",
                                 "Maximum number of iterations to perform when using power method to get eigen vectors. Should be at least over 75."};
    opts.push_back(maxEigenIter);
    SelfDoc::Opt minPercent = {"min-percent", SelfDoc::Opt::Double, ".1", ".1", "0", "1",
                            "Minimum percentage of probes to use for summarization."};
    opts.push_back(minPercent);
    SelfDoc::Opt hardMin = {"hard-min", SelfDoc::Opt::Integer, "4", "4", "1", "NA",
                                 "Hard minimum on number of probes to use for summarization."};
    opts.push_back(hardMin);
    SelfDoc::Opt cutVal = {"cut-val", SelfDoc::Opt::String, "zero", "ncut", "NA", "NA",
                           "How to choose boundary for partition: 'ncut' for best normalized cut or 'zero' to just cut at 0."};
    opts.push_back(cutVal);
    SelfDoc::Opt log2Data = {"log2", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                             "Log 2 transform data before doing selection."};
    opts.push_back(log2Data);
    SelfDoc::Opt normDist = {"normdist", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                             "Should the distance matrix be normalized?"};
    opts.push_back(normDist);
    SelfDoc::Opt metric = {"metric", SelfDoc::Opt::String, "angle", "angle", "NA", "NA",
                             "What distance metric to use: 'angle', 'corr', or 'gauss-radial'"};
    opts.push_back(metric);
    SelfDoc::Opt doRatio = {"ratio", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                             "Should we use the ratio of a probe to its median?"};
    opts.push_back(doRatio);

    SelfDoc::Opt margin = {"margin", SelfDoc::Opt::Double, ".9", ".9", "0", "1",
                               "Percentage of 'good' probes to keep. Sometimes have lower quality probes on edge of cluster boundary."};
    opts.push_back(margin);
    SelfDoc::Opt infoFilter = {"info-criterion", SelfDoc::Opt::String, "aic", "aic", "NA", "NA",
                              "Should we use and information criter ('bic','aic','none') to determine if a strong enough signal was discovered to warrant feature selection?"};
    opts.push_back(infoFilter);
    return opts;
  }

  /** 
   * @brief This static function should be overridden by child classes
   * to return an object of the correct type initialized correctly
   * with the parameters in the string, string map. All objects
   * created this way should be deleted when finished using.
   * 
   * @param param - Map of key/value pairs to initialize the object.
   * 
   * @return Pointer toCreate object, this should be sub casted as necessary.
   */
  static SelfCreate *newObject(std::map<std::string,std::string> &param) {
    SelfDoc doc = explainSelf();
    bool doDebug = false;
    int maxEigIter = 200;
    bool fullEigen = false;
    std::string cutVal = "ncut";
    double minPercent = .1;
    int hardMin = 4;
    bool log2Data = false;
    bool normDist = true;
    std::string metric = "angle";
    bool doRatio = true;
    double margin = 1;
    std::string infoFilter = "bic";
    fillInValue(doDebug, "debug", param, doc);
    fillInValue(maxEigIter, "max-eig-iter", param, doc);
    fillInValue(fullEigen, "full-eigen", param, doc);
    fillInValue(hardMin, "hard-min", param, doc);
    fillInValue(minPercent, "min-percent", param, doc);
    fillInValue(cutVal, "cut-val", param, doc);
    fillInValue(log2Data, "log2", param, doc);
    fillInValue(normDist, "normdist", param, doc);
    fillInValue(metric, "metric", param, doc);
    fillInValue(doRatio, "ratio", param, doc);
    fillInValue(margin, "margin", param, doc);
    fillInValue(infoFilter, "info-criterion", param, doc);
    AnalysisStreamExpSpectSel *stream = new AnalysisStreamExpSpectSel(doDebug, fullEigen, maxEigIter, 
                                                                      hardMin, minPercent, cutVal, log2Data, 
                                                                      .5, metric, normDist, doRatio, margin, infoFilter);
    return stream;
  }


  
  /// Should the data be logged before doing feature selection?
  bool m_DoLog;
  /// Hard cap for how few probes we will use.
  int m_HardMinimum;
  /// Hard cap on minimal proportion of probes to use
  double m_HardProportion;
  /// Hard cap on minimal proportion of probes to use
  double m_HardMaxProportion;
  /// When doing power method for calculating eigen values what is the max number of iterations?
  int m_MaxEigIterations;
  /// Calculate eigenvectors explicitly rather than using power method?
  bool m_FullEigen;
  /// What distance metric are we using?
  enum DistMetric m_Metric;
  /// Where to cut cluster boundaries, DBL_MAX indicates to do it dynamically otherwise 0 is only other valid entry.
  double m_Partition;
  /// File to output some metrics too.
  bool m_Debug;
  /// When doing exponentiation we have a parameter called sigma.
  double m_Sigma;
  /// Should we normalize the distance matrix?
  bool m_NormDist;
  /// Should we use the ratio of a probe to its median?
  bool m_DoRatio;
  /// What percentage of probes that are "good" should we keep. Get a few lower quality hanger ons
  double m_Margin;
  /// Which filter should we use to determine if it worth it to do feature selection?
  enum AnalysisStreamExpPcaSel::InfoFilter m_InfoFilter;
  
  std::ofstream m_Log;
  std::ofstream m_AllProbes;
  std::ofstream m_UsedProbes;
  std::ofstream m_ProbeWeights;
};

#endif /* _ANALYSISSTREAMEXPSPECTSEL_H_ */
