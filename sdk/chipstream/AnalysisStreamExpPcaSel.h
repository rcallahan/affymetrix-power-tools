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
 * @file   AnalysisStreamExpPcaSel.h
 * @author Chuck Sugnet
 * @date   Fri Nov  3 10:18:15 2006
 *
 * @brief An analysis stream that takes the expression data and does a pca
 * selection on the probes before passing them to the quantification
 * method. Inspired by Jim's Corimbia algorithm.
 *
 */
#ifndef _ANALYSISSTREAMEXPPCASEL_H_
#define _ANALYSISSTREAMEXPPCASEL_H_

//
#include "chipstream/AnalysisStreamExpression.h"
#include "chipstream/SketchQuantNormTran.h"
//
#include "newmat.h"
//
#include <iostream>
#include <set>
//

/**
 * @brief An analysis stream that takes the expression data and does a
 * pca selection on the probes before passing them to the
 * quantification method. Inspired by Jim's Corimbia algorithm.
 *
 */
class AnalysisStreamExpPcaSel : public AnalysisStreamExpression {

public:

  enum InfoFilter {
    NoFilter, // Don't do any filtering.
    BIC, // Bayesian information criteria, aka Schwartz criterion.
    AIC // Akaike information criterion
  };

  static enum InfoFilter stringToFilter(const std::string& filter);

  /** Constructor. */
  AnalysisStreamExpPcaSel(bool doLog, bool doCorr, bool doDebug,
                          const std::string &infoFilter,
                          int hardMin, double minProportion, bool qnormOnly);

  ~AnalysisStreamExpPcaSel();

  void setQuantNorm(SketchQuantNormTran *qnorm);

  /**
   * Tell a CelReader what chipstream objects it needs to know about.
   * @param CelReader - Cel file reader than needs to be notified
   * about chipstream objects.
   */
  virtual void registerChipStreamObjs(IntensityReader &reader);

  /**
   * Do the analysis for a particular group of probe sets. First
   * selecting probes that are closest to the principal component of
   * the data.
   *
   * @param psGroup - Collection of probe sets to get probes from.
   * @param iMart - Object containing raw data values for all chips.
   * @param doReport - Should the quantification report object be called?
   *
   * @return true if success, false otherwise.
   */
/*   virtual bool doAnalysis(ProbeSetGroup &psGroup, IntensityMart &iMart, bool doReport); */
  virtual bool doAnalysis(ProbeSetGroup &psGroup, 
                          IntensityMart &iMart, 
                          bool doReport,
                          bool alleleSummaryOnly = false);

  /**
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() {
    std::vector<SelfDoc::Opt> opts;
    SelfDoc::Opt doLog = {"log", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                          "Do log2 transformation on values before doing PCA feature selection."};
    opts.push_back(doLog);
    SelfDoc::Opt doDebug = {"debug", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                            "Print out debugging data files (can be very large)."};
    opts.push_back(doDebug);
    SelfDoc::Opt doCorr = {"corr", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                          "Use correlation rather than scatter matrix for PCA."};
    opts.push_back(doCorr);
    SelfDoc::Opt infoFilter = {"info-criterion", SelfDoc::Opt::String, "aic", "aic", "NA", "NA",
                              "Should we use and information criter ('bic','aic','none') to determine if a strong enough signal was discovered to warrant feature selection?"};
    opts.push_back(infoFilter);
    SelfDoc::Opt hardMin = {"hard-min", SelfDoc::Opt::Integer, "4", "4", "2", "NA",
                                 "Hard minimum on number of probes to use for summarization."};
    opts.push_back(hardMin);
    SelfDoc::Opt minPercent = {"min-percent", SelfDoc::Opt::Double, ".2", ".2", "0", "1",
                            "Minimum percentage of probes to use for summarization."};
    opts.push_back(minPercent);
    SelfDoc::Opt qNormOnly = {"qnorm-only", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                          "Use quantile normalized values rather than regular chipstream for this analysis."};
    opts.push_back(qNormOnly);
    return opts;
  }

  /**
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName("pca-select");
    doc.setDocDescription("Determines PCA for probes and picks probes that are near the principal component as the probes to use for downstream analysis.");
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
    bool doLog = true;
    bool doCorr = false;
    bool doDebug = false;
    bool qnormOnly = true;
    double minPercent = .1;
    int hardMin = 5;
    std::string infoFilter = "bic";
    fillInValue(doLog, "log", param, doc);
    fillInValue(doCorr, "corr", param, doc);
    fillInValue(doDebug, "debug", param, doc);
    fillInValue(infoFilter, "info-criterion", param, doc);
    fillInValue(minPercent, "min-percent", param, doc);
    fillInValue(hardMin, "hard-min", param, doc);
    fillInValue(qnormOnly, "qnorm-only", param, doc);
    AnalysisStreamExpPcaSel *stream = new AnalysisStreamExpPcaSel(doLog, doCorr, doDebug, infoFilter,
                                                                  hardMin, minPercent, qnormOnly);
    return stream;
  }

/*   static void subColAvg(Matrix &M); */
/*   static void MatrixScatter(Matrix &M, Matrix &C); */
/*   static void MatrixCov(Matrix &M, Matrix &C); */
/*   static void MatrixCor(Matrix &M, Matrix &C); */
/*   static void MaxEigen(Matrix &M, double &maxValue, ColumnVector &MaxVec); */

  static void fillInPmData(Matrix &mat, std::vector<probeid_t> &probeIds, ProbeSetGroup &psGroup,
                           std::vector<ChipStream *> &cStream, IntensityMart &iMart,
                           bool doLog);

  static void doFeatureSelection(ColumnVector &W, std::set<probeid_t> &goodIds,
                                 Matrix &PM, std::vector<probeid_t> &probeIds, bool doCorr,
                                 std::ofstream *out=NULL, const char *name=NULL);

  static void log2Matrix(Matrix &M);

  bool fillInSelectProbes(ProbeSetGroup &selectGroup, std::vector<ChipStream *> &cStream,
                          IntensityMart &iMart, ProbeSetGroup &psGroup, std::vector<double> &confVals);


  static void reportProbeLevelData(std::ofstream &out, ProbeSetGroup &psGroup, IntensityMart &iMart,
                                   std::vector<ChipStream *> &cStream);

  static void reportProbesUsed(std::ofstream &out,
                               ProbeSetGroup &used, ProbeSetGroup &orig,
                               std::vector<double> &confVals);

protected:
  bool m_DoLog;
  bool m_DoCorrelation;
  bool m_Debug;
  /// Which filter should we use to determine if it worth it to do feature selection?
  enum InfoFilter m_InfoFilter;
  std::ofstream m_Log;
  std::ofstream m_AllProbes;
  std::ofstream m_UsedProbes;
  std::ofstream m_ProbeWeights;
  SketchQuantNormTran *m_QNorm;
  /// What is the minimum number of probes acceptable after selection?
  int m_HardMinimum;
  /// What is the minimum percent of probes acceptable after selection?
  double m_HardProportion;
  /// Shoule we make our own quantile normalization object or use that from regular analysis stream.
  bool m_QuantNormOnly;
};

#endif /* _ANALYSISSTREAMEXPPCASEL_H_ */
