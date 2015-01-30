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
 * @file   QuantDabg.h
 * @author Chuck Sugnet
 * @date   Thu Dec 29 16:15:28 2005
 * 
 * @brief Class for doing Detected Above BackGround (DABG) analysis. Idea is to
 * use the intensity of probes in a probeset to calculate a p-value for the
 * entire probeset. The individual intensities are compared to a set of
 * background probes.
 */
#ifndef _QUANTDABG_H_
#define _QUANTDABG_H_

//
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityMart.h"
#include "chipstream/PmAdjuster.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantExprMethod.h"
//
#include "dabg/Dabg.h"
#include "util/Err.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/// our standardized name
#define QUANTDABGSTR "dabg"

/**
 * @brief Class for doing Detected Above BackGround (DABG) analysis. Idea is to
 * use the intensity of probes in a probeset to calculate a p-value for the
 * entire probeset. The individual intensities are compared to a set of
 * background probes.
 */
class QuantDabg : public QuantExprMethod {

public:

  /** 
   * What version of the algorithm are we implementing?
   * @return version string.
   */
  std::string getVersion() {
    return "2.0";
  }

  /**
   * Constructor.
   */
  QuantDabg(bool doFisher, bool doPercentile, double percentile, bool reportLog, int adjustment) :
    m_DoFisher(doFisher), m_DoPercentile(doPercentile), 
    m_PercentileEst(percentile), m_ReportLog(reportLog), m_Adjustment(adjustment) { 
    if(m_DoPercentile && m_DoFisher) {
      Err::errAbort("QuantDabg::QuantDabg() - Can't do chisq and percentile at same time.");
    }
    setupSelfDoc(*this);
    m_Type = getDocName();
    m_DabgInit = false;
	setOptValue("chisq", doFisher);
    setOptValue("percentile", ToStr(m_PercentileEst));
	setOptValue("usepercentile", doPercentile);
    setOptValue("neglog10", reportLog);
  }

  /** 
   * @brief Does this quantification method estimate feature effects?
   * @return false as we don't have feature effects.
   */
  inline bool haveFeatureEffects() { return false;}

  /** 
   * @brief Does this quantification method supply residuals?
   * @return false as dabg doesn't have a residual
   */
  inline bool haveResiduals() { return true;}

  inline double getResidual(unsigned int probeIx, unsigned int chipIx) { 
    return m_ProbePvals[chipIx][probeIx];
  }

  /** 
   * @brief Get the feature effect of a particular probe.
   * @param probeIx - Index of probe in probe sets.
   * @return feature (probe) effects.
   */
  inline double getFeatureEffect(unsigned int probeIx) {
    Err::errAbort("Don't do feature effects.");
    return -1;
  } 

  /** 
   * Which probes should be used as the background distribution.
   * @param bgProbes - vector of probes to be used, memory 
   * owned elsewhere.
   */
  void setBgProbes(const std::vector<Probe *> &bgProbes);

  /** 
   * Which probes should be used as the background distribution.
   * @param bgProbes - vector of probe ids to be used, memory 
   * owned elsewhere.
   */
  void setBgProbes(const std::vector<int> &bgProbes);

  /** 
   * How many probes are in the background list?
   * @return - number of probes being used for background.
   */
  uint32_t getBgProbeCount() {
    return (uint32_t) m_BgProbes.size();
  }

  void setParameters(PsBoard &board) {
    DataStore *info = board.getProbeInfo();
    std::vector<int> bgProbes;
    info->getGcControlProbes(bgProbes);
    setBgProbes(bgProbes);
    info->getProbeGc(m_ProbeGcVec);
  }

  /** 
   * @brief Get the p-value that a probeset came from the background
   * distribution for a particular chip.
   *
   * @param chipIx - Index of chip in experiment.
   * @return p-value for probeset in particular chip.
   */
  inline double getSignalEstimate(unsigned int chipIx) {
    assert(chipIx < m_ChipCount);
    return m_SetPvals[chipIx];
  }

  /** 
   * @brief Alias for getSignalEstimate().
   * @param chipIx - Index of chip in experiment.
   * @return p-value for probeset in particular chip.
   */
  inline double getTargetEffect(unsigned int chipIx) { 
    return getSignalEstimate(chipIx);
  }

  /** 
   * Nothing to clear yet. 
   */
  void clear() {}

  /** 
   * @brief Set the PM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   * @param data - Data to be set.
   */
  inline void setPMDataAt(unsigned int probeIx, unsigned int chipIx, double data) {
    assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
    m_PM[chipIx][probeIx] = data;
  }

/** 
   * @brief Get the PM data at probe set probe index and chip index.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   */
  inline double getPMDataAt(unsigned int probeIx, unsigned int chipIx) {
    assert(chipIx < m_ChipCount && probeIx < m_ProbeCount);
	return (double) m_PM[chipIx][probeIx];
  }

  /** 
   * @brief Set the MM data at probe set probe index and chip
   * index. Currently a no-operation as dabg doesn't use mm probes.
   * 
   * @param probeIx - Index of probe in probe set.
   * @param chipIx - Index of chip.
   * @param data - Data to be set.
   */
  inline void setMMDataAt(unsigned int probeIx, unsigned int chipIx, double data) {}

  /** 
   * Set the number of probes and number of chips to be used for estimation.
   * @param numProbes - Probe count.
   * @param numChips - Microarray count.
   */
  void setBounds(unsigned int numProbes, unsigned int numChips);

  /** 
   * What scale are the target effects on? Are they linear (like
   * plier) or log2 (like rma)
   * @return - What scale are the results provided in?
   */
  virtual Scale getScale() { return m_ReportLog ? NegLog10 : Pvalue; }

  /** 
   * What type of quantification method is this? Summarization like
   * plier detection like dabg.
   * @return - Type of quantification method.
   */
  virtual enum QuantType getQuantType() { return Detection; }

  /** 
   * @brief Do the heavy lifting of median polish. If PM-MM is zero
   * threshold at a small positive number to avoid taking logs of
   * numbers that aren't positive.
   */
  void computeEstimate();

  /** 
   * @brief Set up the quantification method given all the data about the probe
   * set, chip  and data.
   * 
   * @param psGroup - Probes to be used for final estimate.
   * @param iMart - Raw data from chips.
   * @param iTrans - Transformations to be applied to data before use.
   * @param pmAdjust - How to estimate background, or MM probe.
   * @return True if setup sucessful, false otherwise.
   */
  bool setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
             std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust);

  /** 
   * @brief Return number of features (probes).
   * @return feature count.
   */
  inline unsigned int getNumFeatures()  { return m_ProbeCount; }

  /** 
   * @brief Return number of targets (experiments or chips).
   * @return target count
   */
  inline unsigned int getNumTargets() { return m_ChipCount; }

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> opts;
	SelfDoc::Opt chisq = {"chisq", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                          "Use Fisher's chi-squared method for combining individual probe p-values."};
    opts.push_back(chisq);
    SelfDoc::Opt usepercentile = {"usepercentile", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                          "Use a individual probe p-value at a particular percentile as the estimator."};
    opts.push_back(usepercentile);
    SelfDoc::Opt percentile = {"percentile", SelfDoc::Opt::Double, ".25", ".25", "0", "1",
                               "Percentile to use for the probe p-value as an estimator for the entire probeset."};
    opts.push_back(percentile);
    SelfDoc::Opt neglog10 = {"neglog10", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                             "Report -1 * log_10(p-value) rather than raw p-value."};
    opts.push_back(neglog10);
	SelfDoc::Opt adjust = {"adjust-dof", SelfDoc::Opt::Integer, "0", "0", "NA", "NA",
                       "Adjustment of the degrees of freedom for the chi-squared test."};
    opts.push_back(adjust);
	SelfDoc::Opt probeMd5 = {"subsetmd5", SelfDoc::Opt::String, "", "", "NA", "NA",
                             "Md5sum of the probe ids being used as background."};
    opts.push_back(probeMd5);
    return opts;
  }

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName(QUANTDABGSTR);
    doc.setDocDescription("Calculates the p-value that the intensities in a probeset could have been observed by chance in a background distribution. Used as a substitute for standard absent/present calls when mismatch probes are not available.");
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
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

  /**
   * @brief What should the residual file be called?
   */
  virtual const std::string &getResidualSuffix() {
    static std::string suffix = ".probe-pvalues.txt";
    return suffix;
  }

  /**
   * Specify the gc count for probes as indexed by position
   */
  void setGcProbes(const std::vector<char> &gcVec) {
    m_ProbeGcVec = gcVec;
  }

  /**
   * set the bgpFile name
   * @param bgpFile
   */
  void setBgpFileName(std::string bgpFile) {
    m_BgpFile = bgpFile;
  }

private:
  /** 
   * Setup the dabg objects with intensity data from background probes. This should
   * be called only once.
   * 
   * @param iMart - Contains data for background probes for each chip.
   * @param iTrans - Chipstream transformations for raw data in iMart.
   * @param pmAdjust - Adjuster for perfect match probes.
   */
  void initializeDabg(const IntensityMart &iMart, 
                      std::vector<ChipStream *> &iTrans, 
                      PmAdjuster &pmAdjust);

  /** 
   * Return the -1 * log_10(val). Val is capped at FLT_MIN to avoid taking lots of values <= 0. 
   * @param val - Regular p-value.
   * @return -1 * log_10(val)
   */
  double negLog10Prob(double val);

  /// Number of chips in experiment
  unsigned int m_ChipCount;
  /// Number of probes (features) in probe sets
  unsigned int m_ProbeCount;
  /// Data matrix for PM probes.
  std::vector<std::vector<double> > m_PM;
  /// G/C count for each probe in probeset being called
  std::vector<int> m_GcProbeCounts;
  /// P-values for probeset in different chips.
  std::vector<double> m_SetPvals;
  /// P-values for individual probes in different chips.
  std::vector<std::vector<double> > m_ProbePvals;
  /// Vector of probes to use as background distribution. Memory owned elsewhere
  std::vector<int> m_BgProbes;
  /// Have we initialized the dabg object with data yet?
  bool m_DabgInit;
  /// Objects implementing algorithm, one per chip.
  std::vector<affx::Dabg> m_Dabgs;
  /// Use Fisher's chi-squared method for combining individual probe p-values.
  bool m_DoFisher;
  /// Use a partcular percentile of probe p-values as an estimator.
  bool m_DoPercentile;
  /// Percentile of probe p-value to use as estimator.
  double m_PercentileEst;
  /// Should we report -1 * log_10(p-value) rather than raw p-value.
  double m_ReportLog;
  /// Adjust the degrees of freedom by one.
  int m_Adjustment;
  /// Md5sum of probe ids used for background
  std::string m_ProbeMd5Sum;
  /// GC count of each probe as indexed by position.
  std::vector<char> m_ProbeGcVec;
  /// Name of bgpFile
  std::string m_BgpFile;
};

#endif /* _QUANTDABG_H_ */
