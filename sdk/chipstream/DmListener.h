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
 * @file   DmListener.h
 * @author Chuck Sugnet
 * @date   Sun Mar 12 12:20:34 2006
 *
 * @brief Class for using DM calls as cell files become available.
 */

#ifndef _DMLISTENER_H_
#define _DMLISTENER_H_

//
#include "chipstream/CelListener.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/GenoSeed.h"
#include "chipstream/ProbeListFactory.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "dm/DM.h"
#include "util/Err.h"
//
#include <map>
#include <vector>
//

/**
 * @brief Class for using DM calls as cell files become available.
 */
class DmListener : public CelListener, public GenoSeed, public ChipSummary
{

public:

  /**
   * Utility class for result of DM algorithm.
   */
  class DmCall
  {
  public:

    DmCall();
    bool good;  ///< Was the resulting call successful?
    affx::GType call; ///< The genotype call: NN(-1),AA(0),AB(1),BB(2)
    float conf;       ///< p-value for DM call.
  };

  /**
   * Set a new layout -- allows us to save some memory and only deal
   * with SNPs in the current iteration. Call this method before
   * starting new interations to set the new layout. You may
   * also want to call clearForIteration.
   * @param layout - the new layout of probes to run DM on
   */
  void setLayout(ChipLayout &layout);

  /**
   * Clear genotype results. Gender results are retained.
   * Call rates (total and chrX) are retained and accumulated
   * over calls to this. The only thing this does is it
   * drops the genotype calls.
   * This allows us to save memory by not keeping around
   * genotypes that are no longer needed.
   */
  void clearForIteration();

  /**
   * Do we need to keep genotypes around for providing
   * seed calls? Default yes. If set to no, then
   * getGenoCalls will fail if called.
   * @param s - true or false
   */
  void setKeepGenotypes(bool s);

  /**
   * Constructor takes a vector of ProbeSets for which calls will be
   * done as cel files become available.
   *
   * @param probesets - Vector of probesets to be calculated.
   * @param hetMult - Factor to add to log likelihood to balance het/hom calls
   *                  0 for no effect.
   * @param maxChips - the maximum number of chips we expect to see
   * @param thresholds - vector of floats to threshold the calls at
   * @param chrXSnps - map of chrX snp probeset names
   */
  DmListener(std::vector<ProbeListPacked> &probesets, float hetMult, int maxChips,
             std::vector<float> &thresholds, std::map<std::string, bool> &chrXSnps,
             std::string callRateName = "");

  /**
   * Virtual destructor.
   */
  virtual ~DmListener();

  /**
   * Will this probeset work for DM? That is does it have MM probes and
   * matched quartets of probes on A and B allele?
   * @param ps - Probeset to check.
   * @return - true if ok for DM, false otherwise.
   */
  static bool okForDm(const ProbeSet *ps);

  /**
   * Make DM genotyping call on a particular probeset given the data in the
   * cel file.
   * @param hetMult - het multiplier to use
   * @param ps - probeset to do genotyping call.
   * @param cel - cel file to use for data.
   * @return Genotype call and confidence.
   */
  static DmCall makeCall(float hetMult, const ProbeSet *ps,
                         affymetrix_fusion_io::FusionCELData *cel);

  /**
   * Fill in a vector of CQuartet containing the cel file info for
   * the given probeset.
   * @param qVec - vector of CQuartets to fill in
   * @param ps - probeset which we are asking for
   * @param cel - the cel file we are processing
   */
  static void fillInQuartet(std::vector<CQuartet> &qVec, const ProbeSet *ps,
                            affymetrix_fusion_io::FusionCELData *cel);

  /**
   * Threshold the call based on the default threshold value. Also
   * keep track of total call rates and ChrX het rates.
   * @param name - probeset name
   * @param call - the dm call
   * @param genoCall - the threshold processed call to fill in
   */
  void assignCall(std::string &name, DmListener::DmCall &call,
                  std::vector<affx::GType> &genoCall);

  /**
   * Process another cel files worth of data.
   * @param cel - Filehandle to open cel file.
   */
  void newChip(affymetrix_fusion_io::FusionCELData *cel);

  /**
   * Get the genotype calls for a particular probeset.
   * @param name - name of probeset to get genotype calls for.
   * @return - A vector of the calls from dm algorithm.
   */
  std::vector<affx::GType> getGenoCalls(const std::string &name);

  /**
   * @brief Check if a call exists for a given probeset name
   * @param name - name of probeset to get genotype calls for.
   * @return boolean result
   */
  bool checkGenoCallsName(const std::string &name);

  /** Accessor function. */
  int getChipCount();

  /**
   * Get the call rate at a particular threshold. Note that we can keep
   * track of the rate at multiple thresholds, but only keep the actual
   * calls for a particular threshold due to space (RAM) concerns.
   *
   * @param callRates - Vector to be filled in with call rates, one for each chip.
   * @param thresholdIndex - Which threshold to use (as defined in constructor).
   */
  void getCallRates(std::vector<float> &callRates, int thresholdIndex);

  /**
   * Get the call rate for a specific chip
   * @param chip - the chip index to pull
   * @param thresholdIndx - which threshold should be applied
   * @return The call rate for that chip
   */
  float getCallRate(int chip, int thresholdIndex);

  /**
   * Get the call rates for chrX snps at a particular threshold. Note that we
   * can keep track of the rate at multiple thresholds, but only keep the actual
   * calls for a particular threshold due to space (RAM) concerns.
   *
   * @param callRates - Vector to be filled in with call rates for chrx snps,
   *                    one for each chip.
   * @param thresholdIndex - Which threshold to use (as defined in constructor.
   */
  void getHetChrXRates(std::vector<float> &hetRates, int thresholdIndex);

  /**
   * Get the chrX het rate for a specific chip
   * @param chip - the chip index to pull
   * @param thresholdIndex - which threshold should be applied
   * @return the chrX het rate for that chip
   */
  float getHetChrXRate(int chip, int thresholdIndex);


  /**
   * Set the index (into the thresholds vector in constructor) that we will
   * use by default when queried for call rates.
   *
   * @param index - index into the call rates at different thresholds
   * we keep track of as defined by the thresholds vector in
   * constructor.
   */
  void setDefaultCallIndex(int index);

  /**
   * Set the index (into the thresholds vector in constructor) that we will
   * use by default when queried for call rates.
   *
   * @param index - index into the call rates at different thresholds
   * we keep track of as defined by the thresholds vector in
   * constructor.
   */
  void setDefaultChrXCallIndex(int index);

  /**
   * Get the call rates for snps for each chip.
   * @param callRates - Vector to be filled in with call rates for chrx snps,
   * one for each chip.
   */
  void getCallRates(std::vector<float> &callRates);

  /**
   * Get the call rate for a specific chip
   * @param chip - chip index to get call rate for
   * @return The call rate
   */
  float getCallRate(int chip);

  /**
   * Get the call rates for chrX snps for each chip.
   * @param callRates - Vector to be filled in with call rates for chrx snps,
   * one for each chip.
   */
  void getHetChrXRates(std::vector<float> &hetRates);

  /**
   * Get the chrX het rate for a specific chip
   * @param chip - the chip index to get call rate for
   * @return the chrC het rate
   */
  float getHetChrXRate(int chip);

private:

  /**
   * Set the full name to use for the call rate metric reported
   * via the chip summary interface.
   * @param label - the full name to use for the call rate metric
   */
  void setCallRateName(std::string label);

  /**
   * Set the prefix to use on chip summary metric labels. Defaults
   * to "dm-listener".
   * @param prefix - the prefix to use
   */
  void setPrefix(std::string prefix);

  /* Quick and dirty. */
  float square(float x);

  /// Chip count, how many have we seen.
  int m_ChipCount;
  /// How many chips are we expecting to see?
  int m_MaxChips;
  /// Vector of probesets to be called;
  std::vector<ProbeListPacked> m_PlVec;
  /// thresholds, first is used as cutoff others just have call rates calculated
  std::vector<float> m_Thresholds;
  /// calls at diferent thresholds.
  std::vector<std::vector<int> > m_PassCalls;
  /// tally of all the possible calls.
  std::vector<int> m_TotalCalls;
  /// chrX het rate at different thresholds
  std::vector<std::vector<int> > m_HetChrXCalls;
  /// tally of all the possible chrX calls
  std::vector<int> m_ChrXCalls;
  /// map of chrX SNPs
  std::map<std::string, bool> m_ChrXSnps;
  ///  Genotype calls indexed by chip.
  std::map<std::string, std::vector<affx::GType> > m_KnownCalls;
  /// Multiplier to help balance hom/het calls
  float m_HetMultiplier;
  /// threshold index to use for call rates by default.
  int m_CallIndex;
  /// threshold index to use for chrX call rates by default.
  int m_ChrXCallIndex;
  /// have we already seen this chip
  std::vector<bool> m_Seen;
  /// Do we keep genotypes around
  bool m_KeepGenotypes;
  /// prefix for dm call rate metric reported
  std::string m_Prefix;
  /// The name to use for the call rate Metric
  std::string m_CallRateName;
  /// Track what snps have been seen over all iterations
  std::map< std::string, std::vector< bool > > m_SeenSNPs;
  /// Gender computed from first iteration
  std::vector<affx::Gender> m_FirstPassGender;
  /// chrX het rate from first iteration
  std::vector<float> m_FirstPassHetRate;
};

#endif /* _DMLISTENER_H_ */
