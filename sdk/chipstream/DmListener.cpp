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
 * @file   DmListener.cpp
 * @author Chuck Sugnet
 * @date   Sun Mar 12 11:55:49 2006
 *
 * @brief Class for using DM calls as cell files become available.
 */

//
#include "chipstream/DmListener.h"
//
#include <cmath>
//
using namespace std;
using namespace affx;

DmListener::DmCall::DmCall()
{
  good = false;
  call = affx::NN; // unknown
  conf = 0;
}

/**
 * Constructor takes a vector of ProbeSets for which calls will be
 * done as cel files become available.
 *
 * @param psVec - Vector of probesets to be calculated.  @param
 * hetMult - Factor to add to log likelihood to balance het/hom calls
 * 0 for no effect.
 *
 * @return
 */


DmListener::DmListener(std::vector<ProbeListPacked> &probesets, float hetMult,
                       int maxChips, std::vector<float> &thresholds,
                       std::map<string, bool> &chrXSnps, std::string callRateName)
{
  std::vector<ProbeSet *>::iterator iter;
  Err::check(!thresholds.empty(), "DmListener::DmListener() - Can't supply no thresholds.");
  m_HetMultiplier = hetMult;
  m_ChipCount = 0;
  m_MaxChips = maxChips;
  m_Thresholds = thresholds;
  m_TotalCalls.resize(m_MaxChips);
  m_PassCalls.resize(m_Thresholds.size());
  m_HetChrXCalls.resize(m_Thresholds.size());
  m_FirstPassGender.resize(maxChips, affx::UnknownGender);
  m_FirstPassHetRate.resize(maxChips, 0.0);
  m_ChrXCalls.resize(maxChips);
  m_ChrXSnps = chrXSnps;
  m_CallIndex = 0;
  m_ChrXCallIndex = 0;
  for (uint32_t i = 0; i < m_PassCalls.size(); i++) {
    m_PassCalls[i].resize(m_MaxChips);
    m_HetChrXCalls[i].resize(m_MaxChips);
  }
  for (uint32_t i = 0; i < probesets.size(); i++) {
    m_PlVec.push_back(probesets[i]);
  }

  m_KeepGenotypes = true;


  m_GenderName = "dm-chrX-het-rate";
  m_GenderDescription = "DM chrX het rate";
  setPrefix("dm-listener");

  if (callRateName != "")
    setCallRateName(callRateName);

  if (m_ChrXSnps.size() > 0) {
    declareMetric(m_GenderName + "_gender", ChipSummary::Metric::String);
    declareMetric(m_GenderName + "_gender_chrX_het_rate", ChipSummary::Metric::Double);
  }
  declareMetric(m_CallRateName, ChipSummary::Metric::Double);
}

DmListener::~DmListener()
{
}


/**
 * Clear genotype results. Gender results are retained.
 * Call rates (total and chrX) are retained and accumulated
 * over calls to this. The only thing this does is it
 * drops the genotype calls.
 * This allows us to save memory by not keeping around
 * genotypes that are no longer needed.
 */

void DmListener::clearForIteration()
{
  m_KnownCalls.clear();
  m_ChipCount = 0;
}

/**
 * Do we need to keep genotypes around for providing
 * seed calls? Default yes. If set to no, then
 * getGenoCalls will fail if called.
 * @param s - true or false
 */
void DmListener::setKeepGenotypes(bool s)
{
  m_KeepGenotypes = s;
}


void DmListener::setLayout(ChipLayout &layout)
{
  m_PlVec.clear();
  for (uint32_t i = 0; i < layout.getProbeSetCount(); i++) {
    m_PlVec.push_back(layout.getProbeListAtIndex(i));
  }
}

/**
 * @brief Check if a call exists for a given probeset name
 * @param name - name of probeset to get genotype calls for.
 * @return boolean result
 */
bool DmListener::checkGenoCallsName(const std::string &name)
{
  return m_KnownCalls.end() != m_KnownCalls.find(name);
}

/** Accessor function. */
int DmListener::getChipCount()
{
  return m_ChipCount;
}

/**
 * Get the call rate at a particular threshold. Note that we can keep
 * track of the rate at multiple thresholds, but only keep the actual
 * calls for a particular threshold due to space (RAM) concerns.
 *
 * @param callRates - Vector to be filled in with call rates, one for each chip.
 * @param thresholdIndex - Which threshold to use (as defined in constructor).
 */
void DmListener::getCallRates(std::vector<float> &callRates, int thresholdIndex)
{
  callRates.clear();
  callRates.resize(m_ChipCount);
  assert(thresholdIndex < m_PassCalls.size());
  assert(m_PassCalls[thresholdIndex].size() == m_ChipCount);
  for (uint32_t i = 0; i < m_ChipCount; i++) {
    assert(m_TotalCalls[i] > 0);
    callRates[i] = (float)m_PassCalls[thresholdIndex][i] / m_TotalCalls[i];
  }
}

/**
 * Get the call rate for a specific chip
 * @param chip - the chip index to pull
 * @param thresholdIndx - which threshold should be applied
 * @return The call rate for that chip
 */
float DmListener::getCallRate(int chip, int thresholdIndex)
{
  if (m_TotalCalls[chip] == 0) {
    return -1;
  } else {
    //Verbose::out(1,"Chip " + ToStr(chip) + " has " + ToStr(m_PassCalls[thresholdIndex][chip]) + " passing calls out of " + ToStr(m_TotalCalls[chip]));
    return (float)m_PassCalls[thresholdIndex][chip] / m_TotalCalls[chip];
  }
}


/**
 * Get the call rates for chrX snps at a particular threshold. Note that we
 * can keep track of the rate at multiple thresholds, but only keep the actual
 * calls for a particular threshold due to space (RAM) concerns.
 *
 * @param callRates - Vector to be filled in with call rates for chrx snps,
 *                    one for each chip.
 * @param thresholdIndex - Which threshold to use (as defined in constructor.
 */
void DmListener::getHetChrXRates(std::vector<float> &hetRates, int thresholdIndex)
{
  hetRates.clear();
  hetRates.resize(m_ChipCount);
  assert(thresholdIndex < m_HetChrXCalls.size());
  assert(m_HetChrXCalls[thresholdIndex].size() == m_ChipCount);
  for (uint32_t i = 0; i < m_ChipCount; i++) {
    if (m_ChrXCalls[i] == 0) {
      hetRates[i] = -1;
    } else {
      hetRates[i] = (float)m_HetChrXCalls[thresholdIndex][i] / m_ChrXCalls[i];
    }
  }
}

/**
 * Get the chrX het rate for a specific chip
 * @param chip - the chip index to pull
 * @param thresholdIndex - which threshold should be applied
 * @return the chrX het rate for that chip
 */
float DmListener::getHetChrXRate(int chip, int thresholdIndex)
{
  if (m_ChrXCalls[chip] == 0) {
    return -1;
  } else {
    return (float)m_HetChrXCalls[thresholdIndex][chip] / m_ChrXCalls[chip];
  }
}


/**
 * Set the index (into the thresholds vector in constructor) that we will
 * use by default when queried for call rates.
 *
 * @param index - index into the call rates at different thresholds
 * we keep track of as defined by the thresholds vector in
 * constructor.
 */
void DmListener::setDefaultCallIndex(int index)
{
  if ((size_t)index >= m_PassCalls.size()) {
    Err::errAbort("DmListener::setDefaultCallIndex() - Can't set default to: " +
                  ToStr(index) + " when there are at most: " + ToStr(m_PassCalls.size()));
  }
  m_CallIndex = index;
}

/**
 * Set the index (into the thresholds vector in constructor) that we will
 * use by default when queried for call rates.
 *
 * @param index - index into the call rates at different thresholds
 * we keep track of as defined by the thresholds vector in
 * constructor.
 */
void DmListener::setDefaultChrXCallIndex(int index)
{
  if ((size_t)index >= m_HetChrXCalls.size()) {
    Err::errAbort("DmListener::setDefaultChrXCallIndex() - Can't set default to: " +
                  ToStr(index) + " when there are at most: " + ToStr(m_HetChrXCalls.size()));
  }
  m_ChrXCallIndex = index;
}

/**
 * Get the call rates for snps for each chip.
 * @param callRates - Vector to be filled in with call rates for chrx snps,
 * one for each chip.
 */
void DmListener::getCallRates(std::vector<float> &callRates)
{
  getCallRates(callRates, m_CallIndex);
}

/**
 * Get the call rate for a specific chip
 * @param chip - chip index to get call rate for
 * @return The call rate
 */
float DmListener::getCallRate(int chip)
{
  return getCallRate(chip, m_CallIndex);
}

/**
 * Get the call rates for chrX snps for each chip.
 * @param callRates - Vector to be filled in with call rates for chrx snps,
 * one for each chip.
 */
void DmListener::getHetChrXRates(std::vector<float> &hetRates)
{
  getHetChrXRates(hetRates, m_ChrXCallIndex);
}

/**
 * Get the chrX het rate for a specific chip
 * @param chip - the chip index to get call rate for
 * @return the chrC het rate
 */
float DmListener::getHetChrXRate(int chip)
{
  return getHetChrXRate(chip, m_ChrXCallIndex);
}

/**
 * Will this probeset work for DM? That is does it have MM probes and
 * matched quartets of probes on A and B allele?
 * @param ps - Probeset to check.
 * @return - true if ok for DM, false otherwise.
 */
bool DmListener::okForDm(const ProbeSet *ps)
{
  bool ok = true;
  if (ps->psType != ProbeSet::GenoType) {
    ok = false;
  } else if (ps->numGroups != 2 && ps->numGroups != 4) {
    ok = false;
  } else if (ps->atomsPerGroup.empty()) {
    Err::errAbort("DmListener::okForDm() - Can't have multiple groups, but atomsPerGroup be emtpy in probeset: '" + ToStr(ps->name) + "'");
    ok = false;
  } else if (!ps->hasMM()) {
    ok = false;
  } else {
    for (int psIx = 0; psIx < ps->numGroups - 1; psIx += 2) {
      if (ps->atomsPerGroup[psIx] != ps->atomsPerGroup[psIx+1])
        ok = false;
    }
  }
  return ok;
}

void DmListener::fillInQuartet(vector<CQuartet> &qVec, const ProbeSet *ps, affymetrix_fusion_io::FusionCELData *cel)
{
  assert(ps);
  if (!okForDm(ps))
    Err::errAbort("Probeset " + ToStr(ps->name) + " not suitable for DM.");

  /* Genotype probe sets can have either 2 or 4 groups, for each
     type the even (0,2) are A alleles and odd (1,3) are B allele.s */
  int atomIndex = 0;
  for (int psIx = 0; psIx < ps->numGroups - 1; psIx += 2) {
    Err::check(ps->atomsPerGroup[psIx] == ps->atomsPerGroup[psIx+1],
               "DmListener::makeCall() - Though groups were symmetrical.");
    int atomsPer = ps->atomsPerGroup[psIx];
    /* Each block will have atomsPer atoms in there. */
    for (int atomIx = 0; atomIx < atomsPer; atomIx++) {
      int aIndex = atomIx + (atomIndex);
      int bIndex = aIndex + atomsPer;
      Atom *A = ps->atoms[aIndex];
      Atom *B = ps->atoms[bIndex];

      // Fill in PM/MM pair
      if (ps->hasMM()) {
        int aPM = -1, aMM = -1;
        assert(A->probes.size() % 2 == 0);
        assert(B->probes.size() % 2 == 0);
        for (int probeIx = 0; probeIx + 1 < A->probes.size(); probeIx += 2) {
          CQuartet q;
          int first = probeIx;
          int next = probeIx + 1;
          /* this is a little akward, need to pull out the a & b allele info,
             but I don't think there is a hard and fast rule about the order
             of the pm and mm probes, so check both for both a & b. */

          /* A allele */
          if (Probe::isPm(*A->probes[first]) && Probe::isMm(*A->probes[next])) {
            aPM = first;
            aMM = next;
          } else if (Probe::isPm(*A->probes[next]) && Probe::isMm(*A->probes[first])) {
            aPM = next;
            aMM = first;
          } else {
            Err::errAbort("Need a pm and mm probe in probeset: " + ToStr(ps->name) + " atom: " + ToStr(atomIx));
          }
          q.SetIntensity(CQuartet::A_Allele, CQuartet::PM, cel->GetIntensity(A->probes[aPM]->id));
          q.SetIntensity(CQuartet::A_Allele, CQuartet::MM, cel->GetIntensity(A->probes[aMM]->id));
          double pm_stdev = cel->GetStdv(A->probes[aPM]->id);
          double mm_stdev = cel->GetStdv(A->probes[aMM]->id);
          
          q.SetVariance(CQuartet::A_Allele, CQuartet::PM, pm_stdev * pm_stdev);
          q.SetVariance(CQuartet::A_Allele, CQuartet::MM, mm_stdev * mm_stdev);
          q.SetPixels(CQuartet::A_Allele, CQuartet::PM, cel->GetPixels(A->probes[aPM]->id));
          q.SetPixels(CQuartet::A_Allele, CQuartet::MM, cel->GetPixels(A->probes[aMM]->id));

          /* B Allele */
          if (Probe::isPm(*B->probes[first]) && Probe::isMm(*B->probes[next])) {
            aPM = first;
            aMM = next;
          } else if (Probe::isPm(*B->probes[next]) && Probe::isMm(*B->probes[first])) {
            aPM = next;
            aMM = first;
          } else {
            Err::errAbort("Need a pm and mm probe in probeset: " + ToStr(ps->name) + " atom: " + ToStr(atomIx));
          }
          q.SetIntensity(CQuartet::B_Allele, CQuartet::PM, cel->GetIntensity(B->probes[aPM]->id));
          q.SetIntensity(CQuartet::B_Allele, CQuartet::MM, cel->GetIntensity(B->probes[aMM]->id));
          pm_stdev = cel->GetStdv(B->probes[aPM]->id);
          mm_stdev = cel->GetStdv(B->probes[aMM]->id);
          q.SetVariance(CQuartet::B_Allele, CQuartet::PM, pm_stdev * pm_stdev);
          q.SetVariance(CQuartet::B_Allele, CQuartet::MM, mm_stdev * mm_stdev);
          q.SetPixels(CQuartet::B_Allele, CQuartet::PM, cel->GetPixels(B->probes[aPM]->id));
          q.SetPixels(CQuartet::B_Allele, CQuartet::MM, cel->GetPixels(B->probes[aMM]->id));
          qVec.push_back(q);
        }
      } else {
        Err::errAbort("Need a pm and mm probe in probeset: " + ToStr(ps->name) + " atom: " + ToStr(atomIx));
      }
    } // for each atom
    atomIndex += (2 * atomsPer);
  } // for each group
}

/**
 * Make DM genotyping call on a particular probeset given the data in the
 * cel file.
 *
 * @param ps - probeset to do genotyping call.
 * @param cel - cel file to use for data.
 *
 * @return Genotype call and confidence.
 */
DmListener::DmCall DmListener::makeCall(float hetMult, const ProbeSet *ps, affymetrix_fusion_io::FusionCELData *cel)
{
  assert(ps);
  // set up the result.
  DmCall result;
  result.call = NN; // -1; unknown
  result.conf = 0;
  // if this probeset isn't compatible with DM then return NN call.
  if (!okForDm(ps)) {
    return result;
  }

  vector<CQuartet> qVec;
  fillInQuartet(qVec, ps, cel);

  std::pair<float, int> call;
  // Here is the actual work.
  if (CallDM(qVec, call, hetMult)) {
    result.good = true;
    result.conf = call.first;
    result.call = GType_from_int(call.second);
  } else {
    Verbose::out(1, "Dm call for probeset: '" + ToStr(ps->name) + "' failed.");
  }

  // sanity check...
  // assert(result.call >= -1 && result.call <= 2);
  assert(GType_valid(result.call));
  if (result.conf < 0 || result.conf > 1) {
    Err::errAbort("DmListener::makeCall() - Expecting p-value between 0 and 1, got: " + ToStr(result.conf));
  }
  return result;
}

void DmListener::assignCall(std::string &name, DmListener::DmCall &call, std::vector<affx::GType> &genoCall)
{
  // keep statistics for all call rates.
  // some SNPs we might see on multiple iterations (ie chrX snps)
  map<string, vector<bool> >::iterator mapIter;

  // Initialize SNP to not seen, if not in map
  mapIter = m_SeenSNPs.find(name);
  if (mapIter == m_SeenSNPs.end()) {
    m_SeenSNPs[name].resize(m_MaxChips, false);
  }

  // If not seen, the compute stats
  if (m_SeenSNPs[name][m_ChipCount] == false) {
    for (unsigned int threshIx = 0; threshIx < m_Thresholds.size(); threshIx++) {
      if (call.good && call.conf < m_Thresholds[threshIx]) {
        m_PassCalls[threshIx][m_ChipCount]++;
      }
    }
    if (call.good) {
      m_TotalCalls[m_ChipCount]++;
    }
    //Verbose::out(1,"Adding Call: " + ToStr(name) + ", " + ToStr(call.conf) + ", " + ToStr(call.good));
    if (m_ChrXSnps.find(name) != m_ChrXSnps.end()) {
      for (unsigned int threshIx = 0; threshIx < m_Thresholds.size(); threshIx++) {
        if (call.good && call.conf < m_Thresholds[threshIx] && call.call == AB) {
          m_HetChrXCalls[threshIx][m_ChipCount]++;
        }
      }
      m_ChrXCalls[m_ChipCount]++;
    }
    m_SeenSNPs[name][m_ChipCount] = true;
  }

  // use index 0 to make actual call, when DM fails we set to NN.
  if (!call.good || call.conf >= m_Thresholds[0]) {
    call.call = NN;
  }
  genoCall[m_ChipCount] = call.call;
}

/**
 * Process another cel files worth of data.
 *
 * @param cel - Filehandle to open cel file.
 */
void DmListener::newChip(affymetrix_fusion_io::FusionCELData *cel)
{
  for (unsigned int psIx = 0; psIx < m_PlVec.size(); psIx++) {
    ProbeSet *ps = ProbeListFactory::asProbeSet(m_PlVec[psIx]);
    std::map<string, std::vector<GType> >::iterator mapIter;
    DmCall call = makeCall(m_HetMultiplier, ps, cel);

    std::string name = ps->name;
    mapIter = m_KnownCalls.find(name);
    if (mapIter == m_KnownCalls.end()) {
      // First time seen.
      std::vector<GType> genoCall(m_MaxChips);
      assignCall(name, call, genoCall);
      if (m_KeepGenotypes)
        m_KnownCalls[name] = genoCall;
    } else {
      // More have been seen
      assignCall(name, call, mapIter->second);
    }
    delete ps;
  }
  if (m_ChipCount == 0) {
    // trim excess capacity.
    std::map<string, std::vector<affx::GType> > tmp(m_KnownCalls);
    m_KnownCalls.swap(tmp);
  }

  // track the gender info over multiple iterations
  if (m_ChipCount >= m_Seen.size())
    m_Seen.resize(m_ChipCount + 1, false);
  if (!m_Seen[m_ChipCount]) {
    std::vector<ChipSummary::Metric> metrics;
    m_SummaryStats.push_back(metrics);
    if (m_ChrXSnps.size() > 0) {
      m_SummaryStats[m_ChipCount].push_back(ChipSummary::Metric(m_GenderName + "_gender",
                                            affx::getGenderString(getGender(m_ChipCount))));
      m_SummaryStats[m_ChipCount].push_back(ChipSummary::Metric(m_GenderName + "_gender_chrX_het_rate",
                                            getHetChrXRate(m_ChipCount)));
      m_FirstPassGender[m_ChipCount] = getGender(m_ChipCount);
      m_FirstPassHetRate[m_ChipCount] = getHetChrXRate(m_ChipCount);
    }
    m_SummaryStats[m_ChipCount].push_back(ChipSummary::Metric(m_CallRateName,
                                          getCallRate(m_ChipCount)));
    m_Seen[m_ChipCount] = true;
  } else {
    // Update call rate
    bool found = false;
    for (int i = 0; i < m_SummaryStats[m_ChipCount].size(); i++) {
      if (m_SummaryStats[m_ChipCount][i].m_Name == m_CallRateName) {
        float callRate = getCallRate(m_ChipCount);
        m_SummaryStats[m_ChipCount][i].m_Double = callRate;
        //Verbose::out(1,"Call rate update on chip " + ToStr(m_ChipCount) + " to " + ToStr(callRate));
        found = true;
      }
    }
    if (!found) {
      Err::errAbort("Unable to update dm-listener call rate for chip " + ToStr(m_ChipCount));
    }

    // We expect no change in gender and chrX het rate -- ie all
    // the chrX snps must be seen on the first pass
    if (m_ChrXSnps.size() > 0) {
      if (fabs(m_FirstPassHetRate[m_ChipCount] - getHetChrXRate(m_ChipCount)) > 0.00001)
        Err::errAbort(ToStr("DmListener chrX het rate changed between iterations!") +
                      " It was " + ToStr(m_FirstPassHetRate[m_ChipCount]) +
                      " but is now " + ToStr(getHetChrXRate(m_ChipCount)));
      if (m_FirstPassGender[m_ChipCount] != getGender(m_ChipCount))
        Err::errAbort(ToStr("DmListener computed gender changed between iterations!") +
                      " It was " + affx::getGenderString(m_FirstPassGender[m_ChipCount]) +
                      " but is now " + affx::getGenderString(getGender(m_ChipCount)));
    }
  }

  setValid(true);
  m_ChipCount++;
}

/**
 * Get the genotype calls for a particular probeset.
 * @param name - name of probeset to get genotype calls for.
 * @return - A vector of the calls from dm algorithm.
 */
std::vector<affx::GType> DmListener::getGenoCalls(const std::string &name)
{
  if (!m_KeepGenotypes)
    Err::errAbort("Genotypes requested from DmListener bug flag to keep them was set to false!");
  std::map<std::string, std::vector<affx::GType> >::iterator mapIter;
  mapIter = m_KnownCalls.find(name);
  if (mapIter == m_KnownCalls.end()) {
    Err::errAbort("DmListener::getGenoCall() - Can't find genotypes for name: " + name);
  }
  return mapIter->second;
}

/**
 * Set the full name to use for the call rate metric reported
 * via the chip summary interface.
 * @param label - the full name to use for the call rate metric
 */
void DmListener::setCallRateName(std::string label) {
  m_CallRateName = label;
}

/**
 * Set the prefix to use on chip summary metric labels. Defaults
 * to "dm-listener".
 * @param prefix - the prefix to use
 */
void DmListener::setPrefix(std::string prefix) {
  m_Prefix = prefix;
  m_CallRateName = m_Prefix + "_call_rate";
}

/* Quick and dirty. */
float DmListener::square(float x) {
  return x * x;
}

