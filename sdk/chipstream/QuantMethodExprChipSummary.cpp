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
 * @file   QuantMethodExprChipSummary.cpp
 * @author Alan Williams
 * 
 * @brief Summarize and output basic statistics on the performance of
 * a particular summarization or detection method across a number of
 * cel files.
 */

#include "chipstream/QuantMethodExprChipSummary.h"
//
#include "file/TsvFile/TsvFile.h"
#include "stats/stats.h"

using namespace std;
using namespace affx;

/**
 * @brief Sort pairs by the first coordinate, or the second
 * coordinate if the first coordinates are equal.
 */
bool QuantMethodExprChipSummary::CoordSort::operator() (const std::pair<float, float>* p1, 
                                                      const std::pair<float, float>* p2) const {
  if (fabs ((double)p1->first - (double)p2->first) > DBL_MIN)
    return p1->first < p2->first;
  else
    return p1->second < p2->second;
}


/** 
 * Constructor. If both the groupNames are empty then
 * all probesets are included in statistics.
 */
QuantMethodExprChipSummary::ProbeSetGroupStats::ProbeSetGroupStats(const std::string &name, uint32_t numChips, 
                                                                 bool cacheData, bool doThreshold,
                                                                 double minThreshold,
                                                                 double maxThreshold,
                                                                 const std::set<std::string> &groupNames) {
  assert(numChips);
  m_Name = name;
  m_ValidNames = groupNames;
  m_DoThreshold = doThreshold;
  m_MinThreshold = minThreshold;
  m_MaxThreshold = maxThreshold;
  m_PsCount = 0;
  m_AtomCount = 0;
  m_CacheData = cacheData;
  m_Summaries.resize(numChips);
  m_RLEs.resize(numChips);
  m_MADs.resize(numChips);
  m_Data.resize(numChips);
}

/** 
 * Check to see if we are including a particular ProbeSetGroup in
 * our statistics.
 * @param psGroup - ProbeSetGroup to check for inclusion.
 * @return bool - true if should be included, false otherwise.
 */
bool QuantMethodExprChipSummary::ProbeSetGroupStats::inGroup(ProbeSetGroup &psGroup) {
  if(m_ValidNames.empty())
    return true;
  if(psGroup.name != NULL && m_ValidNames.find(psGroup.name) != m_ValidNames.end())
    return true;
  return false;
}


/** 
 * Given a probeset group, extract and keep track of summary data
 * statistics.
 *
 * @param qMethod - The quantification method to extract statistics from.
 */
void QuantMethodExprChipSummary::ProbeSetGroupStats::report(ProbeSetGroup &psGroup, QuantExprMethod &qMethod) {
  if(inGroup(psGroup)) {
    uint32_t chipCount = qMethod.getNumTargets();
    std::vector<float> chipEffects;
    chipEffects.reserve(chipCount);
    
    /* Basic count of probes and probesets. */
    m_PsCount++;
    for(uint32_t i = 0; i < qMethod.getNumFeatures(); i++) {
      if(qMethod.featureUsed(i)) 
        m_AtomCount++;
    }
    
    /* Chip level statistics. */
    for(uint32_t chipIx = 0; chipIx < chipCount; chipIx++) {
      /* For present absent calls just do percent. */
      if(m_DoThreshold) {
        float val = qMethod.getTargetEffect(chipIx);
        if(m_CacheData) 
          m_Data[chipIx].push_back(val);
        if(val >= m_MinThreshold && val <= m_MaxThreshold) {
          m_Summaries[chipIx].addData(1);
        }
        else {
          m_Summaries[chipIx].addData(0);
        }
      }
      else {
        /* Basic value for each chip. */
        float chipVal = qMethod.getTargetEffect(chipIx);
        m_Summaries[chipIx].addData(chipVal);
        
        if(m_CacheData) 
          m_Data[chipIx].push_back(chipVal);
        
        /* Save for later as need to do medians. */
        if(qMethod.getScale() != QuantExprMethod::Log2) 
          chipVal = log2(chipVal);
        
        chipEffects.push_back(chipVal);
        std::vector<float> absResiduals;
        /* Do Median Absolute Deviation (MAD) of residuals if exist. */
        if(qMethod.haveResiduals()) {
          uint32_t featCount = qMethod.getNumFeatures();
          absResiduals.reserve(featCount);
          for(uint32_t featIx = 0; featIx < featCount; featIx++) {
            if(qMethod.featureUsed(featIx)) {
              absResiduals.push_back(fabs(qMethod.getResidual(featIx, chipIx)));
            }
          }
          float residMed = median_in_place(absResiduals.begin(), absResiduals.end());
          m_MADs[chipIx].addData(residMed);
        }
      }
    }
    if(!m_DoThreshold) {
      /* Loop through again to calculate median expression and RLE. */
      float chipMed = median(chipEffects.begin(), chipEffects.end());
      for(uint32_t chipIx = 0; chipIx < chipCount; chipIx++) {
        float RleVal = fabs(chipEffects[chipIx] - chipMed);
        m_RLEs[chipIx].addData(RleVal);
      }
    }
  }
}

/** Constructor. */
QuantMethodExprChipSummary::QuantMethodExprChipSummary(const std::vector<std::string> &chipNames,
                                                   bool doThreshold, double minThreshold, double maxThreshold,
                                                   std::string qccFile,
                                                   vector<MetaProbeset *> &metaSets) {
  setValid(false);
  std::set<std::string> groupNames;
  m_ChipNames = chipNames;
  m_posControlIx = -1;
  m_negControlIx = -1;
  m_DoThreshold = doThreshold;
  m_HavePosControls  = false;
  m_HaveNegControls = false;

  m_QcPsResults.resize(chipNames.size());

  set<string> psToRun;
  for(int i=0; i<metaSets.size(); i++) {
      psToRun.insert(ToStr(metaSets[i]->name));
  }

  /* Always add a group for all probesets, leaving both sets empty
     acts as wildcard. */
  addGroupStat("all_probeset", false, groupNames,
               doThreshold, minThreshold, maxThreshold,
               false, false);

  if(qccFile!="")
    readGroupsFile(qccFile, doThreshold, minThreshold, maxThreshold, psToRun);

  // Then AUC
  if (m_HavePosControls && m_HaveNegControls) {
    declareMetric("pos_vs_neg_auc", ChipSummary::Metric::Double);
  } 

  // Add qc probeset signal metrics
  for(map<string,string>::iterator sIx = m_QuantInHeader.begin(); sIx != m_QuantInHeader.end(); sIx++) {
      map<string,string>::iterator gIx = m_QuantInHeaderGroup.find(sIx->first);
      set<string>::iterator qcIx = psToRun.find(sIx->first);
      if(gIx != m_QuantInHeaderGroup.end() && qcIx != psToRun.end()) {
        //Verbose::out(1,"Adding metric '" + gIx->second + "-" + sIx->second + "'");
        declareMetric(gIx->second + "-" + sIx->second, ChipSummary::Metric::Double);
        m_QcPsNameVec.push_back(sIx->first);
      }
  }


}

/** 
 * Load up the groups from a simple tab delimited text file. Need to
 * have a column called 'group_name' and then at least one column
 * named 'probeset_id' or 'probeset_name'. If groups called
 * 'pos_control' and 'neg_control' are both supplied then the AUC
 * from the ROC is estimated and reported as well. Looks something
 * like:
 *
 * \verbatim
 group_name   probeset_name
 pos_control  AFFX-123
 pos_control  AFFX-456
 neg_control  AFFX-789
 \endverbatim
 * @param fileName - path to text file to read from.
 * @param doThreshold - Should we be doing thresholding (detection based)?
 * @param minThreshold - Minimum value to pass threshold.
 * @param maxThreshold - Maximum value to pass threshold.
 */
void QuantMethodExprChipSummary::readGroupsFile(const std::string &fileName, bool doThreshold, 
                                              double minThreshold, double maxThreshold,
                                              std::set<std::string> &psToRun) {
  TsvFile tsv;
  map<string, set<string> > groupSets;
  map<string, set<string> >::iterator mapIx;
  string groupName;
  int psIdInt = -1;
  int quantInHeader = 0;
  string psIdString;
  tsv.bind(0, "group_name", &groupName, TSV_BIND_REQUIRED);
  tsv.bind(0, "probeset_name", &psIdString, TSV_BIND_OPTIONAL);
  tsv.bind(0, "probeset_id", &psIdInt, TSV_BIND_OPTIONAL);
  tsv.bind(0, "quantification_in_header", &quantInHeader, TSV_BIND_OPTIONAL);
  if(tsv.open(fileName) != TSV_OK) 
    Err::errAbort("QuantMethodExprChipSummary::readGroupsFile() - Can't open file '" + fileName + "' to read.");
  if(tsv.cname2cidx(0, "probeset_name") < 0 && tsv.cname2cidx(0, "probeset_id") < 0) 
    Err::errAbort("QuantMethodExprChipSummary::readGroupsFile() - Must have either 'probeset_name' or 'probeset_id' column");
  while(tsv.nextLevel(0) == TSV_OK) {

    string psName;
    if(psIdInt >= 0) 
      psName = ToStr(psIdInt);
    else if(psIdString != "") 
      psName = psIdString;
    else
      Err::errAbort("Must have valid probeset_id or probeset_name");

    // Skip if we are not going to process this probeset
    set<string>::iterator qcIx = psToRun.find(psName);
    if(qcIx == psToRun.end())
        continue;

    mapIx = groupSets.find(groupName);

    /* Check to see if this is a new group. */
    if(mapIx == groupSets.end()) {
      set<string> newSet;
      newSet.insert(psName);
      groupSets[groupName] = newSet;
    }
    else {
      /* Already a group, just add info to the sets. */
      mapIx->second.insert(psName);
    }

    /* Check to see if we keep the probeset info around */
    if(quantInHeader == 1) {
      if(psIdString != "")
        m_QuantInHeader[psName] = psIdString;
      else
        m_QuantInHeader[psName] = psName;
      m_QuantInHeaderGroup[psName] = groupName;
    }
  }
  tsv.close();

  for (mapIx = groupSets.begin(); mapIx != groupSets.end(); ++mapIx) {
    bool isPosControl = false, isNegControl = false, cacheData = false;
    if(mapIx->first == "pos_control") {
      isPosControl = true;
      cacheData = true;
      m_HavePosControls = true;
    }
    else if(mapIx->first == "neg_control") {
      isNegControl = true;
      cacheData = true;
      m_HaveNegControls = true;
    }
    addGroupStat(mapIx->first, cacheData, mapIx->second,
                 doThreshold, minThreshold, maxThreshold, isPosControl, isNegControl);
  }
}

/** 
 * Add a group that we wish to keep statistics for. If a probeset's name is in
 * the name set the results for that probeset will be added to the summary
 * stats. If both the name set is empty it acts as a wildcard and
 * statistics are kept for every probeset.
 * 
 * @param name - Reference name for this group.
 * @param cacheData - Should we keep the data in RAM as we go? (i.e. for ROC curve).
 * @param names - Set containing names of probesets to be kept.
 * @param doThreshold - Should we be doing thresholding (detection based)?
 * @param minThreshold - Minimum value to pass threshold.
 * @param maxThreshold - Maximum value to pass threshold.
 * @param isPosControl - Should this group of stats be considered a positive control
 * @param isNegControl - Should this group of stats be considered a negative control
 */
void QuantMethodExprChipSummary::addGroupStat(const std::string &name, 
                                            bool cacheData,
                                            const std::set<std::string> &names,
                                            bool doThreshold, double minThresh, double maxThresh,
                                            bool isPosControl, bool isNegControl) {

  if(isPosControl && isNegControl) 
    Err::errAbort("QuantMethodExprChipSummary::addGroupStat() - Can't have both neg and pos control be true");
  else if(isPosControl) {
    m_posControlIx = m_Stats.size();
    cacheData = true;
  }
  else if(isNegControl) {
    m_negControlIx = m_Stats.size();
    cacheData = true;
  }
  m_Stats.push_back(ProbeSetGroupStats(name, m_ChipNames.size(), 
                                       cacheData, doThreshold, minThresh, maxThresh, names));

  declareMetric(name + "_probesets", ChipSummary::Metric::Integer);
  declareMetric(name + "_atoms", ChipSummary::Metric::Integer);

  if(doThreshold) {
    declareMetric(name + "_percent_called", ChipSummary::Metric::Double);
  } else {
    declareMetric(name + "_mean", ChipSummary::Metric::Double);
    declareMetric(name + "_stdev", ChipSummary::Metric::Double);
    declareMetric(name + "_mad_residual_mean", ChipSummary::Metric::Double);
    declareMetric(name + "_mad_residual_stdev", ChipSummary::Metric::Double);
    declareMetric(name + "_rle_mean", ChipSummary::Metric::Double);
    declareMetric(name + "_rle_stdev", ChipSummary::Metric::Double);
  }
}

bool QuantMethodExprChipSummary::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
    for(int i=0; i<iMart.getCelFileCount(); i++) {
        std::vector<ChipSummary::Metric> metrics;
        m_SummaryStats.push_back(metrics);
    }
    return true;
}

/** 
 * After every probeset computation this function is called an is an opportunity
 * to query the quantification method for results, residuals, etc.
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * @return true if success, false otherwise.
 */  
bool QuantMethodExprChipSummary::report(ProbeSetGroup &psGroup, QuantMethod &qMethod,
                                      const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                                      PmAdjuster &pmAdjust) {
  QuantExprMethod *qeMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
  if(qeMethod == NULL) {
    Err::errAbort("Can't call QuantMethodExprChipSummary::report() with something other than a QuantExprMethod.");
  }
  std::vector<ProbeSetGroupStats>::iterator statIx;
  for(statIx = m_Stats.begin(); statIx != m_Stats.end(); ++statIx) {
    statIx->report(psGroup, *qeMethod);
  }

  // deal with m_QuantInHeader
  map<string, string>::iterator sIx;
  sIx = m_QuantInHeader.find(psGroup.name);
  if(sIx != m_QuantInHeader.end()) {
    m_QcPsSeen.push_back(sIx->first);
    for(int chip=0; chip<qeMethod->getNumTargets(); chip++) {
        m_QcPsResults[chip].push_back(qeMethod->getSignalEstimate(chip));
    }
  }

  return true;
}

bool QuantMethodExprChipSummary::finish(QuantMethod &qMethod) {

    // Order of metrics:
    // - per catagory metrics
    // - pos vs neg auc
    // - various probeset signals for qc probesets

    // Load up other metrics computed in m_Stats
    for(int chip = 0; chip < m_SummaryStats.size(); chip++) {
        std::vector<ProbeSetGroupStats>::iterator statIx;
        for(statIx = m_Stats.begin(); statIx != m_Stats.end(); ++statIx) {
            m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_probesets", (int)statIx->m_PsCount));
            m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_atoms", (int)statIx->m_AtomCount));
            if(statIx->m_PsCount == 0) {
                // Dummy values because we failed to compute the expected probesets
                if(statIx->m_DoThreshold) {
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_percent_called", -1.0));
                } else {
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_mean", -1.0));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_stdev", -1.0));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_mad_residual_mean", -1.0));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_mad_residual_stdev", -1.0));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_rle_mean", -1.0));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_rle_stdev", -1.0));
                }
            } else {
                if(statIx->m_DoThreshold) {
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_percent_called", 
                                                statIx->m_Summaries[chip].getMean()));
                } else {
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_mean", 
                                                statIx->m_Summaries[chip].getMean()));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_stdev", 
                                                statIx->m_Summaries[chip].getStdev()));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_mad_residual_mean", 
                                                statIx->m_MADs[chip].getMean()));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_mad_residual_stdev", 
                                                statIx->m_MADs[chip].getStdev()));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_rle_mean", 
                                                statIx->m_RLEs[chip].getMean()));
                    m_SummaryStats[chip].push_back(ChipSummary::Metric(statIx->m_Name + "_rle_stdev", 
                                                statIx->m_RLEs[chip].getStdev()));
                }
            }
        }
    }

    // Compute AUC
    if(m_HavePosControls && m_HaveNegControls) {
        if(m_negControlIx >= 0 && m_posControlIx >= 0) {
            for(int chip=0; chip<m_SummaryStats.size(); chip++) {
                m_SummaryStats[chip].push_back(ChipSummary::Metric("pos_vs_neg_auc",
                    calcAuc(m_Stats[m_posControlIx].m_Data[chip], m_Stats[m_negControlIx].m_Data[chip], m_DoThreshold)));
            }
        } else {
            // Dummy values because we failed to compute the expected probesets
            for(int chip=0; chip<m_SummaryStats.size(); chip++) {
                m_SummaryStats[chip].push_back(ChipSummary::Metric("pos_vs_neg_auc",-1.0f));
            }
        }
    }

    // Add qc probeset signal metrics
    for(int qcPsIdx=0; qcPsIdx<m_QcPsNameVec.size(); qcPsIdx++) {
        // Do we have results for this probeset?
        int seenIdx = 0;
        bool found = false;
        for(seenIdx = 0; seenIdx < m_QcPsSeen.size(); seenIdx++) {
            if(m_QcPsSeen[seenIdx] == m_QcPsNameVec[qcPsIdx]) {
                found = true;
                break;
            }
        }

        for(int chip=0; chip<m_SummaryStats.size(); chip++) {
            std::string name = m_QuantInHeaderGroup[m_QcPsNameVec[qcPsIdx]] + "-" + 
                                    m_QuantInHeader[m_QcPsNameVec[qcPsIdx]];
            if(found) {
                m_SummaryStats[chip].push_back(ChipSummary::Metric(name, m_QcPsResults[chip][seenIdx]));
            } else {
                m_SummaryStats[chip].push_back(ChipSummary::Metric(name, -1.0));
            }
        }
    }

    setValid(true);

    return true;
}

/** 
 * Calculate the area under the ROC curve generated by the scores
 * for the positive and negative vectors. This code comes directly from
 * Peter's translation of the original exact perl code
 * @param pos - Vector of scores for positive examples.
 * @param neg - Vector of scores for negative examples.
 * @return - Estimate of area under the ROC curve.
 */
double QuantMethodExprChipSummary::calcAuc(std::vector<float> &posOrig, std::vector<float> &negOrig, bool increasing) {
  vector<float> pos(posOrig), neg(negOrig);
  if(!increasing) {
    vector<float>::iterator posIx, negIx;
    for(posIx = pos.begin(); posIx != pos.end(); ++posIx) {
      *posIx = -1 * (*posIx);
    }
    for(negIx = neg.begin(); negIx != neg.end(); ++negIx) {
      *negIx = -1 * (*negIx);
    }
  }
  sort (pos.begin(), pos.end());
  sort (neg.begin(), neg.end());
  vector<float> all (pos);
  all.insert (all.end(), neg.begin(), neg.end());
  sort (all.begin(), all.end());

  // Fun with C++: roc (x,y) coordinates are a set, sorted by first coordinate,
  // the second coordinate if the first are equal.
  typedef set<pair<float, float>*, CoordSort> coordSet;
  coordSet roc;
  // looks like the set will sometimes drop items that need to be deleted, keep track here
  vector<pair<float, float> *> toFree; 
  pair<float, float>* lowerLeft = new pair<float, float>;
  lowerLeft->first = lowerLeft->second = 0;
  roc.insert (lowerLeft);
  toFree.push_back(lowerLeft);
  pair<float, float>* upperRight = new pair<float, float>;
  upperRight->first = upperRight->second = 1;
  roc.insert (upperRight);
  toFree.push_back(upperRight);
  const unsigned int posSize = pos.size();
  const unsigned int negSize = neg.size();
  const float posSizeDbl = (float)posSize;
  const float negSizeDbl = (float)negSize;
  const unsigned int allSize = all.size();
  
  unsigned int pI = 0;
  unsigned int nI = 0;
  for(unsigned int j = 0; j < allSize; ++j) {
    for(unsigned int k = pI; k < posSize; ++k) {
      if (pos[k] <= all[j])
        pI = k;
      else
        break;
    }
    for(unsigned int k = nI; k < negSize; ++k) {
      if(neg[k] <= all[j])
        nI = k;
      else
        break;
    }
    pair<float, float>* point = new pair<float, float>;
    point->first = (float)(nI + 1) / negSizeDbl;
    point->second = (float)(pI + 1) / posSizeDbl;
    roc.insert (point);
    toFree.push_back(point);
  }
  
  double auc = 0;
  const coordSet::const_iterator rocEnd = roc.end();
  for(coordSet::const_iterator iter = roc.begin(); ; ++iter) {
    coordSet::const_iterator next = iter;
    // Process to the next to last entry.
    if(++next == rocEnd)
      break;
    const double w = (*next)->first - (*iter)->first;
    const double h1 = min ((*iter)->second, (*next)->second);
    const double h2 = max ((*iter)->second, (*next)->second);
    auc += w * h1 + 0.5 * w * (h2 - h1);
  }
  // Clean up.
  vector<pair<float, float> *>::iterator iter;
  for(iter = toFree.begin(); iter != toFree.end(); ++iter) {
    delete *iter;
  }
  return auc;
}
