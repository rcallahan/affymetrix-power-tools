////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
 * @file   QuantMethodGTypeChipSummary.h
 * @author Alan Williams
 * 
 * @brief  Class for accumulating some run statistics from a genotyping run.
 * 
 */

#ifndef QUANTMETHODGTYPECHIPSUMMARY_H
#define QUANTMETHODGTYPECHIPSUMMARY_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/ChipSummary.h"
#include "chipstream/GenderCalls.h"
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethodReport.h"
#include "chipstream/SpecialSnps.h"
#include "chipstream/SummaryStats.h"
//
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <cfloat>
#include <iostream>
#include <vector>
//

class QuantMethodGTypeChipSummary : public QuantMethodReport, public ChipSummary {
private:
  enum CallType {
      total = 0,
      auto_total = 1,
      hom   = 2,
      auto_hom = 3,
      het   = 4,
      auto_het = 5,
      nc    = 6,
      auto_nc = 7,
      pra   = 8,
      cn0   = 9,
      
  };

  static const int CallTypeSize = 10;


public: 

  QuantMethodGTypeChipSummary(const std::string& analysisName,
                              GenderCalls *genders,
                              SpecialSnpMap& SpecialSnps,
                              bool alleleSummariesOnly = false,
                              double minThreshold = -1 * DBL_MAX, 
                              double maxThreshold = DBL_MAX) {
    setValid(false);
    m_AnalysisName = analysisName;
    m_SpecialSnps = SpecialSnps;
    m_MinThreshold = minThreshold;
    m_MaxThreshold = maxThreshold;
    m_Genders = genders;
    m_AlleleSummariesOnly = alleleSummariesOnly;
    // set to reduce debugging warnings.
    setFormat(affx::TsvReport::FMT_TSV);
    setFilename("QuantMethodGTypeChipSummary-DEBUG"); // for debugging, not used.

    if(m_Genders->getGenders().empty() == false) {
        declareMetric("computed_gender",ChipSummary::Metric::String);
    }
    declareMetric("call_rate",ChipSummary::Metric::Double);
    declareMetric("total_call_rate",ChipSummary::Metric::Double);
    declareMetric("het_rate",ChipSummary::Metric::Double);
    declareMetric("total_het_rate",ChipSummary::Metric::Double);
    declareMetric("hom_rate",ChipSummary::Metric::Double);
    declareMetric("total_hom_rate",ChipSummary::Metric::Double);
    declareMetric("cluster_distance_mean",ChipSummary::Metric::Double);
    declareMetric("cluster_distance_stdev",ChipSummary::Metric::Double);
  }

  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param qMethod - Quantification method to be used.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  bool prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
    bool allOk = true;
    int chipCount = iMart.getCelFileCount();
    for(int i=0; i<chipCount; i++) {
        std::vector<ChipSummary::Metric> metrics;
        m_SummaryStats.push_back(metrics);
    }
    m_Calls.resize(CallTypeSize);
    for(unsigned int i = 0; i < m_Calls.size(); i++) {
      m_Calls[i].resize(chipCount);
    }
    m_MeanClosestDistance.resize(chipCount);
    return allOk;
  }

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
  bool report(ProbeSetGroup &psGroup, QuantMethod &qMethod,
              const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
              PmAdjuster &pmAdjust) {
    QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(&qMethod);
    bool realCalls = true;
    if(gMethod == NULL) {
      Err::errAbort("Can only use a QuantMethodGTypeReport with QuantGTypeMethods.");
    }
    Err::check(gMethod->getNumCalls() == m_Calls[0].size(), 
               "Error: QuantMethodGTypeChipSummary::report() - Wrong number of calls.");
    std::vector<double> distanceValues(gMethod->getNumCalls());
    /* Keep track of the calls. */
    const std::string probeSetName = ToStr(psGroup.name);
    bool autosomal = !isHaploid(m_SpecialSnps, probeSetName);
    for(unsigned int i = 0; i < gMethod->getNumCalls(); i++) {
      ///@todo Nice to get away from notion of GTypeCall. This is a hack which works with QuantLabelZMulti
      ///      because we are only interested in if the call is a hom or het.
      affx::GType forcedCall = gMethod->getGTypeForcedCall(i);
      double conf = gMethod->getConfidence(i);
      float dist = gMethod->getConfGenotype(forcedCall, i);
      if(!Util::isFinite(conf) || !Util::isFinite(dist)) {
        Verbose::out(1, "Warning! - Non-finite (Nan/Inf) confidence produced in probeset: '" + probeSetName + 
                     "' at column: " + ToStr(i));
      }
      if(dist == FLT_MAX) 
        realCalls = false;
      distanceValues[i] = dist;


      QuantLabelZMulti *qmm = dynamic_cast<QuantLabelZMulti *>(&qMethod);
      if(qmm != NULL) {
            int call = qmm->getCall(i);
            GenoCallCoder *coder = qmm->m_coder;
            if(coder == NULL) 
                Err::errAbort("GenoCallCoder is NULL");
            if(coder->abstractAlleleToGenotypeCallNum("NoCall") == call) {
                m_Calls[nc][i]++;
                if (autosomal) {
                    m_Calls[auto_nc][i]++;
                }
            } 
            else if(coder->abstractAlleleToGenotypeCallNum("PossibleRareAllele") == call) {
                m_Calls[pra][i]++;
            } 
            else if(coder->abstractAlleleToGenotypeCallNum("ZeroCopyNumber") == call) {
                m_Calls[cn0][i]++;
            } 
            else if(coder->isHom(call)) {
                m_Calls[hom][i]++;
                if (autosomal) {
                    m_Calls[auto_hom][i]++;
                }
            } 
            else if(coder->isHet(call)) {
                m_Calls[het][i]++;
                if (autosomal) {
                    m_Calls[auto_het][i]++;
                }
            } 
            else {
                Verbose::out(1,"Unaccounted for call (" + ToStr(call) + "). Treating as no-call for summary stats.");
                m_Calls[nc][i]++;

            }
      } 
      else {
            affx::GType call = gMethod->getGTypeCall(i);
            if (call == affx::AA) {
                m_Calls[hom][i]++;
                if (autosomal) {
                    m_Calls[auto_hom][i]++;
                }
            } else if(call == affx::AB) {
                m_Calls[het][i]++;
                if (autosomal) {
                    m_Calls[auto_het][i]++;
                }
            } else if(call == affx::BB) {
                m_Calls[hom][i]++;
                if (autosomal) {
                    m_Calls[auto_hom][i]++;
                }
            } else if(call == affx::NN) {
                m_Calls[nc][i]++;
                if (autosomal) {
                    m_Calls[auto_nc][i]++;
                }
            } else {
                Verbose::out(1,"Unaccounted for call (" + ToStr(call) + "). Treating as no-call for summary stats.");
                m_Calls[nc][i]++;
                if (autosomal) {
                    m_Calls[auto_nc][i]++;
                }
            }
      }
      m_Calls[total][i]++;
      if (autosomal) {
        m_Calls[auto_total][i]++;        
      }
    }
    /* if there were real predictions then keep track of them. */
    if(realCalls) {
      m_MeanClosestDistance.addData(distanceValues);
    }

    return true;
  }

  /** 
   * If a probeset fails to compute for whatever reason then this method is
   * called rather than the normal report call above. By default does nothing.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
   bool reportFailure(ProbeSetGroup &psGroup, QuantMethod &qMethod, 
                      const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                      PmAdjuster &pmAdjust) {

     // do nothing.
     return true;
   }

  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * @param qMethod - Quantification method that was used.
   * @return true if success, false otherwise.
   */
  bool finish(QuantMethod &qMethod) {

    vector<affx::Gender> genders(0);
    if(m_Genders != NULL) {
        genders = m_Genders->getGenders();
    }

    // Load up other metrics computed in m_Stats
    for(size_t chip = 0; chip < m_SummaryStats.size(); chip++) {

        /// @todo should this be changed to m_AnalysisName + "_gender_used"? 
        if (genders.empty() == false) {
            m_SummaryStats[chip].push_back(ChipSummary::Metric("computed_gender", affx::getGenderString(genders[chip])));
        }
        // was added -- but now removing to reduce redundancy. analysis can be had from the
        // meta header info (which previously was not provided)
        //m_SummaryStats[chip].push_back(ChipSummary::Metric(m_AnalysisName + "_call_rate", getCallRate(chip)));
        m_SummaryStats[chip].push_back(ChipSummary::Metric("call_rate", getAutoCallRate(chip)));
        m_SummaryStats[chip].push_back(ChipSummary::Metric("total_call_rate", getTotalCallRate(chip)));
        m_SummaryStats[chip].push_back(ChipSummary::Metric("het_rate", getAutoHetCallRate(chip)));
        m_SummaryStats[chip].push_back(ChipSummary::Metric("total_het_rate", getTotalHetCallRate(chip)));
        m_SummaryStats[chip].push_back(ChipSummary::Metric("hom_rate", getAutoHomCallRate(chip)));
        m_SummaryStats[chip].push_back(ChipSummary::Metric("total_hom_rate", getTotalHomCallRate(chip)));
        m_SummaryStats[chip].push_back(ChipSummary::Metric("cluster_distance_mean", m_MeanClosestDistance.getMean(chip)));
        // this little hack is here to suppress warning messages that would
        // probably vex some users.
        double distanceStdDev = (m_AlleleSummariesOnly ? 0 : m_MeanClosestDistance.getStdev(chip));
        m_SummaryStats[chip].push_back(ChipSummary::Metric("cluster_distance_stdev", distanceStdDev));
    }

    // @todo: this shouldnt be needed -jhg
    // zero out the header buffering (We have headers pushed ontous which we never use.
    m_header_buffer.resize(0);
    
    //
    setValid(true);
    return true;
  }

private: 

  /** Get the autosomal call rate */
  double getAutoCallRate(int index) {
    return 100.0 * (m_Calls[auto_total][index] - m_Calls[auto_nc][index]) / m_Calls[auto_total][index];
  }

  /** Get the total call rate */
  double getTotalCallRate(int index) {
    return 100.0 * (m_Calls[total][index] - m_Calls[nc][index]) / m_Calls[total][index];
  }

  /** Get the autosomal het call rate */
  double getAutoHetCallRate(int index) {
    return 100.0 * (m_Calls[auto_het][index]) / (m_Calls[auto_total][index]);
  }

  /** Get the total het call rate */
  double getTotalHetCallRate(int index) {
    return 100.0 * (m_Calls[het][index]) / (m_Calls[total][index]);
  }

  /** Get the autosomal hom call rate */
  double getAutoHomCallRate(int index) {
    return 100.0 * (m_Calls[auto_hom][index]) / (m_Calls[auto_total][index]);
  }

  /** Get the total hom call rate */
  double getTotalHomCallRate(int index) {
    return 100.0 * (m_Calls[hom][index]) / (m_Calls[total][index]);
  }

private:
  /// The name of our analysis.
  std::string m_AnalysisName;
  /// What level of confidence is needed to make a call?
  double m_MinThreshold;
  double m_MaxThreshold;
  /// Keep track of call counts per chip. Indexes are CallType enum
  std::vector<std::vector<int> > m_Calls;
  /// Distance to nearest cluster center.
  SummaryStats m_MeanClosestDistance; 
  /// Gender of sample.
  GenderCalls *m_Genders;
  /// map of haloid snp names
  SpecialSnpMap m_SpecialSnps;
  bool m_AlleleSummariesOnly;

//   1) mean/stdev median absolute deviation
//   2) allele summary.

};

#endif /* QUANTMETHODGTYPECHIPSUMMARY_H */
