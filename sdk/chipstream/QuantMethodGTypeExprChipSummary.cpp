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
 * @file   QuantMethodGTypeExprChipSummary.cpp
 * @author Alan Williams
 *
 * @brief  Class for accumulating some run statistics from a
 *         genotyping run from the expression quant method
 *
 */


#include "chipstream/QuantMethodGTypeExprChipSummary.h"

QuantMethodGTypeExprChipSummary::QuantMethodGTypeExprChipSummary(const std::string& analysisName) {
  setValid(false);
  m_AnalysisName = analysisName;
  // these are set to reduce the debuging warnings.
  setFilename("QuantMethodGTypeExprChipSummary-DEBUG"); // not used.
  setFormat(affx::TsvReport::FMT_TSV); // not used.

  declareMetric("allele_summarization_mean",ChipSummary::Metric::Double);
  declareMetric("allele_summarization_stdev",ChipSummary::Metric::Double);
  declareMetric("allele_deviation_mean",ChipSummary::Metric::Double);
  declareMetric("allele_deviation_stdev",ChipSummary::Metric::Double);
  declareMetric("allele_mad_residuals_mean",ChipSummary::Metric::Double);
  declareMetric("allele_mad_residuals_stdev",ChipSummary::Metric::Double);
}

/** 
 * Get set up for a run of reporting probesets. Often used to open file
 * streams and print headers to files etc.
 * 
 * @param qMethod - Quantification method to be used.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeExprChipSummary::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
  bool allOk = true;
  int chipCount = iMart.getCelFileCount();
  for(int i=0; i<chipCount; i++) {
    std::vector<ChipSummary::Metric> metrics;
    m_SummaryStats.push_back(metrics);
  }
  //    m_IntensitySummary.resize(chipCount);
  m_AlleleSummary.resize(chipCount);
  m_MadResidualSummary.resize(chipCount);
  m_AbsDeviation.resize(chipCount);
  return allOk;
}

/** 
 * After every probeset computation this function is called an is an opportunity
 * to query the quantification method for results, residuals, etc.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeExprChipSummary::report(ProbeSetGroup &psGroup, QuantMethod &qMethod,
                                             const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                                             PmAdjuster &pmAdjust) {
  QuantExprMethod *qeMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
  if(qeMethod == NULL) {
    Err::errAbort("QuantMethodGTypeExprChipSummary::report() Can only use a QuantExprMethod object.");
  }
  addAlleleSummary(*qeMethod);
  addMadResiduals(*qeMethod);
  return true;
}

/** 
 * If a probeset fails to compute for whatever reason then this method is
 * called rather than the normal report call above. By default does nothing.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeExprChipSummary::reportFailure(ProbeSetGroup &psGroup, QuantMethod &qMethod,
                                                    const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                                                    PmAdjuster &pmAdjust) {

  // do nothing.
  return true;
}

/** 
 * Add the summary statistics for the raw intensity of probes in this probeset group.
 *  
 * @param psGroup - Probes to be added.
 * @param iMart - Data to get intensity from.
 */
void QuantMethodGTypeExprChipSummary::addIntensityValues(ProbeSetGroup &psGroup, const IntensityMart &iMart) {
  int chipCount = iMart.getCelFileCount();
  std::vector<double> values(chipCount);
  for(unsigned int groupIx = 0; groupIx < psGroup.probeSets.size(); groupIx++) {
    const ProbeSet *ps = psGroup.probeSets[groupIx];
    for(unsigned int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
      Atom *a = ps->atoms[atomIx];
      for(unsigned int probeIx = 0; probeIx < a->probes.size(); probeIx++) {
        Probe *p = a->probes[probeIx];
        for(unsigned int i = 0; i < chipCount; i++) {
          values[i] = log2(max((double)iMart.getProbeIntensity(p->id, i), DBL_MIN));
        } /* chips. */
        m_IntensitySummary.addData(values);
      } /* probes. */
    } /* atoms. */
  } /* probesets. */
}
  
void QuantMethodGTypeExprChipSummary::addAlleleSummary(QuantExprMethod &qMethod) {
  std::vector<double> values;
  std::vector<double> deviation;
  for(unsigned int i = 0; i < qMethod.getNumTargets(); i++) {
    values.push_back(log2(max(qMethod.getTargetEffect(i), DBL_MIN)));
  }

  m_AlleleSummary.addData(values);
  double med = median(values.begin(), values.end());
  for(unsigned int i = 0; i < values.size(); i++) {
    deviation.push_back(fabs(med - values[i]));
  }
  m_AbsDeviation.addData(deviation);
}

void QuantMethodGTypeExprChipSummary::addMadResiduals(QuantExprMethod &qMethod) {
  static bool warned = false;
  bool useResidual = qMethod.haveResiduals();
  std::vector<double> values(qMethod.getNumTargets());
  std::vector<double> chipVals(qMethod.getNumFeatures());
  if(!useResidual && !warned) {
    Verbose::out(1, "Warning: qMethod: " + ToStr(qMethod.getType()) + " doesn't have residuals, using zero.");
    warned = true;
  }
  for(unsigned int chipIx = 0; chipIx < qMethod.getNumTargets(); chipIx++) {
    for(unsigned int featIx = 0; featIx < qMethod.getNumFeatures(); featIx++) {
      if(useResidual)
        chipVals[featIx] = fabs(qMethod.getResidual(featIx, chipIx));
      else
        chipVals[featIx] = 0;
    }
    values[chipIx] = median(chipVals.begin(), chipVals.end());
  }
  m_MadResidualSummary.addData(values);
}

/** 
 * No more probesets will be processed, this is a chance to finish outputting
 * results and clean up.
 * @param qMethod - Quantification method that was used.
 * @return true if success, false otherwise.
 */
bool QuantMethodGTypeExprChipSummary::finish(QuantMethod &qMethod) {

  // Load up other metrics computed in m_Stats
  for(int chip = 0; chip < m_SummaryStats.size(); chip++) {

    m_SummaryStats[chip].push_back(ChipSummary::Metric("allele_summarization_mean", m_AlleleSummary.getMean(chip)));
    m_SummaryStats[chip].push_back(ChipSummary::Metric("allele_summarization_stdev", m_AlleleSummary.getStdev(chip)));
    m_SummaryStats[chip].push_back(ChipSummary::Metric("allele_deviation_mean", m_AbsDeviation.getMean(chip)));
    m_SummaryStats[chip].push_back(ChipSummary::Metric("allele_deviation_stdev", m_AbsDeviation.getStdev(chip)));
    m_SummaryStats[chip].push_back(ChipSummary::Metric("allele_mad_residuals_mean", m_MadResidualSummary.getMean(chip)));
    m_SummaryStats[chip].push_back(ChipSummary::Metric("allele_mad_residuals_stdev", m_MadResidualSummary.getStdev(chip)));
  }

  setValid(true);
  return true;
}


