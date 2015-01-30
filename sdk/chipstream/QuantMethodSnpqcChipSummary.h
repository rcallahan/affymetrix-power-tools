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

/// @file   QuantMethodSnpqcSummary.h
/// @brief  Class for accumulating some snpqc metrics from a list of probeset names.

#ifndef QUANTMETHODSNPQCSUMMARY_H
#define QUANTMETHODSNPQCSUMMARY_H

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
#include <string>
#include <vector>

class QuantMethodSnpqcChipSummary : public QuantMethodReport, public ChipSummary {

private:
  //
  enum CallType {
      ct_processed_count = 0,
      ct_hom,
      ct_het,
      ct_none,
      ct__size,
  };

public: 

  QuantMethodSnpqcChipSummary(const std::string& analysisName,
                          const std::string& metricPrefix) {
    setValid(false);
    m_AnalysisName = analysisName;
    m_MetricPrefix = metricPrefix;
    //
    setFormat(affx::TsvReport::FMT_TSV);
    setFilename("QuantMethodSnpqcSummary-DEBUG"); // for debugging, not used.

    // the count of probesetids in the input file.
    declareMetric(m_MetricPrefix+"_file_count",ChipSummary::Metric::Double);
    // count of probesets processed this snpqc
    declareMetric(m_MetricPrefix+"_processed_count",ChipSummary::Metric::Double);
    // % of probesets which have a call. (Not No-call)
    declareMetric(m_MetricPrefix+"_call_rate",ChipSummary::Metric::Double);
    // % of probesets which are AB
    declareMetric(m_MetricPrefix+"_het_rate",ChipSummary::Metric::Double);
    // % of probesets which are AA or BB
    declareMetric(m_MetricPrefix+"_hom_rate",ChipSummary::Metric::Double);
    //
    m_Calls.resize(ct__size);
  }

  void setProbesetNames(const std::vector<std::string>& names) {
    m_ProbesetNames=names;
    /// we keep the probeset names sorted for quick lookup
    sort(m_ProbesetNames.begin(),m_ProbesetNames.end());
  }
     
  bool prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
    // resize our m_Calls to match the number of chips.
    int chipCount = iMart.getCelFileCount();
    m_Calls.resize(ct__size);
    for (int i=0; i<m_Calls.size(); i++) {
      m_Calls[i].resize(chipCount);
    }
    return true;
  }

  bool report(ProbeSetGroup &psGroup,
              QuantMethod &qMethod,
              const IntensityMart &iMart,
              std::vector<ChipStream *> &iTrans, 
              PmAdjuster &pmAdjust) 
  {
    // check we have a QuantMethod which does genotyping.
    QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(&qMethod);
    if (gMethod == NULL) {
      Err::errAbort("Can only use a QuantMethodSnpqcReport with QuantGTypeMethods.");
    }

    // Is this probeset group part of our reporting set?
    std::string psName=gMethod->getProbeSetName();
    if (binary_search(m_ProbesetNames.begin(),m_ProbesetNames.end(),psName)==false) {
      // if not, stop processing it
      return true;
    }
    // yep!
    // Verbose::out(1,"QuantMethodSnpqcChipSummary::report()  psName='"+(psName)+"'");

    // accum the metrics.
    for (int i=0; i<gMethod->getNumCalls(); i++) {
      affx::GType call = gMethod->getGTypeCall(i);
      m_Calls[ct_processed_count][i]++;
      //
      if ((call==affx::AA) || (call==affx::BB)){
        m_Calls[ct_hom][i]++;
      }
      else if (call == affx::AB) {
        m_Calls[ct_het][i]++;
      }
      else if (call==affx::NN) {
        m_Calls[ct_none][i]++;
      }
      else {
        Err::errAbort("Not a good good call. (AA,AB,BB,NN)");
      }
    }
    return true;
  }

  bool reportFailure(ProbeSetGroup &psGroup, QuantMethod &qMethod, 
                     const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                     PmAdjuster &pmAdjust) {
    
    // do nothing.
    return true;
  }
  
  bool finish(QuantMethod &qMethod) {
    for (size_t chip = 0; chip < m_Calls[ct_processed_count].size(); chip++) {
      // cast to double for math.
      double processed_count=m_Calls[ct_processed_count][chip];
      //
      setMetric(chip,m_MetricPrefix+"_file_count",
                (double(m_ProbesetNames.size())));
      setMetric(chip,m_MetricPrefix+"_processed_count",
                processed_count);
      setMetric(chip,m_MetricPrefix+"_call_rate",
                (double(100.0*(m_Calls[ct_het][chip]+m_Calls[ct_hom][chip])/processed_count)));
      setMetric(chip,m_MetricPrefix+"_het_rate",
                (double(100.0*m_Calls[ct_het][chip]/processed_count)));
      setMetric(chip,m_MetricPrefix+"_hom_rate",
                (double(100.0*m_Calls[ct_hom][chip]/processed_count)));
    }
    //
    setValid(true);
    return true;
  }

private:
  /// The name of our analysis.
  std::string m_AnalysisName;
  /// The prefix of our qc calls. (SNPQC)
  std::string m_MetricPrefix;
  ///
  std::vector<std::string> m_ProbesetNames;
  /// Keep track of call counts per chip. Indexes are CallType enum
  std::vector<std::vector<int> > m_Calls;
};

#endif /* QUANTMETHODSNPQCSUMMARY_H */
