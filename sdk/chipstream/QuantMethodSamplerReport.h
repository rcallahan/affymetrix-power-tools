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
 * @file   QuantMethodSamplerReport.h
 * @author Chuck Sugnet
 * @date   Wed Feb  1 09:09:36 2006
 * 
 * @brief  Class for reporting qc stats to a file.
 */
#ifndef _QUANTMETHODSAMPLERREPORT_H_
#define _QUANTMETHODSAMPLERREPORT_H_

//
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantMethodReport.h"
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cstring>
#include <string>
//
class QuantMethodSamplerReport : public QuantMethodReport {
public:

  /** 
   * Constructor.
   *  
   * @param reportModulus - Sampling rate for how often a report is
   * generated. For example 9 means that 1 out of every 10 probeset
   * groups will be reported on.
   * @param prefix - Prefix name for report files, usually something like experimentDir/rma
   * @param precision - How many digits to print after the decimal place.
   * @param doRawIntensities - Print a sampling report for raw probe intensities.
   * @param doSummaries - Print a sampling for general summaries.
   * @param doMadResiduals - Print a report for median absolute residuals.
   */
  QuantMethodSamplerReport(int reportModulus, 
                           const std::string& prefix,
                           int precision,
                           bool doRawIntensities,
                           bool doSummaries,
                           bool doMadResiduals) : 
    m_Modulus(reportModulus),
    m_Prefix(prefix),
    m_Precision(precision),
    m_DoRawIntensity(doRawIntensities),
    m_DoSummaries(doSummaries),
    m_DoMadResiduals(doMadResiduals) {
    m_ModCount = 0;
    m_Delim = '\t';
    m_IntensityCount = 0; 
    m_IntensityModulus = m_Modulus;
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
    QuantExprMethod *qeMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
    if(qeMethod == NULL)
      Err::errAbort("Can't have a QuantMethodSamplerReport object called with things other than QuantExprMethod objects");
    if(m_DoRawIntensity)
      Fs::mustOpenToWrite(m_RawOut, m_Prefix + ".intensities.qc.txt");
    if(m_DoSummaries)
      Fs::mustOpenToWrite(m_SummariesOut, m_Prefix + ".summaries.qc.txt");
    if(m_DoMadResiduals && qeMethod->haveResiduals()) {
      Fs::mustOpenToWrite(m_MadOut, m_Prefix + ".mad-residuals.qc.txt");
      m_MadOut.setf(ios::fixed, ios::floatfield);
      m_MadOut.precision(m_Precision);
    }
    else
      m_DoMadResiduals = false;
    
    return true;
  }
  
  /** 
   * Print out the median absolute residual (MAD) for probeset group.
   * 
   * @param out - Stream we are printing to.
   * @param psGroup - Group containing probes to print out.
   * @param qMethod - Quantification method to query for residuals.
   */
  void reportMadResiduals(std::ofstream &out, ProbeSetGroup &psGroup, 
                          QuantExprMethod &qMethod) {
    out << psGroup.name;
    std::vector<float> residuals;
    residuals.reserve(qMethod.getNumTargets());

    for(size_t chipIx = 0; chipIx < qMethod.getNumTargets(); chipIx++) {
      out.put(m_Delim);
      residuals.clear();
      /* Get all the absolute values of the residuals in a vector. */
      for(size_t featureIx = 0; featureIx < qMethod.getNumFeatures(); featureIx++) {
        residuals.push_back(fabs(float(qMethod.getResidual(featureIx, chipIx))));
      }
      /* median residual for a particular chip is what we report. */
      double mad = median(residuals.begin(), residuals.end());
      out << mad;
    }
    out << std::endl;
  }
  
  /** 
   * Report the summary value for a particular probeset group.
   * 
   * @param out - Stream we are printing to.
   * @param psGroup - Group of probesets used to make summary.
   * @param qMethod - Quantification method to query for summaries.
   */
  void reportSummaries(std::ofstream &out, ProbeSetGroup &psGroup, 
                       QuantExprMethod &qMethod) {
    out << psGroup.name;
    for(size_t chipIx = 0; chipIx < qMethod.getNumTargets(); chipIx++) {
      out.put(m_Delim);
      double summary = qMethod.getTargetEffect(chipIx);
      out << summary;
    }
    out << std::endl;
  }

  /** 
   * Print out the log2() intensities of PM probes in psGroup.
   * 
   * @param out - Stream we are printing to.
   * @param psGroup - Group of probesets that have PM probes to report.
   * @param iMart - Table to look up probe intensities in.
   */
  void reportIntensities(std::ofstream &out, ProbeSetGroup &psGroup, const IntensityMart &iMart) {
    unsigned int psIx = 0, atomIx = 0, probeIx = 0, chipIx = 0;
    for(psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
      const ProbeSet *ps = psGroup.probeSets[psIx];
      for(atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
        Atom *a = ps->atoms[atomIx];
	unsigned int channelIx = a->getChannelCode();
        for(probeIx = 0; probeIx < a->probes.size(); probeIx++) {
          Probe *p = a->probes[probeIx];
          if(Probe::isPm(*p)) {
            out << psGroup.name << "." << p->id;
            for(chipIx = 0; chipIx < iMart.getCelFileCount(); chipIx++) {
              out.put(m_Delim);
              out << log2(max(iMart.getProbeIntensity(p->id, chipIx, channelIx),FLT_MIN));
            } 
            out << std::endl;
          } /* if pm probe. */
        } /* probes. */
      } /* atoms. */
    } /* probesets. */
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
  bool report(ProbeSetGroup &psGroup, 
              QuantMethod &qMethod,
              const IntensityMart &iMart, 
              std::vector<ChipStream *> &iTrans, 
              PmAdjuster &pmAdjust) {
    QuantExprMethod *qeMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
    if(qeMethod == NULL)
      Err::errAbort("Can't have a QuantMethodSamplerReport object called with things other than QuantExprMethod objects");
    /* Only report some probesets, specifically total/(m_Modulus+1) of all probesets. */
    if(m_ModCount-- <= 0) {
      if(m_DoMadResiduals) 
        reportMadResiduals(m_MadOut, psGroup, *qeMethod);
      if(m_DoRawIntensity) {
        /* We're not reporting all probes for all sampled probesets, would be too many. */
        if(m_IntensityCount-- <= 0) {
          reportIntensities(m_RawOut, psGroup, iMart);
          m_IntensityCount = m_IntensityModulus;
        }
      }
      if(m_DoSummaries) {
        reportSummaries(m_SummariesOut, psGroup, *qeMethod);
      }
      m_ModCount = m_Modulus;
    }
    return true;
  };
  
  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * 
   * @param qMethod - Quantification method that was used.
   * 
   * @return true if success, false otherwise.
   */
  bool finish(QuantMethod &qMethod) {
    if(m_RawOut.is_open()) 
      m_RawOut.close();
    if(m_SummariesOut.is_open()) 
      m_SummariesOut.close();
    if(m_MadOut.is_open()) 
      m_MadOut.close();
    return true;
  }

private:

  int m_Modulus; ///< How often do we sample from reports given (i.e. 9 means 1 in 10 probesets gets reported.
  std::string m_Prefix; ///< Prefix for filenames.
  int m_Precision; ///< How many digits after the decimal place to output.
  bool m_DoRawIntensity; ///< Output raw intensities for probesets seen?
  bool m_DoSummaries;    ///< Output summarized values for probesets seen?
  bool m_DoMadResiduals; ///< Output median absolute residual for each probeset seen?
  char m_Delim; ///< What delimiter to use for fields, currently a tab.
  int m_ModCount; ///< How many probesets have we seen since last sampling.
  int m_IntensityModulus; ///< How often are we going to print the raw intensities of an entire probe set?
  int m_IntensityCount;   ///< Counter for how many probesets we have seen since last intensity sampling.
  std::vector<std::string> m_ColNames; ///< Names of columns
  std::ofstream m_RawOut;  ///< File for outputting raw values.
  std::ofstream m_SummariesOut; ///< File for outputting summaries.
  std::ofstream m_MadOut;       ///< File for outputting MAD values.
};

#endif /* _QUANTMETHODSAMPLERREPORT_H_ */
