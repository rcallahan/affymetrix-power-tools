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
 * @file   MedNormTran.cpp
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:42:26 2005
 * 
 * @brief  Class for doing median, or average, normalization
 */

//
#include "chipstream/MedNormTran.h"
//
#include "stats/stats.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Verbose.h"

using namespace std;

/** 
 * Constructor
 * 
 * @param target - Value to set median or average value to.
 * @param doAverage - Do average instead of median.
 * @param calcTarget - Calculate target as average or median of all data.
 * @param lowPrecision - Should we pretend like this got wrote to a cel file and 
 *                       read back again (meaning truncation)
 */
MedNormTran::MedNormTran(float target, bool doAverage, bool calcTarget, bool lowPrecision) :
  m_TargetNorm(target), m_DoAverage(doAverage), m_TargetUnset(calcTarget),
  m_LowPrecision(lowPrecision) {
  m_SubProbeCount = 0;
  m_Type = MEDNORMSTR;
  if(calcTarget == true && target > 0.0) {
    Err::errAbort("MedNormTran::MedNormTran() - Can't set target and calcTarget at same time.");
  }
  setupSelfDoc(*this);
  setOptValue("target", ToStr(m_TargetNorm));
  setOptValue("doavg", m_DoAverage);
  setOptValue("lowprecision", m_LowPrecision);
  setOptValue("calctarget", m_TargetUnset);
  }


/** 
 * @brief transform the intensity point supplied coming from a particular
 * probe in a particular microarray.
 * 
 * @param probeIx - Probe index from the cel file.
 * @param chipIx - Set of chip indexes from same sample.
 * @param intensity - Original intensities. 
 * @param return - transformed intensities.
 */
float MedNormTran::transform(int probeIx, int chipIx, float intensity) {
  intensity = float(intensity * m_Scale[chipIx]);
  if (m_LowPrecision) { 
    intensity = Convert::floatLowPrecision(intensity); 
  }
  return intensity;
}

void MedNormTran::transform(int chipIx, std::vector<float>& intensity) {
  for (int probeIx = 0; probeIx < intensity.size(); probeIx++) {
    intensity[probeIx] = transform(probeIx, chipIx, intensity[probeIx]);
  }
}

/** 
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> MedNormTran::getDefaultDocOptions() { 
  std::vector<SelfDoc::Opt> opts;
  SelfDoc::Opt target = {"target", SelfDoc::Opt::Double, "0.0", "0.0", "0", "NA", 
                         "Target intensity to set all chips median (or average) to."};
  opts.push_back(target);
  SelfDoc::Opt doavg = {"doavg", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                        "Set to true to do average rather than median."};
  opts.push_back(doavg);
  SelfDoc::Opt calc = {"calctarget", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                       "Calculate a target from median of chips."};
  opts.push_back(calc);
  SelfDoc::Opt lowprecision = {"lowprecision", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                               "Set to 'true' to truncate values as seen when writing results to a normalized cel file."};
  opts.push_back(lowprecision);
  return opts;
}

void MedNormTran::setupSelfDoc(SelfDoc &doc) {
  doc.setDocName(MEDNORMSTR);
  doc.setDocDescription("Class for doing median normalization. Adjust intensities such that all chips have the same median (or average).");
  doc.setDocOptions(getDefaultDocOptions());
}

/** 
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc MedNormTran::explainSelf() {
  SelfDoc doc;
  setupSelfDoc(doc);
  return doc;
}

/** 
 * @brief Only use a subset of all probes for normalization. For
 * example RMA only uses PM probes or might want to just use control
 * genes.
 * @param subsetProbes - bitmask indicating which probes to use for
 * normalization.
 */
void MedNormTran::setSubProbes(const std::vector<bool> &subsetProbes) {
    unsigned int i = 0;
    m_SubProbes = subsetProbes;
    m_SubProbeCount = 0;
    for(i = 0; i < m_SubProbes.size(); i++) 
      if(m_SubProbes[i] == true) 
        m_SubProbeCount++;
  } 


void MedNormTran::newChip(std::vector<float> &data) {
   initializeData(data);
}

/** 
 * @brief Method for being passed a new cel file worth of data.
 * @param data - cel file vector of data.
 */  
void MedNormTran::newDataSet(IntensityMart* iMart) {
  int dataSetCount = iMart->getCelDataSetCount();
  m_TransformedIMart = iMart;
  // if there aren't any chipstream nodes after this, then don't store
  // all of the intensities
  if (m_Streams.empty()) {
    m_TransformedIMart->setStoreAllCelIntensities(false);
  }

  std::vector<float> data;
  for (int d = 0; d < dataSetCount; d++) {
    data = iMart->getCelData(d);
    newChip(data);
    if(!m_TargetUnset) {
      transform(d, data);
      m_TransformedIMart->setProbeIntensity(d, data);
    }
  }
  finishedChips();
  if(m_TargetUnset) {
    for (int d = 0; d < dataSetCount; d++) {
      data = iMart->getCelData(d);
      transform(d, data);
      m_TransformedIMart->setProbeIntensity(d, data);
    }
  }
  chipStreamPassNewChip(m_TransformedIMart);
}

/** 
 * @brief Method for being passed a new cel file worth of data. Calculates
 * and stores the median (or average) of data supplied.
 * @param data - cel file vector of data.
 */  
void MedNormTran::initializeData(const std::vector<float> &data) {
  vector<float> copy(m_SubProbeCount);
  unsigned int i = 0, currentPm = 0;
  float summary = 0;
  
  /* Sanity check. */
  if(m_SubProbes.size() > 0 && data.size() != m_SubProbes.size()) 
    Err::errAbort("MedNormTran::newChip() - Chip data different size than subset probe vector");
  
  /* If we're doing a subset of the data, put subset into copy and
     normalize from copy. */
  if(m_SubProbes.size() > 0) {
    /* First filter out PM probes, for RMA only use pm probes. */
    for(i = 0; i < m_SubProbes.size(); i++) {
      if(m_SubProbes[i] == true) {
        copy[currentPm++] = data[i];
      }
    }
    /* Get summary from copy. */
    if(m_DoAverage)
      summary = float(average(copy.begin(),copy.end()));
    else {
      summary = median_in_place(copy.begin(),copy.end());
       Verbose::out(2, "Median for chip: " + ToStr(m_Summary.size()) + " is: " + ToStr(summary));
    }
  }
  else {
    /* Get summary from original data. */
    if(m_DoAverage)
      summary = float(average(data.begin(),data.end()));
    else {
      summary = median(data.begin(),data.end());
      Verbose::out(2, "Median for whole chip: " + ToStr(m_Summary.size()) + " is: " + ToStr(summary));
    }
  }
  if(summary <= 0) {
    Err::errAbort("Summary value for chip " + 
                     ToStr( m_Summary.size() + 1) +
                     " is not positive: " + ToStr(summary));
  }
  m_Summary.push_back(summary);
  Verbose::out(3, "Summary value for chip " + ToStr(m_Summary.size()) +
               " is: " + ToStr(summary));
  if (!m_TargetUnset) {
    m_Scale.push_back( m_TargetNorm / summary );
  }
}

/** 
 * @brief Signal that no more data is coming (i.e. newChip() will
 * not be called anymore. Currently the class builds up all the
 * summaries and figures out the normalization factores when this
 * function is called.  This makes it difficult to pass through data
 * to downstream.
 */
void MedNormTran::setTarget() {
  unsigned int chipIx = 0;
 
  ///@todo handle each channel independently

  /* Calculate our target. */
  if(m_TargetUnset) {
    if(m_DoAverage)
      m_TargetNorm = average(m_Summary.begin(), m_Summary.end());
    else
      m_TargetNorm = median(m_Summary.begin(), m_Summary.end());
  }
  
  /* Calculate scaling factor for each chip. */
  for(chipIx = 0; chipIx < m_Summary.size(); chipIx++) {
    m_Scale.push_back( m_TargetNorm / m_Summary[chipIx] );
    Verbose::out(2, "Scale for chip: " + ToStr(chipIx) + " is: " + ToStr(m_Scale[m_Scale.size() -1]));
  }
  
}

/** 
 * @brief Signal that no more data is coming (i.e. newChip() will
 * not be called anymore. Currently the class builds up all the
 * summaries and figures out the normalization factores when this
 * function is called.  This makes it difficult to pass through data
 * to downstream.
 */
void MedNormTran::endDataSet() {  
//   if (m_InDiskMart) {
//     if (m_TargetUnset) {
//       int dataSetCount = m_InDiskMart->getCelDataSetCount();
//       for (int d = 0; d < dataSetCount; d++) {
//         std::vector<float> data;
//         data = m_InDiskMart->getCelData(d);
//         transform(d, data);
//         m_TransformedIMart->setProbeIntensity("", data);
//       }
//     }
//   }
  
  /* Flush data in stream. */
  chipStreamEndChips();
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
SelfCreate *MedNormTran::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  float target = 0;
  bool lowPrecision = false;
  bool doAvg = false, calcTarget=false;

  fillInValue(target, "target", param, doc);
  if(target == 0) 
    calcTarget = true;
  fillInValue(doAvg, "doavg", param, doc);
  fillInValue(lowPrecision, "lowprecision", param, doc);

  MedNormTran *mNorm = new MedNormTran(target, doAvg, calcTarget, lowPrecision);
  return mNorm;
}
