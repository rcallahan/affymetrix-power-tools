////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
 * @file   QuantAvgDiff.h
 * @author Pete Klosterman
 * @date   Mon Jun 26 11:27:19 2006
 *
 * @brief Class for doing MAS 4 average difference signal estimation.
 */
#ifndef _QUANTAVGDIFF_H_
#define _QUANTAVGDIFF_H_

//
#include "chipstream/QuantExprMethod.h"
//
#include <cassert>
//

/// Standardized name.
#define QUANT_AVGDIFF_STR "avgdiff"

/**
 * @brief Class for doing MAS 4 average difference signal estimation.
 * This is the average difference between the pm and mm intensity.
 */
class QuantAvgDiff : public QuantExprMethod
{

public:

  /** 
   * What version of the algorithm are we implementing?
   * @return version string.
   */
  std::string getVersion() {
    return "2.0";
  }

  /**
   * @brief Constructor.
   * @param param Map of key/value pairs to set user-defined parameters.
   */
  QuantAvgDiff (std::map<std::string,std::string> &param);

  /**
   * Set the number of probes and number of chips to be used for estimation.
   * @param numProbes Probe count.
   * @param numChips Microarray count.
   */
  void setBounds (unsigned int numProbes, unsigned int numChips);

  /**
   * @brief Set up for estimation given the data about the probe
   * set, chip layout and intensities.
   *
   * @param psGroup Probes to be used.
   * @param layout Chip layout annotation.
   * @param iMart Raw intensities.
   * @param iTrans Transformations to be applied to data before use.
   * @param pmAdjust Method to estimate background or MM probe.
   * @return True if setup sucessful, false otherwise.
   */
  bool setUp (ProbeSetGroup& psGroup, const IntensityMart& iMart,
    std::vector<ChipStream *>& iTrans, PmAdjuster& pmAdjust);

  /**
   * @brief Return number of features (probes).
   * @return Feature count.
   */
  inline unsigned int getNumFeatures() { return m_ProbeCount; }

  /**
   * @brief Return number of targets (experiments or chips).
   * @return Target count.
   */
  inline unsigned int getNumTargets() { return m_ChipCount; }

  /**
   * @brief Does this quantification method estimate feature effects?
   * @return false This method doesn't have feature effects.
   */
  inline bool haveFeatureEffects() { return false;}

  /**
   * @brief Does this quantification method supply residuals?
   * @return false This method doesn't have a residual (probe-level result).
   */
  inline bool haveResiduals() { return false;}

  /**
   * @brief Get the residual.  None is supplied by this method;
   * errAbort() is called.
   * @return Residuals.
   */
  inline double getResidual (unsigned int probeIx, unsigned int chipIx)
  {
    Err::errAbort ("No residual is available in the avgdiff method.");
    return -1.0;
  }

  /**
   * @brief Get the feature effect of a particular probe. None is
   * supplied by this method; errAbort() is called.
   * @param probeIx Index of probe in probe sets.
   * @return Feature (probe) effects.
   */
  inline double getFeatureEffect (unsigned int probeIx)
  {
    Err::errAbort ("No feature effect is available in the avgdiff method.");
    return -1.0;
  }

  /**
   * @brief Get the signal estimate for a particular chip.
   *
   * @param chipIx Index of chip in experiment.
   * @return Signal estimate.
   */
  inline double getSignalEstimate (unsigned int chipIx)
  {
    assert (chipIx < m_ChipCount);
    return m_SetSignal[chipIx];
  }

  /**
   * @brief Alias for getSignalEstimate().
   * @param chipIx Index of chip in experiment.
   * @return Signal estimate.
   */
  inline double getTargetEffect (unsigned int chipIx)
  {
    return getSignalEstimate (chipIx);
  }

  /**
   * Reset data.
   */
  void clear()
  {
    m_Probes.clear();
  }

  /**
   * @brief Set the PM data at probe set probe index and chip index.
   *
   * @param probeIx Index of probe in probe set.
   * @param chipIx Index of chip.
   * @param data Data to be set.
   */
  inline void setPMDataAt (unsigned int probeIx, unsigned int chipIx, double data)
  {
    assert (chipIx < m_ChipCount && probeIx < m_ProbeCount);
    m_PM[chipIx][probeIx] = (float)data;
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
   * @brief Set the MM data at probe set probe index and chip index.
   *
   * @param probeIx Index of probe in probe set.
   * @param chipIx Index of chip.
   * @param data Data to be set.
   */
  inline void setMMDataAt (unsigned int probeIx, unsigned int chipIx, double data)
  {
    assert (chipIx < m_ChipCount && probeIx < m_ProbeCount);
    m_MM[chipIx][probeIx] = (float)data;
  }

  /**
   * @brief Calculate the measurements for the requested collection of
   * probe sets.
   */
  void computeEstimate();

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName (QUANT_AVGDIFF_STR);
    doc.setDocDescription ("Calculates the average measurement for a probeset using the MAS 4 average difference algorithm, namely the average difference between the pm and mm probe signal.");
    // No parameters.
  }

  /**
   * @brief Supply a brief description of this algorithm.
   * @return SelfDoc object.
   */
  static SelfDoc explainSelf()
  {
    SelfDoc doc;
    setupSelfDoc(doc);
    return doc;
  }

  /**
   * @brief Create a new AvgDiff object.
   *
   * @param param Map of key/value pairs to initialize the object.
   * @return Pointer to SelfCreate object.
   */
  static SelfCreate *newObject (std::map<std::string,std::string> &param);

  /**
   * Custom configuration for this QuantMethod
   */
  virtual void setParameters(PsBoard &board) 
  {
    // no-op here to override the error from parent class
  }

private:

  /// Number of chips in experiment.
  unsigned int m_ChipCount;
  /// Number of probes (features) in probe sets.
  unsigned int m_ProbeCount;
  /// Data matrix for PM probes.
  std::vector<std::vector<float> > m_PM;
  /// Data matrix for MM probes.
  std::vector<std::vector<float> > m_MM;
  /// Signal for probeset in different chips.
  std::vector<double> m_SetSignal;

};

#endif /* _QUANTAVGDIFF_H_ */
