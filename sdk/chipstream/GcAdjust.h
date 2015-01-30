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
 * @file   GcAdjust.h
 * @author Chuck Sugnet
 * @date   Mon Oct 24 09:42:02 2005
 * 
 * @brief Class for doing an adjustment based on the median intensity of probes
 * with similar GC content.
 */
#ifndef _GCADJUST_H_
#define _GCADJUST_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/PmAdjuster.h"
#include "chipstream/QuantMethod.h"
//
#include "stats/stats.h"
#include "util/Verbose.h"
//
#include <map>
#include <vector>
//

/// String describing gc adjust algorithm/module
#define GCADJUSTSTR "pm-gcbg"

/**
 *  Class for doing an adjustment based on the median intensity of probes
 * with similar GC content.
 */
class GcAdjust : public PmAdjuster {

public:

  /** Constructor. */
  GcAdjust(int maxGc=26) : m_MaxGc(maxGc), m_BinsFilled(false) {
    setDocName(GCADJUSTSTR);
    setDocDescription("Do an adjustment based on the median intensity of probes with similar GC content.");
    m_Type = getDocName();
    setDocOptions(getDefaultDocOptions());
  }

  /** 
   * @brief Constructor taking information about both the chip and
   * probes to be used for estimating parameters.
   * 
   * @param layout - Annotation of probes on microarray.
   * @param controlProbes - Probes to be used for estimating parameters.
   */
  void setLayout(ChipLayout &layout, vector<Probe *> &controlProbes);

  /** 
   * @brief Constructor taking information about both the chip and
   * probes to be used for estimating parameters.
   * 
   * @param board - Blackboard containing lots of state.
   */
  void setParameters(PsBoard &board);
  
  /** 
   * @brief Set the probes necessary for estimating GC background to
   * true so they will be loaded.
   * 
   * @param probes - bitmask for probe to be loaded.
   */
  void setProbes(vector<bool> &probes);

  /** 
   * How many probes are in the background list?
   * @return - number of probes being used for background.
   */
  uint32_t getBgProbeCount() {
    return (uint32_t) m_Probes.size();
  }

  /** 
   * @brief Fill in the GC bins using the raw data in intensity mart
   * and processed by iTrans ChipStream.
   * 
   * @param chipIx - Index of chip to be filled in.
   * @param iMart - Raw data to be used.
   * @param iTrans - Chipstreams that will be used to transform the raw data.
   */
  void fillChipBins(int chipIx, const IntensityMart &iMart, std::vector<ChipStream *> &iTrans);

  /** 
   * @brief Calculate the parameters for each GC count.
   * 
   * @param iMart - Raw data from chips.
   * @param iTrans - Transformations to apply to data.
   */
  void calcParams(const IntensityMart &iMart, std::vector<ChipStream *> &iTrans);


  /** 
   * @brief For probe supplied look up the median of background probes with
   * same GC nucleotide count.
   * @param probeIx - Index of probe in cel file data.
   * @param chipIx - Microarray or chip index.
   * @param iMart - IntensityMart which holds raw data.
   * @param iTrans - Vector of transformations that should be performed on raw data.
   * @param pmIintensity - Intensity of perfect match probe to be adjusted, may
   * be modified from original value depending on adjuster
   * @param bgrdAdjust - Background adjustment, if any, recommended (i.e. MM intensity)
   */
  void pmAdjustment(int probeIx, int chipIx, 
                    const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                    float &pmIntensity, float &bgrdAdjust);

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> opts;
    SelfDoc::Opt probeMd5 = {"subsetmd5", SelfDoc::Opt::String, "", "", "NA", "NA",
                             "Md5sum of the probe ids being used as background."};
    opts.push_back(probeMd5);
    return opts;
  }

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf() { 
    GcAdjust doc;

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
  static SelfCreate *newObject(std::map<std::string,std::string> &param) {
    if(param.size() != 0) { Err::errAbort(ToStr("No parameters for ") + ToStr(GCADJUSTSTR)); }
    return new GcAdjust();
  }
  
private:

  /// Maximum GC count probe seen.
  int m_MaxGc;
  /// Have we calculated our parameters yet?
  bool m_BinsFilled;
  /// GC count of each probe as indexed by position.
  std::vector<char> m_ProbeGcVec;
  /// Parameters learned (median of GC background probes).
  std::vector<vector<float> > m_Bins;
  /// Probe to be used for estimating background median
  std::vector<int> m_Probes;
  /// mdsum of probe ids
  std::string m_ProbeMd5Sum;
};

#endif /* _GCADJUST_H_ */
