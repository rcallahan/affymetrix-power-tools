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
 * @file   GcBg.h
 * @author Chuck Sugnet
 * @date   Mon Sep 25 16:26:13 PDT 2006
 * 
 * @brief Class for doing an background adjustment on all probes based
 * on the median intensity of probes with similar GC content.
 */
#ifndef _GCBG_H_
#define _GCBG_H_

#include "chipstream/PsBoard.h"
//
#include "chipstream/ChipStream.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantMethod.h"
//
#include "stats/stats.h"
#include "util/Verbose.h"
//
#include <map>
#include <vector>
//

/// String describing gc adjust algorithm/module
#define GCBG "gc-bg"

/**
 *  Class for doing an adjustment based on the median intensity of probes
 * with similar GC content.
 */
class GcBg : public ChipStream {

public:

  /** Constructor. */
  GcBg(int maxGc=26, bool attenuate=true, float l=0.005, float h=-1) {
    m_MaxGc = maxGc;
    m_Type = GCBG;
    m_ChipCount = 0;
    setupSelfDoc(*this);
    m_Attenuate = attenuate;
    m_L = l;
    m_H = h;
    setOptValue("attenuate", m_Attenuate);
    setOptValue("l", ToStr(m_L));
    setOptValue("h", ToStr(m_H));
  }

  /** 
   * Caluculates medians for different levels (bins) of GC count
   * between 0 and binCount in size and appends to the chipBins. The
   * probes vector contains a list of probes that are thought to be
   * representative of background at different GC counts. The
   * probeGcVec contains the GC count for all probes on the array.
   * 
   * @param data - One microarray's worth of data.
   * @param binCount - Number of bins we are calculating GC content for.
   * @param chipBins - Matrix of chip by gc count containing the
   * median intensity for a collection of gc background
   * probes. Ordering is chipBins[chip][gcCount] and parameters from
   * latest chip in data vector will be appended as last row.
   * @param probes - Probes to use for estimating background for particular GC bings.
   * @param probeGcVec - GC count of every probe on microarray indexed by id.
   * @param warnings - Print warnings if a particular bin is empty?
   * For example: print a warning if there are no probes in the GC
   * count zero bin.
   */
  static void learnParameters(const std::vector<float> &data, 
                              int binCount, 
                              std::vector<vector<float> > &chipBins,
                              std::vector<int> &probes, 
                              std::vector<char> &probeGcVec, 
                              bool warnings);

  /** 
   * Give the background adjusted intensity for a particular probe on a particular chip.
   * @param probeIx - Probe of interest.
   * @param chipIx - Set of Chips of interest.
   * @param intensity - Current set of intensities.
   * @param return - background adjusted set of intensities.
   */
  float transform(int probeIx, int chipIx, float intensity, PsBoard &board) {
    return transform(probeIx, chipIx, intensity);
  }

  float transform(int probeIx, int chipIx, float intensity);

  float transform(int probeIx, int chipIx, float intensity, int probeGc);

  void transform(int chipIx, std::vector<float>& intensity);

  void newChip(std::vector<float> &data);

  /** 
   * @brief Method for being passed a new cel file worth of data.
   * @param data - vector of vectors of cel file data from same sample.
   */
  void newDataSet(IntensityMart* iMart);

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName(GCBG);
    doc.setDocDescription("Subtract bacground based on median intensity of probes with similar GC content.");

    doc.setDocOptions(getDefaultDocOptions());
  }

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> opts;

    SelfDoc::Opt attenuate = {"attenuate", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
                              "Indicate whether or not to attenuate gcbg value when subtracting background. If not set then downstream algorithms must be able to handle negative values."};
    opts.push_back(attenuate);
    SelfDoc::Opt l = {"l", SelfDoc::Opt::Double, "0.005", "0.005", "0", "1",
                      "Tunable parameter for attenuating gc background."};
    opts.push_back(l);
    SelfDoc::Opt h = {"h", SelfDoc::Opt::Double, "-1", "-1", "-1", "NA",
                      "Used fixed constant to attenuate gc bacground. Default results in using 4*PM*GCBG*L."};
    opts.push_back(h);
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
    SelfDoc doc;
    setupSelfDoc(doc);
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
    SelfDoc doc = explainSelf();
    bool attenuate = true;
    float l = 0.005;
    float h = -1;
    fillInValue(attenuate, "attenuate", param, doc);
    fillInValue(l, "l", param, doc);
    fillInValue(h, "h", param, doc);
    GcBg *gcbg = new GcBg(26, attenuate, l, h);
    return gcbg;
  }

  /** 
   * Setting the GC counts for all the probes on the array.
   * @param vec - GC count for every probe on array indexed by id.
   */
  void setProbeGcVec(const std::vector<char> &vec);

  /** 
   * Setting the GC counts for all the probes on the array.
   * @param vec - GC count for every probe on array indexed by id.
   */
  void setProbeGcVec(const std::vector<unsigned char> &vec);

  /** 
   * Set the probes to use for estimating background.
   * @param vec - Vector of probes that are thought to be
   * representative of background binding.
   */
  void setControlProbes(std::vector<Probe *> &vec);

  /** 
   * Set the indexes of the probes to use for estimating background.
   * @param vec - Vector of probe indexes that are thought to be
   * representative of background binding.
   */
  void setControlProbes(std::vector<int> &vec);

  /**
   * Setup the background probes and background gc indexes by reading
   * them from the blackboard
   */
  void setParameters(PsBoard &board);

private:

  /// Parameters learned (median of GC background probes).
  std::vector<vector<float> > m_Bins;
  /// Probe to be used for estimating background median
  std::vector<int> m_BgProbes;
  /// GC count of each probe matching that above
  std::vector<char> m_BgProbeGc;
  /// mdsum of probe ids
  std::string m_ProbeMd5Sum;
  /// Should we use an attenuated mismatch, or simply cap the lower bounds when using non-PM-only adjuster
  bool m_Attenuate;
  /// Tunable parameter for handling background value in non-PM-only adjusters
  float m_L;
  /// Override use of computed H from data with fixed constant in background attenuation
  float m_H;
  /// How many chips have been seen.
  int m_ChipCount;
  /// Maximum GC count allowed.
  int m_MaxGc;
};

#endif /* _GCBG_H_ */
