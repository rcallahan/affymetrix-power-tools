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
 * @file   SketchQuantNormTran.h
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:17:07 2005
 * 
 * @brief Class for doing normalization. Can do sketch and full quantile (just
 * set sketch to chip size) and supports bioconductor compatibility.
 */
#ifndef _SKETCHQUANTNORMTRAN_H_
#define _SKETCHQUANTNORMTRAN_H_

//
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityTransformer.h"
#include "chipstream/PsBoard.h"
//
#include "file5/File5.h"
#include "normalization/normalization.h"
#include "portability/affy-base-types.h"
#include "util/Convert.h"
#include "util/Options.h"
#include "chipstream/DataStore.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <algorithm>
#include <vector>
//

/// String describing quantile norm algorithm
#define SKETCHQUANTNORMSTR "quant-norm"

/**
 * SketchQuantNormTran for doing normalization. Can do sketch and full
 * quantile (just set sketch to chip size) and supports bioconductor
 * compatibility. Unfortunatley this is kind of a bad stream currently as it
 * doesn't really flow any data through...
 */
class SketchQuantNormTran : public ChipStream {

public:

  /** 
   * Constructor.
   *
   * @param sketchSize - Number of data points to use for
   * approximating full quantile normalization.
   * @param biocCompat - Should ties in rank be done in same manner as
   * bioconductor affy package?
   * @param lowPrecision - Should we truncate values as if we had been written to cel file?
   * @param usePmSubset - Should we use just the perfect match probes for normalization?
   * @param target - what should we scale the target sketch to
   * @param doavg - should we use the mean (rather than median)
   */
  SketchQuantNormTran(int sketchSize, bool biocCompat, 
                      bool lowPrecision, bool usePmSubset,
                      float target, bool doavg);

  /**
   * Destructor. 
   */
  ~SketchQuantNormTran();

  /** 
   * @brief Only use a subset of all probes for normalization. For
   * example RMA only uses PM probes or might want to just use control
   * genes.
   * @param subsetProbes - bitmask indicating which probes to use for
   * normalization.
   */
  void setSubProbes(const vector<bool> &subsetProbes);

  /** 
   * @brief Add a stream to the list of those that will receive
   * downstream data.
   * @param stream - ChipStream that wants to be fed our modified
   * data.
   */
  virtual void registerStream(ChipStream *stream) {
    m_Streams.push_back(stream);
  }
  
  /** 
   * @brief register our parent stream that will be passing
   * data to this object.
   * @param stream - who are we getting data from?
   */
  virtual void registerParent(ChipStream *stream) {
    m_ParentStream = stream;
  }

  /** 
   * @brief Get a reference to our parent stream.
   * @return Reference to parent stream.
   */
  virtual const ChipStream &getParent() {
    return *m_ParentStream;
  }

  void setProbeCount(int probeCount) { 
    if(m_SketchSize == -1) {
      m_SketchSize = Max((int)(.01 * probeCount), 50000);
      m_SketchSize = Min(50000, probeCount);
      Verbose::out(2, "Setting sketch size to: " + ToStr(m_SketchSize));
    }
    else if(m_SketchSize == 0) {
      m_SketchSize = probeCount;
      Verbose::out(2, "Setting sketch size to entire chip: " + ToStr(m_SketchSize));
    }
  }

  void transformDataSuppliedSketch(int index, std::vector<float> &data);

  void transformDataSuppliedSketch(std::vector<float> &data);

  void newChip(std::vector<float> &data);

  void finishedChips();

  /** 
   * @brief Method for being passed a new cel file worth of data.
   * @param data - vector of cel file data.
   */
  void newDataSet(IntensityMart* iMart);
   
  /** 
   * @brief Signal that no more data is coming (i.e. newChip() will not be
   * called anymore. Currently the class builds up all the sketches and then
   * does the normalization when this function is called. This makes it
   * difficult to pass through data to downstream. If a precomputed sketch
   * is supplied then the data can be passed through chip by chip.
   */
  void endDataSet();

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions();

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setUpSelfDoc(SelfDoc &doc);

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf() {
    SelfDoc doc;
    setUpSelfDoc(doc);
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
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

  /** 
   * @brief Is the portion of the chip comprising pm probes being used as a subset?
   * @return true if pm subset is being used false otherwise.
   */
  inline bool getUsePmSubset() { return m_UsePmSubset; }

  /** 
   * @brief Set a filename to write out the quantile distrubtion to. This
   * can then be used to normalize new cel files to an old batch without
   * rerunning the entire batch. 
   * 
   * @param fileName - name of file to write distribution sampled to.
   */
  void saveTargetSketch(const std::string& fileName) {
    m_FileOutName = fileName;
  }
  void saveTargetSketch_a5(const std::string& fileName) {
    m_a5_filename = fileName;
  }
  void saveTargetSketch_group_a5(affx::File5_Group* group5) {
    m_a5_shared_group = group5;
  }

  /** 
   * @brief Set a target distribution to normalize all the chips to.
   * Note that this routine will give slightly different values than the 
   * matching routine 
   * @param targetSketch 
   */
  void setTargetSketch(const std::vector<float> &targetSketch) {
    if(targetSketch.size() != m_SketchSize) {
      Err::errAbort("SketchQuantNormTran::setTargetSketch() - " + 
                      ToStr("target sketch (N=") + ToStr(targetSketch.size()) + 
                      ") must be same size as sketch (N=" + ToStr(m_SketchSize));
    }
    m_AverageSketch = targetSketch;
    if(!m_BiocCompat) {
      m_PartialSums.resize(m_AverageSketch.end() - m_AverageSketch.begin() + 1);
      m_PartialSums[0] = 0;
      std::transform(m_AverageSketch.begin(), m_AverageSketch.end(), m_PartialSums.begin() + 1, adder<double>());
    }
    m_UsePrecompSketch = true;
    scaleTargetSketch();
  }

  /** 
   * How big is the sketch being used? If zero then entire dataset.
   * @return - Size of the subset of probes being used for normalization.
   */
  unsigned int getSketchSize() { return m_SketchSize; }

  virtual void setParameters(PsBoard &board);

  /**
   * Open and read a target distribution from a text file. Must have a column
   * called 'intensities' with distribution quantiles sorted from highest to
   * lowest.
   *
   * @param fileName - text file containing column of quantiles.
   * @param targetSketch - vector to load up with our intensities
   */
  static void readTargetSketchFromFile(const std::string& fileName, std::vector<float> &targetSketch);

//protected:
  
  /** 
   * @brief transform the intensity point supplied coming from a particular
   * probe in a particular microarray.
   * 
   * @param probeIx - Probe index from the cel file.
   * @param chipIx - Set of chips from same sample.
   * @param intensity - Original intensities.
   * @param return - transformed intensities.
   */
  // transform a vector of intensities
  void transform(int chipIx, std::vector<float>& intensity);
  // transform a single intensity
  float transform(int probeIx, int chipIx, float intensity);

  /**
   * update the partial sums based on m_AverageSketch
   */
  void updatePartialSums() {
    m_PartialSums.clear();
    if(!m_BiocCompat) {
        m_PartialSums.resize(m_AverageSketch.end() - m_AverageSketch.begin() + 1);
        m_PartialSums[0] = 0;
        std::transform(m_AverageSketch.begin(), m_AverageSketch.end(), m_PartialSums.begin() + 1, adder<double>());
    }
  }

  /**
   * @brief Scale the target sketch
   */
  void scaleTargetSketch();
  
  /** 
   * Save the current m_AverageSketch to a file.
   * @param fileName - where to write the sketch.
   */  
  void saveTargetSketchToFile(const std::string& fileName);
  void saveTargetSketchToFile_a5(const std::string& fileName);
  void saveTargetSketchToGroup_a5(affx::File5_Group* group5);

  /** 
   * Obtain a sketch for a given vector of data.
   * @param data - data to sample sketch from.
   */
  void extractSketch(std::vector<float> &data);

  /** 
   * End reading chips when we have to compute the target
   * (m_AverageSketch) ourselves.
   */
  void computeTargetSketch();

  /** 
   * When we have a precomputed distribution to normalize against
   * we can handle things differently, like passing data through rather
   * than caching it for computing target sketch.
   * @param data - vector of vectors of new chip data from same sample.
   */
  void newChipSuppliedTargetSketch(std::vector<float> &data);
  /// Function object for extracting sketches.
  ExtractSketch< vector<float>::iterator> *m_ExtractSketch;
  /// Vector representing our target (or average) distribution.
  std::vector<float> m_AverageSketch;
  /// Vector with the partial sums of the average sketch for handling ties correctly
  std::vector<double> m_PartialSums;
  /// Collection of sub-vectors that are our sketches.
  std::vector< std::vector<float>::iterator > m_Sketches;
  /// Memory that we have to free when we are done.
  std::vector< std::vector<float> * > m_ToFree;
  /// Bitmask telling us which probes to use for normalization.
  std::vector<bool> m_SubProbes;
  /// How big of a sketch are we using?
  int m_SketchSize;
  /// How many probes are we actually using for normalization.
  int m_SubProbeCount;
  /// Should we resolve ties same way as bioconductor affy package?
  bool m_BiocCompat;
  /// Should we mimic the truncation seen in normalized cel files?
  bool m_LowPrecision;
  /// Should we only use subset of probes?
  bool m_UsePmSubset;
  /// Should we use a pre supplied sketch?
  bool m_UsePrecompSketch;
  /// target to scale the sketch to
  float m_Target;
  /// scale to the average (versus median)
  bool m_DoAverage;
  /// Where to save precomputed sketch.
  std::string m_FileOutName;
  ///
  std::string m_a5_filename;
  ///
  affx::File5_Group* m_a5_shared_group;

  /// Md5sum of probe ids used.
  std::string m_ProbeMd5sum;
};


#endif /* _SKETCHQUANTNORMTRAN_H_ */
