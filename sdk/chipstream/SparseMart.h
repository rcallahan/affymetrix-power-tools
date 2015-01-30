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
 * @file   SparseMart.h
 * @author Chuck Sugnet
 * @date   Fri Oct 21 17:34:46 2005
 *
 * @brief  Class for representing cel file intensity data in memory.
 * Can save a lot of the space if not all of a cel file needs to be
 * remembered for downstream analysis.
 */
#ifndef _SPARSEMART_H_
#define _SPARSEMART_H_

//
#include "chipstream/IntensityMart.h"
//
#include "util/Util.h"
//
#include <algorithm>
#include <iostream>
#include <vector>
//

/**
 * SparseMart Class for representing only a subset of cel file intensity data in
 * memory. Can save a lot of the space if not all of a cel file needs
 * to be remembered for downstream analysis.
 * @image html memory-optimization.png "Only store needed probes (red) in memory. Use a translation vector to avoid having to reset probes ids in memory."
 */
class SparseMart : public IntensityMart
{

public:

  /// Marker for when there is no data supplied.
  static const int m_NoData = -1;

  /**
   * @brief Constructor that takes a boolean mask indicating which
   * probes will be needed for downstream analysis.
   * @param probes - Bitmask where 'true' indicates probe will be
   * needed and 'false' indicates it can be discarded.
   */
  SparseMart(std::vector<bool> &probes);

  /**
   * Constructor that assumes all probes will be loaded, not so sparse
   * I guess.
   * @param size - number of probes on array.
   */
  SparseMart(unsigned int size);

  
  /** 
   * @brief Constructor that takes a vector indicating an ordering of probe ids
   * @param order - vector of probe ids.  Postition in vector indicates ordering.
   */
  SparseMart(const std::vector<probeidx_t> &analysisOrder,
             const std::vector<std::string>& celNames, 
             bool storeAllCelIntensities = false);


  /**
   * @brief Clear stored intensities and intensity counts
   */
  void clear();

  /**
   * @brief Method to create empty sparseMart with same parameters.
   */
  SparseMart* copyMetaDataToEmptyMart() const;

  /**
   * @brief Given the probe index and chip index return the intensity
   * data appropriate for that probe in that chip. Depending on the
   * object that intensity may have been modified from the original
   * found on the chip.
   *
   * @param probeIx - Probe Index number.
   * @param chipIx - Chip Index number.
   * @return double - intensity for that position on array.
   */
  float getProbeIntensity(probeid_t probeIx, chipid_t chipIx, unsigned int channelIx = 0) const;

  /** 
   * @brief Method for getting vector of intensities in original order in CEL file.
   * @param dataSetIx - index of CEL intensity dataset in DiskIntensityMart 
   * @return vector of intensities in CEL file order
   */
  std::vector<float> getCelData(int dataSetIx);

  /** 
   * @brief Method for getting vector of intensities in original order in CEL file.
   * @param celIx - index of input CEL file
   * @param channelIx - index of CEL channel that desired dataSet is on
   * @return vector of intensities in CEL file order
   */
  std::vector<float> getCelData(chipid_t celIx, unsigned int channelIx);

  /**
   * Is a particular probe's data in the sparsemart?
   * @param probeIx - index of probe.
   * @return - true if data is contained false otherwise.
   */
  bool isProbeAvailable(probeid_t probeIx) const;

  /**
   * @brief Given a vector of data use it to fill in all of the datapoints
   * that are going to be needed.
   *
   * @param dataName - Name of the vector of data (usually the cel filename).
   * @param data - cel file intensity data.
   */

  void setProbeIntensity(const int dataIdx, const std::vector<float> &data);

  /**
   * @brief Get the names (cel files) for the various data that has been seen.
   * @return Reference to all of the filenames.
   */
  const std::vector<std::string> &getCelFileNames() const;

  /**
   * @brief Get the total number of probes that this mart can supply.
   * @return int - total number of probes.
   */
  size_t getProbeCount() const;

  /**
   * @brief Get the total number of chips that this mart can supply
   * @return int - total number of chips.
   */

  void reserve(int numProbes, int numChips);

  void setProbes(std::vector<bool> &probes);

  /**
   * @brief store a set of data index groupings , and in particular
   * identify and set datasets grouped by CEL file)
   */
  virtual void setChannelMapping(const IdxGroup &idxGroup);

  /**
   * @brief set object that will translate from chip,channel to cache index
   *
   * @param chip_channel_map - First dimension is chipIx, second is
   * channelIx, value is cache index
   */
   void setChipChannelMap(std::vector<std::vector<unsigned int> >& chip_channel_map);

  int getCelDataSetCount() const;

  int getCelFileCount() const;


private:

  /// Number of columns in m_Data
  int m_ProbeSize;
  /// Maximum size of columns in m_Data
  int m_MaxProbeSize;
  /// Number of columns in m_Data
  int m_NumChips;
  /// Maximum size of columns in m_Data
  int m_MaxNumChips;
  /// Input analysis ordering of probe ids, e.g. from ChipLayout
  std::vector<int> m_AnalysisOrder;
  /// number of unique probe ids in m_AnalysisOrder
  int m_UniqueAnalysisOrderSize;
  /// Ordering of original probe intensities in m_Data (might be an
  /// expanded version of m_AnalysisOrder)
  /// IE: Data[order[probe_id]] = m_Data[probe_id]
  std::vector<probeidx_t> m_Order;
  /// Mapping of probe ids from analysis order to original CEL order
  /// IE: Data[probe_id]=m_Data[m_Map[probe_id]];
  std::vector<int> m_Map;
  std::vector<std::vector<float> > m_Data;
  /// Filenames for each chip. 
  std::vector<std::string> m_FileNames;

};

#endif /* _SPARSEMART_H_ */
