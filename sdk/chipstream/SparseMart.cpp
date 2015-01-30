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
 * @file   SparseMart.cpp
 * @author Chuck Sugnet
 * @date   Fri Oct 21 17:38:28 2005
 *
 * @brief Class for representing cel file intensity data in memory.
 * Can save a lot of the space if not all of a cel file needs to be
 * remembered for downstream analysis.
 */

//
#include "chipstream/SparseMart.h"
//
#include "util/Util.h"

using namespace std;
/**
 * @brief Constructor that takes a boolean mask indicating which
 * probes will be needed for downstream analysis.
 * @param probes - Bitmask where 'true' indicates probe will be
 * needed and 'false' indicates it can be discarded.
 */
SparseMart::SparseMart(std::vector<bool> &probes) {
    unsigned int i = 0;
    unsigned int count = 0;
    m_AnalysisOrder.reserve(probes.size());
    for (i = 0; i < probes.size(); i++) {
        if (probes[i] == true) {
            m_AnalysisOrder.push_back(i);
            count++;
        }
    }
    m_ProbeSize = count;
    m_MaxProbeSize = 0;
    m_NumChips = 0;
    m_MaxNumChips = 0;
    m_UniqueAnalysisOrderSize = 0;
}

/** 
 * @brief Constructor that takes a vector indicating an ordering of probe ids
 * @param order - vector of probe ids.  Postition in vector indicates ordering.
 */
SparseMart::SparseMart(const std::vector<probeidx_t> &analysisOrder, 
                       const std::vector<std::string>& celNames, 
                       bool storeAllCelIntensities) {
    m_AnalysisOrder = analysisOrder;
    m_FileNames = celNames;
    m_StoreAllCelIntensities = storeAllCelIntensities;
    m_MaxProbeSize = 0;
    m_NumChips = 0;
    m_ProbeSize = 0;
    m_MaxNumChips = 0;
    m_UniqueAnalysisOrderSize = 0;
}


/**
 * Constructor that assumes all probes will be loaded, not so sparse
 * I guess.
 * @param size - number of probes on array.
 */
SparseMart::SparseMart(unsigned int size) {
    m_Map.reserve(size);
    for (unsigned int i = 0; i < size; i++) {
        m_Map.push_back(i);
    }
    m_ProbeSize = size;
    m_MaxProbeSize = 0;
    m_NumChips = 0;
    m_MaxNumChips = 0;
    m_UniqueAnalysisOrderSize = 0;
}



/**
 * @brief Method to create empty sparseMart with same parameters.
 */
SparseMart* SparseMart::copyMetaDataToEmptyMart() const {

    SparseMart* newMart = new SparseMart(m_AnalysisOrder, 
                                         m_FileNames, 
                                         m_StoreAllCelIntensities);
    newMart->m_Order = m_Order;
    newMart->m_Map = m_Map;
    newMart->m_UniqueAnalysisOrderSize = m_UniqueAnalysisOrderSize;
    newMart->m_ChipChannelToCacheMap = m_ChipChannelToCacheMap;
    return newMart;
}

/**
 * @brief store a set of data index groupings , and in particular
 * identify and set datasets grouped by CEL file
 */
void SparseMart::setChannelMapping(const IdxGroup &idxGroup) {
    m_CelChannels = idxGroup;
    std::vector<std::vector<unsigned int> > TempChannelGroup = m_CelChannels.getGroupingVec("channels");
    setChipChannelMap(TempChannelGroup);
}

/**
 * @brief set object that will translate from chip,channel to cache index
 * @param chip_channel_map - First dimension is chipIx, second is
 * channelIx, value is cache index
 */
void SparseMart::setChipChannelMap(std::vector<std::vector<unsigned int> >& chipChannelMap) {
    m_ChipChannelToCacheMap = chipChannelMap;
}

/** 
 * @brief Given a vector of data use it to fill in all of the datapoints
 * that are going to be needed.
 *
 * @param dataName - Name of the vector of data (usually the cel filename).
 * @param data - cel file intensity data.
 */
void SparseMart::setProbeIntensity(const int dataIdx, 
                                   const std::vector<float> &data) {
    unsigned int i = 0;
    int probeCount = data.size();
    assert(probeCount > 0);
    if (m_Map.size() == 0) {
        // initialize map and prune AnalysisOrder of duplicate ids
        m_Order.reserve(probeCount);
        m_Map.resize(probeCount);

        fill(m_Map.begin(), m_Map.end(), -1);
        fill(m_Order.begin(), m_Order.end(), -1);

        int indexCount = 0;
        for (i = 0; i < m_AnalysisOrder.size(); i++) {
            int probeIndex = m_AnalysisOrder[i]; 
            if (m_Map[probeIndex] == -1) {
                m_Order.push_back(probeIndex);
                m_Map[probeIndex] = indexCount++;
            }
        }
        // store the size of the analysis order after any duplicate probes
        // ids have been removed.
        m_UniqueAnalysisOrderSize = m_Order.size();
        if (m_StoreAllCelIntensities && probeCount > m_Order.size()) {
            // probe intensities that are not in m_AnalysisOrder are tacked
            // onto the end of m_Order.  in this way, m_AnalysisOrder probe
            // intensities will be cached efficiently, and yet
            // non-m_AnalysisOrder probe intensities will still be available
            // in the diskmart.
            for (i = 0; i < probeCount; i++) {
                if (m_Map[i] == -1) {
                    m_Order.push_back(i);
                    m_Map[i] = indexCount++;
                }
            }
        }
    }

    int write_data_size = data.size();
    if (!m_StoreAllCelIntensities &&
        write_data_size > m_UniqueAnalysisOrderSize) {
        write_data_size = m_UniqueAnalysisOrderSize;
    }

    if (dataIdx >= m_Data.size()) {
//     Verbose::out(2, "SparseMart::setProbeIntensity() - Adding a chip, try reserve to be more memory efficient.");
        // Enlarge m_Data to accomodate new data
        m_Data.resize(dataIdx + 1);
    }
    if (m_Data[dataIdx].empty()) {
        // No data written for this dataIdx yet.  Resize accordingly
        m_Data[dataIdx].resize(write_data_size, 0.0);
        m_NumChips++;
    }

    for (int i = 0; i < write_data_size; i++) {
        if (m_Order[i] >=0) {
//       assert(m_Order[i] < write_data_size);
            m_Data[dataIdx][i] = data[m_Order[i]];
        }
    }
}
/** 
 * @brief Given the probe index and chip index return the intensity
 * data appropriate for that probe in that chip. Depending on the
 * object that intensity may have been modified from the original
 * found on the chip.
 * 
 * @param probeIx - Probe Index number.
 * @param chipIx - Chip Index number.
 * @return float - intensity for that position on array.
 */
float SparseMart::getProbeIntensity(probeid_t probeIx, 
                                    chipid_t chipIx, 
                                    unsigned int channelIx) const {
    assert(probeIx < m_Map.size() && probeIx >= 0);
    assert(m_Map[probeIx] >= 0);
  
    int dataSetIx;
    if (m_ChipChannelToCacheMap.empty()) {
        dataSetIx = chipIx;
    }
    else {
        dataSetIx = m_ChipChannelToCacheMap[chipIx][channelIx];
    }
    assert(dataSetIx < m_Data.size() && dataSetIx >= 0);
    return m_Data[dataSetIx][m_Map[probeIx]];
}


void SparseMart::clear() {
    unsigned int i = 0;
    m_ProbeSize = 0;
    m_NumChips = 0;
    fill(m_Map.begin(), m_Map.end(), -1);
    for (i = 0; i < m_Data.size(); i++) {
        fill(m_Data[i].begin(), m_Data[i].end(), 0.0f);
    }
}

void SparseMart::reserve(int numProbes, int numChips) {
    if (numProbes < m_ProbeSize) {
        Err::errAbort("SparseMart::reserve() - Have to reserve at least enough probes for current size: " + ToStr(m_ProbeSize));
    }
    m_MaxProbeSize = numProbes;
    m_MaxNumChips = numChips;
    m_Data.resize(m_MaxNumChips);
    m_FileNames.resize(m_MaxNumChips);
    Verbose::out(2, 
                 "Reserving with: " + ToStr(numChips) + " chips and " + 
                 ToStr(numProbes) + " probes.");
    for (int i = 0; i < numChips; i++) {
        m_Data[i].resize(m_MaxProbeSize);
        std::fill(m_Data[i].begin(), m_Data[i].end(), 0.0f);
    }
}

void SparseMart::setProbes(std::vector<bool> &probes) {
    unsigned int i = 0;
    unsigned int count = 0;
    m_Map.resize(probes.size());
    /* set up the map. */
    for (i = 0; i < probes.size(); i++) {
        if (probes[i] == true) {
            m_Map[i] = count++;
        } else {
            m_Map[i] = -1;
        }
    }
    m_ProbeSize = count;
    /* If there are more probes now, reallocate our buffers. */
    if (m_ProbeSize > m_MaxProbeSize) {
        Verbose::out(2, 
                     "SparseMart::setProbes() - Current probes: " + 
                     ToStr(m_ProbeSize) + " is greater than allocated: " + 
                     ToStr(m_MaxProbeSize) + " doing reallocation.");
        reserve(m_ProbeSize, m_Data.size());
    }
    /* Always fill with 0's for sanity and debugging. */
    for (i = 0; i < m_Data.size(); i++) {
		float zero = 0.0;
        fill(m_Data[i].begin(), m_Data[i].end(), zero);
    }
}

/**
 * @brief Return all of the intensities for given chipIx
 * @param dataSetIx - Chip Index number.
 * @return celOrderedData - vector of all intensities in the SparseMart for
 * given dataSetIx in CEL order.
 */
std::vector<float> SparseMart::getCelData(int dataSetIx) {
    assert(dataSetIx < m_Data.size() && dataSetIx >= 0);
    std::vector<float> celOrderedData(m_Data[dataSetIx].size());
    for (int i = 0; i < m_Map.size(); i++) {
        celOrderedData[i] = m_Data[dataSetIx][m_Map[i]];
    }
    return celOrderedData;
}


std::vector<float> SparseMart::getCelData(chipid_t celIx, unsigned int channelIx) {
  int dataSetIx = 0;
  // convert celname index to cache index, if a mapping exists
  if (m_ChipChannelToCacheMap.empty()) {
    if (channelIx == 0) {
      dataSetIx = celIx;
    }
    else {
      Err::errAbort("\nSparseMart::getCelData - Accessing CEL data with multi-channel index when no multi-channel data is available.");
    }
  }
  else {
    assert(celIx >= 0 && celIx < m_ChipChannelToCacheMap.size());
    assert(channelIx >= 0 && channelIx < m_ChipChannelToCacheMap[celIx].size());
    dataSetIx = m_ChipChannelToCacheMap[celIx][channelIx];
  }
  std::vector<float> data = getCelData(dataSetIx);
  return data;
}


/**
 * Is a particular probe's data in the sparsemart?
 * @param probeIx - index of probe.
 * @return - true if data is contained false otherwise.
 */
bool SparseMart::isProbeAvailable(probeid_t probeIx) const {
    assert(probeIx < m_Map.size() && probeIx >= 0);
    return m_Map[probeIx] >= 0;
}


/**
 * @brief Get the names (cel files) for the various data that has been seen.
 * @return Reference to all of the filenames.
 */
const std::vector<std::string> &SparseMart::getCelFileNames() const {
    return m_FileNames;
}

/**
 * @brief Get the total number of probes that this mart can supply.
 * @return int - total number of probes.
 */
size_t SparseMart::getProbeCount() const {
    return m_Map.size();
}

/**
 * @brief Get the total number of chips that this mart can supply
 * @return int - total number of chips.
 */
int SparseMart::getCelDataSetCount() const {
    return m_NumChips;
}

int SparseMart::getCelFileCount() const {
    return m_FileNames.size();
}
