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
 * @file   IntensityMart.h
 * @author Chuck Sugnet
 * @date   Fri Sep 23 09:46:10 2005
 * 
 * @brief  Base class for ojbects that will dispense data on 
 *         a [feature,chip] level.
 */

#ifndef _INTENSITYMART_H_
#define _INTENSITYMART_H_

//
#include "chipstream/AptTypes.h"
#include "chipstream/IdxGroup.h"
//
#include "file/CELFileData.h"
#include "portability/affy-base-types.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 *  IntensityMart
 * @brief Base class for all objects that will dispense data on a 
 *        [feature,chip] level.
 */
class IntensityMart {

public:

    /** 
     * @brief Virtual destructor for a virtual class.
     */
    virtual ~IntensityMart() {}

    /**
     * @brief Method to create empty sparseMart with same parameters.
     */
    virtual IntensityMart* copyMetaDataToEmptyMart() const = 0;

    /** 
     * @brief Given the probe index and chip index return the intensity
     * data appropriate for that probe in that chip. Depending on the
     * object that intensity may have been modified from the original
     * found on the chip.
     * 
     * @param probeIx - Probe Index number.
     * @param chipIx - Chip Index number.
     * 
     * @return double - intensity for that position on array.
     */
    virtual float getProbeIntensity(probeid_t probeIx, chipid_t chipIx, unsigned int channelIx = 0) const = 0;

    /** 
     * @brief Return all of the intensities for given chipIx
     * @param chipIx - Chip Index number.
     * @return double - vector of all intensities in the IntensityMart for given chipIx.
     */
    virtual std::vector<float> getCelData(int dataSetIx) = 0;

    /** 
     * @brief Method for getting vector of intensities in original order in CEL file.
     * @param celIx - index of input CEL file
     * @param channelIx - index of CEL channel that desired dataSet is on
     * @return vector of intensities in CEL file order
     */
    virtual std::vector<float> getCelData(chipid_t celIx, 
                                          unsigned int channelIx) = 0;

    /** 
     * @brief Given a vector of data use it to fill in all of the datapoints
     * that are going to be needed.
     * 
     * @param dataName - Name of the vector of data (usually the cel filename).
     * @param data - cel file intensity data. 
     */
    virtual void setProbeIntensity(const int dataIdx, const std::vector<float> &data) = 0;

    /** 
     * @brief Get the names (cel files) for the various data that has been seen.
     * @return Reference to all of the filenames.
     */
    virtual const std::vector<std::string> &getCelFileNames() const = 0;

    /** 
     * @brief Get the total number of probes that this mart can supply.
     * @return int - total number of probes.
     */
    virtual size_t getProbeCount() const = 0;

    /** 
     * @brief Get the total number of CEL channels that this mart can supply
     * @return int - total number of chips.
     */
    virtual int getCelDataSetCount() const = 0;

    /** 
     * @brief Get the total number of CELs/multi-CELs that this mart can supply
     * @return int - total number of chips.
     */
    virtual int getCelFileCount() const = 0;


    /**
     * @brief store a set of data index groupings (e.g. datasets grouped
     * by CEL file) 
     */
    virtual void setChannelMapping(const IdxGroup &idxGroup);

    virtual int getChannelCount() const;

    /** 
     * @brief Method for setting boolean to indicate if diskMart should store all given intesities or only the ones specified in m_Order.
     * @param flag - boolean indicating desired behavior
     */
    void setStoreAllCelIntensities(bool flag);

protected:
  
    /// Object to store information about how the datasets are
    /// grouped, e.g. the "channel" group indicates which datasets
    /// belong to the same CEL file.
    IdxGroup m_CelChannels;

    /// table to lookup translation of chipIx,cel_channel to dataset index
    std::vector<std::vector<unsigned int> > m_ChipChannelToCacheMap;

    /// flag to indicate that all CEL intensities should be stored, not
    /// just the ones listed in the ChipLayout object
    mutable bool m_StoreAllCelIntensities;
};

#endif /* _INTENSITYMART_H_ */
