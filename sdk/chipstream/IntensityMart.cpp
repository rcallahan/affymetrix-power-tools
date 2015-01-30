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
 * @file   IntensityMart.cpp
 * @author Ray Wheeler
 * @date   Tue Nov 10 16:30:47 PST 2009
 * 
 * @brief  Base class for objects that will dispense data on 
 *         a [feature,chip] level.
 */


#include "chipstream/IntensityMart.h"
#include "chipstream/IdxGroup.h"


/**
 * @brief store a set of data index groupings (e.g. datasets grouped
 * by CEL file) 
 */
void IntensityMart::setChannelMapping(const IdxGroup &idxGroup) {
    m_CelChannels = idxGroup;
}

/**
 * @brief get number of CEL channels according to m_ChipChannelToCacheMap
 */
/// Not sure if IntensityMart is the best place to get this data,
/// nor if this is the best way to store it.
int IntensityMart::getChannelCount() const {
    if (m_ChipChannelToCacheMap.empty()) {
        return(1); // 1 channel
    }
    else {
        // number of channels assumed same for all chips
        return(m_ChipChannelToCacheMap[0].size()); 
    }
}

/** 
 * @brief Method for setting boolean to indicate if diskMart should store all given intesities or only the ones specified in m_Order.
 * @param flag - boolean indicating desired behavior
 */
void IntensityMart::setStoreAllCelIntensities(bool flag) {
    m_StoreAllCelIntensities = flag;
}
