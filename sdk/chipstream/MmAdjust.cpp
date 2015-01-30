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
 * @file   MmAdjust.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 09:32:49 2005
 * 
 * @brief  Class for supplying the mismatch intensity as pm adjustment.
 */

//
#include "chipstream/MmAdjust.h"
//
#include "chipstream/QuantMethod.h"
//
#include "util/Err.h"
//

using namespace std; 

/** 
 * @brief Setup ourselves based on data in the blackboard, specifically the
 * the pm/mm pairs.
 * 
 * @param board - blackboard with various state
 */
void MmAdjust::setParameters(PsBoard &board) {
  DataStore *info = board.getProbeInfo();
  info->getProbePmMm(m_Vec);
}

/** 
 * @brief Constructor that takes the list of probe sets to remember
 * the pm/mm pairs.
 * 
 * @param layout - Annotation of probes on microarray.
 */
void MmAdjust::setLayout(ChipLayout &layout) {
  m_Vec.resize(layout.getProbeCount());
  fill(m_Vec.begin(), m_Vec.end(), -1);
  m_Vec = layout.getPmMmVec();
}

/** 
 * Subset of probes that should be loaded into memory representation.
 * @param probes - bitmask of probes to be used.
 */
void MmAdjust::setProbes(std::vector<bool> &probes) {
  // do nothing, probes should be loaded via probeset
}

/** 
 * @brief Given a PM probe intensity, supply the MM probe intensity.
 * @param probeIx - Index of probe in cel file data.
 * @param chipIx - Microarray or chip index.
 * @param iMart - IntensityMart which holds raw data.
 * @param iTrans - Vector of transformations that should be performed on raw data.
 * @param pmIintensity - Intensity of perfect match probe to be adjusted, may
 * be modified from original value depending on adjuster
 * @param bgrdAdjust - Background adjustment, if any, recommended (i.e. MM intensity)
 */
void MmAdjust::pmAdjustment(int probeIx, int chipIx, 
                            const IntensityMart &iMart, std::vector<ChipStream *> &iTrans, 
                            float &pmIntensity, float &bgrdAdjust) {
  int mmIx = 0;
  std::map<int,int>::iterator iter;
  if(m_Vec.empty()) 
    Err::errAbort("MmAdjust::pmAdjustment() - Appears that chip layout has not been set, no mismatch probes.");
  mmIx = m_Vec[probeIx];
  if(mmIx == -1) 
    Err::errAbort("No MM probe for probe with id: " + Convert::toString(probeIx + 1));
  bgrdAdjust = QuantMethod::transformPrimaryData(mmIx, chipIx, iMart, iTrans);
}
