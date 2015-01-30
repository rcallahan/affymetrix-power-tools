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
 * @file   PmSumAdjust.cpp
 * @author Chuck Sugnet
 * @date   Mon May 14 11:45:08 2007
 *
 * @brief Class for doing an additive adjustment based on a probe that
 * may be hybridizing to another allele
 *
 */


//
#include "chipstream/PmSumAdjust.h"
//
#include "chipstream/ChipLayout.h"

/** Constructor. */
PmSumAdjust::PmSumAdjust(bool noMatchOk)
{
  setDocName(PMSUMADJUSTSTR);
  setDocDescription("Add the PM intensities from two alleles.");
  m_NoMatchOk = noMatchOk;
  m_Type = getDocName();
  setDocOptions(getDefaultDocOptions());
}

/**
 * @brief Constructor that takes the list of probe sets to remember
 * the pm/pm pairs.
 *
 * @param layout - Annotation of probes on microarray.
 */
void PmSumAdjust::setLayout(ChipLayout &layout)
{
  m_Vec = layout.getPmAlleleMatchVec();
}

/**
 * @brief Constructor that takes the list of probe sets to remember
 * the pm/pm pairs.
 *
 * @param board - blackboard with various state to access
 */
void PmSumAdjust::setParams(PsBoard &board)
{
  DataStore *info = board.getProbeInfo();
  info->getProbePmAlleleMatch(m_Vec);
}


/**
 * Subset of probes that should be loaded into memory representation.
 * @param probes - bitmask of probes to be used.
 */
void PmSumAdjust::setProbes(std::vector<bool> &probes)
{
  // do nothing, probes should be loaded via probeset
}

/**
 * @brief Given a PM probe intensity add the matching PM intensity
 * for the other allele to the pmIntensity.
 *
 * @param probeIx - Index of probe in cel file data.
 * @param chipIx - Microarray or chip index.
 * @param iMart - IntensityMart which holds raw data.
 * @param iTrans - Vector of transformations that should be performed on raw data.
 * @param pmIintensity - Intensity of perfect match probe to be adjusted, may
 * be modified from original value depending on adjuster
 * @param bgrdAdjust - Background adjustment, if any, recommended (i.e. MM intensity)
 */
void PmSumAdjust::pmAdjustment(int probeIx, int chipIx,
                               const IntensityMart &iMart, std::vector<ChipStream *> &iTrans,
                               float &pmIntensity, float &bgrdAdjust)
{
  int otherPM = 0;
  float otherPmIntensity = 0;
  if (m_Vec.empty())
    Err::errAbort("PmSumAdjust::pmAdjustment() - Appears that chip layout has not been set, no matching probes.");
  otherPM = m_Vec[probeIx];
  if (otherPM == -1 && !m_NoMatchOk)
    Err::errAbort("No  probe for probe with id: " + Convert::toString(probeIx + 1));
  else if (otherPM != -1)
    otherPmIntensity = QuantMethod::transformPrimaryData(otherPM, chipIx, iMart, iTrans);
  pmIntensity += otherPmIntensity;
  bgrdAdjust = 0;
}

/**
 * Fill in the information for Self documentation.
 * @param doc - Self documenter to be filled in.
 */
void PmSumAdjust::setupSelfDoc(SelfDoc &doc)
{
  doc.setDocName(PMSUMADJUSTSTR);
  doc.setDocDescription("Add itensity of PM probe for other allele to PM probes.");
  // No parameters...
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc PmSumAdjust::explainSelf()
{
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
SelfCreate * PmSumAdjust::newObject(std::map<std::string, std::string> &param)
{
  if (param.size() != 0) {
    Err::errAbort(ToStr("No parameters for ") + ToStr(PMSUMADJUSTSTR));
  }
  return new PmSumAdjust(true);
}

