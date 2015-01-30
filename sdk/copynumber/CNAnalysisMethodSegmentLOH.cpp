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
 * @file CNAnalysisMethodSegmentLOH.cpp
 *
 * @brief This file contains the CNAnalysisMethodSegmentLOH class members.
 */
#include "copynumber/CNAnalysisMethodSegmentLOH.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/AffxStatistics.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodSegmentLOH::explainSelf()
{
    CNAnalysisMethodSegmentLOH obj;
    SelfDoc doc;
    doc.setDocName(obj.getType());
    doc.setDocDescription(obj.getDescription());
    doc.setDocOptions(obj.getDefaultDocOptions());
    return doc;
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> CNAnalysisMethodSegmentLOH::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  return opts;
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
SelfCreate* CNAnalysisMethodSegmentLOH::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodSegmentLOH* pMethod = new CNAnalysisMethodSegmentLOH();
    std::string strPrefix = getPrefix();

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodSegmentLOH::CNAnalysisMethodSegmentLOH()
{
	m_bLocalProbeSetsDetermined = false;
}

/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodSegmentLOH::run()
{
    Verbose::out(1, "CNAnalysisMethodSegmentLOH::run(...) start");
    isSetup();
    determineLocalProbeSets();
    m_vSegments.deleteAll();
    newSegments(getSegmentType(), getProbeSets(), m_pEngine->getOptInt("minSegSeparation"));
    calculateSummaryLOH();
    Verbose::out(1, "CNAnalysisMethodSegmentLOH::run(...) end");
}

void CNAnalysisMethodSegmentLOH::determineLocalProbeSets()
{
    if (m_bLocalProbeSetsDetermined) {return;}
    int iNumberOfProbeSets = CNAnalysisMethod::getProbeSets()->getCount();
    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
    for (int iIndex = 0; (iIndex < iNumberOfProbeSets); iIndex++)
    {
		if (isCytoScanHD)
		{
			if (CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAsSNP() )
			{
				getProbeSets()->add( CNAnalysisMethod::getProbeSets()->getAt(iIndex));
			}
		}
		else
		{
			getProbeSets()->add( CNAnalysisMethod::getProbeSets()->getAt(iIndex));
		}
    }
    m_bLocalProbeSetsDetermined = true;
}
