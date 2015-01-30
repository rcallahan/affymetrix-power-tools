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
 * @file CNAnalysisMethodCancer.cpp
 *
 * @brief This file contains the CNAnalysisMethodCancer class members.
 */
#include "copynumber/CNAnalysisMethodCancer.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/AdapterTypeNormTran.h"
#include "chipstream/ChipLayout.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "normalization/normalization.h"
#include "plier/affyplier.h"
#include "util/AffxStatistics.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodCancer::explainSelf()
{
    CNAnalysisMethodCancer obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodCancer::getDefaultDocOptions()
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
SelfCreate* CNAnalysisMethodCancer::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodCancer* pMethod = new CNAnalysisMethodCancer();
    std::string strPrefix = getPrefix();

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodCancer::CNAnalysisMethodCancer() : CNAnalysisMethodChipstream()
{
}

void CNAnalysisMethodCancer::determineSketchProbeSets()
{
    Verbose::out(1, "CNAnalysisMethodCancer::determineSketchProbeSets()");
    for (int iIndex = 0; (iIndex < getProbeSets()->getCount()); iIndex++)
    {
        CNProbeSet* p = getProbeSets()->getAt(iIndex);
        p->setUseForSketch(true);
    }
}
