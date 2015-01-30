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

#ifndef _CNAnalysisMethodCancer_H_
#define _CNAnalysisMethodCancer_H_
/**
 * @file CNAnalysisMethodCancer.h
 *
 * @brief This header contains the CNAnalysisMethodCancer class definition.
 */

#include "copynumber/CNAnalysisMethodChipstream.h"

/**
 * @brief  The Cancer Analysis Method.
 *
 */
class CNAnalysisMethodCancer : public CNAnalysisMethodChipstream
{
public:
    static std::string getType() {return "cancer";}
    static std::string getDescription() {return "Copynumber Cancer";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodCancer();
    virtual ~CNAnalysisMethodCancer() {}

    virtual bool isSegmentTypeAnalysis() {return true;}

    virtual AffxString getName() {return getType();}

protected:
    virtual void determineSketchProbeSets();
};

#endif


