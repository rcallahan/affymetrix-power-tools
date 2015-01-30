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

#ifndef _CNAnalysisMethodCNNeutralLOH_H_
#define _CNAnalysisMethodCNNeutralLOH_H_
/**
 * @file CNAnalysisMethodCNNeutralLOH.h
 *
 * @brief This header contains the CNAnalysisMethodCNNeutralLOH class definition.
 */

#include "copynumber/CNAnalysisMethodSegment.h"

/**
 * @brief  The CNNeutralLOH Analysis Method.
 *
 */
class CNAnalysisMethodCNNeutralLOH : public CNAnalysisMethodSegment
{
private:
    CNProbeSetArray m_vLocalProbeSets;;
    bool m_bLocalProbeSetsDetermined;

public:
    static std::string getType() {return "cn-neutral-loh";}
    static std::string getDescription() {return "Copynumber CNNeutralLOH";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodCNNeutralLOH();
    virtual ~CNAnalysisMethodCNNeutralLOH() {}

    virtual bool isSegmentTypeAnalysis() {return true;}

    virtual AffxString getName() {return getType();}

    virtual void run();
    void determineLocalProbeSets();
    CNProbeSetArray* getProbeSets() {return &m_vLocalProbeSets;}
};

#endif


