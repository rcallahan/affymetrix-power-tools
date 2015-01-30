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

#ifndef _CNAnalysisMethodCNCyto2Gender_H_
#define _CNAnalysisMethodCNCyto2Gender_H_
/**
 * @file CNAnalysisMethodCNCyto2Gender.h
 *
 * @brief This header contains the CNAnalysisMethodCNCyto2Gender
 *        class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

/**
 * @brief  The CNGender analysis method.
 *
 */

class CNAnalysisMethodCNCyto2Gender : public CNAnalysisMethod
{
private:
    double m_dCutoff;
        CNProbeSetArray  m_vLocalProbeSets;
        bool m_bLocalProbeSetsDetermined;
public:
    static std::string getType() { return "cn-cyto2-gender"; }
    static std::string getDescription() { return "Copynumber CNCyto2Gender"; }
    static std::string getVersion() { return "1.0"; }

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

        void determineLocalProbeSets();

        CNProbeSetArray* getProbeSets()
        {
                return &m_vLocalProbeSets;
        }


public:
    CNAnalysisMethodCNCyto2Gender();
    virtual ~CNAnalysisMethodCNCyto2Gender() { }

    virtual bool isSegmentTypeAnalysis() { return false; }

    virtual AffxString getName() { return getType(); }

    virtual void run();
};

#endif


