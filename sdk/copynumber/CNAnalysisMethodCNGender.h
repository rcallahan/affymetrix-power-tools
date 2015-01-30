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

#ifndef _CNAnalysisMethodCNGender_H_
#define _CNAnalysisMethodCNGender_H_
/**
 * @file CNAnalysisMethodCNGender.h
 *
 * @brief This header contains the CNAnalysisMethodCNGender class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

/**
 * @brief  The CNGender analysis method.
 *
 */
class CNAnalysisMethodCNGender : public CNAnalysisMethod
{
protected:

    double m_dMaleChrXLowerThreshold;
    double m_dMaleChrXUpperThreshold;
    double m_dMaleChrYLowerThreshold;
    double m_dMaleChrYUpperThreshold;

    double m_dFemaleChrXLowerThreshold;
    double m_dFemaleChrXUpperThreshold;
    double m_dFemaleChrYLowerThreshold;
    double m_dFemaleChrYUpperThreshold;

    double m_dMAPDThreshold;

public:
    static std::string getType() {return "cn-gender";}
    static std::string getDescription() {return "Copynumber CNGender";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodCNGender();
    virtual ~CNAnalysisMethodCNGender() {}

    virtual bool isSegmentTypeAnalysis() {return false;}

    virtual AffxString getName() {return getType();}

    virtual void run();
};

#endif


