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

#ifndef _CNFamilialAnalysisMethodHemizygousParentOfOrigin_H_
#define _CNFamilialAnalysisMethodHemizygousParentOfOrigin_H_
/**
 * @file CNFamilialAnalysisMethodHemizygousParentOfOrigin.h
 *
 * @brief This header contains the CNFamilialAnalysisMethodHemizygousParentOfOrigin class definition.
 */

#include "copynumber/CNFamilialAnalysisMethod.h"

/**
 * @brief  The HemizygousParentOfOrigin Familial analysis method.
 *
 */
class CNFamilialAnalysisMethodHemizygousParentOfOrigin : public CNFamilialAnalysisMethod
{
private:
    int m_iMarkerCountCutoff;
    int m_iCallCutoff;
    int m_iMinGenomicSpan;
    int m_iXChromosome;
    int m_iYChromosome;

public:
    static std::string getType() {return "hemizygous-parent-of-origin";}
    static std::string getDescription() {return "Familial HemizygousParentOfOrigin";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNFamilialAnalysisMethodHemizygousParentOfOrigin();
    virtual ~CNFamilialAnalysisMethodHemizygousParentOfOrigin() {}

    virtual bool isFamilialAnalysis() {return false;}
    virtual bool isSegmentOverlapAnalysis() {return false;}
    virtual bool isSegmentTypeAnalysis() {return true;}

    virtual AffxString getName() {return getType();}

    virtual void run();
};

#endif


