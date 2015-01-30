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

#ifndef _CNFamilialAnalysisMethodSegmentOverlap_H_
#define _CNFamilialAnalysisMethodSegmentOverlap_H_
/**
 * @file CNFamilialAnalysisMethodSegmentOverlap.h
 *
 * @brief This header contains the CNFamilialAnalysisMethodSegmentOverlap class definition.
 */

#include "copynumber/CNFamilialAnalysisMethod.h"

/**
 * @brief  The SegmentOverlap Familial analysis method.
 *
 */
class CNFamilialAnalysisMethodSegmentOverlap : public CNFamilialAnalysisMethod
{
private:
    float m_fSegmentOverlapThreshold;

public:
    static std::string getType() {return "segment-overlap";}
    static std::string getDescription() {return "Familial SegmentOverlap";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNFamilialAnalysisMethodSegmentOverlap();
    virtual ~CNFamilialAnalysisMethodSegmentOverlap() {}

    virtual bool isFamilialAnalysis() {return false;}
    virtual bool isSegmentOverlapAnalysis() {return true;}
    virtual bool isSegmentTypeAnalysis() {return false;}

    virtual AffxString getName() {return getType();}

    virtual void run();
};

#endif


