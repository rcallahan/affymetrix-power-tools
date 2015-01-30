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

#ifndef _CNAnalysisMethodLOHCytoScan_H_
#define _CNAnalysisMethodLOHCytoScan_H_
/**
 * @file CNAnalysisMethodLOHCytoScan.h
 *
 * @brief This header contains the CNAnalysisMethodLOHCytoScan class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

/**
 * @brief  The LOH analysis method.
 *
 */
class CNAnalysisMethodLOHCytoScan : public CNAnalysisMethod
{
public:
    //LOH
    double m_dLohCNErrorRate;
    double m_dLohCNBeta;
    double m_dLohCNAlpha;
    int m_iLohCNSeparation;
    int m_iLohCNMinimumMarkerCount;
    double m_dLohCNNoCallThreshold;
    int m_iLohCNMinGenomicSpan;

    CNProbeSetArray  m_vLocalProbeSets;
    bool m_bLocalProbeSetsDetermined;


public:
    CNAnalysisMethodLOHCytoScan();
    virtual ~CNAnalysisMethodLOHCytoScan() {}

    static std::string getType() {return "lohCytoScan";}
    static std::string getDescription() {return "CopyNumber LOH CytoScan";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    virtual AffxString getName() {return getType();}
    virtual void run();
    virtual bool isSegmentTypeAnalysis() {return false;}

    CNProbeSetArray* getProbeSets() {return &m_vLocalProbeSets;}
public:
    bool lohPreProcessing(std::vector<char>& vHomHet, int& iMarkerCount, int& iHetCutoff);
    bool lohFind(std::vector<char>& vGenotypeCalls, std::vector<int>& vPositions, int iMarkerCount, int iHetCutoff, int iMinGenomicSpan, std::vector<char>& vLoh);

private:
    void imputeMissingLOHCalls();
    void determineLocalProbeSets();


};

#endif


