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

#ifndef _CNAnalysisMethodMosaicism_H_
#define _CNAnalysisMethodMosaicism_H_
/**
 * @file CNAnalysisMethodMosaicism.h
 *
 * @brief This header contains the CNAnalysisMethodMosaicism class definition.
 */

#include "copynumber/CNAnalysisMethodSegment.h"

/**
 * @brief  The Mosaicism analysis method.
 *
 */
class CNAnalysisMethodMosaicism : public CNAnalysisMethodSegment
{
public:
  // set for debugging use.
  std::string m_current_chr;
  std::string m_current_experiment;

    static std::string getType() {return "mosaicism";}
    static std::string getDescription() {return "Copynumber Mosaicism";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodMosaicism();
    virtual ~CNAnalysisMethodMosaicism() {}

    virtual bool isSegmentTypeAnalysis() {return true;}

    virtual AffxString getName() {return getType();}

    virtual void run();

        void determineLocalProbeSets();

        CNProbeSetArray* getProbeSets()
        {
                return &m_vLocalProbeSets;
        }

  void writeRunningMeans(const std::string& path,
                         int cnt,
                         double* pLog2Ratios,
                         double* pRunMean);

private:
    std::vector<double> m_vGainsBoundries;
    std::vector<double> m_vLossesBoundries;
    int m_iMarkerBandwidth;
    int m_iConfidenceWindow;
    int m_iStartIndex;
    int m_iEndIndex;
    bool m_bRunYChromosome;


        CNProbeSetArray  m_vLocalProbeSets;
        bool m_bLocalProbeSetsDetermined;
private:
    void newSegments(unsigned char cChromosome, double* pRunMean, int iProbeSetCount, std::vector<int>& vPositions, CNSegmentArray& vSegments);
    MosaicClass_t getMosaicClass(double dMostExtreme, std::vector<double>& vBoundries);
    float getCalibratedCN(double dCurrentMax);
    float getCalibratedCN(double dCurrentMax, int chromosome);
    //int getBoundaryAdjustment(double dMosaicClass);
    int getBoundaryAdjustment(MosaicClass_t val); 
    void cleanSegments(CNSegmentArray& vSegments);
    void runmed(double* p, int iCount, int iWindowSize);
    void calculateSegmentConfidence(double* pRunMean, CNSegment* p);
};

#endif
