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

#ifndef _CNAnalysisMethodLog2Ratio_H_
#define _CNAnalysisMethodLog2Ratio_H_
/**
 * @file CNAnalysisMethodLog2Ratio.h
 *
 * @brief This header contains the CNAnalysisMethodLog2Ratio class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

/**
 * @brief  The Log2Ratio analysis method.
 *
 */
class CNAnalysisMethodLog2Ratio : public CNAnalysisMethod
{
private:
    bool m_bLocalProbeSetsDetermined;

protected:
    double m_dYTarget;
    bool m_bGCCorrection;
    bool m_bMedianAutosomeMedianNormalization;
    int m_iMedianSmoothWindowSize;
    double m_dTrimLow;
    double m_dTrimHigh;


    CNProbeSetArray  m_vLocalProbeSets;

    AffxMultiDimensionalArray<float> m_v1;
    AffxMultiDimensionalArray<float> m_v2;
    AffxMultiDimensionalArray<float> m_v3;
    AffxMultiDimensionalArray<float> m_v4;
    AffxMultiDimensionalArray<float> m_v5;

    AffxMultiDimensionalArray<float>& getV1() {initialize(m_v1); return m_v1;}
    AffxMultiDimensionalArray<float>& getV2() {initialize(m_v2); return m_v2;}
    AffxMultiDimensionalArray<float>& getV3() {initialize(m_v3); return m_v3;}
    AffxMultiDimensionalArray<float>& getV4() {initialize(m_v4); return m_v4;}
    AffxMultiDimensionalArray<float>& getV5() {initialize(m_v5); return m_v5;}

    /**
     * @brief Initialize the vector
     * @param AffxMultiDimensionalArray<float>& - The vector to initialize
     */
    void initialize(AffxMultiDimensionalArray<float>& v)
    {
        if (v.getXDimension() == getProbeSets()->getCount()) {v.initialize();}
        else {v.initialize(getProbeSets()->getCount());}
    }

public:
    static std::string getType() {return "log2-ratio";}
    static std::string getDescription() {return "Copynumber Log2Ratio";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodLog2Ratio();
    virtual ~CNAnalysisMethodLog2Ratio() {}

    virtual bool isSegmentTypeAnalysis() {return false;}

    virtual AffxString getName() {return getType();}

    virtual void run();

    void determineLocalProbeSets();

    CNProbeSetArray* getProbeSets()
    {
            return &m_vLocalProbeSets;
    }



protected:
    double calculateIqr();
    void calculateLog2Ratios(bool adjustYTarget = true);
    void GCCorrection(CNExperiment* pobjExperiment);
    void processMarkerGroupGCCorrection(bool bSnp, bool bSty, bool bNsp, std::vector<int>& vBinIndexes);
    inline bool isInGCCorrectionGroup(CNProbeSet& objProbeSet, int iLastAutosomeChromosome, bool bSnp, bool bSty, bool bNsp, int iBinIndex);
    void calculateMedianAutosomeMedian(CNExperiment* pobjExperiment);
    void normalizeLog2Ratios(CNExperiment* pobjExperiment);
    void calculateQCMetrics(CNExperiment* pobjExperiment);
    void trim();
    void applyWaveCorrection();
    void adjustUsingCovariateLRAdjustment();
    void applyYTargetAdjustment();
    void highPassFilterLog2Ratios();
};

#endif
