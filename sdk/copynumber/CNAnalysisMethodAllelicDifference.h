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

#ifndef _CNAnalysisMethodAllelicDifference_H_
#define _CNAnalysisMethodAllelicDifference_H_
/**
 * @file CNAnalysisMethodAllelicDifference.h
 *
 * @brief This header contains the CNAnalysisMethodAllelicDifference class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

/**
 * @brief  The AllelicDifference analysis method.
 *
 */
class CNAnalysisMethodAllelicDifference : public CNAnalysisMethod
{
protected:
    double m_dAllelicDifferenceOutlierTrim;

    CNProbeSetArray m_vLocalProbeSets;
    bool m_bLocalProbeSetsDetermined;

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
        if (v.getXDimension() == CNAnalysisMethod::getProbeSets()->getCount()) {v.initialize();}
        else {v.initialize(CNAnalysisMethod::getProbeSets()->getCount());}
    }

    void determineLocalProbeSets();
    CNProbeSetArray* getProbeSets() {return &m_vLocalProbeSets;}

public:
    static std::string getType() {return "allelic-difference";}
    static std::string getDescription() {return "Copynumber AllelicDifference";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    static float getRequiredValue(CNProbeSet* probeset) {return probeset->getAllelicDifference();}
    static void setRequiredValue(CNProbeSet* probeset, float ad) { probeset->setAllelicDifference(ad); }

public:
    CNAnalysisMethodAllelicDifference();
    virtual ~CNAnalysisMethodAllelicDifference() {}

    virtual bool isSegmentTypeAnalysis() {return false;}

    virtual AffxString getName() {return getType();}

    virtual void run();

protected:
    void calculateAllelicDifferences();

    virtual void filterCutOff(
                    CNProbeSetArray* probeSets,
                    int chrStart,
                    int chrEnd,
                    float cutoff3,
                    float cutoff4,
                    vector<bool>& filteredProbeSets
                    );

    virtual void filterNoMansLand(
                    CNProbeSetArray* probeSets,
                    vector<bool>& filteredProbeSets,
                    const vector<pair<int, int> >& stepWindowBds,
                    const vector<vector<double> >& peaks
                    );

    virtual void filterShrinkTowardPeaks(
                    CNProbeSetArray* probeSets,
                    vector<bool>& filteredProbeSets,
                    const vector<pair<int, int> >& stepWindowBds,
                    const vector<vector<double> >& peaks,
                    float shrinkFactor3,
                    float shrinkFactor4,
                    vector<pair<int, double> >& closestPeak,
                    vector<pair<int, double> >& processedValues,
                    bool saveAllelePeaksFlag
                    );

    virtual void ADOutlierTrim(double& ad);
};

#endif


