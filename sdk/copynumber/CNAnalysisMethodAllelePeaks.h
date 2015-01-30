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

#ifndef _CNAnalysisMethodAllelePeaks_H_
#define _CNAnalysisMethodAllelePeaks_H_
/**
 * @file CNAnalysisMethodAllelePeaks.h
 *
 * @brief This header contains the CNAnalysisMethodAllelePeaks class definition.
 */

#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNReporter.h"
//
#include <cassert>
#include <cmath>
#include <vector>
//

/**
 * @brief  The Allele Peaks reporter method.
 *
 */
class CNAnalysisMethodAllelePeaks : public CNAnalysisMethod
{
private:
    bool m_bSNPAlreadyLoaded;

    CNProbeSetArray  m_vLocalProbeSets;
    bool m_bLocalProbeSetsDetermined;

    void computeSNPQC(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets);
    void determineLocalProbeSets();
    CNProbeSetArray* getProbeSets();
    void fillSnpChrBounds();

public:
    static std::string getType() {return "allele-peaks";}
    static std::string getDescription() {return "Copynumber AllelePeaks";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    static float getRequiredValue(CNProbeSet* probeset) {return probeset->getSCAR();}
    static void setRequiredValue(CNProbeSet* probeset, float scar) {probeset->setSCAR(scar);}

public:
    CNAnalysisMethodAllelePeaks();
    virtual ~CNAnalysisMethodAllelePeaks() {}

    virtual bool isSegmentTypeAnalysis() {return false;}

    virtual AffxString getName() {return getType();}
    virtual void run();

protected:
    int getMaxProbeSetNameLength()
    {
        int iMaxProbeSetNameLength = 0;
        for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
            iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
        }
        return iMaxProbeSetNameLength;
    }

    bool loadSnpReferenceFile(const AffxString& strFileName);

    virtual void filterCutOff(
                        CNProbeSetArray* probeSets,
                        int chrStart, int chrEnd,
                        float cutoff3,
                        float cutoff4,
                        std::vector<bool>& filteredProbeSets
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

};

#endif


