////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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

/**
 * @file   CNAnalysisMethodCovariateSignalAdjuster.h
 *
 * @brief  Class for doing covariate-based signal adjustment
 */
#ifndef _CNAnalysisMethodCovariateSignalAdjuster_H_
#define _CNAnalysisMethodCovariateSignalAdjuster_H_

//#include "chipstream/ChipLayout.h"
//
#include "copynumber/CNAnalysisMethod.h"
//
#include "util/AffxByteArray.h"
#include "util/AffxConv.h"
//#include "util/AffxFile.h"
#include "util/AffxMultiDimensionalArray.h"
#include "util/Convert.h"
#include "util/Util.h"
//
#include <vector>
//


/**
 * @brief CNAnalysisMethodCovariateSignalAdjuster - Class for doing covariate-based signal adjustment
 */
class CNAnalysisMethodCovariateSignalAdjuster : public CNAnalysisMethod
{
public:
    static std::string getType() {return "covariate-signal-adjuster";}
    static std::string getDescription() {return "Covariate Signal Adjuster";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodCovariateSignalAdjuster();
    ~CNAnalysisMethodCovariateSignalAdjuster();

    virtual bool isSegmentTypeAnalysis() {return false;}

    virtual AffxString getName() {return getType();}

    virtual void run();

protected:
    void adjustSignals();

    int determineBins(vector<float>& covValues, vector<int>& binning, int covariateIndex);

    std::vector<int>& getSignalCovariates();

    float getIntensity(CNProbe* probe)
    {
        if (m_alternateIntensities == NULL)
        {
            return probe->getIntensity();
        }
        else
        {
            return m_alternateIntensities->at(probe->getProbeID() - 1);
        }
    }

    void setAdjustedIntensity(CNProbe* probe, float intensity)
    {
        if (m_alternateIntensities == NULL)
        {
            probe->setIntensity(intensity);
        }
        else
        {
            m_alternateIntensities->at(probe->getProbeID() - 1) = intensity;
        }
    }

    bool doesProcessFlagMatchMarkerClass(int markerClass, int processFlag)
    {
        // markerClass 0 = CN, 1 = SNP
        //
        // Marker class: process flag(s)
        // SNP: 2,3,4
        // CN: 1

        if (markerClass == 0 && processFlag == 1)
        {
            // CN
            return true;
        }
        else if (markerClass == 1 && (processFlag == 2 || processFlag == 3 || processFlag == 4))
        {
            // SNP
            return true;
        }
        return false;
    }

public:
    // This method was added divert the intensity adjustment from the CNProbe::setIntensity to the passed in vectors.
    // Assumes probe intensities are in ProbeID order
    void setAlternateIntensities(std::vector<float>& alternateIntensities)
    {
        m_alternateIntensities = &alternateIntensities;
    }

protected:
    std::vector<float>* m_alternateIntensities;
};

#endif  //_CNAnalysisMethodCovariateSignalAdjuster_H_
