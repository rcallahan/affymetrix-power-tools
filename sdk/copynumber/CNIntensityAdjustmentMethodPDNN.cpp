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

/**
 * @file CNIntensityAdjustmentMethodPDNN.cpp
 *
 * @brief This file contains the CNIntensityAdjustmentMethodPDNN class members.
 */
#include "copynumber/CNIntensityAdjustmentMethodPDNN.h"
//
#include "copynumber/CNAnalysisMethodFactory.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/AffxStatistics.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNIntensityAdjustmentMethodPDNN::explainSelf()
{
    CNIntensityAdjustmentMethodPDNN obj;
    SelfDoc doc;
    doc.setDocName(obj.getType());
    doc.setDocDescription(obj.getDescription());
    doc.setDocOptions(obj.getDefaultDocOptions());
    return doc;
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> CNIntensityAdjustmentMethodPDNN::getDefaultDocOptions()
{
    std::vector<SelfDoc::Opt> opts;

    // SelfDoc::Opt(name, type, value, default, min, max, description)
    SelfDoc::Opt PIBinCount = {"predicted-intensity-bin-count", SelfDoc::Opt::Integer, "20", "20", "1", "126", "Predicted Intensity bin count."};
    opts.push_back(PIBinCount);

    SelfDoc::Opt GCBinCount = {"gc-bin-count", SelfDoc::Opt::Integer, "20", "20", "1", "126", "Predicted Intensity bin count."};
    opts.push_back(GCBinCount);

    SelfDoc::Opt ResidualTrim = {"residual-trim", SelfDoc::Opt::Double, "2.0", "2.0", "NA", "NA", "Residual Trim value."};
    opts.push_back(ResidualTrim);

    return opts;
}

/**
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 *
 * @param param - Map of key/value pairs to initialize the object.
 *
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate* CNIntensityAdjustmentMethodPDNN::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNIntensityAdjustmentMethodPDNN* pMethod = new CNIntensityAdjustmentMethodPDNN();
    std::string strPrefix = getPrefix();

    pMethod->m_iPIBinCount = setupIntParameter("predicted-intensity-bin-count", strPrefix, params, doc, opts);
    pMethod->m_iGCBinCount = setupIntParameter("gc-bin-count", strPrefix, params, doc, opts);
    pMethod->m_dResidualTrim = setupDoubleParameter("residual-trim", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNIntensityAdjustmentMethodPDNN::CNIntensityAdjustmentMethodPDNN()
{
    m_iPIBinCount = 0;
    m_iGCBinCount = 0;
    m_dResidualTrim = 0;
    m_pvProbes = NULL;
        m_bLocalProbeSetsDetermined=false;
}

void CNIntensityAdjustmentMethodPDNN::determineLocalProbeSets()
{
        if(m_bLocalProbeSetsDetermined)
        {
                return;
        }

        int iNumberOfProbeSets=CNAnalysisMethod::getProbeSets()->getCount();
        for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
        {
                if( CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAll() )
                {
                        getProbeSets()->add( CNAnalysisMethod::getProbeSets()->getAt(iIndex));
                }
        }
        m_bLocalProbeSetsDetermined=true;
}




/**
 * @brief Run the analysis.
 */
void CNIntensityAdjustmentMethodPDNN::run()
{
    Verbose::out(1, "CNIntensityAdjustmentMethodPDNN::run(...) start");
    isSetup();

//  First filter out the probes so that they all come from probeSets that have process flag 1,2 or 3.
//  Note that this means that they will have been read from the annotations file.

        determineLocalProbeSets();

        std::vector<CNProbe*> pvPDNNProbes;
        for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
        {
                if( CNAnalysisMethod::getProbeSets()->getAt(m_pvProbes->getAt(iIndex)->getProbeSetIndex())->processAll() )
                {
                       pvPDNNProbes.push_back(m_pvProbes->getAt(iIndex));
                }
        }

    std::vector<int> vBinIndexes(pvPDNNProbes.size());
    if (getExperiment()->getIndex() == 0)
    {
        //Predicted Intensity bins
        std::vector<int> vPIBinIndexes(pvPDNNProbes.size());
        std::vector<float> vPredictedIntensities(pvPDNNProbes.size());
        for (int iIndex = 0; iIndex < pvPDNNProbes.size(); iIndex++)
        {
            vPredictedIntensities[iIndex] = pvPDNNProbes[iIndex]->getPredictedIntensity();
        }
        bin(vPredictedIntensities, vPIBinIndexes, m_iPIBinCount);

        // GC bins
        int iGCCount = 0;
        for (int iIndex = 0; (iIndex < getProbeSets()->getCount()); iIndex++)
        {
            if (getProbeSets()->getAt(iIndex)->getGCContent() != 0 ) {iGCCount++;}
        }
        AffxMultiDimensionalArray<float> vGCTemp(iGCCount);
        iGCCount = 0;
        for (int iIndex = 0; (iIndex < getProbeSets()->getCount()); iIndex++)
        {
            if (getProbeSets()->getAt(iIndex)->getGCContent() != 0)
            {
                vGCTemp.set(iGCCount, getProbeSets()->getAt(iIndex)->getGCContent());
                iGCCount++;
            }
        }
        double dGCMedian = vGCTemp.median();
        std::vector<int> vGCBinIndexes(pvPDNNProbes.size());
        std::vector<float> vGCContent(pvPDNNProbes.size());
        for (int iIndex = 0; (iIndex < pvPDNNProbes.size()); iIndex++)
        {
            if (pvPDNNProbes[iIndex]->getProbeSetIndex() == -1)
            {
                vGCContent[iIndex] = dGCMedian;
            }
            else
            {
                vGCContent[iIndex] = CNAnalysisMethod::getProbeSets()->getAt(pvPDNNProbes[iIndex]->getProbeSetIndex()) ->getGCContent();
                if (vGCContent[iIndex] == 0) {vGCContent[iIndex] = dGCMedian;}
            }
        }
        bin(vGCContent, vGCBinIndexes, m_iGCBinCount);

        // Set overall bin index.
        // int iBinCount = m_iPIBinCount * m_iGCBinCount;
        for (int iRowIndex = 0; (iRowIndex < pvPDNNProbes.size()); iRowIndex++)
        {
            vBinIndexes[iRowIndex] = ((vPIBinIndexes[iRowIndex] * m_iGCBinCount) + vGCBinIndexes[iRowIndex]);
            pvPDNNProbes[iRowIndex]->setPDNNBin(vBinIndexes[iRowIndex]);
        }
    }
    else
    {
        for (int iRowIndex = 0; (iRowIndex < pvPDNNProbes.size()); iRowIndex++)
        {
            vBinIndexes[iRowIndex] = pvPDNNProbes[iRowIndex]->getPDNNBin();
        }
    }

    // Apply adjustment.
    int iBinCount = m_iPIBinCount * m_iGCBinCount;
    Verbose::progressBegin(1, "CNAnalysisMethodChipstream::adjustIntensitiesPDNN(...) ", (iBinCount / 100), 1, (iBinCount / 100));
    for (int iBinIndex = 0; (iBinIndex < iBinCount); iBinIndex++)
    {
        if ((iBinIndex % 100) == 0) {Verbose::progressStep(1);}
        // Get autosome count per bin.
        int iAutosomeProbeCount = 0;
        for (int iRowIndex = 0; (iRowIndex < pvPDNNProbes.size()); iRowIndex++)
        {
            if (vBinIndexes[iRowIndex] == iBinIndex)
            {
                if (CNAnalysisMethod::getProbeSets()->getAt(pvPDNNProbes[iRowIndex]->getProbeSetIndex())->getChromosome() < m_iXChromosome) // is autosome.
                {
                    iAutosomeProbeCount++;
                }
            }
        }
        // Calculate autosome median for each bin.
        AffxMultiDimensionalArray<float> vNormalizedIntensities(iAutosomeProbeCount);
        AffxMultiDimensionalArray<float> vMedianNormalizedIntensities(iAutosomeProbeCount);
        int iIndex = 0;
        for (int iRowIndex = 0; (iRowIndex < pvPDNNProbes.size()); iRowIndex++)
        {
            if (vBinIndexes[iRowIndex] == iBinIndex)
            {
                if (CNAnalysisMethod::getProbeSets()->getAt(pvPDNNProbes[iRowIndex]->getProbeSetIndex())->getChromosome() < m_iXChromosome) // is autosome.
                {
                    vNormalizedIntensities.set(iIndex, log(pvPDNNProbes[iRowIndex]->getIntensity()));
                    vMedianNormalizedIntensities.set(iIndex, log(pvPDNNProbes[iRowIndex]->getMedianIntensity()));
                    iIndex++;
                }
            }
        }
        // Perform actual intensity adjustment.
        double dMedianNormalizedIntensities = vNormalizedIntensities.median();
        double dMedianMedianNormalizedIntensities = vMedianNormalizedIntensities.median();
        float fResidual = 0;
        for (int iRowIndex = 0; (iRowIndex < pvPDNNProbes.size()); iRowIndex++)
        {
            if (vBinIndexes[iRowIndex] == iBinIndex)
            {
                if (iAutosomeProbeCount == 0) {pvPDNNProbes[iRowIndex]->setResidual(0);}
                else
                {
                    float fIntensity = log(pvPDNNProbes[iRowIndex]->getIntensity());
                    fIntensity = fIntensity - dMedianNormalizedIntensities + dMedianMedianNormalizedIntensities;
                    fResidual = fIntensity - dMedianMedianNormalizedIntensities;
                    if (fResidual > m_dResidualTrim) {fResidual = (float)m_dResidualTrim;}
                    else if (fResidual < -m_dResidualTrim) {fResidual = (float)-m_dResidualTrim;}
                    pvPDNNProbes[iRowIndex]->setIntensity(exp(fIntensity));
                    pvPDNNProbes[iRowIndex]->setResidual(fResidual);
                }
            }
        }
    }

        if (m_pEngine->getOptBool("keep-intermediate-data"))
        {
                writeIntensities("pdnn");
        }

    Verbose::progressEnd(1, "Done");

/*
    if (m_pEngine->getOptInt("array-size") == 2015) // Full array.
    {
        std::vector<float> vBackground(m_pvProbes->getCount());
        HighPassFilter(*m_pEngine).run(vProbeIDs, vResiduals, vBackground);
        for (int iRowIndex = 0; (iRowIndex < m_pvProbes->getCount()); iRowIndex++)
        {
            CNProbe* pobjProbe = m_pvProbes->getAt(iRowIndex);
            float fIntensity = log(pobjProbe->getIntensity());
            fIntensity -= vBackground[iRowIndex];
            pobjProbe->setIntensity(exp(fIntensity));
        }
    }
*/
    Verbose::out(1, "CNIntensityAdjustmentMethodPDNN::run(...) end");
}
