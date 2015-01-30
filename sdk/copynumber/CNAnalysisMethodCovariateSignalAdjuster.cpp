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
 * @file   CNAnalysisMethodCovariateSignalAdjuster.cpp
 *
 * @brief  Class for doing covariate-based signal adjustment
 */

//
#include "copynumber/CNAnalysisMethodCovariateParams.h"
#include "copynumber/CNAnalysisMethodCovariateSignalAdjuster.h"
#include "copynumber/Annotation.h"
//
//#include "chipstream/MedNormTran.h"
//#include "chipstream/SketchQuantNormTran.h"
//
#include "stats/stats.h"
//#include "util/AffxFile.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Verbose.h"
//

using namespace std;

/**
 * Constructor
 *
 */
CNAnalysisMethodCovariateSignalAdjuster::CNAnalysisMethodCovariateSignalAdjuster()
{
    m_alternateIntensities = NULL;
}

CNAnalysisMethodCovariateSignalAdjuster::~CNAnalysisMethodCovariateSignalAdjuster()
{
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodCovariateSignalAdjuster::explainSelf()
{
    CNAnalysisMethodCovariateSignalAdjuster obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodCovariateSignalAdjuster::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  SelfDoc::Opt CovariateOrder = {"order", SelfDoc::Opt::String, "", "", "", "", "Covariate order"};
  opts.push_back(CovariateOrder);

  SelfDoc::Opt BinType = {"bin-type", SelfDoc::Opt::String, "", "", "", "", "Bin type"};
  opts.push_back(BinType);

  SelfDoc::Opt BinCount = {"bin-count", SelfDoc::Opt::String, "", "", "", "", "Bin count"};
  opts.push_back(BinCount);

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
SelfCreate* CNAnalysisMethodCovariateSignalAdjuster::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodCovariateSignalAdjuster* pMethod = new CNAnalysisMethodCovariateSignalAdjuster();

    //pMethod->Order = setupStringParameter("order", "", params, doc, opts);
    //pMethod->BinType = setupStringParameter("bin-type", "", params, doc, opts);
    //pMethod->BinCount = setupStringParameter("bin-count", "", params, doc, opts);

    return pMethod;
}

void CNAnalysisMethodCovariateSignalAdjuster::run()
{
    // Verbose::out(1, "CNAnalysisMethodCovariateSignalAdjuster::run()...");
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Modifying Probe Intensities with Covariate Signal Adjuster.");
    adjustSignals();
}

void CNAnalysisMethodCovariateSignalAdjuster::adjustSignals()
{
    int iXChromosome = m_pEngine->getOptInt("xChromosome");
    int iYChromosome = m_pEngine->getOptInt("yChromosome");
    int numCovariates = getSignalCovariates().size();
    if (numCovariates == 0)
    {
        return;
    }

    // Check if any covariate values were read from the annotation file or covariates file
    if (getProbeSets()->getCount() > 0 && getProbeSets()->getAt(0)->getNumCovariates() == 0)
    {
        if (m_pEngine->getOptBool("disable-covariates-file-warning") == false)
        {
            Verbose::warn(1, "Covariate adjustment is disabled. No covariate values were loaded. annotation-file does not support covariates.");
            m_pEngine->setOpt("disable-covariates-file-warning", "true");
        }
        return;
    }

    Verbose::out(1, "CNAnalysisMethodCovariateSignalAdjuster::run()...");
    Verbose::progressBegin(3, "Run covariate signal adjustment", numCovariates*2, 1, numCovariates*2);

    // Find out if marker-class is one of the covariates
    int markerClassIndex = CovariateParams::mapCovariateNameToIndex("marker-class");

    if (m_pEngine->getOptBool("keep-intermediate-data")) {
        writeIntensities("beforeCovariateAdjust");
    }

    for (int iIndex = 0; iIndex < numCovariates; ++iIndex)
    {
        int covIndex = getSignalCovariates()[iIndex];    // index into covariates vector in CNProbeSet
        int maxFlagIndex = 2;      // the default number of marker classes (i.e. SNP and CN)
        if (markerClassIndex == covIndex)
        {
            // This is the marker-class covariate, so do not distinguish SNP and CN in what follows
            maxFlagIndex = 1;      // do the loop below only once (i.e. it's one class of markers)
        }
        for (int markerClass = 0; markerClass < maxFlagIndex; ++markerClass)   // markerClass 0 = CN, 1 = SNP
        {
            string intermSuffix;
            if (maxFlagIndex == 1) {
                intermSuffix = "_MARKER_CLASS";
            } else if (markerClass == 0) {
                intermSuffix = "_CN";
            } else if (markerClass == 1) {
                intermSuffix = "_SNP";
            } else {
                intermSuffix = "_INVALID_FLAG";
            }

            vector<float> covValues;                 // corresponding covariates for binning
            vector<CNProbe*> probesToProcess;        // probe with valid iIndex'th covariate value
            vector<int> binning;                     // corresponding bin numbers
    
            Verbose::progressStep(3);
    
            for (int pIndex = 0; pIndex < getProbes()->getCount(); ++pIndex)
            {
                CNProbe* pobjProbe = getProbes()->getAt(pIndex);
                CNProbeSet* pobjProbeSet = getProbeSets()->getAt(pobjProbe->getProbeSetIndex());
                // Collect probe sets with non-NaN covariate values
                float covVal = pobjProbeSet->getCovariateValue(covIndex);
                    if ((maxFlagIndex == 1 || doesProcessFlagMatchMarkerClass(markerClass, pobjProbeSet->getProcessFlag())) &&
                        !isNaN(covVal))
                {
                    probesToProcess.push_back(pobjProbe);
                    covValues.push_back(covVal);
                }
            }
    
            binning.assign(covValues.size(), -1);
            int numBins = determineBins(covValues, binning, iIndex);
            vector<int> histogram(numBins);
    
            // Collect the medians for the bins.
            // NB: medians of autosomes only.
            vector<float> medians;
            vector<float> allData;
            try
            {
                for (int binIndex = 0; binIndex < numBins; ++binIndex)
                {
                    vector<float> bData;
                    for (int iIndex = 0; iIndex < probesToProcess.size(); iIndex++)
                    {
                        if (binIndex == 0) histogram[binning[iIndex]] += 1; // make histogram only while processing the first bin.
    
                        if (binning[iIndex] == binIndex)
                        {
                            CNProbeSet *pset = getProbeSets()->getAt(probesToProcess[iIndex]->getProbeSetIndex());
                            if (pset->getChromosome() != iXChromosome && pset->getChromosome() != iYChromosome)
                            {
                                float intensity = getIntensity(probesToProcess[iIndex]);
                                if (!isNaN(intensity))
                                {
                                    bData.push_back(intensity);
                                    allData.push_back(intensity);
                                }
                            }
                        }
                    }
    
                    if (bData.size() > 0)
                    {
                        vector<float>::iterator iter = bData.begin() + bData.size()/2;
                        std::nth_element(bData.begin(), iter, bData.end());
                        medians.push_back(*iter);
                    }
                    else
                    {
                        medians.push_back(0.0);
                    }
                }
            }
            catch(exception& e)
            {
                Err::errAbort(e.what());
            }
    
            // Compute the grand median - I onlyl have to do this once I should try another way
            vector<float>::iterator iter = allData.begin() + allData.size()/2;
            std::nth_element(allData.begin(), iter, allData.end());
            float grandMedian = *iter;
        
            // Adjust signal (note probe outside validProbes[] are not adjusted)
            for (int pIndex = 0; pIndex < probesToProcess.size(); ++pIndex)
            {
                if (medians[binning[pIndex]] != 0.0)
                {
                    float adjSignal = getIntensity(probesToProcess[pIndex]) * grandMedian / medians[binning[pIndex]];
                    if (!isNaN(adjSignal))
                    {
                        setAdjustedIntensity(probesToProcess[pIndex], adjSignal);
                    }
                }
            }
            if (m_pEngine->getOptBool("keep-intermediate-data"))
            {
                writeIntensities("afterCovariateAdjust_" + Convert::toString(iIndex+1) + intermSuffix);
                writeMediansVectorAndMedian("covariateSignalMedians_" + Convert::toString(iIndex+1) + intermSuffix, medians, grandMedian);
            }
        }
    }
    Verbose::progressEnd(3, "Done.");
}

int CNAnalysisMethodCovariateSignalAdjuster::determineBins(vector<float>& covValues, vector<int>& binning, int covariateIndex)
{
    if (CovariateParams::useEquallyPopulatedSignalBins(covariateIndex))
    {
        int numBins = CovariateParams::getNumSignalBins(covariateIndex);
        CNAnalysisMethod::binEqualNumber(covValues, binning, numBins);
        return numBins;
    }
    else if (CovariateParams::useEquallySpacedSignalBins(covariateIndex))
    {
        int numBins = CovariateParams::getNumSignalBins(covariateIndex);
        CNAnalysisMethod::binEqualSpacing(covValues, binning, numBins);
        return numBins;
    }
    else if (CovariateParams::isSignalCovariateDiscrete(covariateIndex))
    {
        return CNAnalysisMethod::covariateIsBinAssignment(covValues, binning);
    }
    else // unknown binning
    {
        // This is a late error.
        Err::errAbort("Unknown signal covariate binning type");
    }
    return 0;
}

std::vector<int>& CNAnalysisMethodCovariateSignalAdjuster::getSignalCovariates()
{
    return CovariateParams::m_vSignalCovariates;
}
