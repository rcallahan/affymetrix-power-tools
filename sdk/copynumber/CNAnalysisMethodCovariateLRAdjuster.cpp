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
 * @file   CNAnalysisMethodCovariateLRAdjuster.cpp
 *
 * @brief  Class for doing covariate-based log2ratio adjustment
 */

//
#include "copynumber/CNAnalysisMethodCovariateParams.h"
#include "copynumber/CNAnalysisMethodCovariateLRAdjuster.h"
#include "copynumber/Annotation.h"
//
#include "stats/stats.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Verbose.h"
//

using namespace std;

/**
 * Constructor
 *
 */
CNAnalysisMethodCovariateLRAdjuster::CNAnalysisMethodCovariateLRAdjuster()
{
}

CNAnalysisMethodCovariateLRAdjuster::~CNAnalysisMethodCovariateLRAdjuster()
{
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodCovariateLRAdjuster::explainSelf()
{
    CNAnalysisMethodCovariateLRAdjuster obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodCovariateLRAdjuster::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  SelfDoc::Opt CovariateOrder = {"order", SelfDoc::Opt::String, "", "", "", "", "Covariate order"};
  opts.push_back(CovariateOrder);

  SelfDoc::Opt BinType = {"bin-type", SelfDoc::Opt::String, "", "", "", "", "Bin type"};
  opts.push_back(BinType);

  SelfDoc::Opt BinCount = {"bin-count", SelfDoc::Opt::String, "", "", "", "", "Bin count"};
  opts.push_back(BinCount);

  SelfDoc::Opt IQRScaling = {"iqr-scaling", SelfDoc::Opt::String, "", "", "", "", "Apply IQR scaling"};
  opts.push_back(IQRScaling);

  SelfDoc::Opt SubtractFromXY = {"subtract-from-XY", SelfDoc::Opt::String, "", "", "", "", "Subtract median from X and Y."};
  opts.push_back(SubtractFromXY);

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
SelfCreate* CNAnalysisMethodCovariateLRAdjuster::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodCovariateLRAdjuster* pMethod = new CNAnalysisMethodCovariateLRAdjuster();

    //pMethod->Order = setupStringParameter("order", "", params, doc, opts);
    //pMethod->BinType = setupStringParameter("bin-type", "", params, doc, opts);
    //pMethod->BinCount = setupStringParameter("bin-count", "", params, doc, opts);

    return pMethod;
}

// TODO: this needs to be templatized and moved to CNAnalysisMethod.h
// to avoid clonig CNReferenceMethodWaveCorrection::covariateAdjustLog2Ratio()
//
void CNAnalysisMethodCovariateLRAdjuster::run()
{
    int iXChromosome = m_pEngine->getOptInt("xChromosome");
    int iYChromosome = m_pEngine->getOptInt("yChromosome");
    int numCovariates = CovariateParams::m_vLRCovariates.size();
    if (numCovariates == 0)
    {
        return;
    }

    // Check if any covariate values were read from the annotation file or covariates file
    if (getProbeSets()->size() > 0 && getProbeSets()->getAt(0)->getNumCovariates() == 0)
    {
        if (m_pEngine->getOptBool("disable-covariates-file-warning") == false)
        {
            Verbose::warn(1, "Covariate adjustment is disabled. No covariate values were loaded. annotation-file does not support covariates.");
            m_pEngine->setOpt("disable-covariates-file-warning", "true");
        }
        return;
    }

    const float IQRepsilon = 1.0e-10;

    Verbose::out(1, "CNAnalysisMethodCovariateLRAdjuster::run()...");
    Verbose::progressBegin(3, "Run covariate log2 ratio adjustment", numCovariates*2, 1, numCovariates*2);

    char validProcessFlags[] = { 1, 3 };      // CN and SNP, resp.

    // Find out if marker-class is one of the covariates
    int markerClassIndex = CovariateParams::mapCovariateNameToIndex("marker-class");

    if (m_pEngine->getOptBool("keep-intermediate-data")) {
        writeLog2Ratios("beforeCovariateAdjust");
    }

    for (int iCovariateIndex = 0; iCovariateIndex < numCovariates; iCovariateIndex++)
    {
        int iTransIndex = CovariateParams::m_vLRCovariates[iCovariateIndex];
        int maxFlagIndex = 2;      // the default number of marker classes (i.e. SNP and CN)
        if (markerClassIndex == iTransIndex)
        {
            // This is the marker-class covariate, so do not distinguish SNP and CN in what follows
            maxFlagIndex = 1;      // do the loop below only once (i.e. it's one class of markers)
        }
        for (int iFlagIndex = 0; iFlagIndex < maxFlagIndex; iFlagIndex++)
        {
            string intermSuffix;
            if (maxFlagIndex == 1) {
                intermSuffix = "_MARKER_CLASS";
            } else if (validProcessFlags[iFlagIndex] == 1) {
                intermSuffix = "_CN";
            } else if (validProcessFlags[iFlagIndex] == 3) {
                intermSuffix = "_SNP";
            } else {
                intermSuffix = "_INVALID_FLAG";
            }

            vector<float> covariates;
            vector<std::pair<int, int> > binInfo;
            vector<int> probeSetsToProcess;

            Verbose::progressStep(3);

            // NB: local probe sets of CNAnalysisMethodLog2Ratio (process flag = 1 or 3)
            for (int iPSIndex = 0; iPSIndex < getProbeSets()->size(); iPSIndex++)
            {
                CNProbeSet *pset = (*getProbeSets())[iPSIndex];
                float covValue = pset->getCovariateValue(iTransIndex);

                if ((maxFlagIndex == 1 || pset->getProcessFlag() == validProcessFlags[iFlagIndex]) &&
                    covValue == covValue)
                {
                    covariates.push_back(covValue);
                    binInfo.push_back(std::make_pair(0, iPSIndex));
                    probeSetsToProcess.push_back(iPSIndex);
                }
            }
            vector<int> binning(covariates.size());
            int numBins = determineBins(covariates, binning, iCovariateIndex);

            // fill in the binning info
            for (int i = 0; i < binInfo.size(); i++) {
                binInfo[i].first = binning[i];
            }

            // figure out bin membership all at once
            std::sort(binInfo.begin(), binInfo.end(), binCompare());
    
            std::vector<int> binBoundaries(numBins + 1);
            binBoundaries[0] = 0;
            for (int i = 1; i < numBins; i++)
            {
                std::pair<int, int> seek = std::make_pair(i, 0);
                std::vector<std::pair<int, int> >::iterator bdry =
                            std::lower_bound(binInfo.begin(), binInfo.end(), seek, binCompare());
    
                binBoundaries[i] = bdry - binInfo.begin();
            }
            binBoundaries[numBins] = binInfo.size();   // sentinel

            // Adjust log2 ratios
            vector<float> mediansLR(numBins);
            vector<float> IQRanges(numBins);

            for (int iBinIndex = 0; iBinIndex < numBins; iBinIndex++)
            {
                // Find bin median of autosome log2 ratios
                vector<float> log2ratios;

                for (
                    int iIndex = binBoundaries[iBinIndex];
                    iIndex < binBoundaries[iBinIndex + 1];
                    iIndex++
                ){
                    CNProbeSet *pset = (*getProbeSets())[binInfo[iIndex].second];
                    if (pset->getChromosome() != iXChromosome && pset->getChromosome() != iYChromosome)
                    {
                        if (pset->getLog2Ratio() == pset->getLog2Ratio())
                        {
                            log2ratios.push_back(pset->getLog2Ratio());
                        }
                    }
                }

                int binSize = log2ratios.size();
                if (binSize == 0)
                {
                    mediansLR[iBinIndex] = 0.0;
                    IQRanges[iBinIndex] = numeric_limits<float>::quiet_NaN();   // invalid value to flag an empty bin
                    continue;
                }
                std::nth_element(log2ratios.begin(), log2ratios.begin() + binSize/2, log2ratios.end());
                mediansLR[iBinIndex] = log2ratios[binSize/2];

                std::nth_element(log2ratios.begin(), log2ratios.begin() + binSize/4, log2ratios.end());
                float lowerQuartile = log2ratios[binSize/4];

                std::nth_element(log2ratios.begin(), log2ratios.begin() + 3*binSize/4, log2ratios.end());
                IQRanges[iBinIndex] = log2ratios[3*binSize/4] - lowerQuartile;
            }

            bool subtractMedianFromXY = (CovariateParams::m_vLRSubractFromXY[iCovariateIndex] == "on");
    
            // Subtract the medians from markers in this bin (including X and Y)
            for (int iIndex = 0; iIndex < probeSetsToProcess.size(); iIndex++)
            {
                CNProbeSet *pset = (*getProbeSets())[probeSetsToProcess[iIndex]];
                if (subtractMedianFromXY || (pset->getChromosome() != iXChromosome && pset->getChromosome() != iYChromosome))
                {
                    if (pset->getLog2Ratio() == pset->getLog2Ratio())
                    {
                        pset->setLog2Ratio(pset->getLog2Ratio() - mediansLR[binning[iIndex]]);
                    }
                }
            }

            // If requested, perform IQR scaling
            if (CovariateParams::m_vIQRScaling[iCovariateIndex] == "on")
            {
                vector<float> validIQRanges;
                for (int iBinIndex = 0; iBinIndex < numBins; iBinIndex++)
                {
                    if (IQRanges[iBinIndex] == IQRanges[iBinIndex])
                    {
                        validIQRanges.push_back(IQRanges[iBinIndex]);
                    }
                }
                int numValidRanges = validIQRanges.size();
                std::nth_element(validIQRanges.begin(), validIQRanges.begin() + numValidRanges/2, validIQRanges.end());
                float medianIQRange = validIQRanges[numValidRanges/2];

                vector<float> scalingFactors(numBins);
                for (int iBinIndex = 0; iBinIndex < numBins; iBinIndex++)
                {
                    if(IQRanges[iBinIndex] == IQRanges[iBinIndex])
                    {
                        if (IQRanges[iBinIndex] < IQRepsilon)
                        {
                            scalingFactors[iBinIndex] = 1.0;
                        } else {
                            scalingFactors[iBinIndex] = medianIQRange/IQRanges[iBinIndex];
                        }
                    }
                }
    
                for (int iIndex = 0; iIndex < probeSetsToProcess.size(); iIndex++)
                {
                    CNProbeSet *pset = (*getProbeSets())[probeSetsToProcess[iIndex]];
    
                    // Turn off iqr-scaling for X/Y if subtractMedianFromXY == false - per Ben 1/28/2011
                    if (subtractMedianFromXY || (pset->getChromosome() != iXChromosome && pset->getChromosome() != iYChromosome))
                    {
                        if (pset->getLog2Ratio() == pset->getLog2Ratio())
                        {
                            pset->setLog2Ratio(pset->getLog2Ratio() * scalingFactors[binning[iIndex]]);
                        }
                    }
                }
            }
            if (m_pEngine->getOptBool("keep-intermediate-data"))
            {
                writeLog2Ratios("afterCovariateAdjust_" + Convert::toString(iCovariateIndex+1) + intermSuffix);
                writeMediansVector("covariateL2rMedians_" + Convert::toString(iCovariateIndex+1) + intermSuffix, mediansLR);
            }
        }
    }
    Verbose::progressEnd(3, "Done.");
}

int CNAnalysisMethodCovariateLRAdjuster::determineBins(vector<float>& covValues, vector<int>& binning, int covariateIndex)
{
    if (CovariateParams::useEquallyPopulatedLRBins(covariateIndex))
    {
        int numBins = CovariateParams::getNumLRBins(covariateIndex);
        CNAnalysisMethod::binEqualNumber(covValues, binning, numBins);
        return numBins;
    }
    else if (CovariateParams::useEquallySpacedLRBins(covariateIndex))
    {
        int numBins = CovariateParams::getNumLRBins(covariateIndex);
        CNAnalysisMethod::binEqualSpacing(covValues, binning, numBins);
        return numBins;
    }
    else if (CovariateParams::isLRCovariateDiscrete(covariateIndex))
    {
        return CNAnalysisMethod::covariateIsBinAssignment(covValues, binning);
    }
    else // unknown binning
    {
        // This is a late error.
        Err::errAbort("Unknown log2ratio covariate binning type");
    }
    return 0;
}

std::vector<int>& CNAnalysisMethodCovariateLRAdjuster::getLRCovariates()
{
    return CovariateParams::m_vLRCovariates;
}
