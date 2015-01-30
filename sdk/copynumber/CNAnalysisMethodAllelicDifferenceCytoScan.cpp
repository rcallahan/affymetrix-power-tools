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
 * @file CNAnalysisMethodAllelicDifferenceCytoScan.cpp
 *
 * @brief This file contains the CNAnalysisMethodAllelicDifferenceCytoScan class members.
 */
#include "copynumber/CNAnalysisMethodAllelicDifferenceCytoScan.h"
#include "copynumber/CNAnalysisMethodCovariateParams.h"
//
#include "chipstream/HomHiLoCelListener.h"
#include "util/Fs.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodAllelicDifferenceCytoScan::explainSelf()
{
    CNAnalysisMethodAllelicDifferenceCytoScan obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodAllelicDifferenceCytoScan::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  SelfDoc::Opt opt1 = {"outlier-trim", SelfDoc::Opt::Double, "3.0", "3.0", "NA", "NA", "AllelicDifference Outlier Trim"};
  opts.push_back(opt1);

  SelfDoc::Opt Step = {"step", SelfDoc::Opt::Integer, "20", "20", "1", "NA", "AllelicDifference step size"};
  opts.push_back(Step);

  SelfDoc::Opt Window = {"window", SelfDoc::Opt::Integer, "100", "100", "30", "NA", "AllelicDifference window"};
  opts.push_back(Window);

  SelfDoc::Opt PointCount = {"point-count", SelfDoc::Opt::Integer, "128", "128", "30", "1024", "AllelicDifference number of points"};
  opts.push_back(PointCount);

  SelfDoc::Opt Bandwidth = {"bandwidth", SelfDoc::Opt::Double, "0.25", "0.25", "0", "1", "AllelicDifference bandwidth"};
  opts.push_back(Bandwidth);

  SelfDoc::Opt Cutoff = {"cutoff", SelfDoc::Opt::Double, "0.05", "0.05", "0", "0.5", "AllelicDifference cutoff"};
  opts.push_back(Cutoff);

  SelfDoc::Opt CleanThreshold = {"clean-threshold", SelfDoc::Opt::Double, "0.35", "0.35", "0", "1.0", "AllelicDifference clean threshold"};
  opts.push_back(CleanThreshold);

  SelfDoc::Opt Symmetry = {"symmetry", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA", "AllelicDifference SCAR mirror flag"};
  opts.push_back(Symmetry);

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
SelfCreate* CNAnalysisMethodAllelicDifferenceCytoScan::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodAllelicDifferenceCytoScan* pMethod = new CNAnalysisMethodAllelicDifferenceCytoScan();
    std::string strPrefix = getPrefix();

    pMethod->m_dAllelicDifferenceOutlierTrim = setupDoubleParameter("outlier-trim", strPrefix, params, doc, opts);
    pMethod->m_iStep = setupIntParameter("step", strPrefix, params, doc, opts);
    pMethod->m_iWindow = setupIntParameter("window", strPrefix, params, doc, opts);
    pMethod->m_iPointCount = setupIntParameter("point-count", strPrefix, params, doc, opts);
    pMethod->m_dBandwidth = setupDoubleParameter("bandwidth", strPrefix, params, doc, opts);
    pMethod->m_dCutoff = setupDoubleParameter("cutoff", strPrefix, params, doc, opts);
    pMethod->m_dCleanThreshold = setupDoubleParameter("clean-threshold", strPrefix, params, doc, opts);
    pMethod->m_bSymmetry = setupBoolParameter("symmetry", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodAllelicDifferenceCytoScan::run()
{
    Verbose::out(1, "CNAnalysisMethodAllelicDifferenceCytoScan::run(...) start");
    isSetup();
//    Verbose::progressBegin(1, "CNAnalysisMethodAllelicDifference::run(...) ", 1, 1, 1);
//    Verbose::progressStep(1);
    calculateAllelicDifferences();
    removeAllelicDifferencesModulation();
    postprocessAllelicDifferences();
//    Verbose::progressEnd(1, "Done");
    Verbose::out(1, "CNAnalysisMethodAllelicDifferenceCytoScan::run(...) end");
}

void CNAnalysisMethodAllelicDifferenceCytoScan::removeAllelicDifferencesModulation()
{
    getEngine()->setOpt("AD-modulation-no-single-peak", "false");

    if (CovariateParams::m_vAPCovariates.empty()) {
        return;
    }

    Verbose::out(1, "Running allelic differences modulation correction...");

    // for intermediate data
    const bool keepIntermediateData = m_pEngine->getOptBool("keep-intermediate-data");
    const double xiLimit = 2.5;

    bool AD_modulation_no_single_peak_state = false;

    determineLocalProbeSets();
    fillChrBoundsImpl(getProbeSets(), m_chrBounds);

    std::vector<double> weights;
    computeWeights(getProbeSets(), weights);

    int numCovariates = CovariateParams::m_vAPCovariates.size();
    for (int iCovariateIndex = 0; iCovariateIndex < numCovariates; iCovariateIndex++)
    {
        m_vCoarseAllelePeaks.clear();

        std::vector<double> allelicDiffs;

        // find master peaks
        std::vector<double> master_xi;
        std::vector<double> master_density;
        std::vector<int> master_valleyIndex;
        std::vector<int> master_peakIndex;
    
        for (int i = 0; i < getProbeSets()->size(); i++)
        {
            CNProbeSet *pset = (*getProbeSets())[i];
            allelicDiffs.push_back(pset->getAllelicDifference());
        }
        findDensityAndPeaks( allelicDiffs,
                             weights,
                             master_xi,
                             master_density,
                             master_valleyIndex,
                             master_peakIndex,
                             CovariateParams::m_masterPeakPointCount,
                             CovariateParams::m_masterPeakBandwidth
                             );

        // remove junk peaks and valleys at the extremities
        trimToRange(master_valleyIndex, master_xi, xiLimit);
        trimToRange(master_peakIndex, master_xi, xiLimit);

        std::vector<double> intermediate;
        if (keepIntermediateData)
        {
            intermediate = allelicDiffs;      // save for later
            intermDensitiesFileInit(iCovariateIndex);
            intermDensitiesOutput( master_xi,
                                   master_density,
                                   master_peakIndex,
                                   "RawMasterDensity",
                                   "RawMasterPeaks"
                                   );
        }

        // nothing else to do if only one peak (data presumed too noisy)
        if (master_peakIndex.size() <= 1)
        {
            if (keepIntermediateData) {
                outputIntermediateDataAD(intermediate, iCovariateIndex);
                intermDensitiesFileClose();
            }
            break;
        }

        // coarse "allele peaks" processing if requested for this covariate
        if (CovariateParams::m_vCoarseAPAdjust[iCovariateIndex] == "on")
        {
            allelicDiffs.clear();
            weights.clear();
            master_xi.clear();
            master_density.clear();
            master_valleyIndex.clear();
            master_peakIndex.clear();

            CNProbeSetArray coarseProbeSets;
            coarseAPProcess(iCovariateIndex);

            for (int i = 0; i < getProbeSets()->size(); i++)
            {
                double tmp = m_vCoarseAllelePeaks[i];
                if (!isNaN(tmp))
                {
                    coarseProbeSets.push_back((*getProbeSets())[i]);
                    allelicDiffs.push_back(tmp);
                }
            }
            std::vector<double> coarseWeights;
            computeWeights(&coarseProbeSets, coarseWeights);
            findDensityAndPeaks( allelicDiffs,
                                 coarseWeights,
                                 master_xi,
                                 master_density,
                                 master_valleyIndex,
                                 master_peakIndex,
                                 CovariateParams::m_masterPeakPointCount,
                                 CovariateParams::m_masterPeakBandwidth
                                 );

            // remove junk peaks and valleys at the extremities
            trimToRange(master_valleyIndex, master_xi, xiLimit);
            trimToRange(master_peakIndex, master_xi, xiLimit);

            if (keepIntermediateData)
            {
                intermDensitiesOutput( master_xi,
                                       master_density,
                                       master_peakIndex,
                                       "MasterDensity",
                                       "MasterPeaks"
                                       );
            }
        }
 
        // gather the covariates for binning...
        std::vector<float> covariates;
        std::vector<std::pair<int, int> > binInfo;
        int iTransIndex = CovariateParams::m_vAPCovariates[iCovariateIndex];
        for (int iIndex = 0; iIndex < getProbeSets()->size(); iIndex++)
        {
            CNProbeSet *pset = (*getProbeSets())[iIndex];
            float covValue = pset->getCovariateValue(iTransIndex);
            if (!isNaN(covValue)) {
                covariates.push_back(covValue);
                binInfo.push_back(std::make_pair(0, iIndex));
            }
        }

        // all covariates are NaN, don't bother
        if (covariates.empty())
        {
            if (keepIntermediateData) {
                intermDensitiesFileClose();
            }
            continue;
        }

        // ...and bin them
        std::vector<int> binning(covariates.size());
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

        for (int iBinIndex = 0; iBinIndex < numBins; iBinIndex++)
        {
            std::vector<int> probeSetsInBinIdx;
            for (
                int iIndex = binBoundaries[iBinIndex];
                iIndex < binBoundaries[iBinIndex + 1];
                iIndex++
            ){
                probeSetsInBinIdx.push_back( binInfo[iIndex].second );
            }

            // in the coarse peaks case this is a subset of probeSetsInBinIdx[]
            // corresponding to the probe sets after the no man's land filtering.
            // Otherwise, it's the same as probeSetsInBinIdx[].
            CNProbeSetArray probeSetsInBin;
            allelicDiffs.clear();
            if (CovariateParams::m_vCoarseAPAdjust[iCovariateIndex] == "on")
            {
                for (int i = 0; i < probeSetsInBinIdx.size(); i++) {
                    double tmp = m_vCoarseAllelePeaks[probeSetsInBinIdx[i]];
                    // use only filtered values (no man's land in allele peaks)
                    if (!isNaN(tmp)) {
                        allelicDiffs.push_back(tmp); 
                        probeSetsInBin.push_back((*getProbeSets())[probeSetsInBinIdx[i]]);
                    }
                }
            } else {
                for (int i = 0; i < probeSetsInBinIdx.size(); i++) {
                    allelicDiffs.push_back( (*getProbeSets())[probeSetsInBinIdx[i]]->getAllelicDifference() );
                    probeSetsInBin.push_back((*getProbeSets())[probeSetsInBinIdx[i]]);
                }
            }

            // bin too small, cannot calculate the bandwidth, don't bother
            if (probeSetsInBin.size() <= 1)
            {
                continue;
            }

            std::vector<double> covBin_weights;
            computeWeights(&probeSetsInBin, covBin_weights);

            // find covariate bin peaks
            std::vector<double> covBin_xi;
            std::vector<double> covBin_density;
            std::vector<int> covBin_valleyIndex;
            std::vector<int> covBin_peakIndex;

            findDensityAndPeaks( allelicDiffs,
                                 covBin_weights,
                                 covBin_xi,
                                 covBin_density,
                                 covBin_valleyIndex,
                                 covBin_peakIndex,
                                 CovariateParams::m_covariatePeakPointCount,
                                 CovariateParams::m_covariatePeakBandwidth
                                 );

            // remove junk peaks and valleys at the extremities
            trimToRange(covBin_valleyIndex, covBin_xi, xiLimit);
            trimToRange(covBin_peakIndex, covBin_xi, xiLimit);

            // the map is: key=cov. bin peak index, value=master peak index
            std::map<int, int> matchedPeaks;

            matchPeaks( master_xi,
                        master_density,
                        master_peakIndex,
                        covBin_xi,
                        covBin_peakIndex,
                        matchedPeaks
                        );

            if (keepIntermediateData) {
                std::string binIndexStr = Convert::toString(iBinIndex);
                intermDensitiesOutput( covBin_xi,
                                       covBin_density,
                                       covBin_peakIndex,
                                       "Bin_" + binIndexStr + "_Density",
                                       "Bin_" + binIndexStr + "_Peaks"
                                       );

                intermMatchedPeaks( matchedPeaks,
                                    covBin_xi,
                                    master_xi,
                                    "Bin_" + binIndexStr + "_MatchedPeaks"
                                    );
            }

            // do not adjust this bin if cov. bin data is too noisy (one peak or less)
            if (matchedPeaks.size() <= 1)
            {
                continue;
            }

            adjustByMatchedPeaks( getProbeSets(), probeSetsInBinIdx, master_xi, covBin_xi, matchedPeaks );

            // mark the fact that ADs were altered
            AD_modulation_no_single_peak_state = true;
        }

        if (keepIntermediateData) {
            outputIntermediateDataAD(intermediate, iCovariateIndex);
            intermDensitiesFileClose();
        }
    }

    // flag the relevant cychp header parameter telling it whether the correction was applied
    getEngine()->setOpt("AD-modulation-no-single-peak", AD_modulation_no_single_peak_state ? "true" : "false");

    m_vCoarseAllelePeaks.clear();                       // set size=0
    std::vector<double>().swap(m_vCoarseAllelePeaks);   // set capacity=0

    Verbose::out(1, "Done");
}

void CNAnalysisMethodAllelicDifferenceCytoScan::trimToRange(
                                                    std::vector<int>& indices,
                                                    const std::vector<double>& xi,
                                                    const double xiLimit
                                                    )
{
    std::vector<int> indices_tmp(indices);
    indices.clear();

    // ignore indices whose xi values lie outside [-xiLimit, +xiLimit]
    for (int i = 0; i < indices_tmp.size(); i++) {
        if (std::fabs(xi[indices_tmp[i]]) <= xiLimit) {
            indices.push_back(indices_tmp[i]);
        }
    }
}

void CNAnalysisMethodAllelicDifferenceCytoScan::findDensityAndPeaks(
                                                    const std::vector<double>& allelicDiffs,
                                                    const std::vector<double>& weights,
                                                    std::vector<double>& xi,
                                                    std::vector<double>& density,
                                                    std::vector<int>& valleyIndex,
                                                    std::vector<int>& peakIndex,
                                                    int numPoints,
                                                    double bandwidth
                                                    )
{
    const bool symmetrise = true;

    double bw;
    if (symmetrise)
    {
        int N = allelicDiffs.size();
        std::vector<double> symAllelicDiffs(allelicDiffs);
        symAllelicDiffs.resize(2*N);
        for (int i = N; i < 2*N; i++)
        {
            symAllelicDiffs[i] = -allelicDiffs[i - N];
        }
        bw = bwnrd(symAllelicDiffs, bandwidth);
    }
    else {
        bw = bwnrd(allelicDiffs, bandwidth);
    }
    fitDensityCurve(allelicDiffs, weights, xi, density, numPoints, bw, symmetrise);

    std::vector<double> sqDensity(density.size());
    for (int i = 0; i < density.size(); i++) {
        sqDensity[i] = density[i] * density[i];
    }
    double overallSignal = trapzoid(xi, sqDensity);

    findpeaks(valleyIndex, peakIndex, density, 0.046*overallSignal, xi);
}

void CNAnalysisMethodAllelicDifferenceCytoScan::computeWeights(
                                                    CNProbeSetArray *psets,
                                                    std::vector<double>& weights
                                                    )
{
    int N = psets->size();
    weights.assign(N, 1.0/N);
}

void CNAnalysisMethodAllelicDifferenceCytoScan::adjustByMatchedPeaks(
                                                        CNProbeSetArray *psets,
                                                        const std::vector<int>& probeSetsInBinIdx,
                                                        const std::vector<double>& master_xi,
                                                        const std::vector<double>& covBin_xi,
                                                        const std::map<int, int>& matchedPeaks
                                                        )
{
    if (matchedPeaks.empty())
    {
        return;
    }

    for (int ii = 0; ii < probeSetsInBinIdx.size(); ii++)
    {
        int iIndex = probeSetsInBinIdx[ii];
        float oldAD = (*psets)[iIndex]->getAllelicDifference();

        // note we rely on std::map keeping keys sorted at all times
        std::map<int, int>::const_iterator it, prev_it;
        for (it = matchedPeaks.begin(); it != matchedPeaks.end(); ++it) {
            double cov_xi = covBin_xi[(*it).first];
            if (oldAD - cov_xi < 0) {
                break;
            }
            prev_it = it;
        }

        float newAD;
        if (it == matchedPeaks.begin())
        {
            // interpolation left of the leftmost peak
            newAD = oldAD * ( master_xi[(*it).second] / covBin_xi[(*it).first] );
        }
        else if (it != matchedPeaks.end())
        {
            // linear interpolation between peaks
            double tt = (oldAD - covBin_xi[(*prev_it).first]) / (covBin_xi[(*it).first] - covBin_xi[(*prev_it).first]);
            newAD = master_xi[(*prev_it).second] + tt * (master_xi[(*it).second] - master_xi[(*prev_it).second]);
        }
        else
        {
            // interpolation right of the rightmost peak
            newAD = oldAD * ( master_xi[(*prev_it).second] / covBin_xi[(*prev_it).first] );
        }

        (*psets)[iIndex]->setAllelicDifference(newAD);
    }
}

void CNAnalysisMethodAllelicDifferenceCytoScan::matchPeaks(
                                                    const std::vector<double>& master_xi,
                                                    const std::vector<double>& master_density,
                                                    const std::vector<int>& master_peakIndex,
                                                    const std::vector<double>& covBin_xi,
                                                    const std::vector<int>& covBin_peakIndex,
                                                    std::map<int, int>& matchedPeaks
                                                    )
{
    const float peakTrim = 2.5;
    const float significantlyCloserRatio = 0.5;

    // Throw away peaks outside the [-peakTrim, +peakTrim] range and
    // insert master peak index together with corresponding density into a pair.
    // This is done to enable sorting master peaks by density (aka. ranking master peaks)
    // which is needed for matching.
    std::vector<std::pair<int, double> > trimmed_master_peakIndex;
    for (int i = 0; i < master_peakIndex.size(); i++) {
        int j = master_peakIndex[i];
        if (master_xi[j] >= -peakTrim && master_xi[j] <= peakTrim) {
            trimmed_master_peakIndex.push_back( std::make_pair(j, master_density[j]) );
        }
    }

    // sort trimmed_master_peakIndex in the density height (rank) order
    std::sort( trimmed_master_peakIndex.begin(),
               trimmed_master_peakIndex.end(),
               comparePeakRank
               );

    // for covariate bin peaks collect just the peak indices inside [-peakTrim, +peakTrim]
    std::vector<int> trimmed_covBin_peakIndex;
    for (int i = 0; i < covBin_peakIndex.size(); i++) {
        int j = covBin_peakIndex[i];
        if (covBin_xi[j] >= -peakTrim && covBin_xi[j] <= peakTrim) {
            trimmed_covBin_peakIndex.push_back(j);
        }
    }

    if (trimmed_covBin_peakIndex.empty())
    {
        return;
    }

    // peak matching
    for (int iIndex = 0; iIndex < trimmed_master_peakIndex.size(); iIndex++) {
        // find cov. bin peak closest to this master peak
        int closestIndex = 0;
        double closestDist = std::fabs( covBin_xi[trimmed_covBin_peakIndex[0]] - master_xi[trimmed_master_peakIndex[iIndex].first] );
        for (int jIndex = 1; jIndex < trimmed_covBin_peakIndex.size(); jIndex++)
        {
            double curDist = std::fabs( covBin_xi[trimmed_covBin_peakIndex[jIndex]] - master_xi[trimmed_master_peakIndex[iIndex].first] );
            if (closestDist > curDist)
            {
                closestDist = curDist;
                closestIndex = jIndex;
            }
        }

        // see if it was already matched to a master peak
        std::map<int, int>::iterator it = matchedPeaks.find( trimmed_covBin_peakIndex[closestIndex] );
        bool matched = it != matchedPeaks.end();
        if (matched) {
            // now see if it is "significantly closer" to current master peak than
            // to the master peak it's been associated with, and reassign if it is indeed "significantly closer"
            double distAssociated = std::fabs( covBin_xi[trimmed_covBin_peakIndex[closestIndex]] - master_xi[(*it).second] );
            double distCurrent  = std::fabs( covBin_xi[trimmed_covBin_peakIndex[closestIndex]] - master_xi[trimmed_master_peakIndex[iIndex].first] );
            if (distCurrent < significantlyCloserRatio * distAssociated)
            {
                matchedPeaks[trimmed_covBin_peakIndex[closestIndex]] = trimmed_master_peakIndex[iIndex].first;
            }
        } else {
            // match it to current master peak
            matchedPeaks[trimmed_covBin_peakIndex[closestIndex]] = trimmed_master_peakIndex[iIndex].first;
        }
    }
}

bool CNAnalysisMethodAllelicDifferenceCytoScan::comparePeakRank(const std::pair<int, double>& a, const std::pair<int, double>& b)
{
    return a.second > b.second;
}

void CNAnalysisMethodAllelicDifferenceCytoScan::outputIntermediateDataAD(const std::vector<double>& intermediate,
                                                                         int iCovariateIndex
                                                                         )
{
    affx::File5_File file5;
    affx::File5_Group* group5;
    affx::File5_Tsv* tsv5 = NULL;

    std::string analysisString = "analysis";
    std::string covarString = ".covar" + Convert::toString(iCovariateIndex);
    std::string strExperimentName = getExperiment()->getExperimentName();
    std::string fileName = Fs::join(getEngine()->getOpt("out-dir"), analysisString,
                                    strExperimentName + covarString + ".adjustedAD.a5");

    Verbose::out(3, "Writing intermediate data for " + fileName);

    try {
        file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
        group5 = file5.openGroup("AllelicDifferences", affx::FILE5_REPLACE);

        tsv5 = group5->openTsv("AllelicDifferences", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "ProbeSetName", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "OldAllelicDifference", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 4, "NewAllelicDifference", affx::FILE5_DTYPE_DOUBLE);

        for (int i = 0; i < getProbeSets()->size(); i++)
        {
            tsv5->set_string(0, 0, (*getProbeSets())[i]->getProbeSetName());
            tsv5->set_i(0, 1, (*getProbeSets())[i]->getChromosome());
            tsv5->set_i(0, 2, (*getProbeSets())[i]->getPosition());
            tsv5->set_d(0, 3, intermediate[i]);
            tsv5->set_d(0, 4, (*getProbeSets())[i]->getAllelicDifference());
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        delete group5;
        file5.close();
    }
    catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump allelic difference data to file."));}
}

void CNAnalysisMethodAllelicDifferenceCytoScan::intermDensitiesFileInit(
                                                    const int iCovariateIndex
                                                    )
{
    std::string analysisString = "analysis";
    std::string covarString = ".covar" + Convert::toString(iCovariateIndex);
    std::string strExperimentName = getExperiment()->getExperimentName();
    std::string fileName = Fs::join(getEngine()->getOpt("out-dir"), analysisString,
                                    strExperimentName + covarString + ".densitiesPeaks.a5");

    Verbose::out(3, "Writing intermediate data for " + fileName);

    try {
        m_file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
        m_group5 = m_file5.openGroup("DensitiesAndPeaks", affx::FILE5_REPLACE);
    }
    catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump density curves data to file."));}
}


void CNAnalysisMethodAllelicDifferenceCytoScan::intermDensitiesOutput(
                                                    const std::vector<double>& xi,
                                                    const std::vector<double>& density,
                                                    const std::vector<int>& peakIndex,
                                                    const std::string& densityName,
                                                    const std::string& peaksName
                                                    )
{
    affx::File5_Tsv* tsv5 = NULL;

    tsv5 = m_group5->openTsv(densityName.c_str(), affx::FILE5_REPLACE);
    tsv5->defineColumn(0, 0, "xi", affx::FILE5_DTYPE_DOUBLE);
    tsv5->defineColumn(0, 1, "DensityValue", affx::FILE5_DTYPE_DOUBLE);

    for (int i = 0; i < xi.size(); i++)
    {
        tsv5->set_d(0, 0, xi[i]);
        tsv5->set_d(0, 1, density[i]);
        tsv5->writeLevel(0);
    }
    tsv5->close();
    delete tsv5;

    tsv5 = m_group5->openTsv(peaksName.c_str(), affx::FILE5_REPLACE);
    tsv5->defineColumn(0, 0, "PeakIndex", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0, 1, "PeakLocation", affx::FILE5_DTYPE_DOUBLE);

    for (int i = 0; i < peakIndex.size(); i++)
    {
        tsv5->set_i(0, 0, peakIndex[i]);
        tsv5->set_d(0, 1, xi[ peakIndex[i] ]);
        tsv5->writeLevel(0);
    }
    tsv5->close();
    delete tsv5;
}


void CNAnalysisMethodAllelicDifferenceCytoScan::intermMatchedPeaks(
                                                    const std::map<int, int>& matchedPeaks,
                                                    const std::vector<double>& covBin_xi,
                                                    const std::vector<double>& master_xi,
                                                    const std::string& matchedName
                                                    )
{
    affx::File5_Tsv* tsv5 = NULL;

    tsv5 = m_group5->openTsv(matchedName.c_str(), affx::FILE5_REPLACE);
    tsv5->defineColumn(0, 0, "MasterPeakIndex", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0, 1, "BinPeakIndex", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0, 2, "MasterPeakLocation", affx::FILE5_DTYPE_DOUBLE);
    tsv5->defineColumn(0, 3, "BinPeakLocation", affx::FILE5_DTYPE_DOUBLE);

    std::map<int, int>::const_iterator it;
    for (it = matchedPeaks.begin(); it != matchedPeaks.end(); ++it)
    {
        tsv5->set_i(0, 0, (*it).second);
        tsv5->set_i(0, 1, (*it).first);
        tsv5->set_d(0, 2, master_xi[ (*it).second ]);
        tsv5->set_d(0, 3, covBin_xi[ (*it).first ]);
        tsv5->writeLevel(0);
    }
    tsv5->close();
    delete tsv5;
}
                                                    

void CNAnalysisMethodAllelicDifferenceCytoScan::intermDensitiesFileClose()
{
    delete m_group5;
    m_file5.close();
}

void CNAnalysisMethodAllelicDifferenceCytoScan::coarseAPProcess(int iCovariateIndex)
{
    // same parameters as for standard allelic differences postprocessing for now
    const float dThreePeakFLD_X = 0.0;
    const float dThreePeakFLD_Y = 0.0;

    const float dFourPeakFLD_X = 0.0;
    const float dFourPeakFLD_Y = 0.0;

    float threePeakShrink_x[] = { 0.0, 15.0, 20.0 };
    float threePeakShrink_y[] = { 0.4, 0.3, 0.2 };

    float fourPeakShrink_x[] = { 0.0, 15.0, 20.0, 25.0 };
    float fourPeakShrink_y[] = { 0.6, 0.5, 0.45, 0.4 };

    m_vThreePeakFLD_X.push_back(dThreePeakFLD_X);
    m_vThreePeakFLD_Y.push_back(dThreePeakFLD_Y);

    m_vFourPeakFLD_X.push_back(dFourPeakFLD_X);
    m_vFourPeakFLD_Y.push_back(dFourPeakFLD_Y);

    for (int i = 0; i < 3; i++) {
        m_vThreePeakShrink_X.push_back(threePeakShrink_x[i]);
        m_vThreePeakShrink_Y.push_back(threePeakShrink_y[i]);
    }
    for (int i = 0; i < 4; i++) {
        m_vFourPeakShrink_X.push_back(fourPeakShrink_x[i]);
        m_vFourPeakShrink_Y.push_back(fourPeakShrink_y[i]);
    }
 
    CNAnalysisMethod::PeakShrinkOverride paramOverride;

    paramOverride.m_iStep_override           = CovariateParams::m_vCoarseAPAdjustStep          [iCovariateIndex];
    paramOverride.m_iWindow_override         = CovariateParams::m_vCoarseAPAdjustWindow        [iCovariateIndex];
    paramOverride.m_iPointCount_override     = CovariateParams::m_vCoarseAPAdjustPointCount    [iCovariateIndex];
    paramOverride.m_dBandwidth_override      = CovariateParams::m_vCoarseAPAdjustBandwidth     [iCovariateIndex];
    paramOverride.m_dCutoff_override         = CovariateParams::m_vCoarseAPAdjustCutoff        [iCovariateIndex];
    paramOverride.m_dCleanThreshold_override = CovariateParams::m_vCoarseAPAdjustCleanthreshold[iCovariateIndex];
    paramOverride.m_iCovariateIndex          = iCovariateIndex;
    paramOverride.m_bSymmetry_override       = true;

    shrinkToPeaks(getProbeSets(), &paramOverride);
}

void CNAnalysisMethodAllelicDifferenceCytoScan::postprocessAllelicDifferences()
{
    resetAllelePeakInitialValues();
    determineLocalProbeSets();
    fillChrBoundsImpl(getProbeSets(), m_chrBounds);
    calculateSNPQC(getProbeSets());

    const float dThreePeakFLD_X = 0.0;
    const float dThreePeakFLD_Y = 0.0;

    const float dFourPeakFLD_X = 0.0;
    const float dFourPeakFLD_Y = 0.0;

    float threePeakShrink_x[] = { 0.0, 15.0, 20.0 };
    float threePeakShrink_y[] = { 0.4, 0.3, 0.2 };

    float fourPeakShrink_x[] = { 0.0, 15.0, 20.0, 25.0 };
    float fourPeakShrink_y[] = { 0.6, 0.5, 0.45, 0.4 };

    m_vThreePeakFLD_X.push_back(dThreePeakFLD_X);
    m_vThreePeakFLD_Y.push_back(dThreePeakFLD_Y);

    m_vFourPeakFLD_X.push_back(dFourPeakFLD_X);
    m_vFourPeakFLD_Y.push_back(dFourPeakFLD_Y);

    for (int i = 0; i < 3; i++) {
        m_vThreePeakShrink_X.push_back(threePeakShrink_x[i]);
        m_vThreePeakShrink_Y.push_back(threePeakShrink_y[i]);
    }
    for (int i = 0; i < 4; i++) {
        m_vFourPeakShrink_X.push_back(fourPeakShrink_x[i]);
        m_vFourPeakShrink_Y.push_back(fourPeakShrink_y[i]);
    }

    shrinkToPeaks(getProbeSets());
}

void CNAnalysisMethodAllelicDifferenceCytoScan::calculateSNPQC(CNProbeSetArray* vProbeSets)
{
    const double log2 = log(2.0);

        // Make a copy of pointers to probe sets sorted by name
        CNProbeSetArray resortedProbeSets = *vProbeSets;
        resortedProbeSets.quickSort(0);

    AffxArray<AffxString> arSNPIds;
    if ((m_pEngine->isOptDefined("snp-qc-snp-list")) && (m_pEngine->getOpt("snp-qc-snp-list") != ""))
    {
      affx::TsvFile tsv;
      tsv.m_optAutoTrim = true;
      tsv.openTable(m_pEngine->getOpt("snp-qc-snp-list"));
      std::string strSnp;
      while (tsv.nextLevel(0) == affx::TSV_OK) {
        tsv.get(0,0,strSnp);
        arSNPIds.add(new AffxString(strSnp));
      }
      tsv.close();

    }
    arSNPIds.quickSort(0);

    std::vector<double> vData;

        AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
        affx::File5_File file5;
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;

        file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
        group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        tsv5 = group5->openTsv("CN5.snp-posteriors", affx::FILE5_OPEN);

        while (tsv5->nextLine() == affx::FILE5_OK) {
            string probeSetName;
            double clusterAAmean;
            double clusterABmean;
            double clusterBBmean;

            tsv5->get(0, 0, &probeSetName);
            tsv5->get(0, 2, &clusterAAmean);
            tsv5->get(0, 6, &clusterABmean);
            tsv5->get(0, 10, &clusterBBmean);

            if (!arSNPIds.empty()) {
                AffxString str(probeSetName);
                if (arSNPIds.binarySearch(str, 0) == -1) {
                    continue;
                }
            }
            CNProbeSet searchProbeSet;
            searchProbeSet.setProbeSetName(probeSetName);
            pair<CNProbeSetArray::iterator, CNProbeSetArray::iterator> range = equal_range(
                                                                                    resortedProbeSets.begin(),
                                                                                    resortedProbeSets.end(),
                                                                                    &searchProbeSet,
                                                                                    NameCompare()
                                                                                    );
            if (range.first == range.second) {
                continue;
            }

            CNProbeSet* foundProbeSet = *range.first;
            if (foundProbeSet->processAsVisualization())
            {
                double log2ratio = log(foundProbeSet->getAAlleleSignal()/foundProbeSet->getBAlleleSignal())/log2;
                double scaledLog2Ratio = log2ratio - clusterABmean;
                scaledLog2Ratio /= 0.5*(abs(clusterAAmean - clusterABmean) + abs(clusterBBmean - clusterABmean));
                vData.push_back(scaledLog2Ratio);
            }
        }
        tsv5->close();
        delete tsv5;

        group5->close();
        delete group5;
        file5.close();

    arSNPIds.deleteAll();

    for (int i = 0; (i < vData.size()); i++)
    {
        if (vData[i] != vData[i])
        {
            vData[i] = 0;
        }
    }

    m_pobjExperiment->setSNPQC(HomHiLoCelListener::computeSNPQC(vData));
    // If this is run now, it shouldn't be rerun in CNAnalysisEngine::postProcessing()
    m_pobjExperiment->setIsSNPQCset(true);
}

void CNAnalysisMethodAllelicDifferenceCytoScan::ADOutlierTrim(double& ad)
{
    if (ad > m_dAllelicDifferenceOutlierTrim) {
        ad = m_dAllelicDifferenceOutlierTrim;
    }
    else if (ad < -m_dAllelicDifferenceOutlierTrim) {
        ad = -m_dAllelicDifferenceOutlierTrim;
    }
}

int CNAnalysisMethodAllelicDifferenceCytoScan::determineBins(vector<float>& covValues, vector<int>& binning, int covariateIndex)
{
    if (CovariateParams::useEquallyPopulatedAPBins(covariateIndex))
    {   
        int numBins = CovariateParams::getNumAPBins(covariateIndex);
        CNAnalysisMethod::binEqualNumber(covValues, binning, numBins);
        return numBins;
    }       
    else if (CovariateParams::useEquallySpacedAPBins(covariateIndex))
    {       
        int numBins = CovariateParams::getNumAPBins(covariateIndex);
        CNAnalysisMethod::binEqualSpacing(covValues, binning, numBins);
        return numBins;
    }           
    else if (CovariateParams::isAPCovariateDiscrete(covariateIndex))
    {       
        return CNAnalysisMethod::covariateIsBinAssignment(covValues, binning);
    }           
    else // unknown binning
    {       
        // This is a late error.
        Err::errAbort("Unknown allelic difference covariate binning type");
    }
    return 0;
}
