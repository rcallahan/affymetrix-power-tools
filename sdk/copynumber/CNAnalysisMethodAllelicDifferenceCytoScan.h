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

#ifndef _CNAnalysisMethodAllelicDifferenceCytoScan_H_
#define _CNAnalysisMethodAllelicDifferenceCytoScan_H_
/**
 * @file CNAnalysisMethodAllelicDifferenceCytoScan.h
 *
 * @brief This header contains the CNAnalysisMethodAllelicDifferenceCytoScan class definition.
 */

#include "copynumber/CNAnalysisMethodAllelicDifference.h"
#include "file5/File5.h"

/**
 * @brief  The AllelicDifference analysis method - the CytoScan variation.
 *
 */
class CNAnalysisMethodAllelicDifferenceCytoScan : public CNAnalysisMethodAllelicDifference
{
public:
    static std::string getType() {return "allelic-difference-CytoScan";}
    static std::string getDescription() {return "Copynumber AllelicDifference CytoScan";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    virtual ~CNAnalysisMethodAllelicDifferenceCytoScan() {}

    virtual void run();

    int determineBins(vector<float>& covValues, vector<int>& binning, int covariateIndex);

    void postprocessAllelicDifferences();
    void calculateSNPQC(CNProbeSetArray* vProbeSets);
    void removeAllelicDifferencesModulation();
    void coarseAPProcess(int iCovariateIndex);

protected:
    virtual void ADOutlierTrim(double& ad);

    void findDensityAndPeaks(
                          const std::vector<double>& allelicDiffs_ptr,
                          const std::vector<double>& weights,
                          std::vector<double>& xi,
                          std::vector<double>& density,
                          std::vector<int>& valleyIndex,
                          std::vector<int>& peakIndex,
                          int numPoints,
                          double bandwidth
                          );

    void computeWeights( CNProbeSetArray *psets,
                         std::vector<double>& weights
                         );

    void trimToRange(std::vector<int>& indices, const std::vector<double>& xi, const double xiLimit);

    void outputIntermediateDataAD(const std::vector<double>& intermediate, int iCovariateIndex);

    void outputIntermediateDensityPeaks( const std::vector<double>& xi,
                                         const std::vector<double>& density,
                                         const std::vector<int>& peakIndex,
                                         const int iCovariateIndex,
                                         const int iBinIndex
                                         );

    void intermDensitiesFileInit(const int iCovariateIndex);

    void intermDensitiesOutput( const std::vector<double>& xi,
                                const std::vector<double>& density,
                                const std::vector<int>& peakIndex,
                                const std::string& densityName,
                                const std::string& peaksName
                                );

    void intermMatchedPeaks( const std::map<int, int>& matchedPeaks,
                             const std::vector<double>& covBin_xi,
                             const std::vector<double>& master_xi,
                             const std::string& matchedName
                             );

    void intermDensitiesFileClose();

    static void adjustByMatchedPeaks( CNProbeSetArray *psets,
                                      const std::vector<int>& probeSetsInBinIdx,
                                      const std::vector<double>& master_xi,
                                      const std::vector<double>& covBin_xi,
                                      const std::map<int, int>& matchedPeaks
                                      );

    static void matchPeaks( const std::vector<double>& master_xi,
                            const std::vector<double>& master_density,
                            const std::vector<int>& master_peakIndex,
                            const std::vector<double>& covBin_xi,
                            const std::vector<int>& covBin_peakIndex,
                            std::map<int, int>& matchedPeaks
                            );

    static bool comparePeakRank(const std::pair<int, double>& a, const std::pair<int, double>& b);

private:
    affx::File5_File m_file5;
    affx::File5_Group* m_group5;
};

struct NameCompare : std::binary_function<const CNProbeSet*, const CNProbeSet*, bool>
{
    bool operator()(const CNProbeSet* lhs, const CNProbeSet* rhs) const {
        return lhs->getProbeSetName() < rhs->getProbeSetName();
    }
};

#endif
