////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

//
#include "birdseed-dev/GenotypeCaller.h"
//
#include "birdseed-dev/FitSNPGaussiansPriors3.h"
//
#include "broadutil/BroadException.h"
//
#include <algorithm>
#include <cstring>
#include <iostream>
#include <string.h>
//

using namespace std;
using namespace birdseed::dev;

static const double MIN_VARIANCE = 1e-14;

static void floorClusterWeights(Clusters *clusters, int verbosity)
{
    // Make sure all weights are above minimum
    for (size_t i = 0; i < clusters->weights.size(); ++i) {
        clusters->weights[i] = max(clusters->weights[i], final_weight_min);
    }
    if (verbosity >= 2) {
        cout << "After weight min: " << clusters->tostring() << "\n";
    }
}


static CallAndConfidence doCall(const PXI_ZJMatrix::ROW &pxi_zj,
                                bool isDiploid,
                                const IntensityMatrix::ROW &intensities,
                                const Clusters &clusters,
                                const Matrix2x2 *B,
                                const string &priorName)
{
    int biggest = -1;
    int nextBiggest = -1;
    double biggestVal = -INFINITY;
    double nextBiggestVal = -INFINITY;
    for (size_t i = 0; i < pxi_zj.size(); ++i) {
        if (pxi_zj[i] > biggestVal) {
            nextBiggest = biggest;
            nextBiggestVal = biggestVal;
            biggest = i;
            biggestVal = pxi_zj[i];
        } else if (pxi_zj[i] > nextBiggestVal) {
            nextBiggest = i;
            nextBiggestVal = pxi_zj[i];
        }
    }

    double pvalues[MAX_NUM_CLUSTERS];
    memset(pvalues, 0, sizeof(pvalues));

    //
    if ((biggest==-1)||(nextBiggest==-1)) {
      // dump pxi_zj
//      printf("pxi_zj[%d]=[",(int)index);
//      for (size_t i = 0; i < pxi_zj.size(); ++i) {
//        printf("%s%f",((i!=0)?",":""),pxi_zj[i]);
//      }
//      printf("]\n");
      //
      return CallAndConfidence(CallAndConfidence::NOCALL, 0.0, pvalues, isDiploid);
    }

    for (size_t i = 0; i < pxi_zj.size(); ++i) {
        // apt-probeset-genotype want a higher value to mean lower likelihood
        pvalues[i] = 1 - pxi_zj[i];
    }

    if (pxi_zj[biggest] <= 0) {
        cerr << "WARNING: pxi_zj[biggest] is not positive: " << pxi_zj[biggest] << " for " << priorName << endl;
    }
    
    double confidence;
    if (pxi_zj[biggest] > 0) {
        confidence = pxi_zj[nextBiggest] / pxi_zj[biggest];
        double confidenceRelativeToCluster =
            calculateConfidenceRelativeToCluster(B[biggest],
                                                 intensities,
                                                 clusters.means[biggest]);
        confidence = confidence * relative_distance_confidence_weight +
            confidenceRelativeToCluster * (1.0 - relative_distance_confidence_weight);
    } else {
        confidence = 1.0;
    }

    if (isDiploid) {
        switch (biggest) {
        case Priors::AA_INDEX:
            return CallAndConfidence(CallAndConfidence::AA, confidence, pvalues, isDiploid);
        case Priors::AB_INDEX:
            return CallAndConfidence(CallAndConfidence::AB, confidence, pvalues, isDiploid);
        case Priors::BB_INDEX:
            return CallAndConfidence(CallAndConfidence::BB, confidence, pvalues, isDiploid);
        default:
            assert(false);
            return CallAndConfidence(CallAndConfidence::NOCALL, 0.0, pvalues, isDiploid);
        }
    } else {
        switch (biggest) {
        case Priors::A_INDEX:
            return CallAndConfidence(CallAndConfidence::AA, confidence, pvalues, isDiploid);
        case Priors::B_INDEX:
            return CallAndConfidence(CallAndConfidence::BB, confidence, pvalues, isDiploid);
        default:
            assert(false);
            return CallAndConfidence(CallAndConfidence::NOCALL, 0.0, pvalues, isDiploid);
        }
    }
}

static void doCalls(vector<CallAndConfidence> *calls,
                    const Clusters &clusters,
                    const IntensityMatrix &intensities,
                    bool isDiploid,
                    const string &priorName,
                    size_t verbosity)
{
    size_t numClusters = (isDiploid? MAX_NUM_CLUSTERS: MAX_NUM_CLUSTERS - 1);
    PXI_ZJMatrix pxi_zj(intensities.numRows(), numClusters);
    calculatePXI_ZJ(&pxi_zj, clusters, intensities, numClusters);

    Matrix2x2 B[MAX_NUM_CLUSTERS];
    for (size_t i = 0; i < numClusters; ++i) {
        B[i] = makeBMatrix(clusters.vars[i], clusters.covar);
    }

    calls->reserve(intensities.numRows());
    for (size_t i = 0; i < intensities.numRows(); ++i) {
        calls->push_back(doCall(pxi_zj[i], isDiploid, intensities[i], clusters, B, priorName));
    }
}


static void adjustPriors(Priors *priors, const IntensityMatrix &intensities, double correctionFactor)
{
    double squared = correctionFactor * correctionFactor;

    bool seenNonZeroObservations = false;
    for (size_t i = 0; i < priors->getNumPriors(); ++i) {
        priors->getPrior(i).m_mean = priors->getPrior(i).m_mean * correctionFactor;
        priors->getPrior(i).m_covarMatrix = priors->getPrior(i).m_covarMatrix * squared;
        priors->getPrior(i).m_covarMatrix[0][0] =
            max(priors->getPrior(i).m_covarMatrix[0][0], MIN_VARIANCE);
        priors->getPrior(i).m_covarMatrix[1][1] =
            max(priors->getPrior(i).m_covarMatrix[1][1], MIN_VARIANCE);

        if (priors->getPrior(i).m_numObservations > 0) {
            seenNonZeroObservations = true;
        }
    }
    if (!seenNonZeroObservations) {
        for (size_t i = 0; i < priors->getNumPriors(); ++i) {
            // Don't allow numObservations to be zero, to avoid get NaNs from FP arithmetic
            priors->getPrior(i).m_numObservations = 1;
        }
    }
}

BasicGenotypeCaller::BasicGenotypeCaller(const IntensityMatrix &intensities,
                                         const Priors &priors,
                                         double correctionFactor,
                                         const string &priorName,
                                         int verbosity,
                                         std::ostream *clusterOstrm):
    calls(),
    index(0)
{
    if (intensities.numRows() < MIN_SAMPLES_TO_CLUSTER) {
        std::stringstream strm;
        strm << "Not enough samples to clusters.  " << MIN_SAMPLES_TO_CLUSTER << " are needed but there are only " <<
            intensities.numRows() << ".";
        throw BroadException(strm.str().c_str(), __FILE__, __LINE__);
    }
    
    setFitSNPGaussiansVerbosity(verbosity);
    if (correctionFactor < 0) {
        correctionFactor = calculateSNPSpecificCorrectionFactor(intensities, priors);
        if (verbosity >= 2) {
            cout.precision(10);
            cout << "SNP-specific correction factor for " << priorName << "=" << correctionFactor << endl;
        }
    }
    Priors adjustedPriors(priors);
    adjustPriors(&adjustedPriors, intensities, correctionFactor);
    if (verbosity >= 2) {
        cout << "Priors\n";
        for (size_t i = 0; i < adjustedPriors.getNumPriors(); ++i) {
            const Prior &prior = adjustedPriors.getPrior(i);
            cout << prior.m_mean[0] << " " << prior.m_mean[1] << " " << prior.m_covarMatrix[0][0] << " "
                 << " " << prior.m_covarMatrix[1][0] << " " << prior.m_covarMatrix[1][1] << " "
                 << prior.m_numObservations << endl;
        }
    }
    
    Clusters clusters = FitSNPGaussianPriors3(intensities, adjustedPriors);
    if (clusterOstrm != NULL) {
        (*clusterOstrm) << priorName << clusters.standardClusterString() << endl;
    }
    floorClusterWeights(&clusters, verbosity);
    doCalls(&calls, clusters, intensities, priors.isDiploid(), priorName, verbosity);
}

double DummyGenotypeCaller::pvalues[MAX_NUM_CLUSTERS] = {0.0, 0.0, 0.0};

CallAndConfidence DummyGenotypeCaller::nextCall()
{
    return CallAndConfidence(CallAndConfidence::NOCALL, 1.0, pvalues, isDiploid);
}

ForcedClusterGenotypeCaller::ForcedClusterGenotypeCaller(const IntensityMatrix &intensities,
                                                         const Clusters &clusters,
                                                         double correctionFactor,
                                                         const std::string &priorName,
                                                         int verbosity,
                                                         std::ostream *clusterOstrm):
    calls(),
    index(0)
{
    bool isDiploid = clusters.k == MAX_NUM_CLUSTERS;
    Clusters flooredClusters(clusters);
    if (clusterOstrm != NULL) {
        (*clusterOstrm) << priorName << clusters.standardClusterString() << endl;
    }
    floorClusterWeights(&flooredClusters, verbosity);
    doCalls(&calls, flooredClusters, intensities, isDiploid, priorName, verbosity);
}

/******************************************************************/
/**************************[END OF GenotypeCaller.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
