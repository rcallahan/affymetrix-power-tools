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

/*
 * FILE FitSNPGaussiansPriors3.h
 */

#ifndef _FITSNPGAUSSIANSPRIORS3_V1_H_
#define _FITSNPGAUSSIANSPRIORS3_V1_H_

#include "birdseed-dev/Clusters.h"
#include "birdseed-dev/Matrix.h"
#include "birdseed-dev/Prior.h"
#include "broadutil/APTUtil.h"
#include "chipstream/SelfDoc.h"
#include "util/PgOptions.h"
//


namespace birdseed
{
    namespace v1
    {
        using namespace birdseed::dev;

        // 0 = quiet; 1 = some message; 2 = lots of messages
        void setFitSNPGaussiansVerbosity(int verbose);

        // For getting descriptions of parameters that can be controlled via command-line options
        PgOpt **getAlgorithmicOpts();

        // For setting parameters from command-line options
        void setAlgorithmicParameters(CONSTHACK PgOptions &opts);

        // For getting parameter descriptions for chipstream components.  These are appended to
        // opts vector.
        void getSelfDocOptions(std::vector<SelfDoc::Opt> *opts);

        // For setting parameters from chipstream argument.
        void setSelfDocOptions(SelfDoc *doc, CONSTHACK std::map<std::string,std::string> &param);
    
        // Do the clustering
        Clusters FitSNPGaussianPriors3(const IntensityMatrix &intensities, const Priors &priors);

        // Generate the matrix of probabilities of a sample being in a particular cluster (N samples x K clusters)
        void calculatePXI_ZJ(PXI_ZJMatrix *pxi_zj,
                            const birdseed::dev::Clusters &clusters,
                            const IntensityMatrix &intensities, size_t k);

        double calculateSNPSpecificCorrectionFactor(const IntensityMatrix &intensities,
                                                    const Priors &priors);

        // Amount to weight confidence factor determined by comparing most likely cluster to next most likely,
        // as compared to confidence factor determined by distance of sample from cluster mean.
        extern double relative_distance_confidence_weight;

        // After calculating clusters, ensure all weights are >= this
        extern double final_weight_min;
    
    
        // Calculate confidence factor based on distance of sample from cluster mean
        double calculateConfidenceRelativeToCluster(const Matrix2x2 &B,
                                                    const IntensityMatrix::ROW &intensities,
                                                    const IntensityMatrix::ROW &clusterMean);

        // Generate the B matrix that gets passed to calculateConfidenceRelativeToCluster
        Matrix2x2 makeBMatrix(const Matrix2x2::ROW &clusterVars, double covar);
    
    };
};


#endif /* _FITSNPGAUSSIANSPRIORS3_V1_H_ */

/******************************************************************/
/**************************[END OF FitSNPGaussiansPriors3.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
