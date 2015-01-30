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
 * FILE GenotypeCaller.h
 */

#ifndef _GENOTYPECALLER_V1_H
#define _GENOTYPECALLER_V1_H

//
#include "birdseed-dev/Clusters.h"
#include "birdseed-dev/ClustersReader.h"
#include "birdseed-dev/Matrix.h"
#include "birdseed-dev/Prior.h"
#include "birdseed-dev/PriorsReader.h"
#include "birdseed-dev/birdseeddefs.h"
#include "broadutil/BroadException.h"
//
#include <ostream>
#include <vector>
//

namespace birdseed 
{
namespace v1
{
    using namespace birdseed::dev;
    
    // Result for a single sample-SNP combination
    struct CallAndConfidence
    {
        // For haploid, the call will be either AA or BB.  Only the first two positions of pvalues will be set to something meaningful.
        
        // Counter-intuitively, BB=0, and AA=2, despite the positions of AA & BB in the priors and other data structures.
        enum Call {BB = 0, AB = 1, AA = 2, NOCALL = 3};
        Call call;
        double confidence;
        double pvalues[MAX_NUM_CLUSTERS];
        bool isDiploid;
        
        
        CallAndConfidence(Call call, double confidence, const double pxi_zj[MAX_NUM_CLUSTERS], bool isDiploid):
            call(call),
            confidence(confidence),
            pvalues(),
            isDiploid(isDiploid)
        {
            for (size_t i = 0; i < MAX_NUM_CLUSTERS; ++i) {
                pvalues[i] = pxi_zj[i];
            }
        }
    };

    class GenotypeCaller
    {
      public:
        virtual ~GenotypeCaller(){}

        // Iterate through the call for each sample, in the order they appear in the IntensityMatrix.
        virtual CallAndConfidence nextCall() = 0;
    };
    
    
    // Calls all the samples for a single SNP.
    class BasicGenotypeCaller: public GenotypeCaller
    {
      private:
        std::vector<CallAndConfidence> calls;
        size_t index;
    
      public:
        // Do the calling for a SNP.
        BasicGenotypeCaller(const IntensityMatrix &intensities,
                            const Priors &priors,
                            double correctionFactor,
                            const std::string &priorName,
                            int verbosity = 0,
                            std::ostream *clusterOstrm = NULL);

        // Iterate through the call for each sample, in the order they appear in the IntensityMatrix.
        CallAndConfidence nextCall()
        {
            return calls[index++];
        }
    };

    // Returns NO CALL for every sample
    class DummyGenotypeCaller: public GenotypeCaller
    {
      private:
        static double pvalues[MAX_NUM_CLUSTERS];
        bool isDiploid;
      public:
        DummyGenotypeCaller(bool isDiploid):
            isDiploid(isDiploid)
        {}
        
        CallAndConfidence nextCall();
    };
    
    // Make calls based on clusters that are specified rather than calculated
    class ForcedClusterGenotypeCaller: public GenotypeCaller
    {
      private:
        std::vector<CallAndConfidence> calls;
        size_t index;
    
      public:
        // Do the calling for a SNP.
        ForcedClusterGenotypeCaller(const IntensityMatrix &intensities,
                                    const Clusters &clusters,
                                    double correctionFactor,
                                    const std::string &priorName,
                                    int verbosity = 0,
                                    std::ostream *clusterOstrm = NULL);

        // Iterate through the call for each sample, in the order they appear in the IntensityMatrix.
        CallAndConfidence nextCall()
        {
            return calls[index++];
        }
    };

    template<class PRIORS, class PRIORSREADER, class GENOTYPECALLER>
        class GenderAwareGenotypeCallerTemplate: public GenotypeCaller
    {
      public:
        GenderAwareGenotypeCallerTemplate(const IntensityMatrix &intensities,
                                  const std::vector<Gender> &genders,
                                          PRIORSREADER *priorsReader,
                                  const char *snpName,
                                  double correctionFactor,
                                          int verbosity=0,
                                          std::ostream *clusterOstrm = NULL);
        
        // Iterate through the call for each sample, in the order they appear in the IntensityMatrix.
        CallAndConfidence nextCall();
      private:
        std::vector<Gender> genders;
        std::auto_ptr<GenotypeCaller> femaleCaller;
        std::auto_ptr<GenotypeCaller> maleCaller;
        std::vector<Gender>::const_iterator genderIt;
    };

    
    template<class PRIORS, class PRIORSREADER, class GENOTYPECALLER>
    inline GenderAwareGenotypeCallerTemplate<PRIORS, PRIORSREADER, GENOTYPECALLER>::
    GenderAwareGenotypeCallerTemplate(const IntensityMatrix &intensities,
                                      const std::vector<Gender> &genders,
                                      PRIORSREADER *priorsReader,
                                      const char *snpName,
                                      double correctionFactor,
                                      int verbosity,
                                      std::ostream *clusterOstrm):
        genders(genders),
        genderIt(this->genders.begin())
    {
        const PRIORS *priors;
        if (!priorsReader->ploidyDiffersByGender(snpName) || genders.empty()) {
            priors = priorsReader->getPriorsForSNP(snpName);
            if (priors != NULL) {
                std::string priorName;
                priorName = priorsReader->getPriorName(snpName);
                femaleCaller.reset(new GENOTYPECALLER(intensities, *priors, correctionFactor, priorName,
                                                      verbosity, clusterOstrm));
            } else {
                femaleCaller.reset(new DummyGenotypeCaller(true));
            }
        } else {
            if (intensities.numRows() != genders.size()) {
                throw BroadException("Mismatch between intensities.numRows() and genders.size()", __FILE__, __LINE__);
            }

            // Reserve the max # of elements in intensity matrices to avoid extra mallocing & copying
            IntensityMatrix maleIntensities(0, intensities.numRows());
            IntensityMatrix femaleIntensities(0, intensities.numRows());

            for (size_t i = 0; i < intensities.numRows(); ++i) {
                if (genders[i] == GENDER_MALE) {
                    maleIntensities.push_back(intensities[i]);
                } else {
                    assert(genders[i] == GENDER_FEMALE || genders[i] == GENDER_UNKNOWN);
                    femaleIntensities.push_back(intensities[i]);
                }
            }
            if (femaleIntensities.numRows() > 0) {
                priors = priorsReader->getPriorsForSNP(snpName, GENDER_FEMALE);
                if (priors != NULL) {
                    std::string priorName;
                    priorName = priorsReader->getPriorName(snpName, GENDER_FEMALE);
                    femaleCaller.reset(new GENOTYPECALLER(femaleIntensities,
                                                          *priors,
                                                          correctionFactor,
                                                          priorName,
                                                          verbosity,
                                                          clusterOstrm));
                } else {
                    femaleCaller.reset(new DummyGenotypeCaller(true));
                }
            }
            if (maleIntensities.numRows() > 0) {
                priors = priorsReader->getPriorsForSNP(snpName, GENDER_MALE);
                if (priors != NULL) {
                    std::string priorName;
                    priorName = priorsReader->getPriorName(snpName, GENDER_MALE);
                    maleCaller.reset(new GENOTYPECALLER(maleIntensities,
                                                        *priors,
                                                        correctionFactor,
                                                        priorName,
                                                        verbosity,
                                                        clusterOstrm));
                } else {
                    maleCaller.reset(new DummyGenotypeCaller(false));
                }
            }
        }
    }

    template<class PRIORS, class PRIORSREADER, class GENOTYPECALLER>
    inline
    CallAndConfidence GenderAwareGenotypeCallerTemplate<PRIORS, PRIORSREADER, GENOTYPECALLER>::nextCall()
    {
        if (maleCaller.get() == NULL) {
            return femaleCaller->nextCall();
        } else if (*genderIt++ == GENDER_MALE) {
            return maleCaller->nextCall();
        } else {
            return femaleCaller->nextCall();
        }
    }

    typedef GenderAwareGenotypeCallerTemplate<Priors, PriorsReader, BasicGenotypeCaller> GenderAwareGenotypeCaller;
    typedef GenderAwareGenotypeCallerTemplate<Clusters, ClustersReader, ForcedClusterGenotypeCaller>
        GenderAwareForcedClusterGenotypeCaller;
};
};


#endif /* _GENOTYPECALLER_V1_H */

/******************************************************************/
/**************************[END OF GenotypeCaller.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
