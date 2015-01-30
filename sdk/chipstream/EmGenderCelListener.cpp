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

#include "chipstream/EmGenderCelListener.h"
//
#include "chipstream/GenoUtility.h"
#include "chipstream/apt-geno-qc/GenoQC.h"
//
#include "stats/stats.h"
//

using namespace std;
using namespace affx;
using namespace affymetrix_fusion_io;

/** 
 * Constructor.
 * @param chrXpSets - Vector of chrX probesets that are not in the
 * pseudo autosomal region and should be haploid in males and
 * diploid in females.
 */
EmGenderCelListener::EmGenderCelListener(std::vector<ProbeListPacked> &chrXpSets,
                                         double k ,
                                         double emThresh ,
                                         double emCutOff ,
                                         double genderCutOff  
                                         ) : m_ChrXProbeSets(chrXpSets),
                                             m_K(k),
                                             m_EmThresh(emThresh),
                                             m_EmCutOff(emCutOff),
                                             m_GenderCutOff(genderCutOff)
{
  m_GenderName="em-cluster-chrX-het-contrast";
  m_GenderDescription="EM clustering of chrX contrast values from SNP probesets";

  declareMetric(m_GenderName + "_gender",ChipSummary::Metric::String);
  declareMetric(m_GenderName + "_gender_chrX_het_rate",ChipSummary::Metric::Double);
}

void EmGenderCelListener::CallGender(const std::vector<double>& contrast,
                                     affx::Gender &gender,
                                     double &hcr) {
  CallGender(contrast, 0.05f, 0.5f, 0.1f, gender, hcr);
}

/** 
 * Loop through the probesets provided and calculate a contrast value
 * for each one using the median of PM probes for A allele and B
 * allele.
 * @param cel - Cel File to get data from.
 * @param chrXProbeSets - Probesets to process (should be chrX non-pseudo autosomal).
 * @param contrastValues - Contrast values vector to be filled in.
 */
void EmGenderCelListener::fillInContrastValues(FusionCELData *cel, 
                                               std::vector<ProbeListPacked> &chrXProbeSets, 
                                               std::vector<double> &contrastValues, 
                                               double k) {
  contrastValues.clear();
  contrastValues.reserve(chrXProbeSets.size());
  for(int psIx = 0; psIx < chrXProbeSets.size(); psIx++) {
    const ProbeSet *ps =  ProbeListFactory::asProbeSet(chrXProbeSets[psIx]);
    contrastValues.push_back(CalculateEMContrast(cel, ps, k));
    delete ps;
  }
}

/**
 * Calculate contrast values
 */
double EmGenderCelListener::CalculateEMContrast(FusionCELData* cel, const ProbeSet* ps, const double k) {
    float Amedian = -1, Bmedian = -1;
    bool medianOk = alleleMedians(cel, ps, Amedian, Bmedian);
    if(!medianOk || Amedian <= 0 || Bmedian <= 0) {
      Err::errAbort("EmGenderCelListener::fillInContrastValues() - Warning. alleleMedian() failed for probeset: " + ToStr(ps->name));
    }
    double contrast = -2, strength = -2;
    ContrastCentersStretch(k, Amedian, Bmedian, contrast, strength);
    if(fabs(contrast) > 1) {
      Err::errAbort("EmGenderCelListener::fillInContrastValues() - Can't have abs(contrast) > 1 for probeset: " + ToStr(ps->name));
    }
    return contrast;
}

/** 
 * Process another cel files worth of data.
 * @param cel - Filehandle to open cel file.
 */
void EmGenderCelListener::newChip(FusionCELData *cel) {
  double hcr;
  std::vector<ChipSummary::Metric> metrics;
  m_SummaryStats.push_back(metrics);
  m_CelNames.push_back(cel->GetFileName());

  vector<double> contrastValues;
  fillInContrastValues(cel, m_ChrXProbeSets, contrastValues, m_K);
  // use the default thresholds.
  Gender gender = UnknownGender;
  CallGender(contrastValues,m_EmThresh,m_EmCutOff,m_GenderCutOff,gender,hcr);
  m_Genders.push_back(gender);

  m_SummaryStats[m_CelNames.size()-1].push_back(ChipSummary::Metric(m_GenderName + "_gender",affx::getGenderString(gender)));
  m_SummaryStats[m_CelNames.size()-1].push_back(ChipSummary::Metric(m_GenderName + "_gender_chrX_het_rate",hcr));
  setValid(true);
}

/**
* @brief EM-based gender calling algorithm, call gender on a single sample
*/
void EmGenderCelListener::CallGender(const std::vector<double>& contrast,
                                     const double em_thresh,
                                     const double em_cutoff,
                                     const double gender_cutoff,
                                     affx::Gender &gender,
                                     double &hcr) {
    CEMSeed seed;
    seed.setMu(-0.66f, 0.0f, 0.66f);
    seed.setSigma(0.1f, 0.1f, 0.1f);
    seed.setWeight(0.33f, 0.34f, 0.33f);
    seed.setMinMu(-2.0f, -.05f, .25f);
    seed.setMinSigma(0.02f, 0.02f, 0.02f);
    seed.setMaxMu(-.25f, .05f, 2.0f);
    seed.setMaxSigma(.3f, .3f, .3f);

    CPrimeEM em;
    em.setData(contrast);
    em.setThreshold(em_thresh);
    em.EMEstimate(seed);
    CEMEst* pEst = em.getEMEstimates();
    int nAA = count_if(pEst->m_class.begin(),pEst->m_class.end(),bind2nd(equal_to<int>(),0));
    int nBB = count_if(pEst->m_class.begin(),pEst->m_class.end(),bind2nd(equal_to<int>(),2));
    int nAB = count_if(pEst->m_class.begin(),pEst->m_class.end(),bind2nd(equal_to<int>(),1));
    hcr=-1;
    if(nAA+nBB+nBB > 0 && float(nAA+nAB+nBB)/pEst->m_class.size() > em_cutoff) {
        hcr = float(nAB)/float(nAA+nAB+nBB);
        if(hcr < gender_cutoff)
            gender = affx::Male;
        else
            gender = affx::Female;
    } else
        gender = affx::UnknownGender;
}
