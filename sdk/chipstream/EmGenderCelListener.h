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

#ifndef EMGENDERCELLISTENER_H
#define EMGENDERCELLISTENER_H

//
#include "chipstream/BioTypes.h"
#include "chipstream/GenderCelListener.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/ProbeSet.h"
//
#include <cstring>
#include <string>
#include <utility>
#include <vector>
//

/**
 * @brief Class for determining gender from a cel set of cel files.
 */
class EmGenderCelListener : public GenderCelListener {

public:

  /** 
   * Constructor.
   * @param chrXpSets - Vector of chrX probesets that are not in the
   * pseudo autosomal region and should be haploid in males and
   * diploid in females.
   */
  EmGenderCelListener(std::vector<ProbeListPacked> &chrXpSets,
                      double k = 4.0f,
                      double emThresh = 0.05f,
                      double emCutOff = 0.5f,
                      double genderCutOff = 0.1f 
                      );
  
  /**
   * Virtual destructor.
   */
  virtual ~EmGenderCelListener(){}

  /** 
   * Process another cel files worth of data.
   * @param cel - Filehandle to open cel file.
   */
  void newChip(affymetrix_fusion_io::FusionCELData *cel);

  /** 
   * Get the genders for the cel files that have been seen.
   * @return - vector of gender calls, one for each cel file in order seen.
   */
  std::vector<affx::Gender> getGenders() { return m_Genders; }
  
  /** 
   * Get the names for the cel files that have been seen.
   * @return - vector of cel file names, one for each cel file in order seen.
   */
  std::vector<std::string> getCelNames() { return m_CelNames;  }

  /** 
   * Loop through the probesets provided and calculate a contrast value
   * for each one using the median of PM probes for A allele and B
   * allele. (CCS Space)
   * @param cel - Cel File to get data from.
   * @param chrXProbeSets - Probesets to process (should be chrX non-pseudo autosomal).
   * @param contrastValues - Contrast values vector to be filled in.
   */
  static void fillInContrastValues(affymetrix_fusion_io::FusionCELData *cel, 
                                   std::vector<ProbeListPacked> &chrXProbeSets, 
                                   std::vector<double> &contrastValues, 
                                   double k = 4);

  /**
   * Calculate the contrast value for a single probeset
   * @param cel - Cel File to get data from.
   * @param ps - pointer to ProbeSet
   * @param contrastValues - Contrast values vector to be filled in.
   * @return contrast value
   */
  static double CalculateEMContrast(affymetrix_fusion_io::FusionCELData* cel, 
                                    const ProbeSet* ps, const double k);

  /**
   * @brief EM-based gender calling algorithm, call gender on a single
   * sample. 
   * @param contrast - vector of allele contrasts from SNPs on X chromosome
   * @param em_thresh - threshold on EM's confidence of labelling
   * @param em_cutoff - cutoff on total number of SNPs passed em's threshold, below cutoff return unknown
   * @param gender_cutoff - cutoff on percentage of heterozygotes of SNPs which passed the threshold
   *
   * @param gender - gender var to be filled in
   * @param hcr - estimated chrX het rate to be filled in
   */
   static void CallGender(const std::vector<double>& contrast,
                          const double em_thresh,
                          const double em_cutoff,
                          const double gender_cutoff,
                          affx::Gender &gender,
                          double &hcr);
    
   static void CallGender(const std::vector<double>& contrast,
                          affx::Gender &gender,
                          double &hcr);

 protected:

  /// Vector of chrX probesets that are not in the pseudo autosomal
  /// region and should be haploid in males and diploid in females.
  /// Pointer memory is owned elsewhere...
  std::vector<ProbeListPacked> m_ChrXProbeSets;
  /// Name of the cel files that have been called.
  std::vector<std::string> m_CelNames;
  /// Gender calls for chips seen.
  std::vector<affx::Gender> m_Genders;
  /// Our K parameter for the contrast centers transformation.
  double m_K; 
  /// EM Threshold to use
  double m_EmThresh;
  double m_EmCutOff;
  double m_GenderCutOff;
};

#endif /* EMGENDERCELLISTENER_H */
