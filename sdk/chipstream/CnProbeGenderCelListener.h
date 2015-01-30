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

#ifndef CNPROBEGENDERCELLISTENER_H
#define CNPROBEGENDERCELLISTENER_H

#include "chipstream/AptTypes.h"
#include "chipstream/BioTypes.h"
#include "chipstream/GenderCelListener.h"
//
#include <cstring>
#include <string>
#include <vector>
//

/**
 * @brief Class for determining gender from a cel set of cel files.
 */
class CnProbeGenderCelListener : public GenderCelListener {

public:

  /** 
   * Constructor.
   */
  CnProbeGenderCelListener(const std::string &chrXProbeFile, const std::string &chrYProbeFile,
                           double femaleThresh = 0.48, double maleThresh = 0.71, bool chrWZ = false);

  CnProbeGenderCelListener(const std::vector< std::vector<probeid_t> >& chrXProbes,
                           const std::vector< std::vector<probeid_t> >& chrYProbes,
                           double femaleThresh = 0.48, double maleThresh = 0.71, bool chrZW = false);
  /** 
   * Virtual destructor
   */
  virtual ~CnProbeGenderCelListener() {}

  /** 
   * Process another cel files worth of data.
   * @param cel - Filehandle to open cel file.
   */
  void newChip(affymetrix_fusion_io::FusionCELData *cel);

  /** 
   * Get the genders for the cel files that have been seen.
   * @return - vector of gender calls, one for each cel file in order seen.
   */
  std::vector<affx::Gender> getGenders() {
    return m_Genders;
  }

  /** 
   * Get the last gender called.
   * @return - gender call.
   */
  affx::Gender getLastGender() {return m_Genders[m_Genders.size() - 1];}
  float getLastRatio() {return m_vRatios[m_vRatios.size() - 1];}
  
  /** 
   * Get the names for the cel files that have been seen.
   * @return - vector of cel file names, one for each cel file in order seen.
   */
  std::vector<std::string> getCelNames() {
    return m_CelNames;
  }

  void initialize(double femaleThresh, double maleThresh);

protected:
  affx::Gender callGender(affymetrix_fusion_io::FusionCELData *cel, 
                          const std::vector< std::vector<probeid_t> > &chrxProbes, 
                          const std::vector< std::vector<probeid_t> > &chryProbes,
                          float& fRatio);

  /// Name of the cel files that have been called.
  std::vector<std::string> m_CelNames;
  /// Gender calls for chips seen.
  std::vector<affx::Gender> m_Genders;
  std::vector<float> m_vRatios;

  /// Probe ID vectors for ChrX and ChrY
  std::vector< std::vector<probeid_t> > m_ChrXProbes;
  std::vector< std::vector<probeid_t> > m_ChrYProbes;

  // strings for printing chromosome letters
  // For ZW gender calling, chrZ and W are used in place of X and Y
  // respectively.  Using this string enables the messaging to be appropriate
  // to the correct chromosome.
  std::string m_DipSexChrLetter;
  std::string m_HapSexChrLetter;
  bool m_ZWGenderCalling;
private:
    double m_FemaleThresh;
    double m_MaleThresh;
};

#endif /* CNPPROBEGENDERCELLISTENER_H */
