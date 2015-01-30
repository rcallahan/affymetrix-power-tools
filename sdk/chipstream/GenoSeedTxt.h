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
 * @file   GenoSeedTxt.h
 * @author Chuck Sugnet
 * @date   Fri Jun 30 15:13:27 PDT 2006
 * 
 * @brief Implements interface for a genotyping oracle used for
 * genotyping seeds in brlmm. In this case the genotype seed calls are
 * simply read from a text file.
 * 
 */

#ifndef _GENOSEEDTXT_H_
#define _GENOSEEDTXT_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/GenoSeed.h"
//
#include "util/Err.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

/** 
 * @brief Implements interface for a genotyping oracle used for
 * genotyping seeds in brlmm. In this case the genotype seed calls are
 * simply read from a text file.
 */
class GenoSeedTxt : public GenoSeed {

public:

  /** 
   * Open up a text file and read the genotypes from it. The chrX
   * snps will be used to determine the het call rate on chrX. If
   * supplied.
   * 
   * @param fileName - Name of text file to read genotypes from.
   * @param layout - Physical description of chip and probesets.
   * @param chrXSnps - Which snps are on chrX.
   */
  GenoSeedTxt(const std::string& fileName, ChipLayout &layout, 
              std::vector<std::string> &colName,
              std::map<std::string,bool> &chrXSnps);

  /** How many chips does this object know about. */
  int getChipCount() { return m_ChipCount; }

  /** 
   * Get the call rates for snps for each chip.
   * @param callRates - Vector to be filled in with call rates for chrx snps,
   * one for each chip.
   */
  void getCallRates(std::vector<float> &callRates) {
    callRates = m_CallRates;
  }

  /** 
   * Get the call rates for chrX snps for each chip.
   * @param callRates - Vector to be filled in with call rates for chrx snps,
   * one for each chip.
   */
  void getHetChrXRates(std::vector<float> &hetRates) {
    hetRates = m_HetChrXCallRates;
  }

  virtual float getHetChrXRate(int chip) {
      return m_HetChrXCallRates[chip];
  }

      
  /** 
   * Get the genotype calls for a particular probeset.
   * @param name - name of probeset to get genotype calls for.
   * @return - A vector of the calls from dm algorithm.
   */
  std::vector<affx::GType> getGenoCalls(const std::string &name) {
    std::map<std::string, std::vector<affx::GType> >::iterator i = m_KnownCalls.find(name);
    if(i == m_KnownCalls.end()) {
      Err::errAbort("GenoSeedTxt::getGenoCalls() - Don't recognize probeset '" + name + "'");
    }
    return i->second;
  }

  /**
   * @brief Check if a call exists for a given probeset
   * @param name - name of probeset to get genotype calls for.
   * @return boolean result
   */
  inline bool checkGenoCallsName(const std::string &name) {
    return m_KnownCalls.end() != m_KnownCalls.find(name);
  }


  void fillProbesets(std::vector<std::string> &probesets) {
      std::map<std::string, std::vector<affx::GType> >::iterator iter;
      for(iter = m_KnownCalls.begin(); iter != m_KnownCalls.end(); iter++) {
        probesets.push_back(iter->first);
      }
  }

protected:

  /** 
   * Read in the precomputed genotypes from a tab delimited file. Must
   * be a column for each Celfile with same header name as celfile and an row
   * for every snp to be analyzed. There can be extra columns and order of columns
   * is not important as the resulting order will be the same as the colNames
   * vector.
   * 
   * @param fileName - path to text file to read genotypes from.
   * @param colNames - Name of the columns we're expecting to find in this file.
   * @param gtypes - Map to store the genotypes in.
   */
  void readInPrecompGenoTypes(const std::string& fileName, 
                                         std::vector<std::string> &colNames,
                                         std::map<std::string,bool> &chrXSnps);

  std::vector<float> m_CallRates; ///< Call rates for each chip.
  std::vector<float> m_HetChrXCallRates; ///< Het call rates for chrX for each chip.
  std::map<std::string, std::vector<affx::GType> > m_KnownCalls; ///< Calls for each SNP.
  int m_ChipCount;  ///< How many chips do we have calculations for.
};

#endif /* _GENOSEEDTXT_H_ */
