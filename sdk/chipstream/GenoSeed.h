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
 * @file   GenoSeed.h
 * @author Chuck Sugnet
 * @date   Fri Jun 30 14:44:22 2006
 * 
 * @brief Abstract interface for a genotyping oracle used for
 * genotyping seeds in brlmm. Originally seed calls came from DM, but
 * now they are a subject of much research.
 * 
 */

#ifndef _GENOSEED_H_
#define _GENOSEED_H_

#include "chipstream/BioTypes.h"
#include "chipstream/GenderCalls.h"

/** Abstract interface for a genotyping oracle used for genotyping
    seeds in brlmm. */
class GenoSeed : public GenderCalls {

public:

  /** How many chips does this object know about. */
  virtual int getChipCount() = 0;

  /** 
   * Get the call rates for snps for each chip.
   * @param callRates - Vector to be filled in with call rates for chrx snps,
   * one for each chip.
   */
  virtual void getCallRates(std::vector<float> &callRates) = 0;

  /** 
   * Get the call rates for chrX snps for each chip.
   * @param callRates - Vector to be filled in with call rates for chrx snps,
   * one for each chip.
   */
  virtual void getHetChrXRates(std::vector<float> &hetRates) = 0;

  virtual float getHetChrXRate(int chip) = 0;

  /** 
   * Get the genotype calls for a particular probeset.
   * @param name - name of probeset to get genotype calls for.
   * @return - A vector of the calls from dm algorithm.
   */
  virtual std::vector<affx::GType> getGenoCalls(const std::string &name) = 0;

  /**
   * @brief Check if a call exists for a given probeset
   * @param name - name of probeset to get genotype calls for.
   * @return boolean result
   */
  virtual bool checkGenoCallsName(const std::string &name) = 0;

  /** Virtual destructor for virtual class. */
  virtual ~GenoSeed();

  /**
   * Figure out the gender of samples based on the number of het calls on
   * chromosome X. This is all very human centric...
   */
  virtual std::vector<affx::Gender> getGenders();

  virtual affx::Gender getGender(int chip);

};

#endif /* _GENOSEED_H_ */
