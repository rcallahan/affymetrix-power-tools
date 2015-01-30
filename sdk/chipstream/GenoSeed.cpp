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
 * @file   GenoSeed.cpp
 * @author Chuck Sugnet
 * @date   Fri Jun 30 14:44:22 2006
 * 
 * @brief Abstract interface for a genotyping oracle used for
 * genotyping seeds in brlmm. Originally seed calls came from DM, but
 * now they are a subject of much research.
 * 
 */

#include "chipstream/GenoSeed.h"


/** Virtual destructor for virtual class. */
GenoSeed::~GenoSeed()
{
}

/**
 * Figure out the gender of samples based on the number of het calls on
 * chromosome X. This is all very human centric...
 */
std::vector<affx::Gender> GenoSeed::getGenders() {

  std::vector<affx::Gender> genders;

  /* Loop through and decide if we are male or female. */
  for(unsigned int i = 0; i < getChipCount(); i++) {
    genders.push_back(getGender(i));
  }
  return genders;
}

affx::Gender GenoSeed::getGender(int chip) {
  double hetThreshold = .075; // This comes from gtype thresholds.
  float hetRate = getHetChrXRate(chip);

  if(hetRate == -1) {
    return affx::UnknownGender;
  }
  else {
    if(hetRate < hetThreshold) {
      return affx::Male;
    }
    else {
      return affx::Female;
    }
  }
}

