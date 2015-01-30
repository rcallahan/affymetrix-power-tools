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
 * @file   DM.h
 * @author Xiaojun Di
 * @date   March 8 2006
 * 
 * DRAFT: wraps the DM algorithmics for genotyping under SDK
 * 
 */

#ifndef __DM_H__
#define __DM_H__

//
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

/// infinity of floats
#define FLOATINFINITY	numeric_limits<float>::max()
/// number of cells per quartets
#define NUMCELLS		4
/// log of 2 pi
#define LOG2PI			log(6.283185)	

/**
* An object to probe quartet, each quartet contains 4 cells
* each cell contains intensity, variance and number of pixels, store variance instead of stdev
*/
class CQuartet
{
public:
  
  /// alleles
  typedef enum
    {
      A_Allele,
      B_Allele
    } Allele;
  
  /// match types
  typedef enum
    {
      PM,
      MM
    } MatchType;
  
  /// four cells in the order of PM_A, MM_A, PM_B and MM_B
  /// number of pixels per cell
  int pixels[NUMCELLS];
  /// intensities of each cell of the quartet
  float intensity[NUMCELLS];
  /// variances of each cell of the quartet
  float variance[NUMCELLS];
  
  /**
  * @brief get the probe quartet given the allele type and match type
  * @param allele - allele type: A or B
  * @param matchtype - the match type: perfect or mis
  * @return the probe index
  */
  int GetIndex(const int& allele, const int& matchtype);
  
  /**
  * @brief default constructor
  */
  CQuartet();
  
  /**
  * @brief set intensity for allele type=allele and match type = matchtype
  * @param allele - allele type: A or B
  * @param matchtype - the match type: perfect or mis
  * @param val - the value to be set
  * @return void
  */
  inline void SetIntensity(const int& allele, const int& matchtype, const float& val)
  {
    intensity[GetIndex(allele,matchtype)] = val;
  }
  
  /**
  * @brief set variance for allele type=allele and match type = matchtype
  * @param allele - allele type: A or B
  * @param matchtype - the match type: perfect or mis
  * @param val - the value to be set
  * @return void
  */  
  inline void SetVariance(const int& allele, const int& matchtype, const float& val)
  {
    variance[GetIndex(allele,matchtype)] = val;
  }
  
  /**
  * @brief set number of pixels for allele type=allele and match type = matchtype
  * @param allele - allele type: A or B
  * @param matchtype - the match type: perfect or mis
  * @param val - the value to be set
  * @return void
  */  
  inline void SetPixels(const int& allele, const int& matchtype, const int& val)
  {
    pixels[GetIndex(allele,matchtype)] = val;
  }
};

/**
* the global function to invoke DM algorithm
* @param vquartets - a vector of quartets containing not just the intensities but also the variances and pixels
* @param hetMult - het multiplier, hetMult is used to adjust the likelihood of AB by multiplying a constant factor, 
* which is equivalent to adding a constant (hetMult) to its loglikelihood, default to 0 or no adjustment is requested
*
* @param results -- pair of p-value and associated genotype call to be returned
* @return true if no error, false otherwise
*/
bool CallDM(const std::vector<CQuartet> &vquartets, std::pair<float,int>& results, const float& hetMult = 0);

#endif // __DM_H__
