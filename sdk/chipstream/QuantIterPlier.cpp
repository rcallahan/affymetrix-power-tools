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
 * @file   QuantIterPlier.cpp
 * @author Chuck Sugnet
 * @date   Mon Oct 24 11:48:48 2005
 * 
 * @brief Class for doing probe set quantification estimate by iteratively
 * calling PLIER with the probes that best correlate with signal estimate.
 */

//
#include "chipstream/QuantIterPlier.h"
//
#include "util/Verbose.h"

/** 
 * Constructor.  
 */
QuantIterPlier::QuantIterPlier(QuantPlierParams &plierParams)
  : QuantPlierBase(plierParams) {
  initMemory();
  setupSelfDoc(*this);
  m_Type = getDocName();
  // this is bogus -- we dont know how much we need!
  // Initialize memory buffers.
  // allocMemory(100, 50);
  setUsePrecompFeatureEffects(false);
}

/**
 * Destructor.
 */
QuantIterPlier::~QuantIterPlier() {
  freeMemory();
}

void QuantIterPlier::initMemory() {
  m_MaxChips=0;
  m_MaxProbes=0;

  m_MM=NULL;
  m_PM=NULL;
  m_Residuals=NULL;

  m_ChipEffects=NULL;
  m_ProbeEffects=NULL;
  m_ProbesUsed=NULL;
}

/** 
 * @brief Free up all the memory that has bee allocated.
 */
void QuantIterPlier::freeMemory() {
  for(unsigned int i = 0; i < m_MaxChips; i++) {
    delete [] m_PM[i];
    delete [] m_MM[i];
    delete [] m_Residuals[i];
  }
  delete [] m_MM;
  delete [] m_PM;
  delete [] m_Residuals;
  m_MaxChips=0;
  m_MaxProbes=0;
  //
  delete [] m_ChipEffects;
  delete [] m_ProbeEffects;
  delete [] m_ProbesUsed;
}

/** 
 * @brief Allocate enough memory for at least maxChips and
 * maxProbes.
 * 
 * @param maxChips - Number of chips that can be computed.
 * @param maxProbes - Number of probes or features that can be computed.
 */
void QuantIterPlier::allocMemory(unsigned int maxChips, unsigned int maxProbes) {
  // fre what we have.
  freeMemory();

  m_MaxChips = maxChips;
  m_MaxProbes = maxProbes;

  m_ProbeEffects = new double[m_MaxProbes];
  m_ProbesUsed = new  int[m_MaxProbes];
  m_ChipEffects = new double[m_MaxChips];

  m_PM = new double *[m_MaxChips];
  m_MM = new double *[m_MaxChips];
  m_Residuals = new double *[m_MaxChips];
  for(unsigned int i = 0; i < m_MaxChips; i++) {
    m_PM[i] = new double[m_MaxProbes];
    m_MM[i] = new double[m_MaxProbes];
    m_Residuals[i] = new double[m_MaxProbes];
  }
}

/** 
 * @brief Do the heavy lifting of plier.
 */
void QuantIterPlier::computeEstimate() {
  long errorCode = 0;
  unsigned int i = 0;
  int *moreIterations = (int *)malloc(m_Iterations.size()*sizeof(int));

  m_IterPlier.memEnsureSize(m_ChipCount,m_ProbeCount);

  for (i = 0; i < m_Iterations.size(); i++) {
    moreIterations[i] = m_Iterations[i];
  }
  errorCode = m_IterPlier.runCorrIterPlier(m_Plier, m_ProbeCount, m_ChipCount,
                                           m_PM, m_MM, m_ChipEffects, m_ProbeEffects,
                                           m_Residuals, m_ProbesUsed, m_Iterations.size(), 
                                           moreIterations);
  free((void *)moreIterations);
  if(errorCode != 0) 
    Err::errAbort("Problem running plier. Error code: " + ToStr(errorCode));
}
