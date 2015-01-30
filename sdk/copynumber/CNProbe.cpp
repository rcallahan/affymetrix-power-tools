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


#include "copynumber/CNProbe.h"

/**
 * Constructor.
*/
CNProbe::CNProbe()
{
  m_iProbeSetIndex = -1;
  m_cAllele = 0;
  m_uiProbeID = 0;
  m_fIntensity = 0;
  m_fResidual = 0;
  m_dProbeEffect = 0;
  m_bUseForSketch = false;
  m_fMedianIntensity = 0;
  m_fPredictedIntensity = 0;
  m_nPDNNBin = 0;
}

/**
 * Compare function.
 * @param that - A reference to an instance of this class.
 * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
 * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
*/
int CNProbe::compareTo(CNProbe& that, int iCompareCode)
{
  int iCompareResult = 0;
  switch (iCompareCode) {
  case 0:
    iCompareResult = AffxArray<int>::compare(m_iProbeSetIndex, that.m_iProbeSetIndex);
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<char>::compare(m_cAllele, that.m_cAllele);
    }
    break;
  case 1:
    iCompareResult = AffxArray<unsigned int>::compare(m_uiProbeID, that.m_uiProbeID);
    break;
  case 2:
    iCompareResult = AffxArray<int>::compare(m_iProbeSetIndex, that.m_iProbeSetIndex);
    break;
  case 3:
    iCompareResult = AffxArray<int>::compare(m_iProbeSetIndex, that.m_iProbeSetIndex);
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<char>::compare(m_cAllele, that.m_cAllele);
    }
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<unsigned int>::compare(m_uiProbeID, that.m_uiProbeID);
    }
    break;
  }
  return iCompareResult;
}

CNProbeArray::CNProbeArray() {}
CNProbeArray::~CNProbeArray()
{
  deleteAll();
}

unsigned int CNProbeArray::getMaximumNumberProbesPerProbeSet()
{
  // Determine the maximum number of probes for a probe set
  quickSort(0); // by ProbeSet index, by allele
  unsigned int uiMaxProbeCount = 0;
  unsigned int uiProbeCount = 0;
  CNProbe* pPrev = NULL;
  for (int iIndex = 0; (iIndex < getCount()); iIndex++) {
    CNProbe* p = getAt(iIndex);
    if (p->getProbeSetIndex() == -1) {
      continue;
    }
    if ((pPrev == NULL) || (pPrev->compareTo(*p, 0) == 0)) {
      if (pPrev == NULL) {
        pPrev = p;
      }
      uiProbeCount++;
    } else {
      if (uiProbeCount > uiMaxProbeCount) {
        uiMaxProbeCount = uiProbeCount;
      }
      pPrev = p;
      uiProbeCount = 1;
    }
  }
  if (uiProbeCount > uiMaxProbeCount) {
    uiMaxProbeCount = uiProbeCount;
  }
  return uiMaxProbeCount;
}


