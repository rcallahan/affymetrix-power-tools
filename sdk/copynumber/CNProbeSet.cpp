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
 * @file CNProbeSet.cpp
 *
 * @brief This file contains the ProbeSet class definitions.
 */

#include "copynumber/CNProbeSet.h"
//

CNProbeSet::CNProbeSet()
{
  m_pSnpDistribution = NULL;
  clear();
}

CNProbeSet::~CNProbeSet()
{
  if (m_pSnpDistribution != NULL) {
    delete m_pSnpDistribution;
  }
}

void CNProbeSet::clear()
{
  m_bProcess = true;
  m_vWaves.clear();
  m_vAllCovariates.clear();
  m_strProbeSetName.clear();
  m_cChromosome = 0;
  m_iPosition = 0;
  m_fMedianSignal = 0;
  m_fXXMedianSignal = 0;
  m_fYMedianSignal = 0;
  m_fAAMedianSignal = 0;
  m_fABMedianSignal = 0;
  m_fBBMedianSignal = 0;


  m_cProcessFlag = 0;
  m_cStyAdapterCode = -1;
  m_cNspAdapterCode = -1;
  m_fGCContent = 0;
  m_fAMedianIntensity = 0;
  m_fBMedianIntensity = 0;
  m_fAAlleleSignal = 0;
  m_fBAlleleSignal = 0;
  m_fLog2Ratio = 0;
  m_fLog2RatioMedianSmooth = 0;
  m_fAllelicDifference = 0;
  m_fGcAdjustment = 0;
  m_bPseudoAutosomalRegion = false;
  m_iGCBinIndex = -1;
  m_cGenotypeCall = -1;
  m_fGenotypeConfidence = 0;
  m_iCNState = -1;
  m_iImputedCNState = -1;
  m_fCNConfidence = 0;
  m_fSmoothedLog2Ratio = 0;
  m_fLoh = 0;
  m_fSCAR = 0;
  m_fFLD = 0.0;
  m_iMaxPeaks = 0;
  m_bUseForSketch = false;
  if (m_pSnpDistribution != NULL) {
    delete m_pSnpDistribution;
    m_pSnpDistribution = NULL;
  }
  m_fMuAA = 0;
  m_fMuAB = 0;
  m_fMuBB = 0;
  m_fInformation = 0;
  m_iUseInEMAlgorithm = 0;
  m_fMosaicismMixture = numeric_limits<float>::quiet_NaN();
  m_iReplicateCount = 0;
  m_uiAllelePeaks1 = 0;
  m_uiAllelePeaks2 = 0;
  m_bValidSCARExists=false;
  m_bValidFLDExists = false;
  m_bValidHomHetExists=true;
  m_bTrulySNP=false;
  m_bTrulyCN=false;
}

float CNProbeSet::getSignalContrast(double dK)
{
  if ((m_fAAlleleSignal + m_fBAlleleSignal) == 0) {
    throw(Except("Zero median signals found. CEL file may be corrupted."));
  }
  float fContrast = (float)(sinh(dK * (m_fAAlleleSignal - m_fBAlleleSignal) / (m_fAAlleleSignal + m_fBAlleleSignal)) / sinh(dK));
  if (fabs(fContrast) > 1) {
    throw(Except("Can't have abs(contrast) > 1."));
  }
  return fContrast;
}

float CNProbeSet::getSignalContrastMvA()
{
  const double eps = 0.001;
  if (m_fAAlleleSignal < eps)
      m_fAAlleleSignal = eps;
  if (m_fBAlleleSignal < eps)
      m_fBAlleleSignal = eps;

  return log2(m_fBAlleleSignal) - log2(m_fAAlleleSignal);
}

float CNProbeSet::getSignalStrength()
{
  if ((m_fAAlleleSignal + m_fBAlleleSignal) == 0) {
    throw(Except("Zero median signals found. CEL file may be corrupted."));
  }
  return log2(m_fAAlleleSignal + m_fBAlleleSignal);
}

float CNProbeSet::getSignalStrengthMvA()
{
  const double eps = 0.001;
  if (m_fAAlleleSignal < eps)
      m_fAAlleleSignal = eps;
  if (m_fBAlleleSignal < eps)
      m_fBAlleleSignal = eps;

  return (log2(m_fAAlleleSignal) + log2(m_fBAlleleSignal))/2.0;
}

float CNProbeSet::getIntensityContrast(double dK)
{
  if ((m_fAMedianIntensity + m_fBMedianIntensity) == 0) {
    throw(Except("Zero median intensities found. CEL file may be corrupted."));
  }
  float fContrast = (float)(sinh(dK * (m_fAMedianIntensity - m_fBMedianIntensity) / (m_fAMedianIntensity + m_fBMedianIntensity)) / sinh(dK));
  if (fabs(fContrast) > 1) {
    throw(Except("Can't have abs(contrast) > 1."));
  }
  return fContrast;
}

float CNProbeSet::getIntensityStrength()
{
  if ((m_fAMedianIntensity + m_fBMedianIntensity) == 0) {
    throw(Except("Zero median intensities found. CEL file may be corrupted."));
  }
  return log2(m_fAMedianIntensity + m_fBMedianIntensity);
}

float CNProbeSet::getCalibratedLog2Ratio(
  double dAlphaCNCalibrate,
  double dBetaCNCalibrate)
//  Sacrifice a bit of compiler optimization for readability.  Commented out
//  is the original code for reference.
//  {
//  return (float)exp(((getLog2Ratio() / dAlphaCNCalibrate) + dBetaCNCalibrate) * log(2.0));}
{
  double log2ratio = getLog2Ratio();
  log2ratio /= dAlphaCNCalibrate;
  log2ratio += dBetaCNCalibrate;
  return (float)exp(log2ratio*log(2.0));
}

float CNProbeSet::getCalibratedSmoothedLog2Ratio(
  double dAlphaCNCalibrate,
  double dBetaCNCalibrate)
//  Sacrifice a bit of compiler optimization for readability.  Commented out
//  is the original code for reference.
//  {
//  return (float)exp(((getLog2Ratio() / dAlphaCNCalibrate) + dBetaCNCalibrate) * log(2.0));}
{
  double log2ratio = getSmoothedLog2Ratio();
  log2ratio /= dAlphaCNCalibrate;
  log2ratio += dBetaCNCalibrate;
  return (float)exp(log2ratio*log(2.0));
}

/**
 * Compare function.
 * @param that - A reference to an instance of this class.
 * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
 * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
*/
int CNProbeSet::compareTo(CNProbeSet& that, int iCompareCode)
{
  int iCompareResult = 0;
  switch (iCompareCode) {
  case 0:
    //iCompareResult = m_strProbeSetName.compareTo(that.m_strProbeSetName, 0);
    iCompareResult = strcmp(m_strProbeSetName.c_str(), that.m_strProbeSetName.c_str());
    break;
  case 1:
    iCompareResult = AffxArray<char>::compare(m_cChromosome, that.m_cChromosome);
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<int>::compare(m_iPosition, that.m_iPosition);
    }
    //if (iCompareResult == 0) {iCompareResult = m_strProbeSetName.compareTo(that.m_strProbeSetName, 0);}
    if (iCompareResult == 0) {
      iCompareResult = strcmp(m_strProbeSetName.c_str(), that.m_strProbeSetName.c_str());
    }
    break;
  }
  return iCompareResult;
}

/**
 * @brief Return the CNNeutral call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNProbeSet::getCNNeutral(int iSexCode, int iXChromosome, int iYChromosome)
{
//  return ((m_iCNState == 2) ? 1 : 0);
  return (((getCNGain(iSexCode, iXChromosome, iYChromosome) == 0) && (getCNLoss(iSexCode, iXChromosome, iYChromosome) == 0)) ? (char)1 : (char)0);
}

/**
 * @brief Return the CNGain call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNProbeSet::getCNGain(int iSexCode, int iXChromosome, int iYChromosome)
{
//  return ((m_iCNState > 2) ? 1 : 0);
  char cCNGain = 0;
  if ((m_cChromosome == iXChromosome) && (!m_bPseudoAutosomalRegion)) {
    if (iSexCode == affx::Male) {
      if (m_iCNState > 1) {
        cCNGain = 1;
      }
    } else {
      if (m_iCNState > 2) {
        cCNGain = 1;
      }
    }
  } else if ((m_cChromosome == iYChromosome) && (!m_bPseudoAutosomalRegion)) {
    if (iSexCode == affx::Male) {
      if (m_iCNState > 1) {
        cCNGain = 1;
      }
    } else {
      if (m_iCNState > 0) {
        cCNGain = 1;
      }
    }
  } else {
    if (m_iCNState > 2) {
      cCNGain = 1;
    }
  }
  return cCNGain;
}

/**
 * @brief Return the CNLoss call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNProbeSet::getCNLoss(int iSexCode, int iXChromosome, int iYChromosome)
{
//  return ((m_iCNState < 2) ? 1 : 0);
  char cCNLoss = 0;
  if ((m_cChromosome == iXChromosome) && (!m_bPseudoAutosomalRegion)) {
    if (iSexCode == affx::Male) {
      if (m_iCNState < 1) {
        cCNLoss = 1;
      }
    } else {
      if (m_iCNState < 2) {
        cCNLoss = 1;
      }
    }
  } else if ((m_cChromosome == iYChromosome) && (!m_bPseudoAutosomalRegion)) {
    if (iSexCode == affx::Male) {
      if (m_iCNState < 1) {
        cCNLoss = 1;
      }
    } else {
      cCNLoss = 0;
//    if (m_iCNState < 0) {cCNLoss = 1;}
    }
  } else {
    if (m_iCNState < 2) {
      cCNLoss = 1;
    }
  }
  return cCNLoss;
}

/**
 * @brief Return the CNNeutralLOH call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNProbeSet::getCNNeutralLOH(int iSexCode, int iXChromosome, int iYChromosome)
{
  //char cCNNeutral = getCNNeutral(iSexCode, iXChromosome, iYChromosome);
  // Always assume CNNeutral is CN == 2, per Jim and Carl 02/18/2009
  return (((m_iCNState == 2) && (m_fLoh == 1)) ? (char)1 : (char)0);
}

/**
 * @brief Return the NormalDipload call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNProbeSet::getNormalDiploid(int iSexCode, int iXChromosome, int iYChromosome)
{
  //char cCNNeutral = getCNNeutral(iSexCode, iXChromosome, iYChromosome);
  // Always assume CNNeutral is CN == 2, per Jim and Carl 02/18/2009
  return (((m_iCNState == 2) && (m_fLoh == 0)) ? (char)1 : (char)0);
}

/**
 * @brief Return the specified segment type call
 * @param int - The segment type code
 * @param int - The sex code for this sample
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNProbeSet::getSegmentTypeCall(int iSegmentType, int iSexCode, int iXChromosome, int iYChromosome)
{
  char cCall = -1;
  switch (iSegmentType) {
  case 1:
    cCall = getCNState();
    break;
  case 2:
    cCall = ((getLoh() == 1) ? (char)1 : (char)0);
    break;
  case 3:
    cCall = ((getLoh() == 1) ? (char)1 : (char)0);
    break;
  case 4:
    cCall = getCNNeutralLOH(iSexCode, iXChromosome, iYChromosome);
    break;
  case 5:
    cCall = getNormalDiploid(iSexCode, iXChromosome, iYChromosome);
    break;
  }
  return cCall;
}

/**
 * @brief Return the specified segment type confidence
 * @param int - The segment type code
 * @return float - The comnfidence
 */
float CNProbeSet::getSegmentTypeConfidence(int iSegmentType)
{
  float fConfidence = 0;
  switch (iSegmentType) {
  case 1:
    fConfidence = getCNConfidence();
    break;
  }
  return fConfidence;
}

//bool CNProbeSet::lohCalled()
//{
//  return (m_fLoh==0.0 || m_fLoh==1.0) ;
//}

/**
 * @brief Return the affx genotype call code
 * @return char - The call
 */
char CNProbeSet::getGenotypeCallCode()
{
  switch (m_cGenotypeCall) {
  case 0: return ALLELE_A_CALL; break;
  case 1: return ALLELE_AB_CALL; break;
  case 2: return ALLELE_B_CALL; break;
  default: return ALLELE_NO_CALL;
  }
  return ALLELE_NO_CALL;
}

void CNProbeSet::setCovariateValue(int index, float value)
{
    if (m_vAllCovariates.size() <= index)
    {
        int oldSize = m_vAllCovariates.size();

        // swap trick to force the capacity equal EXACTLY to index+1 (which resize() does not do!)
        m_vAllCovariates.resize(index+1);
        std::vector<float>(m_vAllCovariates).swap(m_vAllCovariates);

        // pad the gap (if any) with NaNs
        for (int i = oldSize; i < index; i++)
        {
            m_vAllCovariates[i] = numeric_limits<float>::quiet_NaN();
        }
    }
    m_vAllCovariates[index] = value;
}

CNProbeSetArray::CNProbeSetArray() {}
CNProbeSetArray::~CNProbeSetArray()
{
  nullAll();
}

bool CNProbeSetArray::isSketchSubset()
{
  unsigned int uiCount = 0;
  for (int iIndex = 0; (iIndex < getCount()); iIndex++) {
    CNProbeSet* p = getAt(iIndex);
    if (p->isUseForSketch()) {
      uiCount++;
    }
  }
  return (uiCount == getCount());
}
int CNProbeSetArray::getCNProcessCount()
{
  unsigned int uiCount = 0;
  for (int iIndex = 0; (iIndex < getCount()); iIndex++) {
    CNProbeSet* p = getAt(iIndex);
    if (p->isProcess() && p->processAsCN()) {
      uiCount++;
    }
  }
  return uiCount;
}


int CNProbeSetArray::getProcessCount()
{
  unsigned int uiCount = 0;
  for (int iIndex = 0; (iIndex < getCount()); iIndex++) {
    CNProbeSet* p = getAt(iIndex);
    if (p->isProcess() ) {
      uiCount++;
    }
  }
  return uiCount;
}

void CNProbeSetArray::calculateLog2RatioMedianSmooth(int iWindowSize)
{
  if ((iWindowSize % 2) == 0) {
    iWindowSize++;
  } // Window size should be odd.
  if (getCount() < iWindowSize) {
    return;
  }
  int iHalfWindow = (int)(iWindowSize / 2);
  AffxMultiDimensionalArray<float> v(iWindowSize);
  for (int iIndex = iHalfWindow; (iIndex < (getCount() - iHalfWindow)); iIndex++) {
    int j = 0;
    for (int i = (iIndex - iHalfWindow); (i <= (iIndex + iHalfWindow)); i++) {
      v.set(j, getAt(i)->getLog2Ratio()); j++;
    }
    getAt(iIndex)->setLog2RatioMedianSmooth(v.median());
  }
  for (int iIndex = (iHalfWindow - 1); (iIndex >= 0); iIndex--) {
    getAt(iIndex)->setLog2RatioMedianSmooth(getAt(iIndex + 1)->getLog2RatioMedianSmooth());
  }
  for (int iIndex = (getCount() - iHalfWindow); (iIndex < getCount()); iIndex++) {
    getAt(iIndex)->setLog2RatioMedianSmooth(getAt(iIndex - 1)->getLog2RatioMedianSmooth());
  }
}

void CNProbeSetArray::setupGCCorrectionBins(int iGCBinCount)
{
  AffxMultiDimensionalArray<float> vGcContent(getProcessCount());
  int iIndex = 0;
  for (int iRowIndex = 0; (iRowIndex < getCount()); iRowIndex++) {
    CNProbeSet* p = getAt(iRowIndex);
    if (!p->isProcess() || !( p->processAsCN() || p->processAsSNP() )) {
      continue;
    }
    float fGcContent = p->getGCContent();
    vGcContent.set(iIndex, fGcContent);
    iIndex++;
  }
  AffxMultiDimensionalArray<float> v(iGCBinCount);
  AffxMultiDimensionalArray<int> vCounts(iGCBinCount);
  float fPercentile = 0;
  for (int iBinIndex = 0; (iBinIndex < iGCBinCount); iBinIndex++) {
    fPercentile = vGcContent.percentile((1.0 / (double)iGCBinCount) * (iBinIndex + 1));
    v.set(iBinIndex, fPercentile);
  }
  v.set((iGCBinCount - 1), 1);
  for (int iRowIndex = 0; (iRowIndex < getCount()); iRowIndex++) {
    CNProbeSet* pobjProbeSet = getAt(iRowIndex);
    if (!pobjProbeSet->isProcess() || !( pobjProbeSet->processAsCN() || pobjProbeSet->processAsSNP() )) {
      continue;
    }
    float fGcContent = pobjProbeSet->getGCContent();
    for (int iBinIndex = 0; (iBinIndex < iGCBinCount); iBinIndex++) {
      if (fGcContent <= v.get(iBinIndex)) {
        if (iBinIndex == 0) {
          pobjProbeSet->setGCBinIndex(iBinIndex);
          vCounts.increment(iBinIndex);
          break;
        } else if (fGcContent > v.get(iBinIndex - 1)) {
          pobjProbeSet->setGCBinIndex(iBinIndex);
          vCounts.increment(iBinIndex);
          break;
        }
      }
    }
    if (pobjProbeSet->getGCBinIndex() == -1) {
      Verbose::out(1, pobjProbeSet->getProbeSetName() + " has no gc-bin assignment. " + ::getDouble(fGcContent, 6));
    }
  }
  /*
  Verbose::out(1, "*");
  AffxString str;
  for (int iBinIndex = 0; (iBinIndex < iGCBinCount); iBinIndex++)
  {
   if (iBinIndex == 0)
   {
    Verbose::out(1, str + "GC Correction Bin " + ::getInt(iBinIndex + 1) + "\t= 0       \t<= n <= " + ::getDouble(v.get(iBinIndex), 6) + "\tMarkerCount = " + ::getInt(vCounts.get(iBinIndex)));
   }
   else if (iBinIndex == (iGCBinCount - 1))
   {
    Verbose::out(1, str + "GC Correction Bin " + ::getInt(iBinIndex + 1) + "\t= " + ::getDouble(v.get(iBinIndex - 1), 6) + "\t<  n <= 1       \tMarkerCount = " + ::getInt(vCounts.get(iBinIndex)));
   }
   else
   {
    Verbose::out(1, str + "GC Correction Bin " + ::getInt(iBinIndex + 1) + "\t= " + ::getDouble(v.get(iBinIndex - 1), 6) + "\t<  n <= " + ::getDouble(v.get(iBinIndex), 6) + "\tMarkerCount = " + ::getInt(vCounts.get(iBinIndex)));
   }
  }
  Verbose::out(1, "*");
  */
}


