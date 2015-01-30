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

#ifndef _CNProbe_H_
#define _CNProbe_H_
/**
 * @file CNProbe.h
 *
 * @brief This header contains the CNProbe class definition.
 */

#include "util/AffxArray.h"

/**
 * @brief  A class for storing probe intensity data.
 *
 */
class CNProbe
{
private:
  int m_iProbeSetIndex;
  char m_cAllele;
  unsigned int m_uiProbeID;
  float m_fIntensity;
  float m_fResidual;
  double m_dProbeEffect;
  bool m_bUseForSketch;
  float m_fMedianIntensity;
  float m_fPredictedIntensity;
  short m_nPDNNBin;

public:
  /**
   * Constructor.
  */
  CNProbe();
  void setProbeSetIndex(int i) { m_iProbeSetIndex = i; }
  int getProbeSetIndex() { return m_iProbeSetIndex; }
  void setAllele(char c) { m_cAllele = c; }
  char getAllele() { return m_cAllele; }
  void setProbeID(unsigned int ui) { m_uiProbeID = ui; }
  unsigned int getProbeID() { return m_uiProbeID; }
  void setIntensity(float f) { m_fIntensity = f; }
  float getIntensity() { return m_fIntensity; }
  void setResidual(float f) { m_fResidual = f; }
  float getResidual() { return m_fResidual; }
  void setProbeEffect(double d) { m_dProbeEffect = d; }
  double getProbeEffect() { return m_dProbeEffect; }
  bool isUseForSketch() { return m_bUseForSketch; }
  void setUseForSketch(bool b) { m_bUseForSketch = b; }
  void setMedianIntensity(float f) { m_fMedianIntensity = f; }
  float getMedianIntensity() { return m_fMedianIntensity; }
  void setPredictedIntensity(float f) { m_fPredictedIntensity = f; }
  float getPredictedIntensity() { return m_fPredictedIntensity; }
  void setPDNNBin(short n) { m_nPDNNBin = n; }
  short getPDNNBin() { return m_nPDNNBin; }

  /**
   * Compare function.
   * @param that - A reference to an instance of this class.
   * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
   * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
  */
  int compareTo(CNProbe& that, int iCompareCode);

  template<int k> struct ComparePred {
      bool operator()(const CNProbe* lhs, const CNProbe* rhs) const {
          Err::errAbort("CNProbe: ComparePred instantiated with an invalid compare code = " + ToStr(k));
          return false;
      }
  };
};

template<> struct CNProbe::ComparePred<0> {
    bool operator()(const CNProbe* lhs, const CNProbe* rhs) const {
        int iCompareResult = 0;
        iCompareResult = AffxArray<int>::compare(lhs->m_iProbeSetIndex, rhs->m_iProbeSetIndex);
        if (iCompareResult == 0) {
            iCompareResult = AffxArray<char>::compare(lhs->m_cAllele, rhs->m_cAllele);
        }
        return iCompareResult < 0;
    }
};
template<> struct CNProbe::ComparePred<1> {
    bool operator()(const CNProbe* lhs, const CNProbe* rhs) const {
        return AffxArray<unsigned int>::compare(lhs->m_uiProbeID, rhs->m_uiProbeID) < 0;
    }
};
template<> struct CNProbe::ComparePred<2> {
    bool operator()(const CNProbe* lhs, const CNProbe* rhs) const {
        return AffxArray<int>::compare(lhs->m_iProbeSetIndex, rhs->m_iProbeSetIndex) < 0;
    }
};
template<> struct CNProbe::ComparePred<3> {
    bool operator()(const CNProbe* lhs, const CNProbe* rhs) const {
        int iCompareResult = 0;
        iCompareResult = AffxArray<int>::compare(lhs->m_iProbeSetIndex, rhs->m_iProbeSetIndex);
        if (iCompareResult == 0) {
            iCompareResult = AffxArray<char>::compare(lhs->m_cAllele, rhs->m_cAllele);
        }
        if (iCompareResult == 0) {
            iCompareResult = AffxArray<unsigned int>::compare(lhs->m_uiProbeID, rhs->m_uiProbeID);
        }
        return iCompareResult < 0;
    }
};

/**
 * @brief  A vector of CNProbe pointers
 */
class CNProbeArray : public AffxArray<CNProbe>
{
public:
  CNProbeArray();
  ~CNProbeArray();

  unsigned int getMaximumNumberProbesPerProbeSet();
};

#endif


