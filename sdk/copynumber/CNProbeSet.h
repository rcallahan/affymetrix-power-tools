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

#ifndef _CNProbeSet_H_
#define _CNProbeSet_H_
/**
 * @file CNProbeSet.h
 *
 * @brief This header contains the ProbeSet class definitions.
 */

#include "chipstream/BioTypes.h"
#include "file/CHPFileData.h"
//#include "label/snp.label.h"
#include "util/AffxArray.h"
#include "util/AffxMultiDimensionalArray.h"
//

class snp_distribution;

/**
 * @brief  A class to store probe set data.
 *
 */
class CNProbeSet
{
private:
  bool m_bProcess;
  std::vector<float> m_vWaves;
  std::vector<float> m_vAllCovariates;
  std::string m_strProbeSetName;
  unsigned char m_cChromosome;
  int m_iPosition;
  float m_fMedianSignal;
  float m_fXXMedianSignal;
  float m_fYMedianSignal;
  float m_fAAMedianSignal;
  float m_fABMedianSignal;
  float m_fBBMedianSignal;
  char m_cProcessFlag;
  char m_cStyAdapterCode;
  char m_cNspAdapterCode;
  float m_fGCContent;
  float m_fAMedianIntensity;
  float m_fBMedianIntensity;
  float m_fAAlleleSignal;
  float m_fBAlleleSignal;
  float m_fLog2Ratio;
  float m_fLog2RatioMedianSmooth;
  float m_fAllelicDifference;
  float m_fGcAdjustment;
  bool m_bPseudoAutosomalRegion;
  int m_iGCBinIndex;
  char m_cGenotypeCall;
  float m_fGenotypeConfidence;
  int m_iCNState;
  int m_iImputedCNState;
  float m_fCNConfidence;
  float m_fSmoothedLog2Ratio;
  float m_fLoh;
  float m_fSCAR;
  float m_fFLD;
  int m_iMaxPeaks;
  bool m_bUseForSketch;
  snp_distribution* m_pSnpDistribution;
  float m_fMuAA;
  float m_fMuAB;
  float m_fMuBB;
  float m_fInformation;
  int m_iUseInEMAlgorithm;
  float m_fMosaicismMixture;
  int m_iReplicateCount;
  unsigned int m_uiAllelePeaks1;
  unsigned int m_uiAllelePeaks2;
  bool m_bValidSCARExists;
  bool m_bValidFLDExists;
  bool m_bValidHomHetExists;
  bool m_bTrulySNP;
  bool m_bTrulyCN;

public:
  CNProbeSet();

  ~CNProbeSet();

  CNProbeSet(const CNProbeSet& that);

  const CNProbeSet& operator=(const CNProbeSet& that);

  void clear();

  void setProcess(bool b) { m_bProcess = b; }

  bool isProcess() { return m_bProcess; }

  // AffxString getProbeSetName() {return m_strProbeSetName;}
  std::string getProbeSetName() const { return m_strProbeSetName; }
  void setProbeSetName(const AffxString& str) { m_strProbeSetName = str; }
  unsigned char getChromosome() { return m_cChromosome; }
  void setChromosome(unsigned char c) { m_cChromosome = c; }
  int getPosition() { return m_iPosition; }
  void setPosition(int i) { m_iPosition = i; }
  float getMedianSignal() { return m_fMedianSignal; }
  void setMedianSignal(float f) { m_fMedianSignal = f; }
  float getXXMedianSignal() { return m_fXXMedianSignal; }
  void setXXMedianSignal(float f) { m_fXXMedianSignal = f; }
  float getYMedianSignal() { return m_fYMedianSignal; }
  void setYMedianSignal(float f) { m_fYMedianSignal = f; }
  float getAAMedianSignal() { return m_fAAMedianSignal; }
  void setAAMedianSignal(float f) { m_fAAMedianSignal = f; }
  float getABMedianSignal() { return m_fABMedianSignal; }
  void setABMedianSignal(float f) { m_fABMedianSignal = f; }
  float getBBMedianSignal() { return m_fBBMedianSignal; }
  void setBBMedianSignal(float f) { m_fBBMedianSignal = f; }
  char getProcessFlag() { return m_cProcessFlag; }
  void setProcessFlag(char c) { m_cProcessFlag = c; }
  void setTrulySNP(bool bInput) { m_bTrulySNP = bInput; }
  void setTrulyCN(bool bInput) { m_bTrulyCN = bInput; }
  bool processAsCNNormalize() { return m_bTrulyCN; }
  bool processAsSNPNormalize() { return m_bTrulySNP; }
  bool processAsCN() { return ((m_cProcessFlag == 1) || (m_cProcessFlag == 3)); }
  bool processAsSNP() { return ((m_cProcessFlag == 2) || (m_cProcessFlag == 3)); }
  bool processAsVisualization() { return ((m_cProcessFlag == 2) || (m_cProcessFlag == 3) || (m_cProcessFlag == 4)); }
  bool processAll() { return (m_cProcessFlag > 0); }
  bool isSty() { return m_cStyAdapterCode != -1; }
  bool isNsp() { return m_cNspAdapterCode != -1; }
  char getStyAdapterCode() { return m_cStyAdapterCode; }
  void setStyAdapterCode(char c) { m_cStyAdapterCode = c; }
  char getNspAdapterCode() { return m_cNspAdapterCode; }
  void setNspAdapterCode(char c) { m_cNspAdapterCode = c; }
  float getGCContent() { return m_fGCContent; }
  void setGCContent(float f) { m_fGCContent = f; }
  float getAMedianIntensity() { return m_fAMedianIntensity; }
  void setAMedianIntensity(float f) { m_fAMedianIntensity = f; }
  float getBMedianIntensity() { return m_fBMedianIntensity; }
  void setBMedianIntensity(float f) { m_fBMedianIntensity = f; }
  float getAAlleleSignal() { return m_fAAlleleSignal; }
  void setAAlleleSignal(float f) { m_fAAlleleSignal = f; }
  float getBAlleleSignal() { return m_fBAlleleSignal; }
  void setBAlleleSignal(float f) { m_fBAlleleSignal = f; }
  float getLog2Ratio() { return m_fLog2Ratio; }
  void setLog2Ratio(float f) { m_fLog2Ratio = f; }
  float getLog2RatioMedianSmooth() { return m_fLog2RatioMedianSmooth; }
  void setLog2RatioMedianSmooth(float f) { m_fLog2RatioMedianSmooth = f; }
  float getAllelicDifference() { return m_fAllelicDifference; }
  void setAllelicDifference(float f) { m_fAllelicDifference = f; }
  float getGcAdjustment() { return m_fGcAdjustment; }
  void setGcAdjustment(float f) { m_fGcAdjustment = f; }
  bool isPseudoAutosomalRegion() { return m_bPseudoAutosomalRegion; }
  void setPseudoAutosomalRegion(bool b) { m_bPseudoAutosomalRegion = b; }
  int getGCBinIndex() { return m_iGCBinIndex; }
  void setGCBinIndex(int i) { m_iGCBinIndex = i; }
  char getGenotypeCall() { return m_cGenotypeCall; }
  void setGenotypeCall(char c) { m_cGenotypeCall = c; }
  float getGenotypeConfidence() { return m_fGenotypeConfidence; }
  void setGenotypeConfidence(float f) { m_fGenotypeConfidence = f; }
  int getImputedCNState() { return m_iImputedCNState; }
  void setImputedCNState(int i) { m_iImputedCNState = i; }
  int getCNState() { return m_iCNState; }
  void setCNState(int i) { m_iCNState = i; }
  float getCNConfidence() { return m_fCNConfidence; }
  void setCNConfidence(float f) { m_fCNConfidence = f; }
  float getSmoothedLog2Ratio() { return m_fSmoothedLog2Ratio; }
  void setSmoothedLog2Ratio(float f) { m_fSmoothedLog2Ratio = f; }
  float getLoh() { return m_fLoh; }
  void setLoh(float f) { m_fLoh = f; }
  float getSCAR() { return m_fSCAR; }
  void setSCAR(float f) { m_fSCAR = f; }
  float getFLD() { return m_fFLD; }
  void setFLD(float f) { m_fFLD = f; }
  int getMaxPeaks() { return m_iMaxPeaks; }
  void setMaxPeaks(int i) { m_iMaxPeaks = i; }
  bool isUseForSketch() { return m_bUseForSketch; }
  void setUseForSketch(bool b) { m_bUseForSketch = b; }
  snp_distribution* getSnpDistribution() { return m_pSnpDistribution; }
  void setSnpDistribution(snp_distribution* p) { m_pSnpDistribution = p; }
  float getMosaicismMixture() { return m_fMosaicismMixture; }
  void setMosaicismMixture(float f) { m_fMosaicismMixture = f; }
  int getReplicateCount() { return m_iReplicateCount; }
  void setReplicateCount(int i) { m_iReplicateCount = i; }
  unsigned int getAllelePeaks1() { return m_uiAllelePeaks1; }
  void setAllelePeaks1(unsigned int ui) { m_uiAllelePeaks1 = ui; }
  unsigned int getAllelePeaks2() { return m_uiAllelePeaks2; }
  void setAllelePeaks2(unsigned int ui) { m_uiAllelePeaks2 = ui; }
  std::vector<float>& getWaves() { return m_vWaves; }
  float getWaveValue(int iWaveIndex)
  {
    if ((iWaveIndex < 0) || (iWaveIndex >= m_vWaves.size())) {
      return 0;
    } else {
      return m_vWaves[iWaveIndex];
    }
  }
  float getSignalEstimate() { return (m_fAAlleleSignal + m_fBAlleleSignal); }
  float getMuAA() { return m_fMuAA; }
  void setMuAA(float f) { m_fMuAA = f; }
  float getMuAB() { return m_fMuAB; }
  void setMuAB(float f) { m_fMuAB = f; }
  float getMuBB() { return m_fMuBB; }
  void setMuBB(float f) { m_fMuBB = f; }
  float getInformation() { return m_fInformation; }
  void setInformation(float input) { m_fInformation = input; }
  int getUseInEMAlgorithm() { return m_iUseInEMAlgorithm; }
  void setUseInEMAlgorithm(int input) { m_iUseInEMAlgorithm = input; }
  float getSignalContrast(double dK = 2.0);
  float getSignalContrastMvA();
  float getSignalStrength();
  float getSignalStrengthMvA();
  float getIntensityContrast(double dK = 2.0);
  float getIntensityStrength();
  float getCalibratedLog2Ratio(double dAlphaCNCalibrate, double dBetaCNCalibrate);
  float getCalibratedSmoothedLog2Ratio(double dAlphaCNCalibrate, double dBetaCNCalibrate);
  bool getValidSCARExists() { return m_bValidSCARExists; }
  void setValidSCARExists(bool b) { m_bValidSCARExists = b; }
  bool getValidFLDExists() { return m_bValidFLDExists; }
  void setValidFLDExists(bool b) { m_bValidFLDExists = b; }
  bool getValidHomHetExists() { return m_bValidHomHetExists; }
  void setValidHomHetExists(bool b) { m_bValidHomHetExists = b; }
  bool isUseForAllelePeaks() { return m_bValidSCARExists; }
  bool isUseForCyto2LOH() { return (m_bValidSCARExists && m_bValidHomHetExists) && (m_fInformation != 0.0) ; }
  bool lohCalled() { return (m_fLoh==0.0 || m_fLoh==1.0) ; }
  /**
   * Compare function.
   * @param that - A reference to an instance of this class.
   * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
   * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
  */
  int compareTo(CNProbeSet& that, int iCompareCode);
  /**
   * @brief Return the CNNeutral call (0 or 1)
   * @param int - The sex code for this sample.
   * @param int - The numeric value associated with the X chromosome
   * @param int - The numeric value associated with the Y chromosome
   * @return char - The call
   */
  char getCNNeutral(int iSexCode, int iXChromosome, int iYChromosome);
  /**
   * @brief Return the CNGain call (0 or 1)
   * @param int - The sex code for this sample.
   * @param int - The numeric value associated with the X chromosome
   * @param int - The numeric value associated with the Y chromosome
   * @return char - The call
   */
  char getCNGain(int iSexCode, int iXChromosome, int iYChromosome);
  /**
   * @brief Return the CNLoss call (0 or 1)
   * @param int - The sex code for this sample.
   * @param int - The numeric value associated with the X chromosome
   * @param int - The numeric value associated with the Y chromosome
   * @return char - The call
   */
  char getCNLoss(int iSexCode, int iXChromosome, int iYChromosome);
  /**
   * @brief Return the CNNeutralLOH call (0 or 1)
   * @param int - The sex code for this sample.
   * @param int - The numeric value associated with the X chromosome
   * @param int - The numeric value associated with the Y chromosome
   * @return char - The call
   */
  char getCNNeutralLOH(int iSexCode, int iXChromosome, int iYChromosome);
  /**
   * @brief Return the NormalDipload call (0 or 1)
   * @param int - The sex code for this sample.
   * @param int - The numeric value associated with the X chromosome
   * @param int - The numeric value associated with the Y chromosome
   * @return char - The call
   */
  char getNormalDiploid(int iSexCode, int iXChromosome, int iYChromosome);
  /**
   * @brief Return the specified segment type call
   * @param int - The segment type code
   * @param int - The sex code for this sample
   * @param int - The numeric value associated with the X chromosome
   * @param int - The numeric value associated with the Y chromosome
   * @return char - The call
   */
  char getSegmentTypeCall(int iSegmentType, int iSexCode, int iXChromosome, int iYChromosome);
  /**
   * @brief Return the specified segment type confidence
   * @param int - The segment type code
   * @return float - The comnfidence
   */
  float getSegmentTypeConfidence(int iSegmentType);
  /**
   * @brief Return the affx genotype call code
   * @return char - The call
   */
  char getGenotypeCallCode();

  void setCovariateValue(int index, float value);
  float getCovariateValue(int index) { return m_vAllCovariates[index]; }
  int getNumCovariates() { return m_vAllCovariates.size(); }

  template<int k> struct ComparePred {
      bool operator()(const CNProbeSet* lhs, const CNProbeSet* rhs) const {
          Err::errAbort("CNProbeSet: ComparePred instantiated with an invalid compare code = " + ToStr(k));
          return false;
      }
  };
};

template<> struct CNProbeSet::ComparePred<0> {
    bool operator()(const CNProbeSet* lhs, const CNProbeSet* rhs) const {
        return strcmp(lhs->m_strProbeSetName.c_str(), rhs->m_strProbeSetName.c_str()) < 0;
    }
};
template<> struct CNProbeSet::ComparePred<1> {
    bool operator()(const CNProbeSet* lhs, const CNProbeSet* rhs) const {
        int iCompareResult = AffxArray<char>::compare(lhs->m_cChromosome, rhs->m_cChromosome);
        if (iCompareResult == 0) {
            iCompareResult = AffxArray<int>::compare(lhs->m_iPosition, rhs->m_iPosition);
        }
        if (iCompareResult == 0) {
            iCompareResult = strcmp(lhs->m_strProbeSetName.c_str(), rhs->m_strProbeSetName.c_str());
        }
        return iCompareResult < 0;
    }
};

/**
 * @brief  A vector of ProbeSet Pointers.
 */
class CNProbeSetArray : public AffxArray<CNProbeSet>
{
public:
  CNProbeSetArray();

  ~CNProbeSetArray();

  bool isSketchSubset();

  int getProcessCount();
  int getCNProcessCount();
  void calculateLog2RatioMedianSmooth(int iWindowSize);

  void setupGCCorrectionBins(int iGCBinCount);
};

#endif


