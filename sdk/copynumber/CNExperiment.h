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

#ifndef _CNExperiment_H_
#define _CNExperiment_H_
/**
 * @file CNExperiment.h
 *
 * @brief This header contains the CNExperiment class definition.
 */

#include "chipstream/BioTypes.h"
#include "util/AffxArray.h"
#include "util/AffxString.h"
//
/**
 * @brief  A class for storing CNExperiment data.
 *
 */
class CNExperiment
{
public:
  typedef std::vector<std::pair<int, float> > wavAmplitudes_t;

private:
  int m_iIndex;
  //MG temporary code for Bitao
  AffxString m_strExperimentNormalDiploidFileName;
  //
  AffxString m_strExperimentName;
  affx::Gender m_eCNCallGender;
  affx::Gender m_eRawIntensityRatioGender;
  float m_fCNCallGenderRatio;
  float m_fRawIntensityRatio;
  bool m_bCNCallGenderComputed;
  bool m_bXX;
  bool m_bY;
  float m_fXXRatio;
  float m_fYRatio;
  float m_fMedianAutosomeMedian;
  float m_fMadDiffCN;
  float m_fIqr;
  float m_fMeanAbsRle;
  static int m_iInstanceCount;
  static AffxArray<AffxString> m_arQCMetricColumnNames;
  AffxArray<AffxString> m_arQCMetricColumnValues;
  float m_fGCCorrectionMetric;
  float m_fChrXMean;
  float m_fChrYMean;
  float m_fMedianCnState;
  float m_fHomFrequency;
  float m_fHetFrequency;
  double m_dPError;
  float m_fCNCallGenderConfidence;
  float m_fMedianRawIntensity;
  float m_fCallRate;
  float m_fSNPQC;
  float m_fRawSNPQC;
  bool m_bSNPQCset;
  float m_fAntigenomicRatio;
  std::vector<std::pair<char, float> > m_vChromosomeLOH;
  float m_fGenomeLOH;
  float m_fAutosomeGenomeLOH;
  int m_iNumberOfChromosomes;
  int m_iWavinessSegCountLoss;
  int m_iWavinessSegCountGain;
  float m_fWavinessSd;
  wavAmplitudes_t m_vWavinessAmplitudes;
  int m_iNumberOfWavesUsed;
  std::vector<std::string> m_vCNReferenceHeader;
  int m_iCNCallGenderNotZero;
  int m_iCNCallGenderYCount;
  float m_fL2Gradient;

public:
  /**
   * @brief Constructor
   */
  CNExperiment();

  /**
   * @brief Destructor
   */
  virtual ~CNExperiment();

  int getIndex();
  void setIndex(int i);
  //MG temporary code for Bitao
  AffxString getExperimentNormalDiploidFileName();
  void setExperimentNormalDiploidFileName(const AffxString& str);
  //
  int getCNCallGenderNotZero();
  void setCNCallGenderNotZero(int i);
  int getCNCallGenderYCount();
  void setCNCallGenderYCount(int i);
  bool getCNCallGenderComputed();
  void setCNCallGenderComputed(bool b);
  AffxString getExperimentName();
  void setExperimentName(const AffxString& str);
  bool hasXX();
  void setXX(bool b);
  bool hasY();
  void setY(bool b);
  float getXXRatio();
  void setXXRatio(float f);
  float getYRatio();
  void setYRatio(float f);
  float getMedianAutosomeMedian();
  void setMedianAutosomeMedian(float f);
  float getMadDiffCN();
  void setMadDiffCN(float f);
  float getIqr();
  void setIqr(float f);
  float getMeanAbsRle();
  void setMeanAbsRle(float f);
  float getGCCorrectionMetric();
  void setGCCorrectionMetric(float f);
  static AffxArray<AffxString>* getQCMetricColumnNames();
  AffxArray<AffxString>* getQCMetricColumnValues();
  float getChrXMean();
  void setChrXMean(float f);
  float getChrYMean();
  void setChrYMean(float f);
  float getCNCallGenderRatio();
  void setCNCallgenderRatio(float f);
  float getRawIntensityRatio();
  void setRawIntensityRatio(float f);
  float getMedianCnState();
  void setMedianCnState(float f);
  float getHomFrequency();
  void setHomFrequency(float f);
  float getHetFrequency();
  void setHetFrequency(float f);
  double getPError();
  void setPError(double d);
  float getCNCallGenderConfidence();
  void setCNCallGenderConfidence(float f);
  void setCNCallGenderRatio(float f);
  float getMedianRawIntensity();
  void setMedianRawIntensity(float f);
  float getCallRate();
  void setCallRate(float f);
  float getSNPQC();
  void setSNPQC(float f);
  bool isSNPQCset();
  void setIsSNPQCset(bool b);
  float getRawSNPQC();
  void setRawSNPQC(float f);
  float getAntigenomicRatio();
  void setAntigenomicRatio(float f);
  int getWavinessSegCountLoss();
  void setWavinessSegCountLoss(int i);
  int getWavinessSegCountGain();
  void setWavinessSegCountGain(int i);
  int getWavinessSegCountTotal();
  float getWavinessSd();
  void setWavinessSd(float f);
  wavAmplitudes_t getWavinessAmplitudes();
  void addWavinessAmplitude(int iWaveIndex, float wavAmplitude);
  std::pair<char, float> getChromosomeLOH(int iIndex);
  void addChromosomeSummaryData(std::pair<char, float> data);
  float getGenomeLOH();
  void setGenomeLOH(float fGenomeLOH);
  float getAutosomeGenomeLOH();
  void setAutosomeGenomeLOH(float fAutosomeGenomeLOH);
  int getNumberOfChromosomesToReport();
  void setNumberOfChromosomesToReport(int iNumberOfChromosomes);
  int getNumberOfWavesUsed();
  void setNumberOfWavesUsed(int iNumberOfWavesUsed);
  void setCNReferenceHeader(const std::vector<std::string>& vInput );
  std::vector <std::string>* getCNReferenceHeader();
  float getL2Gradient();
  void setL2Gradient(float value);

  /**
   * @brief Compare function.
   * @param that - A reference to an instance of this class.
   * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
   * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
   */
  int compareTo(CNExperiment& that, int iCompareCode);

  /**
   * @brief Set the gender for the sample using the affx gender code
   * @param affx::Gender - The affx gender code
   */
  void setCNCallGender(affx::Gender eGender);
  void setRawIntensityRatioGender(affx::Gender eGender);

  /**
   * @brief Set the gender for the sample using a string 'male' or 'female'
   * @param AffxString& - The gender string
   */
  void setRawIntensityRatioGenderFromString(const AffxString& strIn);
  void setCNCallGenderFromString(const AffxString& strIn);
  AffxString getRawIntensityRatioGender();
  AffxString getCNCallGender();
  int getCNCallGenderAsInt();
  int getRawIntensityRatioGenderAsInt();

  template<int k> struct ComparePred {
      bool operator()(const CNExperiment* lhs, const CNExperiment* rhs) const {
          Err::errAbort("CNExperiment: ComparePred instantiated with an invalid compare code = " + ToStr(k));
          return false;
      }
  };
};

template<> struct CNExperiment::ComparePred<0> {
    bool operator()(const CNExperiment* lhs, const CNExperiment* rhs) const {
        return lhs->m_strExperimentName.compareTo(rhs->m_strExperimentName, 0) < 0;
    }
};

/**
 * @brief  A vector of CNExperiment pointers
 */
class CNExperimentArray : public AffxArray<CNExperiment>
{
public:
  CNExperimentArray();
  ~CNExperimentArray();

public:
    static void copyExperiments(CNExperimentArray& from, CNExperimentArray& to);
};

#endif


