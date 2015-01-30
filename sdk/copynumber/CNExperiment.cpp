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

#include "copynumber/CNExperiment.h"

int CNExperiment::m_iInstanceCount = 0;
AffxArray<AffxString> CNExperiment::m_arQCMetricColumnNames;


CNExperiment::CNExperiment()
{
  m_iIndex = 0;
  m_iInstanceCount++;
  m_eCNCallGender = affx::UnknownGender;
  m_eRawIntensityRatioGender = affx::UnknownGender;
  m_fCNCallGenderRatio = 0;
  m_fRawIntensityRatio = 0;
  m_bXX = false;
  m_bY = false;
  m_bCNCallGenderComputed = false;
  m_fXXRatio = 0;
  m_fYRatio = 0;
  m_fMedianAutosomeMedian = 0;
  m_fMadDiffCN = 0;
  m_fIqr = 0;
  m_fMeanAbsRle = 0;
  m_fGCCorrectionMetric = 0;
  m_fChrXMean = 0;
  m_fChrYMean = 0;
  m_fMedianCnState = 0;
  m_fHomFrequency = 0;
  m_fHetFrequency = 0;
  m_dPError=0.0;
  m_fCNCallGenderConfidence = 0;
  m_fMedianRawIntensity = 0;
  m_fCallRate = 0.0;
  m_fSNPQC = 0;
  m_fRawSNPQC = 0.0;
  m_bSNPQCset = false;
  m_fAntigenomicRatio = 0;
  m_fGenomeLOH = 0;
  m_fAutosomeGenomeLOH = 0;
  m_iNumberOfChromosomes = 0;
  m_iWavinessSegCountLoss = 0;
  m_iWavinessSegCountGain = 0;
  m_fWavinessSd = 0;
  m_iNumberOfWavesUsed = -1;
  m_iCNCallGenderNotZero = 0;
  m_iCNCallGenderYCount = 0;
  m_fL2Gradient = 0.0;
}

/**
 * @brief Destructor
 */
CNExperiment::~CNExperiment()
{
  m_arQCMetricColumnValues.deleteAll();
  m_iInstanceCount--;
  if (m_iInstanceCount == 0) {
    m_arQCMetricColumnNames.deleteAll();
  }

}
int CNExperiment::getIndex()
{
  return m_iIndex;
}
void CNExperiment::setIndex(int i)
{
  m_iIndex = i;
}

//MG temporary code for Bitao
AffxString CNExperiment::getExperimentNormalDiploidFileName()
{
  return m_strExperimentNormalDiploidFileName;
}
void CNExperiment::setExperimentNormalDiploidFileName(const AffxString& str)
{
  m_strExperimentNormalDiploidFileName = str;
}
//
AffxString CNExperiment::getExperimentName()
{
  return m_strExperimentName;
}
void CNExperiment::setExperimentName(const AffxString& str)
{
  m_strExperimentName = str;
}
int CNExperiment::getCNCallGenderNotZero()
{
  return m_iCNCallGenderNotZero;
}
void CNExperiment::setCNCallGenderNotZero(int i)
{
  m_iCNCallGenderNotZero = i;
}
int CNExperiment::getCNCallGenderYCount()
{
  return m_iCNCallGenderYCount;
}
void CNExperiment::setCNCallGenderYCount(int i)
{
  m_iCNCallGenderYCount = i;
}
bool CNExperiment::hasXX()
{
  return m_bXX;
}
void CNExperiment::setXX(bool b)
{
  m_bXX = b;
}
bool CNExperiment::getCNCallGenderComputed()
{
  return m_bCNCallGenderComputed;
}
void CNExperiment::setCNCallGenderComputed(bool b)
{
  m_bCNCallGenderComputed = b;
}
bool CNExperiment::hasY()
{
  return m_bY;
}
void CNExperiment::setY(bool b)
{
  m_bY = b;
}
float CNExperiment::getXXRatio()
{
  return m_fXXRatio;
}
void CNExperiment::setXXRatio(float f)
{
  m_fXXRatio = f;
}
float CNExperiment::getYRatio()
{
  return m_fYRatio;
}
void CNExperiment::setYRatio(float f)
{
  m_fYRatio = f;
}
float CNExperiment::getMedianAutosomeMedian()
{
  return m_fMedianAutosomeMedian;
}
float CNExperiment::getCNCallGenderRatio()
{
  return m_fCNCallGenderRatio;
}
void CNExperiment::setCNCallGenderRatio(float f)
{
  m_fCNCallGenderRatio = f;
}
float CNExperiment::getRawIntensityRatio()
{
  return m_fRawIntensityRatio;
}
void CNExperiment::setRawIntensityRatio(float f)
{
  m_fRawIntensityRatio = f;
}
void CNExperiment::setMedianAutosomeMedian(float f)
{
  m_fMedianAutosomeMedian = f;
}
float CNExperiment::getMadDiffCN()
{
  return m_fMadDiffCN;
}
void CNExperiment::setMadDiffCN(float f)
{
  m_fMadDiffCN = f;
}
float CNExperiment::getIqr()
{
  return m_fIqr;
}
void CNExperiment::setIqr(float f)
{
  m_fIqr = f;
}
float CNExperiment::getMeanAbsRle()
{
  return m_fMeanAbsRle;
}
void CNExperiment::setMeanAbsRle(float f)
{
  m_fMeanAbsRle = f;
}
float CNExperiment::getGCCorrectionMetric()
{
  return m_fGCCorrectionMetric;
}
void CNExperiment::setGCCorrectionMetric(float f)
{
  m_fGCCorrectionMetric = f;
}
AffxArray<AffxString>* CNExperiment::getQCMetricColumnNames()
{
  return &m_arQCMetricColumnNames;
}
AffxArray<AffxString>* CNExperiment::getQCMetricColumnValues()
{
  return &m_arQCMetricColumnValues;
}
float CNExperiment::getChrXMean()
{
  return m_fChrXMean;
}
void CNExperiment::setChrXMean(float f)
{
  m_fChrXMean = f;
}
float CNExperiment::getChrYMean()
{
  return m_fChrYMean;
}
void CNExperiment::setChrYMean(float f)
{
  m_fChrYMean = f;
}
float CNExperiment::getMedianCnState()
{
  return m_fMedianCnState;
}
void CNExperiment::setMedianCnState(float f)
{
  m_fMedianCnState = f;
}
float CNExperiment::getHomFrequency()
{
  return m_fHomFrequency;
}
void CNExperiment::setHomFrequency(float f)
{
  m_fHomFrequency = f;
}
float CNExperiment::getHetFrequency()
{
  return m_fHetFrequency;
}
void CNExperiment::setHetFrequency(float f)
{
  m_fHetFrequency = f;
}
// m_dPError is a value set on the command line for the CNAnalysisMethodLOH module.  It is also needed in the
// CNAnalysisMethodSetment in order to compute confidences.  Rather than have this value set on the commmand line
// for this latter module we pass it to the experiment object and access it when segmenting.
double CNExperiment::getPError()
{
    return m_dPError;
}
void CNExperiment::setPError( double d)
{
   m_dPError = d;
}
float CNExperiment::getCNCallGenderConfidence()
{
  return m_fCNCallGenderConfidence;
}
void CNExperiment::setCNCallGenderConfidence(float f)
{
  m_fCNCallGenderConfidence = f;
}
float CNExperiment::getMedianRawIntensity()
{
  return m_fMedianRawIntensity;
}
void CNExperiment::setMedianRawIntensity(float f)
{
  m_fMedianRawIntensity = f;
}
float CNExperiment::getCallRate()
{
    return m_fCallRate;
}
void CNExperiment::setCallRate(float f)
{
    m_fCallRate = f;
}
float CNExperiment::getSNPQC()
{
  return m_fSNPQC;
}
void CNExperiment::setSNPQC(float f)
{
  m_fSNPQC = f;
}
bool CNExperiment::isSNPQCset()
{
  return m_bSNPQCset;
}
void CNExperiment::setIsSNPQCset(bool b)
{
  m_bSNPQCset = b;
}
float CNExperiment::getRawSNPQC()
{
    return m_fRawSNPQC;
}
void CNExperiment::setRawSNPQC(float f)
{
    m_fRawSNPQC = f;
}
float CNExperiment::getAntigenomicRatio()
{
  return m_fAntigenomicRatio;
}
void CNExperiment::setAntigenomicRatio(float f)
{
  m_fAntigenomicRatio = f;
}
int CNExperiment::getWavinessSegCountLoss()
{
  return m_iWavinessSegCountLoss;
}
void CNExperiment::setWavinessSegCountLoss(int i)
{
  m_iWavinessSegCountLoss = i;
}
int CNExperiment::getWavinessSegCountGain()
{
  return m_iWavinessSegCountGain;
}
void CNExperiment::setWavinessSegCountGain(int i)
{
  m_iWavinessSegCountGain = i;
}
int CNExperiment::getWavinessSegCountTotal()
{
  if (m_iWavinessSegCountGain == -1 && m_iWavinessSegCountLoss == -1) {
    return -1;
  }
  return m_iWavinessSegCountGain + m_iWavinessSegCountLoss;
}
float CNExperiment::getWavinessSd()
{
  return m_fWavinessSd;
}
void CNExperiment::setWavinessSd(float f)
{
  m_fWavinessSd = f;
}
CNExperiment::wavAmplitudes_t CNExperiment::getWavinessAmplitudes()
{
  return m_vWavinessAmplitudes;
}
void CNExperiment::addWavinessAmplitude(int iWaveIndex, float wavAmplitude)
{
  m_vWavinessAmplitudes.push_back(std::make_pair(iWaveIndex, wavAmplitude));
}
std::pair<char, float> CNExperiment::getChromosomeLOH(int iIndex)
{
  return m_vChromosomeLOH[iIndex];
}
void CNExperiment::addChromosomeSummaryData(std::pair<char, float> data)
{
  m_vChromosomeLOH.push_back(data);
}
float CNExperiment::getGenomeLOH()
{
  return m_fGenomeLOH;
}
void CNExperiment::setGenomeLOH(float fGenomeLOH)
{
  m_fGenomeLOH = fGenomeLOH;
}
float CNExperiment::getAutosomeGenomeLOH()
{
    return m_fAutosomeGenomeLOH;
}
void CNExperiment::setAutosomeGenomeLOH(float fAutosomeGenomeLOH)
{
    m_fAutosomeGenomeLOH = fAutosomeGenomeLOH;
}
int CNExperiment::getNumberOfChromosomesToReport()
{
  return m_iNumberOfChromosomes;
}
void CNExperiment::setNumberOfChromosomesToReport(int iNumberOfChromosomes)
{
  m_iNumberOfChromosomes = iNumberOfChromosomes;
}


int CNExperiment::getNumberOfWavesUsed()
{
  return m_iNumberOfWavesUsed;
}

void CNExperiment::setNumberOfWavesUsed(int iNumberOfWavesUsed)
{
  m_iNumberOfWavesUsed = iNumberOfWavesUsed;
}

void CNExperiment::setCNReferenceHeader( const std::vector<std::string>& vInput)
{
  for(int iIndex=0; iIndex < vInput.size(); iIndex++)
  {
      m_vCNReferenceHeader.push_back(vInput[iIndex]);
  }
}

std::vector<std::string>* CNExperiment::getCNReferenceHeader ()
{
  return &m_vCNReferenceHeader;
}

/**
 * @brief Compare function.
 * @param that - A reference to an instance of this class.
 * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
 * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
 */
int CNExperiment::compareTo(CNExperiment& that, int iCompareCode)
{
  int iCompareResult = 0;
  switch (iCompareCode) {
  case 0:
    iCompareResult = m_strExperimentName.compareTo(that.m_strExperimentName, 0);
    break;
  }
  return iCompareResult;
}

/**
 * @brief Set the gender for the sample using the affx gender code
 * @param affx::Gender - The affx gender code
 */
void CNExperiment::setCNCallGender(affx::Gender eGender)
{
  m_eCNCallGender = eGender;
}
void CNExperiment::setCNCallGenderFromString(const AffxString& strIn)
{
  AffxString str = strIn;
  AffxString strTemp = str.toLowerCase();
  if (strTemp == "male") {
    setCNCallGender(affx::Male);
  } else if (strTemp == "female") {
    setCNCallGender(affx::Female);
  } else {
    setCNCallGender(affx::UnknownGender);
  }
}
void CNExperiment::setRawIntensityRatioGender(affx::Gender eGender)
{
  m_eRawIntensityRatioGender = eGender;
}

/**
 * @brief Set the gender for the sample using a string 'male' or 'female'
 * @param AffxString& - The gender string
 */
void CNExperiment::setRawIntensityRatioGenderFromString(const AffxString& strIn)
{
  AffxString str = strIn;
  AffxString strTemp = str.toLowerCase();
  if (strTemp == "male") {
    setRawIntensityRatioGender(affx::Male);
  } else if (strTemp == "female") {
    setRawIntensityRatioGender(affx::Female);
  } else {
    setRawIntensityRatioGender(affx::UnknownGender);
  }
}

AffxString CNExperiment::getRawIntensityRatioGender()
{
  return affx::getGenderString(m_eRawIntensityRatioGender);
}
AffxString CNExperiment::getCNCallGender()
{
  return affx::getGenderString(m_eCNCallGender);
}
int CNExperiment::getCNCallGenderAsInt()
{
  return (int)m_eCNCallGender;
}
int CNExperiment::getRawIntensityRatioGenderAsInt()
{
  return (int)m_eRawIntensityRatioGender;
}

float CNExperiment::getL2Gradient() { return m_fL2Gradient; }
void CNExperiment::setL2Gradient(float value) { m_fL2Gradient=value; }

CNExperimentArray::CNExperimentArray() {}

CNExperimentArray::~CNExperimentArray()
{
  deleteAll();
}

// This method was deliberately kept out of the CNExperiment class so as not
// to change the copy sematics of the class.
void CNExperimentArray::copyExperiments(CNExperimentArray& from, CNExperimentArray& to)
{
    to.reserve(from.size());
    for (int iexp = 0; iexp < from.getCount(); ++iexp) {
        CNExperiment* pexp = new CNExperiment();
        *pexp = *from.getAt(iexp);

        // deep copy vector of pointers (AffxArray<AffxString>)
        pexp->getQCMetricColumnValues()->resize(0);
        for (int istr = 0; istr < from.getAt(iexp)->getQCMetricColumnValues()->getCount(); ++istr) {
            pexp->getQCMetricColumnValues()->add(new AffxString(*from.getAt(iexp)->getQCMetricColumnValues()->getAt(istr)));
        }
        to.push_back(pexp);
    }
}
