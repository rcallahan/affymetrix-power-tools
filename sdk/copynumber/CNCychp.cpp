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

//
#include "copynumber/CNCychp.h"
//
#include "util/Fs.h"
#include "util/Util.h"
//
#include "../external/xerces/src/xercesc/framework/XMLPScanToken.hpp"
#include "../external/xerces/src/xercesc/parsers/SAXParser.hpp"
#include "../external/xerces/src/xercesc/sax/AttributeList.hpp"
#include "../external/xerces/src/xercesc/sax/HandlerBase.hpp"
#include "../external/xerces/src/xercesc/util/OutOfMemoryException.hpp"
#include "../external/xerces/src/xercesc/util/PlatformUtils.hpp"
#include "../external/xerces/src/xercesc/util/XMLString.hpp"
#include "../external/xerces/src/xercesc/util/XMLUni.hpp"
#include "../external/xerces/src/xercesc/util/XMLUniDefs.hpp"
//
#include <limits>
//
#define XERCES_STATIC_LIBRARY
#define XML_LIBRARY

/**
 * Constructor.
 */
CNCychpHeader::CNCychpHeader()
{
  clear();
}

void CNCychpHeader::clear()
{
  m_strARRFileName.clear();
  m_strARRGuid.clear();
  m_strChpID.clear();
  m_ucXChromosome = 0;
  m_ucYChromosome = 0;
  m_eGender = affx::UnknownGender;
  m_fYTarget = 0;
  m_fMAPD = 0;
  m_fPVQC = 0;
  m_fCQC = -1;
  m_bXX = false;
  m_bY = false;
  m_iWavinessSegCount = -1;
}

void CNCychpHeader::setARRFileName(const AffxString& str)
{
  m_strARRFileName = str;
}
AffxString CNCychpHeader::getARRFileName()
{
  return m_strARRFileName;
}
void CNCychpHeader::setARRGuid(const AffxString& str)
{
  m_strARRGuid = str;
}
AffxString CNCychpHeader::getARRGuid()
{
  return m_strARRGuid;
}
void CNCychpHeader::setChpID(affymetrix_calvin_utilities::AffymetrixGuidType& str)
{
  m_strChpID = str;
}
affymetrix_calvin_utilities::AffymetrixGuidType CNCychpHeader::getChpID()
{
  return m_strChpID;
}
void CNCychpHeader::setXChromosome(unsigned char uc)
{
  m_ucXChromosome = uc;
}
unsigned char CNCychpHeader::getXChromosome()
{
  return m_ucXChromosome;
}
void CNCychpHeader::setYChromosome(unsigned char uc)
{
  m_ucYChromosome = uc;
}
unsigned char CNCychpHeader::getYChromosome()
{
  return m_ucYChromosome;
}
void CNCychpHeader::setGender(affx::Gender eGender)
{
  m_eGender = eGender;
}
affx::Gender CNCychpHeader::getGender()
{
  return m_eGender;
}
void CNCychpHeader::setYTarget(float f)
{
  m_fYTarget = f;
}
float CNCychpHeader::getYTarget()
{
  return m_fYTarget;
}
void CNCychpHeader::setMAPD(float f)
{
  m_fMAPD = f;
}
float CNCychpHeader::getMAPD()
{
  return m_fMAPD;
}
void CNCychpHeader::setPVQC(float f)
{
  m_fPVQC = f;
}
float CNCychpHeader::getPVQC()
{
  return m_fPVQC;
}
void CNCychpHeader::setCQC(float f)
{
  m_fCQC = f;
}
float CNCychpHeader::getCQC()
{
  return m_fCQC;
}
void CNCychpHeader::setXX(bool b)
{
  m_bXX = b;
}
bool CNCychpHeader::getXX()
{
  return m_bXX;
}
void CNCychpHeader::setY(bool b)
{
  m_bY = b;
}
bool CNCychpHeader::getY()
{
  return m_bY;
}
void CNCychpHeader::setWavinessSegCount(int i)
{
  m_iWavinessSegCount = i;
}
int CNCychpHeader::getWavinessSegCount()
{
  return m_iWavinessSegCount;
}
void CNCychpHeader::setCNReferenceFileName(const AffxString& str)
{
  m_strCNReferenceFileName = str;
}
AffxString CNCychpHeader::getCNReferenceFileName()
{
  return m_strCNReferenceFileName;
}

void CNCychpHeader::setGenderFromString(const AffxString& strIn)
{
  AffxString str = strIn;
  AffxString strTemp = str.toLowerCase();
  if (strTemp == "male") {
    setGender(affx::Male);
  } else if (strTemp == "female") {
    setGender(affx::Female);
  } else {
    setGender(affx::UnknownGender);
  }
}

AffxString CNCychpHeader::getGenderAsString()
{
  return affx::getGenderString(m_eGender);
}


/**
 * Constructor.
 */
CNCychpChromosomesSummary::CNCychpChromosomesSummary()
{
  m_ucChromosome = 0;
  m_uiStartIndex = 0;
  m_uiMarkerCount = 0;
  m_fMinSignal = 0;
  m_fMaxSignal = 0;
  m_fMedianCnState = 0;
  m_fHomFrequency = 0;
  m_fHetFrequency = 0;
  m_fMosaicism = 0;
  m_fLOH = 0;
}

void CNCychpChromosomesSummary::setChromosome(unsigned char uc)
{
  m_ucChromosome = uc;
}
unsigned char CNCychpChromosomesSummary::getChromosome()
{
  return m_ucChromosome;
}
void CNCychpChromosomesSummary::setChromosomeDisplay(const AffxString& str)
{
  m_strChromosomeDisplay = str;
}
AffxString CNCychpChromosomesSummary::getChromosomeDisplay()
{
  return m_strChromosomeDisplay;
}
void CNCychpChromosomesSummary::setStartIndex(unsigned int ui)
{
  m_uiStartIndex = ui;
}
unsigned int CNCychpChromosomesSummary::getStartIndex()
{
  return m_uiStartIndex;
}
void CNCychpChromosomesSummary::setMarkerCount(unsigned int ui)
{
  m_uiMarkerCount = ui;
}
unsigned int CNCychpChromosomesSummary::getMarkerCount()
{
  return m_uiMarkerCount;
}
void CNCychpChromosomesSummary::setMinSignal(float f)
{
  m_fMinSignal = f;
}
float CNCychpChromosomesSummary::getMinSignal()
{
  return m_fMinSignal;
}
void CNCychpChromosomesSummary::setMaxSignal(float f)
{
  m_fMaxSignal = f;
}
float CNCychpChromosomesSummary::getMaxSignal()
{
  return m_fMaxSignal;
}
void CNCychpChromosomesSummary::setMedianCnState(float f)
{
  m_fMedianCnState = f;
}
float CNCychpChromosomesSummary::getMedianCnState()
{
  return m_fMedianCnState;
}
void CNCychpChromosomesSummary::setHomFrequency(float f)
{
  m_fHomFrequency = f;
}
float CNCychpChromosomesSummary::getHomFrequency()
{
  return m_fHomFrequency;
}
void CNCychpChromosomesSummary::setHetFrequency(float f)
{
  m_fHetFrequency = f;
}
float CNCychpChromosomesSummary::getHetFrequency()
{
  return m_fHetFrequency;
}
void CNCychpChromosomesSummary::setMosaicism(float f)
{
  m_fMosaicism = f;
}
float CNCychpChromosomesSummary::getMosaicism()
{
  return m_fMosaicism;
}
void CNCychpChromosomesSummary::setLOH(float f)
{
  m_fLOH = f;
}
float CNCychpChromosomesSummary::getLOH()
{
  return m_fLOH;
}

/**
 * Compare function.
 * @param that - A reference to an instance of this class.
 * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
 * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
 */
int CNCychpChromosomesSummary::compareTo(CNCychpChromosomesSummary& that, int iCompareCode)
{
  int iCompareResult = 0;
  switch (iCompareCode) {
  case 0:
    iCompareResult = AffxArray<unsigned char>::compare(m_ucChromosome, that.m_ucChromosome);
    break;
  }
  return iCompareResult;
}

/**
 * Constructor
*/
CNCychpProbeSetsCopyNumber::CNCychpProbeSetsCopyNumber()
{
  m_ucChromosome = 0;
  m_uiPosition = 0;
  m_fLog2Ratio = 0;
  m_fSmoothSignal = 0;
  m_fAAlleleSignal = 0;
  m_fBAlleleSignal = 0;
  m_fSCAR = 0;
  m_cGenotypeCall = -1;
  m_bIsGenotypeCallSet = false;
  m_fGenotypeConfidence = 0;
  m_ucCnCall = 0;
  m_ucLohCall = 0;
  m_cSegmentInput = -1;
  m_cSegmentOutput = -1;
}

CNCychpChromosomesSummaries::~CNCychpChromosomesSummaries()
{
  clear();
}
void CNCychpChromosomesSummaries::clear()
{
  deleteAll();
}
void CNCychpChromosomesSummaries::newCychpChromosomesSummary(affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex)
{
  CNCychpChromosomesSummary* p = new CNCychpChromosomesSummary;
  u_int8_t uc = 0;
  u_int32_t ui = 0;
  float f = 0;
  AffxString str;
  pDataSet->GetData((int)uiRowIndex, 0, uc);
  p->setChromosome(uc);
  pDataSet->GetData((int)uiRowIndex, 1, str);
  p->setChromosomeDisplay(str);
  pDataSet->GetData((int)uiRowIndex, 2, ui);
  p->setStartIndex(ui);
  pDataSet->GetData((int)uiRowIndex, 3, ui);
  p->setMarkerCount(ui);
  pDataSet->GetData((int)uiRowIndex, 4, f);
  p->setMinSignal(f);
  pDataSet->GetData((int)uiRowIndex, 5, f);
  p->setMaxSignal(f);
  pDataSet->GetData((int)uiRowIndex, 6, f);
  p->setMedianCnState(f);
  pDataSet->GetData((int)uiRowIndex, 7, f);
  p->setHomFrequency(f);
  pDataSet->GetData((int)uiRowIndex, 8, f);
  p->setHetFrequency(f);
  pDataSet->GetData((int)uiRowIndex, 9, f);
  p->setMosaicism(f);
  pDataSet->GetData((int)uiRowIndex, 10, f);
  p->setLOH(f);
  add(p);
}

/**
 * @brief Set the allele code ffor the genotype call from the internal representation (-1, 0, 1, 2)
 * @param char - The internal representation
 */
void CNCychpProbeSetsCopyNumber::setGenotypeCallFromCode(char c)
{
  switch (c) {
  case ALLELE_NO_CALL: m_cGenotypeCall = -1; break;
  case ALLELE_A_CALL: m_cGenotypeCall = 0; break;
  case ALLELE_AB_CALL: m_cGenotypeCall = 1; break;
  case ALLELE_B_CALL: m_cGenotypeCall = 2; break;
  default: m_cGenotypeCall = -1;
  }
}

/**
 * Compare function.
 * @param that - A reference to an instance of this class.
 * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
 * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
 */
int CNCychpProbeSetsCopyNumber::compareTo(CNCychpProbeSetsCopyNumber& that, int iCompareCode)
{
  int iCompareResult = 0;
  switch (iCompareCode) {
  case 0:
    iCompareResult = m_strProbeSetName.compareTo(that.m_strProbeSetName, 0);
    break;
  case 1:
    iCompareResult = AffxArray<unsigned char>::compare(m_ucChromosome, that.m_ucChromosome);
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<unsigned int>::compare(m_uiPosition, that.m_uiPosition);
    }
    if (iCompareResult == 0) {
      iCompareResult = m_strProbeSetName.compareTo(that.m_strProbeSetName, 0);
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
char CNCychpProbeSetsCopyNumber::getCNNeutral(int iSexCode, int iXChromosome, int iYChromosome)
{
  return (((getCNGain(iSexCode, iXChromosome, iYChromosome) == 0) && (getCNLoss(iSexCode, iXChromosome, iYChromosome) == 0)) ? (char)1 : (char)0);
}

/**
 * @brief Return the CNGain call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNCychpProbeSetsCopyNumber::getCNGain(int iSexCode, int iXChromosome, int iYChromosome)
{
  char cCnGain = 0;
  if (m_ucChromosome == iXChromosome) {
    if (iSexCode == affx::Male) {
      if (m_ucCnCall > 1) {
        cCnGain = 1;
      }
    } else {
      if (m_ucCnCall > 2) {
        cCnGain = 1;
      }
    }
  } else if (m_ucChromosome == iYChromosome) {
    if (iSexCode == affx::Male) {
      if (m_ucCnCall > 1) {
        cCnGain = 1;
      }
    } else {
      if (m_ucCnCall > 0) {
        cCnGain = 1;
      }
    }
  } else {
    if (m_ucCnCall > 2) {
      cCnGain = 1;
    }
  }
  return cCnGain;
}

/**
 * @brief Return the CNLoss call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNCychpProbeSetsCopyNumber::getCNLoss(int iSexCode, int iXChromosome, int iYChromosome)
{
  char cCnLoss = 0;
  if (m_ucChromosome == iXChromosome) {
    if (iSexCode == affx::Male) {
      if (m_ucCnCall < 1) {
        cCnLoss = 1;
      }
    } else {
      if (m_ucCnCall < 2) {
        cCnLoss = 1;
      }
    }
  } else if (m_ucChromosome == iYChromosome) {
    if (iSexCode == affx::Male) {
      if (m_ucCnCall < 1) {
        cCnLoss = 1;
      }
    } else {
      cCnLoss = 0;
      //                if (m_ucCnCall < 0) {cCnNLoss = 1;}
    }
  } else {
    if (m_ucCnCall < 2) {
      cCnLoss = 1;
    }
  }
  return cCnLoss;
}

/**
 * @brief Return the CNNeutralLOH call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNCychpProbeSetsCopyNumber::getCNNeutralLoh(int iSexCode, int iXChromosome, int iYChromosome)
{
  char cCnNeutral = getCNNeutral(iSexCode, iXChromosome, iYChromosome);
  return (((cCnNeutral == 1) && (m_ucLohCall == 1)) ? (char)1 : (char)0);
}

/**
 * @brief Return the NormalDiploid call (0 or 1)
 * @param int - The sex code for this sample.
 * @param int - The numeric value associated with the X chromosome
 * @param int - The numeric value associated with the Y chromosome
 * @return char - The call
 */
char CNCychpProbeSetsCopyNumber::getNormalDiploid(int iSexCode, int iXChromosome, int iYChromosome)
{
  char cCnNeutral = getCNNeutral(iSexCode, iXChromosome, iYChromosome);
  return (((cCnNeutral == 1) && (m_ucLohCall == 0)) ? (char)1 : (char)0);
}


CNCychpProbeSetsCopyNumbers::~CNCychpProbeSetsCopyNumbers()
{
  clear();
}
void CNCychpProbeSetsCopyNumbers::clear()
{
  deleteAll();
}

void CNCychpProbeSetsCopyNumbers::newCychpProbeSetsCopyNumber(affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex)
{
  CNCychpProbeSetsCopyNumber* p = new CNCychpProbeSetsCopyNumber;
  u_int8_t uc = 0;
  u_int32_t ui = 0;
  float f = 0;
  AffxString str;
  pDataSet->GetData((int)uiRowIndex, 0, str);
  p->setProbeSetName(str);
  pDataSet->GetData((int)uiRowIndex, 1, uc);
  p->setChromosome(uc);
  pDataSet->GetData((int)uiRowIndex, 2, ui);
  p->setPosition(ui);
  pDataSet->GetData((int)uiRowIndex, 3, f);
  p->setLog2Ratio(f);
  pDataSet->GetData((int)uiRowIndex, 5, f);
  p->setSmoothSignal(f);
  add(p);
}

void CNCychpProbeSetsCopyNumbers::setCychpAlgorithmDataMarkerABSignal(affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex)
{
  u_int32_t ui = 0;
  float f = 0;
  pDataSet->GetData((int)uiRowIndex, 0, ui);
  CNCychpProbeSetsCopyNumber* p = getAt(ui);
  pDataSet->GetData((int)uiRowIndex, 1, f);
  p->setBAlleleSignal(f);
  pDataSet->GetData((int)uiRowIndex, 2, f);
  p->setAAlleleSignal(f);
  pDataSet->GetData((int)uiRowIndex, 3, f);
  p->setSCAR(f);
}

void CNCychpProbeSetsCopyNumbers::setCychpGenotypingCalls(affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex)
{
  u_int32_t ui = 0;
  int8_t c = 0;
  pDataSet->GetData((int)uiRowIndex, 0, ui);
  CNCychpProbeSetsCopyNumber* p = getAt(ui);
  pDataSet->GetData((int)uiRowIndex, 3, c);    // ForcedCall column
  p->setGenotypeCallFromCode(c);
  p->setIsGenotypeCallSet(true);
}


/**
 * Constructor
 */
CNCychpSegment::CNCychpSegment()
{
  m_uiDataSetIndex = 0;
  m_uiSegmentID = 0;
  m_ucChromosome = 0;
  m_uiStartPosition = 0;
  m_uiStopPosition = 0;
  m_iMarkerCount = 0;
  m_uiMeanMarkerDistance = 0;
  m_fState = 0;
  m_fConfidence = 0;
  m_fMixture = 0;
}

void CNCychpSegment::setDataSetIndex(unsigned int ui)
{
  m_uiDataSetIndex = ui;
}
unsigned int CNCychpSegment::getDataSetIndex()
{
  return m_uiDataSetIndex;
}
void CNCychpSegment::setSegmentID(unsigned int ui)
{
  m_uiSegmentID = ui;
}
unsigned int CNCychpSegment::getSegmentID()
{
  return m_uiSegmentID;
}
void CNCychpSegment::setChromosome(unsigned char uc)
{
  m_ucChromosome = uc;
}
unsigned char CNCychpSegment::getChromosome()
{
  return m_ucChromosome;
}
void CNCychpSegment::setStartPosition(unsigned int ui)
{
  m_uiStartPosition = ui;
}
unsigned int CNCychpSegment::getStartPosition()
{
  return m_uiStartPosition;
}
void CNCychpSegment::setStopPosition(unsigned int ui)
{
  m_uiStopPosition = ui;
}
unsigned int CNCychpSegment::getStopPosition()
{
  return m_uiStopPosition;
}
void CNCychpSegment::setMarkerCount(int i)
{
  m_iMarkerCount = i;
}
int CNCychpSegment::getMarkerCount()
{
  return m_iMarkerCount;
}
void CNCychpSegment::setMeanMarkerDistance(unsigned int ui)
{
  m_uiMeanMarkerDistance = ui;
}
unsigned int CNCychpSegment::getMeanMarkerDistance()
{
  return m_uiMeanMarkerDistance;
}
void CNCychpSegment::setFState(float f)
{
  m_fState = f;
}
float CNCychpSegment::getFState()
{
  return m_fState;
}
void CNCychpSegment::setConfidence(float f)
{
  m_fConfidence = f;
}
float CNCychpSegment::getConfidence()
{
  return m_fConfidence;
}
void CNCychpSegment::setMixture(float f)
{
  m_fMixture = f;
}
float CNCychpSegment::getMixture()
{
  return m_fMixture;
}

/**
 * @brief Return the string associated with the segment type for this object.
 * @return AffxString - The segment type string.
 */
AffxString CNCychpSegment::getSegmentType()
{
  switch (m_uiDataSetIndex) {
  case 0: return "CN";
  case 1: return "LOH";
  case 2: return "CNNeutralLOH";
  case 3: return "NormalDiploid";
  case 4: return "Mosaicism";
  default: return "Unknown";
  }
  return "Unknown";
}

/**
 * Compare function.
 * @param that - A reference to an instance of this class.
 * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
 * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
 */
int CNCychpSegment::compareTo(CNCychpSegment& that, int iCompareCode)
{
  int iCompareResult = 0;
  switch (iCompareCode) {
  case 0:
    iCompareResult = AffxArray<unsigned int>::compare(m_uiDataSetIndex, that.m_uiDataSetIndex);
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<unsigned int>::compare(m_uiSegmentID, that.m_uiSegmentID);
    }
    break;
  case 1:
    iCompareResult = AffxArray<unsigned int>::compare(m_uiDataSetIndex, that.m_uiDataSetIndex);
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<unsigned char>::compare(m_ucChromosome, that.m_ucChromosome);
    }
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<unsigned int>::compare(m_uiStartPosition, that.m_uiStartPosition);
    }
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<unsigned int>::compare(m_uiStopPosition, that.m_uiStopPosition);
    }
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<unsigned int>::compare(m_uiSegmentID, that.m_uiSegmentID);
    }
    break;
  }
  return iCompareResult;
}


CNCychpSegments::~CNCychpSegments()
{
  clear();
}
void CNCychpSegments::clear()
{
  deleteAll();
}
void CNCychpSegments::newCychpSegmentCN(unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex)
{
  CNCychpSegment* p = new CNCychpSegment;
  u_int8_t uc = 0;
  u_int32_t ui = 0;
  float f = 0;
  int32_t i = 0;
  AffxString str;
  p->setDataSetIndex(uiDataSetIndex);
  pDataSet->GetData((int)uiRowIndex, 0, ui);
  p->setSegmentID(ui);
  pDataSet->GetData((int)uiRowIndex, 1, uc);
  p->setChromosome(uc);
  pDataSet->GetData((int)uiRowIndex, 2, ui);
  p->setStartPosition(ui);
  pDataSet->GetData((int)uiRowIndex, 3, ui);
  p->setStopPosition(ui);
  pDataSet->GetData((int)uiRowIndex, 4, i);
  p->setMarkerCount(i);
  pDataSet->GetData((int)uiRowIndex, 5, ui);
  p->setMeanMarkerDistance(ui);
  pDataSet->GetData((int)uiRowIndex, 6, f);
  p->setFState(f);
  pDataSet->GetData((int)uiRowIndex, 7, f);
  p->setConfidence(f);
  add(p);
}
void CNCychpSegments::newCychpSegment(unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex)
{
  CNCychpSegment* p = new CNCychpSegment;
  u_int8_t uc = 0;
  u_int32_t ui = 0;
  float f = 0;
  int32_t i = 0;
  AffxString str;
  p->setDataSetIndex(uiDataSetIndex);
  pDataSet->GetData((int)uiRowIndex, 0, ui);
  p->setSegmentID(ui);
  pDataSet->GetData((int)uiRowIndex, 1, uc);
  p->setChromosome(uc);
  pDataSet->GetData((int)uiRowIndex, 2, ui);
  p->setStartPosition(ui);
  pDataSet->GetData((int)uiRowIndex, 3, ui);
  p->setStopPosition(ui);
  pDataSet->GetData((int)uiRowIndex, 4, i);
  p->setMarkerCount(i);
  pDataSet->GetData((int)uiRowIndex, 5, ui);
  p->setMeanMarkerDistance(ui);
  pDataSet->GetData((int)uiRowIndex, 6, uc);
  p->setFState(uc);
  pDataSet->GetData((int)uiRowIndex, 7, f);
  p->setConfidence(f);
  add(p);
}
void CNCychpSegments::newCychpSegmentMosaicism(unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex)
{
  CNCychpSegment* p = new CNCychpSegment;
  u_int8_t uc = 0;
  u_int32_t ui = 0;
  float f = 0;
  int32_t i = 0;
  AffxString str;
  p->setDataSetIndex(uiDataSetIndex);
  pDataSet->GetData((int)uiRowIndex, 0, ui);
  p->setSegmentID(ui);
  pDataSet->GetData((int)uiRowIndex, 1, uc);
  p->setChromosome(uc);
  pDataSet->GetData((int)uiRowIndex, 2, ui);
  p->setStartPosition(ui);
  pDataSet->GetData((int)uiRowIndex, 3, ui);
  p->setStopPosition(ui);
  pDataSet->GetData((int)uiRowIndex, 4, i);
  p->setMarkerCount(i);
  pDataSet->GetData((int)uiRowIndex, 5, ui); p->setMeanMarkerDistance(ui);
  //        pDataSet->GetData((int)uiRowIndex, 6, uc); p->setFState(uc);
  pDataSet->GetData((int)uiRowIndex, 7, f);
  p->setConfidence(f);
  pDataSet->GetData((int)uiRowIndex, 8, f);
  p->setFState(f);
  pDataSet->GetData((int)uiRowIndex, 9, f);
  p->setMixture(f);
  add(p);
}


/**
 * For Mendelian Inheritance Error (MIE) results
 */
std::map<unsigned char, unsigned int> CNCychp::m_mMarkerCount;
std::map<unsigned char, int> CNCychp::m_mMIE_Trio;
std::map<unsigned char, int> CNCychp::m_mMIE_Mat;
std::map<unsigned char, int> CNCychp::m_mMIE_Pat;


/**
 * CNCychp constructor.
 */
CNCychp::CNCychp()
{
  clear();
}

/**
 * CNCychp destructor.
 */
CNCychp::~CNCychp()
{
  clear();
}

/**
 * Initialize the CNCychp data and free any memory allocated.
 */
void CNCychp::clear()
{
  m_strFileName.clear();
  m_strFamilialType.clear();
  m_ucFamilialCall = 0;
  m_fFamilialConfidence = 0;

  getCychpHeader().clear();
  getCychpChromosomesSummaries().clear();
  getCychpProbeSetsCopyNumbers().clear();
  getCychpSegments().clear();
}

/**
 * Warning handling for the CNCychp class.
 * @param strMessage - The warning messagae to report.
 */
void CNCychp::warning(const AffxString& strMessage)
{
  Verbose::out(1, "WARNING: " + strMessage);
}

/**
 * Error handling for the CNCychp class.
 * @param strMessage - The error messagae to report.
 */
void CNCychp::error(const AffxString& strMessage)
{
  clear();
  Err::errAbort(strMessage);
}

AffxString CNCychp::getFileName()
{
  return m_strFileName;
}
void CNCychp::setFileName(const AffxString& str)
{
  m_strFileName = Fs::normalizePath(str);
}

AffxString CNCychp::getARRGuid()
{
  return m_strARRGuid;
}
void CNCychp::setARRGuid(const AffxString& str)
{
  m_strARRGuid = str;
}

CNCychpHeader& CNCychp::getCychpHeader()
{
  return m_objCNCychpHeader;
}
CNCychpChromosomesSummaries& CNCychp::getCychpChromosomesSummaries()
{
  return m_vCychpChromosomesSummaries;
}
CNCychpProbeSetsCopyNumbers& CNCychp::getCychpProbeSetsCopyNumbers()
{
  return m_vCychpProbeSetsCopyNumbers;
}
CNCychpSegments& CNCychp::getCychpSegments()
{
  return m_vCNCychpSegments;
}

AffxString& CNCychp::getFamilialType()
{
  return m_strFamilialType;
}
void CNCychp::setFamilialType(const AffxString& str)
{
  m_strFamilialType = str;
}

unsigned char CNCychp::getFamilialCall()
{
  return m_ucFamilialCall;
}
void CNCychp::setFamilialCall(unsigned char uc)
{
  m_ucFamilialCall = uc;
}

float CNCychp::getFamilialConfidence()
{
  return m_fFamilialConfidence;
}
void CNCychp::setFamilialConfidence(float f)
{
  m_fFamilialConfidence = f;
}

unsigned char CNCychp::getFamilialDuoCall()
{
  return m_ucFamilialDuoCall;
}
void CNCychp::setFamilialDuoCall(unsigned char uc)
{
  m_ucFamilialDuoCall = uc;
}

float CNCychp::getFamilialDuoConfidence()
{
  return m_fFamilialDuoConfidence;
}
void CNCychp::setFamilialDuoConfidence(float f)
{
  m_fFamilialDuoConfidence = f;
}

AffxString CNCychp::getSelectedQC()
{
    return m_selectedQC;
}

void CNCychp::setSelectedQC(const AffxString& str)
{
    m_selectedQC = str;
}

/**
 * Read the CNCychp file, and store the data in memory.
 * @param strFileName - The name of the CNCychp file to load data from.
 * @return - bool value. (true if successful)
 */
bool CNCychp::readFile(const AffxString& strFileName, AffxArray<CNCychpProbeSetsCopyNumber>& vProbeSets, bool bLoadProbeSetName, bool bNonNaNOnly)
{
  AffxString strFamilialType = "";
  AffxString strPrompt = "";
  bool bSuccessful = false;
  clear();
  m_strFamilialType = strFamilialType;
  affymetrix_calvin_io::GenericFileReader reader;
  affymetrix_calvin_io::GenericData genericData;
  try {
    m_strFileName = Fs::convertToUncPath(strFileName);
//        Verbose::out(1, "*");
    if (strPrompt == "") {
      Verbose::out(1, "Reading CNCychp file: " + m_strFileName);
    } else {
      Verbose::out(1, "Reading " + strPrompt + ": " + m_strFileName);
    }
    reader.SetFilename(m_strFileName);
    reader.Open(genericData);
    loadHeader(genericData.Header().GetGenericDataHdr());
    WStringVector vDataGroupNames;
    genericData.DataGroupNames(vDataGroupNames);
    for (unsigned int uiDataGroupIndex = 0; (uiDataGroupIndex < vDataGroupNames.size()); uiDataGroupIndex++) {
      AffxString strDataGroupName = StringUtils::ConvertWCSToMBS(vDataGroupNames[uiDataGroupIndex]);
//            Verbose::out(1, " DataGroupName = " + strDataGroupName);
      WStringVector vDataSetNames;
      genericData.DataSetNames(uiDataGroupIndex, vDataSetNames);
      for (unsigned int uiDataSetIndex = 0; (uiDataSetIndex < vDataSetNames.size()); uiDataSetIndex++) {
        affymetrix_calvin_io::DataSet* pDataSet = genericData.DataSet(uiDataGroupIndex, uiDataSetIndex);
        AffxString strDataSetName = StringUtils::ConvertWCSToMBS(vDataSetNames[uiDataSetIndex]);
//                Verbose::out(1, "  DataSetName = " + strDataSetName);
        loadDataSet(strDataGroupName, strDataSetName, uiDataSetIndex, pDataSet, vProbeSets, bLoadProbeSetName, bNonNaNOnly);
      }
      setCnCalls();
      setLohCalls();
//            Verbose::out(1, "*");
    }
//        Verbose::out(1, "*");
    genericData.Clear();
    bSuccessful = true;
  } catch (...) {
    genericData.Clear(); error("Cannot read CNCychp file: " + m_strFileName);
  }
  return bSuccessful;
}

/**
 * Read the CNCychp file, and store the data in memory.
 * @param strFileName - The name of the CNCychp file to load data from.
 * @return - bool value. (true if successful)
 */
bool CNCychp::readFile(const AffxString& strFileName, const AffxString& strFamilialType, const AffxString& strPrompt)
{
  bool bSuccessful = false;
  clear();
  m_strFamilialType = strFamilialType;
  affymetrix_calvin_io::GenericFileReader reader;
  affymetrix_calvin_io::GenericData genericData;
  try {
    m_strFileName = Fs::convertToUncPath(strFileName);
//        Verbose::out(1, "*");
    if (strPrompt == "") {
      Verbose::out(1, "Reading CNCychp file: " + m_strFileName);
    } else {
      Verbose::out(1, "Reading " + strPrompt + ": " + m_strFileName);
    }
    reader.SetFilename(m_strFileName);
    reader.Open(genericData);
    loadHeader(genericData.Header().GetGenericDataHdr());
    WStringVector vDataGroupNames;
    genericData.DataGroupNames(vDataGroupNames);
    for (unsigned int uiDataGroupIndex = 0; (uiDataGroupIndex < vDataGroupNames.size()); uiDataGroupIndex++) {
      AffxString strDataGroupName = StringUtils::ConvertWCSToMBS(vDataGroupNames[uiDataGroupIndex]);
//            Verbose::out(1, " DataGroupName = " + strDataGroupName);
      WStringVector vDataSetNames;
      genericData.DataSetNames(uiDataGroupIndex, vDataSetNames);
      for (unsigned int uiDataSetIndex = 0; (uiDataSetIndex < vDataSetNames.size()); uiDataSetIndex++) {
        affymetrix_calvin_io::DataSet* pDataSet = genericData.DataSet(uiDataGroupIndex, uiDataSetIndex);
        AffxString strDataSetName = StringUtils::ConvertWCSToMBS(vDataSetNames[uiDataSetIndex]);
//                Verbose::out(1, "  DataSetName = " + strDataSetName);
        loadDataSet(strDataGroupName, strDataSetName, uiDataSetIndex, pDataSet);
      }
      setCnCalls();
      setLohCalls();
//            Verbose::out(1, "*");
    }
//        Verbose::out(1, "*");
    genericData.Clear();
    bSuccessful = true;
  } catch (...) {
    genericData.Clear(); error("Cannot read CNCychp file: " + m_strFileName);
  }
  return bSuccessful;
}

/**
 * Selectively read the CNCychp file, and store the data in memory.
 * @param groupDatasets - map containing calvin group nad dataset names to be read
 * @param strFileName - The name of the CNCychp file to load data from.
 * @return - bool value. (true if successful)
 */
bool CNCychp::readFile(const AffxString& strFileName, const groupDatasets_t& groupDatasets, const AffxString& strFamilialType, const AffxString& strPrompt)
{
  bool bSuccessful = false;
  clear();
  m_strFamilialType = strFamilialType;
  affymetrix_calvin_io::GenericFileReader reader;
  affymetrix_calvin_io::GenericData genericData;
  try {
    m_strFileName = Fs::convertToUncPath(strFileName);
//        Verbose::out(1, "*");
    if (strPrompt == "") {
      Verbose::out(1, "Reading CNCychp file: " + m_strFileName);
    } else {
      Verbose::out(1, "Reading " + strPrompt + ": " + m_strFileName);
    }
    reader.SetFilename(m_strFileName);
    reader.Open(genericData);
    loadHeader(genericData.Header().GetGenericDataHdr());
    WStringVector vDataGroupNames;
    genericData.DataGroupNames(vDataGroupNames);
    for (unsigned int uiDataGroupIndex = 0; (uiDataGroupIndex < vDataGroupNames.size()); uiDataGroupIndex++) {
      AffxString strDataGroupName = StringUtils::ConvertWCSToMBS(vDataGroupNames[uiDataGroupIndex]);
      if (!isGroupPresent(groupDatasets, strDataGroupName)) {
        continue;
      }
      if (isGroupRedundant(strDataGroupName)) {
        continue;
      }
      //Verbose::out(1, " DataGroupName = " + strDataGroupName);
      WStringVector vDataSetNames;
      genericData.DataSetNames(uiDataGroupIndex, vDataSetNames);
      for (unsigned int uiDataSetIndex = 0; (uiDataSetIndex < vDataSetNames.size()); uiDataSetIndex++) {
        affymetrix_calvin_io::DataSet* pDataSet = genericData.DataSet(uiDataGroupIndex, uiDataSetIndex);
        AffxString strDataSetName = StringUtils::ConvertWCSToMBS(vDataSetNames[uiDataSetIndex]);
        if (!isGroupDatasetPresent(groupDatasets, strDataGroupName, strDataSetName)) {
          delete pDataSet;
          continue;
        }
        //Verbose::out(1, "  DataSetName = " + strDataSetName);
        loadDataSet(strDataGroupName, strDataSetName, uiDataSetIndex, pDataSet);
      }
      if (false) { // TB replaced by calls to isGroupPresent() and isGroupDatasetPresent() as needed
        setCnCalls();
        setLohCalls();
      }
      //Verbose::out(1, "*");
    }
    //Verbose::out(1, "*");
    genericData.Clear();
    bSuccessful = true;
  } catch (...) {
    genericData.Clear(); error("Cannot read CNCychp file: " + m_strFileName);
  }
  return bSuccessful;
}

/**
 * Helper function: is the given group name present in the given group/dataset map
 * @param groupDatasets - map containing calvin group nad dataset names
 * @param groupName - calvin group name to check for in groupDatasets
 * @return - bool value. (true if groupDatasets contains groupName)
 */
bool CNCychp::isGroupPresent(const groupDatasets_t& groupDatasets, const AffxString& groupName)
{
  return groupDatasets.find(groupName) != groupDatasets.end();
}

bool CNCychp::isGroupRedundant(const AffxString& groupName)
{
  // For future arrays elaborate as needed

  // Do not read AlgorithmData if array type is CytoScanHD
  return groupName == "AlgorithmData" && getCychpHeader().getArrayType() == "CytoScanHD_Array";
}

/**
 * Helper function: is the given group name AND dataset present in the given group/dataset map
 * @param groupDatasets - map containing calvin group nad dataset names
 * @param groupName - calvin group name to check for in groupDatasets
 * @param datasetName - calvin dataset name to check for in groupDatasets for group groupName
 * @return - bool value. (true if groupDatasets contains groupName/datasetName)
 */
bool CNCychp::isGroupDatasetPresent(const groupDatasets_t& groupDatasets, const AffxString& groupName, const AffxString& datasetName)
{
  bool status = false;
  groupDatasets_t::const_iterator pos = groupDatasets.find(groupName);
  if (pos != groupDatasets.end()) {
    if (pos->second.find(datasetName) != pos->second.end()) {
      status = true;
    }
  }
  return status;
}


/**
 * Load the CNCychp header into memory.
 * @param pHeader - A pointer to the affymetrix_calvin_io::GenericDataHeader for the specified CNCychp file.
 */
void CNCychp::loadHeader(affymetrix_calvin_io::GenericDataHeader* pHeader)
{
  CNCychpHeader& header = getCychpHeader();
  affymetrix_calvin_utilities::AffymetrixGuidType guid = pHeader->GetFileId();
  header.setChpID(guid);
  int iParamCount = pHeader->GetNameValParamCnt();
  affymetrix_calvin_parameter::ParameterNameValueType param;
  bool bSelectedQCFound = m_selectedQC.empty();
  std::string paramName = "affymetrix-chipsummary-" + m_selectedQC;
  for (int iParamIndex = 0; (iParamIndex < iParamCount); iParamIndex++) {
    param = pHeader->GetNameValParam(iParamIndex);
    ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
    if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-option-xChromosome") {
      header.setXChromosome((unsigned char)param.GetValueInt32());
    }
    ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
    else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-option-yChromosome") {
      header.setYChromosome((unsigned char)param.GetValueInt32());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-chipsummary-Y-gender-call") {
      header.setGenderFromString(param.GetValueAscii());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-state-yTarget") {
      header.setYTarget(param.GetValueFloat());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-chipsummary-MAPD") {
      header.setMAPD(param.GetValueFloat());
    } else if ((m_selectedQC.empty() || (m_selectedQC == "snp-qc")) && StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-chipsummary-snp-qc") {
      header.setPVQC(param.GetValueFloat());
      if (!m_selectedQC.empty()) {header.setCQC(param.GetValueFloat());}    // special case to make sure "snp-qc" works for CQC.
      bSelectedQCFound = true;
    //} else if (m_selectedQC.empty() && StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-chipsummary-contrast-qc") {
    //  header.setCQC(param.GetValueFloat());
    } else if (m_selectedQC.empty() && StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-chipsummary-contrast-qc-nsp") {
      header.setCQC(param.GetValueFloat());
    } else if (!m_selectedQC.empty() && (StringUtils::ConvertWCSToMBS(param.GetName()) == paramName)) {
      header.setCQC(param.GetValueFloat());
      bSelectedQCFound = true;
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-chipsummary-XX") {
      header.setXX(param.GetValueInt32());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-chipsummary-Y") {
      header.setY(param.GetValueInt32());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-chipsummary-waviness-seg-count") {
      header.setWavinessSegCount(param.GetValueInt32());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-ARR-file") {
      header.setARRFileName(StringUtils::ConvertWCSToMBS(param.GetValueText()));
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-ARR-guid") {
      header.setARRGuid(StringUtils::ConvertWCSToMBS(param.GetValueText()));
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-option-reference-input") {
      header.setCNReferenceFileName(StringUtils::ConvertWCSToMBS(param.GetValueText()));
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-option-annotation-file") {
      header.setAnnotationFileName(StringUtils::ConvertWCSToMBS(param.GetValueText()));
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-array-type") {
      header.setArrayType(StringUtils::ConvertWCSToMBS(param.GetValueText()));
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-dbsnp-date") {
      header.setDbsnpDate(param.GetValueAscii());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-dbsnp-version") {
      header.setDbsnpVersion(param.GetValueAscii());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-netaffx-annotation-date") {
      header.setNetaffxAnnotDate(param.GetValueAscii());
    } else if (StringUtils::ConvertWCSToMBS(param.GetName()) == "affymetrix-algorithm-param-netaffx-build") {
      header.setNetaffxBuild(param.GetValueAscii());
    }
  }
  if ((header.getXChromosome() == 0) || (header.getYChromosome() == 0)) {
    error("Cannot find affymetrix-algorithm-param-xChromosome and/or affymetrix-algorithm-param-yChromosome is CNCychp file header. FileName: " + m_strFileName);
  }
  if (m_selectedQC != "" && bSelectedQCFound == false)
  {
      error("Cannot find the parameter in the cychp file header associated with the selected QC. Parameter: " 
          + paramName  + ", selected-qc: " + m_selectedQC + ", file name: " + m_strFileName);
  }
  pHeader->Clear();
}

/**
 * Load the CNCychp DataSet into memory.
 * @param strDataGroupName - The name of the CNCychp file data group.
 * @param strDataSetName - The name of the CNCychp file data set.
 * @param pDataSet - A pointer to the affymetrix_calvin_io::DataSet object.
 */
void CNCychp::loadDataSet(const AffxString& strDataGroupName, const AffxString& strDataSetName, unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet)
{
  try {
    pDataSet->Open();
    if (strDataGroupName == "Chromosomes") {
      if (strDataSetName == "Summary") {
        getCychpChromosomesSummaries().reserve(pDataSet->Rows());
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpChromosomesSummaries().newCychpChromosomesSummary(pDataSet, uiRowIndex);
        }
      } else {
        pDataSet->Delete(); error("The data set name in the CNCychp file is unexpected. DataGroup: " + strDataGroupName + ", DataSet: " + strDataSetName + ", File: " + m_strFileName);
      }
    } else if (strDataGroupName == "ProbeSets") {
      if (strDataSetName == "CopyNumber") {
        getCychpProbeSetsCopyNumbers().reserve(pDataSet->Rows());
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpProbeSetsCopyNumbers().newCychpProbeSetsCopyNumber(pDataSet, uiRowIndex);
        }
      } else if (strDataSetName == "AllelePeaks") { // AllelePeaks is currrently only used for browser display.
//                getCychpProbeSetsAllelePeaks().reserve(pDataSet->Rows());
//                for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++)
//                {
//                    getCychpProbeSetsAllelePeaks().newCychpProbeSetsAllelePeaks(pDataSet, uiRowIndex);
//                }
      } else {
        pDataSet->Delete(); error("The data set name in the CNCychp file is unexpected. DataGroup: " + strDataGroupName + ", DataSet: " + strDataSetName + ", File: " + m_strFileName);
      }
    } else if (strDataGroupName == "AlgorithmData") {
      if (strDataSetName == "MarkerABSignal") {
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpProbeSetsCopyNumbers().setCychpAlgorithmDataMarkerABSignal(pDataSet, uiRowIndex);
        }
      } else {
        pDataSet->Delete(); error("The data set name in the CNCychp file is unexpected. DataGroup: " + strDataGroupName + ", DataSet: " + strDataSetName + ", File: " + m_strFileName);
      }
    } else if (strDataGroupName == "Segments") {
      if (strDataSetName == "CN") {
        getCychpSegments().reserve(pDataSet->Rows());
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpSegments().newCychpSegmentCN(uiDataSetIndex, pDataSet, uiRowIndex);
        }
      } else if (strDataSetName == "LOH") {
        getCychpSegments().reserve(pDataSet->Rows());
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpSegments().newCychpSegment(uiDataSetIndex, pDataSet, uiRowIndex);
        }
      } else if (strDataSetName == "CNNeutralLOH") {
        getCychpSegments().reserve(pDataSet->Rows());
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpSegments().newCychpSegment(uiDataSetIndex, pDataSet, uiRowIndex);
        }
      } else if (strDataSetName == "NormalDiploid") {
        getCychpSegments().reserve(pDataSet->Rows());
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpSegments().newCychpSegment(uiDataSetIndex, pDataSet, uiRowIndex);
        }
      } else if (strDataSetName == "Mosaicism") {
        getCychpSegments().reserve(pDataSet->Rows());
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpSegments().newCychpSegmentMosaicism(uiDataSetIndex, pDataSet, uiRowIndex);
        }
      } else {
        pDataSet->Delete(); error("The data set name in the CNCychp file is unexpected. DataGroup: " + strDataGroupName + ", DataSet: " + strDataSetName + ", File: " + m_strFileName);
      }
    } else if (strDataGroupName == "Genotyping") {
      if (strDataSetName == "Calls") {
        for (int uiRowIndex = 0; (uiRowIndex < pDataSet->Rows()); uiRowIndex++) {
          getCychpProbeSetsCopyNumbers().setCychpGenotypingCalls(pDataSet, uiRowIndex);
        }
      } else {
        pDataSet->Delete(); error("The data set name in the CNCychp file is unexpected. DataGroup: " + strDataGroupName + ", DataSet: " + strDataSetName + ", File: " + m_strFileName);
      }
    } else {
      pDataSet->Delete(); error("The data group name in the CNCychp file is unexpected. DataGroup: " + strDataGroupName + ", File: " + m_strFileName);
    }
    pDataSet->Delete();
  } catch (...) {
    pDataSet->Delete();
    throw;
  }
}

/**
 * Load the CNCychp DataSet into memory.
 * @param strDataGroupName - The name of the CNCychp file data group.
 * @param strDataSetName - The name of the CNCychp file data set.
 * @param pDataSet - A pointer to the affymetrix_calvin_io::DataSet object.
 */
void CNCychp::loadDataSet(const AffxString& strDataGroupName, const AffxString& strDataSetName, unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet, AffxArray<CNCychpProbeSetsCopyNumber>& vProbeSets, bool bLoadProbeSetName, bool bNonNaNOnly)
{
  try {
    pDataSet->Open();
    if (strDataGroupName == "ProbeSets") {
      if (strDataSetName == "CopyNumber") {
        u_int8_t uc = 0;
        u_int32_t ui = 0;
        float f = 0;
        AffxString str;
        getCychpProbeSetsCopyNumbers().reserve(pDataSet->Rows());
		int uiRowIndex = 0;
        for (int iIndex = 0; (iIndex < pDataSet->Rows()); iIndex++) {
          if (uiRowIndex >= vProbeSets.getCount()) {
            //Err::errAbort("Not enough Probe Sets allocated to hold the entire CYCHP file.");
			break;
          }
          if (bLoadProbeSetName)
          {
              pDataSet->GetData((int)iIndex, 0, str); vProbeSets[uiRowIndex]->setProbeSetName(str);
              pDataSet->GetData((int)iIndex, 1, uc); vProbeSets[uiRowIndex]->setChromosome(uc);
              pDataSet->GetData((int)iIndex, 2, ui); vProbeSets[uiRowIndex]->setPosition(ui);
              pDataSet->GetData((int)iIndex, 3, f); vProbeSets[uiRowIndex]->setLog2Ratio(f);
			  if ((bNonNaNOnly) && (static_cast<double>(f) != static_cast<double>(f))) {continue;}
              pDataSet->GetData((int)iIndex, 5, f); vProbeSets[uiRowIndex]->setSmoothSignal(f);

          }
          else
          {
              pDataSet->GetData((int)iIndex, 1, uc); vProbeSets[uiRowIndex]->setChromosome(uc);
              pDataSet->GetData((int)iIndex, 3, f); vProbeSets[uiRowIndex]->setLog2Ratio(f);
			  if ((bNonNaNOnly) && (static_cast<double>(f) != static_cast<double>(f))) {continue;}
          }
		  uiRowIndex++;
        }
      }
    }
    pDataSet->Delete();
  } catch (...) {
    pDataSet->Delete();
    throw;
  }
}

/**
 * Determine the CN call values from the CycchpSegments data and set the CNCychpProbeSetsCopyNumber accordingly.
 */
void CNCychp::setCnCalls()
{
  CNCychpProbeSetsCopyNumbers& vCychpProbeSetsCopyNumbers = getCychpProbeSetsCopyNumbers();
  CNCychpSegments& vCychpSegments = getCychpSegments();
  for (int iSegmentIndex = 0; (iSegmentIndex < vCychpSegments.getCount()); iSegmentIndex++) {
    CNCychpSegment* pSegment = vCychpSegments.getAt(iSegmentIndex);
    if (pSegment->getDataSetIndex() == 0) { // CN
      for (int iMarkerIndex = 0; (iMarkerIndex < vCychpProbeSetsCopyNumbers.getCount()); iMarkerIndex++) {
        CNCychpProbeSetsCopyNumber* pMarker = vCychpProbeSetsCopyNumbers.getAt(iMarkerIndex);
        if (pMarker->getChromosome() == pSegment->getChromosome()) {
          if (pMarker->getPosition() >= pSegment->getStartPosition()) {
            if (pMarker->getPosition() <= pSegment->getStopPosition()) {
              pMarker->setCnCall((unsigned char)pSegment->getFState());
            }
          }
        }
      }
    }
  }
}

/**
 * Determine the LOH call values from the CNCychpSegments data and set the CNCychpProbeSetsCopyNumber accordingly.
 */
void CNCychp::setLohCalls()
{
  CNCychpProbeSetsCopyNumbers& vCychpProbeSetsCopyNumbers = getCychpProbeSetsCopyNumbers();
  CNCychpSegments& vCychpSegments = getCychpSegments();
  for (int iSegmentIndex = 0; (iSegmentIndex < vCychpSegments.getCount()); iSegmentIndex++) {
    CNCychpSegment* pSegment = vCychpSegments.getAt(iSegmentIndex);
    if (pSegment->getDataSetIndex() == 1) { // LOH
      for (int iMarkerIndex = 0; (iMarkerIndex < vCychpProbeSetsCopyNumbers.getCount()); iMarkerIndex++) {
        CNCychpProbeSetsCopyNumber* pMarker = vCychpProbeSetsCopyNumbers.getAt(iMarkerIndex);
        if (pMarker->getChromosome() == pSegment->getChromosome()) {
          if (pMarker->getPosition() >= pSegment->getStartPosition()) {
            if (pMarker->getPosition() <= pSegment->getStopPosition()) {
              pMarker->setLohCall((unsigned char)pSegment->getFState());
            }
          }
        }
      }
    }
  }
}

/**
 * @brief a Xerces SAX Handler class for parsing XML
 */
class ARRSAXHandler : public XERCES_CPP_NAMESPACE::HandlerBase
{
private:
  CNCychp* m_pCNCychp;
  unsigned int m_uiOptionCount;

public :
  ARRSAXHandler(CNCychp* p) : m_pCNCychp(p), m_uiOptionCount(0) {
    if (m_pCNCychp == NULL) {
      Err::errAbort("ARRSAXHandler must be constructed with a valid CNCychp pointer.");
    }
  }
  ~ARRSAXHandler() {}

  unsigned int getOptionCount() const {
    return m_uiOptionCount;
  }

  void warning(const XERCES_CPP_NAMESPACE::SAXParseException& exc) {
    Verbose::out(1, "WARNING: " + toString(exc.getMessage()));
  }
  void error(const XERCES_CPP_NAMESPACE::SAXParseException& exc) {
    Err::errAbort(toString(exc.getMessage()));
  }
  void fatalError(const XERCES_CPP_NAMESPACE::SAXParseException& exc) {
    Err::errAbort(toString(exc.getMessage()));
  }

  void startDocument() {
    m_uiOptionCount = 0;
  }

  /**
   * @brief Convert the internal representation to a std::string
   * @param const XMLCh* const - The pointer to convert
   * @return std::string - The resulting string value
   */
  static std::string toString(const XMLCh* const in) {
    char* p = XERCES_CPP_NAMESPACE::XMLString::transcode(in);
    std::string str = p;
    XERCES_CPP_NAMESPACE::XMLString::release(&p);
    return str;
  }

  /**
   * @brief Do the XML parsing
   * @param const XMLCh* const - The pointer to the element name
   * @param XERCES_CPP_NAMESPACE::AttributeList& - The attribute list
   */
  void startElement(const XMLCh* const name, XERCES_CPP_NAMESPACE::AttributeList& attributes) {
    try {
      std::string strElementName = toString(name);
      //        Verbose::out(1, "Element=" + strElementName);
      m_uiOptionCount++;
      std::string strGuid;
      for (unsigned int iIndex = 0; (iIndex < attributes.getLength()); iIndex++) {
        std::string strAttributeName = toString(attributes.getName(iIndex));
        std::string strAttributeType = toString(attributes.getType(iIndex));
        std::string strAttributeValue = toString(attributes.getValue(iIndex));
        //            Verbose::out(1, "AttributeName=" + strAttributeName + ", AttributeType=" + strAttributeType + ", AttributeValue=" + strAttributeValue);
        if (strAttributeName == "GUID") {
          strGuid = strAttributeValue;
        }
      }
      if (strElementName == "ArraySetFile") {
        m_pCNCychp->setARRGuid(strGuid);
      }
    } catch (...) {
      Err::errAbort("ARRSAXHandler::startElement() failed.");
    }
  }
};

/**
 * @brief Load GUD from an XML file.
 * @param strFileName - The name of the XML file to load parameters from.
 */
AffxString CNCychp::getGuidFromARRFile(const std::string& strFileName)
{
  m_strARRGuid = "";
//        Verbose::out(1, "*");
//        Verbose::out(1, "Loading guid from file: " + strFileName);
  // Initialize the XML4C system
  try {
    XERCES_CPP_NAMESPACE::XMLPlatformUtils::Initialize();
  }

  catch (const XERCES_CPP_NAMESPACE::XMLException& toCatch) {
    Err::errAbort("CNCychp::getGuidFromXMLFile() failed at XMLPlatformUtils::Initialize(). Msg: " + ARRSAXHandler::toString(toCatch.getMessage()) + " FileName: " + strFileName);
  }

  //
  //  Create a SAX parser object to use and create our SAX event handlers
  //  and plug them in.
  //
  XERCES_CPP_NAMESPACE::SAXParser* parser = new XERCES_CPP_NAMESPACE::SAXParser;
  ARRSAXHandler handler(this);
  parser->setDocumentHandler(&handler);
  parser->setErrorHandler(&handler);
  parser->setValidationScheme(XERCES_CPP_NAMESPACE::SAXParser::Val_Auto);
  parser->setDoNamespaces(false);
  parser->setDoSchema(false);
  parser->setValidationSchemaFullChecking(false);

  //
  //  Ok, lets do the progressive parse loop. On each time around the
  //  loop, we look and see if the handler has found what its looking
  //  for. When it does, we fall out then.
  //
  unsigned long duration;
  int errorCount = 0;
  try {
    // Create a progressive scan token
    XERCES_CPP_NAMESPACE::XMLPScanToken token;
    const unsigned long startMillis = XERCES_CPP_NAMESPACE::XMLPlatformUtils::getCurrentMillis();
    try {
      if (!parser->parseFirst(strFileName.c_str(), token)) {
        Err::errAbort("CNCychp::getGuidFromXMLFile() failed at parser->parseFirst(). FileName: " + strFileName);
      }
    } catch (...) {
      Err::errAbort("CNCychp::getGuidFromXMLFile() failed at parser->parseFirst(). FileName: " + strFileName);
    }

    //
    //  We started ok, so lets call scanNext() until we find what we want
    //  or hit the end.
    //
    bool gotMore = true;
    while (gotMore && !parser->getErrorCount()) {
      gotMore = parser->parseNext(token);
    }

    const unsigned long endMillis = XERCES_CPP_NAMESPACE::XMLPlatformUtils::getCurrentMillis();
    duration = endMillis - startMillis;

    errorCount = parser->getErrorCount();
    //
    //  Reset the parser-> In this simple progrma, since we just exit
    //  now, its not technically required. But, in programs which
    //  would remain open, you should reset after a progressive parse
    //  in case you broke out before the end of the file. This insures
    //  that all opened files, sockets, etc... are closed.
    //
    parser->parseReset(token);
  } catch (const XERCES_CPP_NAMESPACE::OutOfMemoryException&) {
    delete parser;
    XERCES_CPP_NAMESPACE::XMLPlatformUtils::Terminate();
    Err::errAbort("CNCychp::getGuidFromXMLFile() failed with an OutOfMemoryException. FileName: " + strFileName);
  } catch (const XERCES_CPP_NAMESPACE::XMLException& toCatch) {
    delete parser;
    XERCES_CPP_NAMESPACE::XMLPlatformUtils::Terminate();
    Err::errAbort("CNCychp::getGuidFromXMLFile() failed with an XMLException. Msg: " + ARRSAXHandler::toString(toCatch.getMessage()) + " FileName: " + strFileName);
  }

//    Verbose::out(1, "XMLFileName = " + strFileName + ", OptionCount = " + ToStr(handler.getOptionCount()));
//    Verbose::out(1, "*");
  delete parser;
  XERCES_CPP_NAMESPACE::XMLPlatformUtils::Terminate();
  return m_strARRGuid;
}
