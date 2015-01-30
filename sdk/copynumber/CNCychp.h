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

#ifndef _CNCychp_H_
#define _CNCychp_H_
/**
 * @file CNCychp.h
 *
 * @brief This header contains the CNCychp class definition.
 */

#include "calvin_files/data/src/DataGroupHeader.h"
#include "calvin_files/data/src/DataSetHeader.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/GenericFileWriter.h"
#include "chipstream/BioTypes.h"
#include "file/CHPFileData.h"
#include "util/AffxArray.h"
#include "util/AffxString.h"
#include <map>
#include <string>
#include <set>
//

/**
 * @brief  A class for storing Cychp header data.
 *
 */
class CNCychpHeader
{
private:
  AffxString m_strARRFileName;
  AffxString m_strARRGuid;
  affymetrix_calvin_utilities::AffymetrixGuidType m_strChpID;
  unsigned char m_ucXChromosome;
  unsigned char m_ucYChromosome;
  affx::Gender m_eGender;
  float m_fYTarget;
  float m_fMAPD;
  float m_fPVQC;
  float m_fCQC;
  bool m_bXX;
  bool m_bY;
  AffxString m_strCNReferenceFileName;
  AffxString m_strAnnotationFileName;
  AffxString m_strArrayType;
  AffxString m_strDbsnpDate;
  AffxString m_strDbsnpVersion;
  AffxString m_strNetaffxAnnotDate;
  AffxString m_strsetNetaffxBuild;
  int m_iWavinessSegCount;

public:
  /**
   * Constructor.
  */
  CNCychpHeader();

  void clear();

  void setARRFileName(const AffxString& str);
  AffxString getARRFileName();
  void setARRGuid(const AffxString& str);
  AffxString getARRGuid();
  void setChpID(affymetrix_calvin_utilities::AffymetrixGuidType& str);
  affymetrix_calvin_utilities::AffymetrixGuidType getChpID();
  void setXChromosome(unsigned char uc);
  unsigned char getXChromosome();
  void setYChromosome(unsigned char uc);
  unsigned char getYChromosome();
  void setGender(affx::Gender eGender);
  affx::Gender getGender();
  void setYTarget(float f);
  float getYTarget();
  void setMAPD(float f);
  float getMAPD();
  void setPVQC(float f);
  float getPVQC();
  void setCQC(float f);
  float getCQC();
  void setXX(bool b);
  bool getXX();
  void setY(bool b);
  bool getY();
  void setWavinessSegCount(int i);
  int getWavinessSegCount();
  void setCNReferenceFileName(const AffxString& str);
  void setAnnotationFileName(const AffxString& str) { m_strAnnotationFileName = str; }
  void setArrayType(const AffxString& str) { m_strArrayType = str; }
  void setDbsnpDate(const AffxString& str) { m_strDbsnpDate = str; }
  void setDbsnpVersion(const AffxString& str) { m_strDbsnpVersion = str; }
  void setNetaffxAnnotDate(const AffxString& str) { m_strNetaffxAnnotDate = str; }
  void setNetaffxBuild(const AffxString& str) { m_strsetNetaffxBuild = str; }
  AffxString getCNReferenceFileName();
  AffxString getAnnotationFileName() { return m_strAnnotationFileName; }
  AffxString getArrayType() { return m_strArrayType; }
  AffxString getDbsnpDate() { return m_strDbsnpDate; }
  AffxString getDbsnpVersion() { return m_strDbsnpVersion; }
  AffxString getNetaffxAnnotDate() { return m_strNetaffxAnnotDate; }
  AffxString getNetaffxBuild() { return m_strsetNetaffxBuild; }

  void setGenderFromString(const AffxString& strIn);

  AffxString getGenderAsString();
};

/**
 * @brief  A class for storing Cychp chromosome summary data.
 *
 */
class CNCychpChromosomesSummary
{
private:
  unsigned char m_ucChromosome;
  AffxString m_strChromosomeDisplay;
  unsigned int m_uiStartIndex;
  unsigned int m_uiMarkerCount;
  float m_fMinSignal;
  float m_fMaxSignal;
  float m_fMedianCnState;
  float m_fHomFrequency;
  float m_fHetFrequency;
  float m_fMosaicism;
  float m_fLOH;

public:
  /**
   * Constructor.
  */
  CNCychpChromosomesSummary();

  void setChromosome(unsigned char uc);
  unsigned char getChromosome();
  void setChromosomeDisplay(const AffxString& str);
  AffxString getChromosomeDisplay();
  void setStartIndex(unsigned int ui);
  unsigned int getStartIndex();
  void setMarkerCount(unsigned int ui);
  unsigned int getMarkerCount();
  void setMinSignal(float f);
  float getMinSignal();
  void setMaxSignal(float f);
  float getMaxSignal();
  void setMedianCnState(float f);
  float getMedianCnState();
  void setHomFrequency(float f);
  float getHomFrequency();
  void setHetFrequency(float f);
  float getHetFrequency();
  void setMosaicism(float f);
  float getMosaicism();
  void setLOH(float f);
  float getLOH();

  /**
   * Compare function.
   * @param that - A reference to an instance of this class.
   * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
   * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
  */
  int compareTo(CNCychpChromosomesSummary& that, int iCompareCode);
};

/**
 * @brief  A class for storing an array of Cychp Chromosomes Summary data.
 *
 */
class CNCychpChromosomesSummaries : public AffxArray<CNCychpChromosomesSummary>
{
public:
  virtual ~CNCychpChromosomesSummaries();
  void clear();
  void newCychpChromosomesSummary(affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex);
};

/**
 * @brief  A class for storing Cychp ProbeSets CopyNumber data
 *
 */
class CNCychpProbeSetsCopyNumber
{
private:
  AffxString m_strProbeSetName;
  unsigned char m_ucChromosome;
  unsigned int m_uiPosition;
  float m_fLog2Ratio;
  float m_fSmoothSignal;
  float m_fAAlleleSignal;
  float m_fBAlleleSignal;
  float m_fSCAR;
  char m_cGenotypeCall;
  bool m_bIsGenotypeCallSet;
  float m_fGenotypeConfidence;
  unsigned char m_ucCnCall;
  unsigned char m_ucLohCall;
  char m_cSegmentInput;
  char m_cSegmentOutput;

public:
  /**
   * Constructor
  */
  CNCychpProbeSetsCopyNumber();

  void setProbeSetName(const AffxString& str) { m_strProbeSetName = str; }
  AffxString getProbeSetName() { return m_strProbeSetName; }
  void setChromosome(unsigned char uc) { m_ucChromosome = uc; }
  unsigned char getChromosome() { return m_ucChromosome; }
  void setPosition(unsigned int ui) { m_uiPosition = ui; }
  unsigned int getPosition() { return m_uiPosition; }
  void setLog2Ratio(float f) { m_fLog2Ratio = f; }
  float getLog2Ratio() { return m_fLog2Ratio; }
  void setSmoothSignal(float f) { m_fSmoothSignal = f; }
  float getSmoothSignal() { return m_fSmoothSignal; }
  void setAAlleleSignal(float f) { m_fAAlleleSignal = f; }
  float getAAlleleSignal() { return m_fAAlleleSignal; }
  void setBAlleleSignal(float f) { m_fBAlleleSignal = f; }
  float getBAlleleSignal() { return m_fBAlleleSignal; }
  void setSCAR(float f) { m_fSCAR = f; }
  float getSCAR() { return m_fSCAR; }
  void setGenotypeCall(char c) { m_cGenotypeCall = c; }
  char getGenotypeCall() { return m_cGenotypeCall; }
  void setIsGenotypeCallSet(bool b) { m_bIsGenotypeCallSet = b; }
  bool getIsGenotypeCallSet() { return m_bIsGenotypeCallSet; }
  void setGenotypeConfidence(float f) { m_fGenotypeConfidence = f; }
  float getGenotypeConfidence() { return m_fGenotypeConfidence; }
  void setCnCall(unsigned char uc) { m_ucCnCall = uc; }
  unsigned char getCnCall() { return m_ucCnCall; }
  void setLohCall(unsigned char uc) { m_ucLohCall = uc; }
  unsigned char getLohCall() { return m_ucLohCall; }
  void setSegmentInput(char c) { m_cSegmentInput = c; }
  char getSegmentInput() { return m_cSegmentInput; }
  void setSegmentOutput(char c) { m_cSegmentOutput = c; }
  char getSegmentOutput() { return m_cSegmentOutput; }

  /**
   * @brief Set the allele code ffor the genotype call from the internal representation (-1, 0, 1, 2)
   * @param char - The internal representation
   */
  void setGenotypeCallFromCode(char c);

  /**
   * Compare function.
   * @param that - A reference to an instance of this class.
   * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
   * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
  */
  int compareTo(CNCychpProbeSetsCopyNumber& that, int iCompareCode);

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
  char getCNNeutralLoh(int iSexCode, int iXChromosome, int iYChromosome);

  /**
   * @brief Return the NormalDiploid call (0 or 1)
   * @param int - The sex code for this sample.
   * @param int - The numeric value associated with the X chromosome
   * @param int - The numeric value associated with the Y chromosome
   * @return char - The call
   */
  char getNormalDiploid(int iSexCode, int iXChromosome, int iYChromosome);

  template<int k> struct ComparePred {
      bool operator()(const CNCychpProbeSetsCopyNumber* lhs, const CNCychpProbeSetsCopyNumber* rhs) const {
          Err::errAbort("CNCychpProbeSetsCopyNumber: ComparePred instantiated with an invalid compare code = " + ToStr(k));
          return false;
      }
  };
};

template<> struct CNCychpProbeSetsCopyNumber::ComparePred<0> {
    bool operator()(const CNCychpProbeSetsCopyNumber* lhs, const CNCychpProbeSetsCopyNumber* rhs) const {
        return lhs->m_strProbeSetName.compareTo(rhs->m_strProbeSetName, 0) < 0;
    }
};
template<> struct CNCychpProbeSetsCopyNumber::ComparePred<1> {
    bool operator()(const CNCychpProbeSetsCopyNumber* lhs, const CNCychpProbeSetsCopyNumber* rhs) const {
        int iCompareResult = AffxArray<unsigned char>::compare(lhs->m_ucChromosome, rhs->m_ucChromosome);
        if (iCompareResult == 0) {
            iCompareResult = AffxArray<unsigned int>::compare(lhs->m_uiPosition, rhs->m_uiPosition);
        }
        if (iCompareResult == 0) {
            iCompareResult = lhs->m_strProbeSetName.compareTo(rhs->m_strProbeSetName, 0);
        }
        return iCompareResult < 0;
    }
};

/**
 * @brief  A class for storing an array of Cychp Chromosome Summary data.
 */
class CNCychpProbeSetsCopyNumbers : public AffxArray<CNCychpProbeSetsCopyNumber>
{
public:
  virtual ~CNCychpProbeSetsCopyNumbers();
  void clear();
  void newCychpProbeSetsCopyNumber(affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex);

  void setCychpAlgorithmDataMarkerABSignal(affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex);
  void setCychpGenotypingCalls(affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex);
};

/**
 * @brief  A class for storing Cychp Segments data.
 */
class CNCychpSegment
{
private:
  unsigned int m_uiDataSetIndex;
  unsigned int m_uiSegmentID;
  unsigned char m_ucChromosome;
  unsigned int m_uiStartPosition;
  unsigned int m_uiStopPosition;
  int m_iMarkerCount;
  unsigned int m_uiMeanMarkerDistance;
  float m_fState;
  float m_fConfidence;
  float m_fMixture;

public:
  /**
   * Constructor
  */
  CNCychpSegment();

  void setDataSetIndex(unsigned int ui);
  unsigned int getDataSetIndex();
  void setSegmentID(unsigned int ui);
  unsigned int getSegmentID();
  void setChromosome(unsigned char uc);
  unsigned char getChromosome();
  void setStartPosition(unsigned int ui);
  unsigned int getStartPosition();
  void setStopPosition(unsigned int ui);
  unsigned int getStopPosition();
  void setMarkerCount(int i);
  int getMarkerCount();
  void setMeanMarkerDistance(unsigned int ui);
  unsigned int getMeanMarkerDistance();
  void setFState(float f);
  float getFState();
  void setConfidence(float f);
  float getConfidence();
  void setMixture(float f);
  float getMixture();

  /**
   * @brief Return the string associated with the segment type for this object.
   * @return AffxString - The segment type string.
   */
  AffxString getSegmentType();

  /**
   * Compare function.
   * @param that - A reference to an instance of this class.
   * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
   * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
  */
  int compareTo(CNCychpSegment& that, int iCompareCode);
};

/**
 * @brief  A class for storing an array of Cychp Segment data.
 */
class CNCychpSegments : public AffxArray<CNCychpSegment>
{
public:
  virtual ~CNCychpSegments();
  void clear();
  void newCychpSegmentCN(unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex);
  void newCychpSegment(unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex);
  void newCychpSegmentMosaicism(unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet, unsigned int uiRowIndex);
};

typedef std::map<std::string, std::set<std::string> > groupDatasets_t;

/**
 * @brief  A class for accessing and storing Cychp data.
 */
class CNCychp
{
private:
  AffxString m_strFileName;
  CNCychpHeader m_objCNCychpHeader;
  CNCychpChromosomesSummaries m_vCychpChromosomesSummaries;
  CNCychpProbeSetsCopyNumbers m_vCychpProbeSetsCopyNumbers;
  CNCychpSegments m_vCNCychpSegments;
  AffxString m_strARRGuid;
  AffxString m_selectedQC;

  AffxString m_strFamilialType;
  unsigned char m_ucFamilialCall;
  float m_fFamilialConfidence;
  unsigned char m_ucFamilialDuoCall;
  float m_fFamilialDuoConfidence;

public:
  static std::map<unsigned char, unsigned int> m_mMarkerCount;
  static std::map<unsigned char, int> m_mMIE_Trio;
  static std::map<unsigned char, int> m_mMIE_Mat;
  static std::map<unsigned char, int> m_mMIE_Pat;


public:
  CNCychp();

  virtual ~CNCychp();

  AffxString getFileName();

  void setFileName(const AffxString& str);

  AffxString getARRGuid();

  void setARRGuid(const AffxString& str);

  CNCychpHeader& getCychpHeader();

  CNCychpChromosomesSummaries& getCychpChromosomesSummaries();

  CNCychpProbeSetsCopyNumbers& getCychpProbeSetsCopyNumbers();
  CNCychpSegments& getCychpSegments();

  void clear();
  void warning(const AffxString& strMessage);
  void error(const AffxString& strMessage);

  bool readFile(const AffxString& strFileName, AffxArray<CNCychpProbeSetsCopyNumber>& vProbeSets, bool bLoadProbeSetName = true, bool bNonNanOnly = false);
  bool readFile(const AffxString& strFileName, const AffxString& strFamilialType, const AffxString& strPrompt = "");
  bool readFile(const AffxString& strFileName, const groupDatasets_t& groupDatasets, const AffxString& strFamilialType, const AffxString& strPrompt = "");

  AffxString getGuidFromARRFile(const std::string& strFileName);

  AffxString& getFamilialType();
  void setFamilialType(const AffxString& str);

  unsigned char getFamilialCall();
  void setFamilialCall(unsigned char uc);

  float getFamilialConfidence();
  void setFamilialConfidence(float f);

  unsigned char getFamilialDuoCall();
  void setFamilialDuoCall(unsigned char uc);

  float getFamilialDuoConfidence();
  void setFamilialDuoConfidence(float f);

  AffxString getSelectedQC();
  void setSelectedQC(const AffxString& str);

protected:
  void loadHeader(affymetrix_calvin_io::GenericDataHeader* pHeader);
  void loadDataSet(const AffxString& strDataGroupName, const AffxString& strDataSetName, unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet);
  void loadDataSet(const AffxString& strDataGroupName, const AffxString& strDataSetName, unsigned int uiDataSetIndex, affymetrix_calvin_io::DataSet* pDataSet, AffxArray<CNCychpProbeSetsCopyNumber>& vProbeSets, bool bLoadProbeSetName = true, bool bNonNaNOnly = false);
  void setCnCalls();
  void setLohCalls();
  bool isGroupPresent(const groupDatasets_t& groupDatasets, const AffxString& groupName);
  bool isGroupDatasetPresent(const groupDatasets_t& groupDatasets, const AffxString& groupName, const AffxString& datasetName);
  bool isGroupRedundant(const AffxString& groupName);
};

#endif


