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

#ifndef _Segment_H_
#define _Segment_H_
/**
 * @file CNSegment.h
 *
 * @brief This header contains the Segment and SegmentArray class definitions.
 */

#include "util/AffxArray.h"
#include "util/AffxString.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"

//
enum  MosaicClass_t {
  MCLASS_NONE, // unset
  MCLASS_n10,  // -1.0 // 1
  MCLASS_n07,          // 2
  MCLASS_n05,          // 3
  MCLASS_n03,          // 4
  MCLASS_00,           // 5
  MCLASS_p03,          // 6
  MCLASS_p05,          // 7
  MCLASS_p07,          // 8
  MCLASS_p10,  // +1.0
  MCLASS__LAST,
};

/// negate the class (MCLASS_pXX -> MCLASS_nXX)
MosaicClass_t MosaicClassNeg(MosaicClass_t val);
float         MosaicClassToFloat(MosaicClass_t val);

/**
 * @brief  A class for storing Segment data
 */
class CNSegment
{
private:
  AffxString m_strSegmentName;
  short m_nSegmentType;
  unsigned char m_cChromosome;
  int m_iStartPosition;
  int m_iEndPosition;
  char m_cCall;
  float m_fConfidence;
  int m_iMarkerCount;
  //double m_fMixture;
  //float m_fMixture;
  MosaicClass_t m_Mixture;
  float m_fCalibratedCN;
  int m_iFamilialSampleKey;
  float m_fHeterozygosity;
  float m_fHomozygosity;
  float m_fMeanMarkerDistance;
  int m_iStartIndex;
  int m_iEndIndex;
  bool m_bPseudoAutosomalRegion;

public:
  /**
   * @brief Constructor
   */
  CNSegment() ;

  AffxString getSegmentName();
  void setSegmentName(const AffxString& str);
  short getSegmentType();
  void setSegmentType(short n);
  unsigned char getChromosome();
  void setChromosome(unsigned char c);
  int getStartPosition();
  void setStartPosition(int i);
  int getEndPosition();
  void setEndPosition(int i);
  int getStartIndex();
  void setStartIndex(int i);
  int getEndIndex();
  void setEndIndex(int i);
  int getIndexLen(); // getEndIndex-getStartIndex
  char getCall();
  void setCall(char c);
  float getConfidence();
  void setConfidence(float f);
  int getMarkerCount();
  void setMarkerCount(int i);
  int getFamilialSampleKey();
  void setFamilialSampleKey(int i);
  bool isPseudoAutosomalRegion();
  void setPseudoAutosomalRegion(bool b);

  //double getMixture();
  //float getMixture();
  //void setMixture(double f);
  //void setMixture(float f);
  MosaicClass_t getMixture();
  double getMixtureAsDouble();
  void setMixture(MosaicClass_t val);

  float getCalibratedCN();
  void setCalibratedCN(float f);

  float getHeterozygosity();
  void setHeterozygosity(float f);

  float getHomozygosity();
  void setHomozygosity(float f);

  float getMeanMarkerDistance();
  void setMeanMarkerDistance(float f);

  /**
   * Compare function.
   * @param that - A reference to an instance of this class.
   * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
   * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
  */
  int compareTo(CNSegment& that, int iCompareCode);

  static int getCyhpSegmentTypeCount();

  /**
   * @brief Return the segment type associated with the specified name
   * @param const AffxString& - The specified name
   * @return int - The segment type
   */
  static int getSegmentType(const AffxString& strName);

  /**
   * @brief Return the segment type string associated with the specified segment type
   * @param short - The specified segment type
   * @return AffxString - The segment type string
   */
  static AffxString getSegmentTypeString(short nSegmentType);

  AffxString getSegmentTypeString();
};

/**
 * @brief  A vector of Segment Pointers.
 */
class CNSegmentArray : public AffxArray<CNSegment>
{
public:
  CNSegmentArray();
  ~CNSegmentArray() ;
  int getSegmentCount(int iType);

  int writeToTsv(const std::string& tsv_path);
};

/**
 * @brief  A class for storing SegmentOverlap data.
 */
class CNSegmentOverlap
{
public:
  AffxString SegmentType;
  int ReferenceSampleKey;
  AffxString ReferenceSegmentID;
  int FamilialSampleKey;
  AffxString FamilialSegmentID;
};

/**
 * @brief  A vector of SegmentOverlap Pointers.
 */
class CNSegmentOverlapArray : public AffxArray<CNSegmentOverlap>
{
public:
  CNSegmentOverlapArray();
  ~CNSegmentOverlapArray();
};

#endif
