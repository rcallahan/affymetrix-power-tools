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

#include "copynumber/CNSegment.h"
//
#include "file/TsvFile/TsvFile.h"
//
#include "assert.h"

MosaicClass_t MosaicClassNeg(MosaicClass_t val) {
  assert(val<MCLASS__LAST);

  switch (val) {
  case MCLASS_n10:
    return MCLASS_p10;
    break;
  case MCLASS_n07:
    return MCLASS_p07;
    break;
  case MCLASS_n05:
    return MCLASS_p05;
    break;
  case MCLASS_n03:
    return MCLASS_p03;
    break;
  case MCLASS_00:
    return MCLASS_00;
    break;
  case MCLASS_p03:
    return MCLASS_n03;
    break;
  case MCLASS_p05:
    return MCLASS_n05;
    break;
  case MCLASS_p07:
    return MCLASS_n07;
    break;
  case MCLASS_p10:  // +1.0
    return MCLASS_n10;
    break;
  default:
    assert(0);
    return MCLASS_NONE;
    break;
  }
  assert(0);
  return MCLASS_NONE;
}

float MosaicClassToFloat_array[]=
  {0.0, //
   -1.0, -0.7, -0.5, -0.3,
   0.0,
   +0.3, +0.5, +0.7, +1.0};
                                  
float MosaicClassToFloat(MosaicClass_t val) {
#ifdef MOSAICISM_DEBUG
  if (!(val<MCLASS__LAST)) {
    printf("### MosaicClassToFloat: fail! val=%d\n",val);
    return -99.0;
  }
#endif
  assert(val<MCLASS__LAST);
  return MosaicClassToFloat_array[val];
};
  
//////////

CNSegment::CNSegment()
{
  m_nSegmentType = 0;
  m_cChromosome = 0;
  m_iStartPosition = 0;
  m_iEndPosition = 0;
  m_cCall = 0;
  m_fConfidence = 0;
  m_iMarkerCount = 0;
  //m_Mixture = 0.0;
  m_Mixture = MCLASS_NONE;
  m_fCalibratedCN = 0;
  m_iFamilialSampleKey = 0;
  m_fHeterozygosity = 0;
  m_fHomozygosity = 0;
  m_fMeanMarkerDistance = 0;
  m_iStartIndex = 0;
  m_iEndIndex = 0;
  m_bPseudoAutosomalRegion = false;
}

AffxString CNSegment::getSegmentName()
{
  return m_strSegmentName;
}
void CNSegment::setSegmentName(const AffxString& str)
{
  m_strSegmentName = str;
}
short CNSegment::getSegmentType()
{
  return m_nSegmentType;
}
void CNSegment::setSegmentType(short n)
{
  m_nSegmentType = n;
}
unsigned char CNSegment::getChromosome()
{
  return m_cChromosome;
}
void CNSegment::setChromosome(unsigned char c)
{
  m_cChromosome = c;
}
int CNSegment::getStartPosition()
{
  return m_iStartPosition;
}
void CNSegment::setStartPosition(int i)
{
  m_iStartPosition = i;
}
int CNSegment::getEndPosition()
{
  return m_iEndPosition;
}
void CNSegment::setEndPosition(int i)
{
  m_iEndPosition = i;
}
int CNSegment::getStartIndex()
{
  return m_iStartIndex;
}
void CNSegment::setStartIndex(int i)
{
  m_iStartIndex = i;
}
int CNSegment::getEndIndex()
{
  return m_iEndIndex;
}
void CNSegment::setEndIndex(int i)
{
  m_iEndIndex = i;
}
int CNSegment::getIndexLen() {
  return m_iEndIndex-m_iStartIndex;
}
char CNSegment::getCall()
{
  return m_cCall;
}
void CNSegment::setCall(char c)
{
  m_cCall = c;
}
float CNSegment::getConfidence()
{
  return m_fConfidence;
}
void CNSegment::setConfidence(float f)
{
  m_fConfidence = f;
}
int CNSegment::getMarkerCount()
{
  return m_iMarkerCount;
}
void CNSegment::setMarkerCount(int i)
{
  m_iMarkerCount = i;
}
int CNSegment::getFamilialSampleKey()
{
  return m_iFamilialSampleKey;
}
void CNSegment::setFamilialSampleKey(int i)
{
  m_iFamilialSampleKey = i;
}
bool CNSegment::isPseudoAutosomalRegion()
{
  return m_bPseudoAutosomalRegion;
}
void CNSegment::setPseudoAutosomalRegion(bool b)
{
  m_bPseudoAutosomalRegion = b;
}
//double CNSegment::getMixture()
//float CNSegment::getMixture()
MosaicClass_t CNSegment::getMixture()
{
  return m_Mixture;
}
double CNSegment::getMixtureAsDouble() {
  return MosaicClassToFloat(m_Mixture);
}
//void CNSegment::setMixture(double f)
//void CNSegment::setMixture(float f)
void CNSegment::setMixture(MosaicClass_t val)
{
  assert(val<MCLASS__LAST);
  m_Mixture = val;
}

float CNSegment::getCalibratedCN()
{
  return m_fCalibratedCN;
}
void CNSegment::setCalibratedCN(float f)
{
  m_fCalibratedCN = f;
}

float CNSegment::getHeterozygosity()
{
  return m_fHeterozygosity;
}
void CNSegment::setHeterozygosity(float f)
{
  m_fHeterozygosity = f;
}

float CNSegment::getHomozygosity()
{
  return m_fHomozygosity;
}
void CNSegment::setHomozygosity(float f)
{
  m_fHomozygosity = f;
}

float CNSegment::getMeanMarkerDistance()
{
  return m_fMeanMarkerDistance;
}
void CNSegment::setMeanMarkerDistance(float f)
{
  if (f != f) {
    m_fMeanMarkerDistance = 0;
  } else {
    m_fMeanMarkerDistance = f;
  }
}

/**
 * Compare function.
 * @param that - A reference to an instance of this class.
 * @param iCompareCode - The code to switch on when doing compares. (Each code is a different compare.)
 * @return - int value. (-1 if *this < that, 0 if *this == that, 1 if *this > that)
*/
int CNSegment::compareTo(CNSegment& that, int iCompareCode)
{
  int iCompareResult = 0;
  switch (iCompareCode) {
  case 0:
    iCompareResult = m_strSegmentName.compareTo(that.m_strSegmentName, 0);
    break;
  case 1:
    iCompareResult = AffxArray<char>::compare(m_cChromosome, that.m_cChromosome);
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<int>::compare(m_iStartPosition, that.m_iStartPosition);
    }
    if (iCompareResult == 0) {
      iCompareResult = AffxArray<int>::compare(m_iEndPosition, that.m_iEndPosition);
    }
    if (iCompareResult == 0) {
      iCompareResult = m_strSegmentName.compareTo(that.m_strSegmentName, 0);
    }
    break;
  }
  return iCompareResult;
}

int CNSegment::getCyhpSegmentTypeCount()
{
  return 6;
}

/**
 * @brief Return the segment type associated with the specified name
 * @param const AffxString& - The specified name
 * @return int - The segment type
 */
int CNSegment::getSegmentType(const AffxString& strName)
{
  // Non segment analysis (Needed for analysis ordering).
  if (strName == "log2-ratio") {
    return -14;
  } else if (strName == "log2-ratio-cyto2") {
    return -13;
  } else if (strName == "gaussian-smooth") {
    return -12;
  } else if (strName == "kernel-smooth") {
    return -11;
  } else if (strName == "cn-cyto2") {
    return -10;
  } else if (strName == "cn-snp7") {
    return -9;
  } else if (strName == "cn-state") {
    return -8;
  } else if (strName == "cn-gender") {
    return -7;
  } else if (strName == "cn-cyto2-gender") {
    return -6;
  } else if (strName == "genotype") {
    return -5;
  } else if (strName == "allelic-difference-CytoScan") {
    return -4;
  } else if (strName == "allelic-difference") {
    return -3;
  } else if (strName == "lohCytoScan") {
    return -2;
  } else if (strName == "loh") {
    return -1;
  }
  //cychp
  else if (strName == "cn-segment") {
    return 1;
  } else if (strName == "loh-segment") {
    return 2;
  } else if (strName == "loh-cyto2") {
    return 3;
  } else if (strName == "cn-neutral-loh") {
    return 4;
  } else if (strName == "normal-diploid") {
    return 5;
  } else if (strName == "mosaicism") {
    return 6;
  } else if (strName == "no-call") {
    return 7;
  }
  // Non segment analysis (Needed for analysis ordering).
  else if (strName == "allele-peaks") {
    return 8;
  }
  // familial
  else if (strName == "genotype-concordance") {
    return 11;
  } else if (strName == "genotype-discordance") {
    return 12;
  } else if (strName == "cn-neutral-loh-concordance") {
    return 13;
  } else if (strName == "cn-loss-loh-concordance") {
    return 14;
  } else if (strName == "hetero-upd") {
    return 15;
  } else if (strName == "iso-upd") {
    return 16;
  } else if (strName == "hemizygous-parent-of-origin") {
    return 17;
  } else if (strName == "denovo-cn") {
    return 18;
  } else {
    throw(Except("Cannot find segment type for analysis name " + strName));
  }
  return 0;
}

/**
 * @brief Return the segment type string associated with the specified segment type
 * @param short - The specified segment type
 * @return AffxString - The segment type string
 */
AffxString CNSegment::getSegmentTypeString(short nSegmentType)
{
  switch (nSegmentType) {
    // cychp
  case 1: return "CN";
  case 2: return "LOH";
  case 3: return "LOH";
  case 4: return "CNNeutralLOH";
  case 5: return "NormalDiploid";
  case 6: return "Mosaicism";
  case 7: return "NoCall";
    // familial
  case 11: return "GenotypeConcordance";
  case 12: return "GenotypeDiscordance";
  case 13: return "CNNeutralLOHConcordance";
  case 14: return "CNLossLOHConcordance";
  case 15: return "HeteroUPD";
  case 16: return "IsoUPD";
  case 17: return "HemizygousParentOfOrigin";
  case 18: return "DenovoCopyNumber";
  }
  return "Unknown";
}

AffxString CNSegment::getSegmentTypeString()
{
  return getSegmentTypeString(m_nSegmentType);
}


CNSegmentArray::CNSegmentArray() {}
CNSegmentArray::~CNSegmentArray()
{
  deleteAll();
}
int CNSegmentArray::getSegmentCount(int iType)
{
  int iCount = 0;
  for (int iIndex = 0; (iIndex < getCount()); iIndex++) {
    CNSegment* p = getAt(iIndex);
    if (p->getSegmentType() == iType) {
      iCount++;
    }
  }
  return iCount;
}

int CNSegmentArray::writeToTsv(const std::string& tsv_path) 
{
  Verbose::out(1,std::string("CNSegmentArray::writeToTsv: ")+tsv_path);

  affx::TsvFile tsv;
  tsv.defineFile("idx\t"
                 "index-start\t"
                 "index-end\t"
                 "index-size\t"
                 "type\t"
                 //"chr\t" "base-start\t" "base-end\t"
                 "mixture\t"
                 "call\t"
                 "confidence\t"
                 "calibrated-cn");
  // simple headers.
  tsv.writeTsv_v1(tsv_path);

  //
  int i_max=getCount();
  for (int i=0;i<i_max;i++) {
    CNSegment* seg = getAt(i);
    int cidx=0;
    tsv.set(0,cidx++,i);
    tsv.set(0,cidx++,seg->getStartIndex());
    tsv.set(0,cidx++,seg->getEndIndex());
    tsv.set(0,cidx++,seg->getIndexLen());
    tsv.set(0,cidx++,seg->getSegmentType());
    //tsv.set(0,cidx++,seg->getChromosome());
    //tsv.set(0,cidx++,seg->getStartPosition());
    //tsv.set(0,cidx++,seg->getEndPosition());
    tsv.set(0,cidx++,MosaicClassToFloat(seg->getMixture()));
    tsv.set(0,cidx++,seg->getCall());
    tsv.set(0,cidx++,seg->getConfidence());
    tsv.set(0,cidx++,seg->getCalibratedCN());
    //
    tsv.writeLevel(0);
  }
  tsv.close();
  //
  return 0;
}  

CNSegmentOverlapArray::CNSegmentOverlapArray()
{
}

CNSegmentOverlapArray::~CNSegmentOverlapArray() {
  deleteAll();
}
