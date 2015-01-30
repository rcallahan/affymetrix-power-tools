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
 * @file CNReporter.cpp
 *
 * @brief This file contains the CNReporter class members.
 */

#include "copynumber/CNExperiment.h"
#include "copynumber/CNReporter.h"
#include "copynumber/CNSegment.h"
#include "portability/affy-base-types.h"
#include "util/AffxStatistics.h"
#include "util/Fs.h"
#include "util/Guid.h"
//
#include <limits>
//

using namespace std;

CNReporterMethod::CNReporterMethod()
{
  m_pset = NULL;
}
void CNReporterMethod::setDataSetWriter(affymetrix_calvin_io::DataSetWriter& set)
{
  m_pset = &set;
}
int CNReporterMethod::getRowCount()
{
  return 0;
}

/**
 * @brief Constructor
 */
CNReporter::CNReporter()
{
  m_pEngine = NULL;
  m_pobjExperiment = NULL;
  m_pvProbeSets = NULL;

  m_pvMethods = NULL;
}

/**
 * @brief Destructor
 */
CNReporter::~CNReporter()
{
}


void CNReporter::defineOptions(BaseEngine& e)
{
}
void CNReporter::checkOptions(BaseEngine& e)
{
}


CNExperiment* CNReporter::getExperiment()
{
  return m_pobjExperiment;
}
CNProbeSetArray* CNReporter::getProbeSets()
{
  return m_pvProbeSets;
}
AffxArray<CNAnalysisMethod>* CNReporter::getMethods()
{
  return m_pvMethods;
}

/**
 * @brief Setup the reporter to be run
 * @param BaseEngine& - The engine associated with this reporter
 * @param CNExperiment& - The experiment
 * @param CNProbeSetArray& - The probe set vector associated with the experiment
 * @param AffxArray<CNAnalysisMethod>& - The analysis methods to extract data from
 */
void CNReporter::setup(    BaseEngine& engine,
                        CNExperiment& objExperiment,
                        CNProbeSetArray& vProbeSets,
                        AffxArray<CNAnalysisMethod>& vMethods)
{
  m_pEngine = &engine;
  m_pobjExperiment = &objExperiment;
  m_pvProbeSets = &vProbeSets;
  m_pvMethods = &vMethods;
}

/**
 * @brief Is the reporter setup to be rrun
 */
void CNReporter::isSetup()
{
  if ((m_pEngine == NULL) || (m_pobjExperiment == NULL) || (m_pvProbeSets == NULL)) {
    Err::errAbort("CNReport is not setup properly.");
  }
}


void CNReporter::loadParam(const AffxString& strName, PgOpt::PgOptType type, const AffxString& strValue, affymetrix_calvin_parameter::ParameterNameValueType& param)
{
  param.SetName(StringUtils::ConvertMBSToWCS(strName));
  switch (type) {
  case PgOpt::BOOL_OPT: param.SetValueInt8(((strValue == "true") ? (char)1 : (char)0)); break;
  case PgOpt::DOUBLE_OPT: param.SetValueFloat((float)::getDouble(strValue)); break;
  case PgOpt::INT_OPT: param.SetValueInt32(::getInt(strValue)); break;
  case PgOpt::STRING_OPT:
    if ((strName.indexOf("command-line") != -1) || (strName.indexOf("analysis") != -1) || (strName.indexOf("program-cvs-id") != -1) || (strName.indexOf("version-to-report") != -1) || (strName.endsWith("-dir"))) {
      param.SetValueText(StringUtils::ConvertMBSToWCS(strValue));
    } else {
      param.SetValueText(StringUtils::ConvertMBSToWCS(Fs::basename(strValue)));
    }
    break;
  default: throw(Except("Cannot find PgOpt type for: " + strName));
  }
}

void CNReporter::loadBuffer(char* pBuffer, int& iIndex, AffxString& str, int iLength)
{
  if (iLength == -1) {
    iLength = str.length();
  }
  unsigned int ui = htonl(iLength);
  memcpy((void*)(pBuffer + iIndex), (void*)&ui, sizeof(unsigned int));
  iIndex += sizeof(unsigned int);
  memcpy((void*)(pBuffer + iIndex), str.c_str(), iLength); iIndex += iLength;
}

void CNReporter::loadBuffer(char* pBuffer, int& iIndex, unsigned char uc)
{
  memcpy((void*)(pBuffer + iIndex), (void*)&uc, sizeof(unsigned char)); iIndex += sizeof(unsigned char);
}

void CNReporter::loadBuffer(char* pBuffer, int& iIndex, char c)
{
  memcpy((void*)(pBuffer + iIndex), (void*)&c, sizeof(char)); iIndex += sizeof(char);
}

void CNReporter::loadBuffer(char* pBuffer, int& iIndex, unsigned int ui)
{
  ui = htonl(ui); memcpy((void*)(pBuffer + iIndex), (void*)&ui, sizeof(unsigned int)); iIndex += sizeof(unsigned int);
}

void CNReporter::loadBuffer(char* pBuffer, int& iIndex, float f)
{
  type_punned pun;
  pun.v_float = f;
  loadBuffer(pBuffer, iIndex, (unsigned int)pun.v_int32);
}

AffxString CNReporter::prepareAscii(const AffxString& str, int iLength)
{
  AffxString strPrepared = str;
  int iActualLength = (int)strPrepared.length();
  strPrepared.padBlanks(iLength);
  strPrepared.setAt(iActualLength, 0);
  return strPrepared;
}

/**
 * @brief Return the array name
 * @return AffxString - The array name
 */
AffxString CNReporter::getArrayName()
{
  AffxString strArrayName;
  if (m_pEngine->getOpt("array-name") == "") {
    if ((m_pEngine->isOptDefined("cdf-file")) && (m_pEngine->getOpt("cdf-file") != "")) {
      strArrayName = Fs::basename(m_pEngine->getOpt("cdf-file"));
      int iFindIndex = strArrayName.indexOf(".");
      if (iFindIndex != -1) {
        strArrayName = strArrayName.substring(0, iFindIndex);
      }
    } else if ((m_pEngine->isOptDefined("spf-file")) && (m_pEngine->getOpt("spf-file") != "")) {
      strArrayName = Fs::basename(m_pEngine->getOpt("spf-file"));
      int iFindIndex = strArrayName.indexOf(".");
      if (iFindIndex != -1) {
        strArrayName = strArrayName.substring(0, iFindIndex);
      }
    } else {
      if (m_pEngine->isOptDefined("annotation-file")) {
        strArrayName =  Fs::basename(m_pEngine->getOpt("annotation-file"));
        int iFindIndex = strArrayName.indexOf(".");
        if (iFindIndex != -1) {
          strArrayName = strArrayName.substring(0, iFindIndex);
        }
      }
    }
  } else {
    strArrayName = m_pEngine->getOpt("array-name");
  }
  return strArrayName;
}

void CNReporter::wavinessSegCounts(int& segCountLoss, int& segCountGain, float& sd)
{
  int wavinessBlockSize = m_pEngine->getOptInt("waviness-block-size");
  int wavinessGenSpan   = m_pEngine->getOptInt("waviness-genomic-span");

  segCountLoss = 0;
  segCountGain = 0;

  bool segmentTypeMethodSeen = false;
  int iIndex;
  int iSegmentType = CNSegment::getSegmentType("cn-segment");
  int xChromosome = m_pEngine->getOptInt("xChromosome");
  int yChromosome = m_pEngine->getOptInt("yChromosome");
  const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");

  for (iIndex = 0; iIndex < (int)m_pvMethods->size(); iIndex++) {
    CNAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
    if (pMethod->isSegmentTypeAnalysis() && pMethod->getSegmentCount(iSegmentType) > 0) {
      segmentTypeMethodSeen = true;
      for (int iSegmentIndex = 0; iSegmentIndex < pMethod->getSegments().getCount(); iSegmentIndex++) {
        CNSegment* p = pMethod->getSegments().getAt(iSegmentIndex);

        // only autosome segments of genomic span greater than wavinessGenSpan
        // and containing more than wavinessBlockSize markers
        if (p->getChromosome() != xChromosome &&
            p->getChromosome() != yChromosome &&
            p->getEndPosition() - p->getStartPosition() > wavinessGenSpan  &&
            p->getMarkerCount() > wavinessBlockSize) {
          int state = p->getCall();
          if (state < 2)
            segCountLoss++;
          else if (state > 2)
            segCountGain++;
        }
      }
    }
  }
  if (!segmentTypeMethodSeen) {
    segCountLoss = -1;
    segCountGain = -1;
  }
  int probeSetSize = (int)getProbeSets()->size();
  AffxMultiDimensionalArray<float> log2Ratios(probeSetSize);
  AffxMultiDimensionalArray<float> log2RatioDiffs(probeSetSize);

  CNProbeSet* pobjProbeSetPrev = NULL;
  int log2RatioDiffsLen = 0;
  int iLog2RatioIndex = 0;
  for (iIndex = 0; iIndex < probeSetSize; iIndex++) {
    CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
	if (isCytoScanHD && (!pobjProbeSet->processAsCN())) {continue;}
    if (pobjProbeSet->getChromosome() == xChromosome || pobjProbeSet->getChromosome() == yChromosome)
      break;

    float log2Ratio = pobjProbeSet->getLog2Ratio();
    log2Ratios.set(iLog2RatioIndex, log2Ratio);
    if (iLog2RatioIndex % 2 == 1) {
      log2RatioDiffs.set(log2RatioDiffsLen, log2Ratio - pobjProbeSetPrev->getLog2Ratio());
      log2RatioDiffsLen++;
    }
    pobjProbeSetPrev = pobjProbeSet;
	iLog2RatioIndex++;
  }
  double dv = log2Ratios.var(iLog2RatioIndex) - log2RatioDiffs.var(log2RatioDiffsLen) / 2.0;
  sd = 0.0;
  if (dv > 0)
    sd = sqrt(dv);

  // file away for Report.txt
  CNExperiment* pObjExperiment = getExperiment();
  pObjExperiment->setWavinessSegCountLoss(segCountLoss);
  pObjExperiment->setWavinessSegCountGain(segCountGain);
  pObjExperiment->setWavinessSd(sd);
}

std::string CNReporter::wavinessAmplitudes()
{
  CNExperiment::wavAmplitudes_t wavAmpl = getExperiment()->getWavinessAmplitudes();
  std::sort(wavAmpl.begin(), wavAmpl.end(), wavSortCriterion());
  std::string retStr;
  for (int i = 0; i < wavAmpl.size(); i++) {
    retStr += ToStr(wavAmpl[i].second);
    retStr += ",";
  }
  if (!retStr.empty())
    retStr.resize(retStr.size() - 1);

  return retStr;
}
