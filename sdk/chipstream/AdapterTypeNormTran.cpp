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
 * @file   AdapterTypeNormTran.cpp
 *
 * @brief  Class for doing adapter type normalization
 */

//
#include "chipstream/AdapterTypeNormTran.h"
#include "chipstream/ChipStream.h"
#include "chipstream/MedNormTran.h"
#include "chipstream/SketchQuantNormTran.h"
#include "copynumber/Annotation.h"
#include "stats/stats.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Verbose.h"
//

using namespace std;

/**
 * Constructor
 *
 */
AdapterTypeNormTran::AdapterTypeNormTran(double dStyOnlySnpRatio, double dNspOnlySnpRatio, double dNspOnlyCnRatio,
    double dNspBothSnpRatio, double dNspBothCnRatio) :
    m_dStyOnlySnpRatio(dStyOnlySnpRatio), m_dNspOnlySnpRatio(dNspOnlySnpRatio), m_dNspOnlyCnRatio(dNspOnlyCnRatio),
    m_dNspBothSnpRatio(dNspBothSnpRatio), m_dNspBothCnRatio(dNspBothCnRatio)
{
  m_uiCelFileCount = 0;
  m_uiCelFileIndex = 0;
  m_bAnnotationFilesProcessed = false;
  m_Type = ADAPTERTYPENORMSTR;
  setupSelfDoc(*this);

  setOptValue("StyOnlySnpRatio", ::getDouble(m_dStyOnlySnpRatio, 6));
  setOptValue("NspOnlySnpRatio", ::getDouble(m_dNspOnlySnpRatio, 6));
  setOptValue("NspOnlyCnRatio", ::getDouble(m_dNspOnlyCnRatio, 6));
  setOptValue("NspBothSnpRatio", ::getDouble(m_dNspBothSnpRatio, 6));
  setOptValue("NspBothCnRatio", ::getDouble(m_dNspBothCnRatio, 6));
  
}

AdapterTypeNormTran::~AdapterTypeNormTran()
{
}

std::string AdapterTypeNormTran::getDefaultOptions()
{
  return "StyOnlySnpRatio=0.8287.NspOnlySnpRatio=0.7960.NspOnlyCnRatio=1.4218.NspBothSnpRatio=0.9954.NspBothCnRatio=1.6392";
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> AdapterTypeNormTran::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;
  SelfDoc::Opt StyOnlySnpRatio = {"StyOnlySnpRatio", SelfDoc::Opt::Double, "0.8287", "0.8287", "NA", "NA", "StyOnlySnpRatio"};
  SelfDoc::Opt NspOnlySnpRatio = {"NspOnlySnpRatio", SelfDoc::Opt::Double, "0.7960", "0.7960", "NA", "NA", "NspOnlySnpRatio"};
  SelfDoc::Opt NspOnlyCnRatio = {"NspOnlyCnRatio", SelfDoc::Opt::Double, "1.4218", "1.4218", "NA", "NA", "NspOnlyCnRatio"};
  SelfDoc::Opt NspBothSnpRatio = {"NspBothSnpRatio", SelfDoc::Opt::Double, "0.9954", "0.9954", "NA", "NA", "NspBothSnpRatio"};
  SelfDoc::Opt NspBothCnRatio = {"NspBothCnRatio", SelfDoc::Opt::Double, "1.6392", "1.6392", "NA", "NA", "NspBothCnRatio"};

  opts.push_back(StyOnlySnpRatio);
  opts.push_back(NspOnlySnpRatio);
  opts.push_back(NspOnlyCnRatio);
  opts.push_back(NspBothSnpRatio);
  opts.push_back(NspBothCnRatio);
  return opts;
}

void AdapterTypeNormTran::setup(const std::string& strAnnotationFileName, ChipLayout* pobjChipLayout, unsigned int uiCelFileCount)
{
  m_uiCelFileIndex = 0;
  m_uiCelFileCount = uiCelFileCount;
  m_strAnnotationFileName = strAnnotationFileName;
  m_pobjChipLayout = pobjChipLayout;
  if ((m_strAnnotationFileName == "") || (!Fs::isReadable(m_strAnnotationFileName))) {
    Err::errAbort("AdapterTypeNormTran::AdapterTypeNormTran() - An existing annotation-file must be specified.");
  }
  if (m_pobjChipLayout == NULL) {
    Err::errAbort("AdapterTypeNormTran::AdapterTypeNormTran() - ChipLayout pointer cannot be NULL.");
  }
  loadAnnotation();
}

/**
 * @brief Signal that no more data is coming (i.e. newChip() will
 * not be called anymore. Currently the class builds up all the
 * summaries and figures out the normalization factores when this
 * function is called.  This makes it difficult to pass through data
 * to downstream.
 */
void endDataSet();

void AdapterTypeNormTran::setupSelfDoc(SelfDoc &doc)
{
  doc.setDocName(ADAPTERTYPENORMSTR);
  doc.setDocDescription("Class for doing adapter type normalization. Adjust intensities by adapter type.");
  doc.setDocOptions(getDefaultDocOptions());
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc AdapterTypeNormTran::explainSelf()
{
  SelfDoc doc;
  setupSelfDoc(doc);
  return doc;
}
bool AdapterTypeNormTran::isPm(int iIndex)
{
  return m_mxAdapterPMTypes.get(iIndex) == 1;
}
bool AdapterTypeNormTran::isSnp(int iIndex)
{
  return m_mxAdapterProbeTypes.get(iIndex) == 1;
}
bool AdapterTypeNormTran::isCn(int iIndex)
{
  return m_mxAdapterProbeTypes.get(iIndex) == 2;
}
char AdapterTypeNormTran::getNspAdapterType(int iIndex)
{
  return m_mxAdapterNspTypes.get(iIndex);
}
char AdapterTypeNormTran::getStyAdapterType(int iIndex)
{
  return m_mxAdapterStyTypes.get(iIndex);
}

/**
 * @brief transform the intensity point supplied coming from a particular
 * probe in a particular microarray.
 *
 * @param iProbeIndex - Probe index from the cel file.
 * @param iChipIndex - chip index from sample.
 * @param intensity - Original intensity.
 */
inline float AdapterTypeNormTran::transform(int iProbeIndex, int iChipIndex, float intensity)
{
  char cNspAdapterType = getNspAdapterType(iProbeIndex);
  char cStyAdapterType = getStyAdapterType(iProbeIndex);
  if (isSnp(iProbeIndex)) {
    if ((cNspAdapterType == -1) && (cStyAdapterType >= 0)) {
      intensity = (float)(intensity * m_vGrandMedians.get(iChipIndex) / m_mxStyOnlySnpMedians.get(iChipIndex, cStyAdapterType) * m_dStyOnlySnpRatio);
    } else if ((cNspAdapterType >= 0) && (cStyAdapterType == -1)) {
      intensity = (float)(intensity * m_vGrandMedians.get(iChipIndex) / m_mxNspOnlySnpMedians.get(iChipIndex, cNspAdapterType) * m_dNspOnlySnpRatio);
    } else if ((cNspAdapterType >= 0) && (cStyAdapterType >= 0)) {
      intensity = (float)(intensity / m_mxStyBothSnpMedians.get(iChipIndex, cStyAdapterType));
      intensity = (float)(intensity * m_vGrandMedians.get(iChipIndex) / m_mxNspBothSnpMedians.get(iChipIndex, cNspAdapterType) * m_dNspBothSnpRatio);
    }
  } else if (isCn(iProbeIndex)) {
    if ((cNspAdapterType >= 0) && (cStyAdapterType == -1)) {
      intensity = (float)(intensity * m_vGrandMedians.get(iChipIndex) / m_mxNspOnlyCnMedians.get(iChipIndex, cNspAdapterType) * m_dNspOnlyCnRatio);
    } else if ((cNspAdapterType >= 0) && (cStyAdapterType >= 0)) {
      intensity = (float)(intensity / m_mxStyBothCnMedians.get(iChipIndex, cStyAdapterType));
      intensity = (float)(intensity * m_vGrandMedians.get(iChipIndex) / m_mxNspBothCnMedians.get(iChipIndex, cNspAdapterType) * m_dNspBothCnRatio);
    }
  }
  return intensity;
}

void AdapterTypeNormTran::newChip(std::vector<float> &data)
{
  initializeData(data);
}

void AdapterTypeNormTran::newDataSet(IntensityMart* iMart)
{
  std::vector<float> data;
  int dataSetCount = iMart->getCelDataSetCount();
  m_TransformedIMart = iMart;
  // if there aren't any chipstream nodes after this, then don't store
  // all of the intensities
  if (m_Streams.empty()) {
    m_TransformedIMart->setStoreAllCelIntensities(false);
  }
  for (int d = 0; d < dataSetCount; d++) {
    data = iMart->getCelData(d);
    newChip(data);
    for (uint32_t i = 0; i < data.size(); i++) {
      data[i] = transform(i, d, data[i]);
    }
    m_TransformedIMart->setProbeIntensity(d, data);
  }
  chipStreamPassNewChip(m_TransformedIMart);

}

/**
 * @brief Method for being passed a new cel file worth of data. Calculates
 * and stores the median (or average) of data supplied.
 * @param data - Vector of cel file data from same sample.
 */
void AdapterTypeNormTran::initializeData(const std::vector<float>& data)
{
  // Grand Median.
  AffxMultiDimensionalArray<float> vIntensities((int)data.size());
  int iProbeIndex = 0;
  for (unsigned int iIndex = 0; iIndex < data.size(); iIndex++) {
    if (isPm(iIndex)) {
      vIntensities.set(iProbeIndex, data[iIndex]);
      iProbeIndex++;
    }
  }
  m_vGrandMedians.set(m_uiCelFileIndex, vIntensities.median(iProbeIndex));

  for (int y = 0; (y < m_mxNspOnlySnpMedians.getYDimension()); y++) {
    m_mxNspOnlySnpMedians.set(m_uiCelFileIndex, y, 0);
    m_mxNspOnlyCnMedians.set(m_uiCelFileIndex, y, 0);
    m_mxStyOnlySnpMedians.set(m_uiCelFileIndex, y, 0);
    m_mxStyBothSnpMedians.set(m_uiCelFileIndex, y, 0);
    m_mxStyBothCnMedians.set(m_uiCelFileIndex, y, 0);
    m_mxNspBothSnpMedians.set(m_uiCelFileIndex, y, 0);
    m_mxNspBothCnMedians.set(m_uiCelFileIndex, y, 0);
  }

  for (char cAdapterType = 0; (cAdapterType < m_mxNspOnlySnpMedians.getYDimension()); cAdapterType++) {
    m_mxNspOnlySnpMedians.set(m_uiCelFileIndex, cAdapterType, getNspOnlyMedian(data, (char)1, cAdapterType, vIntensities));
    m_mxNspOnlyCnMedians.set(m_uiCelFileIndex, cAdapterType, getNspOnlyMedian(data, (char)2, cAdapterType, vIntensities));
    m_mxStyOnlySnpMedians.set(m_uiCelFileIndex, cAdapterType, getStyOnlyMedian(data, (char)1, cAdapterType, vIntensities));
    m_mxStyBothSnpMedians.set(m_uiCelFileIndex, cAdapterType, getStyBothMedian(data, (char)1, cAdapterType, vIntensities));
    m_mxStyBothCnMedians.set(m_uiCelFileIndex, cAdapterType, getStyBothMedian(data, (char)2, cAdapterType, vIntensities));
  }
// Verbose::progressStep(1);
  std::vector<float> med_norm(data.size());

  for (unsigned int iIndex = 0; iIndex < data.size(); iIndex++) {
    if ((getNspAdapterType(iIndex) >= 0) && (getStyAdapterType(iIndex) >= 0)) {
      if (isSnp(iIndex)) {
        med_norm[iIndex] = data[iIndex] / (float)m_mxStyBothSnpMedians.get(m_uiCelFileIndex, getStyAdapterType(iIndex));
      } else if (isCn(iIndex)) {
        med_norm[iIndex] = data[iIndex] / (float)m_mxStyBothCnMedians.get(m_uiCelFileIndex, getStyAdapterType(iIndex));
      }
    }
  }
// Verbose::progressStep(1);
  for (char cAdapterType = 0; (cAdapterType < m_mxNspBothSnpMedians.getYDimension()); cAdapterType++) {
    m_mxNspBothSnpMedians.set(m_uiCelFileIndex, cAdapterType, getNspBothMedian(med_norm, (char)1, cAdapterType, vIntensities));
    m_mxNspBothCnMedians.set(m_uiCelFileIndex, cAdapterType, getNspBothMedian(med_norm, (char)2, cAdapterType, vIntensities));
  }
  m_uiCelFileIndex++;
}

double AdapterTypeNormTran::getNspOnlyMedian(const std::vector<float>& data, char cProbeSetType, char cAdapterType, AffxMultiDimensionalArray<float>& vIntensities)
{
  int iProbeIndex = 0;
  for (unsigned int iIndex = 0; iIndex < data.size(); iIndex++) {
    if ((m_mxAdapterProbeTypes.get(iIndex) == cProbeSetType) && (getNspAdapterType(iIndex) == cAdapterType) && (getStyAdapterType(iIndex) == -1)) {
      vIntensities.set(iProbeIndex, data[iIndex]);
      iProbeIndex++;
    }
  }
  return vIntensities.median(iProbeIndex);
}

double AdapterTypeNormTran::getStyOnlyMedian(const std::vector<float>& data, char cProbeSetType, char cAdapterType, AffxMultiDimensionalArray<float>& vIntensities)
{
  int iProbeIndex = 0;
  for (unsigned int iIndex = 0; iIndex < data.size(); iIndex++) {
    if ((m_mxAdapterProbeTypes.get(iIndex) == cProbeSetType) && (getNspAdapterType(iIndex) == -1) && (getStyAdapterType(iIndex) == cAdapterType)) {
      vIntensities.set(iProbeIndex, data[iIndex]);
      iProbeIndex++;
    }
  }
  return vIntensities.median(iProbeIndex);
}

double AdapterTypeNormTran::getNspBothMedian(const std::vector<float>& data, char cProbeSetType, char cAdapterType, AffxMultiDimensionalArray<float>& vIntensities)
{
  int iProbeIndex = 0;
  for (unsigned int iIndex = 0; iIndex < data.size(); iIndex++) {
    if ((m_mxAdapterProbeTypes.get(iIndex) == cProbeSetType) && (getNspAdapterType(iIndex) == cAdapterType) && (getStyAdapterType(iIndex) != -1)) {
      vIntensities.set(iProbeIndex, data[iIndex]);
      iProbeIndex++;
    }
  }
  return vIntensities.median(iProbeIndex);
}

double AdapterTypeNormTran::getStyBothMedian(const std::vector<float>& data, char cProbeSetType, char cAdapterType, AffxMultiDimensionalArray<float>& vIntensities)
{
  int iProbeIndex = 0;
  for (unsigned int iIndex = 0; iIndex < data.size(); iIndex++) {
    if ((m_mxAdapterProbeTypes.get(iIndex) == cProbeSetType) && (getNspAdapterType(iIndex) != -1) && (getStyAdapterType(iIndex) == cAdapterType)) {
      vIntensities.set(iProbeIndex, data[iIndex]);
      iProbeIndex++;
    }
  }
  return vIntensities.median(iProbeIndex);
}

/**
 * @brief Signal that no more data is coming (i.e. newChip() will
 * not be called anymore. Currently the class builds up all the
 * summaries and figures out the normalization factores when this
 * function is called.  This makes it difficult to pass through data
 * to downstream.
 */
void AdapterTypeNormTran::endDataSet()
{
// Verbose::out(1, "endDataSet()");
  // Flush data in stream.

  chipStreamEndChips();

  // Reset for new reading cycle.
  m_uiCelFileIndex = 0;
}

/**
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 *
 * @param param - Map of key/value pairs to initialize the object.
 *
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate* AdapterTypeNormTran::newObject(std::map<std::string, std::string> &param)
{
  SelfDoc doc = explainSelf();
  double dStyOnlySnpRatio = 0;
  double dNspOnlySnpRatio = 0;
  double dNspOnlyCnRatio = 0;
  double dNspBothSnpRatio = 0;
  double dNspBothCnRatio = 0;
  fillInValue(dStyOnlySnpRatio, "StyOnlySnpRatio", param, doc);
  fillInValue(dNspOnlySnpRatio, "NspOnlySnpRatio", param, doc);
  fillInValue(dNspOnlyCnRatio, "NspOnlyCnRatio", param, doc);
  fillInValue(dNspBothSnpRatio, "NspBothSnpRatio", param, doc);
  fillInValue(dNspBothCnRatio, "NspBothCnRatio", param, doc);
  AdapterTypeNormTran* p = new AdapterTypeNormTran(dStyOnlySnpRatio, dNspOnlySnpRatio, dNspOnlyCnRatio, dNspBothSnpRatio, dNspBothCnRatio);
  return p;
}

void AdapterTypeNormTran::loadAnnotation()
{
  Verbose::out(1, "AdapterTypeNormTran processing annotation for " + ::getInt(m_pobjChipLayout->getProbeCount()) + " probes.");

  Verbose::progressBegin(1, "AdapterTypeNormTran::processAnnotationFiles(...) ", 2, 1, 2);
  /*
  Verbose::out(1, "StyOnlySnpRatio = " + ::getDouble(m_dStyOnlySnpRatio, 6));
  Verbose::out(1, "NspOnlySnpRatio = " + ::getDouble(m_dNspOnlySnpRatio, 6));
  Verbose::out(1, "NspOnlyCnRatio = " + ::getDouble(m_dNspOnlyCnRatio, 6));
  Verbose::out(1, "NspBothSnpRatio = " + ::getDouble(m_dNspBothSnpRatio, 6));
  Verbose::out(1, "NspBothCnRatio = " + ::getDouble(m_dNspBothCnRatio, 6));
  */
  m_mxAdapterProbeTypes.initialize(m_pobjChipLayout->getProbeCount()); // 0 == ProbeType (1 = Snp, 2 = Cn)
  m_mxAdapterPMTypes.initialize(m_pobjChipLayout->getProbeCount()); // 1 = PM Probe (0 or 1)
  m_mxAdapterNspTypes.initialize(m_pobjChipLayout->getProbeCount());
  m_mxAdapterStyTypes.initialize(m_pobjChipLayout->getProbeCount());
  for (int x = 0; (x < m_pobjChipLayout->getProbeCount()); x++) {
    m_mxAdapterNspTypes.set(x, -1);
    m_mxAdapterStyTypes.set(x, -1);
  }

  AffxString strProbeSetName;
  // Identify pm only probes.
  int iProbeSetCount = m_pobjChipLayout->getProbeSetCount();
  for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++) {
    ProbeListPacked List = m_pobjChipLayout->getProbeListAtIndex(iProbeSetIndex);
    ProbeSet* ps = ProbeListFactory::asProbeSet(List);
    strProbeSetName = ps->name;
    if ((strProbeSetName.startsWith("AFR")) || (strProbeSetName.startsWith("Random")) ||
		(strProbeSetName.startsWith("SNP")) || (strProbeSetName.startsWith("CN")) ||
		(strProbeSetName.startsWith("S-")) || (strProbeSetName.startsWith("C-")) ||
		(strProbeSetName.startsWith("chr"))) {
      for (unsigned int iAtomIndex = 0; (iAtomIndex < ps->atoms.size()); iAtomIndex++) {
        for (unsigned int iProbeIndex = 0; (iProbeIndex < ps->atoms[iAtomIndex]->probes.size()); iProbeIndex++) {
          if (Probe::isPm(*ps->atoms[iAtomIndex]->probes[iProbeIndex])) {
            int iIndex = ps->atoms[iAtomIndex]->probes[iProbeIndex]->id;
            m_mxAdapterPMTypes.set(iIndex, (char)1);
          }
        }
      }
    }
    delete ps;
  }
  Verbose::progressStep(1);
  int iMaxAdapterCode = 0;
  iMaxAdapterCode = Annotation().loadAnnotationForAdapterTypeNormTran(m_strAnnotationFileName, *m_pobjChipLayout, m_mxAdapterProbeTypes, m_mxAdapterNspTypes, m_mxAdapterStyTypes);

  m_vGrandMedians.initialize(m_uiCelFileCount);
  m_mxNspOnlySnpMedians.initialize(m_uiCelFileCount, iMaxAdapterCode + 1);
  m_mxNspOnlyCnMedians.initialize(m_uiCelFileCount, iMaxAdapterCode + 1);
  m_mxStyOnlySnpMedians.initialize(m_uiCelFileCount, iMaxAdapterCode + 1);
  m_mxStyOnlyCnMedians.initialize(m_uiCelFileCount, iMaxAdapterCode + 1);
  m_mxNspBothSnpMedians.initialize(m_uiCelFileCount, iMaxAdapterCode + 1);
  m_mxNspBothCnMedians.initialize(m_uiCelFileCount, iMaxAdapterCode + 1);
  m_mxStyBothSnpMedians.initialize(m_uiCelFileCount, iMaxAdapterCode + 1);
  m_mxStyBothCnMedians.initialize(m_uiCelFileCount, iMaxAdapterCode + 1);
  Verbose::progressEnd(1, "Done");
}

void AdapterTypeNormTran::setCelFileIndex(unsigned int ui) {
  m_uiCelFileIndex = ui;
}
