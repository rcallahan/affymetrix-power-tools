////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
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
 * @file   QuantMethodFactory.cpp
 * @author Chuck Sugnet
 * @date   Wed Jan 11 14:32:27 2006
 * 
 * @brief  
 * 
 * 
 */

//
#include "chipstream/QuantMethodFactory.h"
//
#include "chipstream/QuantAvgDiff.h"
#include "chipstream/QuantBRLMM.h"
#include "chipstream/QuantBirdseedDev.h"
#include "chipstream/QuantBirdseedLegacy.h"
#include "chipstream/QuantBirdseedv1.h"
#include "chipstream/QuantBirdseedv2.h"
#include "chipstream/QuantLabelZ.h"
#include "chipstream/QuantLabelZMulti.h"
#include "chipstream/QuantMedian.h"
#include "chipstream/QuantSea.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/AffxString.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affx;

/** 
 * @brief Constructor. Registers the objects we know how to create and how
 * they describe themselves.
 */    
QuantMethodFactory::QuantMethodFactory(enum QuantType type) {
	m_PrecompFeatureEffects = NULL;
	m_iFeatureEffectSize = 0;
  if(type == Expression) {
    // Plier
    SelfDoc doc = QuantPlier::explainSelf();
    registerMethod(doc, &QuantPlier::newObject);

    // Sea
    doc = QuantSea::explainSelf();
    registerMethod(doc, &QuantSea::newObject);

    // IterPlier
    doc = QuantIterPlier::explainSelf();
    registerMethod(doc, &QuantIterPlier::newObject);

    // Rma
    doc = QuantRma::explainSelf();
    registerMethod(doc, &QuantRma::newObject);

    // Dabg
    doc = QuantDabg::explainSelf();
    registerMethod(doc, &QuantDabg::newObject);

    // AvgDiff
    doc = QuantAvgDiff::explainSelf();
    registerMethod(doc, &QuantAvgDiff::newObject);

    // Median
    doc = QuantMedian::explainSelf();
    registerMethod(doc, &QuantMedian::newObject);


  }
  else if(type == GenoType) {
    // BRLMM
    SelfDoc doc = QuantBRLMM::explainSelf();
    registerMethod(doc, &QuantBRLMM::newObject);

    // BRLMM-P
    doc = QuantLabelZ::explainSelf();
    registerMethod(doc, &QuantLabelZ::newObject);

    // BRLMM-P-MULTI
    doc = QuantLabelZMulti::explainSelf();
    registerMethod(doc, &QuantLabelZMulti::newObject);

    // Original Birdseed (v1)
    doc = QuantBirdseedv1::explainSelf();
    registerMethod(doc, &QuantBirdseedv1::newObject);

    // Broad Birdseed Sandbox
    doc = QuantBirdseedDev::explainSelf();
    registerMethod(doc, &QuantBirdseedDev::newObject);

    // Broad Birdseed v2
    doc = QuantBirdseedv2::explainSelf();
    registerMethod(doc, &QuantBirdseedv2::newObject);

    // Legacy alias of "birdseed" to "birdseed-v1"
    doc = QuantBirdseedLegacy::explainSelf();
    registerMethod(doc, &QuantBirdseedLegacy::newObject);

  }
  else {
    Err::errAbort("QuantMethodFactory() - Don't recognize type: " + ToStr(type));
  }
}

/** 
 * @brief Constructor. Registers the objects we know how to create and how
 * they describe themselves.
 */    
QuantMethodFactory::QuantMethodFactory() {
  m_PrecompFeatureEffects = NULL;
  m_iFeatureEffectSize = 0;
  // Plier
  SelfDoc doc = QuantPlier::explainSelf();
  registerMethod(doc, &QuantPlier::newObject);
  
  // Sea
  doc = QuantSea::explainSelf();
  registerMethod(doc, &QuantSea::newObject);
  
  // IterPlier
  doc = QuantIterPlier::explainSelf();
  registerMethod(doc, &QuantIterPlier::newObject);
  
  // Rma
  doc = QuantRma::explainSelf();
  registerMethod(doc, &QuantRma::newObject);
  
  // Dabg
  doc = QuantDabg::explainSelf();
  registerMethod(doc, &QuantDabg::newObject);
  
  // AvgDiff
  doc = QuantAvgDiff::explainSelf();
  registerMethod(doc, &QuantAvgDiff::newObject);
  
  // Median
  doc = QuantMedian::explainSelf();
  registerMethod(doc, &QuantMedian::newObject);
  
  // BRLMM
  doc = QuantBRLMM::explainSelf();
  registerMethod(doc, &QuantBRLMM::newObject);
  
  // BRLMM-P
  doc = QuantLabelZ::explainSelf();
  registerMethod(doc, &QuantLabelZ::newObject);
  
  // BRLMM-P-MULTI
  doc = QuantLabelZMulti::explainSelf();
  registerMethod(doc, &QuantLabelZMulti::newObject);
  
  // Original Birdseed (v1)
  doc = QuantBirdseedv1::explainSelf();
  registerMethod(doc, &QuantBirdseedv1::newObject);
  
  // Broad Birdseed Sandbox
  doc = QuantBirdseedDev::explainSelf();
  registerMethod(doc, &QuantBirdseedDev::newObject);
  
  // Broad Birdseed v2
  doc = QuantBirdseedv2::explainSelf();
  registerMethod(doc, &QuantBirdseedv2::newObject);
  
  // Legacy alias of "birdseed" to "birdseed-v1"
  doc = QuantBirdseedLegacy::explainSelf();
  registerMethod(doc, &QuantBirdseedLegacy::newObject);
  
}


QuantMethodFactory::~QuantMethodFactory()
{
//	Verbose::out(3, "QuantMethodfactory is being destroyed.");	
	if (m_PrecompFeatureEffects != NULL) {delete[] m_PrecompFeatureEffects; m_PrecompFeatureEffects = NULL;}
}

/** 
 * @brief Create a pointer to a new QuantMethod object as described
 * in the string specification.
 * @param spec - Specification string i.e. quant-norm.sketch=1000000
 * 
 * @return Pointer to new QuantMethod objects, must be deleted when
 * finished.
 */
QuantMethod *QuantMethodFactory::qMethodForString(const std::string &spec) {
  QuantMethod *stream = NULL;
  SelfCreate *create = NULL;
  create = SelfCreate::selfCreateFromString(spec, m_Docs, m_Creators, "QuantMethod");
  /* Check class type. */
  if(InstanceOf(create, QuantMethod)) {
    stream = static_cast<QuantMethod *>(create);
  }
  else {
    Err::errAbort("Class doesn't appear to be of type QuantMethod.");
  }
  return stream;
}


/** 
 * @brief Factory for creating QuantMethod from string representation.
 * 
 * @param spec - Name of QuantMethod.
 * @param layout - Probe and probe set info.
 * 
 * @return Quantification method requested.
 */
QuantMethod *QuantMethodFactory::quantMethodForString(string &spec, ChipLayout &layout, QuantType type) {
  QuantMethodFactory qFactory(type);
  QuantMethod *qMethod = NULL;
  qMethod = qFactory.qMethodForString(spec);
  assert(qMethod);
  /* Some post processing for specific objects. */
  if(InstanceOf(qMethod, QuantIterPlier)) {
    vector<int> iters;
    iters.push_back(22);
    iters.push_back(11);
    static_cast<QuantIterPlier *>(qMethod)->setIterations(iters);
  }
  else if(InstanceOf(qMethod, QuantDabg)) {
    if(m_GcControlProbes.empty()) {
      Err::errAbort("Must specify background probes when using dabg.");
    }
    static_cast<QuantDabg *>(qMethod)->setBgProbes(m_GcControlProbes);
    static_cast<QuantDabg *>(qMethod)->setGcProbes(layout.getGcProbeVec());
  }
  else if(InstanceOf(qMethod, QuantPlier) && m_PrecompFeatureEffects != NULL) {
    Verbose::out(1, "Setting quant plier with " + ToStr(layout.m_PlFactory.getApidMax()) + " feature effects.");
    static_cast<QuantPlier *>(qMethod)->setFeaturePriorEffects(m_PrecompFeatureEffects, m_iFeatureEffectSize);
  }
  else if(InstanceOf(qMethod, QuantRma) && m_PrecompFeatureEffects != NULL) {
    Verbose::out(1, "Setting quant rma with " + ToStr(layout.m_PlFactory.getApidMax()) + " feature effects.");
    static_cast<QuantRma *>(qMethod)->setFeaturePriorEffects(m_PrecompFeatureEffects, m_iFeatureEffectSize);
  }
  else if(!InstanceOf(qMethod, QuantMethod)) {
    Err::errAbort("QuantMethodFactory::quantMethodForString() - spec '" 
                  + spec + "' didn't result in a QuantMethod.");
  }
  return qMethod;
}

/** 
 * @brief Factory for creating QuantMethod from string representation.
 * 
 * @param spec - Name of QuantMethod.
 * @param layout - Probe and probe set info.
 * 
 * @return Quantification method requested.
 */
QuantMethod *QuantMethodFactory::quantMethodForString(const string &spec, PsBoard &board, QuantType type) {
  QuantMethodFactory qFactory(type);
  QuantMethod *qMethod = NULL;
  qMethod = qFactory.qMethodForString(spec);
  assert(qMethod);
  /* Some post processing for specific objects. */
  qMethod->setParameters(board);
  return qMethod;
}

/** 
 * @brief Factory for creating QuantMethod from string representation.
 * 
 * @param spec - Name of QuantMethod.
 * @param layout - Probe and probe set info.
 * 
 * @return Quantification method requested.
 */
QuantMethod *QuantMethodFactory::quantMethodForString(const string &spec, PsBoard &board) {
  QuantMethod *qMethod = NULL;
  qMethod = qMethodForString(spec);
  assert(qMethod);
  /* Some post processing for specific objects. */
  qMethod->setParameters(board);
  return qMethod;
}

/** 
 * @brief Factory for creating QuantMethod from string representation.
 * 
 * @param spec - Name of QuantMethod.
 * @param layout - Probe and probe set info.
 * 
 * @return Quantification method requested.
 */
QuantExprMethod *QuantMethodFactory::quantExprMethodForString(string &spec, ChipLayout &layout, QuantType type) {
  QuantMethod *qMethod = quantMethodForString(spec, layout, type);
  QuantExprMethod *qeMethod = NULL;
  assert(qMethod);
  qeMethod = dynamic_cast<QuantExprMethod *>(qMethod);
  if(qeMethod == NULL) 
    Err::errAbort("QuantMethodFactory::quantExprMethodForString() - " + qMethod->getType() + " doens't appear to be a QuantExprMethod.");
  return qeMethod;
}

/** 
 * @brief Factory for creating genotyping QuantMethod from string
 * representation.
 * 
 * @param spec - Name of QuantMethod.
 * @param layout - Probe and probe set info.
 * 
 * @return Quantification method requested.
 */
QuantGTypeMethod *QuantMethodFactory::quantGTypeMethodForString(string &spec, ChipLayout &layout, QuantType type) {
  QuantMethod *qMethod = quantMethodForString(spec, layout, type);
  QuantGTypeMethod *gtMethod = NULL;
  gtMethod = dynamic_cast<QuantGTypeMethod *>(qMethod);
  if(gtMethod == NULL) 
    Err::errAbort("QuantMethodFactory::quantGTypeMethodForString() - " + qMethod->getType() + " doens't appear to be a QuantGTypeMethod.");
  return gtMethod;
}

/** 
 * Open a file and read the precomputed probe effects to use in plier
 * from it. File must be tab delimited and contain the columns: "probeset_id",
 * "probe_id", and "feature_response"
 * 
 * @param fileName - Path to file containing probe effects.
 * @param effects - Map of probe ids to probe effects to be filled in.
 * @param opts - Options used in filling in.
 */
int QuantMethodFactory::openPrecompFeatureEffects(     const std::string& fileName, 
                                                        double* effects, 
                                                        const ChipLayout &layout) {
  TsvFile tsv;

  std::string probeset_id; 
  int probe_id=0;
  double feature_response=0.0;
  int allele_id=0;
  int context_id=0;
  int channel_id=0;
  int probe_count = 0;

  bool oldStyleFeatureEffectsFile=false;

  tsv.open(fileName);

  if(tsv.cname2cidx(0, "channel_id") == TSV_ERR_NOTFOUND)
  { 
    oldStyleFeatureEffectsFile = true;

    tsv.bind(0,"probeset_id",     &probeset_id,	          affx::TSV_BIND_OPTIONAL);
    tsv.bind(0,"probe_id",        &probe_id,              affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"feature_response",&feature_response,      affx::TSV_BIND_REQUIRED);
   } else
   {  
    tsv.bind(0,"probeset_id",     &probeset_id,	          affx::TSV_BIND_OPTIONAL);
    tsv.bind(0,"probe_id",        &probe_id,              affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"feature_response",&feature_response,      affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"allele_id",       &allele_id,             affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"context_id",      &context_id,            affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"channel_id",      &channel_id,            affx::TSV_BIND_REQUIRED);
   } 
  if (tsv.formatOk()!=TSV_OK) {
    Err::errAbort(ToStr("Feature effect file ") + fileName + 
                  ToStr(" must have both 'probe_id' and 'feature_response' columns."));
  }
  
  while (tsv.nextLevel(0)==TSV_OK) {
    int featureEffectsId=0;

    if(oldStyleFeatureEffectsFile)
    {
      featureEffectsId = layout.findProbeApidByProbeId( probe_id-1 );

    } else
    { 
      featureEffectsId = layout.m_PlFactory.findProbeApidByName(probeset_id, 
                                                                probe_id-1, 
                                                                allele_id, 
                                                                context_id, 
                                                                channel_id);
   }
    if(featureEffectsId!=-1)
    {
      if (effects[featureEffectsId] != 0.0)
      {
        Err::errAbort("Duplicate feature effects for same probe found in feature effects file.");
      } else 
      {
        effects[featureEffectsId]=feature_response;
      }
    } else
    {
      static int warnCount = 0;
      warnCount++;
      if(warnCount < 10) {
        Verbose::out(2, "Some feature effects have been found for probes in the probeset: " + probeset_id + " that do not appear in the spf/cdf file.  Please check the compatibility of the spf/cdf and feature effects files.");
      }
      if(warnCount == 10) {
          Verbose::out(3, "More feature effects found which are not in library files. Only reporting first 10.");
      }
    }    
    probe_count++;
  }
  tsv.close();
  return probe_count;
}

int QuantMethodFactory::openPrecompFeatureEffects_A5(affx::File5_Group* group5,
                                                      const std::string& fileName,
                                                      double* effects,
                                                      const ChipLayout &layout) 
{

  std::string probeset_id;
  int probe_id=0;
  double feature_response=0.0;
  int allele_id=0;
  int context_id=0;
  int channel_id=0;
  int probe_count = 0;

  affx::File5_Tsv* tsv5 = group5->openTsv(fileName,affx::FILE5_OPEN_RO);
  int oldStyleFeatureEffectsFileFlag = tsv5->getColumnIdx(0, "channel_id");
  while (tsv5->nextLevel(0)==affx::FILE5_OK) {
    int fileReadFlag1=0;
    int fileReadFlag2=0;
    if(oldStyleFeatureEffectsFileFlag == -1)
    {

      fileReadFlag1 = tsv5->get(0,"probe_id",        &probe_id);
      fileReadFlag2 = tsv5->get(0,"feature_response",&feature_response);

    } else
    {
      fileReadFlag1 = tsv5->get(0,"probe_id",        &probe_id);
      fileReadFlag2 = tsv5->get(0,"feature_response",&feature_response);
      tsv5->get(0,"probeset_id",      &probeset_id);
      tsv5->get(0,"allele_id",        &allele_id);
      tsv5->get(0,"context_id",       &context_id);
      tsv5->get(0,"channel_id",       &channel_id);
    }

    if (fileReadFlag1 == -1 || fileReadFlag2 == -1) {
      Err::errAbort(ToStr("Feature effect file ") + fileName +
                  ToStr(" must have both 'probe_id' and 'feature_response' columns."));
    }
    int featureEffectsId=0;
    if(oldStyleFeatureEffectsFileFlag == -1)
    {
      featureEffectsId = layout.findProbeApidByProbeId(TO_ZEROBASED(probe_id));
    } else
    {

      featureEffectsId = layout.m_PlFactory.findProbeApidByName(probeset_id, 
                                                                probe_id-1, 
                                                                allele_id, 
                                                                context_id, 
                                                                channel_id);
    } 
    if(featureEffectsId!=-1)
    {
      if (effects[featureEffectsId] != 0.0)
      {
        Err::errAbort("Duplicate feature effects for same probe found in feature effects file.");
      } else
      {
        effects[featureEffectsId]=feature_response;
      }
    } else
    {
      static int warnCount = 0;
      warnCount++;
      if(warnCount < 10) {
        Verbose::out(2, "Some feature effects have been found for probes in the probeset: " + probeset_id + " that do not appear in the spf/cdf file.  Please check the compatibility of the spf/cdf and feature effects files.");
      }
      if(warnCount == 10) {
          Verbose::out(3, "More feature effects found which are not in library files. Only reporting first 10.");
      }

    }
    probe_count++;
  }  // end while

  tsv5->close();
  delete tsv5;
  return probe_count;
}

  /** 
   * Read the feature responses from a file. File should have columns
   * for 'probeset_id', 'probe_id' and 'feature_response'.
   * 
   * @param fileName - text file with feature responses.
   */
const void QuantMethodFactory::readPrecompFeatureEffectsFromFile(	int iProbeCount,
                                                                        const std::string& fileName,
                                                                        const ChipLayout &layout) {
    int probe_count = 0;
    if (m_PrecompFeatureEffects != NULL)
    {
      delete[] m_PrecompFeatureEffects;
      m_PrecompFeatureEffects = NULL;
    }
    iProbeCount = layout.m_PlFactory.getApidMax();
    m_PrecompFeatureEffects = new double[iProbeCount];
    m_iFeatureEffectSize = iProbeCount;
    memset(m_PrecompFeatureEffects, 0, sizeof(double) * iProbeCount);
    probe_count = openPrecompFeatureEffects(fileName, m_PrecompFeatureEffects, layout);
    Verbose::out(1, "Loaded " + ToStr(probe_count) + " feature effects.");
  }


  void QuantMethodFactory::readPrecompFeatureEffectsFrom_A5(	int iProbeCount,
                                                                affx::File5_Group* group5,
                                                                const std::string& name,
                                                                const ChipLayout &layout) {
    int probe_count = 0;
    if (m_PrecompFeatureEffects != NULL)
    {
      delete[] m_PrecompFeatureEffects; 
      m_PrecompFeatureEffects = NULL;
    }
    iProbeCount = layout.m_PlFactory.getApidMax();
    m_PrecompFeatureEffects = new double[iProbeCount];
    m_iFeatureEffectSize = iProbeCount;
    memset(m_PrecompFeatureEffects, 0, sizeof(double) * iProbeCount);
    probe_count = openPrecompFeatureEffects_A5(group5,name,m_PrecompFeatureEffects, layout);
    Verbose::out(1, "Loaded " + ToStr(probe_count) + " feature effects. (A5)");
  }


void QuantMethodFactory::setupQuantMethodFactory(
        QuantMethodFactory &qFactory,
        affx::File5_File * inputFile,
        affx::File5_Group * inputGroup,
        affx::File5_File * outputFile,
        affx::File5_Group * outputGroup,
        uint32_t probeCount,
        const std::string &outDir,
        const std::string &useFeatEff,
        bool a5FeatureEffectsInputGlobal,
        const std::string &a5FeatureEffectsInputFile,
        const std::string &a5FeatureEffectsInputName,
        const std::string &a5FeatureEffectsInputGroup,
        const std::string &a5InputGroup,
        const std::string &a5Group,
        const std::string &setAnalysisName,
        const std::string &quantMethodSpec,
        const ChipLayout &layout
     ) {
  if(useFeatEff != "") {
    qFactory.readPrecompFeatureEffectsFromFile(probeCount, useFeatEff, layout);
  } else if (a5FeatureEffectsInputGlobal || a5FeatureEffectsInputFile != "") {
    string groupName = "/";
    string dataName = "feature-response";
    if(a5FeatureEffectsInputName != "") {
        dataName = a5FeatureEffectsInputName;
    } else if(setAnalysisName != "") {
        dataName = setAnalysisName + ".feature-response";
    } else {
        ///@todo this should really be based on the analysis name which means
        ///      that this functionality needs to be in a QuantMethod instance
        ///      and not in the factory
        dataName = "feature-response";
    }
    if(a5FeatureEffectsInputGroup != "") 
        groupName = a5FeatureEffectsInputGroup;
    else if(a5InputGroup != "") 
        groupName = a5InputGroup;
    else if(a5Group != "") 
        groupName = a5Group;

    if (a5FeatureEffectsInputGlobal) {
        Verbose::out(1,"Loading feat effects from global A5 file, group '" + groupName + "', data '" + dataName + "'");
        if(inputFile == NULL)
            Err::errAbort("--a5-feature-effects-input-global option given, but no global input file. Must specify --a5-global-file.");
        affx::File5_Group *a5group = inputFile->openGroup(groupName,affx::FILE5_OPEN_RO);
        qFactory.readPrecompFeatureEffectsFrom_A5(probeCount, a5group, dataName, layout);
        a5group->close();
        Freez(a5group);
    } else if (a5FeatureEffectsInputFile!="") {
        Verbose::out(1,"Loading feature effects from '" + a5FeatureEffectsInputFile + "' A5 file, group '" + groupName + "', data '" + dataName + "'");
        affx::File5_File  *a5file = new affx::File5_File();
        a5file->open(a5FeatureEffectsInputFile,affx::FILE5_OPEN_RO);
        affx::File5_Group *a5group = a5file->openGroup(groupName,affx::FILE5_OPEN_RO);
        qFactory.readPrecompFeatureEffectsFrom_A5(probeCount, a5group, dataName, layout);
        a5group->close();
        Freez(a5group);
        a5file->close();
        Freez(a5file);
    }
  }
}



