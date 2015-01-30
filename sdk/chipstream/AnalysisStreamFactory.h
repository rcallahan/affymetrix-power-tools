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
 * @file   AnalysisStreamFactory.h
 * @author Chuck Sugnet
 * @date   Wed Jan 11 15:50:22 2006
 *
 * @brief  Object for making new analysis streams.
 */

#ifndef _ANALYSISSTREAMFACTORY_H_
#define _ANALYSISSTREAMFACTORY_H_

//
#include "chipstream/AnalysisStream.h"
#include "chipstream/AnalysisStreamExpPcaSel.h"
#include "chipstream/AnalysisStreamExpSpectSel.h"
#include "chipstream/AnalysisStreamExpression.h"
#include "chipstream/AnalysisStreamGType.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStreamFactory.h"
#include "chipstream/PmAdjusterFactory.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantMethodFactory.h"
//
#include <cstring>
#include <map>
#include <string>
#include <vector>
//

/**
 * @brief Object for making new analysis streams from a text
 * specification. Encapsulates all the work of separately making
 * chipstream, pm adjuster and quantification methods.
 */
class AnalysisStreamFactory {

public:

  /**
   * Constructor.
   * @param What type of Quantification method will be used? I.e.
   * Expression, GenoType, ...
   */
  AnalysisStreamFactory(QuantMethodFactory::QuantType type = QuantMethodFactory::Expression) : m_QuantMethFac(type) {
    m_QuantType = type;
    m_Docs.push_back(AnalysisStreamExpression::explainSelf());
    m_Creators.push_back(&AnalysisStreamExpression::newObject);
    // pca-select probes stream.
    m_Docs.push_back(AnalysisStreamExpPcaSel::explainSelf());
    m_Creators.push_back(&AnalysisStreamExpPcaSel::newObject);
    // spect-select probes stream.
    m_Docs.push_back(AnalysisStreamExpSpectSel::explainSelf());
    m_Creators.push_back(&AnalysisStreamExpSpectSel::newObject);

	m_uiCelFileCount = 0;
  };

  /** Destructor. */
  ~AnalysisStreamFactory() {
    for(unsigned int i = 0; i < m_GcControlProbes.size(); i++) {
      delete m_GcControlProbes[i];
    }
  }

  /**
   * @brief Factory for created AnalysisStreams from a string description
   *
   * @param s - String description.
   * @param layout - Probe and probe set info.
   * @param stdMethods - Aliases for standard methods like RMA.
   *
   * @return analysis stream requested.
   */
  AnalysisStream *constructAnalysisStream(const std::string& s, ChipLayout &layout,
                                          std::map<std::string,std::string> &stdMethods,
                                          std::string analysisName = "");



 /**
   * @brief Factory for creatng mRNA expression AnalysisStream from a string description.
   *
   * @param s - String description.
   * @param layout - Probe and probe set info.
   * @param stdMethods - Aliases for standard methods like RMA.
   *
   * @return analysis stream requested.
   */
  AnalysisStreamExpression *constructExpressionAnalysisStream(const std::string& s, ChipLayout &layout,
                                                              std::map<std::string,std::string> &stdMethods,
                                                              std::string analysisName = "");

  AnalysisStreamExpression *constructExpressionAnalysisStages(const std::string& s, ChipLayout &layout, 
                                                               std::map<std::string,std::string> &stdMethods,
                                                              std::string analysisName);
  
  /**
   * @brief Factory for creatng Genotype SNP chip AnalysisStream from a string description.
   *
   * @param s - String description.
   * @param layout - Probe and probe set info.
   * @param stdMethods - Aliases for standard methods like RMA.
   *
   * @return analysis stream requested.
   */
  AnalysisStreamGType *constructGTypeAnalysisStream(const std::string& s, ChipLayout &layout,
                                                    std::map<std::string,std::string> &stdMethods,
                                                    std::string analysisName = "");


  /**
   * @brief Factory for created AnalysisStreams from a string description
   *
   * @param s - String description.
   * @param layout - Probe and probe set info.
   * @param analysisName -
   * @return analysis stream requested.
   */
  AnalysisStream *constructAnalysisStream(const std::string& s, ChipLayout &layout, std::string analysisName = "") {
    std::map<std::string,std::string> stdMethods;
    return constructAnalysisStream(s, layout, stdMethods, analysisName);
  }

  /**
   * @brief Load up a list of probes from a .bgp file
   *
   * @param fileName - name of bgp file.
   * @param controlProbes - vector to be filled in with probes.
   */
  static void probeListFromBgpFile(const std::string& fileName,
                                   std::vector<Probe *> &controlProbes);

  /**
   * Specify the probes that should be used for gc controls
   * @param gcControlProbes - Which probes are considered background.
   */
  void setGcControlProbes(const std::vector<Probe *> &gcControlProbes) {
    for(unsigned int i = 0; i < m_GcControlProbes.size(); i++) {
      delete m_GcControlProbes[i];
    }
    m_GcControlProbes.clear();
    for(unsigned int i = 0; i < gcControlProbes.size(); i++) {
      Probe *p = Probe::cloneProbe(*gcControlProbes[i]);
      m_GcControlProbes.push_back(p);
    }
    m_QuantMethFac.setGcControlProbes(m_GcControlProbes);
    m_PmAdjustFac.setControlProbes(m_GcControlProbes);
    m_CsFac.setControlProbes(m_GcControlProbes);
  }

  /**
   * Get the vector of probes that are considered to be background.
   * @return vector of the probe pointers that are being used.
   */
  const std::vector<Probe *> &getGcControlProbes() {
    return m_GcControlProbes;
  }

  /**
   * Parse the feature responses from a text file.  File must be tab
   * delimited and contain the columns: "probeset_id", "atom_id",
   * and "feature_response"
   *
   * @param fileName - file containing feature responses.
   */
  void readPrecompFeatureEffectsFromFile(int iProbeCount, const std::string& fileName, const ChipLayout layout) {
    m_QuantMethFac.readPrecompFeatureEffectsFromFile(iProbeCount, fileName, layout);
  }

  void readPrecompFeatureEffectsFrom_A5(int iProbeCount, affx::File5_Group* group5,
                                        const std::string& name, const ChipLayout layout) {
    m_QuantMethFac.readPrecompFeatureEffectsFrom_A5(iProbeCount, group5, name, layout);
  }

  /**
   * Read in the control probes from a text file. File must be tab
   * delimited and contain columns named "probe_id" and "gc_count"
   *
   * @param fileName - file containing the control probes.
   */
  void readControlProbes(const std::string& fileName) {
    probeListFromBgpFile(fileName, m_GcControlProbes);
    m_QuantMethFac.setGcControlProbes(m_GcControlProbes);
    m_PmAdjustFac.setControlProbes(m_GcControlProbes);
    m_CsFac.setControlProbes(m_GcControlProbes);
  }

  /**
   * Toggle whether or not we're writing out the normalization
   * distribution.
   * @param writeSketch - true to write out the distrubtion false
   * otherwise.
   */
  void setWriteSketchDir(const std::string &writeSketchDir) {
    m_CsFac.setWriteSketchDir(writeSketchDir);
  }
  void setWriteSketchDir_a5(const std::string &writeSketchDir) {
    m_CsFac.setWriteSketchDir_a5(writeSketchDir);
  }
  void setWriteSketchGroup_a5(affx::File5_Group* group5) {
    m_CsFac.setWriteSketchGroup_a5(group5);
  }

  /**
   * Profile
   */
  void setWriteProfileDir(const std::string &writeProfileDir){
	m_CsFac.setWriteProfileDir(writeProfileDir);
}

   void setreadReferenceProfile(const std::string &fileName){
	   m_CsFac.setreadReferenceProfile(fileName);
   }


  /**
   * Read in the target sketch from a
   *
   * @param fileName
   */
  void readTargetSketchFromFile(const std::string& fileName) {
    m_CsFac.readTargetSketchFromFile(fileName);
  }
  void readTargetSketchFromFile_a5(const std::string& fileName) {
    m_CsFac.readTargetSketchFromFile_a5(fileName);
  }
  void readTargetSketchFromFile_a5_tsv(affx::File5_Group* group5,const std::string& tsv_name) {
    m_CsFac.readTargetSketchFromFile_a5_tsv(group5,tsv_name);
  }

  /** 
   * Set the annotation-file
   * 
   * @param fileName 
   */
  void setAnnotationFileName(const std::string& str) {m_CsFac.setAnnotationFileName(str);}

  /*
   * Set teh distribution for doing quantile normalizations.
   * @param sketch - distribution to normalize to.
   */
  //void setTargetSketch(const std::vector<float> &sketch) {
  //  m_CsFac.setTargetSketch(sketch);
  //}

    /**
   * Get the vector of documentation objects that this factory is aware of.
   * @return - vector of documentation objects.
   */
  std::vector<SelfDoc> getDocs() {
    return m_Docs;
  }

  void setCelFileCount(unsigned int ui) {m_uiCelFileCount = ui; m_CsFac.setCelFileCount(ui);}

  void static setupAnalysisStreamFactory(
        AnalysisStreamFactory &asFactory,
        affx::File5_File * inputFile,
        affx::File5_Group * inputGroup,
        affx::File5_File * outputFile,
        affx::File5_Group * outputGroup,
        uint32_t celFileCount,
        uint32_t probeCount,
        const std::string &outDir,
        bool writeSketch,
        const std::string &targetSketch,
	bool writeProfile,
	const std::string &referenceProfile,
        bool a5Sketch,
        bool a5SketchUseGlobal,
        bool a5SketchInputGlobal,
        const std::string &a5SketchInputFile,
        const std::string &a5SketchInputName,
        const std::string &a5SketchInputGroup,
        const std::string &useFeatEff,
        bool a5FeatureEffectsInputGlobal,
        const std::string &a5FeatureEffectsInputFile,
        const std::string &a5FeatureEffectsInputName,
        const std::string &a5FeatureEffectsInputGroup,
        const std::string &a5InputGroup,
        const std::string &a5Group,
        const std::string &setAnalysisName,
        const std::string &quantMethodSpec,
        const std::string &bgpFile,
        const std::string &AnnotFile,
        const ChipLayout &layout
     );
  ChipStreamDataTransform *chipstreamStageFromSpec(AnalysisStream *stream,
                                                   ChipLayout &layout);

public:
  /// Module for making chipstreams.
  ChipStreamFactory m_CsFac;
  /// Module for making pmadjusters.
  PmAdjusterFactory m_PmAdjustFac;
  /// Module for making quantification methods.
  QuantMethodFactory m_QuantMethFac;

protected:

  AnalysisStreamExpression *tryToMakeAnalysisStream(const std::string &spec);

  /**
   * Utility function to fill the chipstream and pm adjuster portions of
   * analysisstream. Allows functions like constructAnalysisStream() and
   * constructExpressionAnalysisStream() to share same methods for doing this
   * and then just cap with QuantificationMethod.
   *
   * @param stream - Stream to be filled in.
   * @param layout - Information about layout of probes on chip.
   * @param pmSpec - String description of pm adjuster.
   * @param chipStreamSpec - Vector of strings describing chipstream modules to be created.
   * @param name - Name to be filled in as pm adjusters and chipstream modules are created.
   */
  void fillInAnalysisStream(AnalysisStream *stream, ChipLayout &layout,
                            std::string &pmSpec, std::vector<std::string> &chipStreamSpec,
                            std::string &name);

  void fillInAnalysisStages(AnalysisStream *stream,
                            ChipLayout &layout,
                            const std::string &pmSpec, 
                            std::vector<std::string> &chipStreamSpec, 
                            std::string &name);

  /// Probes to be used for background distributions.
  /// Note that memory is owned elsewhere.
  std::vector<Probe *> m_GcControlProbes;
  /// What type of quantification are we doing? i.e. Expression, GenoType
  QuantMethodFactory::QuantType m_QuantType;
  /// Self documentation
  std::vector<SelfDoc> m_Docs;
  /// Self creation
  std::vector<SelfCreate::selfCreator> m_Creators;

  unsigned int m_uiCelFileCount;
};

#endif /* _ANALYSISSTREAMFACTORY_H_ */
