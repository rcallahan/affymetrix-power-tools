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
 * @file   ChipStreamFactory.h
 * @author Chuck Sugnet
 * @date   Tue Oct 25 11:56:37 2005
 * 
 * @brief Factory class for making chip streams based on a string
 * representation.
 */

#ifndef _CHIPSTREAMFACTORY_H_
#define _CHIPSTREAMFACTORY_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStream.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
//
#include "file5/File5.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include <cstring>
#include <string>
#include <vector>
//



/**
 * Factory class for making chip streams based on a string representation. 
 */
class ChipStreamFactory {

public:

  /** 
   * @brief Constructor. Registers the objects we know how to create and how
   * they describe themselves.
   */  
  ChipStreamFactory();

  /** 
   * @brief Create a pointer to a new ChipStream object as described
   * in the string specification.
   * @param spec - Specification string i.e. quant-norm.sketch=1000000
   * 
   * @return Pointer to new ChipStream objects, must be deleted when
   * finished.
   */
  ChipStream *cStreamForString(const std::string &spec);

  /** 
   * Set the distribution for doing quantile normalizations.
   * @param sketch - distribution to normalize to.
   */
  void setTargetSketch(const std::vector<float> &sketch) {
    m_TargetSketch = sketch;
  }

  /** 
   * Get the skize of the sketch distribution for doing quantile normalizations
   */
  int getTargetSketchSize() {
      return m_TargetSketch.size() > 0 ? int(m_TargetSketch.size()) : -1;
  }

  /** 
   * Set the annotaiton-file
   * @param str - the annotation file name
   */
  void setAnnotationFileName(const std::string& str) {m_strAnnotationFileName = str;}

  /** 
   * Toggle whether or not we're writing out the normalization
   * distribution.
   * @param writeSketch - true to write out the distrubtion false
   * otherwise.
   */
  void setWriteSketchDir(const std::string &writeSketchDir) {
    m_WriteSketchDir = writeSketchDir;
  }
  void setWriteSketchDir_a5(const std::string &writeSketchDir) {
    // printf("### setWriteSketchDir_a5('%s')\n",writeSketchDir.c_str());
    m_WriteSketchDir_a5 = writeSketchDir;
  }
  void setWriteSketchGroup_a5(affx::File5_Group* group5) {
    // printf("### setWriteSketchDir_a5('%s')\n",writeSketchDir.c_str());
    m_WriteSketchGroup_a5 = group5;
  }

  // toggle whether we're writing the normalization dist
  // passed down one level from AnalysisFactory
  void setWriteProfileDir(const std::string &writeProfileDir){
	  m_WriteProfileDir = writeProfileDir;
  }

  void setreadReferenceProfile(const std::string &readProfile){
	  m_ReadProfile= readProfile;
  }

  void setProbesForNorm(const std::vector<int> &probeIds) {
    m_NormProbes = probeIds;
  }

  void readNormProbesFromFile(const std::string& fileName);


  /** 
   * @brief Factory for making ChipStreams from a string description. 
   * 
   * @param spec - Specification in form "chipstream.key=value.key=value"
   * @param layout - Probe and probe set info.
   * @param description - name of chipstream path so far.
   * @return ChipStream requested.
   */
  ChipStream *chipStreamForString(const std::string &spec, ChipLayout &layout, const std::string &description);

  /** 
   * @brief Factory for making ChipStreams from a string description. 
   * 
   * @param board - Blackboard with program state.
   * @param spec - Specification in form "chipstream.key=value.key=value"
   * @param description - name of chipstream path so far.
   * @return ChipStream requested.
   */
   ChipStream *chipStreamForStringBoard(PsBoard &board,
                                        const std::string &spec, 
                                        const std::string &description);
  
  /** 
   * Open and read a target distribution from a text file. Must have a column
   * called 'intensities' with distribution quantiles sorted from highest to 
   * lowest.
   * 
   * @param fileName - text file containing column of quantiles.
   */
  void readTargetSketchFromFile(const std::string& fileName);
  void readTargetSketchFromFile_a5(const std::string& fileName);
  void readTargetSketchFromFile_a5_tsv(affx::File5_Group* group5,
                                       const std::string& name);
  /** 
   * Get the vector of documentation objects that this factory is aware of.
   * @return - vector of documentation objects.
   */
  std::vector<SelfDoc> getDocs() {
    return m_Docs;
  }

  /** 
   * Set the probes to be used as a background distribution.
   * @param controlProbes - probes to use as controls.
   */
  void setControlProbes(const std::vector<Probe *> controlProbes) {
    m_GcControlProbes = controlProbes;
  }
  
  void setCelFileCount(unsigned int ui) {m_uiCelFileCount = ui;}

  /** 
   * Can this factory be used to build from this specification?
   * @param spec - Specification in form "chipstream.key=value.key=value"
   * @return true if this factory knows how to build from this specification, false otherwise
   */
  bool canBuild(const std::string &spec);
 
protected:

  /// Should we write out our normalization distribution.
  std::string m_WriteSketchDir;
  std::string m_WriteSketchDir_a5;
  affx::File5_Group* m_WriteSketchGroup_a5;

  /// Should we write out our normalization distribution.
  // This is passed from the AnalysisFactory level
  std::string m_WriteProfileDir;
  std::string m_WriteProfileDir_a5;
  affx::File5_Group* m_WriteProfileGroup_a5;

  // this level should know nothing about profiles, only files
  // This is passed from the AnalysisFactory level
  std::string m_ReadProfile;

  
  /// Distribution to normalize data to.  
  std::vector<float> m_TargetSketch;      
  /// Probes to use for normalization.
  std::vector<int> m_NormProbes;
  /// Self documentation
  std::vector<SelfDoc> m_Docs;
  /// Self creation
  std::vector<SelfCreate::selfCreator> m_Creators;
  /// Background probes to use. Note that memory is owned elsewhere.
  std::vector<Probe *> m_GcControlProbes;
  /// annotaiton-file
  std::string m_strAnnotationFileName;
  /// intensity-reporter.file
  std::string m_strIntensityReporterFileName;

  unsigned int m_uiCelFileCount;
};

#endif /* _CHIPSTREAMFACTORY_H_ */
