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
 * @file   ArtifactReduction.h
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:17:07 2005
 * 
 * @brief Class for doing normalization. Can do sketch and full quantile (just
 * set sketch to chip size) and supports bioconductor compatibility.
 */
#ifndef _ARTIFACTREDUCTION_H_
#define _ARTIFACTREDUCTION_H_

//
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStream.h"
#include "chipstream/IntensityTransformer.h"
#include "chipstream/ProbeList.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/ChipLayout.h"
//
#include "file5/File5.h"
#include "portability/affy-base-types.h"
//
#include <map>
#include <string>
#include <vector>

/// String describing quantile norm algorithm
#define ARTIFACTREDUCTIONSTR "artifact-reduction"
#define ARTIFACTREDUCTIONTRUSTFILE "artifact-reduction-trust"

/**
 * ArtifactReduction for doing normalization. Can do sketch and full
 * quantile (just set sketch to chip size) and supports bioconductor
 * compatibility. Unfortunatley this is kind of a bad stream currently as it
 * doesn't really flow any data through...
 */
class ArtifactReduction : public ChipStream {

public:

  /** 
   * Constructor.
   */
  ArtifactReduction();

  /**
   * Destructor. 
   */
  ~ArtifactReduction();
  
  /** 
   * @brief Add a stream to the list of those that will receive
   * downstream data.
   * @param stream - ChipStream that wants to be fed our modified
   * data.
   */
  virtual void registerStream(ChipStream *stream) {
    m_Streams.push_back(stream);
  }
  
  /** 
   * @brief register our parent stream that will be passing
   * data to this object.
   * @param stream - who are we getting data from?
   */
  virtual void registerParent(ChipStream *stream) {
    m_ParentStream = stream;
  }

  /** 
   * @brief Get a reference to our parent stream.
   * @return Reference to parent stream.
   */
  virtual const ChipStream &getParent() {
    return *m_ParentStream;
  } 

  /**
   * @brief Placeholder function for transform()
   *
   * @param probeIx - Probe index on chip.
   * @param chipIx - Set of chip indexes from same sample.
   * @param intensity - CEL intensity
   * @param board - Blackboard with various state data
   */
  virtual float transform(int probeIx, int chipIx, float intensity, PsBoard &board);

  /**
   * @brief Placeholder function for transform()
   *
   * @param probeIx - Probe index on chip.
   * @param chipIx - Set of chip indexes from same sample.
   * @param intensity - CEL intensity
   */
  virtual float transform(int probeIx, int chipIx, float intensity);

  void doReduction(IntensityMart* iMart);

  /** 
   * @brief Method for being passed a new cel file worth of data.
   * @param data - vector of cel file data.
   */
  void newDataSet(IntensityMart* iMart);
   
  virtual void newChip(std::vector<float> &data);

  virtual void finishedChips();

  /** 
   * @brief Signal that no more data is coming (i.e. newChip() will not be
   * called anymore. Currently the class builds up all the sketches and then
   * does the normalization when this function is called. This makes it
   * difficult to pass through data to downstream. If a precomputed sketch
   * is supplied then the data can be passed through chip by chip.
   */
  void endDataSet();

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions();

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setUpSelfDoc(SelfDoc &doc);

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf() {
    SelfDoc doc;
    setUpSelfDoc(doc);
    return doc;
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
  static SelfCreate *newObject(std::map<std::string,std::string> &param);

 

  /** 
   * @brief Set a filename to write out the quantile distrubtion to. This
   * can then be used to normalize new cel files to an old batch without
   * rerunning the entire batch. 
   * 
   * @param fileName - name of file to write distribution sampled to.
   */
  void saveReferenceProfile(const std::string& fileName) {
    m_FileOutName = fileName;
  }
  void saveReferenceProfile_a5(const std::string& fileName) {
    m_a5_filename = fileName;
  }
  void saveReferenceProfile_group_a5(affx::File5_Group* group5) {
    m_a5_shared_group = group5;
  }

  // set the file to read
  void setreadReferenceProfile(const std::string &RefPro){
	  m_ReferenceProfileName = RefPro;
}

  void setParameters(PsBoard &board);

  void readReferenceProfileFromFile(const std::string& fileName);
  void readReferenceProfileFromFile_a5(const std::string& fileName);
  void readReferenceProfileFromFile_a5_tsv(affx::File5_Group* group5,
                                                        const std::string& name);

void readMultiChannelReferenceProfileFromFile(const std::string& fileName);

  // either compute or read, depending
  void GetReferenceProfile(IntensityMart* iMart);
  
  void updateReference(std::vector<float> &data);
  void updateAllReference(std::vector<std::vector<float> > &adata);
  float TypeIResidual(std::vector<float> &resid, std::vector<float> &data, std::vector<float> &ref);
  void TypeTwoResidual(std::vector<float> &resid, std::vector<float> &data);
  void NewTypeTwoResidual(std::vector<float> &resid, std::vector<float> &data);
  void Winsorize(std::vector<float> &data, std::vector<float> &ref, float Clip);
  int EraseBlemishToReference(std::vector<float> &data, std::vector<float> &ref, std::vector<int> &blemishes, float mt);
  long ThresholdResiduals(std::vector<int> &tmpBlemish, std::vector<float> &tmpRes, float Clip);
  void EmptyBlemishes(std::vector<int> &tmpBlemish, std::vector<float> &tmpRes);
  int EraseBlemishToUnblemished(std::vector<float> &data, std::vector<int> &blemishes);
  int NewEraseBlemishToUnblemished(std::vector<float> &data, std::vector<int> &blemishes);
  void FixProbeIterationList();
  void	FixBlemish(std::vector<std::vector<float> > &adata, float Clip, int chipNo);
  void FixGradient(std::vector<std::vector<float> > &adata, float Clip, int chipNo);
  int BlemishedByProbeset(int, std::vector<int> &blemishes);
	
  void AddLayout(ChipLayout* pobjChipLayout){
	m_pobjChipLayout = pobjChipLayout;
  }

  static std::string getProbesetTrustFileKey() { return  ARTIFACTREDUCTIONTRUSTFILE; }
  static std::string createProbesetTrustTmpFileName(Options * opt);
  std::string m_ProbesetTrustTmpFileName;

protected:
  
  
  /** 
   * Save the current m_ReferenceProfile to a file.
   * @param fileName - where to write the profile.
   */  
  void saveProfileToFile(const std::string &fileName, std::vector<float> &Profile, const std::string &ProfileName);
  void saveReferenceProfileToFile(const std::string& fileName);
  void saveReferenceProfileToFile_a5(const std::string& fileName);
  void saveReferenceProfileToGroup_a5(affx::File5_Group* group5);
   void saveResidualsToFile(const std::string& fileName, std::vector<float> &Residuals, const std::string &ResidualType);
  void saveBlemishMapToFile(const std::string &fileName, std::vector<int> &Profile, const std::string &ProfileName);
 void saveMultiChannelProfileToFile(const std::string& fileName, std::vector<std::vector<float> > &Profile, const std::string &ProfileType); 

  /// Vector representing our target (or average) distribution.
  std::vector<float> m_ReferenceProfile;
  /// All channels for reference distribution
  std::vector< std::vector <float> > m_AllReferenceProfile;
  // cache indexes for fast looping
  std::vector< std::vector <int> > m_ProbeIterationList;
  /// How many things have we on-line averaged
  int m_CurCount;
  // SafetyZero prevents logs of zero
  float m_SafetyZero;

  float m_Clip; // type I residual threshold for outlier
  int m_Open; // how large snow to remove
  int m_Close; // how much closure to use first
  int m_Fringe; // how much fringe to put around continuous regions
  int m_ResType; // what type of residual should I use when computing blemishes
  int m_MapVerbose; // output map of blemishes in what detail
  int m_CoincidenceCount; // how many channels must agree
  std::vector<std::string> m_CelFiles; // what are the cel files in diskmart
  
  affx::TsvFile m_ProbesetTrustTsv;
  bool  m_TrustCheck; // eliminate probesets in blemished areas

  /// Memory that we have to free when we are done.
  std::vector< std::vector<float> * > m_ToFree;
  
  DataStore *m_DataCache;
   std::vector<std::string> m_CacheCelNames;


	int m_dc_k;  // size of block to average over in dc gradient removal
	int m_Gradient; // do gradient removal?

  /// Where to save profile.
  std::string m_FileOutName;
  ///
  std::string m_a5_filename;
  ///
  affx::File5_Group* m_a5_shared_group;

  std::string m_ReferenceProfileName;

/// Pointer to the ChipLayout object
  ChipLayout* m_pobjChipLayout;

  bool m_FreeLayout;

  int m_NumChannels;

  std::string m_DataStoreTempFile;
};


#endif /* _ARTIFACTREDUCTION_H_ */
