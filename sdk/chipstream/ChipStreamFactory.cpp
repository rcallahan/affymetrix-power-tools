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
 * @file   ChipStreamFactory.cpp
 * @author Chuck Sugnet
 * @date   Wed Jan 11 14:37:09 2006
 *
 * @brief
 */

//
#include "chipstream/ChipStreamFactory.h"
//
#include "chipstream/AdapterTypeNormTran.h"
#include "chipstream/ArtifactReduction.h"
#include "chipstream/GcBg.h"
#include "chipstream/IntensityReporter.h"
#include "chipstream/MedNormTran.h"
#include "chipstream/RmaBgTran.h"
#include "chipstream/SketchQuantNormTran.h"
//
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"

#include "util/Fs.h"


using namespace affx;

/**
 * @brief Constructor. Registers the objects we know how to create and how
 * they describe themselves.
 */
ChipStreamFactory::ChipStreamFactory() {
  // 
  m_WriteSketchGroup_a5=NULL;

  // RmaBgTran
  m_Docs.push_back(RmaBgTran::explainSelf());
  m_Creators.push_back(&RmaBgTran::newObject);

  // Quantile Normalization
  m_Docs.push_back(SketchQuantNormTran::explainSelf());
  m_Creators.push_back(&SketchQuantNormTran::newObject);

  // Generic Transformation
  m_Docs.push_back(ArtifactReduction::explainSelf());
  m_Creators.push_back(&ArtifactReduction::newObject);

  // Median Normalization
  m_Docs.push_back(MedNormTran::explainSelf());
  m_Creators.push_back(&MedNormTran::newObject);

  // AdapterType Normalization
  m_Docs.push_back(AdapterTypeNormTran::explainSelf());
  m_Creators.push_back(&AdapterTypeNormTran::newObject);

  // GcBg
  m_Docs.push_back(GcBg::explainSelf());
  m_Creators.push_back(&GcBg::newObject);

  // IntensityReporter
  m_Docs.push_back(IntensityReporter::explainSelf());
  m_Creators.push_back(&IntensityReporter::newObject);

  // ChipStream -- placeholder chipstream that does no tranformations
  m_Docs.push_back(ChipStream::explainSelf());
  m_Creators.push_back(&ChipStream::newObject);

  m_uiCelFileCount = 0;
  m_WriteProfileGroup_a5 = NULL;
}

/**
 * @brief Create a pointer to a new ChipStream object as described
 * in the string specification.
 * @param spec - Specification string i.e. quant-norm.sketch=1000000
 *
 * @return Pointer to new ChipStream objects, must be deleted when
 * finished.
 */
ChipStream *ChipStreamFactory::cStreamForString(const std::string &spec) {
  ChipStream *stream = NULL;
  SelfCreate *create = NULL;
  create = SelfCreate::selfCreateFromString(spec, m_Docs, m_Creators, "ChipStream");
  /* Check class type. */
  if(InstanceOf(create, ChipStream)) {
    stream = static_cast<ChipStream *>(create);
  }
  else {
    Err::errAbort("Class doesn't appear to be of type ChipStream.");
  }
  return stream;
}

/**
 * @brief Factory for making ChipStreams from a string description.
 *
 * @param spec - Specification in form "chipstream.key=value.key=value"
 * @param layout - Probe and probe set info.
 * @param description - name of chipstream path so far.
 * @return ChipStream requested.
 */
ChipStream *ChipStreamFactory::chipStreamForString(const std::string &spec, 
                                                   ChipLayout &layout,
                                                   const std::string &description) {
  ChipStreamFactory factory;
  ChipStream *stream = NULL;

  stream = cStreamForString(spec);
  assert(stream);
  int pmSize = layout.getPmProbeMask().size();
//   if(pmSize == 0)
//     Err::errAbort("ChipStreamFactory::chipStreamForString() - Must have some perfect match probes on chip.");
  vector<bool> pmProbes(pmSize);
  vector<vector<float> > probes;
  Verbose::out(2, ToStr("Made ChipStream of type: ") + ToStr(stream->getType()));
  /* Some ChipStream specific operations. */
  if(InstanceOf(stream, RmaBgTran)) {
    static_cast<RmaBgTran *>(stream)->setPmProbes(layout.getPmProbeMask());
  }
  if(InstanceOf(stream, MedNormTran)) {
    MedNormTran *s = static_cast<MedNormTran *>(stream);
    if(!m_NormProbes.empty()) {
      vector<bool> subset(layout.getProbeCount(), false);
      for(unsigned int i = 0; i < m_NormProbes.size(); i++)
        subset[m_NormProbes[i]] = true;
      s->setSubProbes(subset);
    }
  }
  if (InstanceOf(stream, AdapterTypeNormTran)) 
  {
    AdapterTypeNormTran *s = static_cast<AdapterTypeNormTran *>(stream);
	s->setup(m_strAnnotationFileName, &layout, m_uiCelFileCount);
  }
  if(InstanceOf(stream, SketchQuantNormTran)) {
    SketchQuantNormTran *s = static_cast<SketchQuantNormTran *>(stream);
    s->setProbeCount(layout.getProbeCount());
    /* Do we have a specific sketch we are normalizing to? */
    if(m_TargetSketch.size() > 0) {
      s->setTargetSketch(m_TargetSketch);
    }
    /* Set pm probes if using them for normalization. */
    if(s->getUsePmSubset()) {
      s->setSubProbes(layout.getPmProbeMask());
    }
    /* Else if we have another specific subset use that. */
    else {
      if(!m_NormProbes.empty()) {
        vector<bool> subset(layout.getProbeCount(), false);
        for(unsigned int i = 0; i < m_NormProbes.size(); i++)
          subset[m_NormProbes[i]] = true;
        s->setSubProbes(subset);
      }
    }
    if(!m_WriteSketchDir.empty()) {
      string outFile = Fs::join(m_WriteSketchDir,description + s->getType() +".normalization-target.txt");
      s->saveTargetSketch(outFile);
    }
    if(!m_WriteSketchDir_a5.empty()) {
      string outFile = Fs::join(m_WriteSketchDir_a5,description + s->getType()+".normalization-target.a5");
      s->saveTargetSketch_a5(outFile);
    }
    if (m_WriteSketchGroup_a5!=NULL) {
      s->saveTargetSketch_group_a5(m_WriteSketchGroup_a5);
    }

  }
// Generic data: boy this is too complex 
    if(InstanceOf(stream, ArtifactReduction)) {
    ArtifactReduction *s = static_cast<ArtifactReduction *>(stream);
    
    // only sets up file names for i/o
    // 
    // reading from somewhere, pass filename down into the ArtifactReduction object
    if (!m_ReadProfile.empty()){
	    s->setreadReferenceProfile(m_ReadProfile);
	}
   
    if(!m_WriteProfileDir.empty()) {
      string outFile = m_WriteProfileDir + description + s->getType() + ToStr(".generic-target.txt");
      s->saveReferenceProfile(outFile);
    }

    // make the object aware of structure
    s->AddLayout(&layout);

    // pass layout to object if needed
	// look to other methods


    
    /*if(!m_WriteProfileDir_a5.empty()) {
      string outFile = m_WriteProfileDir_a5 + description + s->getType() + ToStr(".generic-target.a5");
      s->saveReferenceProfile_a5(outFile);
    }
    if (m_WriteProfileGroup_a5!=NULL) {
      s->saveReferenceProfile_group_a5(m_WriteReferenceProfileGroup_a5);
      }*/

  }
// end generic data
    
  if(InstanceOf(stream, GcBg)) {
    if(m_GcControlProbes.empty()) {
      Err::errAbort("chipStreamForString() - Must set control probes to make GcBg Chipstream object.");
    }
    GcBg *gcbg = static_cast<GcBg *>(stream);
    gcbg->setControlProbes(m_GcControlProbes);
    gcbg->setProbeGcVec(layout.getGcProbeVec());
  }
  return stream;
}

ChipStream *ChipStreamFactory::chipStreamForStringBoard(PsBoard &board,
                                                        const std::string &spec, 
                                                        const std::string &description) {
  ChipStreamFactory factory;
  ChipStream *stream = NULL;

  stream = cStreamForString(spec);
  assert(stream);
  //int pmSize = board.getProbeInfo()->getProbeCount();
  Verbose::out(2, ToStr("Made ChipStream of type: ") + ToStr(stream->getType()));
  stream->setParameters(board);
  return stream;
}

void ChipStreamFactory::readNormProbesFromFile(const std::string& fileName) {
  TsvFile tsv;
  vector<int> probeIds;
  int id;

  if(tsv.open(fileName) != TSV_OK) {
    Err::errAbort("Couldn't open '" + fileName + "' to read.");
  }
  tsv.bind(0, "probe_id", &id, affx::TSV_BIND_REQUIRED);
  while(tsv.nextLevel(0) == TSV_OK) {
    probeIds.push_back(TO_ZEROBASED(id));
  }
  setProbesForNorm(probeIds);
  tsv.close();
}

/**
 * Open and read a target distribution from a text file. Must have a column
 * called 'intensities' with distribution quantiles sorted from highest to
 * lowest.
 *
 * @param fileName - text file containing column of quantiles.
 */
void ChipStreamFactory::readTargetSketchFromFile(const std::string& fileName) {
  SketchQuantNormTran::readTargetSketchFromFile(fileName, m_TargetSketch);
}

void ChipStreamFactory::readTargetSketchFromFile_a5(const std::string& fileName) {
  affx::File5_File* file5=new affx::File5_File();
  file5->open(fileName,affx::FILE5_OPEN_RO);
  
  readTargetSketchFromFile_a5_tsv(file5,"target-sketch");
  
  file5->close();
  delete file5;
}
  
void ChipStreamFactory::readTargetSketchFromFile_a5_tsv(affx::File5_Group* group5,
                                                        const std::string& name)
{
  //printf("### reading A5 sketch from: '%s'\n",name.c_str());
  double intensity;
  affx::File5_Tsv* tsv5=group5->openTsv("target-sketch",affx::FILE5_OPEN_RO);

  m_TargetSketch.clear();
  while (tsv5->nextLine()==affx::FILE5_OK) {
    tsv5->get(0,0,&intensity);
    m_TargetSketch.push_back(intensity);
  }
  //
  tsv5->close();
  delete tsv5;
}

bool ChipStreamFactory::canBuild(const std::string &spec) {
  return SelfCreate::canMake(spec, m_Docs);
}
