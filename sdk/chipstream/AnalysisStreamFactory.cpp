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
 * @file   AnalysisStreamFactory.cpp
 * @author Chuck Sugnet
 * @date   Wed Jan 11 15:59:39 2006
 * 
 * @brief  Object for making new analysis streams.
 */

//
#include "chipstream/AnalysisStreamFactory.h"
//
#include "chipstream/AnalysisStreamExpPcaSel.h" 
#include "chipstream/AnalysisStreamExpSpectSel.h" 
#include "chipstream/SketchQuantNormTran.h"
//
#include "file/TsvFile/BgpFile.h"
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/Util.h"

using namespace std;
using namespace affx;

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
void AnalysisStreamFactory::fillInAnalysisStream(AnalysisStream *stream,
                                                 ChipLayout &layout,
                                                 std::string &pmSpec, 
                                                 std::vector<std::string> &chipStreamSpec, 
                                                 std::string &name)
{
  unsigned int wordIx = 0;
  PmAdjuster *adjuster = NULL;
  adjuster = m_PmAdjustFac.pmAdjusterForString(pmSpec, layout);
  stream->setPmAdjuster(adjuster);

  for(wordIx = 0; wordIx < chipStreamSpec.size(); wordIx++) {
    ChipStream *trans = m_CsFac.chipStreamForString(chipStreamSpec[wordIx], layout, name);
    stream->addChipStream(trans);
    name += trans->getType();
    name += ".";
  }
  name += adjuster->getType();
}

void AnalysisStreamFactory::fillInAnalysisStages(AnalysisStream *stream,
                                                 ChipLayout &layout,
                                                 const std::string &pmSpec, 
                                                 std::vector<std::string> &chipStreamSpec, 
                                                 std::string &name)
{
  unsigned int wordIx = 0;
  PmAdjuster *adjuster = NULL;
  adjuster = m_PmAdjustFac.pmAdjusterForString(pmSpec, layout);
  stream->setPmAdjuster(adjuster);

  for(wordIx = 0; wordIx < chipStreamSpec.size(); wordIx++) {
    ChipStream *trans = m_CsFac.chipStreamForString(chipStreamSpec[wordIx], layout, name);
    ChipStreamDataTransform *csStage = new ChipStreamDataTransform(trans);
    stream->m_CSStages.push_back(csStage);
    name += trans->getType();
    name += ".";
  }
}

/*
ChipStreamDataTransform * AnalysisStreamFactory::chipstreamStageFromSpec(AnalysisStream *stream,
                                                                         ChipLayout &layout)

{
    ChipStream *trans = m_CsFac.chipStreamForString(spec, layout, name);
    ChipStreamDataTransform *csStage = new ChipStreamDataTransform(trans);
    return csStage;
}
*/

/** 
 * @brief Factory for created AnalysisStreams from a string description
 * 
 * @param s - String description.
 * @param layout - Probe and probe set info.
 * @param stdMethods - Aliases for standard methods like RMA.
 *
 * @return analysis stream requested.
 */
AnalysisStream *AnalysisStreamFactory::constructAnalysisStream(const std::string& s, ChipLayout &layout, 
                                                               std::map<std::string,std::string> &stdMethods,
                                                               std::string analysisName) {
  assert(s!="");
  string description;
  bool doingStdMethod = false;
  vector<string> words;
  string name;
  QuantMethod *qMethod = NULL;
  AnalysisStream *stream = new AnalysisStream();
  
  /* Check for one of our std methods. */
  if(stdMethods.find(s) != stdMethods.end()) {
    doingStdMethod = true;
    description = stdMethods[s];
  }
  else
    description = s;

  Util::chopString(description, ',', words);
  if(words.size() < 2) {
    Err::errAbort("Must specify at least a pm adjustment and summary type.");
  }
  vector<string> chipStreamSpec(words.begin(), words.end() - 2);
  fillInAnalysisStream(stream, layout, words[words.size() -2], chipStreamSpec, name);
  /* Set the quantification method. */
  qMethod = m_QuantMethFac.quantMethodForString(words[words.size() - 1], layout, m_QuantType);
  stream->setQuantMethod(qMethod);
  name += ".";
  name += qMethod->getType();
  stream->setName(name);
  /* Make sure we use the alias name for user's sanity. */
  if(doingStdMethod)
    stream->setName(s);
  if(analysisName != "")
    stream->setName(analysisName);
  return stream;
}

AnalysisStreamExpression *
AnalysisStreamFactory::tryToMakeAnalysisStream(const std::string &spec) {
  string s = spec;
  AnalysisStreamExpression *stream = NULL;
  SelfCreate *create = NULL;
  create = SelfCreate::selfCreateFromString(s, m_Docs, m_Creators, "AnalysisStreamExpression", false);
  if(create != NULL) {
    if(InstanceOf(create, AnalysisStreamExpression)) {
      stream = static_cast<AnalysisStreamExpression *>(create);
    }
    else {
      Err::errAbort("Class doesn't appear to be of type AnalysisStreamExpression.");
    }
  }
  return stream;
}

/** 
 * @brief Factory for creatng mRNA expression AnalysisStream from a string description.
 * 
 * @param s - String description.
 * @param layout - Probe and probe set info.
 * @param stdMethods - Aliases for standard methods like RMA.
 *
 * @return analysis stream requested.
 */
AnalysisStreamExpression *
AnalysisStreamFactory::constructExpressionAnalysisStream(const std::string& s, ChipLayout &layout, 
                                                         std::map<std::string,std::string> &stdMethods,
                                                         std::string analysisName) {
  assert(s!="");
  string description;
  bool doingStdMethod = false;
  vector<string> words;
  string name;
  bool appendName = false;
  QuantExprMethod *qeMethod = NULL;
  //  AnalysisStreamExpression *stream = new AnalysisStreamExpression();
  AnalysisStreamExpression *stream = NULL;
  
  /* Check for one of our std methods. */
  if(stdMethods.find(s) != stdMethods.end()) {
    doingStdMethod = true;
    description = stdMethods[s];
  }
  else
    description = s;

  Util::chopString(description, ',', words);
  if(words.size() < 2) {
    Err::errAbort("Must specify at least a pm adjustment and summary type.");
  }
  vector<string>::reverse_iterator rIter = words.rbegin();
  stream = tryToMakeAnalysisStream(*rIter);
  if(stream != NULL) {
    words.pop_back();
    appendName = true;
  }
  else {
    stream = new AnalysisStreamExpression();
  }
  vector<string> chipStreamSpec(words.begin(), words.end() - 2);
  fillInAnalysisStream(stream, layout, words[words.size() -2], chipStreamSpec, name);

  /* Set the quantification method. */
  qeMethod = m_QuantMethFac.quantExprMethodForString(words[words.size() - 1], layout, m_QuantType);
  stream->setQuantMethod(qeMethod);
  name += ".";
  name += qeMethod->getType();
  if(appendName) {
    name += ".";
    name += stream->getDocName();
  }
  stream->setName(name);
  /* Make sure we use the alias name for user's sanity. */
  if(doingStdMethod)
    stream->setName(s);
  if(analysisName != "")
    stream->setName(analysisName);

  // Setup a separate chip stream for getting probe intensities 
  // when doing PCA select -- ie do not use the user specified chip stream
  if(InstanceOf(stream,AnalysisStreamExpPcaSel)) {
    AnalysisStreamExpPcaSel *pcaStream = static_cast<AnalysisStreamExpPcaSel *>(stream);
    string name = stream->getName() + ".self-qnorm";
    string qSpec = "quant-norm.sketch=" + ToStr(m_CsFac.getTargetSketchSize());
    ChipStream *cstream = m_CsFac.chipStreamForString(qSpec, layout, name);
    SketchQuantNormTran *qnorm = dynamic_cast<SketchQuantNormTran *>(cstream);
    pcaStream->setQuantNorm(qnorm);
  }

  return stream;
}

/** 
 * @brief Factory for creatng mRNA expression AnalysisStream from a string description.
 * 
 * @param s - String description.
 * @param layout - Probe and probe set info.
 * @param stdMethods - Aliases for standard methods like RMA.
 *
 * @return analysis stream requested.
 */
AnalysisStreamExpression *
AnalysisStreamFactory::constructExpressionAnalysisStages(const std::string& s, ChipLayout &layout, 
                                                         std::map<std::string,std::string> &stdMethods,
                                                         std::string analysisName) {
  assert(s!="");
  string description;
  bool doingStdMethod = false;
  vector<string> words;
  string name;
  bool appendName = false;
  QuantExprMethod *qeMethod = NULL;
  AnalysisStreamExpression *stream = NULL;
  
  /* Check for one of our std methods. */
  if(stdMethods.find(s) != stdMethods.end()) {
    doingStdMethod = true;
    description = stdMethods[s];
  }
  else
    description = s;

  Util::chopString(description, ',', words);
  if(words.size() < 2) {
    Err::errAbort("Must specify at least a pm adjustment and summary type.");
  }

  /* See if there are parameters for a special type of AnalysisStream, 
     otherwise make a basic version */
  vector<string>::reverse_iterator rIter = words.rbegin();
  stream = tryToMakeAnalysisStream(*rIter);
  if(stream != NULL) {
    words.pop_back();
    appendName = true;
  }
  else {
    stream = new AnalysisStreamExpression();
  }

  vector<string> chipStreamSpec(words.begin(), words.end() - 2);
  fillInAnalysisStages(stream, layout, words[words.size() - 2], chipStreamSpec, name);

  /* Set the quantification method. */
  qeMethod = m_QuantMethFac.quantExprMethodForString(words[words.size() - 1], layout, m_QuantType);
  stream->setQuantMethod(qeMethod);
  name += ".";
  name += qeMethod->getType();
  if(appendName) {
    name += ".";
    name += stream->getDocName();
  }
  stream->setName(name);
  /* Make sure we use the alias name for user's sanity. */
  if(doingStdMethod)
    stream->setName(s);
  if(analysisName != "")
    stream->setName(analysisName);

  // Setup a separate chip stream for getting probe intensities 
  // when doing PCA select -- ie do not use the user specified chip stream
  if(InstanceOf(stream,AnalysisStreamExpPcaSel)) {
    AnalysisStreamExpPcaSel *pcaStream = static_cast<AnalysisStreamExpPcaSel *>(stream);
    string name = stream->getName() + ".self-qnorm";
    ChipStream *cstream = m_CsFac.chipStreamForString("quant-norm.sketch=-1", layout, name);
    SketchQuantNormTran *qnorm = dynamic_cast<SketchQuantNormTran *>(cstream);
    pcaStream->setQuantNorm(qnorm);
  }

  return stream;
}

/** 
 * @brief Factory for creatng Genotype SNP chip AnalysisStream from a string description.
 * 
 * @param s - String description.
 * @param layout - Probe and probe set info.
 * @param stdMethods - Aliases for standard methods like RMA.
 *
 * @return analysis stream requested.
 */
AnalysisStreamGType *
AnalysisStreamFactory::constructGTypeAnalysisStream(const std::string& s, ChipLayout &layout, 
                                                    std::map<std::string,std::string> &stdMethods,
                                                    std::string analysisName) {
  assert(s!="");
  string description;
  bool doingStdMethod = false;
  vector<string> words;
  string name;
  
  QuantGTypeMethod *gtMethod = NULL;
  AnalysisStreamGType *stream = new AnalysisStreamGType();
  
  /* Check for one of our std methods. */
  if(stdMethods.find(s) != stdMethods.end()) {
    doingStdMethod = true;
    description = stdMethods[s];
  }
  else
    description = s;

  Util::chopString(description, ',', words);
  if(words.size() < 2) {
    Err::errAbort("Must specify at least a pm adjustment and summary type.");
  }
  vector<string> chipStreamSpec(words.begin(), words.end() - 2);
  fillInAnalysisStream(stream, layout, words[words.size() -2], chipStreamSpec, name);

  /* Set the quantification method. */
  gtMethod = m_QuantMethFac.quantGTypeMethodForString(words[words.size() - 1], layout, m_QuantType);
  stream->setQuantMethod(gtMethod);
  name += ".";
  name += gtMethod->getType();
  stream->setName(name);
  /* Make sure we use the alias name for user's sanity. */
  if(doingStdMethod)
    stream->setName(s);
  if(analysisName != "")
    stream->setName(analysisName);
  return stream;
}

/** 
 * @brief Load up a list of probes from a .bgp file
 * 
 * @param fileName - name of bgp file.
 * @param controlProbes - vector to be filled in with probes.
 */
void AnalysisStreamFactory::probeListFromBgpFile(const std::string& fileName, vector<Probe *> &controlProbes) {

  affx::BgpFile bgp;
  if (bgp.open(fileName)!=TSV_OK) {
    Err::errAbort("Unable to open "+ToStr(fileName));
  }
  
  while (bgp.next_bgprobe()==TSV_OK) {
    Probe *p = new Probe();
    p->id=bgp.probe_id-1; // NOTE: -1 because this probe_id is 1-based.
    controlProbes.push_back(p);
  }
  bgp.close();
}

void AnalysisStreamFactory::setupAnalysisStreamFactory(
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
     ) {
  asFactory.setCelFileCount(celFileCount);
  if(writeSketch) {
    asFactory.setWriteSketchDir(outDir + Fs::osPathSep());
  } else if (a5Sketch) {
    if (a5SketchUseGlobal) {
        if(outputGroup == NULL)
            Err::errAbort("--a5-sketch-use-global option given, but no global file. Must specify --a5-global-file.");
        asFactory.setWriteSketchGroup_a5(outputGroup);
    }
    else {
      asFactory.setWriteSketchDir_a5(outDir + Fs::osPathSep());
    }
  }

  if (writeProfile){
      asFactory.setWriteProfileDir(outDir + Fs::osPathSep());
  }
  //
  if (targetSketch != "") {
    Verbose::out(1, "Opening target normalization file: " + Fs::basename(targetSketch));
    asFactory.readTargetSketchFromFile(targetSketch);
  } else if(a5SketchInputGlobal || a5SketchInputFile!=""){
      string groupName = "/";
      string dataName = "target-sketch";
      if(a5SketchInputName != "") 
        dataName = a5SketchInputName;
      if(a5SketchInputGroup != "") 
        groupName = a5SketchInputGroup;
      else if(a5InputGroup != "") 
        groupName = a5InputGroup;
      else if(a5Group != "") 
        groupName = a5Group;

      if (a5SketchInputGlobal) {
          Verbose::out(1,"Loading sketch from global A5 file, group '" + groupName + "', data '" + dataName + "'");
          if(inputFile == NULL)
            Err::errAbort("--a5-sketch-input-global option given, but no global input file. Must specify --a5-global-file.");
          affx::File5_Group *a5group = inputFile->openGroup(groupName,affx::FILE5_OPEN_RO);
          asFactory.readTargetSketchFromFile_a5_tsv(a5group, dataName);
          a5group->close();
          Freez(a5group);
      } else if (a5SketchInputFile!="") {
          Verbose::out(1,"Loading sketch from '" + a5SketchInputFile + "' A5 file, group '" + groupName + "', data '" + dataName + "'");
          affx::File5_File  *a5file = new affx::File5_File();
          a5file->open(a5SketchInputFile,affx::FILE5_OPEN_RO);
          affx::File5_Group *a5group = a5file->openGroup(groupName,affx::FILE5_OPEN_RO);
          asFactory.readTargetSketchFromFile_a5_tsv(a5group, dataName);
          a5group->close();
          Freez(a5group);
          a5file->close();
          Freez(a5file);
      }

  }

	// the factory needs to know the file name so it can pass it to the created objects
  if (referenceProfile !=""){
	  asFactory.setreadReferenceProfile(referenceProfile);
}

  if(bgpFile!="") {
    asFactory.readControlProbes(bgpFile);
  }
  if (AnnotFile != "") {
    asFactory.setAnnotationFileName(AnnotFile);
  }

  QuantMethodFactory::setupQuantMethodFactory(
        asFactory.m_QuantMethFac,
        inputFile,
        inputGroup,
        outputFile,
        outputGroup,
        probeCount,
        outDir,
        useFeatEff,
        a5FeatureEffectsInputGlobal,
        a5FeatureEffectsInputFile,
        a5FeatureEffectsInputName,
        a5FeatureEffectsInputGroup,
        a5InputGroup,
        a5Group,
        setAnalysisName,
        quantMethodSpec,
        layout
  );
}

