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
 * @file   SketchQuantNormTran.cpp
 * @author Chuck Sugnet
 * @date   Fri Oct 21 18:16:03 2005
 * 
 * @brief Class for doing normalization. Can do sketch and full quantile (just
 * set sketch to chip size) and supports bioconductor compatibility.
 */

//
#include "chipstream/SketchQuantNormTran.h"
//
#include "chipstream/PsBoard.h"
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Verbose.h"
#include "util/md5sum.h"

using namespace std;
using namespace affx;

/** 
 * Constructor.
 *
 * @param sketchSize - Number of data points to use for
 * approximating full quantile normalization.
 * @param biocCompat - Should ties in rank be done in same manner as
 * bioconductor affy package?
 * @param lowPrecision - Should we truncate values as if we had been written to cel file?
 * @param usePmSubset - Should we use just the perfect match probes for normalization?
 */
SketchQuantNormTran::SketchQuantNormTran(int sketchSize, bool biocCompat, 
                                         bool lowPrecision, bool usePmSubset,
                                         float target, bool doavg) :    
  m_SketchSize(sketchSize), m_BiocCompat(biocCompat), 
  m_LowPrecision(lowPrecision), m_UsePmSubset(usePmSubset),
  m_Target(target), m_DoAverage(doavg) {
  setUpSelfDoc(*this);
  m_Type = getDocName();
  if(sketchSize > 0)
    m_ExtractSketch = new ExtractSketch<std::vector<float>::iterator >(m_SketchSize);
  else 
    m_ExtractSketch = NULL;
  setOptValue("sketch", ToStr(m_SketchSize));
  setOptValue("bioc", m_BiocCompat);
  setOptValue("lowprecision", m_LowPrecision);
  setOptValue("usepm", m_UsePmSubset);
  setOptValue("target", ToStr(m_Target));
  setOptValue("doavg", m_DoAverage);
  m_UsePrecompSketch = false;
  /// @todo Shouldnt there be a for loop here? Otherwise this is constant. -jhg
  md5sum md5;
  md5.final(m_ProbeMd5sum);
  setOptValue("subsetmd5", m_ProbeMd5sum);
  //
  m_a5_filename="";
  m_a5_shared_group=NULL;
  m_SubProbeCount = 0;
}


/**
 * Destructor. 
 */
SketchQuantNormTran::~SketchQuantNormTran() {
  delete m_ExtractSketch;
  for(unsigned int i = 0; i < m_ToFree.size(); i++) {
    delete m_ToFree[i];
  }
}

/** 
 * @brief Only use a subset of all probes for normalization. For
 * example RMA only uses PM probes or might want to just use control
 * genes.
 * @param subsetProbes - bitmask indicating which probes to use for
 * normalization.
 */
void SketchQuantNormTran::setSubProbes(const vector<bool> &subsetProbes) {
  unsigned int i = 0;
  md5sum md5;
  m_SubProbes = subsetProbes;
  m_SubProbeCount = 0;
  for(i = 0; i < m_SubProbes.size(); i++) {
    if(m_SubProbes[i] == true) {
      m_SubProbeCount++;
      // "nbo" put data in network byte order (big endian) for the checksum
      // note that only the index (i) is used, not the actual probesids.
      md5.update_nbo(i);
    }
  }
  md5.final(m_ProbeMd5sum);
  setOptValue("subsetmd5", m_ProbeMd5sum);
  if(m_SketchSize > m_SubProbeCount) 
    m_SketchSize = m_SubProbeCount;
} 

void SketchQuantNormTran::transformDataSuppliedSketch(std::vector<float> &data) {
  transformDataSuppliedSketch(m_Sketches.size() - 1, data);
}

void SketchQuantNormTran::transformDataSuppliedSketch(int index, std::vector<float> &data) {
  if(m_SketchSize != m_AverageSketch.size()) {
    Err::errAbort("SketchQuantNormTran::newChipSuppliedTargetSketch() - Precomputed target sketch (N=" +
                  ToStr(m_AverageSketch.size()) + ") must equal extracted sketch size (N=" +
                  ToStr(m_SketchSize) + ")," +
                  " cel file: " + m_TransformedIMart->getCelFileNames()[index]);
  }
  // @todo - Do we need to keep this around?
  // /* If we're trying to be compatible with exact we have to pretend that we 
  //  * read our sketch from a truncated cel file. */
  // if(m_LowPrecision) {
  //   std::vector<float>::iterator sketch= m_Sketches[index];
  //   for(unsigned int i = 0; i < m_SketchSize; i++) {
  //     *(sketch + i) = Convert::floatLowPrecision( *(sketch + i) );
  //   }
  // }
  
  /* Transform the data to pass it on to the next chipstream. */
  Verbose::out(2, "Passing data on as we are using preset sketch.");
  for(uint32_t i = 0; i < data.size(); i++) {
    data[i] = transform(i, index, data[i]);
  }
}

/** 
 * When we have a precomputed distribution to normalize against
 * we can handle things differently, like passing data through rather
 * than caching it for computing target sketch.
 * @param data - vector of new chip data from same sample.
 */
void SketchQuantNormTran::newChipSuppliedTargetSketch(std::vector<float> &data) {
  unsigned int index = 0;
  //  unsigned int interpIx = 0;


    /* Pull out the probes used to do the normalization. */
    extractSketch(data);
    /* Sanity check for sizes.*/
    if(m_SketchSize != m_AverageSketch.size()) {
      Err::errAbort("SketchQuantNormTran::newChipSuppliedTargetSketch() - Precomputed target sketch (N=" +
                    ToStr(m_AverageSketch.size()) + ") must equal extracted sketch size (N=" +
                    ToStr(m_SketchSize) + ").");
    }
    index = m_Sketches.size() - 1;
    /* If we're trying to be compatible with exact we have to pretend that we 
     * read our sketch from a truncated cel file. */
    if(m_LowPrecision) {
      std::vector<float>::iterator sketch= m_Sketches[index];
      for(unsigned int i = 0; i < m_SketchSize; i++) {
        *(sketch + i) = Convert::floatLowPrecision( *(sketch + i) );
      }
    }

    if (index == 0)
      scaleTargetSketch();

    /* Transform the data to pass it on to the next chipstream. */
    Verbose::out(2, "Passing data on as we are using preset sketch.");
    for(uint32_t i = 0; i < data.size(); i++) {
      /* *sigh* again, if we're trying to be compatible with exact we
       * have to pretend that we read our sketch from a truncated cel
       * file. Sometimes backward compatibility isn't pretty... */
      if(m_LowPrecision) {
        data[i] = Convert::floatLowPrecision(interpolate_qnorm(data[i], 
                        m_Sketches[index], m_Sketches[index]+m_SketchSize, 
                        m_AverageSketch.begin(), m_AverageSketch.end(), 
                        m_PartialSums.begin(), m_PartialSums.end(), 
                        m_BiocCompat));
        if (!Util::isFinite(data[i])) {
            data[i] = Convert::floatLowPrecision(m_AverageSketch[m_AverageSketch.size() - 1]);
        }
      }
      else {
        data[i] = interpolate_qnorm(data[i], 
                        m_Sketches[index], m_Sketches[index]+m_SketchSize, 
                        m_AverageSketch.begin(), m_AverageSketch.end(), 
                        m_PartialSums.begin(), m_PartialSums.end(), 
                        m_BiocCompat);
        if (!Util::isFinite(data[i])) {
          data[i] = m_AverageSketch[m_AverageSketch.size() - 1];
        }
      }
    }
}



void SketchQuantNormTran::extractSketch(std::vector<float> &data) {
  vector<float> copy;
  unsigned int i = 0, currentPm = 0;

  if(m_ExtractSketch == NULL) {
    m_ExtractSketch = new ExtractSketch<vector<float>::iterator >(m_SketchSize);
  }
  /* Sanity checks. */
  if(m_UsePmSubset && m_SubProbes.empty())
    Err::errAbort("Must specify PM subset when using PM probes for normalization.");
  if(m_SubProbes.size() > 0 && data.size() != m_SubProbes.size()) 
    Err::errAbort("SketchQuantNormTran::newChip() - Chip Data different size than Pm Probe Vector");

  // Add a new sketch.
  std::vector<float> *pVec = new vector<float>(m_SketchSize);
  m_ToFree.push_back(pVec);
  m_Sketches.push_back(pVec->begin());
   /* If we're doing a subset of the data, put subset into copy and
     normalize from copy. */
  if(m_SubProbes.size() > 0) {
    copy.resize(m_SubProbeCount);
    /* First filter out PM probes, for RMA only use pm probes. */
    for(i = 0; i < m_SubProbes.size(); i++) {
      if(m_SubProbes[i] == true) {
        copy[currentPm] = data[i];
        currentPm++;
      }
    }
    // Fill in the sketch from the copy
    (*m_ExtractSketch)(copy.begin(), copy.end(), pVec->begin());
  }
  else {
    // Fill in the sketch from the data.
    (*m_ExtractSketch)(data.begin(), data.end(), pVec->begin());
  }
}
  
/** 
 * @brief transform the intensity point supplied coming from a particular
 * probe in a particular microarray.
 * 
 * @param probeIx - Probe index from the cel file.
 * @param chipIx - Set of chips from same sample.
 * @param intensity - Original intensities.
 * @param return - transformed intensities.
 */
void SketchQuantNormTran::transform(int chipIx, std::vector<float>& intensity) {
  if (m_UsePrecompSketch) {
    transformDataSuppliedSketch(chipIx, intensity);
  }
  else {
    /* all the hard work is done just for us by the iterpolater
       function object for this chip. */
    for (int probeIx = 0; probeIx < intensity.size(); probeIx++) {
      intensity[probeIx] = transform(probeIx, chipIx, intensity[probeIx]);
    }
  }
}

float SketchQuantNormTran::transform(int probeIx, int chipIx, float intensity) {
  /* all the hard work is done just for us by the iterpolater
     function object for this chip. */
    
  intensity = interpolate_qnorm(intensity, 
      			  m_Sketches[chipIx], m_Sketches[chipIx]+m_SketchSize, 
      			  m_AverageSketch.begin(), m_AverageSketch.end(), 
      			  m_PartialSums.begin(), m_PartialSums.end(), 
      			  m_BiocCompat, 0.0f, true);
  if(m_LowPrecision) { intensity = Convert::floatLowPrecision(intensity); }
  if(intensity < 0) {
    Err::errAbort("Negative values found when trying to normalize. chip: " + ToStr(chipIx) +
                  " probe: " + ToStr(probeIx + 1) +
                  " cel file: " + m_TransformedIMart->getCelFileNames()[chipIx]);
  }
      if (!Util::isFinite(intensity)) {
    intensity = m_AverageSketch[m_AverageSketch.size() - 1];
  }
  return intensity;
}

void SketchQuantNormTran::newChip(std::vector<float> &data) {
  extractSketch(data);
}

void SketchQuantNormTran::finishedChips() {
  if (!m_UsePrecompSketch) {
    computeTargetSketch();
  }
}

/** 
 * @brief Method for being passed a new cel file worth of data.
 * @param data - Vector of vectors of cel file data.
 */
void SketchQuantNormTran::newDataSet(IntensityMart* iMart) {
  
  if (m_SketchSize < 100)
    Verbose::out(1, "Warning: Are you sure you want sketchsize < 100?");

  int cel_dataset_count = iMart->getCelDataSetCount();
  m_TransformedIMart = iMart;
  // if there aren't any chipstream nodes after this, then don't store
  // all of the intensities
  if (m_Streams.empty()) {
    m_TransformedIMart->setStoreAllCelIntensities(false);
  }

  std::vector<float> data;
  ///@todo Handle multiple sketches for multiple channels
  Verbose::progressBegin(1, "Computing sketch normalization for " + ToStr(cel_dataset_count) + " cel datasets", cel_dataset_count, 0, cel_dataset_count);
  for (int d = 0; d < cel_dataset_count; d++) {
    Verbose::progressStep(1);
    data = iMart->getCelData(d);
    newChip(data);
    if(m_UsePrecompSketch) {
      transformDataSuppliedSketch(data);
      m_TransformedIMart->setProbeIntensity(d, data);
    }
  }
  Verbose::progressEnd(1, "Done.");
  if (!m_UsePrecompSketch) {
    computeTargetSketch();
  }
  if (!m_UsePrecompSketch) {
    Verbose::progressBegin(1, "Applying sketch normalization to " + ToStr(cel_dataset_count) + " cel datasets", cel_dataset_count, 0, cel_dataset_count);
    for (int d = 0; d < cel_dataset_count; d++) {
      Verbose::progressStep(1);
      data = iMart->getCelData(d);
      transform(d, data);
      m_TransformedIMart->setProbeIntensity(d, data);
    }
    Verbose::progressEnd(1, "Done.");
  }

  chipStreamPassNewChip(m_TransformedIMart);
}


/** 
 * Save the current m_AverageSketch to a file.
 * @param fileName - where to write the sketch.
 */  
void SketchQuantNormTran::saveTargetSketchToFile(const std::string& fileName) {
  // printf("### saveTargetSketchToFile('%s')",fileName.c_str());
  //
  affx::TsvFile tsv;
  tsv.addHeader("tsv-file-type","quantile-norm-sketch");
  tsv.defineFile("intensities");
  tsv.setPrecision(0,0,8);
  tsv.writeTsv_v1(fileName);
  //
  for(unsigned int i = 0; i < m_AverageSketch.size(); i++) {
    tsv.set(0,0,m_AverageSketch[i]);
    tsv.writeLevel(0);
  }
  tsv.close();
}

void SketchQuantNormTran::saveTargetSketchToFile_a5(const std::string& fileName) {
  try {
    affx::File5_File* file5=new affx::File5_File();
    file5->open(fileName,affx::FILE5_REPLACE);
    
    saveTargetSketchToGroup_a5(file5);
    
    file5->close();
    delete file5;
  }
  catch(...) {
    Err::errAbort("Cannot open file: " + fileName);
  }
}

void SketchQuantNormTran::saveTargetSketchToGroup_a5(affx::File5_Group* group5) {
  //
  affx::File5_Tsv* tsv5=group5->openTsv("target-sketch",affx::FILE5_CREATE);
  tsv5->defineColumn(0,0,"intensities",affx::FILE5_DTYPE_DOUBLE,0);
  // 
  for(unsigned int i = 0; i < m_AverageSketch.size(); i++) {
    tsv5->set_d(0,0,m_AverageSketch[i]);
    tsv5->writeLevel(0);
  }
  //
  tsv5->close();
  delete tsv5;
}

/** 
 * @brief Signal that no more data is coming (i.e. newChip() will not be
 * called anymore. Currently the class builds up all the sketches and then
 * does the normalization when this function is called. This makes it
 * difficult to pass through data to downstream. If a precomputed sketch
 * is supplied then the data can be passed through chip by chip.
 */
void SketchQuantNormTran::endDataSet() {
  /* Flush data in stream. */
  chipStreamEndChips();

  if(m_FileOutName != "") {
    saveTargetSketchToFile(m_FileOutName);
  }
  if(m_a5_filename!="") {
    saveTargetSketchToFile_a5(m_a5_filename);
  }
  if (m_a5_shared_group!=NULL) {
    saveTargetSketchToGroup_a5(m_a5_shared_group);
  }

  if (m_TransformedIMart != NULL) {
    // intensities have been transformed, so we don't need the
    // sketches anymore
    for(unsigned int i = 0; i < m_ToFree.size(); i++) {
      delete m_ToFree[i];
    }
    m_ToFree.clear();
    m_Sketches.clear();
  }
}

/**
 * @brief Scale the target sketch
 */
void SketchQuantNormTran::scaleTargetSketch() {
  if(m_Target > 0.0) {
    // Compute scaling factor
    double scalingFactor = 1.0; // default to no scaling, reset below
    if(m_DoAverage) {
        double mean = average(m_AverageSketch.begin(), m_AverageSketch.end());
        scalingFactor = m_Target / mean;
        Verbose::out(2,"Mean for quant-norm sketch is " + ToStr(mean));
    } else {
        double med = median(m_AverageSketch.begin(), m_AverageSketch.end());
        scalingFactor = m_Target / med;
        Verbose::out(2,"Median for quant-norm sketch is " + ToStr(med));
    }
    Verbose::out(2,"Scaling factor for quant-norm sketch set to " + ToStr(scalingFactor));

    // scale the target sketch
    for(int i = 0; i < m_AverageSketch.size(); i++) {
        m_AverageSketch[i] *= scalingFactor;
    }

    // update the partial sums
    updatePartialSums();
  }
}
    
/** 
 * Fill in the information for Self documentation.
 * @param doc - Self documenter to be filled in.
 */
void SketchQuantNormTran::setUpSelfDoc(SelfDoc &doc) {
  doc.setDocName(SKETCHQUANTNORMSTR);
  doc.setDocDescription("Class for doing quantile normalization. Can do sketch and full quantile (just set sketch to chip size or zero) and supports bioconductor compatibility.");
  doc.setDocOptions(getDefaultDocOptions());
}

/** 
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> SketchQuantNormTran::getDefaultDocOptions() { 
  std::vector<SelfDoc::Opt> opts;
  SelfDoc::Opt sketch = {"sketch", SelfDoc::Opt::Integer, "-1", "-1", "-1", "NA", 
                         "How many data points from chip to use for normalization (-1 to use default of max of 1% of the chip or 50,000) Use 0 for full quantile normalization."};
  opts.push_back(sketch);
  SelfDoc::Opt bioc = {"bioc", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", 
                       "Set to 'true' for resolving ties in the same order as bioconductor's affy package."};
  opts.push_back(bioc);
  SelfDoc::Opt lowprecision = {"lowprecision", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", 
                               "Set to 'true' to truncate values as seen when writing results to a normalized cel file."};
  opts.push_back(lowprecision);
  SelfDoc::Opt usepm = {"usepm", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                        "Set to true if only using pm probes to do normalization (i.e. rma)"};
  opts.push_back(usepm);
  SelfDoc::Opt target = {"target", SelfDoc::Opt::Double, "0.0", "0.0", "0", "NA", 
                         "Target intensity to set all chips median (or average) to."};
  opts.push_back(target);
  SelfDoc::Opt doavg = {"doavg", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                        "Set to true to do average rather than median."};
  opts.push_back(doavg);
  SelfDoc::Opt probeMd5 = {"subsetmd5", SelfDoc::Opt::String, "", "", "NA", "NA",
                           "Md5sum of the probe ids being used for normalization."};
  opts.push_back(probeMd5);
  return opts;
}


/** 
 * End reading chips when we have to compute the target
 * (m_AverageSketch) ourselves.
 */
void SketchQuantNormTran::computeTargetSketch() {
  vector<vector <float>::iterator> cellDataItVec;
  // Set up function object with coorect types.
  AverageSketch< vector< vector<float>::iterator >::iterator, vector<float>::iterator, double> as;
  m_AverageSketch.clear();
  m_AverageSketch.resize(m_SketchSize);
  // This is the hard work here...
  as(m_Sketches.begin(), m_Sketches.end(), (*m_Sketches.begin()) + m_SketchSize, m_AverageSketch.begin());
 
  updatePartialSums();

  scaleTargetSketch();

//   if(m_Streams.size() > 0) {
//     Err::errAbort("SketchQuantNormTran stream does not currently pass through values unless a target sketch is set.");
//   }
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
SelfCreate *SketchQuantNormTran::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  int sketchSize = 100000;
  bool bioc = false, lowprecision = false, usepm = false, doavg = false;
  float target = 0.0;
  
  fillInValue(sketchSize, "sketch", param, doc);
  fillInValue(bioc, "bioc", param, doc);
  fillInValue(lowprecision, "lowprecision", param, doc);
  fillInValue(usepm, "usepm", param, doc);
  fillInValue(target, "target", param, doc);
  fillInValue(doavg, "doavg", param, doc);
  SketchQuantNormTran *qnorm = new SketchQuantNormTran(sketchSize, bioc, lowprecision, usepm, target, doavg);
  return qnorm;
}

/**
 * Open and read a target distribution from a text file. Must have a column
 * called 'intensities' with distribution quantiles sorted from highest to
 * lowest.
 *
 * @param fileName - text file containing column of quantiles.
 * @param targetSketch - vector to load up with our intensities
 */
void SketchQuantNormTran::readTargetSketchFromFile(const std::string& fileName, std::vector<float> &targetSketch) {
  TsvFile tsv;
  double intensity;
  // file ok
  if (File5_File::isHdf5file(fileName)==1) {
    Err::errAbort("A5 implementation of target-sketch is not yet implemented.");
  }
  int rc = tsv.open(fileName);
  assert(rc==TSV_OK);
  //
  tsv.bind(0,"intensities",&intensity,affx::TSV_BIND_REQUIRED);

  targetSketch.clear();
  while (tsv.nextLevel(0)==TSV_OK) {
    targetSketch.push_back(intensity);
  }

  tsv.close();
}

void SketchQuantNormTran::setParameters(PsBoard &board) { 
  setProbeCount(board.getProbeInfo()->getProbeCount());
  if(getUsePmSubset()) {
    std::vector<bool> pm;
    board.getProbeInfo()->getProbePm(pm);
    setSubProbes(pm);
  }
  if (!board.getOptions()->getOpt("out-dir").empty()) {
    string outFile = board.getOptions()->getOpt("out-dir")  + getType() + ToStr(".normalization-target.txt");
    saveTargetSketch(outFile);
  }
  if (!board.getOptions()->getOpt("a5-sketch").empty()) {
    //    Err::errAbort("What is the best way to support the global a5 via blackboard?");
  }
  if (!board.getOptions()->getOpt("target-sketch").empty()) {
    if (board.getOptions()->getOptBool("a5-sketch")) {
      Err::errAbort("What is the best way to support the global a5 via blackboard?");
    }
    vector<float> sketch;
    string fileName = board.getOptions()->getOpt("target-sketch");
    Verbose::out(1, "Reading sketch from: " + fileName);
    readTargetSketchFromFile(fileName, sketch);
    setTargetSketch(sketch);
  }
  
}
