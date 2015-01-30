////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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
#include "chipstream/SpfReader.h"

using namespace std;
using namespace affx;

/** 
 * Utility function to cut a delimited string into a vector of ints.
 * 
 * @param s - delimited string.
 * @param delim - delimiter to chop on.
 * @param vec - vector to be filled in.
 */
void SpfReader::chopToIntVector(const std::string &s,
                                char delim, 
                                std::vector<int> &vec) {
  vector<string> words;
  Util::chopString(s, delim, words);
  vec.clear();
  vec.resize(words.size());
  for(int i = 0; i < words.size(); i++) {
    vec[i] = Convert::toInt(words[i].c_str());
  }
}


/** 
 * Parse the results from one record into a ProbeList if desired,
 * othewise just keep track of statistics
 * 
 * @param name - Probeset identifier.
 * @param type - Class of probeset (i.e. genotyping, expression, copynumber) as in ProbeSet::Type enumeration
 * @param numBlocks - How many blocks are there in probeset.
 * @param blockSizes - Comma separated size in probes of each block.
 * @param blockAnns - Comma separated annotation for each block as in ProbeSet::BlockAnnotation enumeration
 * @param numMatch - How many matches (1 for pm only, 2 for pm-mm).
 * @param numProbes - Total number of probes.
 * @param probesStr - Comma separated list of probe ids (indexes).
 * @param lineNumber - Line number in file.
 */
ProbeSet * SpfReader::readProbeSetRecord(std::string &name, int type, int numBlocks, std::string &blockSizesStr,
                                         std::string &blockAnnsStr, int numMatch, int numProbes, std::string &probesStr,
                                         int lineNumber) {
  vector<int> blockSizes, blockAnns, probes;

  if(numProbes <= 0) {
    Err::errAbort("Can't have 0 probes for probeset: " + name);
  }
  //  int maxId = m_XCount * m_YCount;
  // convert the block sizes
  chopToIntVector(blockSizesStr, ',', blockSizes);
  if(blockSizes.size() != numBlocks) {
    Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) + 
                  " expecting " + ToStr(numBlocks) + " block sizes got " + ToStr(blockSizes.size()));
  }
  // convert the block annotations
  chopToIntVector(blockAnnsStr, ',', blockAnns);
  if(blockAnns.size() != numBlocks) {
    Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) + 
                  " expecting " + ToStr(numBlocks) + " block annotations got " + ToStr(blockAnns.size()));
  }
  // convert the probes
  chopToIntVector(probesStr, ',', probes);
  if(probes.size() != numProbes) {
    Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) + 
                  " expecting " + ToStr(numProbes) + " probes got " + ToStr(probes.size()));
  }
  if(probes.size() % numMatch != 0) {
    Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) + " must have an integral number * numMatch count of probes.");
  }
  // Subtract 1 of probes as we us R type 1 based index externally
  for(int pIx = 0; pIx < probes.size(); pIx++) {
    if(probes[pIx] != ProbeList::null_probe) {
      probes[pIx]--;
    }
  }
 
  ProbeSet *ps = new ProbeSet();
  ps->name = Util::cloneString(name.c_str());
  ps->psType = (ProbeSet::ProbeSetType)type;
  ps->numGroups = numBlocks;
  ps->atomsPerGroup.resize(numBlocks, 0);
  bool hasMM = numMatch == 2;
  int pmCount = numProbes / numMatch;
  ps->atoms.resize(numBlocks); // make an atom for each block...
  int currentProbeStart = 0;
  for(int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
    Atom *a = new Atom();
    ps->atoms[atomIx] = a;
    int notNullCount = 0;
    int atomProbeEnd = blockSizes[atomIx] + currentProbeStart;
    for(int i = currentProbeStart; i < atomProbeEnd; i++) {
      if(i >= probes.size()) {
        Err::errAbort("Trying to read probe index " + ToStr(i) + " when max is: " + ToStr(probes.size() - 1) + " in probeset " + name);
      }
      if(probes[i] != ProbeList::null_probe) {
        notNullCount++;
      }
      if(hasMM && (i+pmCount >= probes.size())) {
        Err::errAbort("Trying to read probe index " + ToStr(i+pmCount) + " when max is: " + ToStr(probes.size() - 1)  + " in probeset " + name);
      }
      if(hasMM && probes[i+pmCount] != ProbeList::null_probe) {
        notNullCount++;
      }
    } 
    a->probes.resize(notNullCount);
    a->id = atomIx;
    ps->atomsPerGroup[atomIx] = 1; // can do this as there is just one atom per block.
    int probesInserted = 0;
    for(int i = currentProbeStart; i < atomProbeEnd; i++) {
      int probeIndex = i;
      int pId = probes[probeIndex];
      // regular pm-mm pair
      ///@todo I do not think we keep track of st vs at probes -- so we should really
      ///      have a Probe::PM and Probe::MM and set that here
      if(pId != ProbeList::null_probe && hasMM) {
        if(ProbeListFactory::maybeInsertProbe(a, probesInserted, pId, 0, Probe::PMST, pId))
          probesInserted++;
        int mmId = probes[i+pmCount];
        if(ProbeListFactory::maybeInsertProbe(a, probesInserted, mmId, 0, Probe::MMST, pId))
          probesInserted++;
      }
      else if(pId == ProbeList::null_probe && hasMM) {
        int mmId = probes[i+pmCount];
        if(ProbeListFactory::maybeInsertProbe(a, probesInserted, mmId, 0, Probe::MMST, pId))
          probesInserted++;
      }
      else if(pId != ProbeList::null_probe) {
        if(ProbeListFactory::maybeInsertProbe(a, probesInserted, pId, 0, Probe::PMST, pId))
          probesInserted++;
      }
      else {
        Err::errAbort("Null probe with no mismatch in probeset: " + ToStr(ps->name));
      }
      if(probesInserted > a->probes.size()) {
        Err::errAbort("ChipLayout::readProbeSetRecord() - " + ToStr(ps->name) + "More probes than space for.");
      }
    }   
    currentProbeStart += atomProbeEnd - currentProbeStart;
  }
  return ps;
}


ProbeSet * SpfReader::readNextProbeSet() {
  ProbeSet *ps = NULL;
  if(m_IncrementalIn.nextLevel(0) == TSV_OK) {
    ps = readProbeSetRecord(m_IncName, m_IncType, m_IncNumBlocks, m_IncBlockSizesStr, m_IncBlockAnnsStr, 
                            m_IncNumMatch, m_IncNumProbes, m_IncProbesStr, m_IncrementalIn.lineNumber());
  }
  return ps;
}

void SpfReader::openSpf(const std::string &fileName) {
  if(m_IncrementalIn.open(fileName) != TSV_OK)
    Err::errAbort("Couldn't open " + ToStr(fileName) + " to read.");
  m_IncrementalIn.bind(0, "name", &m_IncName, affx::TSV_BIND_REQUIRED);
  m_IncrementalIn.bind(0, "type", &m_IncType, affx::TSV_BIND_REQUIRED);
  m_IncrementalIn.bind(0, "num_blocks", &m_IncNumBlocks, affx::TSV_BIND_REQUIRED);
  m_IncrementalIn.bind(0, "block_sizes", &m_IncBlockSizesStr, affx::TSV_BIND_REQUIRED);
  m_IncrementalIn.bind(0, "block_annotations", &m_IncBlockAnnsStr, affx::TSV_BIND_REQUIRED);
  m_IncrementalIn.bind(0, "num_match", &m_IncNumMatch, affx::TSV_BIND_REQUIRED);
  m_IncrementalIn.bind(0, "num_probes", &m_IncNumProbes, affx::TSV_BIND_REQUIRED);
  m_IncrementalIn.bind(0, "probes", &m_IncProbesStr, affx::TSV_BIND_REQUIRED); 
}
