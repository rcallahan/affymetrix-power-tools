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

//
#include "chipstream/ChipLayout.h"
//
#include "chipstream/EngineUtil.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "file/TsvFile/PgfFile.h"
#include "stats/stats.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/RowFile.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affx;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;

#define MAX_PROBES 65536
#define MAX_ATOMS 65536
#define NODATA -1

/** Clean up memory in hashes. */
ChipLayout::~ChipLayout() {
  /* Used to iterate through hashes, but it was
     very slow as if the memory in ->first is deleted then
     can't iterate to next element and calling begin() lots
     of times and erase() on each element was sloooow. */
  //m_PsNameHash.clear();
}

/** Constructor. */
ChipLayout::ChipLayout() {
  m_NumExpression = 0;
  m_NumGenotyping = 0;
  m_NeedMismatch = true;
  m_NeedGc = true;
  m_numChannels = 1;
  m_NeedPmAlleleMatch = false;
  m_haveKillList = false;
  m_warningCount = 0;
  m_MaxNameLength = 0;
  m_XCount = 0;
  m_YCount = 0;
}

/**
 * Write out a single record to the tsv file supplied.
 *
 * @param tsv - TsvFile opened and bound.
 * @param pL - ProbeList record to write out.
 */
void
ChipLayout::writeSpfProbeSetRecord(SpfFile &tsv,const ProbeSet &ps)
{
  APT_ERR_ASSERT(tsv.m_spf_format==2,"Only format v2 is supported for now.");

  // block counts;
  string block, blockAnn;

  // Setup the block sizes
  int bIx = 0;
  for (bIx = 0; bIx < ps.numGroups - 1; bIx++) {
    block += ToStr(ps.atomsPerGroup[bIx]) + ",";
  }
  block += ToStr(ps.atomsPerGroup[bIx]);

  // Figure out the block annotations based on type and block number number
  vector<int> blockAnnVec(ps.numGroups, ProbeList::NoAnnotation);
  if (ps.numGroups == 4 && 
     (ps.psType == ProbeSet::GenoType)) {
    blockAnnVec[0] =  ProbeList::A_AlleleForward;
    blockAnnVec[1] = ProbeList::B_AlleleForward;
    blockAnnVec[2] = ProbeList::A_AlleleReverse;
    blockAnnVec[3] = ProbeList::B_AlleleReverse;
  }
  else if (ps.numGroups == 2 && 
          (ps.psType == ProbeSet::GenoType)) {
    blockAnnVec[0] =  ProbeList::A_Allele;
    blockAnnVec[1] =  ProbeList::B_Allele;
  }
  for (bIx = 0; bIx < ps.numGroups - 1; bIx++) {
    blockAnn += ToStr(blockAnnVec[bIx]) + ",";
  }
  blockAnn += ToStr(blockAnnVec[bIx]);

  int pmCount = ProbeListFactory::countPmProbes(&ps);
  int probeCount = pmCount;
  bool hasMM = ProbeListFactory::hasNonPm(ps);
  if (hasMM) {
    probeCount = 2*probeCount;
  }

  string probes;
  int pIx = 0;
  vector<int> probeVec(probeCount, ProbeList::null_probe);
  int pCount = 0;
  vector<int> atomProbes(ps.atoms.size(), 0);
  // fill in the probes.
  for (uint32_t atomIx = 0; atomIx < ps.atoms.size(); atomIx++) {
    Atom *a = ps.atoms[atomIx];
    for (uint32_t probeIx = 0; probeIx < a->probes.size(); ++probeIx) {
      Probe *current = a->probes[probeIx];
      Probe *next = NULL;
      atomProbes[atomIx]++;
      if (probeIx +1 < a->probes.size()) 
        next = a->probes[probeIx + 1];
      // if it is an pm-mm set then put pm in front and mm in matching position.
      if (Probe::isPm(*current) && next != NULL && Probe::isMm(*next)) {
        probeVec[pCount] =  current->id;
        probeVec[pCount + pmCount] = next->id;
        probeIx++; // advance probeIx count as we did two.
      }
      // if it is a lone mm then put a "null_probe" in for the pm probe.
      if (!Probe::isPm(*current)) { 
        probeVec[pCount] = ProbeList::null_probe;
        probeVec[pCount + pmCount] =  current->id;
      }
      // lone pm probe
      else {
        probeVec[pCount] =  current->id;
      }
      pCount++;
    }
  }
  for (pIx = 0; pIx < probeVec.size() -1; pIx++) {
    int oneBase = probeVec[pIx] != ProbeList::null_probe ? 1 : 0;
    probes += ToStr(probeVec[pIx] + oneBase) + ",";
  }
  int oneBase = probeVec[pIx] != ProbeList::null_probe ? 1 : 0;
  probes += ToStr(probeVec[pIx] + oneBase);

  // write the SPF record
  int count = 0;
  tsv.set(0, count++, ps.name);
  //
  tsv.set(0, count++, ps.psType);
  tsv.set(0, count++, ps.numGroups);
  tsv.set(0, count++, block);
  tsv.set(0, count++, blockAnn);
  //
  if (hasMM) 
    tsv.set(0, count++, 2);
  else 
    tsv.set(0, count++, 1);
  //
  tsv.set(0, count++, probeCount);  
  tsv.set(0, count++, probes);
  //
  tsv.writeLevel(0);
}

/** 
 * Write out a single record to the tsv file supplied.
 * 
 * @param tsv - TsvFile opened and bound.
 * @param pL - ProbeList record to write out.
 */
void ChipLayout::writeSpfProbeListRecord(SpfFile &tsv,const ProbeListPacked &pL)
{
  APT_ERR_ASSERT(tsv.m_spf_format==2,"Only format v2 is supported for now.");

  // block counts;
  string blockSize, blockAnn;
  int bIx = 0;
  for (bIx = 0; bIx < pL.block_cnt() - 1; bIx++) {
    blockSize += ToStr(pL.get_blockSize(bIx)) + ",";
    blockAnn += ToStr(pL.get_blockAnn(bIx)) + ",";
  }
  blockSize += ToStr(pL.get_blockSize(bIx));
  blockAnn += ToStr(pL.get_blockAnn(bIx));
  string probes;
  int pIx = 0;

  // @todo harley Fix this up with ProbeList::toZero
  // probe index is 1 based rather than 0 based internals so add 1 if not a null probe
  for (pIx = 0; pIx < pL.probe_cnt() - 1; pIx++) {
    int oneBase = pL.get_probeId(pIx) != ProbeList::null_probe ? 1 : 0;
    probes += ToStr(pL.get_probeId(pIx) + oneBase) + ",";
  }
  int oneBase = pL.get_probeId(pIx) != ProbeList::null_probe ? 1 : 0;
  probes += ToStr(pL.get_probeId(pIx) + oneBase);

  // write the SPF record
  // name type num_blocks block_sizes block_annotations num_match num_probes probes
  int count = 0;
  tsv.set(0, count++, pL.get_name_string()); // name
  tsv.set(0, count++, pL.get_type());        // type
  tsv.set(0, count++, pL.block_cnt());       // num_blocks
  tsv.set(0, count++, blockSize);            // block_sizes
  tsv.set(0, count++, blockAnn);             // block_annotations
  tsv.set(0, count++, pL.get_numMatch());    // num_match
  tsv.set(0, count++, pL.probe_cnt());       // num_probes
  tsv.set(0, count++, probes);               // probes
  //
  tsv.writeLevel(0);
}

void ChipLayout::openSpfForWrite(affx::SpfFile& spffile,const std::string& fileName,int spfFormat) {
  // the standard set of columns
  spffile.define_file(spfFormat);
  // stuff what will appear in the headers into the vars...
  spffile.m_header_chipTypes.clear();
  spffile.m_header_chipTypes.push_back(m_Header["chip_type"][0]);
  spffile.m_header_numProbesets=getProbeSetCount();
  spffile.m_header_numCols=getXCount();
  spffile.m_header_numRows=getYCount();
  spffile.m_header_numChannels=numChannels();
  // ...and this will add the headers in the standard order.
  spffile.addStandardHeaders();
  // ready for data.
  spffile.writeSpf(fileName);
}

/** 
 * Write out our probelist information in our simple probe format
 * 
 * @param fileName - path to write file out to
 */
void ChipLayout::writeSpfProbeList(const std::string& fileName,int spfFormat) {
  SpfFile spffile;
  openSpfForWrite(spffile,fileName,spfFormat);
  int i_max=m_PlFactory.getProbeSetCount();
  for (int i = 0; i < i_max; i++) {
    writeSpfProbeListRecord(spffile,m_PlFactory.getProbeListAtIndex(i));
  }
  spffile.close();
}


ProbeListPacked
ChipLayout::makeProbeListPacked(const std::string& name,
                                int type, int numBlocks,
                                const std::vector<int>& blockSizes,
                                const std::vector<int>& blockAnns,
                                const std::vector<int>& blockAlleles,
                                const std::vector<int>& blockContexts,
                                const std::vector<int>& blockChannels,
                                const std::vector<int>& blockRepTypes,
                                int numMatch,
                                int numProbes,
                                const std::vector<int>& probes)
{
  // this shouldnt happen
  APT_ERR_ASSERT((numMatch==1||numMatch==2),"numMatch must be 1 or 2.");
  APT_ERR_ASSERT(numProbes>0,"Most have at least one probe.");

  //
  ProbeListPacked pl = m_PlFactory.add_ProbeList(numBlocks, numProbes, name);
  //
  pl.set_type(type);
  pl.set_numMatch(numMatch);
  //
  APT_ERR_ASSERT(numBlocks==blockSizes.size()   ,"internal error.");
  APT_ERR_ASSERT(numBlocks==blockAnns.size()    ,"internal error.");
  APT_ERR_ASSERT(numBlocks==blockAlleles.size() ,"internal error.");
  APT_ERR_ASSERT(numBlocks==blockContexts.size(),"internal error.");
  APT_ERR_ASSERT(numBlocks==blockChannels.size(),"internal error.");
  APT_ERR_ASSERT(numBlocks==blockRepTypes.size(),"internal error.");
  //
  for (int i = 0; i < numBlocks; i++) {
    pl.set_blockSize(i, blockSizes[i]);
    pl.set_blockAnn(i, blockAnns[i]);
    pl.set_blockAllele(i,blockAlleles[i]);
    pl.set_blockContext(i,blockContexts[i]);
    pl.set_blockChannel(i,blockChannels[i]);
    pl.set_blockRepType(i,blockRepTypes[i]);
  }
  //
  for (int i = 0; i < numProbes; i++) {
    pl.set_probeId(i, probes[i]);
  }
  return pl;
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
 * @param keep - Should this probelist be kept after getting statistics.
 * @param probeSubSet - Bool array with probes to be loaded marked as true.
 */

void ChipLayout::parseProbeListRecord(const std::string& name,
                                      int type,
                                      int numBlocks,
                                      const std::string& blockSizesStr,
                                      const std::string& blockAnnsStr,
                                      bool have_blockAlleleStr,
                                      const std::string& blockAlleleStr,
                                      bool have_blockContextStr,
                                      const std::string& blockContextStr,
                                      bool have_blockChannelStr,
                                      const std::string& blockChannelStr,
                                      bool have_blockRepTypeStr,
                                      const std::string& blockRepTypeStr,
                                      int numMatch,
                                      int numProbes,
                                      const std::string& probesStr,
                                      int lineNumber,
                                      bool keep,
                                      std::vector<bool> &probeSubSet, 
                                      bool keepOrder)
{
  vector<int> blockSizes, blockAnns, blockAlleles, blockContexts, blockChannels, blockRepTypes, probes;
  //int maxId = m_XCount * m_YCount;
  int maxId = getProbeCount();
  
  // convert the block sizes
  Convert::strToIntVec(blockSizesStr, ',', blockSizes);
  if (blockSizes.size() != numBlocks) {
    Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                  " expecting " + ToStr(numBlocks) + " block sizes got " +
                  ToStr(blockSizes.size()));
  }
  // convert the block annotations
  Convert::strToIntVec(blockAnnsStr, ',', blockAnns);
  if (blockAnns.size() != numBlocks) {
    Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                  " expecting " + ToStr(numBlocks) + " block annotations got " +
                  ToStr(blockAnns.size()));
  }
  //
  if (have_blockAlleleStr) {
    Convert::strToIntVec(blockAlleleStr, ',', blockAlleles);
    if (blockAlleles.size() != numBlocks) {
      Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                    " expecting " + ToStr(numBlocks) + " block alleles got " +
                    ToStr(blockAlleles.size())+
                    " from: '"+blockAlleleStr+"'.");
    } 
  }
  else {
    blockAlleles.resize(0);  
    blockAlleles.resize(numBlocks,0);
    if (type == ProbeSet::GenoType && (numBlocks == 2 || numBlocks == 4)) {
      // Unless otherwise specified, alleles for ProbeSet::GenoType
      // blocks are determined by their position in the list of
      // blocks: 1st and 3rd blocks are allele A, 2nd and 4th blocks
      // are allele B.
      for (int i = 1; i < numBlocks; i+=2) {
        blockAlleles[i] = 1;
      }
    }
  }

  //
  if (have_blockContextStr) {
    Convert::strToIntVec(blockContextStr, ',', blockContexts);
    if (blockContexts.size() != numBlocks) {
      Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                    " expecting " + ToStr(numBlocks) + " block contexts got " +
                    ToStr(blockContexts.size()) +
                    " from '"+blockContextStr+"'.");
    }
  }
  else {
    blockContexts.resize(0);
    blockContexts.resize(numBlocks,0);
  }
  //
  if (have_blockChannelStr) {
    Convert::strToIntVec(blockChannelStr, ',', blockChannels);
    if (blockChannels.size() != numBlocks) {
      Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                    " expecting " + ToStr(numBlocks) + " block contexts got " +
                    ToStr(blockChannels.size()) +
                    " from '"+blockChannelStr+"'.");
    }
  }
  else {
    blockChannels.resize(0);
    blockChannels.resize(numBlocks,0);
  }
  //
  if (have_blockRepTypeStr) {
    Convert::strToIntVec(blockRepTypeStr, ',', blockRepTypes);
    if (blockRepTypes.size() != numBlocks) {
      Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                    " expecting " + ToStr(numBlocks) + " block contexts got " +
                    ToStr(blockRepTypes.size()) +
                    " from '"+blockRepTypeStr+"'.");
    }
  }
  else {
    blockRepTypes.resize(0);
    blockRepTypes.resize(numBlocks,0);
  }

  // convert the probes
  Convert::strToIntVec(probesStr, ',', probes);
  if (probes.size() != numProbes) {
    Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                  " expecting " + ToStr(numProbes) + " probes got " + ToStr(probes.size()));
  }
  if (probes.size() % numMatch != 0) {
    Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                  " must have an integral number * numMatch count of probes.");
  }

  // Subtract 1 of probes as we use R type 1 based index externally
  for (int pIx = 0; pIx < probes.size(); pIx++) {
    if (probes[pIx] != ProbeList::null_probe) {
      probes[pIx]--;
    }
    // check to make sure that we got a reasonable probe index
    if (probes[pIx] > maxId ||
       (probes[pIx] != ProbeList::null_probe && probes[pIx] < 0)) {
      Err::errAbort("Probeset: '" + name + "' on line: " + ToStr(lineNumber) +
                    " expecting all probes to be greater than 0 and less than " + ToStr(maxId) +
                    ". got: " + ToStr(probes[pIx]));
    }
  }
  /* Remember the type in case this probeset isn't actually
     loaded. Convert type to expression if a "false" genotype
     probeset. */
  // Ray added this to cover up a later hack.
  /*
    if (m_PsTypes.empty()) {
    if (type == ProbeSet::GenoType && numBlocks != 2 && numBlocks != 4) {
    m_PsTypes.push_back((char)ProbeSet::Expression);
    }
    else {
    m_PsTypes.push_back((char)type);
    }
    }
  */
  bool multAlleles = ! ChipLayout::sameValues(blockAlleles);
  trackGlobalStats(name, type, numBlocks, blockSizes, numMatch, probes, multAlleles);

  if (keepOrder) {
    // mark probes to be loaded.
    for (int pIx = 0; pIx < probes.size(); pIx++) {
      m_ProbeLayoutOrder.push_back(probes[pIx]);
    }
  }
  
  if (keep) {
    for (int pIx = 0; pIx < probes.size(); pIx++) {
      if (probes[pIx] != ProbeList::null_probe) {
        probeSubSet[probes[pIx]] = true;
      }
    }
    // save probeset
    ProbeListPacked pl = makeProbeListPacked(name, type, numBlocks, blockSizes,
                                             blockAnns,
                                             blockAlleles,
                                             blockContexts,
                                             blockChannels,
                                             blockRepTypes,
                                             numMatch, numProbes, probes);
    // we are using the vec in the factory.
    // m_PlVec.push_back(pl);
  }
}

////////////////////

void checkTsvChipType(affx::TsvFile& tsv,const std::string& chipType)
{
  // blank chipType matches all...
  if (chipType=="") {
    return;
  }

  // found a match!
  if (tsv.hasHeaderEqualTo("chip_type",chipType)==affx::TSV_OK) {
    return;
  }

  //// didnt find a match tell the user why...
  //  no header...
  std::string val;
  if (tsv.getHeader("chip_type", val)!=TSV_OK) {
    Err::errAbort("chip_type is a required header in file: '"+tsv.getFileName()+"'");
  }
  // there was a mismatch...
  Verbose::out(1,"This file supports the following:");
  tsv.headersBegin();
  while (tsv.headersFindNext("chip_type",val)==TSV_OK) {
    Verbose::out(1,"   #%chip_type="+val);
  }
  // signal the error
  Err::errAbort("No match for chip_type='"+chipType+"' in '"+tsv.getFileName()+"'");
}

/**
 * Read in probelist data from a simple probe format file.
 *
 * @param fileName - Name of file to be opened.
 * @param probeSetsToLoad - Which probe sets should be loaded?
 * @param probesetNames - If not null filled in with the
 * probesetname for each probeset even if not actually into
 * chiplayout. Useful for chp files where every single probeset must
 * be included and in order.
 * @param probeSubset - Subset of probes to be loaded.
 * @param chipTypeExpected - What sort of chip should this cdf file be for?
 * @param justStats - just read and generate stats, don't load into memory.
 * @param psTypesToLoad - What types of probesets to load into memory.
 */
void ChipLayout::openSpf(const std::string& fileName,
                         const std::set<const char *, 
                         Util::ltstr> &probeSetsToLoad,
                         std::vector<const char *> *probesetNames,
                         std::vector<bool> &probeSubSet,
                         const std::string &chipTypeExpected,
                         bool justStats,
                         std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad) {

  // what format is this spf file?
  std::string spf_format;

  // open and close the file to take a a peek at the format.
  SpfFile tmp_spf;
  if (tmp_spf.openSpf(fileName)!=affx::TSV_OK) {
    APT_ERR_ABORT("Unable to open spf file '"+fileName+"'");
  }
  tmp_spf.close();

  if (tmp_spf.m_spf_format==0) {
    APT_ERR_ABORT("Unknown spf file format. (format="+ToStr(tmp_spf.m_spf_format)+")");
  }

  // v3?
  if (tmp_spf.m_spf_format==4) {
    openSpf_v4(fileName,
               probeSetsToLoad,
               probesetNames,
               probeSubSet,
               chipTypeExpected,
               justStats,
               psTypesToLoad);
    return;
  }
  if (tmp_spf.m_spf_format==3) {
    openSpf_v3(fileName,
               probeSetsToLoad,
               probesetNames,
               probeSubSet,
               chipTypeExpected,
               justStats,
               psTypesToLoad);
    return;
  }

  // anything else (none,"1","2") is considered to be "classic"
  openSpf_v2(fileName,
             probeSetsToLoad,
             probesetNames,
             probeSubSet,
             chipTypeExpected,
             justStats,
             psTypesToLoad);
}

void ChipLayout::openSpf_v2(const std::string& fileName,
                            const std::set<const char *, 
                            Util::ltstr> &probeSetsToLoad,
                            std::vector<const char *> *probesetNames,
                            std::vector<bool> &probeSubSet,
                            const std::string &chipTypeExpected,
                            bool justStats,
                            std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad)
{
  string chipType;
  int expectedProbesets = 0, seenProbesets = 0;
  int tmp_X;
  int tmp_Y;

  //
  m_haveKillList = false;

  if (m_SpfFile.openSpf(fileName) != TSV_OK) {
    Err::errAbort("Couldn't open " + ToStr(fileName) + " to read.");
  }
  APT_ERR_ASSERT(m_SpfFile.m_spf_format==2,"missing required v2 columns.");

  // errors if no match
  checkTsvChipType(m_SpfFile,chipTypeExpected);
  // @todo this data is already in "m_SpfFile.m_header_*"
  //
  if (m_SpfFile.getHeader("num-probesets", expectedProbesets) != TSV_OK) {
    Err::errAbort("num-probesets is a required field in file: '" + ToStr(fileName) + "'");
  }
  // @todo This should come from the spf object, not reading the header again.
  m_SpfFile.getHeader("num-channels", m_numChannels);
  //
  if (m_SpfFile.getHeader("num-cols",tmp_X) != TSV_OK) {
    Err::errAbort("num-cols is a required field in file: '" + ToStr(fileName) + "'");
  }
  if (m_SpfFile.getHeader("num-rows", tmp_Y) != TSV_OK) {
    Err::errAbort("num-rows is a required field in file: '" + ToStr(fileName) + "'");
  }
  setDimensions(tmp_X,tmp_Y);
  //
  resizeStats(getProbeCount());
  //
  if (probeSubSet.empty()) {
    probeSubSet.resize(m_XCount * m_YCount, false);
  }

  string name, blockSizesStr, blockAnnsStr, blockAlleleStr, blockContextStr, blockChannelStr, blockRepTypeStr, probesStr;
  int type = -1, numBlocks = -1, numMatch = -1, numProbes = -1;

  m_SpfFile.bind(0, "name", &name, affx::TSV_BIND_REQUIRED);
  m_SpfFile.bind(0, "type", &type, affx::TSV_BIND_REQUIRED);
  m_SpfFile.bind(0, "num_blocks", &numBlocks, affx::TSV_BIND_REQUIRED);
  m_SpfFile.bind(0, "block_sizes", &blockSizesStr, affx::TSV_BIND_REQUIRED);
  m_SpfFile.bind(0, "block_annotations", &blockAnnsStr, affx::TSV_BIND_REQUIRED);
  m_SpfFile.bind(0, "num_match", &numMatch, affx::TSV_BIND_REQUIRED);
  m_SpfFile.bind(0, "num_probes", &numProbes, affx::TSV_BIND_REQUIRED);
  m_SpfFile.bind(0, "probes", &probesStr, affx::TSV_BIND_REQUIRED);

  // optional columns
  if (m_SpfFile.m_v2_block_alleles_cidx>=0) {
    m_SpfFile.bind(0,m_SpfFile.m_v2_block_alleles_cidx,&blockAlleleStr,affx::TSV_BIND_REQUIRED);
  }
  if (m_SpfFile.m_v2_block_contexts_cidx>=0) {
    m_SpfFile.bind(0,m_SpfFile.m_v2_block_contexts_cidx,&blockContextStr,affx::TSV_BIND_REQUIRED);
  }
  if (m_SpfFile.m_v2_block_channels_cidx>=0) {
    m_SpfFile.bind(0,m_SpfFile.m_v2_block_channels_cidx,&blockChannelStr,affx::TSV_BIND_REQUIRED);
  }
  if (m_SpfFile.m_v2_block_rep_types_cidx>=0) {
    m_SpfFile.bind(0,m_SpfFile.m_v2_block_rep_types_cidx,&blockRepTypeStr,affx::TSV_BIND_REQUIRED);
  }

  std::set<const char *, Util::ltstr>::const_iterator emptyIter = probeSetsToLoad.end();
  vector<string> words;
  unsigned int dotMod = max(expectedProbesets/40,1);
  Verbose::progressBegin(1, "Reading " + ToStr(expectedProbesets) + " probesets", 40, dotMod, expectedProbesets);
  while (m_SpfFile.nextLevel(0) == TSV_OK) {
    Verbose::progressStep(1);
    bool keep = true;
    seenProbesets++;
    // pretend like they aren't even there if have specific types to load
    if (!psTypesToLoad.empty() && 
        psTypesToLoad.find((affxcdf::GeneChipProbeSetType)type) == psTypesToLoad.end()) {
      continue;
    }
    if (probesetNames != NULL) {
      probesetNames->push_back(Util::cloneString(name.c_str()));
    }
    bool keepOrder = true;
    if (justStats)
      keep = false;
    if (!probeSetsToLoad.empty() && 
        probeSetsToLoad.find(name.c_str()) == emptyIter) {
      keep = false;
      keepOrder = false;
    }

    parseProbeListRecord(name, type, numBlocks, 
                         blockSizesStr, 
                         blockAnnsStr,
                         (m_SpfFile.m_v2_block_alleles_cidx>=0),
                         blockAlleleStr,
                         (m_SpfFile.m_v2_block_contexts_cidx>=0),
                         blockContextStr,
                         (m_SpfFile.m_v2_block_channels_cidx>=0),
                         blockChannelStr,
                         (m_SpfFile.m_v2_block_rep_types_cidx>=0),
                         blockRepTypeStr,
                         numMatch, numProbes, 
                         probesStr, 
                         m_SpfFile.lineNumber(), keep, probeSubSet, keepOrder);
  }

  //
  m_SpfFile.close();
  Verbose::progressEnd(1, "Done.");
  
  if (seenProbesets != expectedProbesets) {
    Err::errAbort("Expecting: " + ToStr(expectedProbesets) + " probesets in file: '" +
                  ToStr(fileName) + "' got: " + ToStr(seenProbesets));
  }
  if (!justStats) {
    makePsMaps();
  }
}

void ChipLayout::openSpf_v3(const std::string& fileName,
                            const std::set<const char *, 
                            Util::ltstr> &probeSetsToLoad,
                            std::vector<const char *> *probesetNames,
                            std::vector<bool> &probeSubSet,
                            const std::string& chipTypeExpected,
                            bool justStats,
                            std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad) {
  //
  string chipType;
  int xCount = 0, yCount = 0;
  int expectedProbesets = 0;
  //int seenProbesets = 0;
  //bool chipTypeOk = false;
  bool keep;

  //
  // m_PlVec.clear();
  m_haveKillList = false;

  //
  if (m_SpfFile.openSpf(fileName) != TSV_OK) {
    Err::errAbort("Couldn't open " + ToStr(fileName) + " to read.");
  }
  APT_ERR_ASSERT(m_SpfFile.m_spf_format==3,"missing required v3 columns.");

  //if (m_SpfFile.hasHeaderEqualTo("spf-format","3")!=TSV_OK) {
  //Err::errAbort("header 'spf-format'!=3 in '" + ToStr(fileName) + "'");
  //}

  checkTsvChipType(m_SpfFile,chipType);

  //
  if (m_SpfFile.getHeader("num-cols", xCount)!=TSV_OK) {
    Err::errAbort("num-cols is a required field in file: '" + ToStr(fileName) + "'");
  }
  if (m_SpfFile.getHeader("num-rows", yCount)!=TSV_OK) {
    Err::errAbort("num-rows is a required field in file: '" + ToStr(fileName) + "'");
  }
  setDimensions(xCount,yCount);
  //
  if (m_SpfFile.getHeader("num-probesets", expectedProbesets)!=TSV_OK) {
    Err::errAbort("num-probesets is a required field in file: '" + ToStr(fileName) + "'");
  }
  //
  resizeStats(m_XCount * m_YCount);
  if (probeSubSet.empty()) {
    probeSubSet.resize(m_XCount * m_YCount, false);
  }
  
  std::string i_name;
  int i_type;
  int i_annotation;
  int i_allele_code;
  int i_context_code;
  int i_probe_id;

  std::vector<int> block_size_vec;
  std::vector<int> block_annotation_vec;
  std::vector<int> block_allele_vec;
  std::vector<int> block_context_vec;
  std::vector<int> block_channel_vec;
  std::vector<int> block_rep_type_vec;
  std::vector<int> probe_ids_vec;
  int probe_cnt;

  // 0: name type
  unsigned int dotMod = max(expectedProbesets/40,1);
  unsigned int seenProbesets = 0;
  Verbose::progressBegin(1, "Reading " + ToStr(expectedProbesets) + " probesets", 40, dotMod, expectedProbesets);
  while (m_SpfFile.nextLevel(0)==TSV_OK) {
    Verbose::progressStep(1);
    seenProbesets++;
    // clear accum block info...
    block_size_vec.resize(0);
    block_annotation_vec.resize(0);
    block_allele_vec.resize(0);
    block_context_vec.resize(0);
    block_channel_vec.resize(0);
    block_rep_type_vec.resize(0);
    probe_ids_vec.resize(0);

    //
    if ((m_SpfFile.get(0,"name",i_name)!=TSV_OK) || (i_name=="")) {
      Err::errAbort("No name");
    }
    if (m_SpfFile.get(0,"type",i_type)!=TSV_OK) {
      Err::errAbort("No type");
    }
    // pretend like they aren't even there if we have specific types to load.
    if ((!psTypesToLoad.empty()) && 
        (psTypesToLoad.find((affxcdf::GeneChipProbeSetType)i_type) == psTypesToLoad.end())) {
      continue;
    }
    
    // 1: allele allele_code
    while (m_SpfFile.nextLevel(1)==TSV_OK) {
      m_SpfFile.get(1,"allele_code",i_allele_code);
      
      // 2: context context_code
      while (m_SpfFile.nextLevel(2)==TSV_OK) {
        m_SpfFile.get(2,"context_code",i_context_code);
        m_SpfFile.get(2,"annotation",i_annotation);

        // 3: probe_id
        probe_cnt=0;
        while (m_SpfFile.nextLevel(3)==TSV_OK) {
          m_SpfFile.get(3,"probe_id",i_probe_id);
          probe_ids_vec.push_back(i_probe_id);
          probe_cnt++;
        }
        // record block info
        block_allele_vec.push_back(i_allele_code);
        block_annotation_vec.push_back(i_annotation);
        block_context_vec.push_back(i_context_code);
        block_size_vec.push_back(probe_cnt);
        block_channel_vec.push_back(0);
        block_rep_type_vec.push_back(0);
      }
    }

    // like readProbeListRecord, but we have things in vector form already.

    // adjust the probe_ids from 1-based to 0-based
    for (int i=0;i<probe_ids_vec.size();i++) {
      if (probe_ids_vec[i]!=ProbeList::null_probe) {
        probe_ids_vec[i]--;
      }
      // add it to the list
      m_ProbeLayoutOrder.push_back(probe_ids_vec[i]);
    }
    
    //
    if (probesetNames != NULL) {
      probesetNames->push_back(Util::cloneString(i_name.c_str()));
    }
    bool multAlleles = ! ChipLayout::sameValues(block_allele_vec);
    trackGlobalStats(i_name,i_type, 
                     block_size_vec.size(), block_size_vec, 
                     probe_ids_vec.size(), probe_ids_vec, multAlleles);
    //
    keep=true;
    if (justStats) {
      keep=false;
    }
    else if ((!probeSetsToLoad.empty()) && 
             (probeSetsToLoad.find(i_name.c_str()) == probeSetsToLoad.end())) {
      keep=false;
    }

    if (keep) {
      // mark probes to be loaded.
      for (int pIx = 0; pIx < probe_ids_vec.size(); pIx++) {
        if (probe_ids_vec[pIx] != ProbeList::null_probe) {
          probeSubSet[probe_ids_vec[pIx]] = true;
        }
      }
      
      // process the probeset we just read...
      ProbeListPacked pl = makeProbeListPacked(i_name, i_type, 
                                               block_size_vec.size(), 
                                               block_size_vec,
                                               block_annotation_vec,
                                               block_allele_vec,
                                               block_context_vec,
                                               block_channel_vec,
                                               block_rep_type_vec,
                                               1,
                                               ///@todo AW bug? probe_ids_vec.size(), // numMatch,
                                               probe_ids_vec.size(), // numProbes, 
                                               probe_ids_vec);
      //m_PlVec.push_back(pl);
    }
  } // allele

  //
  m_SpfFile.close();
  Verbose::progressEnd(1, "Done.");

  if (seenProbesets != expectedProbesets) {
    Err::errAbort("Expecting: " + ToStr(expectedProbesets) + " probesets in file: '" +
                  ToStr(fileName) + "' got: " + ToStr(seenProbesets));
  }
    
  if (!justStats) {
    makePsMaps();
  }
}

void ChipLayout::openSpf_v4(const std::string& fileName,
                            const std::set<const char *, 
                            Util::ltstr> &probeSetsToLoad,
                            std::vector<const char *> *probesetNames,
                            std::vector<bool> &probeSubSet,
                            const std::string& chipTypeExpected,
                            bool justStats,
                            std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad)
{
  //
  string chipType;
  int xCount = 0, yCount = 0;
  int expectedProbesets = 0;
  //int seenProbesets = 0;
  //bool chipTypeOk = false;
  bool keep;
  //
  // m_PlVec.clear();
  m_haveKillList = false;

  //
  if (m_SpfFile.openSpf(fileName) != TSV_OK) {
    Err::errAbort("Couldn't open " + ToStr(fileName) + " to read.");
  }
  APT_ERR_ASSERT(m_SpfFile.m_spf_format==4,"missing required v4 columns.");
  //if (m_SpfFile.hasHeaderEqualTo("spf-format","4")!=TSV_OK) {
  //Err::errAbort("header 'spf-format'!=4 in '" + ToStr(fileName) + "'");
  //}

  checkTsvChipType(m_SpfFile,chipType);

  //
  if (m_SpfFile.getHeader("num-cols", xCount)!=TSV_OK) {
    Err::errAbort("num-cols is a required field in file: '" + ToStr(fileName) + "'");
  }
  if (m_SpfFile.getHeader("num-rows", yCount)!=TSV_OK) {
    Err::errAbort("num-rows is a required field in file: '" + ToStr(fileName) + "'");
  }
  setDimensions(xCount,yCount);
  //
  // @todo This should come from the spf object, not reading the header again.
  m_SpfFile.getHeader("num-channels", m_numChannels);
  //
  if (m_SpfFile.getHeader("num-probesets", expectedProbesets)!=TSV_OK) {
    Err::errAbort("num-probesets is a required field in file: '" + ToStr(fileName) + "'");
  }
  //
  resizeStats(m_XCount * m_YCount);
  if (probeSubSet.empty()) {
    probeSubSet.resize(m_XCount * m_YCount, false);
  }
  
  std::string i_name;
  int i_type;
  int i_annotation_code;
  int i_allele_code;
  int i_context_code;
  int i_channel_code;
  int i_rep_type;
  int i_probe_id;

  std::vector<int> block_size_vec;
  std::vector<int> block_annotation_vec;
  std::vector<int> block_allele_vec;
  std::vector<int> block_context_vec;
  std::vector<int> block_channel_vec;
  std::vector<int> block_rep_type_vec;
  std::vector<int> probe_ids_vec;
  int block_probe_cnt;

  // 0: name type
  unsigned int dotMod = max(expectedProbesets/40,1);
  unsigned int seenProbesets = 0;
  Verbose::progressBegin(1, "Reading " + ToStr(expectedProbesets) + " probesets", 40, dotMod, expectedProbesets);

  while (m_SpfFile.nextLevel(0)==TSV_OK) {
    Verbose::progressStep(1);
    seenProbesets++;
    // clear block info...
    i_name="";
    i_type=0;
    block_size_vec.resize(0);
    block_annotation_vec.resize(0);
    block_allele_vec.resize(0);
    block_context_vec.resize(0);
    block_channel_vec.resize(0);
    block_rep_type_vec.resize(0);
    probe_ids_vec.resize(0);

    // ray added this to fill out the vector.
    // we dont need to do this anymore as ChipLayout::getProbeSetType()
    // gets it from the factory and modifies the return value.
    // m_PsTypes.push_back(ProbeSet::Marker);
    
    //
    if ((m_SpfFile.get(0,"name",i_name)!=TSV_OK) || (i_name=="")) {
      Err::errAbort("No name");
    }
    if (m_SpfFile.get(0,"type",i_type)!=TSV_OK) {
      Err::errAbort("No type");
    }

    // 1: block annotation allele
    while (m_SpfFile.nextLevel(1)==TSV_OK) {
      // These should be ProbeList::UNSET
      i_annotation_code=0;
      i_allele_code=0;
      i_context_code=0;
      i_channel_code=0;
      i_rep_type=0;
      //
      m_SpfFile.get(1,m_SpfFile.m_v4_annotation_code_cidx,i_annotation_code);
      m_SpfFile.get(1,m_SpfFile.m_v4_allele_code_cidx,i_allele_code);
      m_SpfFile.get(1,m_SpfFile.m_v4_context_code_cidx,i_context_code);
      m_SpfFile.get(1,m_SpfFile.m_v4_channel_code_cidx,i_channel_code);
      m_SpfFile.get(1,m_SpfFile.m_v4_rep_type_cidx,i_rep_type);
      
      // 2: probe_id
      block_probe_cnt=0;
      while (m_SpfFile.nextLevel(2)==TSV_OK) {
        i_probe_id=-1;
        m_SpfFile.get(2,m_SpfFile.m_v4_probe_id_cidx,i_probe_id);
        probe_ids_vec.push_back(i_probe_id);
        block_probe_cnt++;
      }

      // record block info
      block_size_vec.push_back(block_probe_cnt);
      block_annotation_vec.push_back(i_annotation_code);
      block_allele_vec.push_back(i_allele_code);
      block_context_vec.push_back(i_context_code);
      block_channel_vec.push_back(i_channel_code);
      block_rep_type_vec.push_back(i_rep_type);
    } // while blocks

    // adjust the probe_ids from 1-based to 0-based
    for (int i=0;i<probe_ids_vec.size();i++) {
      if (probe_ids_vec[i]!=ProbeList::null_probe) {
        probe_ids_vec[i]--;
      }
      // add it to the list
      m_ProbeLayoutOrder.push_back(probe_ids_vec[i]);
    }

    // pretend like they aren't even there if we have specific types to load.
    if ((!psTypesToLoad.empty()) && 
        (psTypesToLoad.find((affxcdf::GeneChipProbeSetType)i_type) == psTypesToLoad.end())) {
      continue;
    }

    //
    bool multAlleles = ! ChipLayout::sameValues(block_allele_vec);
    trackGlobalStats(i_name,i_type, 
                     block_size_vec.size(), 
                     block_size_vec,
                     probe_ids_vec.size(), 
                     probe_ids_vec,
                     multAlleles);
    //
    if (probesetNames != NULL) {
      probesetNames->push_back(Util::cloneString(i_name.c_str()));
    }
    
    //
    keep=true;
    if (justStats) {
      keep=false;
    }
    else if ((!probeSetsToLoad.empty()) && 
             (probeSetsToLoad.find(i_name.c_str()) == probeSetsToLoad.end())) {
      keep=false;
    }
    
    if (keep) {
      // mark probes to be loaded.
      for (int pIx = 0; pIx < probe_ids_vec.size(); pIx++) {
        if (probe_ids_vec[pIx] != ProbeList::null_probe) {
          probeSubSet[probe_ids_vec[pIx]] = true;
        }
      }
      
      // process the probeset we just read...
      ProbeListPacked pl = makeProbeListPacked(i_name, i_type, 
                                               block_size_vec.size(), 
                                               block_size_vec,
                                               block_annotation_vec,
                                               block_allele_vec,
                                               block_context_vec,
                                               block_channel_vec,
                                               block_rep_type_vec,
                                               1,
                                               ///@todo AW bug? probe_ids_vec.size(), // numMatch,
                                               probe_ids_vec.size(), // numProbes, 
                                               probe_ids_vec);
    }
  } // while probesets

  //
  m_SpfFile.close();
  Verbose::progressEnd(1, "Done.");

  if (seenProbesets != expectedProbesets) {
    Err::errAbort("Expecting: " + ToStr(expectedProbesets) + " probesets in file: '" +
                  ToStr(fileName) + "' got: " + ToStr(seenProbesets));
  }
  
  if (!justStats) {
    makePsMaps();
  }
}

ProbeListPacked ChipLayout::getProbeListAtIndex(unsigned int index) {
  return m_PlFactory.getProbeListAtIndex(index);
}
//  if ((index>=0)&&(index<=getProbeSetCount())) {
//    return m_PlVec[index];
//  }
//  //
//  return ProbeListPacked(NULL);
// }

/**
 * Find a probe set using the name as a key.
 * @param psName - Name of probe set desired.
 * @return Reference to probe set, NULL if not found.
 */
ProbeListPacked ChipLayout::getProbeListByName(const std::string& name) {
  return m_PlFactory.getProbeListByName(name);
}
//  psNameIter i;
//  if (m_PsNameHash.empty())
//    Err::errAbort("ChipLayout::getProbeListByName() - Name lookup table is empty.");
//  i = m_PsNameHash.find(name.c_str());
//  if (i == m_PsNameHash.end()) {
//    if (!m_haveKillList)
//        Err::errAbort("Can't find probe set with name: " + string(name));
//    else
//      ProbeListPacked(NULL);
//  }
//  return m_PlVec[i->second];
//}

/**
 * Find a probe set index using the name as a key.
 * @param psName - Name of probe set desired.
 * @return probe set index
 */
int ChipLayout::getProbeSetIndexByName(const std::string& name) {
  return m_PlFactory.getProbeListIndexByName(name);
  //  psNameIter i;
  //  if (m_PsNameHash.empty())
  //    Err::errAbort("ChipLayout::getProbeSetIndexByName() - Name lookup table is empty.");
  //  i = m_PsNameHash.find(name.c_str());
  //  if (i == m_PsNameHash.end()) {
  //    if (!m_haveKillList)
  //        Err::errAbort("Can't find probe set with name: " + string(name));
  //    else
  //        return -1;
  //  }
  //  return i->second;
}

/**
 * Make maps of the probe sets in m_PsVec by id and name if possible.
 */
void ChipLayout::makePsMaps() {
  //  /* We dont need to do this any more */
  //  unsigned int i = 0;
  //  m_HaveNameMap = true;
  //  /* Loop through and hash the probe sets by id and name. */
  //  for (i = 0; i < getProbeSetCount(); i++) {
  //    ProbeListPacked pl = m_PlVec[i];
  //    /* if they are duplicate probe sets. */
  //    if (m_PsNameHash.find(pl.get_name_cstr()) != m_PsNameHash.end())
  //      Err::errAbort("Duplicate probe set: '" + pl.get_name_string()
  //                    + "' is already in name hash.");
  //    m_PsNameHash[pl.get_name_cstr()] = i;
  //  }
  //
}

/**
 * Open and parse out the data of a file.
 * @param fileName - Name of PGF file to be opened.
 * @param xSize - How many features in x dimension.
 * @param ySize - How many features in y dimension.
 * @param probeSetsToLoad - Vector of probe set ids to be loaded.
 * @param probeSubset - Subset of probes to be loaded.
 * @param chipType - What sort of chip should this pgf file be for?
 * @return bool - true if successful.
 */
bool ChipLayout::openPgf(const std::string& fileName,
                         unsigned int xSize, unsigned int ySize,
                         const std::set<const char *, 
                         Util::ltstr> &probeSetsToLoad,
                         std::vector<const char *> *probesetNames,
                         std::vector<const char *> *probesetDisplayNames,
                         std::vector<bool> &probeSubset,
                         const std::string &chipType,
                         probeidmap_t &killList,
                         bool justStats, 
                         bool doPList)
{
  // this was a member var (m_Pgf), now it is just a var.
  affx::PgfFile pgf;
  pgf.close();

  int rc = pgf.open(fileName);
  if (rc != TSV_OK) {
    Err::errAbort("ChipLayout::openPgf() - Couldn't open pgf file: '" + fileName + "'");
  }

  if (!killList.empty())
    m_haveKillList = true;

  // Check that we match one of the chip_types in the file.
  /// @todo  move this into PgfFile.h
  bool chipTypeOk = false;
  std::string tmpChipType;

  pgf.m_tsv.headersBegin();
  if (chipType!="") {
    while ((!chipTypeOk)&&
           (pgf.m_tsv.headersFindNext("chip_type",tmpChipType)==TSV_OK)) {
      if (tmpChipType==chipType) {
        chipTypeOk=true;
      }
    }
    if (!chipTypeOk) {
      Err::errAbort("In file '" + ToStr(fileName) + "' does not support chip type: " + chipType);
      return false;
    }
  }

  // Copy the header values to the layout.
  // it will use them while printing.
  std::string key;
  std::string val;
  pgf.m_tsv.headersBegin();
  while (pgf.m_tsv.headersNext(key,val)==TSV_OK) {
    if (key == "probeset_count") {
      m_ProbeCounts.reserve(Convert::toInt(val.c_str()));
    }
    m_Header[key].push_back(val);
  }
  setDimensions(xSize,ySize);
  // do the actual data suck
  unsigned int numProbes = m_XCount * m_YCount;
  //  if (doPList)
  if (m_SpfFileName != "") {
    openSpfForWrite(m_SpfFile,m_SpfFileName);
  }
  readProbeListPgfFileKillList(pgf, numProbes, m_PlFactory, probeSetsToLoad,
                               probesetNames, probesetDisplayNames, probeSubset, killList, justStats);
  pgf.close();
  if (!justStats)
    makePsMaps();
  Verbose::out(2, "Loaded " + ToStr(getProbeSetCount()) + " probe sets.");
  if (m_SpfFile.is_open())
    m_SpfFile.close();
  return true;
}

/**
 * Read the contents of rowFile into a vector of
 * probe sets.
 * @param pgf - Opened Pgf file to read from
 * @param numProbes - How many probes are on this chip?
 * @param psVec - vector to be filled in.
 * @param probeSetsToLoad - Vector of probe set ids to be loaded.
 * @param probeSubset - Subset of probes to be loaded.
 * @param killList - List of probes not to be used.
 */
void ChipLayout::readProbeListPgfFileKillList(affx::PgfFile &pgf,
                                              unsigned int numProbes,
                                              //std::vector<ProbeListPacked> &plVec,
                                              ProbeListFactory &plFactory,
                                              const std::set<const char *, 
                                              Util::ltstr> &probeSetsToLoad,
                                              std::vector<const char *> *probesetNames,
                                              std::vector<const char *> *probesetDisplayNames,
                                              std::vector<bool> &probeSubset,
                                              probeidmap_t &killList,
                                              bool justStats) {
  bool probeSetsToLoad_empty=probeSetsToLoad.empty();
  set<const char *, Util::ltstr>::const_iterator probeSetsToLoad_iter;
  bool keepProbeSet;

  // If probeset count is set, us it to display progress info
  int expectedProbesets = -1;
  pgf.m_tsv.getHeader("probesets",expectedProbesets);

  string val;
  if (pgf.m_tsv.headersFindNext("probesets",val)==TSV_OK){
    expectedProbesets = Convert::toInt(val.c_str());
  }

  // skip to the start of the data
  pgf.rewind();
  // m_PmProbes is a bitmask of which probes are PM.  (pm=true)
  resizeStats(numProbes);

  // Set up some local static storage. Once a probeset is determined
  // to be kept we will allocate memory on the heap for it. Trying to
  // prevent lots of new and deletes as they can fragement memory.
  ProbeSet psConst;
  int psConstMaxName = 256;
  psConst.name = new char[psConstMaxName];
  psConst.name[psConstMaxName -1] = '\0'; // null terminate.
  vector<Atom> *locAtoms = new vector<Atom>(MAX_ATOMS);
  vector<Probe> *locProbes = new vector<Probe>(MAX_PROBES);
  // gobble up probesets
  string tempName;
  unsigned int dotMod = max(expectedProbesets/40,1);
  if (expectedProbesets > 0)
    Verbose::progressBegin(1, "Reading " + ToStr(expectedProbesets) + " probesets", 40, dotMod, expectedProbesets);
  unsigned int seenProbesets = 0;
  while (pgf.next_probeset()==TSV_OK) {
    if (expectedProbesets > 0)
      Verbose::progressStep(1);
    seenProbesets++;
    //clearProbeSet(psConst); // redundant.
    int pmProbeCount = 0;
    int pmProbePossible = 0;
    clearProbeSet(psConst);

    // it is always this (for now I guess -jhg)
    psConst.psType = ProbeSet::Expression;

    // copy the data over
	tempName = ToStr(pgf.probeset_id);
    if (tempName.size() + 1 >= psConstMaxName) {
      delete[] psConst.name;
      psConst.name = new char[tempName.size() + 1];
      psConstMaxName = tempName.size() + 1;
      psConst.name[psConstMaxName -1] = '\0'; // make sure always null terminated
    }
    strncpy(psConst.name, tempName.c_str(), psConstMaxName -1);

    // do we want to keep this probeset?
    keepProbeSet=true; // yes unless...
    if (!probeSetsToLoad_empty) {
      probeSetsToLoad_iter=probeSetsToLoad.find(psConst.name);
      if (probeSetsToLoad_iter==probeSetsToLoad.end()) {
        keepProbeSet=false;
      }
    }
    set<probeid_t> probesSeen;
    set<atomid_t> atomsSeen;
    set<probeid_t>::iterator seenIter;

    int locAtomIndex = 0;
    int locProbeIndex = 0;
    while (pgf.next_atom()==TSV_OK) {
      bool lastPmDropped = false;
      int atomProbeIndex = locProbeIndex;
      if (locAtomIndex + 1 > locAtoms->size()) {
        Verbose::out(2, "Resizing atom size from: " + ToStr(locAtoms->size()) + " to " + ToStr(locAtoms->size() + MAX_ATOMS));
        resizeAtoms(&locAtoms, locAtoms->size() + MAX_ATOMS, &psConst);
      }
      Atom *atom = &(*locAtoms)[locAtomIndex++];
      // copy over data
      atom->id = pgf.atom_id;
      // check to make sure we haven't seen this atom in this probeset before...
      if (atomsSeen.find(atom->id) != atomsSeen.end())
        Err::errAbort("Duplicate atom ids (" + ToStr(atom->id) + ") in probeset : " + ToStr(psConst.name));
      else
        atomsSeen.insert(atom->id);
      while (pgf.next_probe()==TSV_OK) {
        probeidmap_t::iterator iter;

        // If current probes are greater than local storage have to increase.
        if (locProbeIndex + 1 > locProbes->size()) {
          Verbose::out(2, "Resizing probe pm/mm size from: " + ToStr(locProbes->size()) + " to " +
                       ToStr(locProbes->size() + MAX_PROBES));
          resizeProbes(&locProbes, locProbes->size() + MAX_PROBES, locAtoms, locAtomIndex);
        }
        Probe *probe = &(*locProbes)[locProbeIndex++];
        // copy over data
        probe->id=pgf.probe_id;
        // check to make sure we haven't seen this probe in this probeset before...
        if (probesSeen.find(probe->id) != probesSeen.end())
          Err::errAbort("Duplicate probe ids (" + ToStr(probe->id) + ") in probeset : " + ToStr(psConst.name));
        else
          probesSeen.insert(probe->id);
        probe->type=Probe::typeForString(pgf.probe_type.c_str());
        if (pgf.gc_count >= NULLPROBEGC) {
          Verbose::out(1,"GC count, " + ToStr(pgf.gc_count) + ", exceeds software limit of " + ToStr(NULLPROBEGC - 1) + ". Ignoring.");
          probe->gcCount = (char)NULLPROBEGC;
        } else if (pgf.gc_count == -1) {
          probe->gcCount = (char)NULLPROBEGC;
        } else {
          probe->gcCount = (char)pgf.gc_count;
        }
        probe->id--; // Alan created pgf files 1-based, we want zero based.

        if (Probe::isPm(*probe))
          pmProbePossible++;

        if (inKillList((probeid_t)probe->id+1, psConst.name, killList)) {
          Verbose::out(4,"Ignoring probe " + ToStr(probe->id+1) + " in probeset " + ToStr(psConst.name) + " due to kill list");
          if (Probe::isPm(*probe)) {
            lastPmDropped = true;
          }
          clearProbe(*probe);
          locProbeIndex--;
        }
        // Did we exclude this probe's matching PM last time?
        else if (lastPmDropped && Probe::isMm(*probe)) {
          Verbose::out(4,"Ignoring MM probe " + ToStr(probe->id+1) + " in probeset " + ToStr(psConst.name) + " due to PM in kill list");
          clearProbe(*probe);
          locProbeIndex--;
          lastPmDropped = false;
        }
        // treat it like a normal probe.
        else {
          // if this probeset contains a probe to be kept
          // make sure to keep the probeset around
          if (probe->type == Probe::PMAT || probe->type == Probe::PMST) {
            pmProbeCount++;
          }
          //   atom->probes.push_back(probe);
          lastPmDropped = false;
        }
      } // probe
      atom->probes.resize(locProbeIndex - atomProbeIndex);
      for (uint32_t pIx = 0; pIx + atomProbeIndex < locProbeIndex; pIx++) {
        atom->probes[pIx] = &(*locProbes)[pIx + atomProbeIndex];
      }
    } // atom
    psConst.atoms.reserve(locAtomIndex);
    for (uint32_t aIx = 0; aIx < locAtomIndex; aIx++) {
      if ((*locAtoms)[aIx].probes.size() > 0) {
        psConst.atoms.push_back(&(*locAtoms)[aIx]);
      }
    }
    psConst.numGroups = 1; // pgf files only used for expression probesets.
    psConst.atomsPerGroup.resize(1);
    psConst.atomsPerGroup[0] = psConst.atoms.size();
    // if we excluded all the PM probes don't include this probeset.
    // we have some mm-only probesets in pgf lib files hence the check
    // that the probeset had some pm probes to start with.
    if (pmProbeCount == 0 && pmProbePossible > 0)
      keepProbeSet = false;

    int probeCount = psConst.countProbes();

    m_ProbeCounts.push_back(probeCount);
    // Add to chip level pm, mm, gc counts.
    int numBytes = ProbeListPacked::byte_size(psConst.numGroups, probeCount, strlen(psConst.name));
    trackGlobalStats(psConst, numBytes);

    //std::cout << "Probeset: " << psConst << '\n';

    if (keepProbeSet) {
      for (int aIx = 0; aIx < psConst.atoms.size(); aIx++) {
        for (int pIx = 0; pIx < psConst.atoms[aIx]->probes.size(); pIx++) {
          m_ProbeLayoutOrder.push_back(psConst.atoms[aIx]->probes[pIx]->id);
        }
      }
    }

    // Keep it around for processing?
    if (psConst.atoms.empty()) {
      Verbose::warn(1, "Empty atoms vector for probeset '" + ToStr(psConst.name) +
                    "', probeset name '" + psConst.name + "'");
    } else {
      if (m_SpfFile.is_open() && keepProbeSet) {
        writeSpfProbeSetRecord(m_SpfFile, psConst);
      }
      else if (keepProbeSet && !justStats) {
        // ProbeListPacked pList = plFactory.pListFromPSet(psConst);
        // plVec.push_back(pList); // store the pointer
        plFactory.pListFromPSet(psConst);
      }

      // Keep track of ps type
      // m_PsTypes.push_back((char)psConst.psType);

      if (probesetNames != NULL) {
        probesetNames->push_back(Util::cloneString(psConst.name));
      }
	  if (probesetDisplayNames != NULL) {
            probesetDisplayNames->push_back((pgf.probeset_name.empty() == false ? Util::cloneString(pgf.probeset_name.c_str()) : NULL));
	  }
    }

    clearProbeSet(psConst);
  } // probeset
  if (expectedProbesets > 0) {
    Verbose::progressEnd(1, "Done.");
    if (expectedProbesets > 0 && (seenProbesets != expectedProbesets)) {
      Err::errAbort("Expecting: " + ToStr(expectedProbesets) + " probesets,  got: " + ToStr(seenProbesets));
    }
  }
  delete locProbes;
  for (uint32_t aIx = 0; aIx < locAtoms->size(); aIx++) {
    (*locAtoms)[aIx].probes.clear();
  }
  delete locAtoms;
  psConst.atomsPerGroup.clear();
  psConst.atoms.clear();
}

/**
 * Open a file and parse the probe sets.
 * @param fileName - Name of file to be opened.
 * @param probeSetsToLoad - Which probe sets should be loaded?
 * @param probeSubset - Subset of probes to be loaded.
 * @param chipType - What sort of chip should this cdf file be for?
 * @param probesetNames - If not null filled in with the
 * probesetname for each probeset even if not actually into
 * chiplayout. Useful for chp files where every single probeset must
 * be included and in order.
 * @param - killList map of probes that are to be excluded from
 * the probesets. A probe must have its id mapped to the probeset name
 * for this to work.
 * @param justStats - just read and generate stats, don't load into memory.
 * @param psTypesToLoad - What types of probesets to load into memory.
 * @return bool - true if success. Errors out on problems.
 */
bool ChipLayout::openCdf(const std::string& fileName,
                         const std::set<const char *, 
                         Util::ltstr> &probeSetsToLoad,
                         std::vector<const char *> *probesetNames,
                         std::vector<bool> &probeSubset,
                         const std::string &chipType,
                         probeidmap_t &killList,
                         bool justStats,
                         std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad) {
  FusionCDFData *cdf = new FusionCDFData();
  FusionCDFData *cdfFileHeader = new FusionCDFData();
  //m_HaveNames = true;

  std::string cdf_unc_name=Fs::convertToUncPath(fileName);
  cdfFileHeader->SetFileName(cdf_unc_name.c_str());

  if (!killList.empty())
    m_haveKillList = true;

  if (!(cdfFileHeader->Exists())) {
    delete cdf;
    delete cdfFileHeader;
    Err::errAbort("Can't open CDF file " + FS_QUOTE_PATH(cdf_unc_name) + " to read.");
  }
  if (!cdfFileHeader->ReadHeader()) { // All comes in with one big gulp...
    string msg = "Error opening CDF file to read header "+FS_QUOTE_PATH(cdf_unc_name)
      +": " + cdf->GetError();
    delete cdf;
    delete cdfFileHeader;
    Err::errAbort(msg);
  }
  if (chipType != "") {
    vector<std::string> chipTypes = cdfFileHeader->GetChipTypes();
    bool match = false;
    for (int i = 0; i < chipTypes.size(); i++)
      if (chipTypes[i] == chipType) {
        match = 1;
        break;
      }
    if (!match) {
      string msg = "Cdf file: " + fileName + " is of type " + chipTypes[0] + " but chip type: " +
        chipType + " was requested.";
      delete cdf;
      delete cdfFileHeader;
      Err::errAbort(msg);
    }
  }

  // 
  m_Header["chip_type"].push_back(cdfFileHeader->GetChipType());
  m_Header["lib_set_name"].push_back(cdfFileHeader->GetFileName());
  m_Header["lib_set_version"].push_back( cdfFileHeader->GetFileName());
  m_Header["cdf-guid"].push_back( cdfFileHeader->GetGUID() );

  FusionCDFFileHeader &cdfHeader = cdfFileHeader->GetHeader();
  /* there is nothing here...*/
  Verbose::out(2, "There are: " + ToStr(cdfHeader.GetNumProbeSets()) + " probesets.");
  string numProbesets = ToStr(cdfHeader.GetNumProbeSets());
  m_Header["cdf_total_probesets"].push_back(numProbesets);
  setDimensions(cdfHeader.GetCols(), cdfHeader.GetRows());
  if (probeSubset.empty()) {
    probeSubset.resize(cdfHeader.GetCols() * cdfHeader.GetRows(), false);
  }
  cdfFileHeader->Close();
  delete cdfFileHeader;
  
  if (m_SpfFileName != "") {
    openSpfForWrite(m_SpfFile,m_SpfFileName);
  }

  cdf->SetFileName(cdf_unc_name.c_str());
  if (!cdf->Read()) { // All comes in with one big gulp...
    string msg = "Error opening CDF file: "+FS_QUOTE_PATH(cdf_unc_name)+": " + cdf->GetError();
    delete cdf;
    Err::errAbort(msg);
  }

  /* Convert cdf to probe sets. */
  readProbeListCdfFileKillList( *cdf, 
                                m_PlFactory, 
                                probeSetsToLoad,
                                probesetNames, 
                                probeSubset, 
                                killList, 
                                justStats,
                                psTypesToLoad);
  /* Make our lookup maps. */
  if (!justStats) {
    makePsMaps();
  }
  Verbose::out(2, "Loaded " + ToStr(getProbeSetCount()) + " probe sets.");
  if (m_SpfFile.is_open())
    m_SpfFile.close();

  delete cdf;
  return true;
}

/**
 * Write out the contents of probe sets in PGF format.
 * @param fileName - File to write data to.
 * @param vec - Vector of probe sets to be stored.
 */
void ChipLayout::writePgfFile(const std::string& fileName, 
                              const std::vector<ProbeSet *> &vec) {
  affx::PgfFile pgf;

  pgf.defineFilePgf();

  pgf.write(fileName);
  writePgfFile(pgf,vec);
  pgf.close();
}

/**
 * Write the vector of probesets to a PGF file
 * @param pgf - Pgf to write to
 * @param vec -
 */
/// @todo shouldnt this be called "writePgfProbeSetVector"?
void ChipLayout::writePgfFile(affx::PgfFile& pgf, 
                              const std::vector<ProbeSet *> &vec) {

  ProbeSet* probeset;
  Atom*     atom;
  Probe*    probe;

  for (unsigned int probeset_i = 0;probeset_i < vec.size(); probeset_i++) {
    probeset=vec[probeset_i];
    // write a probeset
    pgf.m_tsv.set(0,"probeset_id"  ,ToStr(probeset->name));
    ///@todo type here is not right - should not be probeset ID
    //pgf.m_tsv.set(0,"type"  ,(int)probeset->id);
    //pgf.m_tsv.set(0,"probeset_name",probeset->name);
    pgf.m_tsv.writeLevel(0);

    for (unsigned int atom_i = 0;atom_i < probeset->atoms.size(); atom_i++) {
      atom=probeset->atoms[atom_i];
      // write an atom
      pgf.m_tsv.set(1,"atom_id",(int)atom->id);
      pgf.m_tsv.writeLevel(1);

      for (unsigned int probe_i = 0;probe_i < atom->probes.size(); probe_i++) {
        probe=atom->probes[probe_i];
        // write a probe
        pgf.m_tsv.set(2,"probe_id"              ,(int)probe->id);
        pgf.m_tsv.set(2,"type"                  ,Probe::stringForType((Probe::ProbeType)probe->type));
        pgf.m_tsv.set(2,"gc_count"              ,(int)probe->gcCount);
        /// @todo  these were not defined in ProbeSet.h so we cant write them out.
        //pgf.m_tsv.set(2,"probe_length"          ,25);
        //pgf.m_tsv.set(2,"interrogation_position",13);
        //pgf.m_tsv.set(2,"probe_sequence"        ,"ACAACGACCGTTCCGGAATCGACAT");
        //pgf.m_tsv.set(2,"exon_position"         ,1703);
        pgf.m_tsv.writeLevel(2);
      }
    }
  }
}

/**
 * Turn the x and y cordinates of a probe on the array into an index. Should
 * be exactly the same index as provided by cel files.
 *
 * @param x - x coordinate on array.
 * @param y - y coordinate on array.
 * @return - index (id) of feature in cel file and other arrays.
 */
int ChipLayout::xyToIndex(int x, int y) {
  return FusionCELData::XYToIndex(x,y,m_YCount,m_XCount);
}


/**
 * Convert from numeric representation to string.
 * @param type - numeric representation.
 */
const char *ChipLayout::typeFromCdfType(unsigned short type) {
  switch(type) {
  case affxcdf::UnknownProbeSetType :
    return "unknown";
  case affxcdf::ExpressionProbeSetType :
    return "expression";
  case affxcdf::GenotypingProbeSetType :
    return "genotype";
  case affxcdf::ResequencingProbeSetType :
    return "reseq";
  case affxcdf::TagProbeSetType :
    return "tag";
  case affxcdf::CopyNumberProbeSetType :
    return "copynumber";
  default :
    Err::errAbort("Don't recognize probe set type: " + ToStr(type));
  }
  return NULL;
}

/**
 * Is this probe a perfect match probe?
 * @param cel - probe information from cdf.
 * @return true if probe is a perfect match, false otherwise.
 */
bool ChipLayout::isPm(const FusionCDFProbeInformation &cel) {
  char pBase = 'x', tBase = 'x';
  pBase = tolower(cel.GetPBase());
  tBase = tolower(cel.GetTBase());
  if ((pBase == 'a' && tBase == 't') ||
     (pBase == 't' && tBase == 'a') ||
     (pBase == 'c' && tBase == 'g') ||
     (pBase == 'g' && tBase == 'c')) {
    return true;
  }
  return false;
}

/**
 * Check to see if a charater is a valid DNA base of 'a','t','g', or 'c'
 * @param c - base to check.
 * @return - return true if c is in 'atgcATGC', false otherwise.
 */
bool ChipLayout::validBase(int c) {
  c =  tolower(c);
  if (c == 'a' || c == 't' || c == 'g' || c == 'c')
    return true;
  return false;
}

/**
 * Read cdf probe information into a Probe object.
 *
 * @param p - probe object to be filled in.
 * @param parent - Atom that is parent of this probe.
 * @param numCol - column on chip.
 * @param numRow - row on chip.
 * @param direction - target direction.
 * @param group - group information from cdf.
 * @param cel -probe information from cdf.
 */
void ChipLayout::readCdfProbe(  Probe *p, 
                                Atom *parent, 
                                unsigned int numCol, 
                                unsigned int numRow,
                                unsigned int direction,
                                const FusionCDFProbeGroupInformation &group,
                                const FusionCDFProbeInformation &cel) {
  int x = 0, y = 0;
  char pBase = 'x', tBase = 'x';
  // Do the first probe, often mismatch.
  x = cel.GetX();
  y = cel.GetY();

  pBase = cel.GetPBase();
  tBase = cel.GetTBase();
  p->gcCount = (char)NULLPROBEGC;
  p->id = xyToIndex(x,y);
  /* When pBase is same as tBase it is a mismatch probe
     as not complementary. */
  /* If target is complementary to probe it is a pm probe. */
  if (isPm(cel)) {
    if (direction == affxcdf::SenseDirection)
      p->type = Probe::PMST;
    else if (direction == affxcdf::AntiSenseDirection)
      p->type = Probe::PMAT;
    else
      p->type = Probe::PMAT;
    ///@todo there should be a Probe::PM type where direction is not known
  }
  /* mismatches choise can be constrained by cross hybe on snps. */
  else if (validBase(pBase) && validBase(tBase)) {
    if (direction == affxcdf::SenseDirection)
      p->type = Probe::MMST;
    else if (direction == affxcdf::AntiSenseDirection)
      p->type = Probe::MMAT;
    else
      p->type = Probe::MMAT;
    ///@todo there should be a Probe::PM type where direction is not known
  }
  else {
    Verbose::out(1, "Don't recognize probe type with tbase: '" + ToStr(tBase) + "'"
                 " and pbase: '" + ToStr(pBase) + "' at x,y: " + ToStr(x) + "," + ToStr(y));
  }

}

/**
 * Is a particular probe in the kill list (i.e. not to be used?)
 * @param probeIx - Id of probe.
 * @param psName - Name of probeset that contains probe.
 * @param killList - Map of probes not to be used.
 * @return true if probe id is in kill list, false otherwise.
 */
bool ChipLayout::inKillList(probeid_t probeId, 
                            const std::string& psName, 
                            probeidmap_t &killList) {
  /* simplest case. */
  if (killList.empty())
    return false;
  /* see if this probe is in the kill list. */
  probeidmap_t::iterator iter;
  iter = killList.find(probeId);
  if (iter == killList.end())
    return false;
  /* check probeset */
  std::vector<std::string>::iterator psIter = find(killList[probeId].begin(),killList[probeId].end(),psName);
  if (psIter == killList[probeId].end())
    return false;
  return true;
}




/**
 * Read a vector of ProbeSets from a cdf file representation.
 *
 * @param cdf - cdf file object with probe set information.
 * @param psVec - vector to put newly created probe sets in.
 * @param probeSetsToLoad - which probe sets should be converted.
 * @param probeSubSet - which individual probes were loaded.
 * @param probesetNames - If not null filled in with the
 * probesetname for each probeset even if not actually into
 * chiplayout. Useful for chp files where every single probeset must
 * be included and in order.
 * @param killList - map of probes ids that shouldn't be used. Map should
 * be of probe ids to probeset names. If a probe is in the map it will
 * be excluded as well as matching MM probe (or PM if probe is MM). We check
 * to make sure that the string that the user is mapping the probe to the same
 * probeset as the cdf does.
 */

#define READCDFPROBESET_ERR_MSG "ERROR: readCdfProbeSet(): "

/// @brief     
/// @param     name      
/// @param     set       
/// @param     psIx      
/// @param     numCol    
/// @param     numRow    
/// @param     probeSubSet
/// @param     keepProbeSet
/// @param     killList  
/// @param     ps        
/// @param     psNameSize
/// @param     locAtoms  
/// @param     locProbes 
/// @param     errorCount   Count of errors so far.
void ChipLayout::readCdfProbeSet(const std::string& name,
                                 FusionCDFProbeSetInformation &set,
                                 unsigned int psIx, 
                                 int numCol, int numRow,
                                 std::vector<bool> &probeSubSet, 
                                 bool &keepProbeSet,
                                 probeidmap_t &killList,
                                 ProbeSet &ps, 
                                 int &psNameSize, 
                                 std::vector<Atom> **locAtoms,
                                 std::vector<Probe> **locProbes,
                                 int& errorCount)
{
  // *Warning* This function has evolved into quite a behemoth with a
  // series of feature requests and no time to refactor it, be
  // careful with modifications.
  FusionCDFProbeGroupInformation group;
  FusionCDFProbeInformation cel1, cel2;
  int locAtomIndex = 0, locProbeIndex = 0;
  clearProbeSet(ps);
  int name_size=name.size();
  if (name_size + 1 >= psNameSize) {
    delete [] ps.name;
    ps.name = new char[name_size + 1];
    psNameSize = name_size + 1;
    ps.name[psNameSize -1] = '\0'; // make sure always null terminated
  }
  strncpy(ps.name, name.c_str(), psNameSize -1);
  affxcdf::GeneChipProbeSetType set_probesettype=set.GetProbeSetType();
  if (set_probesettype == affxcdf::ExpressionProbeSetType)
    ps.psType = ProbeSet::Expression;
  else if (set_probesettype == affxcdf::GenotypingProbeSetType)
    ps.psType = ProbeSet::GenoType;
  else if (set_probesettype == affxcdf::CopyNumberProbeSetType)
    ps.psType = ProbeSet::Copynumber;
  else if (set_probesettype==affxcdf::MarkerProbeSetType)
    ps.psType = ProbeSet::Marker;
  else if (set_probesettype==affxcdf::MultichannelMarkerProbeSetType)
    ps.psType = ProbeSet::MultichannelMarker;
  else
    ps.psType = ProbeSet::Unknown;
  int numGroups = set.GetNumGroups();
  if (ps.psType == ProbeSet::GenoType && numGroups != 2 && numGroups != 4) {
    // there is a bug in some cdf files where controls are marked as genotyping...
    Verbose::out(4, "How can probeset: " + ToStr(ps.name) + " not have 2 or 4 groups and still be genotyping?");
  }
  ps.atomsPerGroup.resize(numGroups);
  int totalAtoms = 0;
  for (int groupIx = 0; groupIx < numGroups; groupIx++) {
    set.GetGroupInformation(groupIx, group);
    totalAtoms += group.GetNumLists();
  }
  if (locAtomIndex + totalAtoms >= (*locAtoms)->size()) {
    Verbose::out(2, "Resizing atom size from: " + ToStr((*locAtoms)->size()) +
                 " to " + ToStr((*locAtoms)->size() + totalAtoms + MAX_ATOMS));
    resizeAtoms(locAtoms, (*locAtoms)->size() + totalAtoms + MAX_ATOMS, &ps);
  }
  ps.atoms.reserve(totalAtoms);
  // For each group (block?) read all the atoms. 
  for (int groupIx = 0; groupIx < numGroups; groupIx++) {
    bool expectingMisMatches = false;
    set.GetGroupInformation(groupIx, group);
    unsigned short direction = group.GetDirection();
    int numAtoms = group.GetNumLists();
    int goodAtoms = 0;

    // Have to be two probes per list (atom) to have a pm & mm probe.
    // Don't trust the GetNumCellsPerList() function, sometimes lies...
    if (set.GetNumLists() != 0 &&
       set.GetNumCells() / set.GetNumLists() == 2)
      expectingMisMatches = true;
    // Read in each atom.
    int cellCount = 0; // Keep track of how many cells from this group we've seen already.
    for (int atomIx = 0; atomIx < numAtoms; atomIx++) {
      //
      Atom *a = &(**locAtoms)[locAtomIndex++];
      a->id = atomIx;

      if (ps.psType == ProbeSet::Marker) {
        a->allele_code=group.GetAlleleCode();
        a->context_code=group.GetWobbleSituation();
      }
      else if (ps.psType == ProbeSet::GenoType) {
        if ( groupIx % 2 == 0)
          a->allele_code=0;
        else
          a->allele_code=1;
      } 
      else if (ps.psType == ProbeSet::MultichannelMarker) {
        a->allele_code=group.GetAlleleCode();
        a->context_code=group.GetWobbleSituation();
        a->channel_code=group.GetChannel();
        m_numChannels = max(m_numChannels, a->channel_code+1);
 affxcdf::ReplicationType set_rep_type = group.GetRepType();
 if (set_rep_type == affxcdf::UnknownRepType) {
   a->replication_type=Atom::UnknownRepType;
 }
 else if (set_rep_type == affxcdf::DifferentRepType) {
   a->replication_type=Atom::DifferentRepType;
 }
 else if (set_rep_type == affxcdf::MixedRepType) {
   a->replication_type=Atom::MixedRepType;
 }
 else if (set_rep_type == affxcdf::IdenticalRepType) {
   a->replication_type=Atom::IdenticalRepType;
 }
 else {
   Err::errAbort("In ChipLayout::readCdfProbeSet: unknown Atom::ReplicationType: " + ToStr(set_rep_type));
 }
      }
      int numCells = group.GetNumCells();
      int probesPerAtom = numCells / numAtoms;

      if (numCells % numAtoms != 0) {
        Verbose::out(1,READCDFPROBESET_ERR_MSG
                     "Don't support atoms with varying number of probes.\n"
                     "psIx="+ToStr(psIx)+"\n"
                     "name='"+ToStr(name)+"'\n"
                     "numCells="+ToStr(numCells)+"\n"
                     "numAtoms="+ToStr(numAtoms)+"\n"
                     );
        // add to the count and punt
        errorCount++;
        return;
      }

      if (locProbeIndex + numCells >= (*locProbes)->size()) {
        Verbose::out(2, "Resizing probe pm/mm size from: " + ToStr((*locProbes)->size()) + " to " +
                     ToStr((*locProbes)->size() + MAX_PROBES));
        resizeProbes(locProbes, (*locProbes)->size() + numCells + MAX_PROBES, *locAtoms, locAtomIndex);
      }

      int atomProbeIndex = locProbeIndex;
      for (int cellIx = 0; cellIx < probesPerAtom; cellIx++) {
        if (expectingMisMatches) {

          // Allocate pm & mm probes
          Probe *pm = &(**locProbes)[locProbeIndex++];
          Probe *mm = &(**locProbes)[locProbeIndex++];
          // Look up two possibilities.
          group.GetCell(cellCount+cellIx, cel1);
          group.GetCell(cellCount+cellIx+1, cel2);
          cellCount += 2;

          // Sanity check.
          if (cel1.GetListIndex() != cel2.GetListIndex()) {
            Verbose::out(1,"ERROR: readCdfProbeSet(): "
                         "Expecting pm & mm probes to have same list index in cdf.\n"
                         "probeset '"+ToStr(ps.name)+"' "
                         "psIx="+ToStr(psIx)
                         );

            //
            errorCount++;
            return;
          }
   // @todo? check atom rep_type before trying to read mm?
          // make sure we read pm/mm correctly or die trying...
          if (ChipLayout::isPm(cel1)) {
            readCdfProbe(pm, a, numCol, numRow, direction, group, cel1);
            readCdfProbe(mm, a, numCol, numRow, direction, group, cel2);
          }
          else if (ChipLayout::isPm(cel2)) {
            readCdfProbe(pm, a, numCol, numRow, direction, group, cel2);
            readCdfProbe(mm, a, numCol, numRow, direction, group, cel1);
          }
          else {
            Verbose::out(1,READCDFPROBESET_ERR_MSG
                         "Expecting a mismatch probe. "
                         "probeset '"+ToStr(ps.name)+"' "+
                         "psIx="+ToStr(psIx)
                         );
            //
            errorCount++;
            return;
          }
          /* if either of the probes are on the kill list the get tossed like never
             existed. */
          // @todo inKillList should have zero-based indexes stored.
          bool pmInKillList = inKillList(FROM_ZEROBASED(pm->id), ps.name, killList);
          if (pmInKillList) {
            Verbose::out(4,"Ignoring probe " + ToStr(FROM_ZEROBASED(pm->id)) + " in probeset " + ToStr(ps.name) + " due to kill list");
            // as the mm probe is behind the pm probe have to put the mm probe in the pm probes place.
            *pm = *mm;
            clearProbe(*mm);
            mm = pm;
            pm = NULL;
            locProbeIndex--;
          }
          if (inKillList(FROM_ZEROBASED(mm->id), ps.name, killList)) {
            Verbose::out(4,"Ignoring probe " + ToStr(FROM_ZEROBASED(mm->id)) + " in probeset " + ToStr(ps.name) + 
                         " due to kill list");
            clearProbe(*mm);
            locProbeIndex--;
            mm = NULL;
          } else if (pmInKillList) {
            Verbose::out(4,"Ignoring MM probe " + ToStr(FROM_ZEROBASED(mm->id)) + " in probeset " + ToStr(ps.name) +
                         " due to PM in kill list");
            clearProbe(*mm);
            locProbeIndex--;
            mm = NULL;
          }

          // If these probes are to be kept, then make sure we keep the probeset
          if (!probeSubSet.empty() && ((pm != NULL && probeSubSet[pm->id]) || (mm != NULL && probeSubSet[mm->id])))
            keepProbeSet = true;
          // If we got order of pm & mm probes right put them in atoms.
          if ((pm == NULL || pm->type == Probe::PMST || pm->type == Probe::PMAT) &&
             (mm == NULL || mm->type == Probe::MMST || mm->type == Probe::MMAT)) {
            if (pm != NULL) {
              m_PmProbes[pm->id] = true;
              //              a->probes[cellIx] = pm;
            }
            if (mm != NULL && !m_MmProbes.empty()) {
              m_MmProbes[mm->id] = true;
              //              a->probes[cellIx] = mm;
            }
          }
          // Yikes, what did we get? Aborting...
          else {
            Verbose::out(1,READCDFPROBESET_ERR_MSG
                         "Expected to get a pm and mm probe pair in cdf file, but didn't."
                         "probeset '"+ToStr(ps.name)+"' "+
                         "psIx="+ToStr(psIx)
                         );
            //
            errorCount++;
            return;
          }
          // Increment cellIx by +1 as we read in two probes.
          cellIx++;
        }
        else {
          if (cellCount >= numCells) {
            Verbose::out(1, "WARNING: This probe set appears to be malformed in the CDF:\t" + name + 
                         "\tAtomCount:" + ToStr(numAtoms) + "\tCellCount:" + ToStr(cellCount) +
                         "\tProbesPerAtom:" + ToStr(probesPerAtom)); Err::errAbort("Test");
            continue;
          }
          // Just read a singleton.
          if (locProbeIndex + 1 >= (*locProbes)->size()) {
            Verbose::out(2, "Resizing probe pm only size from: " + ToStr((*locProbes)->size()) +
                         " to " + ToStr((*locProbes)->size() + MAX_PROBES));
            resizeProbes(locProbes, (*locProbes)->size() + MAX_PROBES, *locAtoms, locAtomIndex);
          }
          Probe *p = &(**locProbes)[locProbeIndex++];
          //          p->parentAtom = a;
          group.GetCell(cellCount, cel1);
          cellCount++;
          readCdfProbe(p, a, numCol, numRow, direction, group, cel1);
          // If user has chosen to exclude this probe from probeset don't add it.
          if (inKillList(FROM_ZEROBASED(p->id), ps.name, killList)) {
            Verbose::out(4,"Ignoring probe " + ToStr(p->id+1) + " in probeset " + ToStr(ps.name) + " due to kill list");
            clearProbe(*p);
            locProbeIndex--;
          }
          else {
            // If this probes is to be kept, then make sure we keep the probeset
            if (!probeSubSet.empty() && probeSubSet[p->id])
              keepProbeSet = true;
            if (p->type == Probe::PMST || p->type == Probe::PMAT) {
              m_PmProbes[p->id] = true;
            }
            else if ((p->type == Probe::MMAT || p->type == Probe::MMST) && !(m_MmProbes.empty())) {
              m_MmProbes[p->id] = true;
            }
          }
        }
        a->probes.resize(locProbeIndex - atomProbeIndex);
        for (uint32_t pIx = 0; pIx < locProbeIndex - atomProbeIndex; pIx++) {
          a->probes[pIx] = &(**locProbes)[pIx + atomProbeIndex];
        }
      }
      //      a->parentPSet = &ps;
      if (a->probes.size() > 0) {
        ps.atoms.push_back(a);
        goodAtoms++;
      }
    } // end for (int atomIx = 0; atomIx < numAtoms; atomIx++) 
    ps.numGroups++;
    ps.atomsPerGroup[groupIx] = numAtoms;
    if (numAtoms != goodAtoms) {
      ps.atomsPerGroup[groupIx] = goodAtoms;
    }
  } // end for (int groupIx = 0; groupIx < numGroups; groupIx++) 
}

/**
 * Read a vector of ProbeSets from a cdf file representation.
 *
 * @param cdf - cdf file object with probe set information.
 * @param psVec - vector to put newly created probe sets in.
 * @param probeSetsToLoad - which probe sets should be converted.
 * @param probeSubSet - which individual probes were loaded.
 * @param probesetNames - If not null filled in with the
 * probesetname for each probeset even if not actually into
 * chiplayout. Useful for chp files where every single probeset must
 * be included and in order.
 * @param killList - map of probes ids that shouldn't be used.
 * @param justStats - just read and generate stats, don't load into memory.
 * @param psTypesToLoad - What types of probesets to load into memory.
 */
void 
ChipLayout::readProbeListCdfFileKillList(FusionCDFData &cdf,
                                         ProbeListFactory &plFactory,
                                         const std::set<const char *, 
                                         Util::ltstr> &probeSetsToLoad,
                                         std::vector<const char *> *probesetNames,
                                         std::vector<bool> &probeSubSet,
                                         probeidmap_t &killList,
                                         bool justStats,
                                         std::set<affxcdf::GeneChipProbeSetType> &psTypesToLoad)
{
  int numCol = 0, numRow = 0;
  bool keepProbeSet = false;
  FusionCDFProbeSetInformation set;
  FusionCDFFileHeader &cdfHeader = cdf.GetHeader();
  numCol = cdfHeader.GetCols();
  numRow = cdfHeader.GetRows();
  uint32_t numProbes = numCol * numRow;
  resizeStats(numProbes);
  assert(probeSubSet.size() == 0 || probeSubSet.size() == numProbes);
  int expectedProbesets = cdfHeader.GetNumProbeSets();
  //m_PsTypes.reserve(expectedProbesets);

  // Set up some local static storage. Once a probeset is determined
  // to be kept we will allocate memory on the heap for it. Trying to
  // prevent lots of new and deletes as they can fragement memory.
  ProbeSet psConst;
  int psConstMaxName = 256;
  psConst.name = new char[psConstMaxName];
  psConst.name[psConstMaxName -1] = '\0'; // null terminate.
  vector<Atom> *locAtoms = new vector<Atom>(MAX_ATOMS);
  vector<Probe> *locProbes = new vector<Probe>(MAX_PROBES);
  std::set<const char *, Util::ltstr>::const_iterator emptyIter = probeSetsToLoad.end();
  /* Loop through probesets, reading in each one. */
  unsigned int dotMod = max(expectedProbesets/40,1);
  Verbose::progressBegin(1, "Reading " + ToStr(expectedProbesets) + " probesets", 40, dotMod, expectedProbesets);
  unsigned int seenProbesets = 0, keptProbesets = 0;

  int errorCount=0; // the count of errors from readCdfProbeSet

  for (int psIx = 0; psIx < expectedProbesets; psIx++) {
    Verbose::progressStep(1);
    seenProbesets++;
    keepProbeSet = false;
    cdf.GetProbeSetInformation(psIx, set);
    // If we're not loading this type of probeset ignore it totally...
    if (!psTypesToLoad.empty() && 
        psTypesToLoad.find(set.GetProbeSetType()) == psTypesToLoad.end()) {
      continue;
    }
    enum affxcdf::GeneChipProbeSetType type = set.GetProbeSetType();
    if (type != affxcdf::ExpressionProbeSetType &&
        type != affxcdf::GenotypingProbeSetType &&
        type != affxcdf::CopyNumberProbeSetType &&
        type != affxcdf::MarkerProbeSetType     &&
 type != affxcdf::MultichannelMarkerProbeSetType &&
        type != affxcdf::UnknownProbeSetType) {
      Verbose::out(1,"Only know how to parse expression, genotyping, copynumber, marker, multichannel and unknown probesets. "
                   "Got type: " + ToStr(type) + " for probeset: " + ToStr(cdf.GetProbeSetName(psIx)) + ". "
                   "Treating as unknown.");
    }
    readCdfProbeSet(cdf.GetProbeSetName(psIx), 
      set, 
      psIx,
      numCol, 
      numRow,
      probeSubSet, 
      keepProbeSet, 
      killList,
      psConst, 
      psConstMaxName, 
      &locAtoms, 
      &locProbes,
      errorCount);

    if (probesetNames != NULL) {
      probesetNames->push_back(Util::cloneString(psConst.name));
    }

    /* Remember the type in case this probeset isn't actually
       loaded. Convert type to expression if a "false" genotype
       probeset. */
    /*
    //This is now done in ChipLayout::getprobesettype
    if (psConst.psType == ProbeSet::GenoType && psConst.numGroups != 2 && psConst.numGroups != 4) {
    m_PsTypes.push_back((char)ProbeSet::Expression);
    }
    else {
    m_PsTypes.push_back((char)psConst.psType);
    }
    */

    /* If this probe set is to be loaded do so otherwise free up memory. */
    keepProbeSet = true;
    if (!probeSetsToLoad.empty() && 
        probeSetsToLoad.find(psConst.name) == emptyIter)
      keepProbeSet = false;

    /* If all the probes in this probeset were on the kill list don't include it. */
    int probeCount = psConst.countProbes();
    if (probeCount == 0) {
      keepProbeSet = false;
    }
    /// @todo if we arent keeping the probeset which has an error, then maybe we should ignore the error
    ///       on that probeset and not abort.

    m_ProbeCounts.push_back(probeCount);
    int numBytes = ProbeListPacked::byte_size(psConst.numGroups, 
                                              probeCount, 
                                              strlen(psConst.name));
    trackGlobalStats(psConst, numBytes);

    if (keepProbeSet) {
      keptProbesets++;
      for (int aIx = 0; aIx < psConst.atoms.size(); aIx++) {
        for (int pIx = 0; pIx < psConst.atoms[aIx]->probes.size(); pIx++) {
          m_ProbeLayoutOrder.push_back(psConst.atoms[aIx]->probes[pIx]->id);
        }
      }

      if (m_SpfFile.is_open()) {
        writeSpfProbeSetRecord(m_SpfFile, psConst);
      }
      else if (!justStats) {
        for (unsigned int aIx = 0; aIx < psConst.atoms.size(); aIx++) {
          Atom &a = *(psConst.atoms[aIx]);
          for (unsigned int probeIx = 0; probeIx < a.probes.size(); probeIx++) {
            Probe &p = *(a.probes[probeIx]);
            probeSubSet[p.id] = true;
          }
        }
        /* If somebody has specified all the probes on the kill list we
           can have a situation with no probes... */
        //ProbeListPacked pList = plFactory.pListFromPSet(psConst);
        //plVec.push_back(pList); // store the pointer
        plFactory.pListFromPSet(psConst);
      }
    }
    clearProbeSet(psConst);
  }

  Verbose::progressEnd(1, "Done.");
  if (expectedProbesets > 0 && (seenProbesets != expectedProbesets)) {
    Err::errAbort("Expecting: " + ToStr(expectedProbesets) + " probesets, got: " + ToStr(seenProbesets));
  }
  Verbose::out(1,"Kept " + ToStr(keptProbesets) + " probesets.");

  delete locProbes;
  for (uint32_t aIx = 0; aIx < locAtoms->size(); aIx++) {
    (*locAtoms)[aIx].probes.clear();
  }

  psConst.atoms.clear();
  psConst.atomsPerGroup.clear();
  delete locAtoms;

  // opps! we had one or more errors while reading the CDF file.
  if (errorCount!=0) {
    Err::errAbort("readCdfProbeSet(): "
                  "Encountered "+ToStr(errorCount)+" errors while reading CDF file.");
  }
}

const int ChipLayout::mmId(const int probeIx) const {
  if (m_PmMmVec.empty())
    Err::errAbort("ChipLayout::mmId() - Appears that chip layout has not been set, no mismatch probes.");
  if (probeIx >= m_PmMmVec.size())
    Err::errAbort("ChipLayout::mmId() - Probe id: " + ToStr(probeIx + 1) + " is off the chip.");
  const int mmIx = m_PmMmVec[probeIx];
  if (mmIx == NODATA)
    Err::errAbort("ChipLayout::mmId() - Probe id: " + ToStr(probeIx + 1) + " does not have an MM probe assigned.");
  return mmIx;
}

void ChipLayout::resizeAtoms(std::vector<Atom> **old, 
                             int32_t newSize, 
                             ProbeSet *ps) {
  assert((*old)->size() < newSize);
  std::vector<Atom> *larger = new std::vector<Atom>(newSize);
  for (uint32_t i = 0; i < (*old)->size(); i++) {
    Atom *o = &(**old)[i];
    Atom *n = &(*larger)[i];
    n->id = o->id;
    n->probes.resize(o->probes.size());
    for (uint32_t pIx = 0; pIx < n->probes.size(); pIx++) {
      n->probes[pIx] = o->probes[pIx];
      //      n->probes[pIx]->parentAtom = n;
    }
    o->probes.clear();
    clearAtom(*o);
  }
  for (uint32_t aIx = 0; aIx < ps->atoms.size(); aIx++) {
    ps->atoms[aIx] = &(*larger)[aIx];
  }
  delete *old;
  *old = larger;
}



void ChipLayout::resizeProbes(std::vector<Probe> **old, 
                              int32_t newSize,
                              std::vector<Atom> *currentAtoms, 
                              uint32_t currentAtomCount) {
  assert(old);
  assert(currentAtoms);
  uint32_t origSize = (*old)->size();
  assert( origSize < newSize);
  std::vector<Probe> *larger = new std::vector<Probe>(newSize);
  for (uint32_t i = 0; i < (*old)->size(); i++) {
    Probe *o = &(**old)[i];
    Probe *n = &(*larger)[i];
    (*n) = (*o);
  }
  int currentProbeIx = 0;
  for (uint32_t aIx = 0; aIx < currentAtomCount; aIx++) {
    Atom *a = &(*currentAtoms)[aIx];
    for (uint32_t pIx = 0; pIx < a->probes.size(); pIx++) {
      Probe *o = a->probes[pIx];
      Probe *n = &(*larger)[currentProbeIx];
      if (o->id != n->id) {
        Err::errAbort("ChipLayout::resizeProbes() - Expecting id: " + ToStr(o->id) + " but got id " + ToStr(n->id));
      }
      a->probes[pIx] = n;
      currentProbeIx++;
    }
  }
  for (uint32_t i = 0; i < (*old)->size(); i++) {
    Probe *o = &(**old)[i];
    clearProbe(*o);
  }
  delete *old;
  *old = larger;
}


void ChipLayout::clearProbe(Probe &p) {
  p.id = 0;
  p.type = Probe::BLANK;
  p.gcCount = (char)NULLPROBEGC;
  //    p.parentAtom = NULL;
}

void ChipLayout::clearAtom(Atom &a) {
  vector<Probe *>::iterator pIx;
  a.id = 0;
  for (pIx = a.probes.begin(); pIx != a.probes.end(); ++pIx) {
    clearProbe(**pIx);
  }
  a.probes.resize(0);
}

void ChipLayout::clearProbeSet(ProbeSet &ps) {
  vector<Atom *>::iterator aIx;
  for (aIx = ps.atoms.begin(); aIx != ps.atoms.end(); ++aIx) {
    clearAtom(**aIx);
  }
  if (ps.name != NULL) {
    ps.name[0] = '\0';
  }
  ps.psType = ProbeSet::Unknown;
  ps.numGroups = 0;
  ps.atoms.resize(0);
}

/**
 * Keep track of global stats as we parse through file.
 *
 * @param name - Name of probeset
 * @param type - Class of probeset (as in ProbeSet::Type)
 * @param numBlocks - How many blocks are in the probeset.
 * @param numMatch - How many matches are there (1 for pm only, 2 for pm-mm)
 * @param probes - Vector of the probe ids/indexes in the array.
 */
void
ChipLayout::trackGlobalStats(const std::string& name,
                             int type,
                             int numBlocks,
                             const std::vector<int>& blockSizes,
                             int numMatch,
                             const std::vector<int>& probes,
                             bool multAlleles)
{
  m_MaxNameLength = Max((int)name.length(), m_MaxNameLength);
  int pmStart = 0, pmEnd = probes.size() / numMatch;
  int mmStart = probes.size() / numMatch;
  int numBytes = ProbeListPacked::byte_size(numBlocks, probes.size(), name.size()+1);
  int additionalMem = Util::round(numBytes * .15) + sizeof(ProbeList);
  additionalMem += (12 * sizeof(char *) + sizeof(std::pair<const char *,unsigned int>)); // for name map
  additionalMem += 16; // for safety
  m_ProbesPerProbeset.addData(probes.size());
  m_ProbesetMemSizes.addData(numBytes + additionalMem);
  m_ProbeCounts.push_back(probes.size());
  // make probes that are in a probeset.
  for (int pIx = 0; pIx < probes.size(); pIx++) {
    if (probes[pIx] != ProbeList::null_probe)
      m_ProbesetProbes[probes[pIx]] = true;
  }
  // mark PM probes
  for (int pIx = pmStart; pIx < pmEnd; pIx++) {
    if (probes[pIx] != ProbeList::null_probe)
      m_PmProbes[probes[pIx]] = true;
  }
  // mark MM probes if requested.
  if (!m_MmProbes.empty()) {
    for (int pIx = pmStart; pIx < pmEnd; pIx++) {
      if (probes[pIx] != ProbeList::null_probe)
        m_MmProbes[probes[pIx]] = true;
    }
  }
  // if MM probes and requested PM->MM matching make it
  if (numMatch == 2 && !m_PmMmVec.empty()) {
    for (int pIx = pmStart; pIx < pmEnd; pIx++) {
      if (probes[pIx] != ProbeList::null_probe &&
         probes[pIx + mmStart] != ProbeList::null_probe) {
        m_PmMmVec[probes[pIx]] = probes[pIx+mmStart];
      }
    }
  }
  incrementTypes(type, numBlocks, multAlleles);
}

void ChipLayout::incrementTypes(int type, int numBlocks, bool multipleAlleles) {
  if (type == ProbeSet::Expression || 
      (type == ProbeSet::GenoType && numBlocks != 2 && numBlocks != 4)) {
    m_NumExpression++;
  }
  else if ( multipleAlleles && 
            (type == ProbeSet::GenoType || 
             type == ProbeSet::Marker ||
             type == ProbeSet::MultichannelMarker)) {
    m_NumGenotyping++;
  }
}

void ChipLayout::trackGlobalStats(ProbeSet &ps, int numBytes) {
  Atom *atom = NULL;
  Probe *p = NULL;
  Probe *pm = NULL, *mm = NULL;
  m_MaxNameLength = Max((int)strlen(ps.name), m_MaxNameLength);
  bool multAllele = false;
  if (!ps.atoms.empty()) {
    int first = ps.atoms[0]->allele_code;
    for (int i = 1; i < ps.atoms.size(); i++) {
      if (ps.atoms[i]->allele_code != first) {
        multAllele = true;
      }
    }
  }
  incrementTypes(ps.psType, ps.atomsPerGroup.size(), multAllele);

  // 15% more for heap, pointer to name, name and probelist itself.
  int additionalMem = Util::round(numBytes * .15) + sizeof(ProbeList);
  additionalMem += (12 * sizeof(char *) + sizeof(std::pair<const char *,unsigned int>)); // for name map
  additionalMem += 16; // for safety

  m_ProbesPerProbeset.addData(ps.countProbes());
  m_ProbesetMemSizes.addData(numBytes + additionalMem);

  //m_ProbeCounts handled by caller
  for (uint32_t atomIx = 0; atomIx < ps.atoms.size(); atomIx++) {
    uint32_t probeIx = 0;
    atom = ps.atoms[atomIx];
    /* Do indivitual probe stats. */
    for (probeIx = 0; probeIx < atom->probes.size(); probeIx++) {
      p = atom->probes[probeIx];
      /* Store the gc count for this particular probe. */
      if (!m_ProbeGcVec.empty()) {
        if (p->id>=m_ProbeGcVec.size()) {
          Err::errAbort("trackGlobalStats: Internal error: id ("+ToStr(p->id)+")>=ProbeGcVec.size()");
        }
        m_ProbeGcVec[p->id] = p->gcCount;
      }
      /* Keep track of pm and mm probes. */
      if (p->type == Probe::PMAT || p->type == Probe::PMST) {
        m_PmProbes[p->id] = true;
      }
      else if ((p->type == Probe::MMAT || p->type == Probe::MMST) && 
               !m_MmProbes.empty()) {
        m_MmProbes[p->id] = true;
      }
      if (ps.psType != ProbeSet::Unknown) {
        m_ProbesetProbes[p->id] = true;
      }
    }

    /* Do the pm->mm mapping. */
    for (probeIx = 0; probeIx < atom->probes.size() && probeIx < atom->probes.size() - 1; probeIx++) {
      pm = atom->probes[probeIx];
      mm = atom->probes[probeIx+1];
      if (pm->type == Probe::PMST && mm->type == Probe::MMST && !m_PmMmVec.empty())
        m_PmMmVec[pm->id] = mm->id;
      if (pm->type == Probe::PMAT && mm->type == Probe::MMAT && !m_PmMmVec.empty())
        m_PmMmVec[pm->id] = mm->id;
    }
  }

  if (ps.psType == ProbeSet::GenoType && !m_PmAlleleMatchVec.empty()) {
    int atomIndex = 0;
    /* Genotype probe sets can have either 2 or 4 groups, for each
       type the even (0,2) are A alleles and odd (1,3) are B allele.s */
    for (int psIx = 0; psIx < ps.numGroups - 1; psIx += 2) {
      Err::check(ps.atomsPerGroup[psIx] == ps.atomsPerGroup[psIx+1],
                 "ChipLayout::trackGlobalStats() - Must have matched A and B allele probes when using this analysis. Probeset '" + ToStr(ps.name) + "' has unmatched probe count.");
      int atomsPer = ps.atomsPerGroup[psIx];
      /* Each block will have atomsPer atoms in there. */
      for (int atomIx = 0; atomIx < atomsPer; atomIx++) {
        int aIndex = atomIx + (atomIndex);
        int bIndex = aIndex + atomsPer;
        Atom *A = ps.atoms[aIndex];
        Atom *B = ps.atoms[bIndex];
        if (A->probes.size() != B->probes.size()) {
          Err::errAbort("For probeset '" + ToStr(ps.name) + " A and B alleles must have same number of probes.");
        }
        for (int probeIx = 0; probeIx < A->probes.size(); probeIx++) {
          m_PmAlleleMatchVec[A->probes[probeIx]->id] = B->probes[probeIx]->id;
          m_PmAlleleMatchVec[B->probes[probeIx]->id] = A->probes[probeIx]->id;
        }
      }
    }
  }
}

/**
 * Reads in probes to kill.
 *
 * @param fileName - file containing BadProbes
 * @param killList - things to kill, determined from file
 * @param rows - number of rows in cel file
 * @param cols - number of cols in cel file
 */
void ChipLayout::fillInKillList(const std::string& fileName, 
                                probeidmap_t &killList, int rows, int cols) {

  assert(fileName!="");
  assert(rows > 0);
  assert(cols > 0);

  killList.clear();
  affx::TsvFile tsv;
  int count = 0;
  probeidmap_t::iterator mapIter;
  vector<string>::iterator psIter;

  bool haveProbeId = false;
  bool haveXY = false;

  string psName;
  int probeId;
  int x,y;

  if (tsv.open(fileName) != TSV_OK)
    Err::errAbort("Couldn't open file: " + ToStr(fileName) + " to read.");

  if (tsv.cname2cidx(0, "probeset_id") == TSV_ERR_NOTFOUND)
    Err::errAbort("Probe kill list file, " + ToStr(fileName) + ", does not have probeset_id column.");
  else
    tsv.bind(0, "probeset_id", &psName, TSV_BIND_REQUIRED);

  // We expect either a probe_id column or an x/y column.
  if (tsv.cname2cidx(0, "probe_id") != TSV_ERR_NOTFOUND) {
    tsv.bind(0, "probe_id", &probeId, TSV_BIND_REQUIRED);
    haveProbeId = true;
  }
  if ((tsv.cname2cidx(0, "x") != TSV_ERR_NOTFOUND) && (tsv.cname2cidx(0, "y") != TSV_ERR_NOTFOUND)) {
    tsv.bind(0, "x", &x, TSV_BIND_REQUIRED);
    tsv.bind(0, "y", &y, TSV_BIND_REQUIRED);
    haveXY = true;
  }

  if (!haveProbeId && !haveXY)
    Err::errAbort("Kill list file must have either a probe_id column or an x and y column.");

  while (tsv.nextLevel(0) == TSV_OK) {

    int probeIdFromXY = -1;
    if (haveXY) {
      int rv=affymetrix_fusion_io::FusionCELData::XYToIndex(x,y,rows,cols);
      probeIdFromXY = FROM_ZEROBASED(rv);
    }

    if (haveProbeId) {
      if (haveXY && probeIdFromXY != probeId) {
        Err::errAbort("ProbeID value (" + ToStr(probeId) + ") does not match x,y (" +
                      ToStr(x) + "," + ToStr(y) + ") based ID (" + ToStr(probeIdFromXY) + ").");
      }
    } else if (haveXY) {
      probeId = probeIdFromXY;
    } else {
      Err::errAbort("In unexpected part of code. No XY or ProbeID in kill list.");
    }

    mapIter = killList.find((uint32_t)probeId);
    if (mapIter != killList.end()) {
      psIter = find(killList[(uint32_t)probeId].begin(),killList[(uint32_t)probeId].end(),psName);
      if (psIter != killList[(uint32_t)probeId].end()) {
        Err::errAbort("Probe id: " + ToStr(probeId) + " specified twice in file: " + fileName);
      } else {
        killList[(uint32_t)probeId].push_back(psName);
        count++;
      }
    } else {
      killList[(uint32_t)probeId].push_back(psName);
      count++;
    }
  }
  Verbose::out(1, "Loaded: " + ToStr(count) + " probe/probeset pairs to be ignored from: " + Fs::basename(fileName));
  tsv.close();
}

/**
 * @brief compare probeset names according to respective probeset ids
 * @param a, b - probeset names to be ordered
 */
// static bool ChipLayout::compareProbesetNameByProbesetIdx(std::string a, std::string b) {
//   return (getProbeSetIndexByName(a) < getProbeSetIndexByName(b));
// }

void ChipLayout::sortProbesetNameByProbesetIdx(std::vector<std::string>::iterator begin, std::vector<std::string>::iterator end) {
  if (begin < end - 1) {
    std::vector<std::string>::iterator right = end - 1;
    std::vector<std::string>::iterator left = begin;
    std::string temp;
    int pivot = getProbeSetIndexByName(*right);
    
    while (left < right) {
      while (left < right && getProbeSetIndexByName(*left) <= pivot) {
 left++;
      }
      while (left < right && pivot < getProbeSetIndexByName(*right)) {
 right--;
      }
      
      if (left < right) {
 temp = *left;
 *left = *right;
 *right = temp;
      }
    }
    
    sortProbesetNameByProbesetIdx(begin, right);
    sortProbesetNameByProbesetIdx(left, end);
  }
}
