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

#include "chipstream/AnalysisStreamExpression.h"

AnalysisStreamExpression::AnalysisStreamExpression(bool doGenoTypes, bool doStrand, bool aOnly) {
    setupSelfDoc(*this);
    m_DoGenotypes = doGenoTypes;
    m_DoStrand = doStrand;
    m_Aonly = aOnly;
    setOptValue("genotype", m_DoGenotypes);
    setOptValue("strand", m_DoStrand);
    setOptValue("allele-a", m_Aonly);
}

/** 
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> AnalysisStreamExpression::getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> opts;
    SelfDoc::Opt doGenoType = {"genotype", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                          "Should genotype probesets also be summarized. [experimental] "};
    opts.push_back(doGenoType);

    SelfDoc::Opt doStrand = {"strand", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                          "Should genotype probesets be separated by strand or just A and B alleles. [experimental]"};
    opts.push_back(doStrand);
    SelfDoc::Opt aOnly = {"allele-a", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                          "Should only the A allele be summarized? (Useful with pm-sum where A and B give same results). [experimental]"};
    opts.push_back(aOnly);
    return opts;
}

/** 
 * Fill in the information for Self documentation.
 * @param doc - Self documenter to be filled in.
 */
void AnalysisStreamExpression::setupSelfDoc(SelfDoc &doc) {
    doc.setDocName("expr");
    doc.setDocDescription("Does expression summarization on probesets.");
    doc.setDocOptions(getDefaultDocOptions());
}

QuantExprMethod *AnalysisStreamExpression::getQuantExprMethod() {
    return m_QExprMethod; 
}

/** 
 * Set the object for summarizing an individual probe set.
 * @param qMethod - Quantification object.
 */
void AnalysisStreamExpression::setQuantMethod(QuantMethod *qMethod) {
    QuantExprMethod *qeMethod = dynamic_cast<QuantExprMethod *>(qMethod);
    if(qeMethod == NULL)
      Err::errAbort("AnalysisStreamExpression::setQuantMethod() - Can only set QuantExprMethod's in AnalysisStreamExpression.");
    m_QMethod = qeMethod;
    m_QExprMethod = qeMethod;
}

/** 
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc AnalysisStreamExpression::explainSelf() { 
    SelfDoc doc;
    setupSelfDoc(doc);
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
SelfCreate *AnalysisStreamExpression::newObject(std::map<std::string,std::string> &param) {
    SelfDoc doc = explainSelf();
    bool doGenotype = false;
    bool doStrand = false;
    bool aOnly = false;
    fillInValue(doGenotype, "genotype", param, doc);
    fillInValue(doStrand, "strand", param, doc);
    fillInValue(aOnly, "allele-a", param, doc);
    AnalysisStreamExpression *stream = new AnalysisStreamExpression(doGenotype, doStrand, aOnly);
    return stream;
}

/** 
 * Fill in psNewGroup with probes from psGroup based on the block
 * indexes in the blocks vector. To divide a genotyping probeset
 * based on Allele (A or B) you would use the 1,3 blocks for A and
 * 2,4 blocks for B. If you just want one strand of a particular
 * allele then it would be 1 for A+, 2 for B+, 3 for A- and 4 for
 * B-. This is sensitive to the vagaries of the genotyping probeset
 * structure which is poorly documented and ever
 * "innovating"... Eventually this should be moved to the ChipLayout
 * class or ProbeSet class.
 * 
 * @param psNewGroup 
 * @param psGroup 
 * @param blocks 
 */
void AnalysisStreamExpression::fillInProbeSetFromBlocks(ProbeSetGroup &psNewGroup, ProbeSetGroup &psGroup,
                              const std::vector<short> &blocks) {
    std::set<probeid_t> toKeep;
    const ProbeSet *ps = psGroup.probeSets[0];
    for(uint32_t blockIx = 0; blockIx < blocks.size(); blockIx++) {
      int start = 0;
      // figure out where to start for this block
      for(uint32_t groupIx = 0; groupIx < blocks[blockIx]; groupIx++) {
        start += ps->atomsPerGroup[groupIx];
      }
      for(uint32_t bIx = start; bIx < ps->atomsPerGroup[blocks[blockIx]] + start; bIx++) {
        Atom *a = ps->atoms[bIx];
        for(uint32_t pIx = 0; pIx < a->probes.size(); pIx++) {
          toKeep.insert(a->probes[pIx]->id);
        }
      }
    }
    makeSelectProbesetGroup(psNewGroup, toKeep, psGroup);
}

void AnalysisStreamExpression::fillInGroup(ProbeSetGroup &psAllele, ProbeSetGroup &psGroup, int blockIx, const char *basename, const char *suffix) {
    const ProbeSet *ps = psGroup.probeSets[0];
    std::vector<short> blocks;
    blocks.push_back(blockIx);
    fillInProbeSetFromBlocks(psAllele, psGroup, blocks);
    std::string name = ps->name + ToStr(suffix);
    FreezArray(psAllele.probeSets[0]->name);
    psAllele.name = psAllele.probeSets[0]->name = Util::cloneString(name.c_str());
    psAllele.probeSets[0]->psType = ProbeSet::Expression;
}

/** 
 * Do the analysis for a particular group of probe sets. If we see a
 * genotype probeset and user wishes it try to break it up into
 * expression probesets. The code for dividing up a genotyping
 * probeset should probably live in ChipLayout eventually, but for
 * now the initial implementation is here.
 * 
 * @param psGroup - Collection of probe sets to get probes from.
 * @param layout - How probes/probesets are laid out on chip.
 * @param iMart - Object containing raw data values for all chips.
 * @param doReport - Should the quantification report object be called?
 * @param alleleSummaryOnly - this is a parameter that makes sense in the base class AnalysisStream::doAnalysis, but not here.  Included here only to make the inheritance work.  Feel free to ignore.
 * 
 * @return true if success, false otherwise.
 */
bool AnalysisStreamExpression::doAnalysis(ProbeSetGroup &psGroup, IntensityMart &iMart, bool doReport, bool alleleSummaryOnly) {
    bool success = true;
    std::vector<ProbeSetGroup *> toRun;
    std::vector<ProbeSetGroup> psAlleles;
    if(psGroup.probeSets.size() <= 0) {
      Err::errAbort("Can't have probeset groups with no probesets.");
    }
    ///@todo handle Marker type probesets
    if(psGroup.probeSets[0]->psType == ProbeSet::GenoType && m_DoGenotypes &&
       (psGroup.probeSets[0]->numGroups == 2 || psGroup.probeSets[0]->numGroups == 4)) {
      const ProbeSet *ps = psGroup.probeSets[0];

      if(m_DoStrand && ps->numGroups == 4) {

        psAlleles.resize(4);
        short blockIx = 0;
        fillInGroup(psAlleles[blockIx], psGroup, blockIx, ps->name, "-A-F");
        toRun.push_back(&psAlleles[blockIx++]);

        if(!m_Aonly) {
          fillInGroup(psAlleles[blockIx], psGroup, blockIx, ps->name, "-B-F");
          toRun.push_back(&psAlleles[blockIx++]);
        }
        fillInGroup(psAlleles[blockIx], psGroup, blockIx, ps->name, "-A-R");
        toRun.push_back(&psAlleles[blockIx++]);
        if(!m_Aonly) {
          fillInGroup(psAlleles[blockIx], psGroup, blockIx, ps->name, "-B-R");
          toRun.push_back(&psAlleles[blockIx++]);
        }
      }
      else {
        std::string name;
        psAlleles.resize(2);
        std::vector<short> blocksA, blocksB;
        short blockIx = 0;
        blocksA.push_back(blockIx++);
        blocksB.push_back(blockIx++);
        if(ps->numGroups == 4) {
          blocksA.push_back(blockIx++);
          blocksB.push_back(blockIx++);
        }

        fillInProbeSetFromBlocks(psAlleles[0], psGroup, blocksA);
        blockIx = 0;
        name = ps->name + ToStr("-A");
        //
        FreezArray(psAlleles[blockIx].name);
        psAlleles[blockIx].name = Util::cloneString(name.c_str());
        //
        FreezArray(psAlleles[blockIx].probeSets[0]->name);
        psAlleles[blockIx].probeSets[0]->name = Util::cloneString(name.c_str());
        psAlleles[blockIx].probeSets[0]->psType = ProbeSet::Expression;
        //
        toRun.push_back(&psAlleles[blockIx]);

        fillInProbeSetFromBlocks(psAlleles[1], psGroup, blocksB);
        if(!m_Aonly) {
          blockIx = 1;
          name = ps->name + ToStr("-B");
          //
          FreezArray(psAlleles[blockIx].name);
          psAlleles[blockIx].name = Util::cloneString(name.c_str());
          //
          FreezArray(psAlleles[blockIx].probeSets[0]->name);
          psAlleles[blockIx].probeSets[0]->name = Util::cloneString(name.c_str());
          psAlleles[blockIx].probeSets[0]->psType = ProbeSet::Expression;
          //
          toRun.push_back(&psAlleles[blockIx]);
        }
      }
    }
    else if((psGroup.probeSets[0]->psType == ProbeSet::Marker || psGroup.probeSets[0]->psType == ProbeSet::MultichannelMarker) && m_DoGenotypes) {
      if(m_DoStrand) {
          psAlleles.resize(psGroup.probeSets[0]->atoms.size());
          for(int blockIx = 0; blockIx < psGroup.probeSets[0]->atoms.size(); blockIx++) {
              std::string suffix = "-" + ToStr(psGroup.probeSets[0]->atoms[blockIx]->getAlleleCode()) +
                  "-" + ToStr(psGroup.probeSets[0]->atoms[blockIx]->getContextCode());
              fillInGroup(psAlleles[blockIx], psGroup, blockIx, psGroup.probeSets[0]->name, suffix.c_str());
              toRun.push_back(&psAlleles[blockIx]);
          }
      } else {
          std::string name;
          std::vector< std::vector<short> > blocks;

          for(short blockIx = 0; blockIx < psGroup.probeSets[0]->atoms.size(); blockIx++) {
              int allele = psGroup.probeSets[0]->atoms[blockIx]->getAlleleCode();
              if(blocks.size() <= allele)
                blocks.resize(allele+1);
              blocks[allele].push_back(blockIx);
          }

          psAlleles.resize(blocks.size());
          for(int allele=0; allele<psAlleles.size(); allele++) {
              if(blocks[allele].size() > 0) {
                fillInProbeSetFromBlocks(psAlleles[allele], psGroup, blocks[allele]);
                name = ToStr(psGroup.probeSets[0]->name) + "-" + ToStr(allele);
                FreezArray(psAlleles[allele].name);
                psAlleles[allele].name = Util::cloneString(name.c_str());
                FreezArray(psAlleles[allele].probeSets[0]->name);
                psAlleles[allele].probeSets[0]->name = Util::cloneString(name.c_str());
                psAlleles[allele].probeSets[0]->psType = ProbeSet::Expression;
                toRun.push_back(&psAlleles[allele]);
              }
          }
      }
    }
    else {
      toRun.push_back(&psGroup);
    }
    for(uint32_t gIx = 0; gIx < toRun.size(); gIx++) {
      ProbeSetGroup &group = *toRun[gIx];
      //std::cout << *(group.probeSets[0]) << '\n';
      if(m_QMethod->setUp(group, iMart, m_CStreams, *m_PmAdjust)) {
        m_QMethod->computeEstimate();
        if(doReport) {
          for(unsigned int i = 0; i < m_Reporters.size(); i++) {
            m_Reporters[i]->report(group, *m_QMethod, iMart, m_CStreams, *m_PmAdjust);
          }
        }
      }
      else {
        Verbose::out(5, "Warning setup failed for name: " + ToStr(group.name));
        success = false;
      }
      if(!success && doReport) {
        for(unsigned int i = 0; i < m_Reporters.size(); i++) {
          m_Reporters[i]->reportFailure(group, *m_QMethod, iMart, m_CStreams, *m_PmAdjust);
        }
      }
    }
    for(uint32_t gIx = 0; gIx < psAlleles.size(); gIx++) {
      ProbeSetGroup &pG = psAlleles[gIx];
      pG.ownMem = false;
      FreezArray(psAlleles[gIx].name);
      for(uint32_t psIx = 0; psIx < pG.probeSets.size(); psIx++) {
        delete pG.probeSets[psIx];
      }
    }
    return success;
}

/** 
 * Make a new probeset group which is a subset of the original based on the
 * probeset ids contained in the goodIds set.
 * 
 * @param selectGroup - New probeset group to be filled in.
 * @param goodIds - Set with ids of probes to keep.
 * @param psGroup - Original probeset to get probes and structure from.
 */
void AnalysisStreamExpression::makeSelectProbesetGroup(ProbeSetGroup &selectGroup, 
                                      std::set<probeid_t> &goodIds,
                                      ProbeSetGroup &psGroup) {
    selectGroup.ownMem = true;
    selectGroup.name = Util::cloneString(psGroup.name);
    selectGroup.probeSets.resize(psGroup.probeSets.size());
    selectGroup.probeSetNames.resize(psGroup.probeSetNames.size());
    for(uint32_t nameIx = 0; nameIx < psGroup.probeSetNames.size(); nameIx++) {
      selectGroup.probeSetNames[nameIx] = Util::cloneString(psGroup.probeSetNames[nameIx]);
    }
    for(uint32_t psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
      const ProbeSet *ps = psGroup.probeSets[psIx];
      ProbeSet *newPs = new ProbeSet();
      newPs->name = Util::cloneString(ps->name);
      newPs->psType = ps->psType;
      newPs->numGroups = ps->numGroups;
      newPs->atomsPerGroup.resize(ps->atomsPerGroup.size());
      std::copy(ps->atomsPerGroup.begin(), ps->atomsPerGroup.end(), newPs->atomsPerGroup.begin());
      newPs->atoms.resize(ps->atoms.size());
      for(uint32_t atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
        Atom *a = ps->atoms[atomIx];
        Atom *newAtom = new Atom();
        std::vector<Probe *> goodProbes;
        goodProbes.reserve(a->probes.size());
        newAtom->id = a->id;
        bool expectingMismatches = false;
        // if we have at two probes, a multiple of two number of probes
        // and the first two probes are a pm & mm we think there should
        // be mm probes.
        if(a->probes.size() > 1 && a->probes.size() % 2 == 0 && 
           (Probe::isPm(*(a->probes[0])) && Probe::isMm(*(a->probes[1])))) {
          expectingMismatches = true;
        }
        for(uint32_t probeIx = 0; probeIx < a->probes.size(); probeIx++) {
          Probe *p = a->probes[probeIx];
          if(goodIds.find(p->id) != goodIds.end()) {
            Probe *newProbe = Probe::cloneProbe(*p);
            goodProbes.push_back(newProbe);
            // If we're expecting mismatches add the next probe as it should be MM of PM/MM pair.
            if(expectingMismatches) {
              probeIx++;
              if(probeIx >= a->probes.size()) {
                Err::errAbort("makeSelectProbesetGroup() - Expecting a mismatch for probeset group: " + 
                              ToStr(psGroup.name) + " but didn't have enough probes.");
              }
              Probe *newMMProbe = Probe::cloneProbe(*(a->probes[probeIx]));
              if(!Probe::isMm(*newMMProbe)) {
                Err::errAbort("makeSelectProbesetGroup() - Expecting a mismatch for probeset group: " + 
                              ToStr(psGroup.name) + " got non-MM probe instead.");
              }
              goodProbes.push_back(newMMProbe);
            }
          }
        }
        newAtom->probes.resize(goodProbes.size());
        std::copy(goodProbes.begin(), goodProbes.end(), newAtom->probes.begin());
        newPs->atoms[atomIx] = newAtom;
      }
      selectGroup.probeSets[psIx] = newPs;
    }
}
