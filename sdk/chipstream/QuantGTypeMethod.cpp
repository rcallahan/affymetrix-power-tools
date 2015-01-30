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
 * @file   QuantGTypeMethod.cpp
 * @author Earl Hubbell
 * @date   Fri Feb 24 11:02:09 2006
 *
 *  Placeholder: should isolate genotyping probe set summary functions
 *
 */

//
#include "chipstream/QuantMethodFactory.h"
//
#include "file/TsvFile/TsvFile.h"
#include "util/Convert.h"
#include "util/Util.h"
#include "util/Verbose.h"

using namespace affx;

QuantGTypeMethod::~QuantGTypeMethod() {}

/** 
 * Get the genotype call at specified index (sample).
 * @param index - sample of interest.
 * @return - Genotyping call made.
 */

affx::GType QuantGTypeMethod::getGTypeCall(unsigned int index) {
  affx::GType call = getGTypeForcedCall(index);
  double conf = getConfidence(index);
  if((conf >= getMinThresh()) && (conf <= getMaxThresh()))
    return call;
  else
    return affx::NN;
}

affx::GType QuantGTypeMethod::getGTypeForcedCall(unsigned int index) {
  if(index >= m_Calls.size()) {
    Err::errAbort("Asking for call at index " + ToStr(index) + " when Probeset " + 
                  m_ProbesetName + " has only " + ToStr(m_Calls.size()) + " calls.");
  }
  return m_Calls[index];
}

/** 
 * Get the genotype call at specified index (sample).
 * @param index - sample of interest.
 * @return - Genotyping call made.
 */

int QuantGTypeMethod::getForcedCall(unsigned int index) {
  return affx::GType_to_int(getGTypeForcedCall(index));
}

/** 
 * Get the genotype call at specified index (sample).
 * @param index - sample of interest.
 * @return - Genotyping call made.
 */

int QuantGTypeMethod::getCall(unsigned int index) {
  return affx::GType_to_int(getGTypeCall(index));
}

/** 
 * Get a string representation for a particular SNP.
 * @return - a string version of a SNP
 */
std::string QuantGTypeMethod::getModelString() {
  return "";
}

/** summarize this allele - no reports */
bool QuantGTypeMethod::summarizeAllele(ProbeSet *pSet,
                                       std::vector<double> &summaryValues,
                                       const IntensityMart &iMart,
                                       std::vector<ChipStream *> &iTrans, 
                                       PmAdjuster &pmAdjust, 
                                       QuantExprMethod *quantMethod,
                                       bool lowPrecision) {
  std::vector<QuantMethodReport *> nullReporters;
  return summarizeAllele(pSet, summaryValues, iMart, iTrans, 
                         pmAdjust, quantMethod, lowPrecision, 
                         false, nullReporters);
}

/** clear a probeset structure */
void QuantGTypeMethod::clearProbeSet(ProbeSet &ps) {
    delete [] ps.name;
    ps.name = NULL;
    ps.atoms.clear();
};

/** summarize an allele using expression summary methods specified*/
bool QuantGTypeMethod::summarizeAllele(ProbeSet *pSet,
                     vector<double> &summaryValues,
                     const IntensityMart &iMart,
                     std::vector<ChipStream *> &iTrans,
                     PmAdjuster &pmAdjust,
                     QuantExprMethod *quantMethod,
                     bool lowPrecision,
                     bool doReport,
                     std::vector<QuantMethodReport *> reporters) {
  bool success = true;
  summaryValues.clear();
  /* This is a little interesting as we are spoofing a probeSetGroup for calling
     functions but don't want the destructor called while that group thinks it
     owns those probesets. */
  ProbeSetGroup group(pSet);
  if(!quantMethod->setUp(group, iMart, iTrans, pmAdjust))
    success = false;
  else {
    quantMethod->computeEstimate();

    for(unsigned int i = 0; i < quantMethod->getNumTargets(); i++) {
      double val = quantMethod->getTargetEffect(i);
      if(lowPrecision) {
        val = Util::round(val);
      }
      summaryValues.push_back(val);
    }
    success = true;
  }
  /* print out a report if necessary. */
  if(success && doReport) {
    for(unsigned int i = 0; i < reporters.size(); i++) {
      reporters[i]->report(group, *quantMethod, iMart, iTrans, pmAdjust);
    }
  }
  /* clear out probesets to prevent deleting them. */
  group.probeSets.clear();
  clearProbeSet(*pSet);
  return success;
}

/** fill in both probe sets with intensities for both alleles */
bool QuantGTypeMethod::fillInAlleleProbeSets(const ProbeSet &gtPs, ProbeSet &aAllele, ProbeSet &bAllele) {
  /* First set up the probe sets to do the summarization of a and b alleles. */
  aAllele.psType = ProbeSet::Expression;
  bAllele.psType = ProbeSet::Expression;
  // Set the names for our spoofed probesets with suffix A and B
  string name = gtPs.name;
  if (gtPs.psType == ProbeSet::Copynumber)
  {
		aAllele.name = Util::cloneString(name.c_str());
		int groupOffset = 0; // atoms are stored in blocks for each allele, keep track of which block we're in.
		for(unsigned int groupIx = 0; groupIx < gtPs.numGroups; groupIx++) 
		{
			for(unsigned int atomIx = groupOffset; atomIx < gtPs.atomsPerGroup[groupIx] + groupOffset; atomIx++) 
			{
				aAllele.atoms.push_back(gtPs.atoms[atomIx]);
			}
			groupOffset += gtPs.atomsPerGroup[groupIx];
		}
  }
  else if (gtPs.psType == ProbeSet::GenoType || gtPs.psType == ProbeSet::Marker || gtPs.psType == ProbeSet::MultichannelMarker)
  {
	  name += "-A";
	  aAllele.name = Util::cloneString(name.c_str());
	  name = gtPs.name;
	  name += "-B";
	  bAllele.name = Util::cloneString(name.c_str());
	  ProbeSet* alleleToGetAtom = NULL;
	  int groupOffset = 0; // atoms are stored in blocks for each allele, keep track of which block we're in.
	  for(unsigned int groupIx = 0; groupIx < gtPs.numGroups; groupIx++) {
	    /* Probe belongs to the A allele, put atoms in that probeset. */
	    if (gtPs.atoms[groupIx]->getAlleleCode() == 0) {
	      alleleToGetAtom = &aAllele;
	    }
	    /* Probe belongs to the B allele, put atoms in that probeset. */
	    else if (gtPs.atoms[groupIx]->getAlleleCode() == 1) {
	      alleleToGetAtom = &bAllele;
	    }
	    /* Uh-oh, if it doesn't belong to A and it doesn't belong to B we're going to complain. */
	    else {
	      return false;
	    }

	    for(unsigned int atomIx = groupOffset; atomIx < gtPs.atomsPerGroup[groupIx] + groupOffset; atomIx++) {
	      alleleToGetAtom->atoms.push_back(gtPs.atoms[atomIx]);
	    }
	    groupOffset += gtPs.atomsPerGroup[groupIx];
	  }
  }
  else {
      return false;
  }
  return true;
}


/** extract the probes for an allele specified*/
bool QuantGTypeMethod::extractProbes(vector <vector <double> > &probeIntensity,
                   ProbeSet *pSet,
                   vector< unsigned int > &probeIds,
                   const IntensityMart &iMart, std::vector<ChipStream *> &iTrans,
                   PmAdjuster &pmAdjust, bool lowPrecision,
                   QuantExprMethod *quantMethod) {
  bool success = true;
  unsigned int i,j,NTarg,NFeat;
  double val;
  /* This is a little interesting as we are spoofing a probeSetGroup for calling
     functions but don't want the destructor called while that group thinks it
     owns those probesets. */
  ProbeSetGroup group(pSet);

  if(!quantMethod->setUp(group, iMart, iTrans, pmAdjust))
    success = false;
  else {
    NTarg = quantMethod->getNumTargets();
    NFeat = quantMethod->getNumFeatures();

    ///@todo need to populate probe IDs
    probeIds.assign(NFeat,0);

    probeIntensity.clear(); // empty the register
    // organize this by probe across experiments
    // because that's how we'll use this.
    if(probeIntensity.empty()) {
      while(probeIntensity.size() < NFeat) {
        probeIntensity.push_back(vector<double>(NTarg));
      }
    }
    for(i = 0; i < NTarg; i++) {
      for (j=0; j<NFeat; j++) {
        val = quantMethod->getPMDataAt(j,i);
        if(lowPrecision)
          val = Util::round(val);
        probeIntensity[j][i] = val;
        const Probe *p = quantMethod->getFeature(j);
        probeIds[j] = p->id;
      }
    }
    success = true;
  }
  /* clear out probesets to prevent deleting them. */
  group.probeSets.clear();
  clearProbeSet(*pSet);
  return success;
}

/**
 * Can the ProbeSet set up the intensities and summaries
 * @param const ProbeSet* - The probe set to be checked.
 * @return bool - true if can be setup else false. 
 */
bool QuantGTypeMethod::canSetUpProbeSet(const ProbeSet* gtPs)
{
	if (gtPs->psType != ProbeSet::Copynumber)
	{
	        if (gtPs->psType != ProbeSet::GenoType && gtPs->psType != ProbeSet::Marker && gtPs->psType != ProbeSet::MultichannelMarker) 
		{
			Verbose::out(4, "Expecting to get genotyping probesets only. (probeset: name=" + ToStr(gtPs->name) + ", type=" + ProbeSet::typeEnumToString(gtPs->psType) + ")");
			return false; // We only process genotyping probesets.
		}
                if (gtPs->psType == ProbeSet::GenoType && gtPs->numGroups != 4 && gtPs->numGroups != 2) 
		{
			Verbose::out(4, "Expecting to get two or four groups in genotyping probeset. (probeset: name=" + ToStr(gtPs->name) + ", numgroups=" + ToStr((int)gtPs->numGroups) + ")");
			return false;
		}
	}
	return true;
}
