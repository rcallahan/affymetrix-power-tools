////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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
#include "chipstream/QuantBirdseed.h"
//
#include "birdseed-dev/PriorsReader.h"
#include "broadutil/BroadUtil.h"
//
#include "chipstream/QuantMethodFactory.h"
#include "chipstream/QuantMethodReport.h"
#include "chipstream/GenotypeInfo.h"
//
#include "util/Fs.h"


/** Constructor, currently creates prior estimates from pre-calculated data from R. */
QuantBirdseed::QuantBirdseed(double confidenceThreshold, double correctionFactor):
  m_ConfidenceThreshold(confidenceThreshold),
  m_CorrectionFactor(correctionFactor),
  verbosity(0)
{
  m_QuantMethod = NULL;
}


QuantBirdseed::~QuantBirdseed()
{
	for (std::vector<QuantMethodReport *>::iterator it = m_Reporters.begin(); it != m_Reporters.end(); it++)
    {
		delete *it;
	}
	m_Reporters.clear();
    if (m_ClusterOstrm.get() != NULL) {
        m_ClusterOstrm->close();
    }
}


// Copied from QuantBRLMM
void QuantBirdseed::blankSelf() {
    clearProbeSet(m_Aallele);
    clearProbeSet(m_Ballele);
    m_AValues.clear();
    m_BValues.clear();
    //m_InitialCalls.clear();
    m_Calls.clear();
    m_Confidences.clear();
    m_Distances.clear();
};

void QuantBirdseed::setClusterOutFile(const std::string& clusterOutPath) {
  m_ClusterOutFile = clusterOutPath;
  m_ClusterOstrm.reset(new ofstream(clusterOutPath.c_str()));

}

bool QuantBirdseed::prepare(const IntensityMart &iMart) {
    bool result = true;
    for(int i=0; i<m_Reporters.size(); i++)
        result = result && m_Reporters[i]->prepare(*this, iMart);

  if (m_ClusterOstrm.get() != NULL) {
    (*m_ClusterOstrm) << "#%affymetrix-algorithm-param-" << "probeset_count" << "=" << m_Info.m_NumProbeSets << endl;
    std::vector<std::string>::const_iterator keyIx, paramIx;
    // The birdseed code from broad doesn't know about our headers so just put them at top of file
    for(keyIx = m_Info.m_ParamNames.begin(), paramIx = m_Info.m_ParamValues.begin();
      keyIx != m_Info.m_ParamNames.end() && paramIx != m_Info.m_ParamValues.end();
	++keyIx, ++paramIx) {
      // @todo should we using the AffymetrixParameterConsts.h #defined values?
      (*m_ClusterOstrm) << "#%affymetrix-algorithm-param-" << *keyIx << "=" << *paramIx << endl;
    }
    for(keyIx = m_Info.m_ClientInfoNames.begin(), paramIx = m_Info.m_ClientInfoValues.begin();
	keyIx != m_Info.m_ClientInfoNames.end() && paramIx != m_Info.m_ClientInfoValues.end();
      ++keyIx, ++paramIx) {
      // @todo should we using the AffymetrixParameterConsts.h #defined values?
      (*m_ClusterOstrm) << "#%affymetrix-application-meta-data-info-" << *keyIx << "=" << *paramIx << endl;
      
    }
  }
  return result;
}

// Copied from QuantBRLMM
/**
 * @brief Set up the quantification method given all the data about the probe
 * set, chip layout and data.
 *
 * @param psGroup - Probes to be used for final estimate.
 * @param layout - Chip layout annotation.
 * @param iMart - Raw data from chips.
 * @param iTrans - Transformations to be applied to data before use.
 * @param pmAdjust - How to estimate background, or MM probe.
 * @return True if setup sucessful, false otherwise.
 */
bool QuantBirdseed::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
                          std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust)
{
#if 0
    // Not doing this anymore.  If correction factor is negative, a SNP-specific correction
    // factor is calculated by birdseed.
    if (m_CorrectionFactor <= 0) {
        m_CorrectionFactor = averageIntensities(iMart);
    }
#endif
    const ProbeSet *gtPs = NULL; // Genotyping probeset
    blankSelf(); // Make sure we clear out the past analysis before starting this one.
    bool success = true;
    /* Sanity checks about probesets. */
    if(psGroup.probeSets.empty())
        Err::errAbort("Zero probesets in ProbeSetGroup (group: " + ToStr(psGroup.name) + ").");
    /* Remember this probeset. */
    gtPs = psGroup.probeSets[0];
    m_GtProbeSet = gtPs;
    m_ProbesetName = gtPs->name;
    if (!canSetUpProbeSet(gtPs)) {return false;}
    if(psGroup.probeSets.size() > 1)
        Err::errAbort("Can't have multiple probesets in a genotyping ProbeSetGroup (group: " + ToStr(psGroup.name) + ").");
 
    if(m_QuantMethod == NULL) {
        QuantMethodFactory factory(QuantMethodFactory::Expression);
        string spec = "plier.optmethod=1";
        ChipLayout layout;
        m_QuantMethod = factory.quantExprMethodForString(spec, layout, QuantMethodFactory::Expression);
    }

    fillInAlleleProbeSets(*gtPs, m_Aallele, m_Ballele);

    /* Get summaries for each allele. */
    success &= summarizeAllele(&m_Aallele, m_AValues, iMart, iTrans, pmAdjust, m_QuantMethod, false, true, m_Reporters);
	if (gtPs->psType == ProbeSet::Copynumber) {blankSelf(); return false;}
    success &= summarizeAllele(&m_Ballele, m_BValues, iMart, iTrans, pmAdjust, m_QuantMethod, false, true, m_Reporters);
    if(!success)
        blankSelf();
    return success;
}

void QuantBirdseed::setParameters(PsBoard &board) {
  Options *o = board.getOptions();
  vector<string> celFiles = o->getOptVector("cels");
  string outDir = o->getOpt("out-dir");
  GenotypeInfo *gInfo = board.getGenotypeInfo();
  if(gInfo == NULL) {
    Err::errAbort("QuantBirdseed::setParameters() - GenotypeInfo is NULL");
  }
  setGenders(gInfo->m_Genders->getGenders());
  if (o->getOpt("read-models-birdseed") != "") {
      string modelFile = o->getOpt("read-models-birdseed");
      m_PriorsReader.reset(new birdseed::dev::PriorsReader(*gInfo->m_SpecialSnps,
                                                           new birdseed::dev::TsvPriorsReader(modelFile, NULL)));
  }

  // @todo refactor - put current analysis path name on bboard
  string outfile = Fs::join(outDir,getType() + ".snp-posteriors");
  bool writeModels = o->getOptBool("write-models");
  if(writeModels) {
      setClusterOutFile(outfile);
      string modelFile = outfile + ".txt";
      board.set("model-file-written", modelFile);
  }
  else if (o->getOptBool("a5-write-models")) {
      ///@todo implement models in A5 for birdseed
      Err::errAbort("--a5-write-models for birdseed not implemented");
  }
  // Add quantification method.
  QuantMethodFactory qFactory;
  string qMethodSpec = o->getOpt("qmethod-spec");
  QuantMethod *qMethod = qFactory.quantMethodForString(qMethodSpec, board);

  // @todo - handle feature effects!
  QuantExprMethod *eMethod = NULL;
  eMethod= static_cast<QuantExprMethod *>(qMethod);
  if (eMethod == NULL) {
    Err::errAbort("Couldn't make a QuantExprMethod for specification: " + qMethodSpec);
  }
  setQuantExprMethod(eMethod);


}
