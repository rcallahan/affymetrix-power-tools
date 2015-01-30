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
 * @file   QuantLabelZ.cpp
 * @author Earl Hubbell
 * @date   Fri Feb 24 11:02:09 2006
 */

//
#include "chipstream/ArtifactReduction.h"
#include "chipstream/QuantLabelZ.h"
//
#include "chipstream/BioTypes.h"
#include "chipstream/MatrixUtil.h"
#include "chipstream/QuantLabelZIO.h"
#include "chipstream/QuantMethodFactory.h"
#include "chipstream/GenotypeInfo.h"
#include "util/Fs.h"

using namespace affx;

//////////

/**
 *  Utility function for returning subset values to place in whole sample list
 *
 * @param Target - destination for values
 * @param src - values defined for subset
 * @param subset - definition of subset of samples
 */
template <typename T1>
void t_place_subset(vector<T1> &Target, vector<T1> &src, vector<int> &subset)
{
  int ti,i;
  ti=0;
  for (i=0; i<subset.size(); i++)
  if (subset[i]>0) {
    Target[i] = src[ti];
    ti++;
  }
}

/**
 *  Utility function to grab a set of samples (like male/female)
 *
 * @param target - destination for values
 * @param source - values defined for subset
 * @param subset - definition of subset of samples
 */
template <typename T1>
void t_extract_subset(std::vector<T1> &target, std::vector<T1> &source, std::vector<int> &subset)
{
  int i,ti;

  ti=0;
  for (i=0; i<source.size(); i++) {
    if (subset[i]>0) {
      target[ti] = source[i];
      ti++;
    }
  }
}

//////////

//////////

QuantLabelZ::QuantLabelZ(enum Transformation transform,
                         double K,
                         bool lowPrecision) 
{
  if (transform > CCS || transform < MvA) {
    Err::errAbort("QuantLabelZ() - Expecting transform to be between: " + ToStr(MvA) + " (MvA) and: " + 
                  ToStr(CCS) + " (CCS), but got: " + ToStr(transform));
  }
  m_bCopyNumber = false;
  m_Param.m_LowPrecision = lowPrecision;
  m_Param.m_Transform = transform;
  m_Param.m_K = K;
  em_thresh = 0.05;
  em_cutoff = 0.5;
  gender_cutoff = 0.1;
  m_QuantMethod = NULL;
  m_SelectionOverride = 0; // do not override values for selection
  sp.Initialize();
  m_KnownGenoTypes = NULL;
  m_SnpPosteriorTsv_ver=-1;
  m_ProbeSetsToReport = NULL;
  m_ZWGenderCalling = false;

  setupSelfDoc(*this);
  setOptValue("transform", stringForTransformation(m_Param.m_Transform));
  setOptValue("K", ToStr(m_Param.m_K));
  setOptValue("lowprecision", m_Param.m_LowPrecision);
  setOptValue("copynumber", m_bCopyNumber);
  // set some bogus filenames for debugging.
  m_NormSummaryTsv.setFilename("m_NormSummaryTsv-DEBUG");
  m_SnpPosteriorTsv.setFilename("m_SnpPosteriorTsv-DEBUG");
}
  
/** Destructor. */
QuantLabelZ::~QuantLabelZ()
{
  if (m_SequentialModelTsv.is_open()) 
    m_SequentialModelTsv.close();
  if (m_NormSummaryTsv.is_open())
    m_NormSummaryTsv.close();
  if (m_SnpPosteriorTsv.is_open())
    m_SnpPosteriorTsv.close();
  if (m_SnpProbeTsv.is_open())
    m_SnpProbeTsv.close();
  //
  for (std::vector<QuantMethodReport *>::iterator it = m_Reporters.begin(); it != m_Reporters.end(); it++)
    {
      delete *it;
    }
  m_Reporters.clear();
  m_vectSnpPriors.deleteAll();
  if (m_QuantMethod != NULL) {delete m_QuantMethod; m_QuantMethod = NULL;}
}

bool QuantLabelZ::prepare(const IntensityMart &iMart)
{
  bool result = true;
  for (int i=0; i<m_Reporters.size(); i++)
    result = result && m_Reporters[i]->prepare(*this, iMart);

  //
  if (!m_SnpPosteriorName.empty()) {
    std::vector<std::string>::const_iterator keyIx, paramIx;
    //
    Verbose::out(1, "Got " + ToStr(m_Info.m_ParamValues.size()) + " headers to write.");
    m_SnpPosteriorTsv.addHeader("affymetrix-algorithm-param-probeset_count", ToStr(m_Info.m_NumProbeSets));
    for (keyIx = m_Info.m_ParamNames.begin(), paramIx = m_Info.m_ParamValues.begin();
	keyIx != m_Info.m_ParamNames.end() && paramIx != m_Info.m_ParamValues.end();
	++keyIx, ++paramIx) {
      // @todo should we using the AffymetrixParameterConsts.h #defined values?
      m_SnpPosteriorTsv.addHeader("affymetrix-algorithm-param-" + *keyIx, *paramIx);
    }
    for (keyIx = m_Info.m_ClientInfoNames.begin(), paramIx = m_Info.m_ClientInfoValues.begin();
	keyIx != m_Info.m_ClientInfoNames.end() && paramIx != m_Info.m_ClientInfoValues.end();
      ++keyIx, ++paramIx) {
      // @todo should we using the AffymetrixParameterConsts.h #defined values?
      m_SnpPosteriorTsv.addHeader("affymetrix-application-meta-data-info-" + *keyIx, *paramIx);
    }
    QuantLabelZ__writeSnpPosteriorTsv(m_SnpPosteriorTsv,m_SnpPosteriorName,m_PosteriorFormat,m_SnpPosteriorTsv_ver); 
  }
  return result;
}

bool QuantLabelZ::finish() 
{
  bool result = true;
  m_SnpPosteriorTsv.close();
  for (int i=0; i<m_Reporters.size(); i++)
    result = result && m_Reporters[i]->finish(*this);
  return result;
}

/** Fill in our self documentation from current state. */
void QuantLabelZ::setupSelfDoc(SelfDoc &doc) {
  doc.setDocName(QUANTBRLMMP);
  doc.setDocDescription("Do genotyping calls with BRLMM-P (perfect match) algorithm.");
  doc.setDocOptions(getDefaultDocOptions());
}

void QuantLabelZ::blankSelf() 
{
  clearProbeSet(m_Aallele);
  clearProbeSet(m_Ballele);
  m_AValues.clear();
  m_BValues.clear();
  m_Calls.clear();
  m_Confidences.clear();
  m_Distances.clear();
  m_InitialCalls.clear();
}

//////////

/**
 * Sets up the intensities and summaries for a probe set
 * Separate function, because shared between several routines (gender, normalization, labeling)
 * Ideally, would be inherited from global genotyping class
 *
 * @param gtPs - probe set to be set up
 * @param layout - probes located on chip
 * @param iMart - serves up raw intensities
 * @param iTrans - transforms intensities
 * @param pmAdjust - background adjust
 * @param doReport - throw summaries out to reporters or not
 */
bool QuantLabelZ::SetUpProbeSet(const ProbeSet *gtPs,
                                const IntensityMart &iMart, std::vector<ChipStream *> &iTrans,
                                PmAdjuster &pmAdjust, bool doReport){

  bool success=true;

  if (!canSetUpProbeSet(gtPs)) {return false;}
  m_ProbesetName = gtPs->name;
  m_GtProbeSet = gtPs;
  m_AValues.clear();
  m_BValues.clear();

  bool report = doReport;
  if (m_ProbeSetsToReport != NULL) {
    if (m_ProbeSetsToReport->find(m_ProbesetName.c_str())==m_ProbeSetsToReport->end()) {
      report = false;
    }
  }
  fillInAlleleProbeSets(*gtPs, m_Aallele, m_Ballele);
  // if by some mischance we have failed to specify method
  if (m_QuantMethod == NULL) {
    QuantMethodFactory factory(QuantMethodFactory::Expression);
    string spec = "plier.optmethod=1";
    ChipLayout layout;
    m_QuantMethod = factory.quantExprMethodForString(spec, layout, QuantMethodFactory::Expression);
  }
  // load up the data
  success &= summarizeAllele(&m_Aallele, m_AValues, iMart, iTrans, pmAdjust, m_QuantMethod, m_Param.m_LowPrecision, report,m_Reporters);
  if (gtPs->psType == ProbeSet::Copynumber) {return false;}
  if (m_Ballele.atoms.size() < 1) {
    return false;
  }
  success &= summarizeAllele(&m_Ballele, m_BValues, iMart, iTrans, pmAdjust, m_QuantMethod, m_Param.m_LowPrecision, report,m_Reporters);

  // grab probe values if doing probe selection
  if (m_SnpProbeTsv.is_open()){
    //do this again, because summarizeAllele clears the probeset m_Aallele and m_Ballele
    fillInAlleleProbeSets(*gtPs, m_Aallele, m_Ballele);
    extractProbes(m_PM_A, &m_Aallele, m_A_probeId, iMart,iTrans,pmAdjust,m_Param.m_LowPrecision,m_QuantMethod);
    extractProbes(m_PM_B, &m_Ballele, m_B_probeId, iMart,iTrans,pmAdjust,m_Param.m_LowPrecision,m_QuantMethod);
  }

  // transform the data as specified
  transform();
  loadInitialCallsForSnp(gtPs->name);
  if ( m_bProbeSetTrustOpt && m_ProbeSetTrust.empty() ) {
    loadProbeSetTrustFromTmpFile();
  }
  return(success);
}

/**
 * Load up the initial calls for SNP
 */
void QuantLabelZ::loadInitialCallsForSnp(const string &probeset) {
  // setup hints for genotype
  m_InitialCalls.clear();
  if (m_KnownGenoTypes != NULL) {
    KnownGenoTypes_iter_t mapIter = m_KnownGenoTypes->find(probeset.c_str());
    //  std::map<string, vector<affx::GType> >::iterator mapIter = m_KnownGenoTypes.find(gtPs->name);
    if (mapIter == m_KnownGenoTypes->end()) {
        // unlike BRLMM, this is not >necessary< for analysis
        Verbose::out(4,"Did not find known genotype for probeset " + probeset);
    }
    else{
        vector<affx::GType> &gType = mapIter->second;
        for (unsigned int i = 0; i < gType.size(); i++) {
        m_InitialCalls.push_back(gType[i]);
        }
    }
  }
}

void QuantLabelZ::loadProbeSetTrustFromTmpFile(){
  affx::TsvFile tsv;
  std::string chpId;
  std::string probesetId;
  int trustedCount;
  
  // Aborts on error. 
  if ( tsv.open(m_ProbeSetTrustTmpFileName) != affx::TSV_OK ) {
    APT_ERR_ABORT(m_ProbeSetTrustTmpFileName + " file not found or permission denied");
  }

  tsv.bind(0, "chp_id", &chpId,   affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "probeset_id", &probesetId,   affx::TSV_BIND_REQUIRED);
  tsv.bind(0, "trust_count", &trustedCount, affx::TSV_BIND_REQUIRED);

  m_ProbeSetTrust.clear();
  
  while (tsv.nextLevel(0)==TSV_OK) {
    std::string trustKey = chpId + ":" + probesetId;
    APT_ERR_ASSERT( m_ProbeSetTrust.count(trustKey) == 0,
                    "Invalid artifact reduction trusted probeset file, "
                    + m_ProbeSetTrustTmpFileName
                    + ", has duplicate id: "
                    + trustKey);
    
    m_ProbeSetTrust[trustKey] = trustedCount;
  }

  tsv.close();
}

/**
 * SelectGoodProbes selects probe pairs that contribute most to classification
 */
void QuantLabelZ::SelectGoodProbes() {
  // rescan the data to find the probes that contributed most to the output
  // compare genotypes to input data

  // step 1: obtain all the raw probe input data from ExprMethod
  // step 1 is actually done in SetUpProbeSet
  // step 2: Iterate over all pairs of probes
  // step 2a: construct transformed data for each pair of probes
  // setp 2b: compute genotyping regression/classification against data
  // step 2c: evaluate goodness of this probe pair
  // step 3: dump best probe(s) to summary file with evaluation
  // step 3a: dump all probes to summary file with evaluation

  vector<double> DummyA, DummyB,Predictor;
  vector< vector<double> > Covariate;
  vector<int> RealCalls,TmpCalls;
  vector<double> Summaries;
  double AIC;
  unsigned int i,j,k,ValidCallNum;
  vector<double> PriorM;
  vector<double> PriorW;

  PriorM.resize(3);
  PriorW.resize(3);
  // use the fitted posterior
  PriorM[0] = sp.posterior.aa.m;
  PriorM[1] = sp.posterior.ab.m;
  PriorM[2] = sp.posterior.bb.m;
  PriorW[0]=PriorW[1]=PriorW[2] = 0.01; // 1 percent data point if no obs

  RealCalls.resize(m_Calls.size());
  for (i=0; i<RealCalls.size(); i++) {
     // make sure this is "int" worthy and 0-2 range
      RealCalls[i] = affx::GType_to_int(m_Calls[i]);
   }
  for (i=0; i<m_InitialCalls.size(); i++) {
      if (m_InitialCalls[i]!=affx::NN && m_SelectionOverride)
        // make sure this is "int" worthy and 0-2 range
        RealCalls[i] = affx::GType_to_int(m_InitialCalls[i]);
   }

  // now I have a set of calls, some of which may be nonexistent
  // filter the data by valid calls
  ValidCallNum=0;
  for (i=0; i<RealCalls.size(); i++) {
    if (RealCalls[i]>-1) {
        ValidCallNum++;
    }
  }
  if (ValidCallNum < 2) {
    Err::errAbort("At least two valid calls are needed for each probeset in order to use --select-probes.  Probeset " + m_ProbesetName + " has one or zero valid calls.");
  }
  TmpCalls.resize(ValidCallNum); // how many calls Valid
  k=0;
  for (i=0; i<RealCalls.size(); i++)
    if (RealCalls[i]>-1){
        TmpCalls[k] = RealCalls[i];
        k++;
    }

  DummyA.resize(ValidCallNum);
  DummyB.resize(ValidCallNum);

  // loop over all probes
  for (j=0; j<m_PM_A.size() && (ValidCallNum>0); j++) {
    k=0;
    for (i=0; i<m_PM_A[j].size(); i++) {
        // filter the data by valid calls
        if (RealCalls[i]>-1) {
            DummyA[k] = m_PM_A[j][i];
            DummyB[k] = m_PM_B[j][i];
            k++;
        }
    }
    GenoUtility_transformData(DummyA, DummyB, m_Param.m_Transform,m_Param.m_K);
    Covariate.clear();
    //
    // take out 'no calls'
    Covariate.push_back(DummyA);
    Summaries.clear();
    //cout << endl;
    //cout << m_ProbesetName << "\t"; // dummy test
    AIC = LogisticRegression(Summaries,Predictor,TmpCalls,Covariate,PriorM,PriorW,100,.001);

    writeSnpProbeValue(m_SnpProbeTsv,j,m_A_probeId[j],m_B_probeId[j],
                        Summaries,Predictor,Predictor,m_ProbesetName);
  }
}

/*
 * @brief normalizes the data or other modifications
 */
void QuantLabelZ::PreprocessValues() {
  // transform data, do any normalization
  // transform the data
  // transform()
  // generate calls for tranformed data

  // here, we apply normalization prediction across everything

  // AValues = Contrast
  // BValues = Strength
  vector<double> Cov; // to be covariates, offset/cubed as required
  double tmpf, strength, contrast;
  unsigned int sIx; // sample iterator
  for (sIx=0; sIx<normMap.size(); sIx++) {
    // isolate this, because not all transforms alike
    strength = m_BValues[sIx];
    contrast = m_AValues[sIx];
    // make offset/cubed covariates
    normMap[sIx].MakeCov(Cov,strength);
    // generate prediction
    tmpf = normMap[sIx].MakePrediction(contrast,Cov);
    // place back in appropriate register
    m_AValues[sIx]=tmpf;
  }

  // write out processed values
  if (m_NormSummaryTsv.is_open()) {
    writeNormSummaryValue(m_ProbesetName,m_AValues,m_BValues);
  }

}

void QuantLabelZ::EliminateUntrustedSamples(vector<int> & SnpSampleCopyNumber){
  // remove any sample who doesn't have trustworthy probes
  // by setting "copynumber" =0 in such cases
  for (int sIx=0; sIx<m_AValues.size(); sIx++){
    if (!isTrusted(sIx)){
      SnpSampleCopyNumber[sIx]=0;
    }
  }
}

// just checks to see if this sample has any trusted probes
bool QuantLabelZ::isTrusted(int cIx){
  if ( m_ProbeSetTrust.empty() ) {
    return true;
  }
  map<string,int>::iterator iter;
  bool retval=true;  // return trusted if not found
  
  // only get zero block for testing
  std::string trustKey = ToStr(cIx)+":"+m_ProbesetName; 
  iter = m_ProbeSetTrust.find(trustKey);
  if (iter!=m_ProbeSetTrust.end() &&  (iter->second<1) ) {
    retval=false;
  }
  return(retval);
}

// void QuantLabelZ::getPrior(string &TmpName, int copyIx, snp_labeled_distribution &objSearch, snp_param &tsp) {
//   // first search for name+copy number
//   objSearch.probeset_id = TmpName;
//   int iSearchIndex = m_vectSnpPriors.binarySearch(objSearch, 0);
//   if (iSearchIndex == -1)
//     {
//       // first fail to generic - avoid command line needed!
//       if (copyIx<2)
//         {
//           objSearch.probeset_id = "GENERIC:" + ToStr(copyIx);
//               }
//       else
//         objSearch.probeset_id = "GENERIC";
//       iSearchIndex = m_vectSnpPriors.binarySearch(objSearch,0);
//       if (iSearchIndex==-1)
//         {
          
//           // Tell the user we tried to search for GENERIC as well but didn't find either
//           Err::errAbort("Can't find model for SNP: " + TmpName + " or " + objSearch.probeset_id);
//         }
//     }
//   // if success, grab & copy
//   if (iSearchIndex>-1)
//     {
//       snp_labeled_distribution* p = m_vectSnpPriors.getAt(iSearchIndex);
//       // this is really one copy operation!!!
//       tsp.prior.Copy(p->Dist);
//       /*
//         tsp.prior.aa = p->Dist.aa;
//         tsp.prior.ab = p->Dist.ab;
//         tsp.prior.bb = p->Dist.bb;
//         tsp.prior.xah = p->Dist.xah;
//         tsp.prior.xab = p->Dist.xab;
//         tsp.prior.xhb = p->Dist.xhb;
//         tsp.prior.yah = p->Dist.yah;
//         tsp.prior.yab = p->Dist.yab;
//         tsp.prior.yhb = p->Dist.yhb;
//         tsp.prior.xyah = p->Dist.xyah;
//         tsp.prior.xyab = p->Dist.xyab;
//         tsp.prior.xyhb = p->Dist.xyhb;
//         tsp.prior.yxah = p->Dist.yxah;
//         tsp.prior.yxab = p->Dist.yxab;
//         tsp.prior.yxhb = p->Dist.yxhb;
//       */
//     }
// }

void QuantLabelZ::getSequentialPrior(std::string &name, snp_param &tsp, snp_labeled_distribution &sDist) {
  if (m_SequentialModelTsv.nextLevel(0) == affx::TSV_OK) {
    QuantLabelZ__SnpPriorFromStrings(sDist, m_BB, m_AB, m_AA, m_CV);
    sDist.probeset_id = m_Id;
    tsp.prior.Copy(sDist.Dist);
    if (name != m_Id) {
      APT_ERR_ABORT("Expecting: " + name + " but got: " + m_Id);
    }
  }
  else {
    APT_ERR_ABORT("Expecting model for snp: " + name);
  }
}

void QuantLabelZ::getPrior(std::string &TmpName, int presentCopyNumber, snp_param &tsp, snp_labeled_distribution &sDist) {
  if (m_SequentialModelTsv.is_open()) {
    getSequentialPrior(TmpName, tsp, sDist);
  }
  else if (m_vectSnpPriors.getCount() > 0) {
    snp_labeled_distribution objSearch;
    snp_labeled_distribution* p=NULL; 
    objSearch.probeset_id = TmpName;
    int iSearchIndex = m_vectSnpPriors.binarySearch(objSearch, 0);
    if (iSearchIndex == -1){
      if (presentCopyNumber<2){
        objSearch.probeset_id = "GENERIC:" + ToStr(presentCopyNumber);
      } else {
        objSearch.probeset_id = "GENERIC";
      }
      iSearchIndex = m_vectSnpPriors.binarySearch(objSearch,0);
    } 
    if (iSearchIndex==-1) {
      // should this fail to global prior?
      Err::errAbort("Can't find model for SNP: " + TmpName);
    } else {
      p = m_vectSnpPriors.getAt(iSearchIndex);
      sDist = *p;
      // this is really one copy operation!!!
      tsp.prior.Copy(p->Dist);
      /*
        tsp.prior.aa = p->Dist.aa;
        tsp.prior.ab = p->Dist.ab;
        tsp.prior.bb = p->Dist.bb;
        tsp.prior.xah = p->Dist.xah;
        tsp.prior.xab = p->Dist.xab;
        tsp.prior.xhb = p->Dist.xhb;
        tsp.prior.yah = p->Dist.yah;
        tsp.prior.yab = p->Dist.yab;
        tsp.prior.yhb = p->Dist.yhb;
        tsp.prior.xyah = p->Dist.xyah;
        tsp.prior.xyab = p->Dist.xyab;
        tsp.prior.xyhb = p->Dist.xyhb;
        tsp.prior.yxah = p->Dist.yxah;
        tsp.prior.yxab = p->Dist.yxab;
        tsp.prior.yxhb = p->Dist.yxhb;
      */
    }
  }
}

std::string QuantLabelZ::getCallProbabilities(int idx) {
  return ToStr(m_Txb[idx])+","+ToStr(m_Txh[idx])+","+ToStr(m_Txa[idx])+","+ToStr(1 - m_Txb[idx] - m_Txh[idx] - m_Txa[idx]);
}

void QuantLabelZ::SubsetSetup(vector<int> &samplesubset, int &subsetlength, int copyIx, vector<int>& SnpSampleCopyNumber) {
  int sIx;

   for (sIx=0; sIx<samplesubset.size(); sIx++) {
      //logic
      if (SnpSampleCopyNumber[sIx]==copyIx)
        samplesubset[sIx]=1;
      else
        samplesubset[sIx]=0;
    }

    // setup input data
    subsetlength = 0;
    for (int i=0; i<samplesubset.size(); i++)
      subsetlength += samplesubset[i];
}

/**
 * @brief Do the heavy lifting of estimation.
 */
void QuantLabelZ::computeEstimate() {


  static int call_count = 0;
  call_count++;
  
  snp_labeled_distribution objSearch;
  // first, do anything necessary
  PreprocessValues();

  vector<affx::GType> tcall;
  vector<double> tconf, tx, ty;
  vector<double> txa, txh, txb;

  // set up snp specific prior if available
  snp_param tsp;

  // numeric method does its magic here
  // call the numerics
  int length;
  length = m_AValues.size();
  m_Calls.assign(length,affx::NN);
  m_Confidences.assign(length,0);
  if (m_OutputProbabilities) {
    m_Txa.assign(length, 0);
    m_Txh.assign(length, 0);
    m_Txb.assign(length, 0);
  }

  // standardized distance to each cluster
  m_Distances.resize(3);
  for (int m_Di=0; m_Di<3; m_Di++)
    m_Distances[m_Di].assign(length,0);

  vector<int> samplesubset;
  vector<int> GenoHint, GenoRef;
  vector<double> SubsetInbredPenalty; 
  
  samplesubset.assign(length,1);
  int subsetlength;
  int copyIx,sIx;
  string TmpName;

  // put together the reference vector of seeds, if it exists
  GenoRef.assign(length,affx::NN);
  for (unsigned int qIx=0; qIx<m_InitialCalls.size(); qIx++)
    GenoRef[qIx] = affx::GType_to_int(m_InitialCalls[qIx]);

  // look up the right copy number for this snp
  int MaleCopy,FemaleCopy;
  map<string,pair<int, int> >::iterator test;
  test=m_SpecialSnps.find(m_ProbesetName);
  if (test==m_SpecialSnps.end()) {
    MaleCopy=2;
    FemaleCopy=2;
  }
  else {
    MaleCopy=test->second.first;
    FemaleCopy=test->second.second;
  }

  // assign sample copy numbers for this snp
  vector<int> SnpSampleCopyNumber;
  SnpSampleCopyNumber.assign(length,2); // default for all regular snps
  for (sIx=0; sIx<length; sIx++) {
    if (m_Genders[sIx] == affx::Female)
      SnpSampleCopyNumber[sIx] = FemaleCopy;
    else if (m_Genders[sIx] == affx::Male)
      SnpSampleCopyNumber[sIx] = MaleCopy;
    else if (m_Genders[sIx] == affx::UnknownGender) {
      if (m_ZWGenderCalling) {
        SnpSampleCopyNumber[sIx] = MaleCopy;
      }
      else {
        SnpSampleCopyNumber[sIx] = FemaleCopy;
      }
    }
    else
      Err::errAbort("Unknown gender seen by genotype caller: " + ToStr(m_Genders[sIx]));
  }

  map<string,vector<int> >::iterator test2;
  test2=m_SpecialSampleSnps.find(m_ProbesetName);
  if (test2!=m_SpecialSampleSnps.end()) {
    for (sIx=0; sIx<length; sIx++) {
      SnpSampleCopyNumber[sIx] = test2->second[sIx];
    }
  }

  if (m_bProbeSetTrustOpt)
    EliminateUntrustedSamples(SnpSampleCopyNumber);

  // Do subsets by copy number
  for (copyIx=2; copyIx>0; copyIx--) {
    
    // assign appropriate subset to current copy number
    SubsetSetup(samplesubset,subsetlength, copyIx, SnpSampleCopyNumber);

    // if there's anything to do
    if (subsetlength>0) {
      // if there's something to do, find the SNP specific prior

      tsp.copy(sp); // duplicate current parameters to keep commands
      tsp.copynumber=copyIx; // what copy for this snp so know right model centers

      // look up this snp and/or prepare for output
      TmpName = m_ProbesetName;
      if (copyIx<2)
        TmpName = m_ProbesetName + ":" +ToStr(copyIx);
          // if any to look up
      snp_labeled_distribution sDist;
      if (m_vectSnpPriors.getCount() > 0 || m_SequentialModelTsv.is_open()) {
        getPrior(TmpName, copyIx, tsp, sDist);
      }

      tcall.resize(subsetlength);
      tconf.resize(subsetlength);
      tx.resize(subsetlength);
      ty.resize(subsetlength);
      GenoHint.assign(subsetlength,-1);
      SubsetInbredPenalty.assign(subsetlength,0); // Null penalty

      t_extract_subset(tx,m_AValues,samplesubset);
      t_extract_subset(ty,m_BValues,samplesubset);
      t_extract_subset(GenoHint,GenoRef,samplesubset);
      t_extract_subset(SubsetInbredPenalty,m_InbredHetPenalty,samplesubset);
      // am I getting anything useful here
     
      // note: using temporary copy of main parameters tsp
      bayes_label(tcall,tconf,tx,ty,GenoHint,SubsetInbredPenalty,tsp,txa,txh,txb,m_OutputProbabilities);

      // transfer data to calls
      t_place_subset(m_Calls,tcall,samplesubset);
      t_place_subset(m_Confidences,tconf,samplesubset);
      if (m_OutputProbabilities) {
        t_place_subset(m_Txa,txa,samplesubset);
        t_place_subset(m_Txh,txh,samplesubset);
        t_place_subset(m_Txb,txb,samplesubset);
      }

      // note: bayes_label uses probability, not distance
      // computing standardized distance to be compatible QC

      vector<double> tsd;
      int tsdi;
      tsd.assign(subsetlength,0);
      for (tsdi=0; tsdi<subsetlength; tsdi++) {
        tsd[tsdi] = abs((tsp.posterior.aa.m-tx[tsdi])/sqrt(tsp.posterior.aa.ss));
      }
      t_place_subset(m_Distances[0],tsd,samplesubset);
      for (tsdi=0; tsdi<subsetlength; tsdi++) {
        tsd[tsdi] = abs((tsp.posterior.ab.m-tx[tsdi])/sqrt(tsp.posterior.ab.ss));
      }
      t_place_subset(m_Distances[1],tsd,samplesubset);
      for (tsdi=0; tsdi<subsetlength; tsdi++) {
        tsd[tsdi] = abs((tsp.posterior.bb.m-tx[tsdi])/sqrt(tsp.posterior.bb.ss));
      }
      t_place_subset(m_Distances[2],tsd,samplesubset);

      // now possible output of snp-specific posterior
      // that is, now we've updated the prior with the data
      // and we want to save what we learned
      if (m_SnpPosteriorTsv.is_open()) {
        writeSnpPosteriorValue(TmpName, tsp);
      }

    } // done processing this copy number
  } // back into the loop for more copy number tries

  // made calls for all samples
  // select probe that best mimics
  // hack: chrX not yet handled differently
  // but should work in regression framework?

  if (m_SnpProbeTsv.is_open()){
    SelectGoodProbes();
  }
  
}

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
bool QuantLabelZ::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
                        std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust) {
	const ProbeSet *gtPs = NULL; //		Genotyping probeset
  blankSelf(); // Make sure we clear out the past analysis before starting this one.
  bool success = true;
  /* Sanity checks about probesets. */
  if (psGroup.probeSets.empty())
    Err::errAbort("Zero probesets in ProbeSetGroup (group: " + ToStr(psGroup.name) + ").");
  if (psGroup.probeSets.size() > 1)
    Err::errAbort("Can't have multiple probesets in a genotyping ProbeSetGroup (group: " + ToStr(psGroup.name) + ").");

  /* Remember this probeset. */
  gtPs = psGroup.probeSets[0];
  success &= SetUpProbeSet(gtPs, iMart,iTrans,pmAdjust, true);
  if (!success)
    blankSelf();
  return success;
}

/**
 * Construct a contrast normalization function from a collection of probesets. *
 * @param probeSets - ProbeSets to be used in computing normalizer.
 * @param layout - description and chip and probesets.
 * @param iMart - Cel file data for probesets.
 * @param iTrans - Transformers to operate on raw data (i.e. normalization and bg subtraction).
 * @param pmAdjust - Individual probe adjustment like pm-mm or pm-only.
 *
 * @return - returns a normalization function set.
 */
vector<NormalizationPredictor> QuantLabelZ::makeNorm(std::vector<std::string> probeSets, ChipLayout &layout,
                                   const IntensityMart &iMart, std::vector<ChipStream *> &iTrans,
                                   PmAdjuster &pmAdjust) {
  //normMap Qmap;
  vector<vector<double> > SampleContrast;
  vector<vector<double> > SampleStrength;
  bool success=true;
  unsigned int i;
  unsigned int nSample; // how do I detect this again?
  nSample=0;

  Verbose::out(1,"making normalization map");
  if (probeSets.empty()) {
     Err::errAbort("bogus snp sample");
  }

  // for each probeset in my sample
  for (unsigned int psIx = 0; psIx < probeSets.size(); psIx++) {
    // Have to skip probesets that don't have 2 or 4 groups as it
    // appears that some are mislabelled in the cdf file.
    ProbeListPacked pList = layout.getProbeListByName(probeSets[psIx]);
    ///@todo deal with missing probes from kill list
    if (pList.isNull()) {
      Err::errAbort("QuantLabelZ::makeNorm(): Cannot deal with missing probes from kill list");
    }
    ProbeSet *gtPs = ProbeListFactory::asProbeSet(pList);
    success = SetUpProbeSet(gtPs, iMart,iTrans,pmAdjust, false);
    if (!success) {
      delete gtPs;
      continue;
    }
    // add snpStats values
    // currently contrast/strength
    // push into the slots by sample
    // set up the slots to put everything in
    if (nSample<1) {
      nSample=m_AValues.size();
      SampleContrast.resize(nSample);
      SampleStrength.resize(nSample);
    }
    for (i=0; i<nSample; i++) {
      SampleContrast[i].push_back(m_AValues[i]);
      SampleStrength[i].push_back(m_BValues[i]);
    }
    // look up parameters
    //cout << "SNP:\t" << psIx << endl;
    delete gtPs;
  }
  NormalizationFrame TmpFrame;
  vector<NormalizationPredictor> Qmap;
  Qmap.resize(nSample);
  for (i=0; i<nSample; i++) {
    // for each sample
    TmpFrame.FillInContrast(SampleContrast[i]);
    TmpFrame.FillInCovariates(SampleStrength[i]);
    TmpFrame.FitEM();
    //store in Qmap
    Qmap[i].Copy(TmpFrame.ThisFit);
    /*cout << "Sample:\t" << i << endl;
    td(Qmap[i].fitAA);
    td(Qmap[i].fitAB);
    td(Qmap[i].fitBB);
    td(Qmap[i].sigma);*/
  }
  // compute normalization function here
  // fun for the whole family
  return Qmap;
}

/**
 * Calls gender using EM.
 * This should not be anywhere near LabelZ, but code structure & time make this easiest.
 *
 * @param contrast - list of contrast values for chrX snps in one sample
 * @param em_thresh -  when are we done?
 * @param em_cutoff - enough data to make decision?
 * @param gender_cutoff - make decision as to gender
 */
Gender t_CallEMGender(const std::vector<double>& contrast,
                      double& em_thresh,
                      double& em_cutoff,
                      double& gender_cutoff) {
  Gender retval;
  CEMSeed seed;
  seed.setMu(-0.66, 0, 0.66);
  seed.setSigma(0.1f, 0.1f, 0.1f);
  seed.setWeight(0.33f, 0.34f, 0.33f);
  seed.setMinMu(-2, -.05, .25);
  seed.setMinSigma(0.02f, 0.02f, 0.02f);
  seed.setMaxMu(-.25, .05, 2);
  seed.setMaxSigma(.3, .3, .3);

  CPrimeEM em;
  em.setData(contrast);
  em.setThreshold(em_thresh);
  em.EMEstimate(seed);
  CEMEst* pEst = em.getEMEstimates();
  //int nNC = count_if (pEst->m_class.begin(),pEst->m_class.end(),bind2nd(equal_to<int>(),-1));
  int nAA = count_if (pEst->m_class.begin(),pEst->m_class.end(),bind2nd(equal_to<int>(),0));
  int nBB = count_if (pEst->m_class.begin(),pEst->m_class.end(),bind2nd(equal_to<int>(),2));
  int nAB = count_if (pEst->m_class.begin(),pEst->m_class.end(),bind2nd(equal_to<int>(),1));
  double HCR=-1;
  if (nAA+nBB+nBB > 0 && double(nAA+nAB+nBB)/pEst->m_class.size() > em_cutoff) {
      HCR = double(nAB)/double(nAA+nAB+nBB);
      if (HCR < gender_cutoff)
        retval = affx::Male; // male
      else
        retval = affx::Female; // female
  }
  else
    retval = affx::UnknownGender; // unknown
  return retval;
}

/**
 * Find genders of sample - uses EM, not labeling as you might think!
 * Again, this function would be elsewhere if we had a good structure
 *
 * @param genders - to be determined - what are the genders of the sample
 * @param probeSets - use these snps to find genders
 * @param layout - where are probes on array
 * @param iMart - serve up raw intensities
 * @param iTrans - adjust intensities (normalize, etc)
 * @param pmAdjust - adjust intensities for background
 */
void QuantLabelZ::FindGender(std::vector<affx::Gender> &genders,
                             std::vector<ProbeSet *> probeSets,
                             ChipLayout &layout,
                             const IntensityMart &iMart,
                             std::vector<ChipStream *> &iTrans,
                             PmAdjuster &pmAdjust) {

  // probeSets contains X-chromosome snps
  // so we're going to grab all the summary values
  // and call some genders
  vector<vector<double> > SampleContrast;
  // dead code?
  //vector<vector<double> > SampleStrength;
  bool success=true;
  unsigned int i;
  unsigned int nSample; // how do I detect this again?
  nSample=0;

  //cout << "Preparing to call gender" << endl;
  if (!probeSets[0]) {
    Err::errAbort("bogus snp sample");
  }

  // for each probeset in my sample
  for (unsigned int psIx = 0; psIx < probeSets.size(); psIx++) {
    const ProbeSet *gtPs = probeSets[psIx];
    // Have to skip probesets that don't have 2 or 4 groups as it
    // appears that some are mislabelled in the cdf file.

    success = SetUpProbeSet(gtPs, iMart,iTrans,pmAdjust, false);
    if (!success)
      continue;
    // add snpStats values
    // currently contrast/strength
    // push into the slots by sample
    // set up the slots to put everything in
    //cout << "Set up pset: " << psIx << endl;
    if (nSample<1) {
      nSample=m_AValues.size();
      SampleContrast.resize(nSample);
      // dead code?
      // SampleStrength.resize(nSample);
    }
    for (i=0; i<nSample; i++) {
      SampleContrast[i].push_back(m_AValues[i]);
      // dead code? -jhg
      //SampleStrength[i].push_back(m_BValues[i]);
    }
  }

  // same here as normalization, now call gender
  // okay have setup the big matrix of contrasts
  // for each sample, call the gender
  //cout << "setup contrast: call gender " << nSample << "\t" << genders.size() << endl;
  std::vector<affx::Gender> Gmap;
  Gmap.resize(nSample);
  for (i=0; i<nSample; i++) {
    if (true) {
      Gmap[i] = t_CallEMGender(SampleContrast[i],em_thresh,em_cutoff,gender_cutoff);
      //cout << i << "\t"<< Gmap[i] << "\t" << genders[i] << endl;
      genders[i] = Gmap[i];
    }
  }
  //cout << endl;
  //return Gmap;
}

//////////

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> QuantLabelZ::getDefaultDocOptions() {
  std::vector<SelfDoc::Opt> opts;

  SelfDoc::Opt k = {"K", SelfDoc::Opt::Double, "4.0", "4.0", "0", "NA",
                    "Scale parameter used used in CCS and CES transformations."};
  opts.push_back(k);
  SelfDoc::Opt transform = {"transform", SelfDoc::Opt::String, "CCS", "CCS", "NA", "NA",
                            "Transformation of initial data are we feeding into the classifier? {'CCS', 'CES', 'MvA','RvT'}"};
  opts.push_back(transform);
  SelfDoc::Opt lowprecision = {"lowprecision", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
                               "R prototype uses summary values rounded to first decimal place. Use this flag to be simulate behavior."};
  opts.push_back(lowprecision);

  SelfDoc::Opt kx = {"KX",SelfDoc::Opt::Float,"4.0","4.0","0.0001","NA",
                     "Prior strength for Homs"};
  opts.push_back(kx);
  SelfDoc::Opt kh = {"KH",SelfDoc::Opt::Float,"0.2","0.2","0.0001","NA",
                     "Prior strength for Hets"};
  opts.push_back(kh);

  SelfDoc::Opt v = {"V",SelfDoc::Opt::Float,"10","10","0.0001","NA",
                    "Prior strength for variances"};
  opts.push_back(v);

  SelfDoc::Opt aam = {"AAM",SelfDoc::Opt::Float,"0.66","0.66","NA","NA",
                      "Prior location of AA mean"};
  opts.push_back(aam);
  SelfDoc::Opt abm = {"ABM",SelfDoc::Opt::Float,"0","0","NA","NA",
                      "Prior location of AB mean"};
  opts.push_back(abm);
  SelfDoc::Opt bbm = {"BBM",SelfDoc::Opt::Float,"-.66","-.66","NA","NA",
                      "Prior location of BB mean"};
  opts.push_back(bbm);
  SelfDoc::Opt bby = {"BBY",SelfDoc::Opt::Float,"10.0","10.0","NA","NA",
                      "Prior location of BB Y mean"};
  opts.push_back(bby);
  SelfDoc::Opt aby = {"ABY",SelfDoc::Opt::Float,"10.0","10.0","NA","NA",
                      "Prior location of AB Y mean"};
  opts.push_back(aby);
  SelfDoc::Opt aay = {"AAY",SelfDoc::Opt::Float,"10.0","10.0","NA","NA",
                      "Prior location of AA Y mean"};
  opts.push_back(aay);

  SelfDoc::Opt aav = {"AAV",SelfDoc::Opt::Float,".005",".005","0","NA",
                      "Prior variance of AA"};
  opts.push_back(aav);
  SelfDoc::Opt abv = {"ABV",SelfDoc::Opt::Float,".010",".010","0","NA",
                      "Prior variance of AB"};
  opts.push_back(abv);
  SelfDoc::Opt bbv = {"BBV",SelfDoc::Opt::Float,".005",".005","0","NA",
                      "Prior variance of BB"};
  opts.push_back(bbv);

  // within cluster variances/covariances
 SelfDoc::Opt aayv = {"AAYV",SelfDoc::Opt::Float,"0.1","0.1","0","NA",
                      "Prior Variance of AA, Y coordinate"};
  opts.push_back(aayv);
  SelfDoc::Opt abyv = {"ABYV",SelfDoc::Opt::Float,"0.1","0.1","0","NA",
                       "Prior Variance of AB, Y coordinate"};
  opts.push_back(abyv);
  SelfDoc::Opt bbyv = {"BBYV",SelfDoc::Opt::Float,"0.1","0.1","0","NA",
                       "Prior Variance of BB, Y coordinate"};
  opts.push_back(bbyv);
  // within cluster covariances
  SelfDoc::Opt aaxy = {"AAXY",SelfDoc::Opt::Float,"0","0","NA","NA",
                       "Prior CoVariance of AA, XY"};
  opts.push_back(aaxy);
  SelfDoc::Opt abxy = {"ABXY",SelfDoc::Opt::Float,"0","0","NA","NA",
                       "Prior CoVariance of AB, XY"};
  opts.push_back(abxy);
  SelfDoc::Opt bbxy = {"BBXY",SelfDoc::Opt::Float,"0","0","NA","NA",
                       "Prior CoVariance of BB, XY "};
  opts.push_back(bbxy);

  // between cluster covariances X to X
  SelfDoc::Opt kxx = {"KXX", SelfDoc::Opt::Float,"0","0","NA","NA",
                      "Prior strength for hom covariance"};
  opts.push_back(kxx);
  SelfDoc::Opt kah = {"KAH", SelfDoc::Opt::Float,"0","0","NA","NA",
                      "Prior strength for A-H covariance"};
  opts.push_back(kah);
  SelfDoc::Opt khb = {"KHB", SelfDoc::Opt::Float, "0","0", "NA","NA",
                      "Prior strength for H-B covariance"};
  opts.push_back(khb);
  // between cluster covariances Y to Y
  SelfDoc::Opt kyah = {"KYAH",SelfDoc::Opt::Float,"0","0","NA","NA",
                       "Prior strength for Y-Y A-H cluster"};
  opts.push_back(kyah);
  SelfDoc::Opt kyab = {"KYAB",SelfDoc::Opt::Float,"0","0","NA","NA",
                       "Prior Strength for Y-Y A-B cluster"};
  opts.push_back(kyab);
  SelfDoc::Opt kyhb = {"KYHB",SelfDoc::Opt::Float,"0","0","NA","NA",
                       "Prior Strength for Y-Y H-B cluster"};
  opts.push_back(kyhb);

  SelfDoc::Opt comvar = {"COMVAR",SelfDoc::Opt::Integer,"1","1","NA","NA",
                         "Flag: common variance"};
  opts.push_back(comvar);
  SelfDoc::Opt hard = {"HARD",SelfDoc::Opt::Integer,"2","2","NA","NA",
                       "Flag: type of hard shell"};
  opts.push_back(hard);
  SelfDoc::Opt sb   = {"SB",SelfDoc::Opt::Float,".05",".05","0","NA",
                       "Size of shell barrier"};
  opts.push_back(sb);
  SelfDoc::Opt cm   = {"CM",SelfDoc::Opt::Integer,"0","0","NA","NA",
                       "Type of call method, CM=1 for posterior, CM=2 for single-sample"};
  opts.push_back(cm);
  SelfDoc::Opt maxscore = {"MS",SelfDoc::Opt::Float,".2",".2","0","2",
                           "Threshold for no-calls"};
  opts.push_back(maxscore);
  SelfDoc::Opt bins = {"bins",SelfDoc::Opt::Integer,"0","0","NA","NA",
                       "Use efficient binning to speed up labeling"};
  opts.push_back(bins);
  SelfDoc::Opt hints = {"hints",SelfDoc::Opt::Integer,"0","0","NA","NA",
                        "Use reference genotype data to indicate clusters"};
  opts.push_back(hints);

  SelfDoc::Opt override = {"override",SelfDoc::Opt::Integer,"0","0","NA","NA",
                           "Use reference genotype data to select probes"};
  opts.push_back(override);


  SelfDoc::Opt mix = {"mix", SelfDoc::Opt::Integer,"0","0","NA","NA",
                      "Apply mixture frequency penalty to clusters"};
  opts.push_back(mix);
  SelfDoc::Opt bic = {"bic", SelfDoc::Opt::Float, "0","0", "NA","NA",
                      "BIC penalty for clusters"};
  opts.push_back(bic);

  SelfDoc::Opt CSepPen = {"CSepPen", SelfDoc::Opt::Float, "0","0", "NA","NA",
                      "Penalty for ClusterSep/Stdev too small"};
  opts.push_back(CSepPen);
  SelfDoc::Opt CSepThr = {"CSepThr",SelfDoc::Opt::Float,"16","16","0.000000001","NA",
                         "ClusterSep/Stdev threshold discounting by Geman-McClure"};
  opts.push_back(CSepThr);

  SelfDoc::Opt IsoHetY = {"IsoHetY",SelfDoc::Opt::Float,"0","0","NA","NA",
	  			"Het mean must lie above the line connecting hom means"};
  opts.push_back(IsoHetY);


  SelfDoc::Opt lambda = {"lambda",SelfDoc::Opt::Float,"1","1","0","1",
                         "Controls mixing of common variances"};
  opts.push_back(lambda);
  SelfDoc::Opt wobble = {"wobble",SelfDoc::Opt::Float,".0001",".0001","0","NA",
                         "Limits prior pseudo-observations to 1/wobble"};
  opts.push_back(wobble);
  SelfDoc::Opt copyqc = {"copyqc",SelfDoc::Opt::Float,"0","0","0","1",
                         "Test for outlier size values (CNV or errors)"};
  opts.push_back(copyqc);
  SelfDoc::Opt copytype = {"copytype",SelfDoc::Opt::Integer,"0","0","NA","NA",
                           "Flag: Method for handling outlier data points(CNV/errors)"};
  opts.push_back(copytype);

  SelfDoc::Opt ocean = {"ocean",SelfDoc::Opt::Float,"0","0","0","1",
                        "Test datapoints against uniform ocean probability"};
  opts.push_back(ocean);

  SelfDoc::Opt inflatePRA = {"inflatePRA",SelfDoc::Opt::Float,"0","0","0","NA",
                             "Make calls adding uncertainty in mean to observed variance"};
  opts.push_back(inflatePRA);

  SelfDoc::Opt clustertype = {"clustertype",SelfDoc::Opt::Integer,"1","1","NA","NA",
                              "Flag: type of cluster (1-d, etc)"};
  opts.push_back(clustertype);

  SelfDoc::Opt cp = {"CP", SelfDoc::Opt::Float,"16","16", "0", "NA",
                     "Penalty for contradicting reference genotype"};
  opts.push_back(cp);
  SelfDoc::Opt hok = {"Hok", SelfDoc::Opt::Integer,"0","0", "NA", "NA",
                      "Allow Hints to be flipped in genotype"};
  opts.push_back(hok);
  SelfDoc::Opt x_em_thresh = {"em_thresh",SelfDoc::Opt::Float,"0.05","0.05","0","NA",
                              "set threshold for em gender routine"};
  opts.push_back(x_em_thresh);
  SelfDoc::Opt x_em_cutoff = {"em_cutoff",SelfDoc::Opt::Float,"0.5","0.5","0","NA",
                              "set cutoff for em gender routine"};
  opts.push_back(x_em_cutoff);
  SelfDoc::Opt x_gender_cutoff = {"gender_cutoff",SelfDoc::Opt::Float,"0.1","0.1","0","NA",
                                  "set cutoff for which gender in em gender call"};
  opts.push_back(x_gender_cutoff);
  SelfDoc::Opt copynumber = {"copynumber",SelfDoc::Opt::Boolean,"false","false","NA","NA",
                                  "Is this method begin called as part of a copynumber run?"};
  opts.push_back(copynumber);
  SelfDoc::Opt trustcheck = {"TrustCheck",SelfDoc::Opt::Boolean,"false","false","NA","NA",
            "Completely untrusted, blimeshed probesets in samples are no-calls"};
  opts.push_back(trustcheck);

  return opts;
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc QuantLabelZ::explainSelf() {
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
SelfCreate *QuantLabelZ::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  bool lowprecision = false;
  enum Transformation transform = CCS;
  double K=4;
  string transformStr;
  map<string,string>::iterator iter;

  fillInValue(K, "K", param, doc);

  fillInValue(transformStr, "transform", param, doc);
  transform = transformationForString(transformStr);
  QuantLabelZ *labelz = new QuantLabelZ(transform, K, lowprecision);

  snp_param ts;
  ts.Initialize();
  // prior strength
  fillInValue(ts.prior.aa.k,"KX",param,doc);
  ts.prior.bb.k=ts.prior.aa.k;
  fillInValue(ts.prior.ab.k,"KH",param,doc);

  fillInValue(ts.prior.ab.v,"V",param,doc);
  ts.prior.aa.v=ts.prior.ab.v;
  ts.prior.bb.v=ts.prior.ab.v;
  // prior centers
  fillInValue(ts.prior.bb.m,"BBM",param,doc);
  fillInValue(ts.prior.ab.m,"ABM",param,doc);
  fillInValue(ts.prior.aa.m,"AAM",param,doc);
  // prior y coordinates
  fillInValue(ts.prior.bb.ym,"BBY",param,doc);
  fillInValue(ts.prior.ab.ym,"ABY",param,doc);
  fillInValue(ts.prior.aa.ym,"AAY",param,doc);
  // prior variances
  fillInValue(ts.prior.aa.ss,"AAV",param,doc);
  fillInValue(ts.prior.bb.ss,"BBV",param,doc);
  fillInValue(ts.prior.ab.ss,"ABV",param,doc);
  // prior y variances
  fillInValue(ts.prior.aa.yss,"AAYV",param,doc);
  fillInValue(ts.prior.ab.yss,"ABYV",param,doc);
  fillInValue(ts.prior.bb.yss,"BBYV",param,doc);
  // prior xy covariances
  fillInValue(ts.prior.aa.xyss,"AAXY",param,doc);
  fillInValue(ts.prior.ab.xyss,"ABXY",param,doc);
  fillInValue(ts.prior.bb.xyss,"BBXY",param,doc);
  // X-X covariances between clusters
   fillInValue(ts.prior.xab,"KXX",param,doc);
  fillInValue(ts.prior.xah,"KAH",param,doc);
  fillInValue(ts.prior.xhb,"KHB",param,doc);
  // Y-Y covariances between clusters
  fillInValue(ts.prior.yah,"KYAH",param,doc);
  fillInValue(ts.prior.yab,"KYAB",param,doc);
  fillInValue(ts.prior.yhb,"KYHB",param,doc);
  // wacky parameters
  fillInValue(ts.comvar,"COMVAR",param,doc);
  fillInValue(ts.hardshell,"HARD",param,doc);
  fillInValue(ts.shellbarrier,"SB",param,doc);
  fillInValue(ts.callmethod,"CM",param,doc);
  fillInValue(ts.bins,"bins",param,doc);
  fillInValue(ts.hints,"hints",param,doc);
  fillInValue(ts.contradictionpenalty, "CP",param,doc);
  fillInValue(ts.hok,"Hok",param,doc);
  fillInValue(ts.mix,"mix",param,doc);
  fillInValue(ts.bic,"bic",param,doc);
  fillInValue(ts.CSepPen,"CSepPen",param,doc);
  fillInValue(ts.CSepThr,"CSepThr",param,doc);
  fillInValue(ts.lambda,"lambda",param,doc);
  fillInValue(ts.wobble,"wobble",param,doc);
  fillInValue(ts.copyqc,"copyqc",param,doc);
  fillInValue(ts.copytype,"copytype",param,doc);
  fillInValue(ts.clustertype,"clustertype",param,doc);
  fillInValue(ts.ocean,"ocean",param,doc);
  fillInValue(ts.inflatePRA,"inflatePRA",param,doc);
  fillInValue(ts.IsoHetY,"IsoHetY",param,doc);
  //printf("BBM: %f AAM: %f CM: %d\n",ts.bb.m,ts.aa.m,ts.callmethod);
  //ts.aa.Dump();
  //ts.ab.Dump();
  //ts.bb.Dump();
  labelz->setLabelParam(ts);
  double t_em_cutoff,t_em_thresh,t_gender_cutoff;
  fillInValue(t_em_cutoff,"em_cutoff",param,doc);
  fillInValue(t_em_thresh,"em_thresh",param,doc);
  fillInValue(t_gender_cutoff,"gender_cutoff",param,doc);
  labelz->setGenderParam(t_em_thresh,t_em_cutoff,t_gender_cutoff);
  double score;
  fillInValue(score,"MS",param,doc);
  labelz->setMaxScore(score);
  int overrideflag;
  fillInValue(overrideflag,"override",param,doc);
  labelz->setOverride(overrideflag);
  bool bCopyNumber = false;
  fillInValue(bCopyNumber,"copynumber",param,doc);
  labelz->setCopyNumber(bCopyNumber);
  bool trust_check = false;
  fillInValue(trust_check,"TrustCheck",param,doc);
  labelz->setTrustCheck(trust_check);
  return labelz;
}

//////////////////////////////
//////////////////////////////
//
//
// QuantLableZ IO functions.
//


//////////

/**
 * Writes normalization functions to a file
 *
 * @param nMap - normalization functions
 * @param fileName - the file to write
 */
void QuantLabelZ::writeNormMap(std::vector<NormalizationPredictor> &nMap,
                               const std::string& fileName)  {
  affx::TsvFile tsv;
  //
  tsv.defineColumn(0,0,"Sample");
  tsv.defineColumn(0,1,"FitBB");
  tsv.defineColumn(0,2,"FitAB");
  tsv.defineColumn(0,3,"FitAA");
  tsv.defineColumn(0,4,"sigma");
  tsv.defineColumn(0,5,"Cov");
  tsv.defineColumn(0,6,"Targ");
  //
  tsv.writeTsv_v1(fileName);

  //
  for (int i=0;i<nMap.size(); i++) {
    //
    tsv.set(0,0,i);
    // these are all vectors values.
    tsv.set(0,1,nMap[i].fitBB);
    tsv.set(0,2,nMap[i].fitAB);
    tsv.set(0,3,nMap[i].fitAA);
    tsv.set(0,4,nMap[i].sigma);
    tsv.set(0,5,nMap[i].CenterCovariates);
    tsv.set(0,6,nMap[i].TargetCenter);
    //
    tsv.writeLevel(0);
  }
  tsv.close();
}

//////////

void QuantLabelZ::writeNormSummaryTsv(const std::string& fileName, 
                                      std::vector<std::string> &headers,
                                      affx::TsvReport::TsvReportFmt_t format) {
  m_NormSummaryTsv.setFilename(fileName);
  m_NormSummaryTsv.setFormat(format);
  //
  // @todo: Is this the best format for the report? -jhg
  // we could skip the "Delta" or "Size" by using more columns or seperate reports.
  m_NormSummaryTsv.addHeaderComment("minimal file for normalized values");
  // the name as a string
  m_NormSummaryTsv.defineStringColumn(0,0,"probeset_id",TSVREPORT_PROBESET_STRLEN);
  // "Delta" or "Size"
  m_NormSummaryTsv.defineStringColumn(0,1,"datatype",TSVREPORT_PROBESET_STRLEN);
  // 
  m_NormSummaryTsv.defineColumns(Fs::basename(headers),affx::FILE5_DTYPE_DOUBLE,0);
  //
  m_NormSummaryTsv.writeTsv_v1();
}

/**
 * Writes a summary of the transformed values to a file - i.e as set up for genotyping, not expression
 *
 * @param m_NormSummaryOut - output stream
 * @param m_AValues - first coordinate of transformed data - always "contrast"-like values
 * @param m_BValues - second coordinate of transformed data - always "strength"-like values
 * @param m_ProbeSetName - what is our SNP?
 */
void QuantLabelZ::writeNormSummaryValue(const std::string& m_ProbeSetName,
                                        const vector<double>& m_AValues,
                                        const vector<double>& m_BValues) {
  // delta line first.
  m_NormSummaryTsv.set_string(0,0,m_ProbeSetName);
  m_NormSummaryTsv.set_string(0,1,"Delta");
  for (int i=0; i<m_AValues.size(); i++) {
    m_NormSummaryTsv.set_d(0,2+i,m_AValues[i]);
  }
  m_NormSummaryTsv.writeLevel(0);

  // Then the "Size" line.
  m_NormSummaryTsv.set_string(0,0,m_ProbeSetName);
  m_NormSummaryTsv.set_string(0,1,"Size");
  for (int i=0; i<m_BValues.size(); i++) {
    m_NormSummaryTsv.set_d(0,2+i,m_BValues[i]);
  }
  m_NormSummaryTsv.writeLevel(0);
}

//////////

/**
 * Write header to the selection file so people don't have to look up everything
 *
 * @param m_Out - the file for a header
 */
void QuantLabelZ::writeSnpProbeTsv(const std::string& fileName,
                                   affx::TsvReport::TsvReportFmt_t format) {
  //
  m_SnpProbeTsv.setFilename(fileName);
  m_SnpProbeTsv.setFormat(format);
  // write out a header in the file
  // describing what we actually are doing
  int c=0;
  m_SnpProbeTsv.defineStringColumn(0,c++,"ProbeSet",TSVREPORT_PROBESET_STRLEN);
  m_SnpProbeTsv.defineColumn(0,c++,"A",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"ProbeNum",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"ProbeAindex",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"ProbeBindex",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"HomLogisticCorrect",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"HetLogisticCorrect",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"LogisticCorrect",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"MeanAA",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"MeanAB",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"MeanBB",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"nAA",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"nAB",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"nBB",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"Stdev",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"FLDAH",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"FLDHB",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"FLDAB",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"Entropy",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"AIC",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"B0",affx::FILE5_DTYPE_DOUBLE);
  m_SnpProbeTsv.defineColumn(0,c++,"B1",affx::FILE5_DTYPE_DOUBLE);
  //
  m_SnpProbeTsv.writeTsv_v1();
}

/**
 * Writeprobe writes the summaries to the file
 *
 * @param m_Out - stream required
 * @param probenum - what probe we are
 * @param Aid - probe identity for A
 * @param Bid - probe identity for B
 * @param Summaries - summary information for probe fit
 * @param m_Ax - A data (fitted)
 * @param m_Bx - B data (fitted - not used)
 * @param m_ProbeSetName
 */
void QuantLabelZ::writeSnpProbeValue(affx::TsvReport& tsv,
                                     int probenum,
                                     unsigned int Aid,
                                     unsigned int Bid,
                                     const vector<double> &Summaries,
                                     const vector<double> &m_Ax,
                                     const vector<double> &m_Bx,
                                     string& m_ProbeSetName) {
  // todo add tsvdata here.
  int c=0;
  tsv.set_string(0,c++,m_ProbeSetName);
  tsv.set_string(0,c++,"A"); // This is defined as double above.
  tsv.set_d(0,c++,probenum);
  tsv.set_d(0,c++,Aid);
  tsv.set_d(0,c++,Bid);
  //
  for (int i=0; i<Summaries.size(); i++) {
    tsv.set_d(0,c++,Summaries[i]);
  }
  for (int i=0;i<m_Ax.size();i++) {
    tsv.set_d(0,c++,m_Ax[i]);
  }
  //
  tsv.writeLevel(0);
}

//////////

/// @brief Pick the output version based on the file format.  (txt=>v1,a5=>v2)
void QuantLabelZ::writeSnpPosteriorTsv(const std::string& fileName,
                                       affx::TsvReport::TsvReportFmt_t format_file) {
  int format_ver = 1;
  if (format_file==TsvReport::FMT_TSV) {
    format_ver=1;
  }
  else if (format_file==TsvReport::FMT_A5) {
	  if (m_bCopyNumber) {
      format_ver = 2;
    }
	  else {
      format_ver = 3;
    }
  }
  else {
    Err::errAbort("writeSnpPosteriorTsv");
  }
  //
  writeSnpPosteriorTsv(fileName,format_file,format_ver);
}

void QuantLabelZ::writeSnpPosteriorTsv(const std::string& fileName,
                                       affx::TsvReport::TsvReportFmt_t format_file,
                                       int format_ver) {
  // remember the version being written (1 or 2)
  m_SnpPosteriorTsv_ver=format_ver;
  m_SnpPosteriorName = fileName;
  m_PosteriorFormat = format_file;
}
//   //
//     std::vector<std::string>::const_iterator keyIx, paramIx;
//   //
//     Verbose::out(1, "Got " + ToStr(m_Info.m_ParamValues.size()) + " headers to write.");
//   for (keyIx = m_Info.m_ParamNames.begin(), paramIx = m_Info.m_ParamValues.begin();
//       keyIx != m_Info.m_ParamNames.end() && paramIx != m_Info.m_ParamValues.end();
//       ++keyIx, ++paramIx) {
//     // @todo should we using the AffymetrixParameterConsts.h #defined values?
//     m_SnpPosteriorTsv.addHeader("affymetrix-algorithm-param-" + *keyIx, *paramIx);
//   }
//   for (keyIx = m_Info.m_ClientInfoNames.begin(), paramIx = m_Info.m_ClientInfoValues.begin();
//       keyIx != m_Info.m_ClientInfoNames.end() && paramIx != m_Info.m_ClientInfoValues.end();
//       ++keyIx, ++paramIx) {
//     // @todo should we using the AffymetrixParameterConsts.h #defined values?
//     m_SnpPosteriorTsv.addHeader("affymetrix-application-meta-data-info-" + *keyIx, *paramIx);
//   }
//   QuantLabelZ__writeSnpPosteriorTsv(m_SnpPosteriorTsv,m_SnpPosteriorName,format_file,format_ver);
//}

/**
 *  Writes one specific posterior to a file along with snp name
 *
 * @param tsp - snp specific parameters, including posterior distribution of cluster centers
 * @param out - output stream
 * @param id - snp id
 * @param fieldDelim - always tab
 * @param sep - always comma
 */
void QuantLabelZ::writeSnpPosteriorValue(const std::string& probeset_id,
                                         snp_param &tsp) {
  QuantLabelZ__writeSnpPosteriorValue(m_SnpPosteriorTsv,m_SnpPosteriorTsv_ver,
                                      probeset_id,tsp);
}

void QuantLabelZ::readSnpPriorMap(const std::string& fileName) {
  QuantLabelZ__readSnpPriorMap(m_vectSnpPriors,fileName);
}
void QuantLabelZ::readSnpPriorMap_tsv5(affx::File5_Tsv* tsv5) {
  QuantLabelZ__readSnpPriorMap_tsv5(m_vectSnpPriors,tsv5);
}

void QuantLabelZ::setParameters(PsBoard &board) {
  Options *o = board.getOptions();
  vector<string> celFiles = o->getOptVector("cels");
  string outDir = o->getOpt("out-dir");
  if (o->getOpt("read-models-brlmmp") != "") {
      if (!o->getOptBool("brlmmp-models-non-sequential")) {
          m_SequentialModelFile = o->getOpt("read-models-brlmmp");
          m_SequentialModelTsv.open(m_SequentialModelFile);
          m_SequentialModelTsv.bind(0,"id", &m_Id, affx::TSV_BIND_REQUIRED);
          m_SequentialModelTsv.bind(0,"BB", &m_BB, affx::TSV_BIND_REQUIRED);
          m_SequentialModelTsv.bind(0,"AB", &m_AB, affx::TSV_BIND_REQUIRED);
          m_SequentialModelTsv.bind(0,"AA", &m_AA, affx::TSV_BIND_REQUIRED);
          
          int fc_cidx =  m_SequentialModelTsv.cname2cidx(0,"CV");
          if (fc_cidx!=affx::TSV_ERR_NOTFOUND) {
              m_SequentialModelTsv.bind(0, "CV", &m_CV, affx::TSV_BIND_REQUIRED);
          }
          else {
              m_CV = "0,0,0";
          }
          
      }
      else {
          readSnpPriorMap(o->getOpt("read-models-brlmmp"));
      }
  }
  
  // @todo refactor - put current analysis path name on bboard
  string outfile = Fs::join(outDir,getType() + ".snp-posteriors");
  bool writeModels = o->getOptBool("write-models");
  if (writeModels) {
    writeSnpPosteriorTsv(outfile,affx::TsvReport::FMT_TSV);
    string modelFile = outfile + ".txt";
    board.set("model-file-written", modelFile);
  }
  // @todo refactor - put current analysis path name on bboard
  outfile =Fs::join(outDir,getType() + ".normalized-summary");
  if (o->getOptBool("summaries")) {
    writeNormSummaryTsv(outfile, celFiles,affx::TsvReport::FMT_TSV);
  }
  if (o->getOptBool("select-probes")) {
    // @todo refactor - put current analysis path name on bboard
    outfile =Fs::join(outDir,getType() + ".select-probes");
    writeSnpProbeTsv(outfile,affx::TsvReport::FMT_TSV);
  }

  // Add quantification method.
  QuantMethodFactory qFactory;
  string qMethodSpec = o->getOpt("qmethod-spec");
  QuantMethod *qMethod = qFactory.quantMethodForString(qMethodSpec, board);
  QuantExprMethod *eMethod = NULL;
  eMethod= static_cast<QuantExprMethod *>(qMethod);
  if (eMethod == NULL) {
    Err::errAbort("Couldn't make a QuantExprMethod for specification: " + qMethodSpec);
  }
  setQuantExprMethod(eMethod);
  
  GenotypeInfo *gInfo = board.getGenotypeInfo();
  if (gInfo == NULL) {
    Err::errAbort("QuantLabelZ::setParameters() - GenotypeInfo is NULL");
  }
  if (!gInfo->m_SpecialSnps->empty())
    setSpecialSnps(*(gInfo->m_SpecialSnps));
  if (!gInfo->m_HaploidSnps->empty())
    setHaploidSnps(*(gInfo->m_HaploidSnps));
  setGenders(gInfo->m_Genders->getGenders());
  // add sample covariate inbreeding penalty here
  // this is highly indirect at this point, but perhaps there's a general model for sample covariates
  setInbredHetPenalty(gInfo->m_Inbred->getInbredStatus());
  if (o->isOptDefined("output-probabilities")) {
    m_OutputProbabilities = (o->getOpt("output-probabilities")!="");
  }
  try {
    board.get(ArtifactReduction::getProbesetTrustFileKey(), &m_ProbeSetTrustTmpFileName);
  }
  catch(...) {
    // HACK: board.get takes abortOnError as a parameter but that parameter is ignored. 
    // just ignore it.
  }
  
}
