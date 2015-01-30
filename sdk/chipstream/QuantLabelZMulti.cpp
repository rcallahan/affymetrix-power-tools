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
 * @file   QuantLabelZMulti.cpp
 * @author Martin Gilchrist 
 * @date   Tue Jun 24 1:50:00 2008
 *
 *
 */

// The following typedefs are here for quick reference, avoids having to bring up .h file.
/*
typedef vector<double> SampleSummaries;
typedef map<int, SampleSummaries> ContextMap;
typedef map<int, ContextMap> AlleleMap;

typedef pair<int, int> AlleleContext;
typedef pair<AlleleContext, double> ASummary;
typedef vector<ASummary> AlleleVector;
typedef vector<ASummary>  MaxContextVector;
typedef vector<MaxContextVector> SampleMaxs;
*/
#include "chipstream/QuantLabelZMulti.h"
#include "chipstream/QuantLabelZIO.h"
//
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantMethodFactory.h"
//
#include "calvin_files/utils/src/GenoCallCoder.h"
#include "util/AffxString.h"
//
#include <iostream>
//

/**
  @brief  Gets the integer valued (not affx::GType) forced call for the probeset. 
*/
int QuantLabelZMulti::getForcedCall(unsigned int index) {
  if(index >= m_CallsMulti.size()) {
    Err::errAbort("Asking for call at index " + ToStr(index) + " when Probeset " +
                m_ProbesetName + " has only " + ToStr(m_CallsMulti.size()) + " calls.");
  }
  return m_CallsMulti[index];
}


/**
  @brief  Gets the integer valued (not affx::GType) call for the probeset. 
*/
int QuantLabelZMulti::getCall(unsigned int index) {
  int call = getForcedCall(index);
  double conf = getConfidence(index);

  if(m_CopyNumbers[index]==0)
    return (m_coder->abstractAlleleToGenotypeCallNum("ZeroCopyNumber"));
  else if (m_PriorObs[index] < m_PraThresh)
      return(m_coder->abstractAlleleToGenotypeCallNum("PossibleRareAllele"));
  else if((conf >= getMinThresh()) && (conf <= getMaxThresh()))
      return call;
  else
      return(m_coder->abstractAlleleToGenotypeCallNum("NoCall"));
}


/** 
 * Get the genotype forced call at specified index (sample).
 * @param index - sample of interest.
 * @return - Genotyping call made.
 */
affx::GType QuantLabelZMulti::getGTypeForcedCall(unsigned int index) {
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
**/
affx::GType QuantLabelZMulti::getGTypeCall(unsigned int index) {
  affx::GType call = getGTypeForcedCall(index);
  if(m_CopyNumbers[index]==0)
    ///@todo  This only works since Summary statistics only need hom vs. het information.
    return affx::AA;
  double conf = getConfidence(index);
  if((conf >= getMinThresh()) && (conf <= getMaxThresh()))
      return call;
  else
      return affx::NN;
}


/** Fill in our self documentation from current state. */
void QuantLabelZMulti::setupSelfDoc(SelfDoc &doc) {
  doc.setDocName(QUANTBRLMMPMULTI);
  doc.setDocDescription("Do genotyping calls with BRLMM-P-MULTI (perfect match) algorithm.");
  std::vector<SelfDoc::Opt> opts = QuantLabelZ::getDefaultDocOptions();

  SelfDoc::Opt praThresh = {"pra-thresh",SelfDoc::Opt::Integer,"0","0","0","NA",
                         "Threshold on cluster mean strength below which a PRA is called."};
  opts.push_back(praThresh);

  SelfDoc::Opt ccAlleles = {"cc-alleles",SelfDoc::Opt::Integer,"6","6","0","NA",
                         "The max number of alleles for call encoding/decoding."};
  opts.push_back(ccAlleles);

  SelfDoc::Opt ccType = {"cc-type",SelfDoc::Opt::String,"UCHAR","UCHAR","0","NA",
                         "Call encoding/decoding data size."};
  opts.push_back(ccType);

  SelfDoc::Opt ccVersion = {"cc-version",SelfDoc::Opt::String,"1.0","1.0","0","NA",
                         "Call encoding/decoding version."};
  opts.push_back(ccVersion);

  doc.setDocOptions(opts);
}



/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc QuantLabelZMulti::explainSelf() {
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
SelfCreate *QuantLabelZMulti::newObject(std::map<std::string,std::string> &param) {
  SelfDoc doc = explainSelf();
  bool lowprecision = false;
  enum Transformation transform = CCS;
  double K=4;
  int praThresh = 3;
  string transformStr;
  int ccAlleles = 6;
  string ccType = "UCHAR";
  string ccVersion = "1.0";
  map<string,string>::iterator iter;

  fillInValue(K, "K", param, doc);
  fillInValue(praThresh, "pra-thresh", param, doc);
  fillInValue(ccAlleles, "cc-alleles", param, doc);
  fillInValue(ccType, "cc-type", param, doc);
  fillInValue(ccVersion, "cc-version", param, doc);

  fillInValue(transformStr, "transform", param, doc);
  transform = transformationForString(transformStr);
  QuantLabelZMulti *labelz = new QuantLabelZMulti(transform, K, lowprecision, praThresh, ccAlleles, ccType, ccVersion);

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

  return labelz;
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
bool QuantLabelZMulti::setUp(ProbeSetGroup &psGroup, const IntensityMart &iMart, 
                             std::vector<ChipStream *> &iTrans, PmAdjuster &pmAdjust) {

  const ProbeSet *gtPs = NULL; // Genotyping probeset

  bool success = true;
  /* Sanity checks about probesets. */
  if(psGroup.probeSets.empty())
    Err::errAbort("Zero probesets in ProbeSetGroup (group: " + ToStr(psGroup.name) + ").");
  if(psGroup.probeSets.size() > 1)
    Err::errAbort("Can't have multiple probesets in a genotyping ProbeSetGroup (group: " + ToStr(psGroup.name) + ").");

  /* Remember this probeset. */
  gtPs = psGroup.probeSets[0];
  success &= SetUpProbeSet(gtPs, iMart, iTrans, pmAdjust, true);

  return success;
}

void QuantLabelZMulti::amalgamateAtoms( const ProbeSet* gtPs, 
                                        ProbeSet &amalgamatedProbeset){


  // This will the set of indices of the Atoms in gtPs that we amalgamate.  This vector is filled several times,
  // once for each allele/context combination that is seen in the Atoms of gtPs. 
  vector<int> toBeAmalgamatedThisTime;

  int atomCount=0;
  int allele=0;
  int context=0;
  int numberOfAtoms=gtPs->atoms.size();
  vector<int> alreadyAmalgamated(numberOfAtoms, 0);
  int insertedAtomCount=0; 
  while(atomCount < numberOfAtoms){
    if( alreadyAmalgamated[atomCount]){
      atomCount++;
      continue;
    }
    //  We have found an Atom indexed by atomCount that has not been amalgamated.  We label it amalgamated.
    toBeAmalgamatedThisTime.erase(toBeAmalgamatedThisTime.begin(), toBeAmalgamatedThisTime.end());
    alreadyAmalgamated[atomCount]=1;
    toBeAmalgamatedThisTime.push_back(atomCount);

    //  We find the allele and context values for the Atom in gtPs indexed by "atomCount". 
    allele=gtPs->atoms[atomCount]->getAlleleCode();
    context=gtPs->atoms[atomCount]->getContextCode();

    //  We determine the indices of all the atoms that have the same allele and context values; 
    //  as the Atom indexed by "atomCount".
    for(int atomIndex=atomCount+1; atomIndex < numberOfAtoms; atomIndex++){

      //  If Atom indexed by atomIndex has already been amalgamated we go on to the next Atom.
      if( alreadyAmalgamated[atomIndex])
        continue;

      //  Determine whether the Atom with index "atomIndex" matches our Atom indexed by "atomCount".  
      if( allele==gtPs->atoms[atomIndex]->getAlleleCode() &&
          context==gtPs->atoms[atomIndex]->getContextCode() ){
            alreadyAmalgamated[atomIndex]=1;
            toBeAmalgamatedThisTime.push_back(atomIndex);
      }//  end determination of status of Atom indexed by "atomIndex".   
    }// end determination of the indices of all Atoms having the same allele and context  
     // values as the Atom indexed by "atomCount".


    // We amalgamate all the atoms in the vector "toBeAmalgamatedThisTime". 
    Atom* atomPointer = amalgamatedProbeset.atoms[insertedAtomCount];
    atomPointer->allele_code = allele;
    atomPointer->context_code = context;
    atomPointer->probes.clear();
    vector<int>::iterator begin=toBeAmalgamatedThisTime.begin();
    vector<int>::iterator end=toBeAmalgamatedThisTime.end();
    for( ; begin!=end; begin++){
      vector<Probe*>::iterator probeIterBegin = ((gtPs->atoms[*begin])->probes).begin();
      vector<Probe*>::iterator probeIterEnd = ((gtPs->atoms[*begin])->probes).end();
      for( ; probeIterBegin!=probeIterEnd; probeIterBegin++){
        atomPointer->probes.push_back( *probeIterBegin );
      }
    } // end amalgamation process for a single allele/context value. 

    insertedAtomCount++;
  } // end of while loop
  amalgamatedProbeset.atoms.resize(insertedAtomCount);

}

/*
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
bool QuantLabelZMulti::SetUpProbeSet(   const ProbeSet *gtPs, 
                                        const IntensityMart &iMart, 
                                        std::vector<ChipStream *> &iTrans, 
                                        PmAdjuster &pmAdjust, 
                                        bool doReport){ 

  m_Summaries.erase(m_Summaries.begin(), m_Summaries.end());
  bool success=true;

  //  We store the marker name for subsequent use in computeEstimate.
  m_ProbesetName = gtPs->name;

  bool report = doReport;
  if(m_ProbeSetsToReport != NULL)
      if(m_ProbeSetsToReport->find(m_ProbesetName.c_str())==m_ProbeSetsToReport->end())
          report = false;


  // Each atom in the probeset gtPs contains the probes for a single allele/context/sense.  Since we need to have
  // probes grouped in an atom without regard to the sense parameter, we amalgamate atoms over the sense parameter.
  ProbeSet  amalgamatedProbeset;

  // We need to add a sufficiently large set of "empty" Atoms to the amalgamatedProbeset due to the fact that a Probeset
  // maintains a vector of POINTERS to Atoms.  Thus if we define Atoms within the function "amalgamateAtoms" and then take
  // their addresses and insert the pointers in the vector, we lose the Atom itself when the function "amalgamateAtoms" pops.
  // The method is thus to create some "empty" atoms and then have the function "amalgamateAtoms" fill them with Probes. 
  // (Actually they are filled with pointer to Probes.  The management of the memory for these pointers to Probes is assumed
  // to be properly maintained by the rest of the program.  This is dangerous, but is something we have lived with and will have
  // to live with until the Probeset hierarchy is redone. 

  std::vector<Atom> atomVector(gtPs->atoms.size());
  for(int i=0; i< gtPs->atoms.size(); i++){
    amalgamatedProbeset.atoms.push_back(&(atomVector[i]));
  }  

  amalgamateAtoms(gtPs, amalgamatedProbeset);

  vector<double> summaryValues;
  ProbeSet probeSetToBeSummarized;
  probeSetToBeSummarized.psType = ProbeSet::Expression;
  string name;

  // Here we loop over the atoms in the probeset which has "amalgamated atoms"  and turn them into 
  // probesets for summarization.
  map<int,int> alleles;
  int numberOfAmalgamatedAtoms = amalgamatedProbeset.atoms.size();
  for(int i=0; i<numberOfAmalgamatedAtoms; i++) {

    // Get an atom and create a probeSet.    
    probeSetToBeSummarized.atoms.push_back(amalgamatedProbeset.atoms[i]);
    probeSetToBeSummarized.psType = ProbeSet::Expression;

    // Arrange a correct name for the probeSet.  That is marker-allele-context.
    int suffix1 = amalgamatedProbeset.atoms[i]->getAlleleCode();  
    int suffix2 = amalgamatedProbeset.atoms[i]->getContextCode();  
    name = m_ProbesetName;
    name += "-";
    name += AffxString::intToString(suffix1, false);
    name += "-";
    name += AffxString::intToString(suffix2, false);

    alleles[suffix1]=1;

    probeSetToBeSummarized.name = Util::cloneString(name.c_str());

    // Summarize.
    success &= summarizeAllele(&probeSetToBeSummarized, 
                    summaryValues, iMart, iTrans, 
                    pmAdjust, m_QuantMethod, m_Param.m_LowPrecision, 
                    report, m_Reporters);

    //  We store the number of samples to make the algorithms in computeEstimate easier to follow.
    m_numberOfSamples = summaryValues.size();
    m_Calls.resize(m_numberOfSamples);
    m_CallsMulti.resize(m_numberOfSamples);
    m_PriorObs.resize(m_numberOfSamples);
    m_Context.resize(m_numberOfSamples);
    m_Confidences.resize(m_numberOfSamples);
    m_Distances.resize(3);
    m_Distances[0].resize(m_numberOfSamples);
    m_Distances[1].resize(m_numberOfSamples);
    m_Distances[2].resize(m_numberOfSamples);

    ContextMap newContextMap;

    // Insert our newly computed summaryValues into our member data structure m_Summaries.
    // This data structure is a map of allele values to a map of context values to a vector of summary values.
    // Insertion begin.
    AlleleMap::iterator alleleIterator =  m_Summaries.find(suffix1); 
    if(alleleIterator != m_Summaries.end()){
      // We have already seen this allele.
      ContextMap::iterator contextIterator = (alleleIterator->second).find(suffix2);
      if(contextIterator != (alleleIterator->second).end()){
        // We have already seen this context. Should never get here. Error.
      } else{
        // We put the summary values into the ContextMap for the present allele.
        (alleleIterator->second).insert(ContextMap::value_type(suffix2,summaryValues));
      }
    } else{
      // We have a new allele not yet seen.
      newContextMap.clear();
      newContextMap.insert(ContextMap::value_type(suffix2, summaryValues)); 
      m_Summaries.insert(AlleleMap::value_type(suffix1, newContextMap));              
    }
    // Insertion complete.

    // Clear out the data structures we've used for preparation for the next atom.
    summaryValues.clear();
    name = "";
    delete [] probeSetToBeSummarized.name;
    probeSetToBeSummarized.atoms.clear();
  } // end Atom loop.

  // We need to clear out the amalgamatedProbeset's collection of Atoms, otherwise, when the function pops, the pointers to
  // the Probesets within the atoms get deleted.  The memory management for these pointers is elsewhere. Again, we'll have to
  // live with this until the ProbeSet hierarchy gets modified.
  for(int i=0; i < amalgamatedProbeset.atoms.size(); i++){
    amalgamatedProbeset.atoms[i]->probes.clear();
  }
  amalgamatedProbeset.atoms.clear(); 

  if(alleles.size() < 2)
      return false;

  return(success);
}  // end  SetUpProbeSet


SampleMaxs QuantLabelZMulti::findMaxContextForEachSample(int numberOfSamples){
  SampleMaxs ourSampleMaxs;
  for(int i=0; i<numberOfSamples; i++){

    AlleleMap::iterator begin = m_Summaries.begin();
    AlleleMap::iterator end = m_Summaries.end();
    AlleleVector maxVector;

    // Iterate over all alleles and find the max context for each. 
    for(; begin!=end; begin++){
      int allele = (*begin).first;
      ContextMap presentContext = (*begin).second;
      ContextMap::iterator contextBegin = presentContext.begin();
      ContextMap::iterator contextEnd = presentContext.end();
      // Iterate over all contexts adding all <<allele,context>, presentSummaries[i]> to presentAlleleVector;
      AlleleVector presentAlleleVector;
      for(; contextBegin!=contextEnd; contextBegin++){
        int context = (*contextBegin).first;
        vector<double> presentSummaries = (*contextBegin).second;
        AlleleContext data(allele,context);
        ASummary presentSummary(data, presentSummaries[i]);
        presentAlleleVector.push_back(presentSummary);
      }

      // After we've entered all the summary values for all the contexts for the given allele we sort 
      // on the size of the summary value.  Then put the largest <<allele,context>,summary> into maxVector, 
      // then go on to the next allele and do the same thing. 
      sort(presentAlleleVector.begin(), presentAlleleVector.end(), LessThan());
      maxVector.push_back(presentAlleleVector[0]);
      // Now empty presentAlleleVector so that we can use it for the next allele.
      presentAlleleVector.erase(presentAlleleVector.begin(), presentAlleleVector.end());
    }

    // After we've found the max context for each allele sort and return the two with the largest summary values.
    sort(maxVector.begin(), maxVector.end(), LessThan());
    MaxContextVector sampleContextVector;
    ASummary summ1 = maxVector[0];
    ASummary summ2 = maxVector[1];
    sampleContextVector.push_back(summ1);
    sampleContextVector.push_back(summ2);
    //  We order the two allele context values in the sampleContextVector so as not to do two separate genotyping
    //  calls.  ie.  Don't want to genotype AB and BA separately. 
    sort(sampleContextVector.begin(), sampleContextVector.end(), AlleleLessThan()); 
    ourSampleMaxs.push_back( sampleContextVector );
    // Now empty maxVector preparing to use it for the next sample. 
    maxVector.erase(maxVector.begin(), maxVector.end());
  }
  return ourSampleMaxs;
}

void QuantLabelZMulti::determineCopyNumber(){
  // This function will take into account all the copy number determination for the samples and probeset 
  // being genotyped.  This info is 1) Input from the copyNumberMap which is usually the place where
  // externally computed CN information is held.  2) Gender calculation of CN.

  //  The default copy number is 2 so if none of the mechanisms below change the CN we return 2.
  m_CopyNumbers.resize(m_numberOfSamples);
  for(int i=0; i<m_numberOfSamples; i++){
    m_CopyNumbers[i] = 2;
  }
  // If the CN has been determined by some external method we accept that CN value 
  map<std::string, std::vector<int> >::iterator iter=m_CopyNumberMap.find(m_ProbesetName); 
  if(iter != m_CopyNumberMap.end()){
    m_CopyNumbers = iter->second;
  }else{ 
  // If not we use the information from the Special Snps file and the gender determination method to fix
  // the CN. 
    int MaleCopy, FemaleCopy; 
    map<string,pair<int, int> >::iterator test;
    test=m_SpecialSnps.find(m_ProbesetName);
    if(test==m_SpecialSnps.end()) {
      MaleCopy=2;
      FemaleCopy=2;
    }else{
      MaleCopy=test->second.first;
      FemaleCopy=test->second.second;
    }
    for (int i=0; i<m_numberOfSamples; i++){
      if(m_Genders[i] == affx::Female || m_Genders[i] == affx::UnknownGender)
        m_CopyNumbers[i] = FemaleCopy;
      else if(m_Genders[i] == affx::Male)
        m_CopyNumbers[i] = MaleCopy;
      else
        Err::errAbort("Unknown gender seen while partitioning samples for brlmmp." );
    }

  }

  // We now have the CN predictions and need to translate these into the correct genotyping method.
  for(int i=0; i<m_numberOfSamples; i++){
    int copyNumberValue=m_CopyNumbers[i]; 
    switch(copyNumberValue){
      case -3:
        // CN has been determined to be >=1, so we want to use diploid model. 
        m_CopyNumbers[i] = 2;
        break;

      case -2:
        // CN was not available (issue of consent perhaps?) so we want to use a diploid model.
        m_CopyNumbers[i] = 2;
        break;

      case -1:
        // CN was not successfully determined so we want to use a diploid model.
        m_CopyNumbers[i] = 2;
        break;

      case 0:
        // CN was determined to be 0 so we assert this fact.
        m_CopyNumbers[i] = 0;
        break;

      case 1:
        // CN was determined to be 1 so we assert this fact.
        m_CopyNumbers[i] = 1;
        break;

      case 2:
        // CN was determined to be 2 so we assert this fact.
        m_CopyNumbers[i] = 2;
        break;

       default:     
        // Unknown CN state reported. This means we need to extend this switch.
          Err::errAbort("Unknown copy number state found while partitioning samples.");
    }// end of switch

  }// end of for 

}// end of determineCopyNumber



int QuantLabelZMulti::translateBrlmmpOutputToCall(int copyNumber, 
        int brlmmpCall, int firstAllele, int secondAllele){

  ///@todo  Since this function was written Ray has added functionality to the encode which will eliminate the switches.
  string translateThisString;
  string firstAlleleString;
  string secondAlleleString;
 
  switch(firstAllele){
    case 0:
      firstAlleleString = "A";
      break;

    case 1:
      firstAlleleString = "B";
      break;

    case 2:
      firstAlleleString = "C";
      break;

    case 3:
      firstAlleleString = "D";
      break;

    case 4:
      firstAlleleString = "E";
      break;

    case 5:
      firstAlleleString = "F";
      break;

    case 6:
      firstAlleleString = "G";
      break;
    case 7:
      firstAlleleString = "H";
      break;
    case 8:
      firstAlleleString = "I";
      break; 
   
  }
  
  switch(secondAllele){
    case 0:
      secondAlleleString = "A";
      break;

    case 1:
      secondAlleleString = "B";
      break;

    case 2:
      secondAlleleString = "C";
      break;

    case 3:
      secondAlleleString = "D";
      break;

    case 4:
      secondAlleleString = "E";
      break;

    case 5:
      secondAlleleString = "F";
      break;

    case 6:
      secondAlleleString = "G";
      break;
    case 7:
      secondAlleleString = "H";
      break;
    case 8:
      secondAlleleString = "I";
      break;

  }

  if(copyNumber == 0){
    translateThisString = "ZeroCopyNumber";
  }

  if(copyNumber==1){
    if(brlmmpCall == 3){
      translateThisString = "NoCall"; 
    }
    if(brlmmpCall == 0){
      translateThisString = firstAlleleString;
    }
    if(brlmmpCall == 1){
      Err::errAbort("Unexpected result -- het call for cn=1!");
    }
    if(brlmmpCall == 2){
      translateThisString = secondAlleleString;
    }
  }

  if(copyNumber==2){
    if(brlmmpCall == 3){
      translateThisString = "NoCall";
    } 
    if(brlmmpCall == 0){
      translateThisString = firstAlleleString + firstAlleleString;
    } 
    if(brlmmpCall == 1){
      //  We have to have the following check because the encoder won't accept BA,  only AB
      if(firstAllele < secondAllele){
        translateThisString = firstAlleleString + secondAlleleString;
      } else {
        translateThisString = secondAlleleString + firstAlleleString;
      }

    } 
    if(brlmmpCall == 2){
      translateThisString = secondAlleleString + secondAlleleString;
    }
  }

  unsigned char code;

  code = m_coder->abstractAlleleToGenotypeCallNum(translateThisString);   
  return code; 
} 

void QuantLabelZMulti::getSequentialPrior(std::string &name, snp_param &tsp, snp_labeled_distribution &sDist) {
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

void QuantLabelZMulti::getPrior(std::string &TmpName, int presentCopyNumber, snp_param &tsp, snp_labeled_distribution &sDist) {
  if (m_SequentialModelTsv.is_open()) {
    getSequentialPrior(TmpName, tsp, sDist);
  }
  else if (m_vectSnpPriors.getCount() > 0) {
    snp_labeled_distribution objSearch;
    snp_labeled_distribution* p=NULL; 
    objSearch.probeset_id = TmpName;
    int iSearchIndex = m_vectSnpPriors.binarySearch(objSearch, 0);
    if (iSearchIndex == -1){
      if(presentCopyNumber<2){
        objSearch.probeset_id = "GENERIC:" + ToStr(presentCopyNumber);
      } else {
        objSearch.probeset_id = "GENERIC";
        iSearchIndex = m_vectSnpPriors.binarySearch(objSearch,0);
      }
      if(iSearchIndex==-1){
        // should this fail to global prior?
        Err::errAbort("Can't find model for SNP: " + TmpName);
      }
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

void QuantLabelZMulti::computeEstimate(){


  // The first task is to determine the copy number applied to the samples and probeset genotyped.
  // This information is used in three places.  First for the partitioning of the samples for submission
  // to brlmm-p, the choice of priors for this computation and the computation of the "call-code" which
  // is the coding of the actual genotyping. 
  determineCopyNumber();

  // Input and Output vectors for bayes_label. 
  vector<affx::GType> tcall;
  vector<double> tconf, tx, ty;

  // For each sample find the two alleles and the corresponding contex values which 
  // have the largest summary values. 
  SampleMaxs ourSampleMaxs = findMaxContextForEachSample(m_numberOfSamples);

  //  Create the the input vectors of sample summaries for the strength/contrast computation.
  vector<double> firstAlleleVector, secondAlleleVector;
  for(int i=0; i<m_numberOfSamples; i++){
    firstAlleleVector.push_back( ((ourSampleMaxs[i])[0]).second ); 
    secondAlleleVector.push_back( ((ourSampleMaxs[i])[1]).second ); 
    m_Context[i] = 
            ToStr(((ourSampleMaxs[i][0]).first).first) + "," +  // A Allele: allele code
            ToStr(((ourSampleMaxs[i][0]).first).second) + "," + // A Allele: context code
            ToStr(((ourSampleMaxs[i][1]).first).first) + "," +  // B Allele: allele code
            ToStr(((ourSampleMaxs[i][1]).first).second)          // B Allele: context code
            ;
  }
  // Compute the strength/contrast values and store in original input vectors of summary values.  
  GenoUtility_transformData(firstAlleleVector, secondAlleleVector, m_Param.m_Transform,m_Param.m_K); 

  // Now begin the main section of computeEstimate where we partition the set of samples according to allele pair values
  // and the copy number of each particular sample. We then prepare for submission to the function bayes_label for genotyping. 
  vector<int> alreadyGenotyped(m_numberOfSamples,0);
  int allele1, allele2;
  int presentCopyNumber=2;

  // This will the set of indices of samples that we submit for genotyping to bayes_label. 
  vector<int> toBeGenotypedThisTime; 


  //  The method of partitioning the samples is to start with the first sample, find its allele pair and copy
  //  values.  Then pass through all the other samples, picking out those with the same characteristics, save
  //  their index in "toBeGenotypedThisTime" and mark them as "alreadyGenotyped". Then pass to the first 
  //  ungenotyped sample.
  int count=0;
  while(count < m_numberOfSamples){

    //  If already genotyped we go on to the next sample.
    if( alreadyGenotyped[count] ){
      count++;
      continue;
    }
    //  We have found a sample indexed by count that has not been genotyped.  We label it genotyped.
    toBeGenotypedThisTime.erase(toBeGenotypedThisTime.begin(), toBeGenotypedThisTime.end()); 
    alreadyGenotyped[count]=1;
    toBeGenotypedThisTime.push_back(count);

    //  We find the allele pairs and copy number of the sample indexed by "count". 
    allele1=((ourSampleMaxs[count])[0]).first.first; 
    allele2=((ourSampleMaxs[count])[1]).first.first; 
    presentCopyNumber=m_CopyNumbers[count];


    //  We determine the index of all samples that have the same allele pair values and copy number
    //  as the sample indexed by "count".
    for(int sampleIndex=count+1; sampleIndex < m_numberOfSamples; sampleIndex++){

      //  If sample indexed by sampleIndex has already been genotyped we go on to the next sample.
      if( alreadyGenotyped[sampleIndex])
        continue;

      //  Determine whether the sample with index "sampleIndex" matches our sample indexed by "count".  
      if( allele1==((ourSampleMaxs[sampleIndex])[0]).first.first && 
          allele2==((ourSampleMaxs[sampleIndex])[1]).first.first && 
          presentCopyNumber==m_CopyNumbers[sampleIndex]){ 
     
            toBeGenotypedThisTime.push_back(sampleIndex);    
            alreadyGenotyped[sampleIndex]=1;  

      }//  end determination of status of sample indexed by "sampleIndex".   
    }// end determination of the indices of all samples having the same allele pair and copy 
     // number value as the sample indexed by "count".



    //  We now find priors, if they exist, and insert them into the parameters used by bayes_label. 
    //  We submit all those samples whose index is in the vector "toBeGenotypedThisTime".
    //  We will use priors for the allele pair values of allele1, allele2 and copy number value of presentCopyNumber. 

    int numberToBeGenotyped = toBeGenotypedThisTime.size();
    snp_param tsp;
    tsp.copy(sp);
    tsp.copynumber = presentCopyNumber;


    ///@todo  Put prior determination in a function within QuantLabelZ to enforce conformity between QLZ and QLZM.  

    // The following is basically a cut and past from QuantLabelZ, but there are differences.
    // The names for the priors are assumed to take the following form: marker-allele1-allele2,   
    // with the addition of :1 if the copynumber is 1. If no priors can be found
    // in that form we look for "Generic:1" in case copynumber is 1.  If CN=0 we use priors for CN=2, with 
    // the same name as for the priors for CN=2. 

    string TmpName; 
    TmpName = m_ProbesetName + "-" + AffxString::intToString(allele1,false) + "-" 
        + AffxString::intToString(allele2,false);

    if (presentCopyNumber==1)
      TmpName = TmpName + ":" +ToStr(presentCopyNumber);

    // Search for the prior which has the correct name.
    snp_labeled_distribution sDist;
    getPrior(TmpName, presentCopyNumber, tsp, sDist);
    //  Prior determination is now complete.

    //  Here is where we make the final preparation for submission to bayes_label. Data structures are cleaned.

    // This is not used in QuantLabelZMulti, but we need the dummy input to satisfy the function signature. 
    vector<int> GenoHint;
    GenoHint.assign(m_numberOfSamples, affx::NN);

    // This is not supplied in QuantLabelZMulti, but we need the dummy input to satisfy the function signature.
    // These should be made compatible, because there is no reason for these capabilities not to be needed in general
    vector<double> SubsetInbredPenalty;
    // null vector, no need to do anything
    // length zero vectors are ignored
    SubsetInbredPenalty.clear();

    // We erase the input and output from the previous genotyping call.
    // Output    
    tcall.erase(tcall.begin(), tcall.end());
    tconf.erase(tconf.begin(), tconf.end());
    tcall.resize(numberToBeGenotyped);
    tconf.resize(numberToBeGenotyped);
    // Input
    tx.erase(tx.begin(), tx.end());
    ty.erase(ty.begin(), ty.end());

    // We transfer the strength/contrast values of the samples we are genotyping into the vectors tx, ty.  
    vector<int>::iterator vecIterBegin = toBeGenotypedThisTime.begin(); 
    vector<int>::iterator vecIterEnd = toBeGenotypedThisTime.end(); 
    for(; vecIterBegin != vecIterEnd; vecIterBegin++){
      tx.push_back(firstAlleleVector[*vecIterBegin]);       
      ty.push_back(secondAlleleVector[*vecIterBegin]);   
    }
 
    // Do the genotyping call.
    bayes_label(tcall,tconf,tx,ty,GenoHint,SubsetInbredPenalty,tsp);

    // Now amalgamate the tcall tconf and tsp values into the data structure m_CallsMulti m_Confidences and 
    // m_Distance for reporting. 
    int tempIndex=0; 
    int sampleCopyNumber=0;
    int brlmmpCall=0; 
    int call=0; 
    vecIterBegin = toBeGenotypedThisTime.begin(); 
    vecIterEnd = toBeGenotypedThisTime.end(); 
    for(; vecIterBegin!=vecIterEnd; vecIterBegin++){
      // To  translate the brlmmp output to a call we need the copy number of the sample and the 
      // brlmmp call itself. 
      sampleCopyNumber = m_CopyNumbers[*vecIterBegin];
      brlmmpCall = tcall[tempIndex];
      // For HET on CN=1, make it a no-call and set the confidence to 1
      if(brlmmpCall == affx::AB && sampleCopyNumber==1){
          tconf[tempIndex] = 1.0;
          brlmmpCall = affx::NN;
      }
      call = translateBrlmmpOutputToCall(sampleCopyNumber, brlmmpCall, allele1, allele2);  

      m_Calls[*vecIterBegin] = brlmmpCall;
      m_CallsMulti[*vecIterBegin] = call;
      m_Confidences[*vecIterBegin] = tconf[tempIndex];

      if(brlmmpCall == affx::AA)
        m_PriorObs[*vecIterBegin] = sDist.Dist.aa.k;
      else if(brlmmpCall == affx::AB)
        m_PriorObs[*vecIterBegin] = sDist.Dist.ab.k;
      else if(brlmmpCall == affx::BB)
        m_PriorObs[*vecIterBegin] = sDist.Dist.bb.k;
      else 
        m_PriorObs[*vecIterBegin] = 0;

      m_Distances[0][*vecIterBegin] =  abs((tsp.posterior.aa.m-tx[tempIndex])/sqrt(tsp.posterior.aa.ss));
      m_Distances[1][*vecIterBegin] =  abs((tsp.posterior.ab.m-tx[tempIndex])/sqrt(tsp.posterior.ab.ss)); 
      m_Distances[2][*vecIterBegin] =  abs((tsp.posterior.bb.m-tx[tempIndex])/sqrt(tsp.posterior.bb.ss));

      tempIndex++;
    }

    // now possible output of snp-specific posterior
    // that is, now we've updated the prior with the data
    // and we want to save what we learned
    if (m_SnpPosteriorTsv.is_open()) {
      bool report = true;
      if (m_ProbeSetsToReport != NULL) {
        if(m_ProbeSetsToReport->find(m_ProbesetName.c_str())==m_ProbeSetsToReport->end()) {
          report = false;
        }
        if (report) {
          writeSnpPosteriorValue(TmpName, tsp);
        }
      }
    }

    // Now increment the index of the sample we just genotyped.  If this sample has already been genotyped the continue
    // statement at the beginning of the while loop will ensure we set count to the smallest ungenotyped sample if
    // it exists.  If no ungenotyped sample exists the loop will exit. 
    count++;
  } // end while  
 
}
