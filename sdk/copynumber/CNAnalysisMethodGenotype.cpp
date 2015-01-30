///////////////////////////////////////////////////////////////
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
 * @file CNAnalysisMethodGenotype.cpp
 *
 * @brief This file contains the CNAnalysisMethodGenotype class members.
 */
#include "copynumber/CNAnalysisMethodGenotype.h"
//
#include "util/AffxStatistics.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodGenotype::explainSelf()
{
    CNAnalysisMethodGenotype obj;
    SelfDoc doc;
    doc.setDocName(obj.getType());
    doc.setDocDescription(obj.getDescription());
    doc.setDocOptions(obj.getDefaultDocOptions());
    return doc;
}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> CNAnalysisMethodGenotype::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)
/*
  example: 
  SelfDoc::Opt lohCN_errorrate = {"lohCN_errorrate", SelfDoc::Opt::Double, "0.05", "0.05", "NA", "NA", "LOH CN Error Rate"};
  opts.push_back(lohCN_errorrate);


*/
  return opts;
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
SelfCreate* CNAnalysisMethodGenotype::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodGenotype* pMethod = new CNAnalysisMethodGenotype();
    std::string strPrefix = getPrefix();
/* example
    pMethod->m_dLohCNErrorRate = setupDoubleParameter("lohCN_errorrate", strPrefix, params, doc, opts);
*/
    return pMethod;
}

//Constructor
CNAnalysisMethodGenotype::CNAnalysisMethodGenotype()
{
    m_iPrevGenderCode = -1;
    m_iGenderCode = -1;
    m_bLocalProbeSetsDetermined = false;
    m_iSeparation=1000000;
}



/**
 * @brief run the analysis
 */
void CNAnalysisMethodGenotype::run()
{
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Calculating SNP Data.");
    Verbose::out(1, "CNAnalysisMethodGenotype::run(...) start");
    isSetup();

    imputeUncalledCNstateForSNPs(getProbeSets());  // Assumes that probesets are Chromosome/Position ordered.

    getProbeSets()->quickSort(0);   // Sorted by ProbeSet Name.  Needed to load priors. 
    calculateGenotypes();

    if (m_pEngine->getEngineName() == "CNCytoEngine" && !m_pEngine->getOptBool("cyto2"))
    {
        getExperiment()->setCallRate(calculateCallRate());
    }

    getProbeSets()->quickSort(1);   // sorted by Chromosome then position

    Verbose::out(1, "CNAnalysisMethodGenotype::run(...) end");
}


float CNAnalysisMethodGenotype::calculateCallRate()
{
    unsigned int totalCount = 0;
    unsigned int callCount = 0;
    for (int iProbeSetIndex = 0; iProbeSetIndex < getProbeSets()->getCount(); iProbeSetIndex++)
    {
        // the three if clauses here (isProcess, SNP, confidence)
        // should sync up with if clauses in CNReporterCychp::writeGenInfoDataSet()
        CNProbeSet* p = (*getProbeSets())[iProbeSetIndex];
        if (!p->isProcess()) {
            continue;
        }
        if (p->processAsSNP())
        {
            totalCount++;
            float fConfidenceThreshold = CNAnalysisMethod::getConfidenceThreshold(m_pEngine->getOpt("brlmmp-parameters"));
            if (p->getGenotypeConfidence() >= fConfidenceThreshold) {
                // count as ALLELE_NO_CALL
                continue;
            }
            if (p->getGenotypeCallCode() != (char)ALLELE_NO_CALL)
            {
                callCount++;
            }
        }
    }
    return callCount/(float)totalCount;
}


void CNAnalysisMethodGenotype::imputeUncalledCNstateForSNPs(CNProbeSetArray* pvProbeSets)
{

    //  The algorithm is to advance up the genome until we hit the X chromosome.  
    //  We assume that the PAR region is a contiguous initial segment of the X chromosome 
    //  and we impute a CN call of 2 for all markers in this region.
    //


    //  Here we advance up to the X chromosome.
    CNProbeSet* pobjProbeSet = NULL;
    int iXChrBeginIndex=-1;

    int iIndex=0;
    for (; iIndex < pvProbeSets->getCount(); iIndex++)  
    {
        pobjProbeSet=pvProbeSets->getAt(iIndex);
        if(pobjProbeSet->getChromosome()==m_iXChromosome)
        {
            iXChrBeginIndex = iIndex;
            break;
        }
        
    }
   

    // Here we set the imputed CN state to 2 for all markers in the PAR region. 
    iIndex = iXChrBeginIndex;
    pobjProbeSet = pvProbeSets->getAt(iIndex);

    if(!(pobjProbeSet->isPseudoAutosomalRegion()))
    {
        Err::errAbort("PAR region is not at beginning of X chromosome as assumed by genotyping algorithm.");
    }

    while(pobjProbeSet->isPseudoAutosomalRegion())
    {
        pobjProbeSet->setImputedCNState(2);
        iIndex++;
        pobjProbeSet = pvProbeSets->getAt(iIndex);
    }

    //  Here iBeginIndex points to the first marker outside of the PAR region.
    int iBeginIndex = iIndex;


    // Here we advance up the genome and check whether the separation of adjacent markers is greater that the
    // distance we expect across the centromere.  If it is we we pause and impute CN states to markers up to iIndex.
    // Note that iIndex points to the last element in a block, so iBeginIndex is advanced one beyond this point.
    // We hit an end when iIndex is advanced to one beyond the X chromosome.  We then impute this last block outside the loop. 
    // One important thing to note is that we impute the upper and lower bounds of each block.
    while(pvProbeSets->getAt(iIndex)->getChromosome() == m_iXChromosome )
    {
        if ( (pvProbeSets->getAt(iIndex+1)->getPosition() - pvProbeSets->getAt(iIndex)->getPosition()) > m_iSeparation)
        {
              impute(iBeginIndex, iIndex, pvProbeSets); 
              iBeginIndex = iIndex+1;
        }    
        iIndex++;
    } 
    impute(iBeginIndex,iIndex-1, pvProbeSets);
}

void CNAnalysisMethodGenotype::impute(int iBeginIndex, int iEndIndex, CNProbeSetArray* pvProbeSets )
{



    //  We first advance up to the marker with the first valid CN call.
    int iIndex = iBeginIndex;
    while(iIndex <= iEndIndex  && pvProbeSets->getAt(iIndex)->getCNState()==-1)
    {
       iIndex++;
    }
    
    // iIndex now points to the first probeset that has a valid called CN state. 
    int iFirstValidCNCall=pvProbeSets->getAt(iIndex)->getCNState();

    //  Here we set all the markers up to and including the above first probeset to have an imputed
    //  call equal to the called state of that first marker 
    for(int iNewIndex=iBeginIndex; iNewIndex <= iIndex; iNewIndex++)
    {
        pvProbeSets->getAt(iNewIndex)->setImputedCNState(iFirstValidCNCall);
    }
   


    //  iPresentBeginIndex points to the marker with a valid CN call.
    int iPresentCNState= iFirstValidCNCall;
    int iPresentBeginIndex= iIndex;
        iIndex++;

    //  The idea of the algorithm is to keep advancing until we find a marker with a valid CN call. This is the second "if" below.
    //  We then impute a CN state to all markers up to this new marker with the last CN state. This is the "else" below.
    //  Note that at that point iIndex points to the new marker with a valid CN call, so we don't impute a call to it.
    //  At the end, we will advance until iIndex points to the last marker in the block we want to impute a CN call to.
    //  There are two cases.  Either it has an valid CN state or not.  If it does it is imputed with this value.  If not it
    //  get the imputed state equal to the previous valid CN state called.
    while(iIndex<=iEndIndex)
    {
        if(iIndex==iEndIndex)
        {
            for(int iNewIndex=iPresentBeginIndex; iNewIndex<= iEndIndex-1; iNewIndex++)
            {
               pvProbeSets->getAt(iNewIndex)->setImputedCNState(iPresentCNState);         
            }
            if(pvProbeSets->getAt(iIndex)->getCNState() != -1)
            {
                pvProbeSets->getAt(iIndex)->setImputedCNState(iPresentCNState);
            }
            break;
        } 
           
        int iNewCNState=pvProbeSets->getAt(iIndex)->getCNState();
        if(iNewCNState==-1 )
        {
            iIndex++;
            continue;
        }
        else
        {
            for(int iNewIndex=iPresentBeginIndex; iNewIndex<=iIndex-1; iNewIndex++)
            {
               pvProbeSets->getAt(iNewIndex)->setImputedCNState(iPresentCNState);         
            }
            iPresentCNState=pvProbeSets->getAt(iIndex)->getCNState(); 
            iPresentBeginIndex=iIndex;
            iIndex++; 
               
        }
    }
}


void CNAnalysisMethodGenotype::calculateGenotypes()
{
    //set("brlmmp-parameters", "CM=2.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.wobble=0.05.MS=1")

    bool bCyto2 = ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")));
    snp_param sp;
    initializeSnpParameters(sp, 2);
    setBrlmmpParameters(sp);

    loadSnpDistributions();

    vector<affx::GType> tcall(1);
    vector<double> tconf(1);
    vector<double> tx(1);
    vector<double> ty(1);
    vector<int> GenoHint(1);
    vector<double> SubsetInbredPenalty(1);
    for (int iProbeSetIndex = 0; (iProbeSetIndex < getProbeSets()->getCount()); iProbeSetIndex++)
    {
        CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
        if ((bCyto2) && (p->getGenotypeCall() != -1)) {
            continue;
        }
        p->setGenotypeCall(-1);
        p->setGenotypeConfidence(0);
        if (p->getSnpDistribution() != NULL)
        {
            snp_param tsp;
            tsp.copy(sp); // duplicate current parameters to keep commands
            tsp.copynumber = 2; // what copy for this snp so know right model centers
            if (getExperiment()->getRawIntensityRatioGenderAsInt() == affx::Male)
            {
                if ((p->getChromosome() == m_iXChromosome) && (!p->isPseudoAutosomalRegion() && p->getImputedCNState() == 1)) 
                {
                    tsp.copynumber = 1;
                }
                if ((p->getChromosome() == m_iXChromosome) && (!p->isPseudoAutosomalRegion() && p->getImputedCNState() >= 2 )) 
                {
                    tsp.copynumber = 2;
                }
                if ((p->getChromosome() == m_iXChromosome) && (!p->isPseudoAutosomalRegion() && p->getImputedCNState() <  1)) 
                {
                    tsp.copynumber = 2;
                }
                else // Y Chromosome, Autosome or Par 
                {
                        if (p->getChromosome() == m_iYChromosome) 
                        {
                                tsp.copynumber = 1;
                        }
                }
            }
            else  // Female
            {
                if (p->getChromosome() == m_iYChromosome) {
                    continue;
                }
            }
            tsp.prior.Copy(*p->getSnpDistribution());
            // note: using temporary copy of main parameters tsp
            tcall[0] = affx::NN;
            tconf[0] = 0;
            tx[0] = p->getAAlleleSignal();
            ty[0] = p->getBAlleleSignal();
            GenoHint[0] = -1;
            SubsetInbredPenalty[0] = 0; // null parameter
            GenoUtility_transformData(tx, ty, MvA, 1);
            bayes_label(tcall,tconf,tx,ty,GenoHint,SubsetInbredPenalty,tsp);
            p->setGenotypeCall(affx::GType_to_int(tcall[0]));
            p->setGenotypeConfidence((float)tconf[0]);
        }
    }
}


void CNAnalysisMethodGenotype::initializeSnpParameters(snp_param& sp, int iCallMethod)
{
    sp.Initialize();
    sp.callmethod = iCallMethod;
    sp.bins = 100;
    sp.mix = 1;
    sp.bic = 2;
    sp.hardshell = 3;
    sp.shellbarrier = 0.45;

    sp.prior.aa.k = 1;
    sp.prior.ab.k = 1.5;
    sp.prior.bb.k = 1; // bb.k always should = aa.k

    sp.prior.aa.m = 2;
    sp.prior.ab.m = 0;
    sp.prior.bb.m = -2;

    sp.prior.aa.ss = 0.06;
    sp.prior.ab.ss = 0.06;
    sp.prior.bb.ss = 0.06;

    sp.prior.aa.v = 10;
    sp.prior.ab.v = 10;
    sp.prior.bb.v = 10;

    sp.prior.xah = -0.6;
    sp.prior.xab = 0.5;
    sp.prior.xhb = -0.6;

    sp.copyqc = 0;
    sp.wobble = 0.05;
}

void CNAnalysisMethodGenotype::setBrlmmpParameters(snp_param& sp)
{
    string spec = m_pEngine->getOpt("brlmmp-parameters");
    if (spec.empty()) {
        return;
    }
    string name;
    map<string, string> param;
    SelfCreate::fillInNameParam(spec, name, param);

    map<string, string>::iterator iter;
    if ((iter = param.find("bins")) != param.end()) {
        sp.bins = getInt(iter->second);
    }
    if ((iter = param.find("mix")) != param.end()) {
        sp.mix = getInt(iter->second);
    }
    if ((iter = param.find("bic")) != param.end()) {
        sp.bic = getDouble(iter->second);
    }
    if ((iter = param.find("HARD")) != param.end()) {
        sp.hardshell = getInt(iter->second);
    }
    if ((iter = param.find("SB")) != param.end()) {
        sp.shellbarrier = getDouble(iter->second);
    }

    if ((iter = param.find("KX")) != param.end()) {
        sp.prior.aa.k = getDouble(iter->second);
    }
    if ((iter = param.find("KH")) != param.end()) {
        sp.prior.ab.k = getDouble(iter->second);
    }
    sp.prior.bb.k = sp.prior.aa.k; // bb.k always should = aa.k

    if ((iter = param.find("AAM")) != param.end()) {
        sp.prior.aa.m = getDouble(iter->second);
    }
    //sp.prior.ab.m = 0;    // set by initializeSnpParameters() already
    if ((iter = param.find("BBM")) != param.end()) {
        sp.prior.bb.m = getDouble(iter->second);
    }

    if ((iter = param.find("AAV")) != param.end()) {
        sp.prior.aa.ss = getDouble(iter->second);
    }
    if ((iter = param.find("ABV")) != param.end()) {
        sp.prior.ab.ss = getDouble(iter->second);
    }
    if ((iter = param.find("BBV")) != param.end()) {
        sp.prior.bb.ss = getDouble(iter->second);
    }
    // done by initializeSnpParameters() already
    //sp.prior.aa.v = 10.0;
    //sp.prior.ab.v = 10.0;
    //sp.prior.bb.v = 10.0;

    if ((iter = param.find("KAH")) != param.end()) {
        sp.prior.xah = getDouble(iter->second);
    }
    if ((iter = param.find("KXX")) != param.end()) {
        sp.prior.xab = getDouble(iter->second);
    }
    if ((iter = param.find("KHB")) != param.end()) {
        sp.prior.xhb = getDouble(iter->second);
    }

    if ((iter = param.find("copyqc")) != param.end()) {
        sp.copyqc = getDouble(iter->second);
    }
    if ((iter = param.find("wobble")) != param.end()) {
        sp.wobble = getDouble(iter->second);
    }
}
void CNAnalysisMethodGenotype::loadSnpDistributions()
{
        AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
        affx::File5_File file5;
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;

        file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);

        if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
        {
            group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            tsv5 = group5->openTsv("SNPPosteriors", affx::FILE5_OPEN);
        }
        else
        {
            group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
            tsv5 = group5->openTsv("CN5.snp-posteriors", affx::FILE5_OPEN);
        }
        QuantLabelZ__readSnpPriorMap_tsv5_v2(m_iXChromosome, m_iYChromosome, *getExperiment(), *getProbeSets(), tsv5);
        tsv5->close();
        delete tsv5;
        group5->close();
        delete group5;
        file5.close();
}

