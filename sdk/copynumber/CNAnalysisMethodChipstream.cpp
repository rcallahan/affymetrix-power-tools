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
 * @file CNAnalysisMethodChipstream.cpp
 *
 * @brief This file contains the CNAnalysisMethodChipstream class members.
 */

#include "copynumber/CNAnalysisMethodChipstream.h"
//
#include "copynumber/CNAnalysisMethodCovariateSignalAdjuster.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNAnalysisMethodReference.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/BioTypes.h"
#include "chipstream/CelReader.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStreamFactory.h"
#include "chipstream/HomHiLoCelListener.h"
#include "chipstream/QuantLabelZIO.h"
#include "chipstream/SketchQuantNormTran.h"
#include "chipstream/SparseMart.h"
#include "file/TsvFile/TsvFile.h"
#include "normalization/normalization.h"
#include "util/AffxStatistics.h"
#include "util/RowFile.h"
//
/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodChipstream::explainSelf()
{
    CNAnalysisMethodChipstream obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodChipstream::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)
  SelfDoc::Opt chipstreamNDBandwidth = {"NDBandwidth", SelfDoc::Opt::Float, "10.0", "10.0", "10.0", "20.0", "Normal-Diploid Bandwidth"};
  opts.push_back(chipstreamNDBandwidth);



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
SelfCreate* CNAnalysisMethodChipstream::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodChipstream* pMethod = new CNAnalysisMethodChipstream();
    std::string strPrefix = getPrefix();

    pMethod->m_fNDBandwidth = setupFloatParameter("NDBandwidth", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodChipstream::CNAnalysisMethodChipstream()
{
    m_iPrevGenderCode = -1;
    m_iGenderCode = -1;

    m_uiMaxProbeCount = 0;
    m_pdProbeEffects = NULL;
    m_pdChipEffects = NULL;
    m_ppdPM = NULL;
    m_ppdMM = NULL;
    m_ppdResiduals = NULL;
    m_pPlier = NULL;
    m_iPlierChipCount = 0;

    m_pvProbes = NULL;
    m_pProbes = NULL;
    m_iInstanceCount = 0;
    m_pobjAdapterTypeNormTran = NULL;

    m_pvSNPReferenceSketch = NULL;
    m_pvCNReferenceSketch = NULL;
    m_pvChipSketch = NULL;
    m_pvSNPData = NULL;
    m_pvCNData = NULL;

    m_dXXRatio = 0;
    m_dYRatio = 0;

    m_fNDBandwidth = 0.0;
}

/**
 * @brief Destructor
 */
CNAnalysisMethodChipstream::~CNAnalysisMethodChipstream()
{
    if (m_pvSNPReferenceSketch != NULL) {deleteSketch();}
    if (m_pProbes != NULL) {deleteProbes();}
    if (m_pPlier != NULL) {deletePlier();}
    if (m_pvSNPData != NULL) {delete m_pvSNPData;}
    if (m_pvCNData != NULL) {delete m_pvCNData;}
}


/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodChipstream::run()
{
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Running Probe Intensity Normalization.");

    Verbose::out(1, "CNAnalysisMethodChipstream::run(...) start");
    isSetup();
    getProbeSets()->quickSort(0); // by ProbeSetName

//    loadGenotypeCallsFromFile();
    if (m_pProbes == NULL) {newProbes();}
    if (m_pvSNPReferenceSketch == NULL) {newSketch();}
    if (m_pPlier == NULL) {newPlier(1);}
    determineSketchProbeSets();
    normalizeIntensities();
    if (m_pEngine->getOptBool("keep-intermediate-data"))
        {
                writeIntensities("single-sample-normalized");
        }

    memory("CNAnalysisMethodChipstream::run(...)");

    if (        m_pEngine->isOptDefined("local-gc-background-intensity-adjustment-method")      && 
                m_pEngine->getOpt("local-gc-background-intensity-adjustment-method") != "none"  &&
                m_pEngine->getOpt("local-gc-background-intensity-adjustment-method") != "")
    {
        adjustIntensitiesPDNN();
    }

    adjustIntensitiesHighPassFilter();

    calculateSignalEstimates(); 
    calculateXY(); 

    if (m_pEngine->getOptBool("keep-intermediate-data"))
    {
            writeSignals("chipstream");
            writeSignalsTsv("chipstream", getProbeSets());
    }

    getProbeSets()->quickSort(1); // by Chromosome, Position
    Verbose::out(1, "CNAnalysisMethodChipstream::run(...) end");
}


/*
void CNAnalysisMethodChipstream::loadSnpProbeSets(std::vector<CNProbeSet*> &vSnpProbeSets)
{
        vSnpProbeSets.clear();
        for (int i = 0; (i < m_pvProbeSets->getCount()); i++)
        {
                CNProbeSet* p = m_pvProbeSets->getAt(i);
                if ((p->isProcess()) && (p->processAsSNP()) && (p->getInformation() != 0))
                {
                        vSnpProbeSets.push_back(p);
                }
        }
}
*/


void CNAnalysisMethodChipstream::determineSketchProbeSets()
{
//MG temporary code for Bitao
        AffxArray<AffxString> ar;
        if ((m_pEngine->getOpt("normal-diploid-files-file") != "") && (m_pEngine->getOpt("normal-diploid-files-file") != "none"))
        {
          affx::TsvFile tsv;
          tsv.m_optAutoTrim = true;

          std::string diploid;
          AffxString normalDiploidFileName = getExperiment()->getExperimentNormalDiploidFileName();
          if (tsv.open(normalDiploidFileName) == affx::TSV_OK)  {
            while ( tsv.nextLevel(0) == affx::TSV_OK ) {
              tsv.get(0,0, diploid);
              ar.add(new AffxString(diploid));
            }
            tsv.clear();
          } else {Err::errAbort("Cannot open file: " + normalDiploidFileName);}
          ar.quickSort(0);
        }
//
        for (int iIndex = 0; (iIndex < getProbeSets()->getCount()); iIndex++)
        {
                CNProbeSet* p = getProbeSets()->getAt(iIndex);
                p->setUseForSketch(false);
                if (p->getChromosome() < m_iXChromosome)
                {
                        if (ar.getCount() > 0)
                        {
                                // Check if probe set exists in normal diploid list.
                                AffxString str = p->getProbeSetName();
                                p->setUseForSketch(ar.binarySearch(str, 0) != -1);
                        }
                        else
                        {
                                // No list, assume all probe sets are normal diploid.
                                p->setUseForSketch(true);
                        }
                }
        }
        ar.deleteAll();
}

void CNAnalysisMethodChipstream::extractSketch(std::vector<float>& vData, std::vector<float>& vSketch)
{
    AffxMultiDimensionalArray<float> v(vData.size());
    for (int i = 0; (i < vData.size()); i++)
    {
        v.set(i, vData[i]);
    }
    v.quickSort();
    for (int i = 0; (i < vSketch.size()); i++)
    {
        vSketch[i] = (float)v.percentile(((double)i / (double)(vSketch.size() - 1)), false);
    }
}

float CNAnalysisMethodChipstream::interpolate(    float fIntensity,
                                                AffxMultiDimensionalArray<float>& vReferenceSketch,
                                                AffxMultiDimensionalArray<float>& vChipSketch,
                                                double dRatio)
{
    unsigned int uiSketchSize = vReferenceSketch.length();
    if (fIntensity < vChipSketch.get(0))
    {
        fIntensity = vReferenceSketch.get(0) * fIntensity / vChipSketch.get(0);
    }
    else if (fIntensity >= vChipSketch.get(uiSketchSize - 1))
    {
        fIntensity = vReferenceSketch.get(uiSketchSize - 1) + (fIntensity - vChipSketch.get(uiSketchSize - 1)) * dRatio;
    }
    else
    {
        unsigned int uiLowIndex = 0;
        unsigned int uiHighIndex = 0;
        vChipSketch.binarySearch(fIntensity, uiLowIndex, uiHighIndex); // Produces <= n <=
        if (fIntensity == vChipSketch.get(uiHighIndex)) {uiLowIndex = uiHighIndex; uiHighIndex++;} // Want <= n <
        if (vChipSketch.get(uiLowIndex) == vChipSketch.get(uiHighIndex))
        {
            fIntensity = vReferenceSketch.get(uiLowIndex);
        }
        else
        {
            fIntensity = vReferenceSketch.get(uiLowIndex) + (fIntensity - vChipSketch.get(uiLowIndex)) * ((vReferenceSketch.get(uiHighIndex) - vReferenceSketch.get(uiLowIndex)) / (vChipSketch.get(uiHighIndex) - vChipSketch.get(uiLowIndex)));
        }
    }
    if (!Util::isFinite(fIntensity)) {fIntensity = vReferenceSketch.get(vReferenceSketch.length() - 1);}
    return fIntensity;
}

void CNAnalysisMethodChipstream::normalizeIntensities()
{
    // Check to see if a subset was selected.
    loadIntensities();

    if (!(m_pEngine->getOptBool("doDualNormalization")))
    {        
        unsigned int uiSketchSize = 0;
        if(m_pEngine->isOptDefined("sketch-size"))
        {
            uiSketchSize = m_pEngine->getOptInt("sketch-size");
        }
        else
        {
            uiSketchSize = 50000;
        }
        if (uiSketchSize > m_pvProbes->getCount()) {uiSketchSize = m_pvProbes->getCount();}
        SketchQuantNormTran *qnorm = new SketchQuantNormTran(uiSketchSize, false, false, false, 0, false);
        vector<float> sketch;
        AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
        affx::File5_File file5;
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;
        file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
        group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        double d = 0;
        tsv5 = group5->openTsv("target-sketch", affx::FILE5_OPEN);
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            tsv5->get(0, 0, &d); sketch.push_back((float)d);
        }
        tsv5->close();
        delete tsv5;
        group5->close();
        delete group5;
        file5.close();

        qnorm->setTargetSketch(sketch);

        qnorm->newChip(*m_pvSNPData);
        for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
        {
            CNProbe* p = m_pvProbes->getAt(iIndex);
            p->setIntensity(qnorm->transform(iIndex, 0, p->getIntensity()));
        }
        delete qnorm;
        return;
    }

    extractSketch(*m_pvSNPData, *m_pvChipSketch);

    AffxMultiDimensionalArray<float> vReferenceSketch(m_pvSNPReferenceSketch->size());
    AffxMultiDimensionalArray<float> vChipSketch(m_pvChipSketch->size());
    for (int iIndex = 0; (iIndex < m_pvChipSketch->size()); iIndex++)
    {
        vReferenceSketch.set(iIndex, m_pvSNPReferenceSketch->at(iIndex));
        vChipSketch.set(iIndex, m_pvChipSketch->at(iIndex));
    }
    double dRatio = (vReferenceSketch.median() / vChipSketch.median());

    // Normalize the snp probe intensities.
    for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
    {
        CNProbe* p = m_pvProbes->getAt(iIndex);
        if (p->getProbeSetIndex() == -1) {continue;}
        if (getProbeSets()->getAt(p->getProbeSetIndex())->processAsSNPNormalize())
        {
            p->setIntensity(interpolate(p->getIntensity(), vReferenceSketch, vChipSketch, dRatio));
        }
    }

    extractSketch(*m_pvCNData, *m_pvChipSketch);

    for (int iIndex = 0; (iIndex < m_pvChipSketch->size()); iIndex++)
    {
        vReferenceSketch.set(iIndex, m_pvCNReferenceSketch->at(iIndex));
        vChipSketch.set(iIndex, m_pvChipSketch->at(iIndex));
    }
    dRatio = (vReferenceSketch.median() / vChipSketch.median());

    // Normalize the cn probe intensities.
    for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
    {
        CNProbe* p = m_pvProbes->getAt(iIndex);
        if (p->getProbeSetIndex() == -1) {continue;}
        if (getProbeSets()->getAt(p->getProbeSetIndex())->processAsCNNormalize())
        {
            p->setIntensity(interpolate(p->getIntensity(), vReferenceSketch, vChipSketch, dRatio));
        }
    }

        bool testFlag = false;
        if(testFlag)
        {
                affx::TsvFile *tsv = new affx::TsvFile;
                tsv->writeTsv(getExperiment()->getExperimentName() + ".intensities.txt");
                tsv->defineColumn(    0,0,"probe_id", affx::TSV_TYPE_INT);
                tsv->defineColumn(    0,1,"intensity", affx::TSV_TYPE_DOUBLE);

            for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
            {
                CNProbe* p = m_pvProbes->getAt(iIndex);
                tsv->set(0,0,p->getProbeID());
                tsv->set(0,1,p->getIntensity());
                tsv->writeLevel(0);
            }
            delete tsv;
        }
}

void CNAnalysisMethodChipstream::calculateSignalEstimates()
{
//    Verbose::out(1, "CNAnalysisMethodChipstream::calculateSignalEstimates()");
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Running Signal Summarization.");
    m_pvProbes->quickSort(0); // by ProbeSetIndex, by Allele
    AffxMultiDimensionalArray<float> vIntensities(m_uiMaxProbeCount);

    long errorCode = 0;
    CNProbe* pPrev = NULL;
    unsigned int uiProbeCount = 0;
    for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
    {
        CNProbe* p = m_pvProbes->getAt(iIndex);
        if ((pPrev == NULL) || (pPrev->compareTo(*p, 0) == 0))
        {
            if (pPrev == NULL) {pPrev = p;}
            m_pdProbeEffects[uiProbeCount] = p->getProbeEffect();
            m_ppdPM[0][uiProbeCount] = p->getIntensity();
            vIntensities.set(uiProbeCount, p->getIntensity());
            m_ppdMM[0][uiProbeCount] = 0;
            m_ppdResiduals[0][uiProbeCount] = 0;
            uiProbeCount++;
        }
        else
        {
            // Run plier
            errorCode = 0;
            m_pdChipEffects[0] = 0;
            m_pPlier->setNumFeature(uiProbeCount);
            m_pPlier->run(&errorCode);
            if (errorCode != 0) {Err::errAbort("Problem running plier. Error code: " + ToStr(errorCode));}
            setProbeSetSignal(pPrev, m_pdChipEffects[0], vIntensities.median(uiProbeCount));

            // Setup for next plier run
            pPrev = p;
            m_pdProbeEffects[0] = p->getProbeEffect();
            m_ppdPM[0][0] = p->getIntensity();
            vIntensities.set(0, p->getIntensity());
            m_ppdMM[0][0] = 0;
            m_ppdResiduals[0][0] = 0;
            uiProbeCount = 1;
        }
    }
    if (pPrev != NULL)
    {
        // Run plier
        errorCode = 0;
        m_pdChipEffects[0] = 0;
        m_pPlier->setNumFeature(uiProbeCount);
        m_pPlier->run(&errorCode);
        if (errorCode != 0) {Err::errAbort("Problem running plier. Error code: " + ToStr(errorCode));}
        setProbeSetSignal(pPrev, m_pdChipEffects[0], vIntensities.median(uiProbeCount));
    }
}

void CNAnalysisMethodChipstream::setProbeSetSignal(CNProbe* p, double dSignal, double dMedianIntensity)
{
    if (p->getAllele() == 0)
    {
        float fSignal = (float)dSignal;
        getProbeSets()->getAt(p->getProbeSetIndex())->setAAlleleSignal((float)(fSignal / 2.0));
        getProbeSets()->getAt(p->getProbeSetIndex())->setBAlleleSignal((float)(fSignal / 2.0));
        getProbeSets()->getAt(p->getProbeSetIndex())->setAMedianIntensity((float)(dMedianIntensity / 2.0));
        getProbeSets()->getAt(p->getProbeSetIndex())->setBMedianIntensity((float)(dMedianIntensity / 2.0));
    }
    else
    {
        if (p->getAllele() == 'A')
        {
            getProbeSets()->getAt(p->getProbeSetIndex())->setAAlleleSignal((float)dSignal);
            getProbeSets()->getAt(p->getProbeSetIndex())->setAMedianIntensity((float)dMedianIntensity);
        }
        else if (p->getAllele() == 'B')
        {
            getProbeSets()->getAt(p->getProbeSetIndex())->setBAlleleSignal((float)dSignal);
            getProbeSets()->getAt(p->getProbeSetIndex())->setBMedianIntensity((float)dMedianIntensity);
        }
    }
}

void CNAnalysisMethodChipstream::calculateGenotypes()
{
    //set("brlmmp-parameters", "CM=2.bins=100.mix=1.bic=2.HARD=3.SB=0.45.KX=1.KH=1.5.KXX=0.5.KAH=-0.6.KHB=-0.6.transform=MVA.AAM=2.0.BBM=-2.0.AAV=0.06.BBV=0.06.ABV=0.06.copyqc=0.wobble=0.05.MS=1")

    bool bCyto2 = ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")));
    snp_param sp;
    initializeSnpParameters(sp, 2);
    setBrlmmpParameters(sp);
    if (!bCyto2) {loadSnpDistributions();}
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
                if ((p->getChromosome() == m_iXChromosome) && (!p->isPseudoAutosomalRegion())) {tsp.copynumber = 1;}
                else if (p->getChromosome() == m_iYChromosome) {tsp.copynumber = 1;}
            }
            else
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

void CNAnalysisMethodChipstream::newProbes()
{
    AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;

    int iCount = 0;
    if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
    {
        iCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "ProbeEffects");
    }
    else
    {
        iCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CN5", "CN5.plier.feature-response");
    }
    if (iCount <= 0) {throw(Except("Cannot find any probes in the CN reference file to process."));}
    file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
    if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
    {
        group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
    }
    else
    {
        group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
    }

    m_pProbes = new CNProbe[iCount];
    int i = 0;
    double d = 0;
    float f = 0;
    AffxString str;
    CNProbeSet ps;
    if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
    {
        tsv5 = group5->openTsv("ProbeEffects", affx::FILE5_OPEN);
        if (tsv5->getColumnCount(0) < 5) {Err::errAbort("The CN Reference file is missing data required by CNAnalysisMethodChipstream::newProbes() ");}
    }
    else
    {
        tsv5 = group5->openTsv("CN5.plier.feature-response", affx::FILE5_OPEN);
    }
    unsigned int uiIndex = 0;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        tsv5->get(0, 0, &i); m_pProbes[uiIndex].setProbeID(i);
        tsv5->get(0, 1, &d); m_pProbes[uiIndex].setProbeEffect(d);

        if (tsv5->getColumnDtype(0,2) == affx::FILE5_DTYPE_INT)
        {
            tsv5->get(0, 2, &i); m_pProbes[uiIndex].setProbeSetIndex(i);
            tsv5->get(0, 3, &i); m_pProbes[uiIndex].setAllele((char)i);
            tsv5->get(0, 4, &f); m_pProbes[uiIndex].setMedianIntensity(f);
            uiIndex++;
        }
        else
        {
            tsv5->get(0, 2, &str);
            if ((str.endsWith("-A")) || (str.endsWith("-B")))
            {
                str = str.substring(0, str.length() - 2);
            }
            ps.setProbeSetName(str);
            int iSearchIndex = m_pvProbeSets->binarySearch(ps, 0);
            if (iSearchIndex != -1)
            {
                m_pProbes[uiIndex].setProbeSetIndex(iSearchIndex);
                uiIndex++;
            } else {iCount--;}
        }
    }
    tsv5->close();
    delete tsv5;

    group5->close();
    delete group5;
    file5.close();

    m_pvProbes = new CNProbeArray;
    m_pvProbes->reserve(iCount);
    for (int iProbeIndex = 0; (iProbeIndex < iCount); iProbeIndex++)
    {
        m_pvProbes->add(&m_pProbes[iProbeIndex]);
    }

    // Initialize the adapter type normalization
    if (!m_pEngine->getOptBool("cyto2"))
    {
        ChipLayout* pLayout = new ChipLayout;
        CNAnalysisMethodReference::loadLayout(pLayout,m_pEngine);
        for (int iProbeIndex = 0; (iProbeIndex < iCount); iProbeIndex++)
        {
            CNProbe* probe = m_pvProbes->at(iProbeIndex);
            if (probe->getProbeSetIndex() == -1) {continue;}
            AffxString strProbeSetName = m_pvProbeSets->at(probe->getProbeSetIndex())->getProbeSetName();

            ProbeListPacked pList = pLayout->getProbeListByName(strProbeSetName);
            if (!pList.isNull())
            {
                ProbeSet* ps = ProbeListFactory::asProbeSet(pList);
                if ((ps->psType == ProbeSet::GenoType) || (ps->psType == ProbeSet::Copynumber) || (strProbeSetName.startsWith("anti")))
                {
                    int iAtomCount = 0;
                    for (unsigned int iGroupIndex = 0; (iGroupIndex < ps->numGroups); iGroupIndex++)
                    {
                        for (unsigned int iAtomIndex = 0; (iAtomIndex < ps->atomsPerGroup[iGroupIndex]); iAtomIndex++)
                        {
                            for (unsigned int iProbeIndex = 0; (iProbeIndex < ps->atoms[iAtomCount]->probes.size()); iProbeIndex++)
                            {
                                if (Probe::isPm(*ps->atoms[iAtomCount]->probes[iProbeIndex]))
                                {
                                    int uiIndex = ps->atoms[iAtomCount]->probes[iProbeIndex]->id;
                                    if (probe->getProbeID() == (uiIndex + 1))
                                    {
                                        probe->setAllele(0);
                                        if (ps->psType == ProbeSet::GenoType)
                                        {
                                            if ((iGroupIndex % 2) == 0)
                                            {
                                                probe->setAllele('A');
                                            }
                                            else
                                            {
                                                probe->setAllele('B');
                                            }
                                        }
                                    }
                                }
                            }
                            iAtomCount++;
                        }
                    }
                }
                delete ps;
            }
        }
        if ((m_pEngine->isOptDefined("adapter-type-normalization")) && (m_pEngine->getOptBool("adapter-type-normalization")))
        {
            m_pobjAdapterTypeNormTran = new AdapterTypeNormTran(0.8287, 0.7960, 1.4218, 0.9954, 1.63920);
            m_pobjAdapterTypeNormTran->setup(m_pEngine->getOpt("annotation-file"), pLayout, 1);
        }
        calculateRawSNPQCs(pLayout);

        delete pLayout;
    }
}

void CNAnalysisMethodChipstream::calculateRawSNPQCs(ChipLayout* pLayout)
{
    AffxArray<AffxString> arSNPIds;
    if ((m_pEngine->isOptDefined("snp-qc-snp-list")) && (m_pEngine->getOpt("snp-qc-snp-list") != ""))
    {
      affx::TsvFile tsv;
      tsv.m_optAutoTrim = true;
      if ( tsv.openTable(m_pEngine->getOpt("snp-qc-snp-list")) == affx::TSV_OK) {
        std::string snp;
        while (tsv.nextLevel(0) == affx::TSV_OK) {
          tsv.get(0,0,snp);
          arSNPIds.add(new AffxString(snp));
        }
        tsv.clear();
      }
    }
    arSNPIds.quickSort(0);

    vector<string> probeSetNames;
    if (arSNPIds.empty()) {
        for (int iIndex = 0; iIndex < m_pvProbeSets->getCount(); iIndex++) {
            CNProbeSet* pobjProbeSet = m_pvProbeSets->at(iIndex);
            if (pobjProbeSet->processAsVisualization())
            {
                probeSetNames.push_back(pobjProbeSet->getProbeSetName());
            }
        }
    } else {
        for (int iIndex = 0; iIndex < m_pvProbeSets->getCount(); iIndex++) {
            CNProbeSet* pobjProbeSet = m_pvProbeSets->at(iIndex);
            if (pobjProbeSet->processAsVisualization())
            {
                AffxString str = pobjProbeSet->getProbeSetName();
                if (arSNPIds.binarySearch(str, 0) != -1)
                {
                    probeSetNames.push_back(pobjProbeSet->getProbeSetName());
                }
            }
        }
    }
    arSNPIds.deleteAll();

    const double log2 = log(2.0);
    for (int iIndex = 0; iIndex < getExperiments()->size(); iIndex++)
    {
        CNExperiment* pExperiment = getExperiments()->getAt(iIndex);
        vector<double> vData;

        vector<string> celFileNames;
        string cel = (m_pEngine->getOptVector("cels"))[pExperiment->getIndex()];

        celFileNames.push_back(cel);
        vector<int> dummy_order(pLayout->getProbeCount());
        for (int i = 0; i < pLayout->getProbeCount(); i++) {
            dummy_order[i] = i;
        }
        ChipStream *chipStream = 0;
        string junk;
        ChipStreamFactory chipStreamFactory;
        chipStream = chipStreamFactory.chipStreamForString("no-trans", *pLayout, junk);
        SparseMart *sparseMart = new SparseMart(dummy_order, celFileNames, true);
        dummy_order.clear();
        vector<int>().swap(dummy_order);   // release memory

        CelReader celReader;
        celReader.setFiles(celFileNames);
        celReader.registerStream(chipStream);
        celReader.registerIntensityMart(sparseMart);

        celReader.readFiles();

        vector<int> probeIds;
        vector<double> probeIntensities;
        for (vector<string>::iterator it = probeSetNames.begin(); it != probeSetNames.end(); ++it) {
            probeIds.clear();
            probeIntensities.clear();
            ProbeListPacked probeset = pLayout->getProbeListByName(*it);
            if (probeset.isNull()) {
                Verbose::warn(0, "Empty probe set '" + *it + "' found in cel file " + cel +
                                 " while calculating SNPQC. This probe set will be ignored.");
                continue;
            }
            probeset.get_probeIds(probeIds);

            for (int idx=0; idx < probeIds.size(); idx++) {
                probeIntensities.push_back(chipStream->getTransformedIntensity(probeIds[idx], 0));
            }
            vector<double>::iterator iterMedianA = probeIntensities.begin() + probeIntensities.size()/4;
            vector<double>::iterator iterEndA    = probeIntensities.begin() + probeIntensities.size()/2;

            vector<double>::iterator iterMedianB = probeIntensities.begin() + (3*probeIntensities.size())/4;
            vector<double>::iterator iterStartB  = probeIntensities.begin() + probeIntensities.size()/2;

            std::nth_element(probeIntensities.begin(), iterMedianA, iterEndA);
            std::nth_element(iterStartB, iterMedianB, probeIntensities.end());

            if (*iterMedianA < 0.0 || *iterMedianB < 0.0) {
                Verbose::warn(0, "Negative intensity found in CEL file " + cel +
                                 ". Probe set '" + *it + "' ignored in SNPQC calculation.");
                continue;
            }
            vData.push_back(log((*iterMedianA)/(*iterMedianB))/log2);
        }
        delete chipStream;
        delete sparseMart;

        pExperiment->setRawSNPQC(HomHiLoCelListener::computeSNPQC(vData));
    }
}

void CNAnalysisMethodChipstream::newPlier(int iChipCount)
{
    m_iPlierChipCount = iChipCount;
    // Allocate memory for PLIER
    m_uiMaxProbeCount = m_pvProbes->getMaximumNumberProbesPerProbeSet();
    m_pdProbeEffects = new double[m_uiMaxProbeCount];
    m_pdChipEffects = new double[iChipCount];

    m_ppdPM = new double *[iChipCount];
    m_ppdMM = new double *[iChipCount];
    m_ppdResiduals = new double *[iChipCount];
    for(int i = 0; i < iChipCount; i++)
    {
        m_ppdPM[i] = new double[m_uiMaxProbeCount];
        m_ppdMM[i] = new double[m_uiMaxProbeCount];
        m_ppdResiduals[i] = new double[m_uiMaxProbeCount];
    }
    m_pPlier = new caffyplier;
    if (iChipCount == 1)
    {
        m_pPlier->setPlierFitFeatureResponse(false);
        m_pPlier->setPlierUseInputModel(true);
        m_pPlier->setSafetyZero(0);
    }
    else
    {
        m_pPlier->setPlierFitFeatureResponse(true);
        m_pPlier->setPlierUseInputModel(false);
        m_pPlier->setFixFeatureEffect(true);
    }
    m_pPlier->setPM(m_ppdPM);
    m_pPlier->setMM(m_ppdMM);
    m_pPlier->setResiduals(m_ppdResiduals);
    m_pPlier->setTargetResponse(m_pdChipEffects);
    m_pPlier->setFeatureResponse(m_pdProbeEffects);
    m_pPlier->setPlierOptOptimizationMethod(1);
    m_pPlier->setNumExp(iChipCount);
}

void CNAnalysisMethodChipstream::deletePlier()
{
    // Deallocate memory for plier
    if (m_pPlier != NULL)
    {
        delete [] m_pdProbeEffects;
        delete [] m_pdChipEffects;
        for (int i = 0; i < m_iPlierChipCount; i++)
        {
            delete [] m_ppdPM[i];
            delete [] m_ppdMM[i];
            delete [] m_ppdResiduals[i];
        }
        delete [] m_ppdResiduals;
        delete [] m_ppdPM;
        delete [] m_ppdMM;

        delete m_pPlier;
    }
}

void CNAnalysisMethodChipstream::newSketch()
{
    AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;

    // load the sketch from the reference
    file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
    if (m_pEngine->getOptBool("doDualNormalization") )
    {
        if(m_pEngine->getOptBool("cyto2"))
        {
            group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            tsv5 = group5->openTsv("SketchSNP", affx::FILE5_OPEN);
        }
        else
        {
            group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            tsv5 = group5->openTsv("SketchSNP", affx::FILE5_OPEN);
        }
    }
    else
    {
        group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        tsv5 = group5->openTsv("target-sketch", affx::FILE5_OPEN);
    }
    unsigned int uiSketchSize = 0;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        uiSketchSize++;
    }
    tsv5->close();
    delete tsv5;
    m_pvChipSketch = new std::vector<float>(uiSketchSize);
    m_pvSNPReferenceSketch = new std::vector<float>(uiSketchSize);
    m_pvCNReferenceSketch = new std::vector<float>(uiSketchSize);
    if (m_pEngine->getOptBool("doDualNormalization"))
    {
        tsv5 = group5->openTsv("SketchSNP", affx::FILE5_OPEN);
    }
    else
    {
        tsv5 = group5->openTsv("target-sketch", affx::FILE5_OPEN);
    }
    double d = 0;
    unsigned int uiSketchIndex = 0;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        tsv5->get(0, 0, &d); m_pvSNPReferenceSketch->at(uiSketchIndex) = (float)d;
        uiSketchIndex++;
    }
    tsv5->close();
    delete tsv5;
    if (m_pEngine->getOptBool("doDualNormalization"))
    {
        tsv5 = group5->openTsv("SketchCN", affx::FILE5_OPEN);
    }
    else
    {
        tsv5 = group5->openTsv("target-sketch", affx::FILE5_OPEN);
    }
    uiSketchIndex = 0;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        tsv5->get(0, 0, &d); m_pvCNReferenceSketch->at(uiSketchIndex) = (float)d;
        uiSketchIndex++;
    }
    tsv5->close();
    delete tsv5;
    if (group5 != NULL)
    {
        group5->close();
        delete group5;
    }
    file5.close();
}

void CNAnalysisMethodChipstream::deleteSketch()
{
    if (m_pvSNPReferenceSketch != NULL)
    {
        delete m_pvSNPReferenceSketch; m_pvSNPReferenceSketch = NULL;
        delete m_pvCNReferenceSketch; m_pvCNReferenceSketch = NULL;
        delete m_pvChipSketch; m_pvChipSketch = NULL;
    }
}

void CNAnalysisMethodChipstream::deleteProbes()
{
    if (m_pProbes != NULL)
    {
        if (m_pobjAdapterTypeNormTran != NULL) {delete m_pobjAdapterTypeNormTran; m_pobjAdapterTypeNormTran = NULL;}
        m_pvProbes->nullAll(); delete m_pvProbes; m_pvProbes = NULL;
        delete[] m_pProbes; m_pProbes = NULL;
    }
}

void CNAnalysisMethodChipstream::loadIntensities()
{
    m_pvProbes->quickSort(1); // by ProbeID

    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");

    // Load the raw intensities to use for the sketch
    if (!m_pEngine->isOptDefined("cels")) {throw(Except("Cannot determine CEL file for loading raw intensities."));}
    AffxString strCelFileName = m_pEngine->getOptVector("cels")[getExperiment()->getIndex()];
    affymetrix_fusion_io::FusionCELData cel;
    cel.SetFileName(strCelFileName.c_str());
    if (!cel.Read())
    {
        throw(Except("\nCan't read cel file: " + cel.GetFileName() + "\nMessage: " + StringUtils::ConvertWCSToMBS(cel.GetError())));
    }

    int numCells = cel.GetRows() * cel.GetCols();
    AffxMultiDimensionalArray<float> v(numCells);
    for (int k = 0; (k < numCells); k++)
    {
        v.set_unchecked(k, cel.GetIntensity(k));
    }
    getExperiment()->setMedianRawIntensity(v.median());
    v.initialize(1);

    int iSNPDataCount = 0;
    int iCNDataCount = 0;
    int iIndex = 0;
    for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
    {
        if ((iIndex < m_pvProbes->getCount()) && (m_pvProbes->getAt(iIndex)->getProbeID() == (iProbeIndex + 1)))
        {
            CNProbe* p = m_pvProbes->getAt(iIndex);
            CNProbeSet* pProbeSet = getProbeSets()->getAt(p->getProbeSetIndex());
            p->setIntensity(cel.GetIntensity(iProbeIndex));
            if (pProbeSet->isUseForSketch())
            {
                if (pProbeSet->processAsSNPNormalize())
                {
                    iSNPDataCount++;
                }
                if (pProbeSet->processAsCNNormalize())
                {
                    iCNDataCount++;
                }
            }
            iIndex++;
        }
    }
    if (!m_pEngine->getOptBool("doDualNormalization")) {iSNPDataCount = numCells; iCNDataCount = 0;}

    // Get ratio of antigenomic to genomic probes.
    {
        AffxString strFileName = m_pEngine->getOpt("antigenomic-probe-file");
        AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
        int iAntigenomicCount = 0;
        affx::File5_File file5;
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;
        affx::TsvFile tsv;
        tsv.m_optAutoTrim = true;
        
        if (strFileName != "")
        {
          std::string col;
          if (tsv.open(strFileName) != affx::TSV_OK) {Err::errAbort("Cannot open file: " + strFileName);}
          while (tsv.nextLevel(0) == affx::TSV_OK) {
            tsv.get(0,0,col);
            if (!col.empty()) {
              iAntigenomicCount++;
            }
          }
          tsv.clear();
        }
        else
        {
            if (m_pEngine->getOptBool("cyto2"))
            {
                iAntigenomicCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "AntigenomicProbes");
                file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
                group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            }
            else
            {
                iAntigenomicCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CN5", "AntigenomicProbes");
                file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
                group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
            }
        }
        if (iAntigenomicCount > 0)
        {
            AffxMultiDimensionalArray<int> v(iAntigenomicCount);
            AffxMultiDimensionalArray<float> vIntensities(iAntigenomicCount);
            
            int i = 0;
            int iIndex = 0;
            if (strFileName != "")
            {
              AffxString probe;
              if (tsv.open(strFileName) != affx::TSV_OK ) {Err::errAbort("Cannot open file: " + strFileName);}
              while (tsv.nextLevel(0) == affx::TSV_OK)  {
                tsv.get(0,0,probe);
                if (!probe.empty()) {
                  i = AffxByteArray(probe).parseInt();
                  v.set(iIndex, i - 1);
                  iIndex++;
                }
              }
              tsv.clear();
            }
            else
            {
                tsv5 = group5->openTsv("AntigenomicProbes", affx::FILE5_OPEN);
                while (tsv5->nextLine() == affx::FILE5_OK)
                {
                    tsv5->get(0, 0, &i);
                    v.set(iIndex, (i - 1));
                    iIndex++;
                }
                tsv5->close();
                delete tsv5;
            }
            v.quickSort();
            iIndex = 0;
            for (int k = 0; (k < numCells); k++)
            {
                if (v.binarySearch(k) != -1)
                {
                    vIntensities.set(iIndex, cel.GetIntensity(k));
                    iIndex++;
                }
            }
            double dMedianAntigenomicProbes = vIntensities.median();
            unsigned int uiCount = 0;
            for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
            {
                CNProbe* p = m_pvProbes->getAt(iIndex);
                CNProbeSet* pProbeSet = getProbeSets()->getAt(p->getProbeSetIndex());
                if (isCytoScanHD && ((pProbeSet->getProcessFlag() != 1) || (!pProbeSet->processAsCNNormalize()))) {continue;}
                if (pProbeSet->processAsCNNormalize())
                {
                    uiCount++;
                }
            }
            vIntensities.initialize(uiCount);
            i = 0;
            for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
            {
                CNProbe* p = m_pvProbes->getAt(iIndex);
                CNProbeSet* pProbeSet = getProbeSets()->getAt(p->getProbeSetIndex());
                if (isCytoScanHD && ((pProbeSet->getProcessFlag() != 1) || (!pProbeSet->processAsCNNormalize()))) {continue;}
                if (pProbeSet->processAsCNNormalize())
                {
                    vIntensities.set(i, p->getIntensity());
                    i++;
                }
            }
            double dMedianCopynumberProbes = vIntensities.median();
            if (dMedianCopynumberProbes != 0)
            {
                getExperiment()->setAntigenomicRatio((float)(dMedianAntigenomicProbes / dMedianCopynumberProbes));
            }
        }
        if (strFileName == "")
        {
            group5->close();
            delete group5;
            file5.close();
        }
    }

    adapterTypeNormalization(cel);
    adjustIntensitiesUsingCovariateSignalAdjustment();

    if (m_pEngine->getOptBool("keep-intermediate-data"))
    {
        if (InstanceOf(this, CNAnalysisMethodReference) == true)
        {
            writeIntensities("reference-before-normalized");
        }
        else
        {
            writeIntensities("single-sample-before-normalized");
        }
    }

    if (m_pvSNPData != NULL)
    {
        delete m_pvSNPData;
    }
    m_pvSNPData = new std::vector<float>(iSNPDataCount);

    if (m_pvCNData != NULL)
    {
        delete m_pvCNData;
    }
    m_pvCNData = new std::vector<float>(iCNDataCount);

    int iSNPIndex = 0;
    int iCNIndex = 0;
    if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
    {
    // Working on 2.7M or 310K chips
        for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
        {
            CNProbe* p = m_pvProbes->getAt(iIndex);
            CNProbeSet* pProbeSet = getProbeSets()->getAt(p->getProbeSetIndex());
            if (pProbeSet->isUseForSketch())
            {
                if (pProbeSet->processAsSNPNormalize())
                {
                    m_pvSNPData->at(iSNPIndex) = p->getIntensity();
                    iSNPIndex++;
                }
                if (pProbeSet->processAsCNNormalize())
                {
                    m_pvCNData->at(iCNIndex) = p->getIntensity();
                    iCNIndex++;
                }
            }
        }
    }
    else
    {
        //MG Working on snp6 or snp7 chips i.e. cyto2 = false 
        if (m_pEngine->getOptBool("doDualNormalization"))
        {
            // Working on snp6/7 chips with dual normalization.
            if ((m_pEngine->isOptDefined("adapter-type-normalization")) && (m_pEngine->getOptBool("adapter-type-normalization")))
            {
                // Working on snp6/7 chips with dual and adapter type normalization.
                for (int iProbeIndex = 0; iProbeIndex < m_pvProbes->getCount(); iProbeIndex++)
                { 
                    CNProbe* p = m_pvProbes->getAt(iProbeIndex);
                    CNProbeSet* pProbeSet = getProbeSets()->getAt(p->getProbeSetIndex());
                    if (pProbeSet->isUseForSketch())
                    {
                        {
                            if (pProbeSet->processAsSNPNormalize())
                            {
                                m_pvSNPData->at(iSNPIndex) = m_pobjAdapterTypeNormTran->transform(p->getProbeID()-1, 0, cel.GetIntensity(p->getProbeID()-1)); 
                                iSNPIndex++;
                            }
                            if (pProbeSet->processAsCNNormalize())
                            {
                                m_pvCNData->at(iCNIndex) = m_pobjAdapterTypeNormTran->transform(p->getProbeID()-1, 0, cel.GetIntensity(p->getProbeID()-1)); 
                                iCNIndex++;
                            }
                        }
                   }
                }
            }
            else
            {
                // Working on snp6/7 chips dual but no adapter type normalization.
                for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
                {
                    CNProbe* p = m_pvProbes->getAt(iIndex);
                    CNProbeSet* pProbeSet = getProbeSets()->getAt(p->getProbeSetIndex());
                    if (pProbeSet->isUseForSketch())
                    {
                        if (pProbeSet->processAsSNPNormalize())
                        {
                            m_pvSNPData->at(iSNPIndex) = p->getIntensity();
                            iSNPIndex++;
                        }
                        if (pProbeSet->processAsCNNormalize())
                        {
                            m_pvCNData->at(iCNIndex) = p->getIntensity();
                            iCNIndex++;
                        }
                   }
                }

            }

        }
        else
        {
            // Working on snp6/7 chips with single normalization.
            if ((m_pEngine->isOptDefined("adapter-type-normalization")) && (m_pEngine->getOptBool("adapter-type-normalization")))
            {
                // Working on snp6/7 chips with single and adapter type normalization. 
                for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
                {
                    // Note that when doing single normalization we load all intensities into the m_pvSNPData vector.
                    m_pvSNPData->at(iSNPIndex) = m_pobjAdapterTypeNormTran->transform(iProbeIndex, 0, cel.GetIntensity(iProbeIndex));
                    iSNPIndex++;
                }


/*  This code seems more logical, but for backwards compatibility with workflow/cyto on snp6 chips we do quantile (i.e.single)
    normalization   on all probesets on the chip, not just the union of SNP and CN probesets. 

                for (int iProbeIndex = 0; iProbeIndex < m_pvProbes->getCount(); iProbeIndex++)
                {
                    CNProbe* p = m_pvProbes->getAt(iProbeIndex);
                    CNProbeSet* pProbeSet = getProbeSets()->getAt(p->getProbeSetIndex());
                    if (pProbeSet->isUseForSketch())
                    {
                        m_pvSNPData->at(iSNPIndex) = m_pobjAdapterTypeNormTran->transform(p->getProbeID()-1, 0, cel.GetIntensity(p->getProbeID()-1));
                        iSNPIndex++;
                    }
                } 
*/
            }
            else
            {
                // Working on snp6/7 chips with single and no adapter type normalization. 
                for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
                {
                    m_pvSNPData->at(iSNPIndex) = cel.GetIntensity(iProbeIndex);
                    iSNPIndex++;
                }
            }

        }
    }

    cel.Close();
}

void CNAnalysisMethodChipstream::adapterTypeNormalization(affymetrix_fusion_io::FusionCELData& cel)
{

  // adapter type normalization
    if ((m_pEngine->isOptDefined("adapter-type-normalization")) && (m_pEngine->getOptBool("adapter-type-normalization")))
    {
        int numCells = cel.GetRows() * cel.GetCols();
        std::vector<float> data(numCells);
        for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
        {
            data[iProbeIndex] = cel.GetIntensity(iProbeIndex);
        }
        //m_pobjAdapterTypeNormTran->newChip(data);
        m_pobjAdapterTypeNormTran->setCelFileIndex(0);
        m_pobjAdapterTypeNormTran->initializeData(data);
        for (int iProbeIndex = 0; (iProbeIndex < m_pvProbes->getCount()); iProbeIndex++)
        {
            CNProbe* p = m_pvProbes->getAt(iProbeIndex);
            p->setIntensity(m_pobjAdapterTypeNormTran->transform((p->getProbeID() - 1), 0, p->getIntensity()));
        }
    }
}

void CNAnalysisMethodChipstream::initializeSnpParameters(snp_param& sp, int iCallMethod)
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

void CNAnalysisMethodChipstream::setBrlmmpParameters(snp_param& sp)
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
    //

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

void CNAnalysisMethodChipstream::loadSnpDistributions()
{
    int m_iGenderCode = getExperiment()->getRawIntensityRatioGenderAsInt(); 
    if (m_iGenderCode != m_iPrevGenderCode)
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
        m_iPrevGenderCode = m_iGenderCode;
    }
}

bool CNAnalysisMethodChipstream::loadGenotypeCallsFromFile()
{
    bool bSuccessful = false;
    AffxString strFileName = m_pEngine->getOpt("genotype-call-override-file");
    if (strFileName != "")
    {
        for (int iProbeSetIndex = 0; (iProbeSetIndex < getProbeSets()->getCount()); iProbeSetIndex++)
        {
            CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
            if (p->processAsSNP())
            {
                p->setGenotypeCall((char)-1);
                p->setGenotypeConfidence(0);
            }
        }
        AffxString strSampleId = getSampleId(getExperiment()->getExperimentName());
        if (strSampleId != "")
        {
          // DO NOT USE TSV.
          // This section was changed from TSV to RowFile for performance reasons.
          // Further investigation needs to be done as to the root cause for
          // the TSV penalty.
          /*
            affx::TsvFile tsv;
            std::string col, colName;
            AffxByteArray baColumn;
            AffxString strProbeSetName;
            CNProbeSet objSearch;
            tsv.m_optAutoTrim = true;
            if (tsv.open(strFileName) == affx::TSV_OK ) {
              for (int iColumnNumber = 1; (iColumnNumber < tsv.getColumnCount(0)); iColumnNumber++) {
                tsv.cidx2cname(0, iColumnNumber, colName);
                if (strSampleId == getSampleId(colName)) {
                  bSuccessful = true;
                  Verbose::progressBegin(1, "Loading override genotype calls for sample\t" + strSampleId + " ", 1, 1, 1);
                  Verbose::progressStep(1);
                  while (tsv.nextLevel(0) == affx::TSV_OK)  {
                    tsv.get(0,0,strProbeSetName);
                    objSearch.setProbeSetName(strProbeSetName);
                    int iSearchIndex = getProbeSets()->binarySearch(objSearch, 0);
                    if (iSearchIndex != -1) {
                      CNProbeSet* p = getProbeSets()->getAt(iSearchIndex);
                      tsv.get(0,iColumnNumber, col);
                      p->setGenotypeCall((char)AffxByteArray(col).parseInt());
                      p->setGenotypeConfidence(0);
                    }
                  }
                  Verbose::progressEnd(1, "Done");
                  break;
                }
              }
              tsv.clear();
            }
          */
            RowFile rowFile;
            std::vector< std::string > row;
            AffxString col;
            AffxString strProbeSetName;
            CNProbeSet objSearch;
            bool okOpen = false;
            try {
              rowFile.open(strFileName);
              okOpen = true;
            }
            catch (...) {
              // Just return false below, callee will abort. 
            }
            if ( okOpen ) {
              okOpen = rowFile.nextRow(row);
              for ( int iColumnNumber = 1; !bSuccessful && okOpen &&
                      (iColumnNumber < row.size()); iColumnNumber++)  {
                col = row[iColumnNumber];
                if (strSampleId == getSampleId(col.trim()))  {
                  bSuccessful = true;
                  Verbose::progressBegin(1, "Loading override genotype calls for sample\t" + strSampleId + " ", 1, 1, 1);
                  Verbose::progressStep(1);
                  const std::string *line;
                  while ( (line = rowFile.nextLine() ) ) {
                    col = line->substr(0,line->find('\t'));
                    strProbeSetName = col.trim();
                    objSearch.setProbeSetName(strProbeSetName);
                    int iSearchIndex = getProbeSets()->binarySearch(objSearch, 0);
                    if (iSearchIndex != -1)
                      {
                        CNProbeSet* p = getProbeSets()->getAt(iSearchIndex);
                        col = *line;
                        for ( int index = 0; index < iColumnNumber; index++ ) {
                          col = col.substr(col.find('\t') + 1);
                        }
                        col = col.substr(0,col.find('\t'));
                        p->setGenotypeCall((char)AffxByteArray(col.trim()).parseInt());
                        p->setGenotypeConfidence(0);
                      }
                  }
                  Verbose::progressEnd(1, "Done");
                }
              }
              rowFile.close();
            }
        }
    }
    return bSuccessful;
}

AffxString CNAnalysisMethodChipstream::getSampleId(const AffxString& strCelFileName)
{
    // Use the entire file name converted to lower case, minus ".cel" suffix
    AffxString strSampleId = Fs::basename(strCelFileName);
    transform(    strSampleId.begin(),
                strSampleId.end(),
                strSampleId.begin(),
                static_cast<int(*)(int)>(tolower));
    AffxString::size_type iFindIndex = strSampleId.rfind(".cel");
    if (iFindIndex != AffxString::npos) {
        strSampleId = strSampleId.substring(0, iFindIndex);
    }
    return strSampleId;
}

void CNAnalysisMethodChipstream::calculateXY()
{
    AffxArray<AffxString> vYProbeSetsToProcess;

    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
    
    if (tsv.open(m_pEngine->getOpt("chrY-probes")) == affx::TSV_OK)
    {
        while (tsv.nextLevel(0) == affx::TSV_OK)
        {
          AffxString strProbeSetName;
          tsv.get(0,1, strProbeSetName);
          if (strProbeSetName != "")
            {
              vYProbeSetsToProcess.add(new AffxString(strProbeSetName));
            }
        }
        tsv.clear();
    } else {throw(Except("Cannot open file: " + m_pEngine->getOpt("chrY-probes")));}
    vYProbeSetsToProcess.quickSort(0);

    CNExperiment* pobjExperiment = getExperiment();
    pobjExperiment->setXX(false);
    pobjExperiment->setY(false);
    // Allocate the data vectors.
    int iChrRefCount = 0;
    int iChrXCount = 0;
    int iChrYCount = 0;
    for (int uiRowIndex = 0; (uiRowIndex < getProbeSets()->getCount()); uiRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(uiRowIndex);
        AffxString strProbeSetName = pobjProbeSet->getProbeSetName();
        if ((pobjProbeSet->getChromosome() == m_iYChromosome) && (vYProbeSetsToProcess.getCount() > 0) && (vYProbeSetsToProcess.binarySearch(strProbeSetName, 0) == -1)) {continue;}
        if (pobjProbeSet->getChromosome() == m_pEngine->getOptInt("reference-chromosome")) {iChrRefCount++;}
        else if (pobjProbeSet->getChromosome() == m_iXChromosome) {iChrXCount++;}
        else if (pobjProbeSet->getChromosome() == m_iYChromosome) {iChrYCount++;}
    }
    AffxMultiDimensionalArray<float> vChrRefMedian(iChrRefCount);
    AffxMultiDimensionalArray<float> vChrXSignal(iChrXCount);
    AffxMultiDimensionalArray<float> vChrYSignal(iChrYCount);

    int iChrRefIndex = 0;
    int iChrXIndex = 0;
    int iChrYIndex = 0;
    // For each probe set load the relevant data into the vectors so we can calculate the median.
    for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
        AffxString strProbeSetName = pobjProbeSet->getProbeSetName();
        if ((pobjProbeSet->getChromosome() == m_iYChromosome) && (vYProbeSetsToProcess.getCount() > 0) && (vYProbeSetsToProcess.binarySearch(strProbeSetName, 0) == -1)) {continue;}
        if (pobjProbeSet->getChromosome() == m_pEngine->getOptInt("reference-chromosome"))
        {
            vChrRefMedian.set(iChrRefIndex, pobjProbeSet->getMedianSignal());
            iChrRefIndex++;
        }
        else if (pobjProbeSet->getChromosome() == m_iXChromosome)
        {
            vChrXSignal.set(iChrXIndex, pobjProbeSet->getSignalEstimate());
            iChrXIndex++;
        }
        else if (pobjProbeSet->getChromosome() == m_iYChromosome)
        {
            vChrYSignal.set(iChrYIndex, pobjProbeSet->getSignalEstimate());
            iChrYIndex++;
        }
    }
    vYProbeSetsToProcess.deleteAll();
    double dChrRefMedian = vChrRefMedian.median();

    // Set the experiment hasXX value.
    m_dXXRatio = (vChrXSignal.median() / dChrRefMedian);
    pobjExperiment->setXXRatio((float)m_dXXRatio);
    pobjExperiment->setXX((m_dXXRatio > m_pEngine->getOptDouble("xx-cutoff")) && (m_dXXRatio < m_pEngine->getOptDouble("xx-cutoff-high")));

    // Set the experiment hasY value.
    m_dYRatio = (vChrYSignal.median() / dChrRefMedian);
    pobjExperiment->setYRatio((float)m_dYRatio);
    pobjExperiment->setY(m_dYRatio > m_pEngine->getOptDouble("y-cutoff"));
    
    Verbose::out(1, "CNAnalysisMethodChipstream::calculateXY()\t" + getExperiment()->getExperimentName() + "\tXXRatio\t" + ::getDouble(m_dXXRatio, 6) + "\tYRatio\t" + ::getDouble(m_dYRatio, 6) + "\t" + (getExperiment()->hasXX() ? "XX" : "X") + (getExperiment()->hasY() ? "Y" : ""));

    if ((!m_pEngine->getOptBool("cyto2")) && (!m_pEngine->getOptBool("cytoscan-hd")))
    {
        getExperiment()->setXX(getExperiment()->getRawIntensityRatioGenderAsInt() == affx::Female);
        getExperiment()->setY(getExperiment()->getRawIntensityRatioGenderAsInt() == affx::Male);
        Verbose::out(1, "CNAnalysisMethodChipstream::calculateXY(), (cyto2 = false) override\t" + getExperiment()->getExperimentName() + "\tXXRatio\t" + ::getDouble(m_dXXRatio, 6) + "\tYRatio\t" + ::getDouble(m_dYRatio, 6) + "\t" + (getExperiment()->hasXX() ? "XX" : "X") + (getExperiment()->hasY() ? "Y" : ""));
    }
}

void CNAnalysisMethodChipstream::adjustIntensitiesPDNN()
{

        AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
        affx::File5_File file5;
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;

        int iCount = 0;
        if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
        {
            iCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "ProbeEffects");
        }
        else
        {
            iCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CN5", "CN5.plier.feature-response");
        }
        if (iCount <= 0) {throw(Except("Cannot find any probes in the CN reference file to process."));}

        file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
        if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
        {
            group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
        }
        else
        {
            group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        }

        if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
        {
            tsv5 = group5->openTsv("ProbeEffects", affx::FILE5_OPEN);
            if (tsv5->getColumnCount(0) < 6) {throw(Except("The CN Reference file does not contain Predicted Intensities.  PDNN cannot be computed."));}
        }
        else
        {
            tsv5 = group5->openTsv("CN5.plier.feature-response", affx::FILE5_OPEN);
            if (tsv5->getColumnCount(0) < 6) {throw(Except("The CN Reference file does not contain Predicted Intensities.  PDNN cannot be computed."));}
        }

        unsigned int uiIndex = 0;
        float fTempValue = 0.0;
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            tsv5->get(0, 5, &fTempValue); 
            m_pProbes[uiIndex].setPredictedIntensity(fTempValue);
            uiIndex++;
        }
        tsv5->close();
        delete tsv5;

        group5->close();
        delete group5;
        file5.close();


        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(m_pEngine->getOpt("local-gc-background-intensity-adjustment-method"));
        am->setEngine(m_pEngine);
        am->setup(*getExperiment(), *getProbeSets());
        am->setProbes(*m_pvProbes);
        am->run();
        delete am;
}

void CNAnalysisMethodChipstream::adjustIntensitiesHighPassFilter()
{
    if (!m_pEngine->isOptDefined("image-correction-intensity-adjustment-method")) {return;}
    if (m_pEngine->getOpt("image-correction-intensity-adjustment-method") == "none") {return;}
    if (m_pEngine->getOpt("image-correction-intensity-adjustment-method") != "")
    {
        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(m_pEngine->getOpt("image-correction-intensity-adjustment-method"));
        am->setEngine(m_pEngine);
        am->setup(*getExperiment(), *getProbeSets());
        am->setProbes(*m_pvProbes);
        am->run();
        delete am;
    }
}

void CNAnalysisMethodChipstream::adjustIntensitiesUsingCovariateSignalAdjustment()
{
    if (!m_pEngine->isOptDefined("signal-adjustment-covariates")) {return;}
    if (m_pEngine->getOpt("signal-adjustment-covariates") == "none") {return;}
    if (m_pEngine->getOpt("signal-adjustment-covariates") != "")
    {
        CNAnalysisMethodFactory amFactory;
        CNAnalysisMethod* am = NULL;
        am = amFactory.CNAnalysisMethodForString(m_pEngine->getOpt("signal-adjustment-covariates"));
        am->setEngine(m_pEngine);
        am->setup(*getExperiment(), *getProbeSets());
        am->setProbes(*m_pvProbes);
        am->run();
        delete am;
    }
}

void CNAnalysisMethodChipstream::adjustAlternateIntensitiesUsingCovariateSignalAdjustment(affymetrix_fusion_io::FusionCELData& cel, 
                                                                                          std::vector<float>& adjustedIntensities)
{
    if (!m_pEngine->isOptDefined("signal-adjustment-covariates")) {return;}
    if (m_pEngine->getOpt("signal-adjustment-covariates") == "none") {return;}
    if (m_pEngine->getOpt("signal-adjustment-covariates") != "")
    {
        int numCells = cel.GetRows() * cel.GetCols();
        adjustedIntensities.resize(numCells);
        for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
        {
            adjustedIntensities[iProbeIndex] = cel.GetIntensity(iProbeIndex);
        }

        CNAnalysisMethodCovariateSignalAdjuster* pMethod = new CNAnalysisMethodCovariateSignalAdjuster();

        pMethod->setEngine(m_pEngine);
        pMethod->setup(*getExperiment(), *getProbeSets());
        pMethod->setProbes(*m_pvProbes);
        pMethod->setAlternateIntensities(adjustedIntensities);
        pMethod->run();
        delete pMethod;
    }
}
