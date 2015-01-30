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
#include "copynumber/CNReferenceMethodWaveCorrection.h"
//
#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNCychp.h"
#include "copynumber/CNLog2RatioData.h"
#include "copynumber/CNAnalysisMethodCovariateParams.h"
#include "copynumber/CNAnalysisMethodCovariateLRAdjuster.h"
//
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "portability/affy-base-types.h"
#include "util/AffxMultiDimensionalArray.h"
#include "util/Fs.h"
#include "util/TmpFileFactory.h"
#include "util/Verbose.h"
//
#include <algorithm>
#include <limits>
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNReferenceMethodWaveCorrection::explainSelf()
    {
    CNReferenceMethodWaveCorrection obj;
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
std::vector<SelfDoc::Opt> CNReferenceMethodWaveCorrection::getDefaultDocOptions()
    {
    std::vector<SelfDoc::Opt> opts;

    // SelfDoc::Opt(name, type, value, default, min, max, description)

    SelfDoc::Opt dTrim = {"trim", SelfDoc::Opt::Double, "2.0", "2.0", "NA", "NA", "Log2Ratio Trim value."};
    opts.push_back(dTrim);

    SelfDoc::Opt dPercentile = {"percentile", SelfDoc::Opt::Double, "0.75", "0.75", "NA", "NA", "High Percentile value."};
    opts.push_back(dPercentile);

    SelfDoc::Opt iWaveCount = {"wave-count", SelfDoc::Opt::Integer, "6", "6", "0", "NA", "Number of waves to add to the reference."};
    opts.push_back(iWaveCount);

    SelfDoc::Opt bDemean = {"demean", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", "Demean the input to the SVD."};
    opts.push_back(bDemean);

//    SelfDoc::Opt dCutoff = {"wave-cutoff", SelfDoc::Opt::Double, "0.000001", "0.000001", "0", "NA", "The cutoff for wave calculation."};
//    opts.push_back(dCutoff);

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
SelfCreate * CNReferenceMethodWaveCorrection::newObject(
    std::map<std::string,std::string>& params)
    {
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNReferenceMethodWaveCorrection* pMethod = new CNReferenceMethodWaveCorrection();
    std::string strPrefix = getPrefix();

    pMethod->m_dTrim = setupDoubleParameter("trim", strPrefix, params, doc, opts);
    pMethod->m_dPercentile = setupDoubleParameter("percentile", strPrefix, params, doc, opts);
    pMethod->m_iWaveCount = setupIntParameter("wave-count", strPrefix, params, doc, opts);
    pMethod->m_bDemean = setupBoolParameter("demean", strPrefix, params, doc, opts);
//    pMethod->m_dWaveCutoff = setupDoubleParameter("wave-cutoff", strPrefix, params, doc, opts);

    return pMethod;
    }

/**
 * @brief Constructor
 */
CNReferenceMethodWaveCorrection::CNReferenceMethodWaveCorrection()
{
    m_iLastAutosomeChromosome = 0;
    m_iAutosomeCount = 0;
    m_dTrim = 0;
    m_dPercentile = 0;
    m_iWaveCount = 0;
    m_bDemean = false;
    m_dWaveCutoff = 0;
    m_iRefWaveCount = 0;
    m_dCNQCCutoff = 0;
    m_dSNPQCCutoff = 0;
    m_uiPassedQCCount = 0;
    m_bKeepTempData = false;
    m_bForce = false;
    m_iWavinessSegCountCutoff = 0;
    m_bUseHighWavinessSegCount = true;
    m_bCyto2 = false;
}

/**
 * @brief Destructor
 */
CNReferenceMethodWaveCorrection::~CNReferenceMethodWaveCorrection()
{
}

float CNWaveCorrectionProbeSet::getCovariateValue(int iIndex)
{
    return AllCovariates[iIndex];
}
unsigned char CNWaveCorrectionProbeSet::getProcessFlag()
{
    return ProcessFlag;
}
unsigned char CNWaveCorrectionProbeSet::getChromosome()
{
    return Chromosome;
}
double CNWaveCorrectionProbeSet::getLog2Ratio()
{
    return Log2Ratio;
}
void CNWaveCorrectionProbeSet::setLog2Ratio(double f)
{
    Log2Ratio = f;
}

/**
 * @brief Run the calculations.
 */
void CNReferenceMethodWaveCorrection::run()
{
    Verbose::out(1, "CNReferenceMethodWaveCorrection::run(...) start");
    bool bCNReferenceRun = true;
    if (m_pEngine->isOptDefined("cn-reference-input"))
    {
        bCNReferenceRun = false;
        m_strReferenceFileName = m_pEngine->getOpt("cn-reference-input");
        // Create temp file
        TmpFileFactory tempFactory;
        tempFactory.setTmpdir(m_pEngine->getOpt("temp-dir"));
        m_strTempFileName = tempFactory.genFilename_basic( "CNReference_", ".a5.tmp");
        m_pEngine->setOpt("temp-reference-file", m_strTempFileName);
    }
    else
    {
        m_strReferenceFileName = m_pEngine->getOpt("reference-output");
        m_strTempFileName = m_pEngine->getOpt("temp-reference-file");
    }
    m_bCyto2 = true;
    m_strGroup = "Cyto2";
    m_strTsv = "MedianSignals";
    int iProbeSetCount = affx::File5_Tsv::getFileTsvLineCount(m_strReferenceFileName, m_strGroup, m_strTsv);
    if (iProbeSetCount <= 0)
    {
        m_bCyto2 = false;
        m_strGroup = "CopyNumber";
        m_strTsv = "Reference";
        iProbeSetCount = affx::File5_Tsv::getFileTsvLineCount(m_strReferenceFileName, m_strGroup, m_strTsv);
    }
    if (iProbeSetCount <= 0) {throw(Except("Cannot find any probe sets to process in the reference."));}
    CNWaveCorrectionProbeSet* pProbeSets = new CNWaveCorrectionProbeSet[iProbeSetCount];
    AffxArray<CNWaveCorrectionProbeSet> vProbeSets;
    vProbeSets.reserve(iProbeSetCount);
    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;
    affx::File5_Tsv* tsv5Waves = NULL;
    affx::File5_Tsv* tsv5Covariates = NULL;
    file5.open(m_strReferenceFileName, affx::FILE5_OPEN_RO);
    group5 = file5.openGroup(m_strGroup, affx::FILE5_OPEN_RO);
    AffxString str;
    int ii = 0;
    double d = 0;
    m_iRefWaveCount = 0;
    int iNumCovariates = 0;
    bool useCovariateTsv = m_pEngine->getEngineName() == "CNCytoEngine"  &&
                           group5->name_exists("Covariates");

    tsv5 = group5->openTsv(m_strTsv, affx::FILE5_OPEN_RO);
    if (bCNReferenceRun)
    {
        if (useCovariateTsv)
        {
            tsv5Covariates = group5->openTsv("Covariates", affx::FILE5_OPEN_RO);
            iNumCovariates = tsv5Covariates->getColumnCount(0);
        }
    }
    else {
        tsv5Waves = group5->openTsv("WaveCorrection", affx::FILE5_OPEN_RO);
        m_iRefWaveCount = tsv5Waves->getColumnCount(0);
    }
    int iColumnCount = tsv5->getColumnCount(0);
    unsigned int uiIndex = 0;
    int iOrder = 0;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        if (m_iRefWaveCount > 0) {tsv5Waves->nextLine();}
        if (iNumCovariates > 0) {tsv5Covariates->nextLine();}
        CNWaveCorrectionProbeSet* p = &pProbeSets[uiIndex];
        uiIndex++;
        p->Order = iOrder; iOrder++;
        tsv5->get(0, 0, &str); p->ProbeSetName = str;
        if (iColumnCount > 7)
        {
            tsv5->get(0, 7, &ii); p->Chromosome = (unsigned char)ii;
            tsv5->get(0, 8, &ii); p->Position = (unsigned int)ii;
        }
        tsv5->get(0, 1, &d); p->MedianSignal = (float)d;
        tsv5->get(0, 2, &d); p->XXMedianSignal = (float)d;
        tsv5->get(0, 3, &d); p->YMedianSignal = (float)d;
        p->Waves.resize(m_iRefWaveCount + m_iWaveCount);
        for (int iWaveIndex = 0; (iWaveIndex < m_iRefWaveCount); iWaveIndex++)
        {
            tsv5Waves->get(0, iWaveIndex, &d);
            p->Waves[iWaveIndex] = d;
        }
        if (iColumnCount > 9)
        {
            tsv5->get(0, 9, &ii);
            if ((ii != 1) && (ii != 3)) {
                continue;
            }
        }
        if (iNumCovariates > 0)
        {
            p->ProcessFlag = (unsigned char)ii;

            // Assign iNumCovariates elements _exactly_
            p->AllCovariates.resize(iNumCovariates);
            std::vector<float>(p->AllCovariates).swap(p->AllCovariates);

            for (int iCovariateIndex = 0; iCovariateIndex < iNumCovariates; iCovariateIndex++)
            {
                tsv5Covariates->get(0, iCovariateIndex, &d);
                p->AllCovariates[iCovariateIndex] = d;
            }
        }
        vProbeSets.add(p);
    }
    tsv5->close();
    delete tsv5;
    if (m_iRefWaveCount > 0)
    {
        tsv5Waves->close();
        delete tsv5Waves;
    }
    if (iNumCovariates > 0)
    {
        tsv5Covariates->close();
        delete tsv5Covariates;
    }

    group5->close();
    delete group5;
    file5.close();

    vProbeSets.quickSort(1);

    unsigned int uiCelCount = 0;
    if (m_pEngine->isOptDefined("cn-reference-input"))
    {
        uiCelCount = m_pEngine->getOptVector("cychps").size();
    }
    else
    {
        uiCelCount = m_pEngine->getOptVector("cels").size();
    }
    
    if (bCNReferenceRun == false && m_strSelectedQC != "") {Verbose::out(1, "Using selected-qc=" + m_strSelectedQC);}

    m_vAutosomeCV.resize(uiCelCount);
    calculateChromosomeStatistics(vProbeSets);
    for (int iWaveIndex = 0; (iWaveIndex < m_iWaveCount); iWaveIndex++)
    {
        loadPCALikeVector(m_strTempFileName, vProbeSets, iWaveIndex);
    }

    updateReference(vProbeSets, iProbeSetCount);
    if ((!bCNReferenceRun) && (!m_bKeepTempData))
    {
        Fs::rm(m_pEngine->getOpt("temp-reference-file"), false);
    }

    delete[] pProbeSets;
    Verbose::out(1, "CNReferenceMethodWaveCorrection::run(...) end");
}

void writeWaveCorrectionLog2Ratios(const string& dir, const string& name, AffxArray<CNWaveCorrectionProbeSet>& vProbeSets)
{
  affx::File5_File file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5 = NULL;

  std::string analysisString = "analysis";
  std::string fileName = Fs::join(dir, analysisString, name);

  Verbose::out(3, "Writing Log2Ratios for " + fileName);

  try{

    file5.open(fileName, affx::FILE5_OPEN_CREATE);
    group5 = file5.openGroup("Signals", affx::FILE5_REPLACE);
    tsv5 = group5->openTsv("Signals", affx::FILE5_REPLACE);
    tsv5->defineColumn(0,0,"ProbeSetName", affx::FILE5_DTYPE_STRING,40);
    tsv5->defineColumn(0,1,"Chromosome", affx::FILE5_DTYPE_CHAR);
    tsv5->defineColumn(0,2,"Position", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0,3,"Log2Ratio", affx::FILE5_DTYPE_FLOAT);

    CNWaveCorrectionProbeSet* pProbeSet;
    for(int iProbeSetIndex=0;  iProbeSetIndex < vProbeSets.getCount(); iProbeSetIndex++)
    {
        pProbeSet = vProbeSets.getAt(iProbeSetIndex);
        {
          tsv5->set_string(0,0, pProbeSet->ProbeSetName);
          tsv5->set_c(0,1, pProbeSet->getChromosome());
          tsv5->set_i(0,2, pProbeSet->Position);
          tsv5->set_f(0,3, pProbeSet->getLog2Ratio());
          tsv5->writeLevel(0);
        }
    }

    tsv5->close();
    delete tsv5;
    delete group5;
    file5.close();
  }
  catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump intensities to file."));}
}

void CNReferenceMethodWaveCorrection::loadPCALikeVector(    const std::string& strTempFileName,
                                                            AffxArray<CNWaveCorrectionProbeSet>& vProbeSets,
                                                            int iWaveIndex)
{
    Verbose::out(1, "CNReferenceMethodWaveCorrection::loadPCALikeVector(...) for wave " + ToStr(iWaveIndex));
    int iProbeSetCount = vProbeSets.size();
    unsigned int uiCelCount = 0;
    bool bCNReferenceRun = true;
    if (m_pEngine->isOptDefined("cn-reference-input"))
    {
        bCNReferenceRun = false;
        uiCelCount = m_pEngine->getOptVector("cychps").size();
    }
    else
    {
        uiCelCount = m_pEngine->getOptVector("cels").size();
    }
    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;
    if (!bCNReferenceRun)
    {
        if (iWaveIndex == 0)
        {
            file5.open(strTempFileName, affx::FILE5_CREATE);
            group5 = file5.openGroup("Signals", affx::FILE5_CREATE);
        }
        else
        {
            file5.open(strTempFileName, affx::FILE5_OPEN);
            group5 = file5.openGroup("Signals", affx::FILE5_OPEN);
        }
    }
    else
    {
        file5.open(strTempFileName, affx::FILE5_OPEN);
        group5 = file5.openGroup("Intensities", affx::FILE5_OPEN);
    }

    std::vector<int> vAutosomeProbeSetIndexes;
    int iAutosomeProcessCount = Min(m_iAutosomeCount,(int)2e4);
    int iAutosomeStep = (int)::roundDouble((double)m_iAutosomeCount / (double)iAutosomeProcessCount);

    for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex += iAutosomeStep)
    {
        CNWaveCorrectionProbeSet* p = vProbeSets[iProbeSetIndex];
        if (p->Chromosome <= m_iLastAutosomeChromosome)
        {
            vAutosomeProbeSetIndexes.push_back(iProbeSetIndex);
        } else {break;}
    }
    iAutosomeProcessCount = (int)vAutosomeProbeSetIndexes.size();
    if (iWaveIndex == 0)
    {
        m_vHasXX.resize(uiCelCount);
        m_vHasY.resize(uiCelCount);
        m_vYTargets.resize(uiCelCount);
        m_vPassedQC.resize(uiCelCount);
    }
    CNCychp cychp;
    CNCychpProbeSetsCopyNumber* pv = NULL;
    AffxArray<CNCychpProbeSetsCopyNumber> v;
    AffxArray<CNWaveCorrectionProbeSet> vProbeSetsNameOrder;
    if (iWaveIndex == 0)
    {
        pv = new CNCychpProbeSetsCopyNumber[vProbeSets.getCount()];
        v.reserve(iProbeSetCount);
        vProbeSetsNameOrder.reserve(iProbeSetCount);
        for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
        {
            v.add(&pv[iProbeSetIndex]);
            vProbeSetsNameOrder.add(vProbeSets[iProbeSetIndex]);
        }
        vProbeSetsNameOrder.quickSort(0);
        m_uiPassedQCCount = uiCelCount;
    }
    if (bCNReferenceRun == false && m_strSelectedQC != "")
    {
        cychp.setSelectedQC(m_strSelectedQC);
    }
    Matrix mxAutosomes(iAutosomeProcessCount, m_uiPassedQCCount);
    unsigned int uiSampleIndex = 0;
    AffxString strCNReferenceInputFileName = "";
    if (m_pEngine->isOptDefined("cn-reference-input")) {strCNReferenceInputFileName = Fs::basename(m_pEngine->getOpt("cn-reference-input"));}
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        Verbose::out(1, "CNReferenceMethodWaveCorrection::loadPCALikeVector(...) Step " + ToStr(iWaveIndex) + " Loading data " +ToStr(uiCelIndex + 1) + " of " + ToStr(uiCelCount));
        if (iWaveIndex == 0)
        {
            // Load up the signals, and calculate log2 ratios.
            if (bCNReferenceRun)
            {
                double dA = 0;
                double dB = 0;
                int iColumnIndex = 0;
                int rowIndex = -1;
                tsv5 = group5->openTsv("Signals-" + ToStr(uiCelIndex), affx::FILE5_OPEN_RO);
                for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
                {
                    CNWaveCorrectionProbeSet* p = vProbeSetsNameOrder[iProbeSetIndex];
                    tsv5->nextLine();
                    while (p->Order != ++rowIndex)
                    {
                        if (p->Order < rowIndex) {Err::errAbort("Logical error in processing Signals-" + ToStr(uiCelIndex) + " - rowIndex > p->Order.");}
                        if (tsv5->nextLine() != affx::FILE5_OK) {Err::errAbort("Logical error in processing Signals-" + ToStr(uiCelIndex) + " - it has too few rows.");}
                    }

                    iColumnIndex = (uiCelIndex * 2);
                    tsv5->get(0, 0, &dA);
                    tsv5->get(0, 1, &dB);
                    double dSignal = (dA + dB);
                    double dLog2Ratio = 0;
                    if (p->Chromosome == m_iXChromosome)
                    {
                        dLog2Ratio = log2(dSignal / p->XXMedianSignal);
                    }
                    else if (p->Chromosome == m_iYChromosome)
                    {
                        dLog2Ratio = log2(dSignal / p->YMedianSignal);
                    }
                    else // Autosome
                    {
                        dLog2Ratio = log2(dSignal / p->MedianSignal);
                    }
                    p->Signal = dSignal;
                    p->Log2Ratio = dLog2Ratio;
                }
                tsv5->close();
                delete tsv5; tsv5 = NULL;

                covariateAdjustLog2Ratio(vProbeSets, uiCelIndex);

                calibrateMedianAutosomeMedian(vProbeSets);

                if (getEngine()->getOptBool("keep-intermediate-data"))
                {
                    string fileName = Fs::basename(getEngine()->getOptVector("cels")[uiCelIndex]) + ".wave" + ::ToStr(iWaveIndex) + ".a5";
                    ::writeWaveCorrectionLog2Ratios(getEngine()->getOpt("out-dir"), fileName, vProbeSets);
                }

                bool bXX = false;
                bool bY = false;
                calculateXY(vProbeSets, iProbeSetCount, bXX, bY);
                m_vHasXX[uiCelIndex] = bXX;
                m_vHasY[uiCelIndex] = bY;
                m_vYTargets[uiCelIndex] = 0;
                m_vPassedQC[uiCelIndex] = true;
            }
            else
            {
                cychp.readFile(m_pEngine->getOptVector("cychps")[uiCelIndex], v, true, true);
                if ((!m_bForce) && (Fs::basename(cychp.getCychpHeader().getCNReferenceFileName()) != strCNReferenceInputFileName))
                {
                    group5->close();
                    delete group5;
                    file5.close();
                    Fs::rm(m_pEngine->getOpt("temp-reference-file"), false);
                    Err::errAbort("The CN reference used to generate the CYCHP file does not match the cn-reference-input file name.\tCYCHP CN reference = " + Fs::basename(cychp.getCychpHeader().getCNReferenceFileName()));
                }
                if ((!m_bCyto2) && (uiCelIndex == 0))
                {
                    vProbeSets.quickSort(0);
                    v.quickSort(0);
                }
                Verbose::out(1, "Gender\t" + cychp.getCychpHeader().getGenderAsString());
                for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
                {
                    CNWaveCorrectionProbeSet* p1 = vProbeSets[iProbeSetIndex];
                    CNCychpProbeSetsCopyNumber* p2 = v[iProbeSetIndex];
                    if (p1->ProbeSetName != p2->getProbeSetName())
                    {
                        group5->close();
                        delete group5;
                        file5.close();
                        Fs::rm(m_pEngine->getOpt("temp-reference-file"), false);
                        Err::errAbort("Probe Set sets mismatch between CN reference and CYCHP file.");
                    }
                    p1->Log2Ratio = p2->getLog2Ratio();
                    p1->Chromosome = p2->getChromosome();
                    p1->Position = p2->getPosition();
                }
                if ((!m_bCyto2) && (uiCelIndex == 0))
                {
                    vProbeSets.quickSort(1);
                    v.quickSort(1);
                }

                m_vHasXX[uiCelIndex] = (cychp.getCychpHeader().getXX());
                m_vHasY[uiCelIndex] = (cychp.getCychpHeader().getY());
                m_vYTargets[uiCelIndex] = cychp.getCychpHeader().getYTarget();
                m_vPassedQC[uiCelIndex] = true;
                if (m_bCyto2)
                {
                    if ((cychp.getCychpHeader().getMAPD() > m_dCNQCCutoff) || (cychp.getCychpHeader().getPVQC() < m_dSNPQCCutoff))
                    {
                        m_vPassedQC[uiCelIndex] = false;
                        Verbose::out(1, "CNQC = " + ::getDouble(cychp.getCychpHeader().getMAPD(), 6) + ", SNPQC = " + ::getDouble(cychp.getCychpHeader().getPVQC(), 6) + ", CYCHP file failed QC.");
                        m_uiPassedQCCount--;
                    }
                    else
                    {
                        Verbose::out(1, "CNQC = " + ::getDouble(cychp.getCychpHeader().getMAPD(), 6) + ", SNPQC = " + ::getDouble(cychp.getCychpHeader().getPVQC(), 6) + ", CYCHP file passed QC.");

                        if (m_bUseHighWavinessSegCount)
                        {
                            if ((cychp.getCychpHeader().getWavinessSegCount() <= m_iWavinessSegCountCutoff) && (cychp.getCychpHeader().getWavinessSegCount() != -1))
                            {
                                m_vPassedQC[uiCelIndex] = false;
                                Verbose::out(1, "WavinessSegCount " + ::getInt(cychp.getCychpHeader().getWavinessSegCount()) + " is <= " + ::getInt(m_iWavinessSegCountCutoff) + ", expecting wavy samples only, this sample is skipped because it is non-wavy.");
                                m_uiPassedQCCount--;

                            }
                            else
                            {
                                Verbose::out(1, "WavinessSegCount " + ::getInt(cychp.getCychpHeader().getWavinessSegCount()) + " is > " + ::getInt(m_iWavinessSegCountCutoff) + ".");
                            }
                        }
                        else
                        {
                            if ((cychp.getCychpHeader().getWavinessSegCount() > m_iWavinessSegCountCutoff) && (cychp.getCychpHeader().getWavinessSegCount() != -1))
                            {
                                m_vPassedQC[uiCelIndex] = false;
                                Verbose::out(1, "WavinessSegCount " + ::getInt(cychp.getCychpHeader().getWavinessSegCount()) + " is > " + ::getInt(m_iWavinessSegCountCutoff) + ", expecting non-wavy samples only, this sample is skipped because it is wavy.");
                                m_uiPassedQCCount--;

                            }
                            else
                            {
                                Verbose::out(1, "WavinessSegCount " + ::getInt(cychp.getCychpHeader().getWavinessSegCount()) + " is <= " + ::getInt(m_iWavinessSegCountCutoff) + ".");
                            }
                        }
                    }
                }
                else
                {
                    if ((cychp.getCychpHeader().getMAPD() > m_dCNQCCutoff) || ((cychp.getCychpHeader().getCQC() != -1) && (cychp.getCychpHeader().getCQC() < m_dSNPQCCutoff)))
                    {
                        m_vPassedQC[uiCelIndex] = false;
                        Verbose::out(1, "CNQC = " + ::getDouble(cychp.getCychpHeader().getMAPD(), 6) + ", SNPQC = " + ::getDouble(cychp.getCychpHeader().getCQC(), 6) + ", CYCHP file failed QC.");
                        m_uiPassedQCCount--;
                    }
                    else
                    {
                        Verbose::out(1, "CNQC = " + ::getDouble(cychp.getCychpHeader().getMAPD(), 6) + ", SNPQC = " + ::getDouble(cychp.getCychpHeader().getCQC(), 6) + ", CYCHP file passed QC.");

                        if (m_bUseHighWavinessSegCount)
                        {
                            if ((cychp.getCychpHeader().getWavinessSegCount() <= m_iWavinessSegCountCutoff) && (cychp.getCychpHeader().getWavinessSegCount() != -1))
                            {
                                m_vPassedQC[uiCelIndex] = false;
                                Verbose::out(1, "WavinessSegCount " + ::getInt(cychp.getCychpHeader().getWavinessSegCount()) + " is <= " + ::getInt(m_iWavinessSegCountCutoff) + ", expecting wavy samples only, this sample is skipped because it is non-wavy.");
                                m_uiPassedQCCount--;

                            }
                            else
                            {
                                Verbose::out(1, "WavinessSegCount " + ::getInt(cychp.getCychpHeader().getWavinessSegCount()) + " is > " + ::getInt(m_iWavinessSegCountCutoff) + ".");
                            }
                        }
                        else
                        {
                            if ((cychp.getCychpHeader().getWavinessSegCount() > m_iWavinessSegCountCutoff) && (cychp.getCychpHeader().getWavinessSegCount() != -1))
                            {
                                m_vPassedQC[uiCelIndex] = false;
                                Verbose::out(1, "WavinessSegCount " + ::getInt(cychp.getCychpHeader().getWavinessSegCount()) + " is > " + ::getInt(m_iWavinessSegCountCutoff) + ", expecting non-wavy samples only, this sample is skipped because it is wavy..");
                                m_uiPassedQCCount--;

                            }
                            else
                            {
                                Verbose::out(1, "WavinessSegCount " + ::getInt(cychp.getCychpHeader().getWavinessSegCount()) + " is <= " + ::getInt(m_iWavinessSegCountCutoff) + ".");
                            }
                        }
                    }
                }
                // Subtract the y-target (first wave only)
                for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
                {
                    CNWaveCorrectionProbeSet* p = vProbeSets[iProbeSetIndex];
                    if (p->Chromosome == m_iYChromosome && m_vPassedQC[uiCelIndex] && m_vHasY[uiCelIndex])
                    {
                        p->Log2Ratio -= m_vYTargets[uiCelIndex];
                    }
                }
                m_vYTargets[uiCelIndex] = 0;
            }
        }
        else
        {
            float f = 0;
            affx::File5_Tsv* ptsv5 = group5->openTsv("Log2Ratios-" + ::getInt(uiCelIndex), affx::FILE5_OPEN_RO);
            for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
            {
                CNWaveCorrectionProbeSet* p = vProbeSets[iProbeSetIndex];
                ptsv5->nextLine();
                ptsv5->get(0, 0, &f);
                p->Log2Ratio = f;
            }
            ptsv5->close();
            delete ptsv5;

            adjustLog2Ratios(vProbeSets, iProbeSetCount, iWaveIndex, uiCelIndex);
            calibrateMedianAutosomeMedian(vProbeSets);
        }

        tsv5 = group5->openTsv("Log2Ratios-" + ::getInt(uiCelIndex), affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "Value", affx::FILE5_DTYPE_FLOAT);
        float f = 0;
        for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
        {
            f = vProbeSets[iProbeSetIndex]->Log2Ratio;
            tsv5->set_f(0, 0, f);
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        if (iWaveIndex > 0)
        {
            if (m_vPassedQC[uiCelIndex])
            {
                loadMatrix(vProbeSets, vAutosomeProbeSetIndexes, uiSampleIndex, mxAutosomes);
                uiSampleIndex++;
            }
            else
            {
                Verbose::out(1, "Skipping this file because it failed QC, or because of the waviness-seg-count.");
            }
        }
    }

    if (iWaveIndex == 0)
    {
        Verbose::out(1, "Using " + ::getUnsignedInt(m_uiPassedQCCount) + " out of " + ::getUnsignedInt(uiCelCount) + " samples for this run.");
        if (((double)m_uiPassedQCCount / (double)uiCelCount) < 0.5)
        {
            Verbose::warn(1, "The percentage of samples failing due to QC or waviness-seg-count is >= 50%.");
        }
        v.nullAll();
        vProbeSetsNameOrder.nullAll();
        delete[] pv;
        if (m_uiPassedQCCount == 0) {Err::errAbort("No samples are valid for this run.");}

        mxAutosomes.ReSize(iAutosomeProcessCount, m_uiPassedQCCount);
        unsigned int uiSampleIndex = 0;
        for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
        {
            if (m_vPassedQC[uiCelIndex])
            {
                float f = 0;
                affx::File5_Tsv* ptsv5 = group5->openTsv("Log2Ratios-" + ::getInt(uiCelIndex), affx::FILE5_OPEN_RO);
                for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
                {
                    CNWaveCorrectionProbeSet* p = vProbeSets[iProbeSetIndex];
                    ptsv5->nextLine();
                    ptsv5->get(0, 0, &f);
                    p->Log2Ratio = f;
                }
                ptsv5->close();
                delete ptsv5;
                loadMatrix(vProbeSets, vAutosomeProbeSetIndexes, uiSampleIndex, mxAutosomes);
                uiSampleIndex++;
            }
        }
    }

    if (m_bDemean)
    {
        for (int y = 0; (y < mxAutosomes.Nrows()); y++)
        {
            double dAutosomeMean = 0;
            for (int x = 0; (x < mxAutosomes.Ncols()); x++)
            {
                dAutosomeMean += mxAutosomes.element(y, x);
            }
            dAutosomeMean /= (double)mxAutosomes.Ncols();
            for (int x = 0; (x < mxAutosomes.Ncols()); x++)
            {
                mxAutosomes.element(y, x) -= dAutosomeMean;
            }
        }
    }

    Matrix mxAutosomePrincipalComponent(vAutosomeProbeSetIndexes.size(), 1);
    Matrix mxAutosomeLog2Ratios(vAutosomeProbeSetIndexes.size(), 1);
    Matrix mxU; DiagonalMatrix mxS;

    SVD(mxAutosomes, mxS, mxU);

    CNAnalysisMethod::memory("CNReferenceMethodWaveCorrection::svd(...)");
    for (int y = 0; (y < mxAutosomes.Nrows()); y++)
    {
        mxAutosomePrincipalComponent.element(y, 0) = mxU.element(y, 0);
    }
//    Verbose::out(1, "*");
    uiSampleIndex = 0;
    for (int x = 0; (x < uiCelCount); x++)
    {
        if (m_vPassedQC[x])
        {
            for (int y = 0; (y < mxAutosomes.Nrows()); y++)
            {
                mxAutosomeLog2Ratios.element(y, 0) = mxAutosomes.element(y, uiSampleIndex);
            }
            m_vAutosomeCV[x] = corr(mxAutosomePrincipalComponent, mxAutosomeLog2Ratios);
    //        Verbose::out(1, ::getInt(x + 1) + "\t" + ::getDouble(m_vAutosomeCV[x], 6));
            uiSampleIndex++;
        }
    }
//    Verbose::out(1, "*");
    Verbose::out(1, "CNReferenceMethodWaveCorrection::loadPCALikeVector(...) Step " + ToStr(iWaveIndex) + " Calculating wave removal");
    affx::File5_Tsv** pptsv5 = new affx::File5_Tsv*[uiCelCount];
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5[uiCelIndex] = group5->openTsv("Log2Ratios-" + ::getInt(uiCelIndex), affx::FILE5_OPEN_RO);
    }
    Matrix mx(iProbeSetCount, 1);
    float f = 0;
    double* pValues = new double[uiCelCount];
    for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
    {
        CNWaveCorrectionProbeSet* p = vProbeSets[iProbeSetIndex];
        int iCount = 0;
        if (p->Chromosome == m_iXChromosome)
        {
            for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
            {
                if (m_vPassedQC[uiCelIndex])
                {
                    pptsv5[uiCelIndex]->nextLine();
                    pptsv5[uiCelIndex]->get(0, 0, &f);
                    if (m_vAutosomeCV[uiCelIndex] < 0) {f *= -1;}
                    if (m_vHasXX[uiCelIndex]) {pValues[iCount] = f; iCount++;}
                }
            }
        }
        else if (p->Chromosome == m_iYChromosome)
        {
            for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
            {
                if (m_vPassedQC[uiCelIndex])
                {
                    pptsv5[uiCelIndex]->nextLine();
                    pptsv5[uiCelIndex]->get(0, 0, &f);
                    if (m_vHasY[uiCelIndex])
                    {
                        f -= m_vYTargets[uiCelIndex];
                        if (m_vAutosomeCV[uiCelIndex] < 0) {f *= -1;}
                        pValues[iCount] = f; iCount++;
                    }
                }
            }
        }
        else // Autosome
        {
            for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
            {
                if (m_vPassedQC[uiCelIndex])
                {
                    pptsv5[uiCelIndex]->nextLine();
                    pptsv5[uiCelIndex]->get(0, 0, &f);
                    if (m_vAutosomeCV[uiCelIndex] < 0) {f *= -1;}
                    pValues[iCount] = f; iCount++;
                }
            }
        }
        if (iCount == 0) {mx.element(iProbeSetIndex, 0) = 0;}
        else
        {
            double dLowPercentile = percentile((1.0 - m_dPercentile) * 100, pValues, iCount);
            double dHighPercentile = percentile(m_dPercentile * 100, pValues, iCount);

            double dDiff = fabs(fabs(dLowPercentile) - fabs(dHighPercentile));
            if ((dDiff < m_dWaveCutoff) || (fabs(dHighPercentile) > fabs(dLowPercentile)))
            {
                mx.element(iProbeSetIndex, 0) = dHighPercentile;
            }
            else
            {
                mx.element(iProbeSetIndex, 0) = dLowPercentile;
            }
            if (mx.element(iProbeSetIndex, 0) != mx.element(iProbeSetIndex, 0)) {mx.element(iProbeSetIndex, 0) = 0;}
        }
    }
    delete[] pValues;
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5[uiCelIndex]->close();
        delete pptsv5[uiCelIndex];
    }
    delete[] pptsv5;

    group5->close();
    delete group5;
    file5.close();

    // Normalize
    double dNorm = norm(mx);
    mx /= dNorm;

    // Orthogonalize
    int iTotalWaveIndex = m_iRefWaveCount + iWaveIndex;
    std::vector<double> vSums;
    vSums.resize(iProbeSetCount);
    for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
    {
        vSums[iProbeSetIndex] = 0;
    }
    for (int iWave = 0; (iWave < iTotalWaveIndex); iWave++)
    {
        double dSum = 0;
        for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
        {
            dSum += vProbeSets[iProbeSetIndex]->Waves[iWave] * mx.element(iProbeSetIndex, 0);
        }
        for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
        {
            vSums[iProbeSetIndex] += dSum * vProbeSets[iProbeSetIndex]->Waves[iWave];
        }
    }
    for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
    {
        mx.element(iProbeSetIndex, 0) -= vSums[iProbeSetIndex];
    }
    // Normalize
    dNorm = norm(mx);
    mx /= dNorm;

    for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
    {
        vProbeSets[iProbeSetIndex]->Waves[iTotalWaveIndex] = mx.element(iProbeSetIndex, 0);
    }
}

// TODO: This needs to be templatized and moved to CNAnalysisMethod.h
// to avoid cloning CNAnalysisMethodCovariateLRAdjuster::run()
//
void CNReferenceMethodWaveCorrection::covariateAdjustLog2Ratio(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, unsigned int uiCelIndex)
{
    if (!m_pEngine->isOptDefined("lr-adjustment-covariates")) {return;}
    if (m_pEngine->getOpt("lr-adjustment-covariates") == "none") {return;}
    if (m_pEngine->getOpt("lr-adjustment-covariates") == "") {return;}
    int numCovariates = CovariateParams::m_vLRCovariates.size();
    if (numCovariates == 0)
    {
        return;
    }

    // Check if any covariate values were read from the annotation file or covariates file
    if (vProbeSets.size() > 0 && vProbeSets[0]->AllCovariates.size() == 0)
    {
        if (m_pEngine->getOptBool("disable-covariates-file-warning") == false)
        {
            Verbose::warn(1, "Covariate adjustment is disabled. No covariate values were loaded. annotation-file does not support covariates.");
            m_pEngine->setOpt("disable-covariates-file-warning", "true");
        }
        return;
    }

    const float IQRepsilon = 1.0e-10;

    Verbose::out(1, "CNReferenceMethodWaveCorrection::covariateAdjustLog2Ratio()...");
    if (m_pvExperiments != NULL) m_pobjExperiment = m_pvExperiments->getAt(uiCelIndex);
    int iXChromosome = m_pEngine->getOptInt("xChromosome");
    int iYChromosome = m_pEngine->getOptInt("yChromosome");
    int validProcessFlags[] = { 1, 3 };      // CN and SNP, resp.

    // Find out if marker-class is one of the covariates
    int markerClassIndex = CovariateParams::mapCovariateNameToIndex("marker-class");

    if (m_pEngine->getOptBool("keep-intermediate-data")) {
        writeWaveCorrectionLog2Ratios("WCBeforeCovariateAdjust", vProbeSets);
    }

    for (int iCovariateIndex = 0; iCovariateIndex < numCovariates; iCovariateIndex++)
    {
        int iTransIndex = CovariateParams::m_vLRCovariates[iCovariateIndex];
        int maxFlagIndex = 2;      // the default number of marker classes (i.e. SNP and CN)
        if (markerClassIndex == iTransIndex)
        {
            // This is the marker-class covariate, so do not distinguish SNP and CN in what follows
            maxFlagIndex = 1;      // do the loop below only once (i.e. it's one class of markers)
        }
        for (int iFlagIndex = 0; iFlagIndex < maxFlagIndex; iFlagIndex++)
        {
            string intermSuffix;
            if (maxFlagIndex == 1) {
                intermSuffix = "_MARKER_CLASS";
            } else if (validProcessFlags[iFlagIndex] == 1) {
                intermSuffix = "_CN";
            } else if (validProcessFlags[iFlagIndex] == 3) {
                intermSuffix = "_SNP";
            } else {
                intermSuffix = "_INVALID_FLAG";
            }

            vector<float> covariates;
            vector<int> probeSetsToProcess;
    
            // NB: local probe sets of CNAnalysisMethodLog2Ratio (process flag = 1 or 3)
            for (int iPSIndex = 0; iPSIndex < vProbeSets.size(); iPSIndex++)
            {
                CNWaveCorrectionProbeSet *pset = vProbeSets[iPSIndex];
                float covValue = pset->getCovariateValue(iTransIndex);
    
                if ((maxFlagIndex == 1 || pset->getProcessFlag() == validProcessFlags[iFlagIndex]) &&
                    covValue == covValue)
                {
                    covariates.push_back(covValue);
                    probeSetsToProcess.push_back(iPSIndex);
                }
            }
            vector<int> binning(covariates.size());
            int numBins = CNAnalysisMethodCovariateLRAdjuster::determineBins(covariates, binning, iCovariateIndex);
    
            // Adjust log2 ratios
            vector<float> mediansLR(numBins);
            vector<float> IQRanges(numBins);
            for (int iBinIndex = 0; iBinIndex < numBins; iBinIndex++)
            {
                // Find bin median of autosome log2 ratios
                vector<float> log2ratios;
                for (int iIndex = 0; iIndex < probeSetsToProcess.size(); iIndex++)
                {
                    if (binning[iIndex] == iBinIndex)
                    {
                        CNWaveCorrectionProbeSet *pset = vProbeSets[probeSetsToProcess[iIndex]];
                        if (pset->getChromosome() != iXChromosome && pset->getChromosome() != iYChromosome)
                        {
                            if (pset->getLog2Ratio() == pset->getLog2Ratio())
                            {
                                log2ratios.push_back(pset->getLog2Ratio());
                            }
                        }
                    }                        
                }
                int binSize = log2ratios.size();
                if (binSize == 0)
                {
                    mediansLR[iBinIndex] = 0.0;
                    IQRanges[iBinIndex] = numeric_limits<float>::quiet_NaN();   // invalid value to flag an empty bin
                    continue;
                }
                std::nth_element(log2ratios.begin(), log2ratios.begin() + binSize/2, log2ratios.end());
                mediansLR[iBinIndex] = log2ratios[binSize/2];
    
                std::nth_element(log2ratios.begin(), log2ratios.begin() + binSize/4, log2ratios.end());
                float lowerQuartile = log2ratios[binSize/4];
    
                std::nth_element(log2ratios.begin(), log2ratios.begin() + 3*binSize/4, log2ratios.end());
                IQRanges[iBinIndex] = log2ratios[3*binSize/4] - lowerQuartile;
            }
        
            bool subtractMedianFromXY = (CovariateParams::m_vLRSubractFromXY[iCovariateIndex] == "on");
    
            // Subtract the medians from markers in this bin (including X and Y)
            for (int iIndex = 0; iIndex < probeSetsToProcess.size(); iIndex++)
            {
                CNWaveCorrectionProbeSet *pset = vProbeSets[probeSetsToProcess[iIndex]];
                if (subtractMedianFromXY || (pset->getChromosome() != iXChromosome && pset->getChromosome() != iYChromosome))
                {
                    if (pset->getLog2Ratio() == pset->getLog2Ratio())
                    {
                        pset->setLog2Ratio(pset->getLog2Ratio() - mediansLR[binning[iIndex]]);
                    }
                }
            }
    
            // If requested, perform IQR scaling
            if (CovariateParams::m_vIQRScaling[iCovariateIndex] == "on")
            {
                vector<float> validIQRanges;
                for (int iBinIndex = 0; iBinIndex < numBins; iBinIndex++)
                {
                    if (IQRanges[iBinIndex] == IQRanges[iBinIndex])
                    {
                        validIQRanges.push_back(IQRanges[iBinIndex]);
                    }
                }
                int numValidRanges = validIQRanges.size();
                std::nth_element(validIQRanges.begin(), validIQRanges.begin() + numValidRanges/2, validIQRanges.end());
                float medianIQRange = validIQRanges[numValidRanges/2];
    
                vector<float> scalingFactors(numBins);
                for (int iBinIndex = 0; iBinIndex < numBins; iBinIndex++)
                {
                    if(IQRanges[iBinIndex] == IQRanges[iBinIndex])
                    {
                        if (IQRanges[iBinIndex] < IQRepsilon)
                        {
                            scalingFactors[iBinIndex] = 1.0;
                        }
                        else {
                            scalingFactors[iBinIndex] = medianIQRange/IQRanges[iBinIndex];
                        }
                    }
                }
    
                for (int iIndex = 0; iIndex < probeSetsToProcess.size(); iIndex++)
                {
                    CNWaveCorrectionProbeSet *pset = vProbeSets[probeSetsToProcess[iIndex]];
    
                    // Turn off iqr-scaling for X/Y if subtractMedianFromXY == false - per Ben 1/28/2011
                    if (subtractMedianFromXY || (pset->getChromosome() != iXChromosome && pset->getChromosome() != iYChromosome))
                    {
                        if (pset->getLog2Ratio() == pset->getLog2Ratio())
                        {
                            pset->setLog2Ratio(pset->getLog2Ratio() * scalingFactors[binning[iIndex]]);
                        }
                    }
                }
            }
            if (m_pEngine->getOptBool("keep-intermediate-data"))
            {
                writeWaveCorrectionLog2Ratios("WCAfterCovariateAdjust_" + Convert::toString(iCovariateIndex+1) + intermSuffix, vProbeSets);
                writeMediansVector("WCCovariateL2rMedians_" + Convert::toString(iCovariateIndex+1) + intermSuffix, mediansLR);
            }
        }
    }
}

// Subtract residuals from Log2 Ratios.
void CNReferenceMethodWaveCorrection::adjustLog2Ratios(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, int iProbeSetCount, int iWaveIndex, unsigned int uiCelIndex)
{
    double dSum = 0;
    double dTempNorm = 0;
    if (iWaveIndex > 0)
    {
        Matrix mxTemp(m_iAutosomeCount, 1);
        int iTempIndex = 0;
        for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
        {
            CNWaveCorrectionProbeSet* p = vProbeSets[iProbeSetIndex];
            if (p->Chromosome <= m_iLastAutosomeChromosome)
            {
                dSum += (p->Waves[m_iRefWaveCount + iWaveIndex - 1] * p->Log2Ratio);
                mxTemp.element(iTempIndex, 0) = p->Waves[m_iRefWaveCount + iWaveIndex - 1]; iTempIndex++;
            }
        }
        dTempNorm = norm(mxTemp);
        for (int iProbeSetIndex = 0; (iProbeSetIndex < iProbeSetCount); iProbeSetIndex++)
        {
                vProbeSets[iProbeSetIndex]->Log2Ratio -= (dSum * vProbeSets[iProbeSetIndex]->Waves[m_iRefWaveCount + iWaveIndex - 1] / dTempNorm);
        }

        if (m_pEngine->isOptDefined("wave-correction-log2ratio-adjustment-method")     &&
            m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") != "none" &&
            m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") != ""     &&
            getExperiment())
        {
            getExperiment()->addWavinessAmplitude(iWaveIndex, dSum/dTempNorm);
        }
    }
}

/**
 * Get the last autosome chromosome. As defined by being the last chromosome before the X chromosome value.
 * @return int - The chromosome number
 */
void CNReferenceMethodWaveCorrection::calculateChromosomeStatistics(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets)
{
    m_iXChromosomeCount = 0;
    m_iYChromosomeCount = 0;
    m_iAutosomeCount = 0;
    m_iLastAutosomeChromosome = 1;
    for (int iIndex = 0; (iIndex < vProbeSets.size()); iIndex++)
    {
        CNWaveCorrectionProbeSet* p = vProbeSets[iIndex];
        if (p->Chromosome == m_iXChromosome) {m_iXChromosomeCount++;}
        else if (p->Chromosome == m_iYChromosome) {m_iYChromosomeCount++;}
        else if (p->Chromosome < m_iXChromosome)
        {
            m_iAutosomeCount++;
            m_iLastAutosomeChromosome = Max(m_iLastAutosomeChromosome, (int)p->Chromosome);
        }
    }
}
/**
 * Calcuate the median autosome median. This is used to normalize the log2 ratios.
 */
void CNReferenceMethodWaveCorrection::calibrateMedianAutosomeMedian(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets)
{
    int iProbeSetCount = vProbeSets.size();
    AffxMultiDimensionalArray<int> vChrCounts(m_iLastAutosomeChromosome);
    AffxMultiDimensionalArray<int> vChrIndexes(m_iLastAutosomeChromosome);
    AffxMultiDimensionalArray<float> vChrMedians(m_iLastAutosomeChromosome);
    AffxMultiDimensionalArray<float>* pvar = new AffxMultiDimensionalArray<float>[m_iLastAutosomeChromosome];
    for (int iRowIndex = 0; (iRowIndex < iProbeSetCount); iRowIndex++)
    {
        CNWaveCorrectionProbeSet* pobjProbeSet = vProbeSets[iRowIndex];
        if (pobjProbeSet->Chromosome <= 0) {continue;}
        if (pobjProbeSet->Chromosome <= m_iLastAutosomeChromosome)
        {
            vChrCounts.increment(pobjProbeSet->Chromosome - 1);
        }
    }
    for (int iIndex = 0; (iIndex < m_iLastAutosomeChromosome); iIndex++)
    {
        pvar[iIndex].initialize(vChrCounts.get(iIndex));
    }
    vChrIndexes.initialize();
    for (int iRowIndex = 0; (iRowIndex < iProbeSetCount); iRowIndex++)
    {
        CNWaveCorrectionProbeSet* pobjProbeSet = vProbeSets[iRowIndex];
        if (pobjProbeSet->Chromosome <= 0) {continue;}
        if (pobjProbeSet->Chromosome <= m_iLastAutosomeChromosome)
        {
            pvar[pobjProbeSet->Chromosome - 1].set(vChrIndexes.get(pobjProbeSet->Chromosome - 1), pobjProbeSet->Log2Ratio);
            vChrIndexes.increment(pobjProbeSet->Chromosome - 1);
        }
    }
    for (int iIndex = 0; (iIndex < m_iLastAutosomeChromosome); iIndex++)
    {
        vChrMedians.set(iIndex, (float)pvar[iIndex].median());
    }
    double dMedianAutosomeMedian = vChrMedians.finiteMedian();
    for (int iRowIndex = 0; (iRowIndex < iProbeSetCount); iRowIndex++)
    {
        CNWaveCorrectionProbeSet* pobjProbeSet = vProbeSets[iRowIndex];
        pobjProbeSet->Log2Ratio = (float)(pobjProbeSet->Log2Ratio - dMedianAutosomeMedian);
    }
    delete[] pvar;
}

void CNReferenceMethodWaveCorrection::calculateXY(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, int iProbeSetCount, bool& bXX, bool& bY)
{
    AffxArray<AffxString> vYProbeSetsToProcess;

    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;

    if ( tsv.open(m_pEngine->getOpt("chrY-probes")) == affx::TSV_OK ) {
        AffxString strProbeSetName;
        while (tsv.nextLevel(0) == affx::TSV_OK) {
            tsv.get(0,1, strProbeSetName);
            if (!strProbeSetName.empty()) {
                vYProbeSetsToProcess.add(new AffxString(strProbeSetName));
            }
        }
        tsv.clear();
    }

    vYProbeSetsToProcess.quickSort(0);

    bXX = false;
    bY = false;
    // Allocate the data vectors.
    int iChrRefCount = 0;
    int iChrXCount = 0;
    int iChrYCount = 0;
    for (int uiRowIndex = 0; (uiRowIndex < iProbeSetCount); uiRowIndex++)
    {
        CNWaveCorrectionProbeSet* pobjProbeSet = vProbeSets[uiRowIndex];
        if ((pobjProbeSet->Chromosome == m_iYChromosome) && (vYProbeSetsToProcess.getCount() > 0) && (vYProbeSetsToProcess.binarySearch(pobjProbeSet->ProbeSetName, 0) == -1)) {continue;}
        if (pobjProbeSet->Chromosome == m_pEngine->getOptInt("reference-chromosome")) {iChrRefCount++;}
        else if (pobjProbeSet->Chromosome == m_iXChromosome) {iChrXCount++;}
        else if (pobjProbeSet->Chromosome == m_iYChromosome) {iChrYCount++;}
    }
    AffxMultiDimensionalArray<float> vChrRefMedian(iChrRefCount);
    AffxMultiDimensionalArray<float> vChrXSignal(iChrXCount);
    AffxMultiDimensionalArray<float> vChrYSignal(iChrYCount);

    int iChrRefIndex = 0;
    int iChrXIndex = 0;
    int iChrYIndex = 0;
    // For each probe set load the relevant data into the vectors so we can calculate the median.
    for (int iRowIndex = 0; (iRowIndex < iProbeSetCount); iRowIndex++)
    {
        CNWaveCorrectionProbeSet* pobjProbeSet = vProbeSets[iRowIndex];
        if ((pobjProbeSet->Chromosome == m_iYChromosome) && (vYProbeSetsToProcess.getCount() > 0) && (vYProbeSetsToProcess.binarySearch(pobjProbeSet->ProbeSetName, 0) == -1)) {continue;}
        if (pobjProbeSet->Chromosome == m_pEngine->getOptInt("reference-chromosome"))
        {
            vChrRefMedian.set(iChrRefIndex, pobjProbeSet->MedianSignal);
            iChrRefIndex++;
        }
        else if (pobjProbeSet->Chromosome == m_iXChromosome)
        {
            vChrXSignal.set(iChrXIndex, pobjProbeSet->Signal);
            iChrXIndex++;
        }
        else if (pobjProbeSet->Chromosome == m_iYChromosome)
        {
            vChrYSignal.set(iChrYIndex, pobjProbeSet->Signal);
            iChrYIndex++;
        }
    }
    vYProbeSetsToProcess.deleteAll();
    double dChrRefMedian = vChrRefMedian.median();

    // Set the experiment hasXX value.
    double dXXRatio = (vChrXSignal.median() / dChrRefMedian);
    bXX = ((dXXRatio > m_pEngine->getOptDouble("xx-cutoff")) && (dXXRatio < m_pEngine->getOptDouble("xx-cutoff-high")));

    // Set the experiment hasY value.
    double dYRatio = (vChrYSignal.median() / dChrRefMedian);
    bY = (dYRatio > m_pEngine->getOptDouble("y-cutoff"));

    Verbose::out(1, "CNReferenceMethodWaveCorrection::calculateXY(...)\tXXRatio\t" + ::getDouble(dXXRatio, 6) + "\tYRatio\t" + ::getDouble(dYRatio, 6) + "\t" + (bXX ? "XX" : "X") + (bY ? "Y" : ""));
}

void CNReferenceMethodWaveCorrection::loadMatrix(       AffxArray<CNWaveCorrectionProbeSet>& vProbeSets,
                                                        std::vector<int>& vProbeSetIndexes,
                                                        unsigned int uiSampleIndex,
                                                        Matrix& mx)
{
    for (int iIndex = 0; (iIndex < vProbeSetIndexes.size()); iIndex++)
    {
        float fLog2Ratio = vProbeSets[vProbeSetIndexes[iIndex]]->Log2Ratio;
        if (fLog2Ratio > m_dTrim) {fLog2Ratio = (float)m_dTrim;}
        else if (fLog2Ratio < -m_dTrim) {fLog2Ratio = (float)-m_dTrim;}
        if (fLog2Ratio != fLog2Ratio)
        {
            throw(Except("NaN Log2Ratios found in data."));
        }
        mx.element(iIndex, uiSampleIndex) = fLog2Ratio;
    }
}

double CNReferenceMethodWaveCorrection::corr(Matrix& mx1, Matrix& mx2)
{
    ColumnVector v1 = mx1.AsColumn();
    ColumnVector v2 = mx2.AsColumn();
    if (v1.Nrows() != v2.Nrows()) {return std::numeric_limits<double>::quiet_NaN();}
    int iLength = v1.Nrows();
    double dSum = 0;
    for (int iElementIndex = 0; (iElementIndex < iLength); iElementIndex++)
    {
        dSum += v1.element(iElementIndex);
    }
    double m1 = (dSum / (double)iLength);
    dSum = 0;
    for (int iElementIndex = 0; (iElementIndex < iLength); iElementIndex++)
    {
        dSum += v2.element(iElementIndex);
    }
    double m2 = (dSum / (double)iLength);
    double dNumerator = 0;
    double d1 = 0;;
    double d2 = 0;
    for (int iElementIndex = 0; (iElementIndex < iLength); iElementIndex++)
    {
        dNumerator += (v1.element(iElementIndex) - m1) * (v2.element(iElementIndex) - m2);
        d1 += (v1.element(iElementIndex) - m1) * (v1.element(iElementIndex) - m1) ;
        d2 += (v2.element(iElementIndex) - m2) * (v2.element(iElementIndex) - m2) ;
    }
    double dDenominator = sqrt(d1 * d2);
    return dNumerator / dDenominator;
}
/*
double CNReferenceMethodWaveCorrection::norm(Matrix& v)
{
    DiagonalMatrix D;
    SVD(v, D);
    return D.Maximum();
}
*/
void CNReferenceMethodWaveCorrection::updateReference(AffxArray<CNWaveCorrectionProbeSet>& vProbeSets, int iTotalProbeSetCount)
{
    int iWaveCount = m_iRefWaveCount + m_iWaveCount;
    if (m_pEngine->isOptDefined("cn-reference-input"))
    {
        Fs::fileCopy(m_pEngine->getOpt("cn-reference-input"), m_pEngine->getOpt("cn-reference-output"));
    }

    vProbeSets.quickSort(0);

    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;
    if (m_pEngine->isOptDefined("cn-reference-input"))
    {
        file5.open(m_pEngine->getOpt("cn-reference-output"), affx::FILE5_OPEN);
    }
    else
    {
        file5.open(m_pEngine->getOpt("reference-output"), affx::FILE5_OPEN);
    }
    group5 = file5.openGroup(m_strGroup, affx::FILE5_OPEN);

    tsv5 = group5->openTsv("WaveCorrection", affx::FILE5_REPLACE);
    for (int iWaveIndex = 0; (iWaveIndex < iWaveCount); iWaveIndex++)
    {
        tsv5->defineColumn(0, iWaveIndex, "X" + ToStr(iWaveIndex + 1), affx::FILE5_DTYPE_DOUBLE);
    }
    int iIndex = 0;
    CNWaveCorrectionProbeSet* p = NULL;
    for (int iFullIndex = 0; (iFullIndex < iTotalProbeSetCount); iFullIndex++)
    {
        if (iIndex < vProbeSets.size()) {p = vProbeSets[iIndex];}
        if ((p != NULL) && (p->Order == iFullIndex))
        {
            for (int iWaveIndex = 0; (iWaveIndex < (int)p->Waves.size()); iWaveIndex++)
            {
                tsv5->set_d(0, iWaveIndex, p->Waves[iWaveIndex]);
            }
            tsv5->writeLevel(0);
            iIndex++;
        }
        else
        {
            for (int iWaveIndex = 0; (iWaveIndex < iWaveCount); iWaveIndex++)
            {
                tsv5->set_d(0, iWaveIndex, std::numeric_limits<double>::quiet_NaN());
            }
            tsv5->writeLevel(0);
        }
    }

    tsv5->close();
    delete tsv5;

    if (m_pEngine->isOptDefined("cn-reference-input"))
    {
        tsv5 = group5->openTsv("AdditionalWavesParameters", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "Parameter", affx::FILE5_DTYPE_STRING, 1024);

        // dump options.
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-xChromosome=" + ::getInt(m_iXChromosome)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-yChromosome=" + ::getInt(m_iYChromosome)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-analysis=" + m_pEngine->getOpt("analysis")); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-cn-reference-input=" + Fs::basename(m_pEngine->getOpt("cn-reference-input"))); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-state-existing-wave-count=" + ::getInt(m_iRefWaveCount)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-cn-reference-output=" + Fs::basename(m_pEngine->getOpt("cn-reference-output"))); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-additional-wave-count=" + ::getInt(m_iWaveCount)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-trim=" + ::getDouble(m_dTrim, 6)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-percentile=" + ::getDouble(m_dPercentile, 6)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-demean=" + std::string(((m_bDemean) ? "true" : "false"))); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-cn-qc-cutoff=" + ::getDouble(m_dCNQCCutoff, 6)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-snp-qc-cutoff=" + ::getDouble(m_dSNPQCCutoff, 6)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-waviness-seg-count-cutoff=" + ::getInt(m_iWavinessSegCountCutoff)); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-use-high-waviness-seg-count=" + std::string(((m_bUseHighWavinessSegCount) ? "true" : "false"))); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-force=" + std::string(((m_bForce) ? "true" : "false"))); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-keep-temp-data=" + std::string(((m_bKeepTempData) ? "true" : "false"))); tsv5->writeLevel(0);
        for (unsigned int ui = 0; (ui < m_pEngine->getOptVector("cychps").size()); ui++)
        {
            if (m_vPassedQC[ui])
            {
                tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-cychps-" + ::getUnsignedInt(ui + 1) + "=" + Fs::basename(m_pEngine->getOptVector("cychps")[ui])); tsv5->writeLevel(0);
            }
        }

        tsv5->close();
        delete tsv5;
    }

    group5->close();
    delete group5;
    file5.close();
}

void CNReferenceMethodWaveCorrection::writeWaveCorrectionLog2Ratios(std::string l2rFileInfix, AffxArray<CNWaveCorrectionProbeSet>& vProbeSets)
{
  affx::File5_File file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5 = NULL;

  // Self preservation
  // TODO: Fix this.
  if (getExperiment() == NULL) return;

  std::string analysisString = "analysis";
  std::string strExperimentName = getExperiment()->getExperimentName();
  std::string fileName = Fs::join(getEngine()->getOpt("out-dir"),
                                  analysisString,
                                  strExperimentName + "." + l2rFileInfix + ".l2r.a5");
  
  std::string dirName = Fs::join(getEngine()->getOpt("out-dir"), analysisString);  


  Verbose::out(3, "Writing Log2Ratios for " + fileName);

  try{


    if(!Fs::dirExists(dirName))
    {
        Fs::mkdir(dirName, false);
    }

    file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    group5 = file5.openGroup("Signals", affx::FILE5_REPLACE);
    tsv5 = group5->openTsv("Signals", affx::FILE5_REPLACE);
    tsv5->defineColumn(0,0,"ProbeSetName", affx::FILE5_DTYPE_STRING,40);
    tsv5->defineColumn(0,1,"Chromosome", affx::FILE5_DTYPE_CHAR);
    tsv5->defineColumn(0,2,"Position", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0,3,"Log2Ratio", affx::FILE5_DTYPE_FLOAT);

    CNWaveCorrectionProbeSet* pProbeSet;
    for(int iProbeSetIndex=0;  iProbeSetIndex < vProbeSets.getCount(); iProbeSetIndex++)
    {
      pProbeSet = vProbeSets.getAt(iProbeSetIndex);
      tsv5->set_string(0,0, pProbeSet->ProbeSetName);
      tsv5->set_c(0,1, pProbeSet->Chromosome);
      tsv5->set_i(0,2, pProbeSet->Position);
      tsv5->set_f(0,3, (float)pProbeSet->getLog2Ratio());
      tsv5->writeLevel(0);
    }

    tsv5->close();
    delete tsv5;
    delete group5;
    file5.close();

  }
  catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump intensities to file."));}
}
