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
 * @file CNAnalysisMethodReference.cpp
 *
 * @brief This file contains the CNAnalysisMethodReference class members.
 */

//
#include "copynumber/CNAnalysisMethodReference.h"
//
#include "copynumber/Annotation.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNAnalysisMethodCovariateParams.h"
#include "copynumber/CNAnalysisMethodCovariateSignalAdjuster.h"
//
#include "chipstream/SketchQuantNormTran.h"
#include "chipstream/TsvReport.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/TmpFileFactory.h"
#include "util/Util.h"

//
using namespace affx;
/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodReference::explainSelf()
{
    CNAnalysisMethodReference obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodReference::getDefaultDocOptions()
{
    std::vector<SelfDoc::Opt> opts;

//    SelfDoc::Opt(name, type, value, default, min, max, description)

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
SelfCreate* CNAnalysisMethodReference::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodReference* pMethod = new CNAnalysisMethodReference();
    std::string strPrefix = getPrefix();

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodReference::CNAnalysisMethodReference()
{
    m_pvExperiments = NULL;

    m_dAASum=0.0;
    m_dABSum=0.0;
    m_dBBSum=0.0;
    m_dAASumSquares=0.0;
    m_dABSumSquares=0.0;
    m_dBBSumSquares=0.0;
    m_iAAIndex=0;
    m_iABIndex=0;
    m_iBBIndex=0;
    m_iWarningCount = 0;
}

void CNAnalysisMethodReference::runPart1()
{
    Verbose::out(1, "CNAnalysisMethodReference::runPart1(...) start");
    isSetup();

    affx::File5_File file5;
    affx::File5_Group* group5;
    AffxString strFileName = m_pEngine->getOpt("reference-file");
    if ( !Fs::dirExists(Fs::dirname(strFileName)) ) {
      Fs::mkdirPath(Fs::dirname(strFileName), false);
    }
    Fs::rmIfExists(strFileName,false);
    file5.open(strFileName, affx::FILE5_CREATE);
    if (m_pEngine->getOptBool("cyto2"))
    {
        group5 = file5.openGroup("Cyto2", affx::FILE5_CREATE);
    }
    else
    {
        group5 = file5.openGroup("CopyNumber", affx::FILE5_CREATE);
        group5->close();
        delete group5;
        group5 = file5.openGroup("CN5", affx::FILE5_CREATE);
    }
    group5->close();
    delete group5;
    file5.close();

    determineSketchProbeSets();

    // Note that the invocation of newProbes is the  newProbes in Reference, not newProbes in Chipstream.
    // All probes in the probesets in m_pProbeSets are loaded.
    newProbes();
    calculateReferenceSketch();
    newSketch();

     // Create temp file
    TmpFileFactory tempFactory;
    tempFactory.setTmpdir(m_pEngine->getOpt("temp-dir"));
    m_strTempFileName = tempFactory.genFilename_basic( "CNReference_", ".a5.tmp");

    // Verbose::out(1, ::getInt(m_pvProbes->getCount()) + " Probes.");
    file5.open(m_strTempFileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    group5 = file5.openGroup("Intensities", affx::FILE5_REPLACE);
    group5->close();
    delete group5;
    file5.close();

    // process experiments
    for (int iExperimentIndex = 0; (iExperimentIndex < getExperiments()->getCount()); iExperimentIndex++)
    {
        Verbose::out(1, "CNAnalysisMethodReference::normalizeIntensities(...)" );
        m_pobjExperiment = getExperiments()->getAt(iExperimentIndex);
        normalizeIntensities();
        if (m_pEngine->getOptBool("keep-intermediate-data"))
        {
                writeIntensities("reference-normalized");
        }

        // write out intensities
        m_pvProbes->quickSort(3); // by ProbeSetIndex, Allele, ProbeID
        file5.open(m_strTempFileName, affx::FILE5_OPEN);
        group5 = file5.openGroup("Intensities", affx::FILE5_OPEN);
        affx::File5_Tsv* tsv5 = group5->openTsv("Intensities-" + ::getInt(getExperiment()->getIndex()), affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "ProbeSetIndex", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 1, "Allele", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "ProbeID", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "Intensity", affx::FILE5_DTYPE_FLOAT);
        for (int iIndex = 0; (iIndex < m_pvProbes->getCount()); iIndex++)
        {
            CNProbe* p = m_pvProbes->getAt(iIndex);
            tsv5->set_i(0, 0, p->getProbeSetIndex());
            tsv5->set_i(0, 1, p->getAllele());
            tsv5->set_i(0, 2, p->getProbeID());
            tsv5->set_f(0, 3, p->getIntensity());
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;
        group5->close();
        delete group5;
        file5.close();
    }

    // Read back out the intensities, calculate median intensity and write out probes to cn reference
    affx::File5_Tsv* tsv5ProbeEffects = NULL;
    file5.open(strFileName, affx::FILE5_OPEN);
    if (m_pEngine->getOptBool("cyto2"))
    {
        group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
        tsv5ProbeEffects = group5->openTsv("ProbeEffects", affx::FILE5_OPEN|affx::FILE5_CREATE);
    }
    else
    {
        group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        tsv5ProbeEffects = group5->openTsv("CN5.plier.feature-response", affx::FILE5_OPEN|affx::FILE5_CREATE);
    }

    tsv5ProbeEffects->defineColumn(0, 0, "probe_id", affx::FILE5_DTYPE_INT);
    tsv5ProbeEffects->defineColumn(0, 1, "feature_response", affx::FILE5_DTYPE_DOUBLE);
    tsv5ProbeEffects->defineColumn(0, 2, "ProbeSetIndex", affx::FILE5_DTYPE_INT);
    tsv5ProbeEffects->defineColumn(0, 3, "AlleleCode", affx::FILE5_DTYPE_INT);
    tsv5ProbeEffects->defineColumn(0, 4, "MedianNormalizedIntensity", affx::FILE5_DTYPE_FLOAT);

    unsigned int uiCelCount = getExperiments()->getCount();
    unsigned int uiProbeCount = m_pvProbes->getCount();
    affx::File5_File file5Input;
    affx::File5_Group* pGroup5Input;
    affx::File5_Tsv** pptsv5 = new affx::File5_Tsv*[uiCelCount];
    file5Input.open(m_strTempFileName, affx::FILE5_OPEN);
    pGroup5Input = file5Input.openGroup("Intensities", affx::FILE5_OPEN);
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5[uiCelIndex] = pGroup5Input->openTsv("Intensities-" + ::getInt(uiCelIndex), affx::FILE5_OPEN);
    }
    int iProbeID = 0;
    CNProbe objSearch;
    AffxMultiDimensionalArray<float> vIntensities(uiCelCount);
    m_pvProbes->quickSort(1); // ProbeID

    float f = 0;
    for (unsigned int uiProbeIndex = 0; (uiProbeIndex < uiProbeCount); uiProbeIndex++)
    {
        for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
        {
            pptsv5[uiCelIndex]->nextLine();
            pptsv5[uiCelIndex]->get(0, 2, &iProbeID);
            pptsv5[uiCelIndex]->get(0, 3, &f); vIntensities.set(uiCelIndex, f);
        }
        double dMedianIntensity = vIntensities.median();
        objSearch.setProbeID(iProbeID);
        int iSearchIndex = m_pvProbes->binarySearch(objSearch, 1);
        if (iSearchIndex != -1)
        {
            CNProbe* pobjProbe = m_pvProbes->getAt(iSearchIndex);
            pobjProbe->setMedianIntensity((float)dMedianIntensity);
            tsv5ProbeEffects->set_i(0, 0, pobjProbe->getProbeID());
            tsv5ProbeEffects->set_d(0, 1, 0); // Place holder for probe effects
            tsv5ProbeEffects->set_i(0, 2, pobjProbe->getProbeSetIndex());
            tsv5ProbeEffects->set_i(0, 3, (int)pobjProbe->getAllele());
            tsv5ProbeEffects->set_f(0, 4, (float)dMedianIntensity);
            tsv5ProbeEffects->writeLevel(0);
        }
    }
    // Close the input file
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5[uiCelIndex]->close();
        delete pptsv5[uiCelIndex];
    }
    delete[] pptsv5;
    pGroup5Input->close();
    delete pGroup5Input;
    file5Input.close();

    // close the output file
    tsv5ProbeEffects->close();
    delete tsv5ProbeEffects; tsv5ProbeEffects = NULL;
    group5->close();
    delete group5;
    file5.close();

    Verbose::out(1, "CNAnalysisMethodReference::runPart1(...) end");
}

// run pdnn here (well, back in cyto engine really)

void CNAnalysisMethodReference::runPart2()
{
    Verbose::out(1, "CNAnalysisMethodReference::runPart2(...) start");
    if ((m_pEngine->getOpt("snp-reference-output-file") != "") || (!m_pEngine->getOptBool("cyto2")))
    {
        m_mxGenotypeCalls.initialize(getExperiments()->getCount(), getProbeSets()->getCount());
        for (int x = 0; (x < getExperiments()->getCount()); x++)
        {
            for (int y = 0; (y < getProbeSets()->getCount()); y++)
            {
                m_mxGenotypeCalls.set(x, y, -2);
            }
        }
    }
    CNAnalysisMethodChipstream::newProbes();
    Verbose::out(1, ::getInt(m_pvProbes->getCount()) + " Probes.");
    newPlier(getExperiments()->getCount());

    m_pvProbes->quickSort(1); // ProbeID
    unsigned int uiCelCount = getExperiments()->getCount();
    unsigned int uiProbeCount = m_pvProbes->getCount();

    affx::File5_File file5Input;
    affx::File5_Group* pGroup5Input;
    affx::File5_Tsv* ptsv5;
    file5Input.open(m_strTempFileName, affx::FILE5_OPEN);
    pGroup5Input = file5Input.openGroup("Intensities", affx::FILE5_OPEN);
    int iProbeID = 0;
    float fIntensity = 0;
    CNProbe objSearch;
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        m_pobjExperiment = getExperiments()->getAt(uiCelIndex);

        if (m_pEngine->getOpt("snp-reference-output-file") != "")
        {
            if (loadGenotypeCallsFromFile())
            {
                for (int y = 0; (y < getProbeSets()->getCount()); y++)
                {
                    m_mxGenotypeCalls.set(getExperiment()->getIndex(), y, getProbeSets()->getAt(y)->getGenotypeCall());
                }
            }
            else
            {
                warn(1, "No external genotyping data is available for sample " + m_pvExperiments->at(uiCelIndex)->getExperimentName() + ".  This will produce an invalid snp-reference file.");
            }
        }

        m_pvProbes->quickSort(1); // ProbeID
        ptsv5 = pGroup5Input->openTsv("Intensities-" + ::getInt(uiCelIndex), affx::FILE5_OPEN);
        for (unsigned int uiProbeIndex = 0; (uiProbeIndex < uiProbeCount); uiProbeIndex++)
        {
            ptsv5->nextLine();
            ptsv5->get(0, 2, &iProbeID);
            ptsv5->get(0, 3, &fIntensity);
            objSearch.setProbeID(iProbeID);
            int iSearchIndex = m_pvProbes->binarySearch(objSearch, 1); // ProbeID
            if (iSearchIndex != -1)
            {
                m_pvProbes->getAt(iSearchIndex)->setIntensity(fIntensity);
            }
        }
        ptsv5->close();
        delete ptsv5;

        if (        m_pEngine->isOptDefined("local-gc-background-correction-reference-method")      &&
                    m_pEngine->getOpt("local-gc-background-correction-reference-method") != "none"  &&
                    m_pEngine->getOpt("local-gc-background-correction-reference-method") != "")
        {
            adjustIntensitiesPDNN();
        }

        adjustIntensitiesHighPassFilter();

        m_pvProbes->quickSort(3); // by ProbeSetIndex, Allele, ProbeID
        ptsv5 = pGroup5Input->openTsv("Intensities-" + ::getInt(uiCelIndex), affx::FILE5_REPLACE);
        ptsv5->defineColumn(0, 0, "ProbeSetIndex", affx::FILE5_DTYPE_INT);
        ptsv5->defineColumn(0, 1, "Allele", affx::FILE5_DTYPE_INT);
        ptsv5->defineColumn(0, 2, "ProbeID", affx::FILE5_DTYPE_INT);
        ptsv5->defineColumn(0, 3, "Intensity", affx::FILE5_DTYPE_FLOAT);
        for (unsigned int uiProbeIndex = 0; (uiProbeIndex < uiProbeCount); uiProbeIndex++)
        {
            CNProbe* p = m_pvProbes->getAt(uiProbeIndex);
            ptsv5->set_i(0, 0, p->getProbeSetIndex());
            ptsv5->set_i(0, 1, p->getAllele());
            ptsv5->set_i(0, 2, p->getProbeID());
            ptsv5->set_f(0, 3, p->getIntensity());
            ptsv5->writeLevel(0);
        }
        ptsv5->close();
        delete ptsv5;
    }
    pGroup5Input->close();
    delete pGroup5Input;
    file5Input.close();

    processIntensities();
    writeCopyNumberReference();

        if (m_pEngine->getOptBool("create-snp-reference"))
        {
                loadSnpReferenceFile(m_pEngine->getOpt("snp-reference-input-file"));
                writeSnpReference();
        }

    Verbose::out(1, "CNAnalysisMethodReference::runPart2(...) end");
}


void CNAnalysisMethodReference::calculateReferenceSketch()
{
    m_pvProbes->quickSort(1); // By ProbeID.

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
    unsigned int uiCelCount =  m_pEngine->getOptVector("cels").size();
    std::vector<float> vChipSketch(uiSketchSize);
    std::vector<double> vSNPReferenceSketch(uiSketchSize);
    std::vector<double> vCNReferenceSketch(uiSketchSize);
    std::vector<float> vSNPData;
    std::vector<float> vCNData;
    std::vector<float> vAllData;
    int iSNPDataCount = 0;
    int iCNDataCount = 0;
    int iOtherDataCount = 0;
    std::vector<char> vDataType;
    CNProbe objSearch;
    SketchQuantNormTran *qnorm = NULL;
    if (!((m_pEngine->getOptBool("doDualNormalization"))))
    {
        qnorm = new SketchQuantNormTran(uiSketchSize, false, false, false, 0, false);
    }
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        m_pobjExperiment = (*getExperiments())[uiCelIndex];
        AffxString strCelFileName = m_pEngine->getOptVector("cels")[uiCelIndex];
        affymetrix_fusion_io::FusionCELData cel;
        cel.SetFileName(strCelFileName.c_str());
        if (!cel.Read())
        {
            throw(Except("\nCan't read cel file: " + cel.GetFileName() + "\nMessage: " + StringUtils::ConvertWCSToMBS(cel.GetError())));
        }

        int numCells = cel.GetRows() * cel.GetCols();
        vDataType.resize(numCells);
        if (uiCelIndex == 0)
        {
            int iIndex = 0;
            for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
            {
                if ((iIndex < m_pvProbes->getCount()) && (m_pvProbes->getAt(iIndex)->getProbeID() == (iProbeIndex + 1)))
                {
                    if (getProbeSets()->getAt(m_pvProbes->getAt(iIndex)->getProbeSetIndex())->isUseForSketch())
                    {
                        if (getProbeSets()->getAt(m_pvProbes->getAt(iIndex)->getProbeSetIndex())->processAsSNPNormalize())
                        {
                            vDataType[iProbeIndex] = 1;
                            iSNPDataCount++;
                        }
                        else if (getProbeSets()->getAt(m_pvProbes->getAt(iIndex)->getProbeSetIndex())->processAsCNNormalize())
                        {
                            vDataType[iProbeIndex] = 2;
                            iCNDataCount++;
                        }
                    }
                    iIndex++;
                }
                else
                {
                    vDataType[iProbeIndex] = 3;
                    iOtherDataCount++;
                }
            }
        }


        if (!m_pEngine->getOptBool("doDualNormalization"))
        {
            vAllData.resize(numCells);
            if ((m_pEngine->isOptDefined("adapter-type-normalization")) && (m_pEngine->getOptBool("adapter-type-normalization")))
            {
                adapterTypeNormalization(cel);
                for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
                {
                    vAllData[iProbeIndex] = m_pobjAdapterTypeNormTran->transform(iProbeIndex, 0, cel.GetIntensity(iProbeIndex));
                }
            }
             else if ((m_pEngine->isOptDefined("signal-adjustment-covariates")) && 
                     (m_pEngine->getOpt("signal-adjustment-covariates") != "") &&
                     (m_pEngine->getOpt("signal-adjustment-covariates") != "none"))
            {
                adjustAlternateIntensitiesUsingCovariateSignalAdjustment(cel, vAllData);
            }
            else
            {
                for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
                {
                    vAllData[iProbeIndex] = cel.GetIntensity(iProbeIndex);
                }
            }
        }
        else
        {
            vSNPData.resize(iSNPDataCount);
            vCNData.resize(iCNDataCount);
            int iSNPIndex = 0;
            int iCNIndex = 0;

            if ((m_pEngine->isOptDefined("adapter-type-normalization")) && (m_pEngine->getOptBool("adapter-type-normalization")))
            {

              adapterTypeNormalization(cel);
                for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
                {
                    if (vDataType[iProbeIndex] == 1)
                    {
                        vSNPData[iSNPIndex] = m_pobjAdapterTypeNormTran->transform(iProbeIndex, 0, cel.GetIntensity(iProbeIndex));
                        iSNPIndex++;
                    }
                    else if (vDataType[iProbeIndex] == 2)
                    {
                        vCNData[iCNIndex] = m_pobjAdapterTypeNormTran->transform(iProbeIndex, 0, cel.GetIntensity(iProbeIndex));
                        iCNIndex++;
                    }
                }
            }
            else if ((m_pEngine->isOptDefined("signal-adjustment-covariates")) && 
                     (m_pEngine->getOpt("signal-adjustment-covariates") != "") &&
                     (m_pEngine->getOpt("signal-adjustment-covariates") != "none"))
            {
               vAllData.resize(numCells);
               adjustAlternateIntensitiesUsingCovariateSignalAdjustment(cel, vAllData);
                for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
                {
                    if (vDataType[iProbeIndex] == 1)
                    {
                        vSNPData[iSNPIndex] = vAllData[iProbeIndex];
                        iSNPIndex++;
                    }
                    else if (vDataType[iProbeIndex] == 2)
                    {
                        vCNData[iCNIndex] = vAllData[iProbeIndex];
                        iCNIndex++;
                    }
                }
            }
            else
            {

                for (int iProbeIndex = 0; iProbeIndex < numCells; iProbeIndex++)
                {
                    if (vDataType[iProbeIndex] == 1)
                    {
                        vSNPData[iSNPIndex] = cel.GetIntensity(iProbeIndex);
                        iSNPIndex++;
                    }
                    else if (vDataType[iProbeIndex] == 2)
                    {
                        vCNData[iCNIndex] = cel.GetIntensity(iProbeIndex);
                        iCNIndex++;
                    }
                }
            }
        }
        cel.Close();
        if (!(m_pEngine->getOptBool("doDualNormalization")))
        {
            qnorm->newChip(vAllData);
        }
        else
        {
            extractSketch(vSNPData, vChipSketch);
            for (int iIndex = 0; (iIndex < uiSketchSize); iIndex++)
            {
                vSNPReferenceSketch[iIndex] += vChipSketch[iIndex];
            }
            extractSketch(vCNData, vChipSketch);
            for (int iIndex = 0; (iIndex < uiSketchSize); iIndex++)
            {
                vCNReferenceSketch[iIndex] += vChipSketch[iIndex];
            }
        }
    }
    if (!(m_pEngine->getOptBool("doDualNormalization")))
    {
        qnorm->finishedChips();
        AffxString strFileName = m_pEngine->getOpt("reference-file");
        affx::File5_File file5;
        affx::File5_Group* group5;
        file5.open(strFileName, affx::FILE5_OPEN);
        group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        qnorm->saveTargetSketchToGroup_a5(group5);
        group5->close();
        delete group5;
        file5.close();
        delete qnorm;
    }
    else
    {
        for (int iIndex = 0; (iIndex < uiSketchSize); iIndex++)
        {
            vSNPReferenceSketch[iIndex] /= (double)uiCelCount;
            vCNReferenceSketch[iIndex] /= (double)uiCelCount;
        }

        writeReferenceSketch("SketchSNP", vSNPReferenceSketch);
        writeReferenceSketch("SketchCN", vCNReferenceSketch);
    }
}

bool CNAnalysisMethodReference::writeReferenceSketch(const std::string& strName, std::vector<double>& vReferenceSketch)
{
    unsigned int uiSketchSize = vReferenceSketch.size();
    AffxString strFileName = m_pEngine->getOpt("reference-file");
    affx::File5_File file5;
    affx::File5_Group* group5;
    affx::File5_Tsv* tsv5 = NULL;

    bool bSuccessful = false;

    if(m_pEngine->getOptBool("cyto2"))
    {
        try
        {
            file5.open(strFileName, affx::FILE5_OPEN);
            group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            tsv5 = group5->openTsv(strName, affx::FILE5_REPLACE);
            tsv5->defineColumn(0, 0, "Sketch", affx::FILE5_DTYPE_DOUBLE);
            for (int iIndex = 0; (iIndex < uiSketchSize); iIndex++)
            {
                tsv5->set_d(0, 0, vReferenceSketch[iIndex]);
                tsv5->writeLevel(0);
            }
            tsv5->close();
            delete tsv5;
            group5->close();
            delete group5;
            file5.close();
            bSuccessful = true;
        } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    }
    else
    {
        try
        {
            file5.open(strFileName, affx::FILE5_OPEN);
            group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            tsv5 = group5->openTsv(strName, affx::FILE5_REPLACE);
            tsv5->defineColumn(0, 0, "Sketch", affx::FILE5_DTYPE_DOUBLE);
            for (int iIndex = 0; (iIndex < uiSketchSize); iIndex++)
            {
                tsv5->set_d(0, 0, vReferenceSketch[iIndex]);
                tsv5->writeLevel(0);
            }
            tsv5->close();
            delete tsv5;
            group5->close();
            delete group5;
            file5.close();
            bSuccessful = true;
        } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    }

    return bSuccessful;
}

void CNAnalysisMethodReference::newProbes()
{
    int* piProbeSetIndexes = NULL;
    char* pcAlleles = NULL;
    ChipLayout* pLayout = new ChipLayout;
    CNAnalysisMethodReference::loadLayout(pLayout, m_pEngine);

    // Initialize the adapter type normalization
    if ((m_pEngine->isOptDefined("adapter-type-normalization")) && (m_pEngine->getOptBool("adapter-type-normalization")))
    {
        m_pobjAdapterTypeNormTran = new AdapterTypeNormTran(0.8287, 0.7960, 1.4218, 0.9954, 1.63920);
        m_pobjAdapterTypeNormTran->setup(m_pEngine->getOpt("annotation-file"), pLayout, 1);
    }

    unsigned int uiProbeCount = pLayout->getProbeCount();
    piProbeSetIndexes = new int[uiProbeCount];
    pcAlleles = new char[uiProbeCount];
    for (int i = 0; (i < uiProbeCount); i++)
    {
        piProbeSetIndexes[i] = -1;
        pcAlleles[i] = -1;
    }
    // Write out antigenomic probes
    AffxString strFileName = m_pEngine->getOpt("reference-file");
    affx::File5_File file5;
    affx::File5_Group* group5;
    affx::File5_Tsv* tsv5 = NULL;

    bool bSuccessful = false;
    try
    {
        file5.open(strFileName, affx::FILE5_OPEN);
        if (!((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2"))))
        {
            group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        }
        else
        {
            group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
        }
        tsv5 = group5->openTsv("AntigenomicProbes", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "ProbeID", affx::FILE5_DTYPE_INT);
        unsigned int uiCount = pLayout->getProbeSetCount();
        for (unsigned int uiIndex = 0; (uiIndex < uiCount); uiIndex++)
        {
            ProbeListPacked List = pLayout->getProbeListAtIndex(uiIndex);
            ProbeSet* ps = ProbeListFactory::asProbeSet(List);
            AffxString strProbeSetName = ps->name;
            strProbeSetName = strProbeSetName.toLowerCase();
            if ((strProbeSetName.startsWith("anti")) || (strProbeSetName.startsWith("randomgc")))
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
                                int iProbeID = (ps->atoms[iAtomCount]->probes[iProbeIndex]->id + 1);
                                tsv5->set_i(0, 0, iProbeID);
                                tsv5->writeLevel(0);
                            }
                        }
                        iAtomCount++;
                    }
                }
            }
            delete ps;
        }
        tsv5->close();
        delete tsv5;
        group5->close();
        delete group5;
        file5.close();
        bSuccessful = true;
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}

    unsigned int uiCount = 0;
    for (int iProbeSetIndex = 0; (iProbeSetIndex < getProbeSets()->getCount()); iProbeSetIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iProbeSetIndex);
        if (pLayout->containsProbeSet(pobjProbeSet->getProbeSetName()))
        {
            ProbeListPacked pList = pLayout->getProbeListByName(pobjProbeSet->getProbeSetName());
            if (!pList.isNull())
            {
                ProbeSet* ps = ProbeListFactory::asProbeSet(pList);
                AffxString strProbeSetName = ps->name;
                strProbeSetName = strProbeSetName.toLowerCase();
                if ((ps->psType == ProbeSet::GenoType) || (ps->psType == ProbeSet::Copynumber) || (strProbeSetName.startsWith("anti")) || (strProbeSetName.startsWith("randomgc")))
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
                                    piProbeSetIndexes[uiIndex] = iProbeSetIndex;
                                    pcAlleles[uiIndex] = 0;
                                    if (ps->psType == ProbeSet::GenoType)
                                    {
                                        if ((iGroupIndex % 2) == 0)
                                        {
                                            pcAlleles[uiIndex] = 'A';
                                        }
                                        else
                                        {
                                            pcAlleles[uiIndex] = 'B';
                                        }
                                    }
                                    uiCount++;
                                }
                            }
                            iAtomCount++;
                        }
                    }
                }
                delete ps;
            }
        }
    }
    memory("CNAnalysisMethodReference::newProbes(...)");
    delete pLayout;

    m_pProbes = new CNProbe[uiCount];
    int iProbeIndex = 0;
    for (int i = 0; (i < uiProbeCount); i++)
    {
        if (piProbeSetIndexes[i] != -1)
        {
            m_pProbes[iProbeIndex].setProbeID(i + 1);
            m_pProbes[iProbeIndex].setProbeSetIndex(piProbeSetIndexes[i]);
            m_pProbes[iProbeIndex].setAllele(pcAlleles[i]);
            iProbeIndex++;
        }
    }
    delete[] piProbeSetIndexes;
    delete[] pcAlleles;

    m_pvProbes = new CNProbeArray;
    m_pvProbes->reserve(uiCount);
    for (int iProbeIndex = 0; (iProbeIndex < uiCount); iProbeIndex++)
    {
        m_pvProbes->add(&m_pProbes[iProbeIndex]);
    }
}

void CNAnalysisMethodReference::processIntensities()
{
    AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
    const bool isCreateSnpRef = m_pEngine->getOptBool("create-snp-reference");
    affx::File5_File file5;
    affx::File5_Group* group5;
    affx::File5_Tsv* tsv5ProbeEffects = NULL;
    affx::File5_Tsv* tsv5SnpPosteriors = NULL;

    unsigned int uiCelCount =  m_pEngine->getOptVector("cels").size();
    unsigned int uiProbeCount = m_pvProbes->getCount();
    int iProbeSetIndex = 0;
    char cAllele = 0;
    unsigned int uiProbeID = 0;
    int iPrevProbeSetIndex = -1;
    char cPrevAllele = -1;
    float* pfIntensities = new float[uiCelCount];
    int iIndex = 0;

    std::vector<int> vProbeIDs(m_pvProbes->getMaximumNumberProbesPerProbeSet());
    m_pvProbes->quickSort(1); // ProbeID
    double* pdAAlleleSignals = new double[uiCelCount];
    double* pdBAlleleSignals = new double[uiCelCount];
    vector<affx::GType> tcall(uiCelCount);
    vector<float> vMedianIntensities(100);

    file5.open(strReferenceFileName, affx::FILE5_OPEN);
    if (!((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2"))))
    {
        group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
        tsv5ProbeEffects = group5->openTsv("CN5.plier.feature-response", affx::FILE5_REPLACE);
        tsv5SnpPosteriors = group5->openTsv("CN5.snp-posteriors", affx::FILE5_OPEN|affx::FILE5_CREATE);
    }
    else
    {
        group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
        tsv5ProbeEffects = group5->openTsv("ProbeEffects", affx::FILE5_REPLACE);
        tsv5SnpPosteriors = group5->openTsv("SNPPosteriors", affx::FILE5_OPEN|affx::FILE5_CREATE);
    }

    tsv5ProbeEffects->defineColumn(0, 0, "probe_id", affx::FILE5_DTYPE_INT);
    tsv5ProbeEffects->defineColumn(0, 1, "feature_response", affx::FILE5_DTYPE_DOUBLE);
    tsv5ProbeEffects->defineColumn(0, 2, "ProbeSetIndex", affx::FILE5_DTYPE_INT);
    tsv5ProbeEffects->defineColumn(0, 3, "AlleleCode", affx::FILE5_DTYPE_INT);
    tsv5ProbeEffects->defineColumn(0, 4, "MedianNormalizedIntensity", affx::FILE5_DTYPE_FLOAT);
    if(m_pEngine->getOptBool("pdnn-reference-values-exist"))
    {
        tsv5ProbeEffects->defineColumn(0, 5, "PredictedIntensity", affx::FILE5_DTYPE_FLOAT);
    }

    tsv5SnpPosteriors->defineColumn(0, 0, "id", affx::FILE5_DTYPE_STRING, 30);
    tsv5SnpPosteriors->defineColumn(0, 1, "Cluster_AA_MeanStrength", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 2, "Cluster_AA_Mean", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 3, "Cluster_AA_Variance", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 4, "Cluster_AA_VarianceStrength", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 5, "Cluster_AB_MeanStrength", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 6, "Cluster_AB_Mean", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 7, "Cluster_AB_Variance", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 8, "Cluster_AB_VarianceStrength", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 9, "Cluster_BB_MeanStrength", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 10, "Cluster_BB_Mean", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 11, "Cluster_BB_Variance", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 12, "Cluster_BB_VarianceStrength", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 13, "Cross-Term_PAH", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 14, "Cross-Term_PAB", affx::FILE5_DTYPE_DOUBLE);
    tsv5SnpPosteriors->defineColumn(0, 15, "Cross-Term_PHB", affx::FILE5_DTYPE_DOUBLE);


    affx::File5_File file5Input;
    affx::File5_Group* pGroup5Input;
    affx::File5_Tsv** pptsv5 = new affx::File5_Tsv*[uiCelCount];
    affx::File5_Tsv** pptsv5Signals = new affx::File5_Tsv*[uiCelCount];

    // Determine Feature Effects, and SNP Posteriors.
    file5Input.open(m_strTempFileName, affx::FILE5_OPEN);
    pGroup5Input = file5Input.openGroup("Intensities", affx::FILE5_OPEN);
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5Signals[uiCelIndex] = pGroup5Input->openTsv("Signals-" + ToStr(uiCelIndex), affx::FILE5_OPEN|affx::FILE5_CREATE);
        pptsv5Signals[uiCelIndex]->defineColumn(0, 0, getExperiments()->getAt(uiCelIndex)->getExperimentName() + "_A", affx::FILE5_DTYPE_DOUBLE);
        pptsv5Signals[uiCelIndex]->defineColumn(0, 1, getExperiments()->getAt(uiCelIndex)->getExperimentName() + "_B", affx::FILE5_DTYPE_DOUBLE);
    }

    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5[uiCelIndex] = pGroup5Input->openTsv("Intensities-" + ::getInt(uiCelIndex), affx::FILE5_OPEN);
    }
    int i = 0;
    AffxMultiDimensionalArray<float> vIntensities(uiCelCount);
    for (unsigned int uiProbeIndex = 0; (uiProbeIndex < uiProbeCount); uiProbeIndex++)
    {
        for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
        {
            pptsv5[uiCelIndex]->nextLine();
            pptsv5[uiCelIndex]->get(0, 0, &iProbeSetIndex);
            pptsv5[uiCelIndex]->get(0, 1, &i); cAllele = (char)i;
            pptsv5[uiCelIndex]->get(0, 2, &i); uiProbeID = (unsigned int)i;
            pptsv5[uiCelIndex]->get(0, 3, &pfIntensities[uiCelIndex]);
            tcall[uiCelIndex] = affx::NN;
        }

        if ((iPrevProbeSetIndex == -1) || ((iPrevProbeSetIndex != iProbeSetIndex) || (cPrevAllele != cAllele)))
        {
            if (iPrevProbeSetIndex != -1)
            {
                callPlier(      iPrevProbeSetIndex,
                                cPrevAllele,
                                vProbeIDs,
                                iIndex,
                                tsv5ProbeEffects,
                                vMedianIntensities);
                CNProbeSet* p = getProbeSets()->getAt(iPrevProbeSetIndex);
                AffxString strProbeSetName=p->getProbeSetName();
                setSignals(cPrevAllele, pdAAlleleSignals, pdBAlleleSignals, pptsv5Signals,strProbeSetName);
                if (cPrevAllele == 'B')
                {
                    callBrlmmp(iPrevProbeSetIndex, iIndex, pdAAlleleSignals, pdBAlleleSignals, tcall, tsv5SnpPosteriors);
                }
                if ((cPrevAllele == 0) || (cPrevAllele == 'B'))
                {
                    calculateMedians(iPrevProbeSetIndex, pdAAlleleSignals, pdBAlleleSignals, tcall);
                    if (isCreateSnpRef)
                    {
                        CNProbeSet* p = getProbeSets()->getAt(iPrevProbeSetIndex);
                        if (p->getChromosome() != m_iYChromosome && p->processAsSNP())
                        {
                            calculateSnpReference(iPrevProbeSetIndex, pdAAlleleSignals, pdBAlleleSignals, tcall);
                        }
                    }
                }
            }
            iPrevProbeSetIndex = iProbeSetIndex;
            cPrevAllele = cAllele;
            iIndex = 0;
        }
        vProbeIDs[iIndex] = uiProbeID;
        m_pdProbeEffects[iIndex] = 0;
        for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
        {
            m_pdChipEffects[uiCelIndex] = 0;
            m_ppdPM[uiCelIndex][iIndex] = pfIntensities[uiCelIndex];
            m_ppdMM[uiCelIndex][iIndex] = 0;
            m_ppdResiduals[uiCelIndex][iIndex] = 0;
            vIntensities.set(uiCelIndex, pfIntensities[uiCelIndex]);
        }
        if (iIndex >= vMedianIntensities.size()) {vMedianIntensities.push_back(0);}
        vMedianIntensities[iIndex] = vIntensities.median();
        iIndex++;
    }
    if (iPrevProbeSetIndex != -1)
    {
        callPlier(      iPrevProbeSetIndex,
                        cPrevAllele,
                        vProbeIDs,
                        iIndex,
                        tsv5ProbeEffects,
                        vMedianIntensities);

        CNProbeSet* p = getProbeSets()->getAt(iPrevProbeSetIndex);
        AffxString strProbeSetName=p->getProbeSetName();
        setSignals(cPrevAllele, pdAAlleleSignals, pdBAlleleSignals, pptsv5Signals, strProbeSetName);
        if (cPrevAllele == 'B')
        {
            callBrlmmp(iPrevProbeSetIndex, iIndex, pdAAlleleSignals, pdBAlleleSignals, tcall, tsv5SnpPosteriors);
        }
        if ((cPrevAllele == 0) || (cPrevAllele == 'B'))
        {
            calculateMedians(iPrevProbeSetIndex, pdAAlleleSignals, pdBAlleleSignals, tcall);
            if (m_pEngine->getOptBool("create-snp-reference"))
            {
                CNProbeSet* p = getProbeSets()->getAt(iPrevProbeSetIndex);
                if ( p->getChromosome() != m_iYChromosome )
                {
                    calculateSnpReference(iPrevProbeSetIndex, pdAAlleleSignals, pdBAlleleSignals, tcall);
                }
            }
        }
    }
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5Signals[uiCelIndex]->close();
        delete pptsv5Signals[uiCelIndex];
    }

    tsv5SnpPosteriors->close();
    delete tsv5SnpPosteriors; tsv5SnpPosteriors = NULL;

    tsv5ProbeEffects->close();
    delete tsv5ProbeEffects; tsv5ProbeEffects = NULL;

    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5[uiCelIndex]->close();
        delete pptsv5[uiCelIndex];
        group5->deleteTsv("Intensities-" + ::getInt(uiCelIndex));
    }

    // Determine hasXX and hasY.
    Verbose::out(1, "*");
    CNAnalysisMethodChipstream chipstream;
    chipstream.setEngine(m_pEngine);
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        // Load up the signals.
        double d = 0;
        affx::File5_Tsv* tsv5Signals = pGroup5Input->openTsv("Signals-" + ToStr(uiCelIndex), affx::FILE5_OPEN_RO);
        for (int iProbeSetIndex = 0; (iProbeSetIndex < getProbeSets()->getCount()); iProbeSetIndex++)
        {
            CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
            tsv5Signals->nextLine();
            tsv5Signals->get(0, 0, &d);
            p->setAAlleleSignal(d);
            tsv5Signals->get(0, 1, &d);
            p->setBAlleleSignal(d);
        }
        tsv5Signals->close();
        delete tsv5Signals; tsv5Signals = NULL;

        // Calculate XY.
        chipstream.setup(*m_pvExperiments->getAt(uiCelIndex), *getProbeSets());
        chipstream.calculateXY();
    }
    Verbose::out(1, "*");

    // Now load up the signals by Experiment and calculate the XX and Y medians.
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5Signals[uiCelIndex] = pGroup5Input->openTsv("Signals-" + ToStr(uiCelIndex), affx::FILE5_OPEN_RO);
    }
    for (int iProbeSetIndex = 0; (iProbeSetIndex < getProbeSets()->getCount()); iProbeSetIndex++)
    {
        double d = 0;
        for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
        {
            pptsv5Signals[uiCelIndex]->nextLine();
            pptsv5Signals[uiCelIndex]->get(0, 0, &d);
            pdAAlleleSignals[uiCelIndex] = d;
            pptsv5Signals[uiCelIndex]->get(0, 1, &d);
            pdBAlleleSignals[uiCelIndex] = d;
        }
        calculateXYMedians(iProbeSetIndex, pdAAlleleSignals, pdBAlleleSignals);
    }
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        pptsv5Signals[uiCelIndex]->close();
        delete pptsv5Signals[uiCelIndex]; pptsv5Signals[uiCelIndex] = NULL;
    }

    pGroup5Input->close();
    delete pGroup5Input;
    file5Input.close();

    group5->close();
    delete group5;
    file5.close();

    delete[] pdAAlleleSignals;
    delete[] pdBAlleleSignals;
    delete[] pfIntensities;
    delete[] pptsv5;
    delete[] pptsv5Signals;
}

void CNAnalysisMethodReference::setSignals(    char cAllele,
                                                double* pdAAlleleSignals,
                                                double* pdBAlleleSignals,
                                                affx::File5_Tsv** pptsv5Signals,
                                                AffxString strProbeSetName )
{
    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
    const bool isCyto2 = m_pEngine->getOptBool("cyto2");

    unsigned int uiCelCount =  m_pEngine->getOptVector("cels").size();
    if ((cAllele == 0) || (cAllele == 'A'))
    {
        for (int i = 0; (i < uiCelCount); i++)
        {
            pdAAlleleSignals[i] = m_pdChipEffects[i];
            pdBAlleleSignals[i] = 0;
            if (cAllele == 0)
            {
                if (isCyto2 || isCytoScanHD)
                {
                    pptsv5Signals[i]->set_d(0, 0, pdAAlleleSignals[i]);
                }
                else // Must be SNP6
                {
                    pptsv5Signals[i]->set_d(0, 0, (float)pdAAlleleSignals[i]);
                }
                pptsv5Signals[i]->set_d(0, 1, 0);
                pptsv5Signals[i]->writeLevel(0);
            }
            else
            {
                if (isCyto2 || isCytoScanHD)
                {
                    pptsv5Signals[i]->set_d(0, 0, pdAAlleleSignals[i]);
                }
                else // Must be SNP6
                {
                    pptsv5Signals[i]->set_d(0, 0, (float)pdAAlleleSignals[i]);
                }
            }
        }
    }
    else if (cAllele == 'B')
    {
        for (int i = 0; (i < uiCelCount); i++)
        {
            pdBAlleleSignals[i] = m_pdChipEffects[i];
                if (isCyto2 || isCytoScanHD)
            {
                pptsv5Signals[i]->set_d(0, 1, pdBAlleleSignals[i]);
            }
            else // Must be SNP6
            {
                pptsv5Signals[i]->set_d(0, 1, (float)pdBAlleleSignals[i]);
            }
            pptsv5Signals[i]->writeLevel(0);
        }
    }
}

void CNAnalysisMethodReference::callPlier(    int iProbeSetIndex,
                                                char cAllele,
                                                std::vector<int>& vProbeIDs,
                                                int iProbeCount,
                                                affx::File5_Tsv* tsv5,
                                                std::vector<float>& vMedianIntensities)
{
    const bool pdnnReferenceValuesExist = m_pEngine->getOptBool("pdnn-reference-values-exist");

    long errorCode = 0;
    m_pPlier->setNumFeature(iProbeCount);
    m_pPlier->run(&errorCode);
    if (errorCode != 0) {Err::errAbort("Problem running plier. Error code: " + ToStr(errorCode));}
    CNProbe objSearch;
    if (tsv5 != NULL)
    {
        for (int i = 0; (i < iProbeCount); i++)
        {
            float fPredictedIntensity = 0;
            float fMedianIntensity = 0;
            objSearch.setProbeID(vProbeIDs[i]);
            int iSearchIndex = m_pvProbes->binarySearch(objSearch, 1); // ProbeID
            if (iSearchIndex != -1)
            {
                if(pdnnReferenceValuesExist)
                {
                    fPredictedIntensity = m_pvProbes->getAt(iSearchIndex)->getPredictedIntensity();
                }
                fMedianIntensity = m_pvProbes->getAt(iSearchIndex)->getMedianIntensity();
            }

           // 09/14/2009 Walt Short - This is the fix for the 'Known' bug in v2.0 that was released.
            // It will insure that the adjusted PDNN intensity number will match between the cn reference build run
            // and the single sample runs.
            //tsv5->set_f(0, 4, vMedianIntensities[i]);
            tsv5->set_i(0, 0, vProbeIDs[i]);
            tsv5->set_d(0, 1, m_pdProbeEffects[i]);
            tsv5->set_i(0, 2, iProbeSetIndex);
            tsv5->set_i(0, 3, (int)cAllele);
            tsv5->set_f(0, 4, fMedianIntensity);
            if(pdnnReferenceValuesExist)
            {
                tsv5->set_f(0, 5, fPredictedIntensity);
            }
            tsv5->writeLevel(0);
        }
    }
}

void CNAnalysisMethodReference::callBrlmmp(    int iProbeSetIndex,
                                                int iProbeCount,
                                                double* pdAAlleleSignals,
                                                double* pdBAlleleSignals,
                                                std::vector<affx::GType>& tcalls,
                                                affx::File5_Tsv* tsv5)
{
    if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2"))) {return;}
    unsigned int uiCelCount =  m_pEngine->getOptVector("cels").size();
    CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
    snp_param sp;
    initializeSnpParameters(sp, 1);
    setBrlmmpParameters(sp);
    vector<affx::GType> tcall(uiCelCount);
    vector<double> tconf(uiCelCount);
    vector<double> tx(uiCelCount);
    vector<double> ty(uiCelCount);
    vector<int> GenoHint(uiCelCount);
    vector<double> SubsetInbred(uiCelCount);
    snp_param tsp;
    for (int iCopyNumber = 2; (iCopyNumber > 0); iCopyNumber--)
    {
        int iProcessCount = 0;
        tsp.copy(sp); // duplicate current parameters to keep commands
        for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
        {
            tsp.copynumber = 2; // what copy for this snp so know right model centers
            if (m_pvExperiments->at(uiCelIndex)->getRawIntensityRatioGenderAsInt() == affx::Male)
            {
                if ((p->getChromosome() == m_iXChromosome) && (!p->isPseudoAutosomalRegion())) {tsp.copynumber = 1;}
                else if (p->getChromosome() == m_iYChromosome) {tsp.copynumber = 1;}
            }
            else
            {
                if (p->getChromosome() == m_iYChromosome) {tsp.copynumber = 0;}
            }
            if (tsp.copynumber == iCopyNumber) {iProcessCount++;}
        }
        if (iProcessCount > 0)
        {
            tcall.resize(iProcessCount);
            tconf.resize(iProcessCount);
            tx.resize(iProcessCount);
            ty.resize(iProcessCount);
            GenoHint.resize(iProcessCount);
            SubsetInbred.resize(iProcessCount);
            int iProcessIndex = 0;
            for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
            {
                tsp.copynumber = 2; // what copy for this snp so know right model centers
                if (m_pvExperiments->at(uiCelIndex)->getRawIntensityRatioGenderAsInt() == affx::Male)
                {
                    if ((p->getChromosome() == m_iXChromosome) && (!p->isPseudoAutosomalRegion())) {tsp.copynumber = 1;}
                    else if (p->getChromosome() == m_iYChromosome) {tsp.copynumber = 1;}
                }
                else
                {
                    if (p->getChromosome() == m_iYChromosome) {tsp.copynumber = 0;}
                }
                if (iCopyNumber == tsp.copynumber)
                {
                    tcall[iProcessIndex] = affx::NN;
                    tconf[iProcessIndex] = 0;
                    tx[iProcessIndex] = pdAAlleleSignals[uiCelIndex];
                    ty[iProcessIndex] = pdBAlleleSignals[uiCelIndex];
                    GenoHint[iProcessIndex] = affx::NN;
                    SubsetInbred[iProcessIndex] = 0; // null value
                    iProcessIndex++;
                }
            }
            tsp.copynumber = iCopyNumber;
            GenoUtility_transformData(tx, ty, MvA, 1);
            bayes_label(tcall,tconf,tx,ty,GenoHint,SubsetInbred,tsp);
            iProcessIndex = 0;
            for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
            {
                tsp.copynumber = 2; // what copy for this snp so know right model centers
                if (m_pvExperiments->at(uiCelIndex)->getRawIntensityRatioGenderAsInt() == affx::Male)
                {
                    if ((p->getChromosome() == m_iXChromosome) && (!p->isPseudoAutosomalRegion())) {tsp.copynumber = 1;}
                    else if (p->getChromosome() == m_iYChromosome) {tsp.copynumber = 1;}
                }
                else
                {
                    if (p->getChromosome() == m_iYChromosome) {tsp.copynumber = 0;}
                }
                if (iCopyNumber == tsp.copynumber)
                {
                    tcalls[uiCelIndex] = tcall[iProcessIndex];
                    iProcessIndex++;
                }
            }
            if (tsv5 != NULL)
            {
                if (iCopyNumber < 2)
                {
                    tsv5->set_string(0, 0, p->getProbeSetName() + ":" + ::getInt(iCopyNumber));
                }
                else
                {
                    tsv5->set_string(0, 0, p->getProbeSetName());
                }
                tsv5->set_d(0, 1, tsp.posterior.aa.k);
                tsv5->set_d(0, 2, tsp.posterior.aa.m);
                tsv5->set_d(0, 3, tsp.posterior.aa.ss);
                tsv5->set_d(0, 4, tsp.posterior.aa.v);
                tsv5->set_d(0, 5, tsp.posterior.ab.k);
                tsv5->set_d(0, 6, tsp.posterior.ab.m);
                tsv5->set_d(0, 7, tsp.posterior.ab.ss);
                tsv5->set_d(0, 8, tsp.posterior.ab.v);
                tsv5->set_d(0, 9, tsp.posterior.bb.k);
                tsv5->set_d(0, 10, tsp.posterior.bb.m);
                tsv5->set_d(0, 11, tsp.posterior.bb.ss);
                tsv5->set_d(0, 12, tsp.posterior.bb.v);
                tsv5->set_d(0, 13, tsp.posterior.xah);
                tsv5->set_d(0, 14, tsp.posterior.xab);
                tsv5->set_d(0, 15, tsp.posterior.xhb);
                tsv5->writeLevel(0);
            }
        }
    }
}

/*
float CNAnalysisMethodReference::computeSCAR(CNProbeSet* pobjProbeSet)
{
    float fSCAR=0.0;
    fSCAR = (
        2.0*( log2(pobjProbeSet->getAAlleleSignal()/pobjProbeSet->getBAlleleSignal()) - pobjProbeSet->getMuAB() )
        )
        /
        (
        fabs( pobjProbeSet->getMuAA() - pobjProbeSet->getMuAB() )
        +
        fabs( pobjProbeSet->getMuBB() - pobjProbeSet->getMuAB() )
        );
    if(fSCAR < -4.0) fSCAR = -4.0;
    if(fSCAR > 4.0) fSCAR = 4.0;

    return fSCAR;
}
*/

void CNAnalysisMethodReference::calculateSnpReference(    int iProbeSetIndex,
                            double* pdAAlleleSignals,
                            double* pdBAlleleSignals,
                            std::vector<affx::GType>& vCalls)
{
    CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
    unsigned int uiCelCount =  m_pEngine->getOptVector("cels").size();
    AffxMultiDimensionalArray<float> vAASignals(uiCelCount);
    AffxMultiDimensionalArray<float> vABSignals(uiCelCount);
    AffxMultiDimensionalArray<float> vBBSignals(uiCelCount);
    int iAAIndex = 0;
    int iABIndex = 0;
    int iBBIndex = 0;
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        if ( (p->getChromosome() == m_iXChromosome)  && (m_pvExperiments->at(uiCelIndex)->getRawIntensityRatioGenderAsInt() == affx::Male) ) continue;

        float fASignal = pdAAlleleSignals[uiCelIndex];
        float fBSignal = pdBAlleleSignals[uiCelIndex];
        float flogRatio = log2(fASignal/fBSignal);

        char cGenotypeCall = -2;
        if (m_pEngine->getOpt("snp-reference-output-file") != "")
        {
            cGenotypeCall = m_mxGenotypeCalls.get(uiCelIndex, iProbeSetIndex);
        }
        if ((cGenotypeCall != affx::AA) && (cGenotypeCall != affx::AB) && (cGenotypeCall != affx::BB))
                {
                        warn(3, "A no-call genotype value is being used for " + p->getProbeSetName() + " and sample " + m_pvExperiments->at(uiCelIndex)->getExperimentName() );
        }
        switch(cGenotypeCall)
        {
        case affx::AA:
            vAASignals.set(iAAIndex, flogRatio);
            iAAIndex++;
            break;
        case affx::AB:
            vABSignals.set(iABIndex, flogRatio);
            iABIndex++;
            break;
        case affx::BB:
            vBBSignals.set(iBBIndex, flogRatio);
            iBBIndex++;
            break;
        }
    }

        if(       (iAAIndex < 1 && iABIndex < 1)
                  ||
                  (iABIndex < 1 && iBBIndex < 1)
                  ||
                  (iAAIndex < 1 && iBBIndex < 1)
          )
        {
                p->setMuAA((float)::roundDouble(INVALID_DOUBLE, 10));
                p->setMuAB((float)::roundDouble(INVALID_DOUBLE, 10));
                p->setMuBB((float)::roundDouble(INVALID_DOUBLE, 10));

                warn(1, "Within the samples provided, missing genotypes were found for the probe set " + p->getProbeSetName() +". No snp-reference values will be calculated for the probe set." );

                return;
        }

        if(iAAIndex < 1)
        {
                double dAB = vABSignals.mean(iABIndex);
                double dBB = vBBSignals.mean(iBBIndex);
                double dAA = dAB - (dBB-dAB);
                p->setMuAA((float)::roundDouble(dAA, 10));
                p->setMuAB((float)::roundDouble(dAB, 10));
                p->setMuBB((float)::roundDouble(dBB, 10));
                return;
        }
        if(iABIndex < 1)
        {
                double dAA = vAASignals.mean(iAAIndex);
                double dBB = vBBSignals.mean(iBBIndex);
                double dAB = (dAA + dBB)/2.0;
                p->setMuAA((float)::roundDouble(dAA, 10));
                p->setMuAB((float)::roundDouble(dAB, 10));
                p->setMuBB((float)::roundDouble(dBB, 10));
                return;
        }
        if(iBBIndex < 1)
        {
                double dAA = vAASignals.mean(iAAIndex);
                double dAB = vABSignals.mean(iABIndex);
                double dBB = dAB - (dAA - dAB);
                p->setMuAA((float)::roundDouble(dAA, 10));
                p->setMuAB((float)::roundDouble(dAB, 10));
                p->setMuBB((float)::roundDouble(dBB, 10));
                return;
        }

        // If we get here we have representatives of all genotypes so all means can be calculated.
        p->setMuAA((float)::roundDouble(vAASignals.mean(iAAIndex), 10));
        p->setMuAB((float)::roundDouble(vABSignals.mean(iABIndex), 10));
        p->setMuBB((float)::roundDouble(vBBSignals.mean(iBBIndex), 10));
        return;
}

void CNAnalysisMethodReference::calculateMedians(int iProbeSetIndex, double* pdAAlleleSignals, double* pdBAlleleSignals, std::vector<affx::GType>& vCalls)
{
    CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
    unsigned int uiCelCount =  m_pEngine->getOptVector("cels").size();
    AffxMultiDimensionalArray<float> vSignals(uiCelCount);
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        float fASignal = pdAAlleleSignals[uiCelIndex];
        float fBSignal = pdBAlleleSignals[uiCelIndex];
        double dSignal = fASignal + fBSignal;

        vSignals.set(uiCelIndex, dSignal);
        char cGenotypeCall = -2;
        if (m_pEngine->getOpt("snp-reference-output-file") != "")
        {
            cGenotypeCall = m_mxGenotypeCalls.get(uiCelIndex, iProbeSetIndex);
        }
        if (cGenotypeCall == -2)
        {
            cGenotypeCall = vCalls[uiCelIndex];
            if ((m_pEngine->getOpt("snp-reference-output-file") != "") || (!m_pEngine->getOptBool("cyto2")))
            {
                m_mxGenotypeCalls.set(uiCelIndex, iProbeSetIndex, cGenotypeCall);
            }
        }
    }
    p->setMedianSignal((float)::roundDouble(vSignals.median(), 10));
}

void CNAnalysisMethodReference::calculateXYMedians(int iProbeSetIndex, double* pdAAlleleSignals, double* pdBAlleleSignals)
{
    CNProbeSet* p = getProbeSets()->getAt(iProbeSetIndex);
    unsigned int uiCelCount =  m_pEngine->getOptVector("cels").size();
    AffxMultiDimensionalArray<float> vXXSignals(uiCelCount);
    AffxMultiDimensionalArray<float> vYSignals(uiCelCount);
    AffxMultiDimensionalArray<float> vAASignals(uiCelCount);
    AffxMultiDimensionalArray<float> vABSignals(uiCelCount);
    AffxMultiDimensionalArray<float> vBBSignals(uiCelCount);
    int iXXIndex = 0;
    int iYIndex = 0;
    int iAAIndex = 0;
    int iABIndex = 0;
    int iBBIndex = 0;
    for (unsigned int uiCelIndex = 0; (uiCelIndex < uiCelCount); uiCelIndex++)
    {
        float fASignal = pdAAlleleSignals[uiCelIndex];
        float fBSignal = pdBAlleleSignals[uiCelIndex];
        double dSignal = fASignal + fBSignal;
        double dDifference = fASignal - fBSignal;

        if (m_pvExperiments->at(uiCelIndex)->hasXX())
        {
            vXXSignals.set(iXXIndex, dSignal);
            iXXIndex++;
        }
        if (m_pvExperiments->at(uiCelIndex)->hasY())
        {
            vYSignals.set(iYIndex, dSignal);
            iYIndex++;
        }
        if ((p->getChromosome() != m_iXChromosome) || (m_pvExperiments->at(uiCelIndex)->hasXX()))
        {
            if ((m_pEngine->getOpt("snp-reference-output-file") != "") || (!m_pEngine->getOptBool("cyto2")))
            {
                char cGenotypeCall = m_mxGenotypeCalls.get(uiCelIndex, iProbeSetIndex);
                switch(cGenotypeCall)
                {
                case affx::AA:
                    vAASignals.set(iAAIndex, dDifference);
                    iAAIndex++;
                    break;
                case affx::AB:
                    vABSignals.set(iABIndex, dDifference);
                    iABIndex++;
                    break;
                case affx::BB:
                    vBBSignals.set(iBBIndex, dDifference);
                    iBBIndex++;
                    break;
                }
            }
        }
    }
    p->setXXMedianSignal((float)::roundDouble(vXXSignals.median(iXXIndex), 10));
    p->setYMedianSignal((float)::roundDouble(vYSignals.median(iYIndex), 10));
    p->setAAMedianSignal((float)::roundDouble(vAASignals.median(iAAIndex), 10));
    p->setABMedianSignal((float)::roundDouble(vABSignals.median(iABIndex), 10));
    p->setBBMedianSignal((float)::roundDouble(vBBSignals.median(iBBIndex), 10));
}

bool CNAnalysisMethodReference::loadSnpReferenceFile(const AffxString& strFileName)
{
        Verbose::out(3, "CNAnalysisMethodReference::loadSnpReferenceFile");
        bool bSuccessful = false;
        if (affx::File5_File::isHdf5file(strFileName))
        {
                try
                {
                        affx::File5_File file5;
                        file5.open(strFileName, affx::FILE5_OPEN_RO);
                        affx::File5_Group* group5 = file5.openGroup("SNPReference", affx::FILE5_OPEN);
                        affx::File5_Tsv* tsv5;


                        tsv5 = group5->openTsv("EMSNParameters", affx::FILE5_OPEN);
                        double dInputDouble=0.0;
                        int iInputInteger=0;
                        while (tsv5->nextLine() == affx::FILE5_OK)
                        {
                                tsv5->get(0, 0, &dInputDouble);  m_fLambdaSpacing=(float)dInputDouble;
                                tsv5->get(0, 1, &dInputDouble);  m_fEMSNSpacing=(float)dInputDouble;
                                tsv5->get(0, 2, &dInputDouble);  m_fLambdaStartingValue=(float)dInputDouble;
                                tsv5->get(0, 3, &dInputDouble);  m_fEMSNStartingValue=(float)dInputDouble;
                                tsv5->get(0, 4, &iInputInteger); m_iNumberOfLambdaPoints=iInputInteger;
                                tsv5->get(0, 5, &iInputInteger); m_iNumberOfEMSNPoints=iInputInteger;
                        }
                        tsv5->close();
                        delete tsv5;

                        m_dEMSNArray.initialize(m_iNumberOfEMSNPoints, m_iNumberOfLambdaPoints);

                        tsv5 = group5->openTsv("EMSNArray", affx::FILE5_OPEN);

                        dInputDouble=0.0;
                        int iRowCount=0;
                        while (tsv5->nextLine() == affx::FILE5_OK)
                        {
                                for(int iIndex=0; iIndex < m_iNumberOfLambdaPoints; iIndex++)
                                {
                                        tsv5->get(0, iIndex, &dInputDouble);
                                        m_dEMSNArray.set(iRowCount, iIndex, dInputDouble);
                                }
                                iRowCount++;
                        }
                        tsv5->close();
                        delete tsv5;

                        if(group5->name_exists("VIRanges") && group5->name_exists("VIValues"))
                        {
                                tsv5 = group5->openTsv("VIRanges", affx::FILE5_OPEN);
                                if(tsv5)
                                {
                                        dInputDouble=0.0;
                                        while (tsv5->nextLine() == affx::FILE5_OK)
                                        {
                                                tsv5->get(0, 0, &dInputDouble);
                                                m_vInformationRanges.push_back(dInputDouble);
                                        }
                                }
                                tsv5->close();
                                delete tsv5;

                                tsv5 = group5->openTsv("VIValues", affx::FILE5_OPEN);
                                if(tsv5)
                                {
                                        dInputDouble=0.0;
                                        while (tsv5->nextLine() == affx::FILE5_OK)
                                        {
                                                tsv5->get(0, 0, &dInputDouble);
                                                m_vVarianceInflationValues.push_back(dInputDouble);
                                        }
                                tsv5->close();
                                delete tsv5;
                                }
                        }
                        else
                        {
                                m_vVarianceInflationValues.push_back(1.0);
                        }

                        tsv5 = group5->openTsv("FLDCutoff3peaks", affx::FILE5_OPEN);
                        float fThreePeakFLD_X;
                        float fThreePeakFLD_Y;
                        while (tsv5->nextLine() == affx::FILE5_OK) {
                            tsv5->get(0, 0, &fThreePeakFLD_X);  m_vThreePeakFLD_X.push_back(fThreePeakFLD_X);
                            tsv5->get(0, 1, &fThreePeakFLD_Y);  m_vThreePeakFLD_Y.push_back(fThreePeakFLD_Y);
                        }
                        tsv5->close();
                        delete tsv5;

                        tsv5 = group5->openTsv("FLDCutoff4peaks", affx::FILE5_OPEN);
                        float fFourPeakFLD_X;
                        float fFourPeakFLD_Y;
                        while (tsv5->nextLine() == affx::FILE5_OK) {
                            tsv5->get(0, 0, &fFourPeakFLD_X);  m_vFourPeakFLD_X.push_back(fFourPeakFLD_X);
                            tsv5->get(0, 1, &fFourPeakFLD_Y);  m_vFourPeakFLD_Y.push_back(fFourPeakFLD_Y);
                        }
                        tsv5->close();
                        delete tsv5;

                        tsv5 = group5->openTsv("ShrinkFactor3peaks", affx::FILE5_OPEN);
                        float fThreePeakShrink_X;
                        float fThreePeakShrink_Y;
                        while (tsv5->nextLine() == affx::FILE5_OK) {
                            tsv5->get(0, 0, &fThreePeakShrink_X);  m_vThreePeakShrink_X.push_back(fThreePeakShrink_X);
                            tsv5->get(0, 1, &fThreePeakShrink_Y);  m_vThreePeakShrink_Y.push_back(fThreePeakShrink_Y);
                        }
                        tsv5->close();
                        delete tsv5;

                        tsv5 = group5->openTsv("ShrinkFactor4peaks", affx::FILE5_OPEN);
                        float fFourPeakShrink_X;
                        float fFourPeakShrink_Y;
                        while (tsv5->nextLine() == affx::FILE5_OK) {
                            tsv5->get(0, 0, &fFourPeakShrink_X);  m_vFourPeakShrink_X.push_back(fFourPeakShrink_X);
                            tsv5->get(0, 1, &fFourPeakShrink_Y);  m_vFourPeakShrink_Y.push_back(fFourPeakShrink_Y);
                        }
                        tsv5->close();
                        delete tsv5;

                        group5->close();
                        delete group5;
                        file5.close();
                        bSuccessful = true;

                } catch(...) {throw(Except("Cannot open file: " + strFileName));}
        }
        return bSuccessful;
}


void CNAnalysisMethodReference::writeSnpReference()
{
    AffxString strFileName = m_pEngine->getOpt("snp-reference-output-file");
    try
    {
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_REPLACE);
        affx::File5_Group* group5 = file5.openGroup("SNPReference", affx::FILE5_REPLACE);
        affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "Parameter", affx::FILE5_DTYPE_STRING, 1024);
        tsv5->set_string(0, 0, "#%guid=" + affxutil::Guid::GenerateNewGuid()); tsv5->writeLevel(0);

                vector<string> optionNames;
                m_pEngine->getOptionNames(optionNames,1);
        for(int i=0; i< optionNames.size(); i++) {
            std::string name = optionNames[i];
                        std::vector<std::string> vals = m_pEngine->getOptVector(name,1);
            if (vals.size() > 1)
            {
                for(int j=0; j< vals.size(); j++)
                {
                                        std::string val = vals[j];
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-" + name + "=" + val); tsv5->writeLevel(0);
                }
            }
            else
            {
                                std::string val = m_pEngine->getOpt(name,1);
                tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-" + name + "=" + val); tsv5->writeLevel(0);
            }
        }
        // dump state.
                std::vector<std::string> stateNames;
                m_pEngine->getOptionNames(stateNames);
        for(int i=0; i< stateNames.size(); i++) {
            std::string name = stateNames[i];
                        std::vector<std::string> vals = m_pEngine->getOptVector(name);
            if (vals.size() > 1)
            {
                for(int j=0; j< vals.size(); j++)
                {
                                        std::string val = vals[j];
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-state-" + name + "=" + val); tsv5->writeLevel(0);
                }
            }
            else
            {
                                std::string val = m_pEngine->getOpt(name);
                tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-state-" + name + "=" + val); tsv5->writeLevel(0);
            }
        }

        tsv5->close();
        delete tsv5;

        int iMaxProbeSetNameLength = 0;
        for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
            iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
        }

        tsv5 = group5->openTsv("SNPReference", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "probeset_id", affx::FILE5_DTYPE_STRING, iMaxProbeSetNameLength);
        tsv5->defineColumn(0, 1, "muAB", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 2, "muAA", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 3, "muBB", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 4, "information", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 5, "useInEMAlgorithm", affx::FILE5_DTYPE_INT);
        for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
                        if(pobjProbeSet->processAsSNP()){

                    tsv5->set_string(0, 0, pobjProbeSet->getProbeSetName());
                    tsv5->set_d(0, 1, (double)pobjProbeSet->getMuAB());
                    tsv5->set_d(0, 2, pobjProbeSet->getMuAA());
                    tsv5->set_d(0, 3, pobjProbeSet->getMuBB());
                    tsv5->set_d(0, 4, pobjProbeSet->getInformation());
                    tsv5->set_i(0, 5, pobjProbeSet->getUseInEMAlgorithm());
                    tsv5->writeLevel(0);
                        }
        }
        tsv5->close();
        delete tsv5;


        tsv5 = group5->openTsv("EMSNParameters", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "LambdaSpacing", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 1, "EMSNSpacing", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 2, "LambdaStartingValue", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 3, "EMSNStartingValue", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 4, "NumberOfLambdaPoints", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 5, "NumberOfEMSNPoints", affx::FILE5_DTYPE_INT);
        tsv5->set_d(0, 0, m_fLambdaSpacing);
        tsv5->set_d(0, 1, m_fEMSNSpacing);
        tsv5->set_d(0, 2, m_fLambdaStartingValue);
        tsv5->set_d(0, 3, m_fEMSNStartingValue);
        tsv5->set_i(0, 4, m_iNumberOfLambdaPoints);
        tsv5->set_i(0, 5, m_iNumberOfEMSNPoints);
        tsv5->writeLevel(0);
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("EMSNArray", affx::FILE5_REPLACE);
                for(int iIndex=0; iIndex < m_iNumberOfLambdaPoints; iIndex++)
                {
                        tsv5->defineColumn(0, iIndex, "col_"+AffxString::intToString(iIndex+1,false), affx::FILE5_DTYPE_DOUBLE);
                }
                double dInsertionValue=0.0;
                for(int iEMSNIndex=0; iEMSNIndex < m_iNumberOfEMSNPoints; iEMSNIndex++)
                {
                        for(int iLambdaIndex=0; iLambdaIndex<m_iNumberOfLambdaPoints; iLambdaIndex++)
                        {
                                dInsertionValue = m_dEMSNArray.get(iEMSNIndex, iLambdaIndex);
                                tsv5->set_d(0, iLambdaIndex, dInsertionValue);
                        }
                tsv5->writeLevel(0);
                }
        tsv5->close();
        delete tsv5;



        tsv5 = group5->openTsv("VIRanges", affx::FILE5_REPLACE);
                tsv5->defineColumn(0, 0, "VIRanges", affx::FILE5_DTYPE_DOUBLE);
                for(int iRangeIndex=0; iRangeIndex<m_vInformationRanges.size(); iRangeIndex++)
                {
                        dInsertionValue = m_vInformationRanges[iRangeIndex];
                        tsv5->set_d(0, 0, dInsertionValue);
                        tsv5->writeLevel(0);
                }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("VIValues", affx::FILE5_REPLACE);
                tsv5->defineColumn(0, 0, "VIValues", affx::FILE5_DTYPE_DOUBLE);
                for(int iValueIndex=0; iValueIndex<m_vVarianceInflationValues.size(); iValueIndex++)
                {
                        dInsertionValue = m_vVarianceInflationValues[iValueIndex];
                        tsv5->set_d(0, 0, dInsertionValue);
                        tsv5->writeLevel(0);
                }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("FLDCutoff3peaks", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "SNPQC", affx::FILE5_DTYPE_FLOAT);
        tsv5->defineColumn(0, 1, "Cutoff", affx::FILE5_DTYPE_FLOAT);
        for (int i = 0; i < m_vThreePeakFLD_X.size(); i++) {
            tsv5->set_f(0, 0, m_vThreePeakFLD_X[i]);
            tsv5->set_f(0, 1, m_vThreePeakFLD_Y[i]);
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("FLDCutoff4peaks", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "SNPQC", affx::FILE5_DTYPE_FLOAT);
        tsv5->defineColumn(0, 1, "Cutoff", affx::FILE5_DTYPE_FLOAT);
        for (int i = 0; i < m_vFourPeakFLD_X.size(); i++) {
            tsv5->set_f(0, 0, m_vFourPeakFLD_X[i]);
            tsv5->set_f(0, 1, m_vFourPeakFLD_Y[i]);
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("ShrinkFactor3peaks", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "SNPQC", affx::FILE5_DTYPE_FLOAT);
        tsv5->defineColumn(0, 1, "Shrinkage", affx::FILE5_DTYPE_FLOAT);
        for (int i = 0; i < m_vThreePeakShrink_X.size(); i++) {
            tsv5->set_f(0, 0, m_vThreePeakShrink_X[i]);
            tsv5->set_f(0, 1, m_vThreePeakShrink_Y[i]);
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("ShrinkFactor4peaks", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "SNPQC", affx::FILE5_DTYPE_FLOAT);
        tsv5->defineColumn(0, 1, "Shrinkage", affx::FILE5_DTYPE_FLOAT);
        for (int i = 0; i < m_vFourPeakShrink_X.size(); i++) {
            tsv5->set_f(0, 0, m_vFourPeakShrink_X[i]);
            tsv5->set_f(0, 1, m_vFourPeakShrink_Y[i]);
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("FLDValues", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "probeset_id", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "FLD", affx::FILE5_DTYPE_FLOAT);
        for (int i = 0; i < getProbeSets()->size(); i++) {
            CNProbeSet* p = getProbeSets()->getAt(i);
            if (p->getValidFLDExists()) {
                tsv5->set_string(0, 0, p->getProbeSetName());
                tsv5->set_f(0, 1, p->getFLD());
                tsv5->writeLevel(0);
            }
        }
        tsv5->close();
        delete tsv5;

        group5->close();
        delete group5;
        file5.close();
    } catch(...) {
   throw(Except("Cannot open file: " + strFileName));}
}

void CNAnalysisMethodReference::writeCopyNumberReference()
{
    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");

    AffxString strArrayName;
    if (m_pEngine->getOpt("array-name") == "")
    {
        if ((m_pEngine->isOptDefined("cdf-file")) && (m_pEngine->getOpt("cdf-file") != ""))
        {
            strArrayName = Fs::basename(m_pEngine->getOpt("cdf-file"));
            int iFindIndex = strArrayName.indexOf(".");
            if (iFindIndex != -1) {strArrayName = strArrayName.substring(0, iFindIndex);}
        }
    } else
        {
                strArrayName = m_pEngine->getOpt("array-name");
        }
        AffxString strFileName = m_pEngine->getOpt("reference-file");
    try
    {
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_OPEN);
        affx::File5_Group* group5 = NULL;
        if (!((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2"))))
        {
            group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
        }
        else
        {
            group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
        }
        affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "Parameter", affx::FILE5_DTYPE_STRING, 1024);
        tsv5->set_string(0, 0, "#%guid=" + affxutil::Guid::GenerateNewGuid()); tsv5->writeLevel(0);
        tsv5->set_string(0, 0, "#%affymetrix-array-type=" + strArrayName); tsv5->writeLevel(0);

        vector<string> optionNames;
        m_pEngine->getOptionNames(optionNames,1);
        for(int i=0; i< optionNames.size(); i++) {
            std::string name = optionNames[i];
            std::vector<std::string> vals = m_pEngine->getOptVector(name,1);
            if (vals.size() > 1)
            {
                for(int j=0; j< vals.size(); j++)
                {
                    std::string val = vals[j];
                    if (m_pEngine->getOptBool("cyto2"))
                    {
                        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-" + name + "=" + val); tsv5->writeLevel(0);
                    }
                    else
                    {
                        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-opt-" + name + "=" + val); tsv5->writeLevel(0);
                    }
                }
            }
            else
            {
                std::string val = m_pEngine->getOpt(name,1);
                if (m_pEngine->getOptBool("cyto2"))
                {
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-option-" + name + "=" + val); tsv5->writeLevel(0);
                }
                else
                {
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-opt-" + name + "=" + val); tsv5->writeLevel(0);
                }
            }
        }
        // dump state.
        std::vector<std::string> stateNames;
        m_pEngine->getOptionNames(stateNames);
        for(int i=0; i< stateNames.size(); i++) {
            std::string name = stateNames[i];
            std::vector<std::string> vals = m_pEngine->getOptVector(name);
            if (vals.size() > 1)
            {
                for(int j=0; j< vals.size(); j++)
                {
                    std::string val = vals[j];
                    if (m_pEngine->getOptBool("cyto2"))
                    {
                        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-state-" + name + "=" + val); tsv5->writeLevel(0);
                    }
                    else
                    {
                        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-state-" + name + "=" + val); tsv5->writeLevel(0);
                    }
                }
            }
            else
            {
                std::string val = m_pEngine->getOpt(name);
                if (m_pEngine->getOptBool("cyto2"))
                {
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-state-" + name + "=" + val); tsv5->writeLevel(0);
                }
                else
                {
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-state-" + name + "=" + val); tsv5->writeLevel(0);
                }
            }
        }

        // Write the information entered using the ---meta-data-info option.
        vector<pair<string, string> > metaData = m_pEngine->getMetaDataDescription();
        for (int i = 0; i < metaData.size(); i++)
        {
             pair<string,string> p = metaData[i];
             string strMetaDataString = "#%affymetrix-application-meta-data-info-" + p.first + "=" + p.second;
             tsv5->set_string(0,0, strMetaDataString);
             tsv5->writeLevel(0);
        }

        // Copy over headers from the annotation file.
        if (Annotation::getParams()->size() > 0)
        {
            for (unsigned int uiParamIndex = 0; (uiParamIndex < Annotation::getParams()->size()); uiParamIndex++)
            {
                tsv5->set_string(0, 0, "#%" + StringUtils::ConvertWCSToMBS(Annotation::getParams()->at(uiParamIndex).GetName()) + "=" + StringUtils::ConvertWCSToMBS(Annotation::getParams()->at(uiParamIndex).ToString())); tsv5->writeLevel(0);
            }
        }

        // Copy covariate parameters
        bool covariateValuesLoaded = (getProbeSets()->getCount() > 0 && getProbeSets()->getAt(0)->getNumCovariates() > 0);

        if (!(m_pEngine->isOptDefined("cyto2") && m_pEngine->getOptBool("cyto2")))
        {
            if (!CovariateParams::m_allCovariateMap.empty() && covariateValuesLoaded)
            {
                string name, val;

                name = "signal-command";
                val = m_pEngine->getOpt("signal-adjustment-covariates");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "log2ratio-command";
                val = m_pEngine->getOpt("lr-adjustment-covariates");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "allele-peaks-command";
                val = m_pEngine->getOpt("allele-peaks-adjustment-covariates");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "all";
                val = CovariateParams::translateToParamString(CovariateParams::m_allCovariateMap);
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "signal";
                val = Convert::intVecToString(CovariateParams::m_vSignalCovariates, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "log2ratio";
                val = Convert::intVecToString(CovariateParams::m_vLRCovariates, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "allele-peaks";
                val = Convert::intVecToString(CovariateParams::m_vAPCovariates, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "signal-binning-types";
                val = Convert::strVecToString(CovariateParams::m_vSignalBinningType, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "signal-bin-sizes";
                val = Convert::intVecToString(CovariateParams::m_vSignalBinSizes, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "log2ratio-binning-types";
                val = Convert::strVecToString(CovariateParams::m_vLRBinningType, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "log2ratio-bin-sizes";
                val = Convert::intVecToString(CovariateParams::m_vLRBinSizes, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "log2ratio-iqr-scaling";
                val = Convert::strVecToString(CovariateParams::m_vIQRScaling, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "log2ratio-subtract-from-XY";
                val = Convert::strVecToString(CovariateParams::m_vLRSubractFromXY, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "allele-peaks-binning-types";
                val = Convert::strVecToString(CovariateParams::m_vAPBinningType, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "allele-peaks-bin-sizes";
                val = Convert::intVecToString(CovariateParams::m_vAPBinSizes, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "allele-peaks-coarse-adjustment";
                val = Convert::strVecToString(CovariateParams::m_vCoarseAPAdjust, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "coarse-allele-peak-adjustment-step";
                val = Convert::intVecToString(CovariateParams::m_vCoarseAPAdjustStep, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "coarse-allele-peak-adjustment-window";
                val = Convert::intVecToString(CovariateParams::m_vCoarseAPAdjustWindow, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "coarse-allele-peak-adjustment-point-count";
                val = Convert::intVecToString(CovariateParams::m_vCoarseAPAdjustPointCount, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "coarse-allele-peak-adjustment-bandwidth";
                val = Convert::doubleVecToString(CovariateParams::m_vCoarseAPAdjustBandwidth, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "coarse-allele-peak-adjustment-cutoff";
                val = Convert::doubleVecToString(CovariateParams::m_vCoarseAPAdjustCutoff, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "coarse-allele-peak-adjustment-clean-threshold";
                val = Convert::doubleVecToString(CovariateParams::m_vCoarseAPAdjustCleanthreshold, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "coarse-allele-peak-adjustment-outlier-trim";
                val = Convert::doubleVecToString(CovariateParams::m_vCoarseAPAdjustOutliertrim, ",");
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "master-peaks-point-count";
                val = Convert::toString(CovariateParams::m_masterPeakPointCount);
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "master-peaks-bandwidth";
                val = Convert::toString(CovariateParams::m_masterPeakBandwidth);
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "covariate-bin-peaks-point-count";
                val = Convert::toString(CovariateParams::m_covariatePeakPointCount);
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "covariate-bin-peaks-bandwidth";
                val = Convert::toString(CovariateParams::m_covariatePeakBandwidth);
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);

                name = "kernel-function-selection";
                val = CovariateParams::m_kernelFunctionSelection;
                tsv5->set_string(0, 0, CovariateParams::m_CNReferencePrefix + name + "=" + val);
                tsv5->writeLevel(0);
            }
        }
        tsv5->close();
        delete tsv5;

        int iMaxProbeSetNameLength = 0;
        for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
            iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
        }
        if (!((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2"))))
        {
            tsv5 = group5->openTsv("Reference", affx::FILE5_REPLACE);
        }
        else
        {
            tsv5 = group5->openTsv("MedianSignals", affx::FILE5_REPLACE);
        }
        tsv5->defineColumn(0, 0, "probeset_id", affx::FILE5_DTYPE_STRING, iMaxProbeSetNameLength);
        tsv5->defineColumn(0, 1, "MedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 2, "XXMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 3, "YMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 4, "AAMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 5, "ABMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 6, "BBMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 7, "Chromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 8, "Position", affx::FILE5_DTYPE_INT);
        if (isCytoScanHD)
        {
            tsv5->defineColumn(0, 9, "ProcessFlag", affx::FILE5_DTYPE_INT);
        }
        for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
            tsv5->set_string(0, 0, pobjProbeSet->getProbeSetName());
            tsv5->set_d(0, 1, (double)pobjProbeSet->getMedianSignal());
            tsv5->set_d(0, 2, pobjProbeSet->getXXMedianSignal());
            tsv5->set_d(0, 3, pobjProbeSet->getYMedianSignal());
            tsv5->set_d(0, 4, pobjProbeSet->getAAMedianSignal());
            tsv5->set_d(0, 5, pobjProbeSet->getABMedianSignal());
            tsv5->set_d(0, 6, pobjProbeSet->getBBMedianSignal());
            tsv5->set_i(0, 7, (int)pobjProbeSet->getChromosome());
            tsv5->set_i(0, 8, pobjProbeSet->getPosition());
            if (isCytoScanHD)
            {
                tsv5->set_i(0, 9, pobjProbeSet->getProcessFlag());
            }
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        int medianSignalIndex = CovariateParams::mapCovariateNameToIndex("median-signal");

        // Write covariate parameter values
        if (!(m_pEngine->isOptDefined("cyto2") && m_pEngine->getOptBool("cyto2")))
        {
            if (!CovariateParams::m_allCovariateMap.empty() && covariateValuesLoaded)
            {
                tsv5 = group5->openTsv("Covariates", affx::FILE5_REPLACE);
                int numCovariates = CovariateParams::m_allCovariateMap.size();
                for (int i = 0; i < numCovariates; i++)
                {
                    tsv5->defineColumn(0, i, "C" + Convert::toString(i), affx::FILE5_DTYPE_DOUBLE);
                }
                for (int iRowIndex = 0; (iRowIndex < getProbeSets()->getCount()); iRowIndex++)
                {
                    CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iRowIndex);
                    for (int i = 0; i < numCovariates; i++)
                    {
                        if (i == medianSignalIndex)
                        {
                            tsv5->set_d(0, i, pobjProbeSet->getMedianSignal());
                            pobjProbeSet->setCovariateValue(i, pobjProbeSet->getMedianSignal());
                        }
                        else {
                            tsv5->set_d(0, i, pobjProbeSet->getCovariateValue(i));
                        }
                    }
                    tsv5->writeLevel(0);
                }
                tsv5->close();
                delete tsv5;
            }
        }
        group5->close();
        delete group5;
        file5.close();
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
}

void CNAnalysisMethodReference::warn(int iVerbosity, const std::string& strMessage)
{
    if (m_iWarningCount <= m_pEngine->getOptInt("warning-message-limit"))
    {
        Verbose::warn(iVerbosity, strMessage);
        m_iWarningCount++;
    }
}

/**
 * Utility function for reading a ChipLayout based on the options passed
 * to the base engine (i.e. cdf-file or spf-file)
 */
void CNAnalysisMethodReference::loadLayout(ChipLayout *pLayout, BaseEngine *engine)
{
  string cdfFile = engine->getOpt("cdf-file");
  string spfFile = engine->getOpt("spf-file");
  if(!cdfFile.empty())
    pLayout->openCdfAll(cdfFile);
  else if(!spfFile.empty())
    pLayout->openSpfAll(spfFile);
  else
    Err::errAbort("No cdf-file or spf-file specified");
}
