////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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
 * @file CNReporterCychp.cpp
 *
 * @brief This file contains the CNReporterCychp class members.
 */

//
#include "copynumber/CNReporterCychp.h"
//
#include "copynumber/Annotation.h"
#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNCychp.h"
//
#include "calvin_files/data/src/DataGroupHeader.h"
#include "calvin_files/data/src/DataSetHeader.h"
#include "calvin_files/data/src/GenericData.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/GenericFileWriter.h"
#include "chipstream/BioTypes.h"
#include "util/AffxStatistics.h"
#include "util/AffxString.h"
#include "util/BaseEngine.h"
#include "util/Fs.h"
#include "util/Util.h"
#include "util/AptVersionInfo.h"
//
#include <cassert>
#include <cmath>
#include <vector>
//

void CNReporterCychp::defineOptions(BaseEngine& e)
{
    e.defineOption("","cychp-output", PgOpt::BOOL_OPT,
                    "Report CYCHP files",
                    "false");
}

void CNReporterCychp::checkOptions(BaseEngine& e)
{
}

void CNReporterCychp::run()
{
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Writing cychp output file.");
    isSetup();
    if (!m_pEngine->getOptBool("cychp-output")) {return;}
    Verbose::out(1, "CNReporterCychp::run(...) Start");
    m_uiSegmentID = 0;

    AffxString strFileName = Fs::convertToUncPath(getFileName());
    Verbose::out(1, "Writing file: " + strFileName);
    try
    {
        writeCychpFile(strFileName);

        // Convert from calvin to text format?
        if (m_pEngine->getOptBool("text-output"))
        {
            CalvinToText converter;
            converter.run(strFileName, strFileName + ".txt", true, true, true);
        }
    } catch(...) {throw(Except("CYCHP write failed for file: " + strFileName));}
    Verbose::out(1, "CNReporterCychp::run(...) End");
}

AffxString CNReporterCychp::getFileName()
{
    if ((m_pEngine->isOptDefined("result-files")) && (m_pEngine->getOptVector("result-files").size() > 0))
    {
        AffxString strFileName = m_pEngine->getOptVector("result-files")[getExperiment()->getIndex()];
        if (strFileName != "") {return strFileName;}
    }
    AffxString strExperimentName = getExperiment()->getExperimentName().toLowerCase();
    if (strExperimentName.endsWith(".cel"))
    {
        strExperimentName = getExperiment()->getExperimentName().substring(0, (unsigned int)strExperimentName.length() - 4);
    }
    else {strExperimentName = getExperiment()->getExperimentName();}

    AffxString strOutputDirectory = m_pEngine->getOpt("out-dir");
    AffxString strAnalysisName    = m_pEngine->getOpt("set-analysis-name");
    AffxString strPrefix          = m_pEngine->getOpt("file-name-prefix");
    AffxString strSuffix          = m_pEngine->getOpt("file-name-suffix");
    AffxString strExt             = m_pEngine->getOpt("file-name-ext");

    std::string prefixSeparator   = ".";
    std::string suffixSeparator   = ".";
    std::string analysisSeparator = ".";

    if (strPrefix.empty() || strPrefix == "none") {
        strPrefix = "";
        prefixSeparator = "";
    }
    if (strSuffix.empty() || strSuffix == "none") {
        strSuffix = "";
        suffixSeparator = "";
    }
    if (strAnalysisName.empty()) {
        analysisSeparator = "";
    }
    if (strExt.empty() || strExt == "none") {
        strExt = "cychp";
    }

    AffxString strFileName = Fs::join(strOutputDirectory,
                                      strPrefix          + prefixSeparator   +
                                      strExperimentName  + "."               +
                                      strSuffix          + suffixSeparator   +
                                      strAnalysisName    + analysisSeparator +
                                      strExt);
    return strFileName;
}

AffxString CNReporterCychp::getCelFileName()
{
    AffxString strCelFileName = getExperiment()->getExperimentName();

    if (m_pEngine->isOptDefined("cels"))
    {
        std::vector<std::string> vCelFiles = m_pEngine->getOptVector("cels");
        for (int iIndex = 0; (iIndex < (int)vCelFiles.size()); iIndex++)
        {
            AffxString str = vCelFiles[iIndex];
            if (str.startsWith(strCelFileName))
            {
                strCelFileName = str;
                break;
            }
            else
            {
              if (AffxString(Fs::basename(str)).startsWith(strCelFileName))
                {
                    strCelFileName = str;
                    break;
                }
            }
        }
    }
    m_strARRFileName = "";
    if ((m_pEngine->isOptDefined("arrs")) && (m_pEngine->getOptVector("arrs").size() > 0))
    {
        m_strARRFileName = m_pEngine->getOptVector("arrs")[getExperiment()->getIndex()];
    }
    if (m_strARRFileName == "")
    {
        AffxString strTempFileName = strCelFileName;
        if (strTempFileName.toLowerCase().endsWith(".cel"))
        {
            m_strARRFileName = strCelFileName.substring(0, (unsigned int)strCelFileName.length() - 4) + ".ARR";
        }
        else
        {
            m_strARRFileName = strCelFileName + ".ARR";
        }
    }

    return strCelFileName;
}

void CNReporterCychp::writeCychpFile(const AffxString& strFileName)
{
    const bool addGenotypingTable = m_pEngine->getEngineName() == "CNCytoEngine" && !m_pEngine->getOptBool("cyto2");

    int iSegmentTypeCount = 0;
    for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
    {
        CNAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
        if (pMethod->isSegmentTypeAnalysis())
        {
            for (int iSegmentType = 1; (iSegmentType <= 7); iSegmentType++)
            {
                if (pMethod->getSegmentCount(iSegmentType) > 0)
                {
                    iSegmentTypeCount++;
                }
            }
        }
    }
    // Omitting "Segments" table if there are no segments to report
    const int genotypingTableIndex = iSegmentTypeCount > 0 ? 4 : 3;

    AffxString strCelFileName = getCelFileName();
    affymetrix_calvin_io::GenericData data;
    data.Header().SetFilename(strFileName);
    data.Header().AddDataGroupHdr(affymetrix_calvin_io::DataGroupHeader(L"Chromosomes"));
    data.Header().AddDataGroupHdr(affymetrix_calvin_io::DataGroupHeader(L"ProbeSets"));
    data.Header().AddDataGroupHdr(affymetrix_calvin_io::DataGroupHeader(L"AlgorithmData"));
    if (iSegmentTypeCount > 0) {
        data.Header().AddDataGroupHdr(affymetrix_calvin_io::DataGroupHeader(L"Segments"));
    }
    if (addGenotypingTable) {
        data.Header().AddDataGroupHdr(affymetrix_calvin_io::DataGroupHeader(L"Genotyping"));
    }

    affymetrix_calvin_io::GenericDataHeader* pHeader = data.Header().GetGenericDataHdr();
    pHeader->SetFileTypeId("affymetrix-multi-data-type-analysis");
    pHeader->SetFileCreationTime(DateTime::GetCurrentDateTime().ToString());
    pHeader->SetFileId(affxutil::Guid::GenerateNewGuid());

    addParentHeader(pHeader, strCelFileName);
    Verbose::progressBegin(1, "CNReporterCychp::report(...)", 7, 1, 7);
    addHeaderParams(pHeader);
    addCNReferenceSelectedInfo(pHeader);
    addHeaderAlgorithmParams(pHeader);
    addHeaderClientMetaInfo(pHeader);
    addHeaderChipSummaryParams(pHeader);

    AffxArray<CNChromosome> arChromosomes;
    loadChromosomes(arChromosomes);
    addChromosomesSummaryDataSetHeader(data.Header().GetDataGroup(0), arChromosomes.getCount());
    addProbeSetsCopyNumberDataSetHeader(data.Header().GetDataGroup(1), m_pvProbeSets->getCount());
    int iRowCount = 0;
    for (int i = 0; (i < m_pvProbeSets->getCount()); i++)
    {
        CNProbeSet* p = m_pvProbeSets->getAt(i);
        if ((p->getAllelePeaks1() > 0) || (p->getAllelePeaks2() > 0)) {iRowCount++;}
    }
    addProbeSetsAllelePeaksDataSetHeader(data.Header().GetDataGroup(1), iRowCount);
    int iSnpCount = 0;
    int iVisualizationCount = 0;
    for (int i = 0; (i < m_pvProbeSets->getCount()); i++)
    {
        CNProbeSet* p = m_pvProbeSets->getAt(i);
        if ((p->isProcess()) && (p->processAsSNP()))
        {
            iSnpCount++;
        }
        if ((p->isProcess()) && (p->processAsVisualization()))
        {
            iVisualizationCount++;
        }
    }
    addAlgorithmDataMarkerABSignalDataSetHeader(data.Header().GetDataGroup(2), iVisualizationCount);
    if (iSegmentTypeCount > 0) {
        addSegmentsDataSets(data.Header().GetDataGroup(3));
    }
    if (addGenotypingTable) {
        addGenotypingDataSetHeader(data.Header().GetDataGroup(genotypingTableIndex), iSnpCount);
    }

    affymetrix_calvin_io::GenericFileWriter* writer = new affymetrix_calvin_io::GenericFileWriter(&data.Header());
    writer->WriteHeader();

    std::vector<int> vChromosomesOffsets(1);
    writeDataGroupHeader(writer, 0, vChromosomesOffsets);

    std::vector<int> vProbeSetsOffsets(2);
    writeDataGroupHeader(writer, 1, vProbeSetsOffsets);
    std::vector<int> vAlgorithmDataOffsets(1);
    writeDataGroupHeader(writer, 2, vAlgorithmDataOffsets);

    std::vector<int> vSegmentsOffsets;
    if (iSegmentTypeCount > 0) {
        Verbose::progressStep(1);
        vSegmentsOffsets.resize(iSegmentTypeCount);
        writeDataGroupHeader(writer, 3, vSegmentsOffsets);
    }

    std::vector<int> vGenotypeOffsets(1);
    if (addGenotypingTable) {
        writeDataGroupHeader(writer, genotypingTableIndex, vGenotypeOffsets);
    }

    Verbose::progressStep(1);
    writeChromosomesSummaryDataSet(writer, 0, 0, vChromosomesOffsets, arChromosomes);
    Verbose::progressStep(1);
    writeProbeSetsCopyNumberDataSet(writer, 1, 0, vProbeSetsOffsets, *m_pvProbeSets);
    Verbose::progressStep(1);
    writeProbeSetsAllelePeaksDataSet(writer, 1, 1, vProbeSetsOffsets, *m_pvProbeSets);
    Verbose::progressStep(1);
    writeAlgorthmDataMarkerABSignalDataSet(writer, 2, 0, vAlgorithmDataOffsets, *m_pvProbeSets);
    Verbose::progressStep(1);
    if (iSegmentTypeCount > 0)
    {
        int iDataSetIndex = 0;
        for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
        {
            CNAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
            if (pMethod->isSegmentTypeAnalysis())
            {
                for (int iSegmentType = 1; (iSegmentType <= 7); iSegmentType++)
                {
                    if (pMethod->getSegmentCount(iSegmentType) > 0)
                    {
                        writeSegmentsDataSet(writer, 3, iDataSetIndex, vSegmentsOffsets, pMethod->getSegments(), (short)iSegmentType);
                        iDataSetIndex++;
                    }
                }
            }
        }
    }
    if (addGenotypingTable) {
        Verbose::progressStep(1);
        writeGenInfoDataSet(writer, genotypingTableIndex, 0, vGenotypeOffsets, *m_pvProbeSets);
    }
    Verbose::progressStep(1);

    delete writer;
    arChromosomes.deleteAll();
    Verbose::progressEnd(1, "Done.");
}

void CNReporterCychp::addParentHeader(  affymetrix_calvin_io::GenericDataHeader* pHeader,
                                        AffxString& strFileName)
{
    // Copy over headers from the cel file.
    try
    {
        affymetrix_calvin_io::GenericFileReader reader;
        affymetrix_calvin_io::GenericData genericData;
        reader.SetFilename(strFileName.c_str());
        reader.Open(genericData);
        pHeader->AddParent(*genericData.Header().GetGenericDataHdr());
        genericData.Clear();
    } catch(...) {Verbose::out(1, "WARNING: Load parent header failed for file: " + strFileName);}
}

void CNReporterCychp::addHeaderParams(affymetrix_calvin_io::GenericDataHeader* pHeader)
{
    affymetrix_calvin_parameter::ParameterNameValueType param;

    param.SetName(L"affymetrix-algorithm-name");
    param.SetValueText(L"CYTO2");
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-version");
    
    param.SetValueText(StringUtils::ConvertMBSToWCS(AptVersionInfo::reportRelease()));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-array-type");
    param.SetValueText(StringUtils::ConvertMBSToWCS(getArrayName()));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-param-ArraySet");
    param.SetValueText(StringUtils::ConvertMBSToWCS(getArrayName()));
    pHeader->AddNameValParam(param);

    param.SetName(L"program-name");
    param.SetValueText(StringUtils::ConvertMBSToWCS(m_pEngine->getOpt("program-name")));
    pHeader->AddNameValParam(param);

    param.SetName(L"program-version");
    param.SetValueText(StringUtils::ConvertMBSToWCS(m_pEngine->getOpt("program-version")));
    pHeader->AddNameValParam(param);

    param.SetName(L"program-company");
    param.SetValueText(StringUtils::ConvertMBSToWCS("Affymetrix"));
    pHeader->AddNameValParam(param);

    param.SetName(L"create-date");
    param.SetValueText(StringUtils::ConvertMBSToWCS(Util::getTimeStamp()));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-meta-data-version");
    param.SetValueInt32(1);
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-param-file-name");
    param.SetValueText(StringUtils::ConvertMBSToWCS(Fs::basename(m_pEngine->getXMLParameterFileName())));
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-algorithm-param-file-guid");
    param.SetValueText(StringUtils::ConvertMBSToWCS(m_pEngine->getXMLParameterFileGuid()));
    pHeader->AddNameValParam(param);

}

void CNReporterCychp::addHeaderClientMetaInfo(affymetrix_calvin_io::GenericDataHeader* pHeader)
{
  ParameterNameValueType param;
  vector<pair<string, string> > metaData = m_pEngine->getMetaDataDescription();
  for (int i = 0; i < metaData.size(); i++) {
    pair<string,string> p = metaData[i];
    string key = "affymetrix-application-meta-data-info-" + p.first;
    param.SetName(StringUtils::ConvertMBSToWCS(key));
    param.SetValueText(StringUtils::ConvertMBSToWCS(p.second));
    pHeader->AddNameValParam(param);
  }

}


void CNReporterCychp::addHeaderAlgorithmParams(affymetrix_calvin_io::GenericDataHeader* pHeader)
{
    affymetrix_calvin_parameter::ParameterNameValueType param;

    param.SetName(L"affymetrix-algorithm-param-cychp-file-version");
    param.SetValueInt32(1);
    pHeader->AddNameValParam(param);

    // dump options
    vector<string> optionNames;
    m_pEngine->getOptionNames(optionNames,1);
    std::vector<PgOpt::PgOptType> optionTypes;
    m_pEngine->getOptionTypes(optionTypes,1);
    for(int i=0; i< optionNames.size(); i++) {
        std::string name = optionNames[i];
        std::vector<std::string> vals = m_pEngine->getOptVector(name,1);
        if (vals.size() > 1)
        {
            for(int j=0; j< vals.size(); j++)
            {
                std::string val = vals[j];
                loadParam("affymetrix-algorithm-param-option-" + name + "-" + ::getInt(j+1), optionTypes[i], val, param);
                pHeader->AddNameValParam(param);
            }
        }
        else
        {
            std::string val = m_pEngine->getOpt(name,1);
            loadParam("affymetrix-algorithm-param-option-" + name, optionTypes[i], val, param);
            pHeader->AddNameValParam(param);
        }
    }
    // dump state.
    std::vector<std::string> stateNames;
    m_pEngine->getOptionNames(stateNames);
    std::vector<PgOpt::PgOptType> stateTypes;
    m_pEngine->getOptionTypes(stateTypes);
    for(int i=0; i< stateNames.size(); i++) {
        std::string name = stateNames[i];
        std::vector<std::string> vals = m_pEngine->getOptVector(name);
        if (vals.size() > 1)
        {
            for(int j=0; j< vals.size(); j++)
            {
                std::string val = vals[j];
                loadParam("affymetrix-algorithm-param-state-" + name + "-" + ::getInt(j+1), stateTypes[i], val, param);
                pHeader->AddNameValParam(param);
            }
        }
        else
        {
            std::string val = m_pEngine->getOpt(name);
            loadParam("affymetrix-algorithm-param-state-" + name, stateTypes[i], val, param);
            pHeader->AddNameValParam(param);
        }
    }

    for (int iIndex = 0; (iIndex < (int)CNAnalysisMethod::getParams()->size()); iIndex++)
    {
        pHeader->AddNameValParam(CNAnalysisMethod::getParams()->at(iIndex));
    }
    // Copy over cel file names.
    if (m_pEngine->isOptDefined("create-reference"))
    {
        if (m_pEngine->getOptBool("create-reference"))
        {
            if (CNAnalysisMethod::getCelFileParams()->size() > 0)
            {
                for (unsigned int uiParamIndex = 0; (uiParamIndex < CNAnalysisMethod::getCelFileParams()->size()); uiParamIndex++)
                {
                    pHeader->AddNameValParam(CNAnalysisMethod::getCelFileParams()->at(uiParamIndex));
                }
            }
        }
    }

    // Copy over headers from the annotation file.
    if (Annotation::getParams()->size() > 0)
    {
        for (unsigned int uiParamIndex = 0; (uiParamIndex < Annotation::getParams()->size()); uiParamIndex++)
        {
            pHeader->AddNameValParam(Annotation::getParams()->at(uiParamIndex));
        }
    }
}

void CNReporterCychp::addCNReferenceHeader(affymetrix_calvin_io::GenericDataHeader* pHeader)
{

    affymetrix_calvin_parameter::ParameterNameValueType param;
    int iNumberOfHeaders = m_pobjExperiment->getCNReferenceHeader()->size();

    for(int iIndex=0; iIndex<iNumberOfHeaders; iIndex++)
    {
        AffxString strParameterName = "affymetrix-CNReferenceHeader-";
        AffxString strHeaderLine =( *(m_pobjExperiment->getCNReferenceHeader()))[iIndex];

        AffxString strTempString1  = strHeaderLine.getKey();
        AffxString strTempString2  = strTempString1.removeCruft();
        strParameterName += strTempString2;
        param.SetName(StringUtils::ConvertMBSToWCS(strParameterName));

        AffxString strValue = strHeaderLine.getValue();
        param.SetValueAscii(strValue);

        pHeader->AddNameValParam(param);
    }
}

void CNReporterCychp::addCNReferenceSelectedInfo(affymetrix_calvin_io::GenericDataHeader* pHeader)
{
    affymetrix_calvin_parameter::ParameterNameValueType param;
    int iNumberOfHeaders = m_pobjExperiment->getCNReferenceHeader()->size();

    // Options from the CN reference header to show in this CYCHP header
    std::set<std::string> allowedHeaderItems;
    allowedHeaderItems.insert("guid");
    allowedHeaderItems.insert("algorithm-param-apt-state-reference-file");

    for(int iIndex=0; iIndex<iNumberOfHeaders; iIndex++)
    {
        AffxString strParameterName = "affymetrix-CNReferenceHeader-";
        AffxString strHeaderLine =( *(m_pobjExperiment->getCNReferenceHeader()))[iIndex];

        AffxString strTempString1  = strHeaderLine.getKey();
        AffxString strTempString2  = strTempString1.removeCruft();
        if (allowedHeaderItems.count(strTempString2))
        {
            strParameterName += strTempString2;
            param.SetName(StringUtils::ConvertMBSToWCS(strParameterName));
    
            AffxString strValue = strHeaderLine.getValue();
            param.SetValueAscii(strValue);
    
            pHeader->AddNameValParam(param);
        }
    }
}

void CNReporterCychp::addHeaderChipSummaryParams(affymetrix_calvin_io::GenericDataHeader* pHeader)
{
    affymetrix_calvin_parameter::ParameterNameValueType param;

    param.SetName(L"affymetrix-chipsummary-raw-intensity-ratio-gender");
    param.SetValueAscii(m_pobjExperiment->getRawIntensityRatioGender());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-raw-intensity-ratio");
    param.SetValueFloat(m_pobjExperiment->getRawIntensityRatio());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-XX-ratio");
    param.SetValueFloat(m_pobjExperiment->getXXRatio());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-Y-ratio");
    param.SetValueFloat(m_pobjExperiment->getYRatio());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-XX");
    param.SetValueInt32(m_pobjExperiment->hasXX() ? 1 : 0);
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-Y");
    param.SetValueInt32(m_pobjExperiment->hasY() ? 1 : 0);
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-Y-gender-call-computed");
    param.SetValueInt32(m_pobjExperiment->getCNCallGenderComputed() ? 1 : 0);
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-Y-gender-call");
    param.SetValueAscii(m_pobjExperiment->getCNCallGender());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-Y-gender-not-zero");
    param.SetValueInt32(m_pobjExperiment->getCNCallGenderNotZero());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-Y-gender-Y-count");
    param.SetValueInt32(m_pobjExperiment->getCNCallGenderYCount());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-Y-gender-confidence");
    param.SetValueFloat(m_pobjExperiment->getCNCallGenderConfidence());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-Y-gender-ratio");
    param.SetValueFloat(m_pobjExperiment->getCNCallGenderRatio());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-MedianAutosomeMedian");
    param.SetValueFloat(m_pobjExperiment->getMedianAutosomeMedian());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-MAPD");
    param.SetValueFloat(m_pobjExperiment->getMadDiffCN());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-iqr");
    param.SetValueFloat(m_pobjExperiment->getIqr());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-all_probeset_rle_mean");
    param.SetValueFloat(m_pobjExperiment->getMeanAbsRle());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-gc-correction-size");
    param.SetValueFloat(m_pobjExperiment->getGCCorrectionMetric());
    pHeader->AddNameValParam(param);

    param.SetName(L"affymetrix-chipsummary-call-rate");
    param.SetValueFloat(m_pobjExperiment->getCallRate());
    pHeader->AddNameValParam(param);

    if (CNExperiment::getQCMetricColumnNames()->getCount() != m_pobjExperiment->getQCMetricColumnValues()->getCount())
    {
        Err::errAbort("Mismatch occurred between QC Metric column names count and QC Metric column values count.");
    }
    for (int iMetricIndex = 0; (iMetricIndex < CNExperiment::getQCMetricColumnNames()->getCount()); iMetricIndex++)
    {
        if ((*CNExperiment::getQCMetricColumnNames()->getAt(iMetricIndex)).endsWith("_gender"))
        {
            param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-" + *CNExperiment::getQCMetricColumnNames()->getAt(iMetricIndex)));
            param.SetValueAscii(*m_pobjExperiment->getQCMetricColumnValues()->getAt(iMetricIndex));
            pHeader->AddNameValParam(param);
        }
        else
        {
            param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-" + *CNExperiment::getQCMetricColumnNames()->getAt(iMetricIndex)));
            param.SetValueFloat((float)::getDouble(*m_pobjExperiment->getQCMetricColumnValues()->getAt(iMetricIndex)));
            pHeader->AddNameValParam(param);
        }
    }

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-wave-count"));
    param.SetValueFloat(m_pobjExperiment->getNumberOfWavesUsed());
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-sample-median-cn-state"));
    param.SetValueFloat(m_pobjExperiment->getMedianCnState());
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-sample-hom-frequency"));
    param.SetValueFloat(m_pobjExperiment->getHomFrequency());
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-sample-het-frequency"));
    param.SetValueFloat(m_pobjExperiment->getHetFrequency());
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-median-raw-intensity"));
    param.SetValueFloat(m_pobjExperiment->getMedianRawIntensity());
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-antigenomic-ratio"));
    param.SetValueFloat(m_pobjExperiment->getAntigenomicRatio());
    pHeader->AddNameValParam(param);

    if (m_pEngine->getOptBool("cyto2"))
    {
        if (m_pEngine->isOptDefined("snp-qc-use-contrast"))
        {
            param.SetName(L"affymetrix-chipsummary-snp-qc-algorithm");
            if (m_pEngine->getOptBool("snp-qc-use-contrast"))
            {
                param.SetValueAscii("CQC");
            }
            else
            {
                param.SetValueAscii("PVQC");
            }
            pHeader->AddNameValParam(param);
        }
        param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-snp-qc"));
        param.SetValueFloat(m_pobjExperiment->getSNPQC());
        pHeader->AddNameValParam(param);
    }
    else
    {
        param.SetName(L"affymetrix-chipsummary-snp-qc-algorithm");
        param.SetValueAscii("EMQC");
        pHeader->AddNameValParam(param);

        param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-snp-qc"));
        param.SetValueFloat(m_pobjExperiment->getSNPQC());
        pHeader->AddNameValParam(param);

        param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-raw-snp-qc"));
        param.SetValueFloat(m_pobjExperiment->getRawSNPQC());
        pHeader->AddNameValParam(param);
    }

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-genome-LOH"));
    param.SetValueFloat(m_pobjExperiment->getGenomeLOH());
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-autosome-genome-LOH"));
    param.SetValueFloat(m_pobjExperiment->getAutosomeGenomeLOH());
    pHeader->AddNameValParam(param);

    // CN waviness parameters
    int segCountLoss, segCountGain;
    float sd;
    wavinessSegCounts(segCountLoss, segCountGain, sd);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-waviness-seg-count-loss"));
    param.SetValueInt32(segCountLoss);
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-waviness-seg-count-gain"));
    param.SetValueInt32(segCountGain);
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-waviness-seg-count"));
    if (segCountLoss == -1 && segCountGain == -1) {
        param.SetValueInt32(-1);
    } else {
        param.SetValueInt32(segCountLoss + segCountGain);
    }
    pHeader->AddNameValParam(param);

    param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-waviness-sd"));
    param.SetValueFloat(sd);
    pHeader->AddNameValParam(param);

    if (m_pEngine->isOptDefined("wave-correction-log2ratio-adjustment-method")     &&
        m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") != "none" &&
        m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") != "")
    {
        param.SetName(StringUtils::ConvertMBSToWCS("affymetrix-chipsummary-waviness-amplitudes"));
        param.SetValueAscii(wavinessAmplitudes());
        pHeader->AddNameValParam(param);
    }

    param.SetName(L"affymetrix-chipsummary-log2-ratio-gradient");
    param.SetValueFloat(m_pobjExperiment->getL2Gradient());
    pHeader->AddNameValParam(param);
}

void CNReporterCychp::loadChromosomes(AffxArray<CNChromosome>& arChromosomes)
{
    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
    const double alhpaCNCalibrate = m_pEngine->getOptDouble("alpha-cn-calibrate");
    const double alhpa_X_CNCalibrate = m_pEngine->getOptDouble("alpha-X-cn-calibrate");
    const double alhpa_Y_CNCalibrate = m_pEngine->getOptDouble("alpha-Y-cn-calibrate");
    const double betaCNCalibrate = m_pEngine->getOptDouble("beta-cn-calibrate");
    const double beta_X_CNCalibrate = m_pEngine->getOptDouble("beta-X-cn-calibrate");
    const double beta_Y_CNCalibrate = m_pEngine->getOptDouble("beta-Y-cn-calibrate");
    int iXChromosome = m_pEngine->getOptInt("xChromosome");
    int iYChromosome = m_pEngine->getOptInt("yChromosome");

    int iMaxChromosome = 0;
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        if (!pobjProbeSet->isProcess()) {continue;}
        iMaxChromosome = Max(iMaxChromosome, (int)pobjProbeSet->getChromosome());
    }
    arChromosomes.deleteAll();
    for (int iChromosome = 0; iChromosome <= iMaxChromosome; iChromosome++)
    {
        float fMinLog2Ratio = 1000000;
        float fMaxLog2Ratio = -1000000;
        std::vector<float> vLog2Ratios;
        vLog2Ratios.clear();
        int iMarkerCount = 0;
        int iTotalMarkerCount = 0;
        int iSnpCount = 0;
        int iHomCount = 0;
        int iHetCount = 0;
        int iStartIndex = -1;
        float fConfidenceThreshold = CNAnalysisMethod::getConfidenceThreshold(m_pEngine->getOpt("brlmmp-parameters"));
        for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
            if (!pobjProbeSet->isProcess() ) {continue;}
            if (pobjProbeSet->getChromosome() == iChromosome)
            {
                if(pobjProbeSet->processAsCN())
                {
                    vLog2Ratios.push_back(pobjProbeSet->getLog2Ratio());
                    fMinLog2Ratio = Min(fMinLog2Ratio, pobjProbeSet->getLog2Ratio());
                    fMaxLog2Ratio = Max(fMaxLog2Ratio, pobjProbeSet->getLog2Ratio());
                    iMarkerCount++;
                }
                iTotalMarkerCount++;
                if (iStartIndex == -1) {iStartIndex = iIndex;}
                if (pobjProbeSet->processAsSNP())
                {
                    if (isCytoScanHD && (pobjProbeSet->getGenotypeConfidence() >= fConfidenceThreshold)) {continue;}
                    if ((pobjProbeSet->getGenotypeCall() == 0) || (pobjProbeSet->getGenotypeCall() == 2))
                    {
                        iHomCount++;
                    }
                    if (pobjProbeSet->getGenotypeCall() == 1)
                    {
                        iHetCount++;
                    }
                    iSnpCount++;
                }
            }
        }
        if (iMarkerCount > 0)
        {
            AffxMultiDimensionalArray<float> vCytoScanLog2Ratios(iMarkerCount);
            AffxMultiDimensionalArray<float> vCnState(iMarkerCount);
            AffxMultiDimensionalArray<float> vMosaicismMixture(iMarkerCount);
            int iMarkerIndex = 0;
            for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
                if (!pobjProbeSet->isProcess() || !pobjProbeSet->processAsCN()) {continue;}
                if (pobjProbeSet->getChromosome() == iChromosome)
                {
                    vCytoScanLog2Ratios.set(iMarkerIndex, pobjProbeSet->getLog2Ratio());
                    float fCalibratedLog2Ratio = 0.0;
                    if (iChromosome < iXChromosome) {
                        fCalibratedLog2Ratio = (float)exp(((pobjProbeSet->getLog2Ratio() / alhpaCNCalibrate) + betaCNCalibrate) * log(2.0));
                    }
                    else if (iChromosome == iXChromosome) {
                        fCalibratedLog2Ratio = (float)exp(((pobjProbeSet->getLog2Ratio() / alhpa_X_CNCalibrate) + beta_X_CNCalibrate) * log(2.0));
                    }
                    else if (iChromosome == iYChromosome) {
                        fCalibratedLog2Ratio = (float)exp(((pobjProbeSet->getLog2Ratio() / alhpa_Y_CNCalibrate) + beta_Y_CNCalibrate) * log(2.0));
                    }
                    vCnState.set(iMarkerIndex, fCalibratedLog2Ratio);
                    if (!isCytoScanHD && (pobjProbeSet->getMosaicismMixture() == numeric_limits<float>::quiet_NaN()))
                    {
                        vMosaicismMixture.set(iMarkerIndex, 0);
                    }
                    else
                    {
                        vMosaicismMixture.set(iMarkerIndex, pobjProbeSet->getMosaicismMixture());
                    }
                    iMarkerIndex++;
                }
            }
            CNChromosome* p = new CNChromosome;
            p->m_ucChromosome = (unsigned char)iChromosome;
            if (iChromosome == m_pEngine->getOptInt("xChromosome")) {p->m_strDisplay = "X";}
            else if (iChromosome == m_pEngine->getOptInt("yChromosome")) {p->m_strDisplay = "Y";}
            else {p->m_strDisplay = ::getInt(iChromosome);}
            p->m_uiStartIndex = iStartIndex;
            p->m_uiMarkerCount = iTotalMarkerCount;
            p->m_fMinSignal = fMinLog2Ratio;
            p->m_fMaxSignal = fMaxLog2Ratio;
            if (isCytoScanHD) {p->m_fMedianSignal = vCytoScanLog2Ratios.median();}
            else
            {
                std::nth_element(vLog2Ratios.begin(), vLog2Ratios.begin()+iMarkerCount/2, vLog2Ratios.end());
                p->m_fMedianSignal = vLog2Ratios[iMarkerCount/2];
            }
            p->m_fMedianCnState = (float)vCnState.median();
            p->m_fHomFrequency = (float)((double)iHomCount / (double)iSnpCount);
            p->m_fHetFrequency = (float)((double)iHetCount / (double)iSnpCount);
            p->m_fMosaicism = (float)vMosaicismMixture.median();
            p->m_fChromosomeLOH = numeric_limits<float>::quiet_NaN();
            for(int iIndex = 0; (iIndex < getExperiment()->getNumberOfChromosomesToReport()); iIndex++)
            {
                std::pair<char, float> pairLOH;
                try {pairLOH = getExperiment()->getChromosomeLOH(iIndex);} catch(...) {continue;}
                if(pairLOH.first == iChromosome)
                {
                    p->m_fChromosomeLOH = pairLOH.second;
                    break;
                }
                else
                {
                    if (isCytoScanHD) {p->m_fChromosomeLOH = numeric_limits<float>::quiet_NaN();}
                    else {p->m_fChromosomeLOH = 0.0;}
                }
            }

            arChromosomes.add(p);
        }
    }
}

void CNReporterCychp::addChromosomesSummaryDataSetHeader(       affymetrix_calvin_io::DataGroupHeader& group,
                                                                int iRowCount)
{
    affymetrix_calvin_io::DataSetHeader dsh;
    dsh.SetRowCnt(iRowCount);
    dsh.SetName(L"Summary");
    dsh.AddUByteColumn(L"Chromosome");
    dsh.AddAsciiColumn(L"Display", 4);
    dsh.AddUIntColumn(L"StartIndex");
    dsh.AddUIntColumn(L"MarkerCount");
    dsh.AddFloatColumn(L"MinSignal");
    dsh.AddFloatColumn(L"MaxSignal");
    dsh.AddFloatColumn(L"MedianCnState");
    dsh.AddFloatColumn(L"HomFrequency");
    dsh.AddFloatColumn(L"HetFrequency");
    dsh.AddFloatColumn(L"Mosaicism");
    dsh.AddFloatColumn(L"LOH");
    dsh.AddFloatColumn(L"MedianSignal");
    group.AddDataSetHdr(dsh);
}

int CNReporterCychp::getMaxProbeSetNameLength()
{
    int iMaxProbeSetNameLength = 0;
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
    }
    return iMaxProbeSetNameLength;
}

void CNReporterCychp::addProbeSetsCopyNumberDataSetHeader(      affymetrix_calvin_io::DataGroupHeader& group,
                                                                int iRowCount)
{
    int iMaxProbeSetNameLength = getMaxProbeSetNameLength();
    affymetrix_calvin_io::DataSetHeader dsh;
    dsh.SetRowCnt(iRowCount);
    dsh.SetName(L"CopyNumber");
    dsh.AddAsciiColumn(L"ProbeSetName", iMaxProbeSetNameLength + 1);
    dsh.AddUByteColumn(L"Chromosome");
    dsh.AddUIntColumn(L"Position");
    dsh.AddFloatColumn(L"Log2Ratio");
    dsh.AddFloatColumn(L"WeightedLog2Ratio");
    dsh.AddFloatColumn(L"SmoothSignal");
    group.AddDataSetHdr(dsh);
}

void CNReporterCychp::addProbeSetsAllelePeaksDataSetHeader(     affymetrix_calvin_io::DataGroupHeader& group,
                                                                int iRowCount)
{
    int iMaxProbeSetNameLength = getMaxProbeSetNameLength();
    affymetrix_calvin_io::DataSetHeader dsh;
    dsh.SetRowCnt(iRowCount);
    dsh.SetName(L"AllelePeaks");
    dsh.AddAsciiColumn(L"ProbeSetName", iMaxProbeSetNameLength + 1);
    dsh.AddUByteColumn(L"Chromosome");
    dsh.AddUIntColumn(L"Position");
    dsh.AddUIntColumn(L"AllelePeaks0");
    dsh.AddUIntColumn(L"AllelePeaks1");
    group.AddDataSetHdr(dsh);
}

void CNReporterCychp::addAlgorithmDataMarkerABSignalDataSetHeader(      affymetrix_calvin_io::DataGroupHeader& group,
                                                                        int iRowCount)
{
    const bool cyto2Run = m_pEngine->getOptBool("cyto2");

    affymetrix_calvin_io::DataSetHeader dsh;
    dsh.SetRowCnt(iRowCount);
    dsh.SetName(L"MarkerABSignal");
    dsh.AddUIntColumn(L"Index");
    if (cyto2Run) {
        dsh.AddFloatColumn(L"ASignal");
        dsh.AddFloatColumn(L"BSignal");
    }
    dsh.AddFloatColumn(L"SCAR");
    group.AddDataSetHdr(dsh);
}

void CNReporterCychp::addGenotypingDataSetHeader(      affymetrix_calvin_io::DataGroupHeader& group,
                                                       int iRowCount)
{
    affymetrix_calvin_io::DataSetHeader dsh;
    dsh.SetRowCnt(iRowCount);
    dsh.SetName(L"Calls");
    dsh.AddUIntColumn(L"Index");
    dsh.AddUByteColumn(L"Call");
    dsh.AddFloatColumn(L"Confidence");
    dsh.AddUByteColumn(L"ForcedCall");
    dsh.AddFloatColumn(L"ASignal");
    dsh.AddFloatColumn(L"BSignal");
    dsh.AddFloatColumn(L"SignalStrength");
    dsh.AddFloatColumn(L"Contrast");
    group.AddDataSetHdr(dsh);
}

void CNReporterCychp::addSegmentsDataSets(affymetrix_calvin_io::DataGroupHeader& group)
{
    for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++)
    {
        CNAnalysisMethod* pMethod = m_pvMethods->at(iIndex);
        if (pMethod->isSegmentTypeAnalysis())
        {
            for (int iSegmentType = 1; (iSegmentType <= 7); iSegmentType++)
            {
                int iCount = pMethod->getSegmentCount(iSegmentType);
                if (iCount > 0)
                {
                    setupSegmentsDataSetHeader(group, iCount, (short)iSegmentType);
                }
            }
        }
    }
}

void CNReporterCychp::setupSegmentsDataSetHeader(    affymetrix_calvin_io::DataGroupHeader& group,
                                                     int iRowCount,
                                                     short nSegmentType)
{
    AffxString strName = CNSegment::getSegmentTypeString(nSegmentType);
    affymetrix_calvin_io::DataSetHeader dsh;
    dsh.SetRowCnt(iRowCount);
    dsh.SetName(StringUtils::ConvertMBSToWCS(strName));
    dsh.AddUIntColumn(L"SegmentID");
    dsh.AddUByteColumn(L"Chromosome");
    dsh.AddUIntColumn(L"StartPosition");
    dsh.AddUIntColumn(L"StopPosition");
    dsh.AddIntColumn(L"MarkerCount");
    dsh.AddUIntColumn(L"MeanMarkerDistance");
    if (strName != "CN")
    {
        dsh.AddUByteColumn(StringUtils::ConvertMBSToWCS(strName));
        dsh.AddFloatColumn(L"Confidence");
    }
    if (strName == "CN")
    {
        dsh.AddFloatColumn(L"State");
        dsh.AddFloatColumn(L"Confidence");
    }
    if (strName == "Mosaicism")
    {
        dsh.AddFloatColumn(L"State");
        dsh.AddFloatColumn(L"Mixture");
    }
    group.AddDataSetHdr(dsh);
}

void CNReporterCychp::writeDataGroupHeader(     affymetrix_calvin_io::GenericFileWriter* pWriter,
                                                int iGroupIndex,
                                                std::vector<int>& vOffsets)
{
    affymetrix_calvin_io::DataGroupWriter& group = pWriter->GetDataGroupWriter(iGroupIndex);
    group.WriteHeader();
    int iDataSetCount = (int)vOffsets.size();
    for (int iDataSetIndex = 0; (iDataSetIndex < iDataSetCount); iDataSetIndex++)
    {
        group.GetDataSetWriter(iDataSetIndex).WriteHeader();
        affymetrix_calvin_io::DataSetWriter& set = group.GetDataSetWriter(iDataSetIndex);
        vOffsets[iDataSetIndex] = pWriter->GetFilePos();
        pWriter->SeekFromCurrentPos(set.GetDataSetSize() + 1);
        set.UpdateNextDataSetOffset();
        group.UpdateNextDataGroupPos();
    }
}

void CNReporterCychp::writeChromosomesSummaryDataSet(  affymetrix_calvin_io::GenericFileWriter* pWriter,
                                                        int iGroupIndex,
                                                        int iSetIndex,
                                                        std::vector<int>& vOffsets,
                                                        AffxArray<CNChromosome>& arChromosomes)
{
    affymetrix_calvin_io::DataSetWriter& set = pWriter->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    pWriter->SeekFromBeginPos(vOffsets[iSetIndex]);
    int iSize = 0;
    iSize += sizeof(unsigned char);
    iSize += sizeof(unsigned int);
    iSize += 4;
    iSize += sizeof(unsigned int);
    iSize += sizeof(unsigned int);
    iSize += sizeof(float);
    iSize += sizeof(float);
    iSize += sizeof(float);
    iSize += sizeof(float);
    iSize += sizeof(float);
    iSize += sizeof(float);
    iSize += sizeof(float);
    int iStructCount = 100000;
    char* pBuffer = new char[iSize * iStructCount];
    int iIndex = 0;
    int iCount = 0;
    AffxString str;
    for (int iChromosomeIndex = 0; (iChromosomeIndex < arChromosomes.getCount()); iChromosomeIndex++)
    {
        CNChromosome* p = arChromosomes.getAt(iChromosomeIndex);
        loadBuffer(pBuffer, iIndex, (unsigned char)p->m_ucChromosome);
        str = prepareAscii(p->m_strDisplay, 4);
        loadBuffer(pBuffer, iIndex, str, 4);
        loadBuffer(pBuffer, iIndex, (unsigned int)p->m_uiStartIndex);
        loadBuffer(pBuffer, iIndex, (unsigned int)p->m_uiMarkerCount);
        loadBuffer(pBuffer, iIndex, p->m_fMinSignal);
        loadBuffer(pBuffer, iIndex, p->m_fMaxSignal);
        loadBuffer(pBuffer, iIndex, p->m_fMedianCnState);
        loadBuffer(pBuffer, iIndex, p->m_fHomFrequency);
        loadBuffer(pBuffer, iIndex, p->m_fHetFrequency);
        loadBuffer(pBuffer, iIndex, p->m_fMosaicism);
        loadBuffer(pBuffer, iIndex, p->m_fChromosomeLOH);
        loadBuffer(pBuffer, iIndex, p->m_fMedianSignal);

        iCount++;
        if (iCount >= iStructCount)
        {
            set.WriteBuffer(pBuffer, iIndex);
            iCount = 0;
            iIndex = 0;
        }
    }
    if (iCount > 0) {set.WriteBuffer(pBuffer, iIndex);}
    delete[] pBuffer;
}

void CNReporterCychp::writeProbeSetsCopyNumberDataSet( affymetrix_calvin_io::GenericFileWriter* pWriter,
                                                        int iGroupIndex,
                                                        int iSetIndex,
                                                        std::vector<int>& vOffsets,
                                                        CNProbeSetArray& arProbeSets)
{
    int iXChromosome = m_pEngine->getOptInt("xChromosome");
    int iYChromosome = m_pEngine->getOptInt("yChromosome");
    int iMaxProbeSetNameLength = getMaxProbeSetNameLength();
    affymetrix_calvin_io::DataSetWriter& set = pWriter->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    pWriter->SeekFromBeginPos(vOffsets[iSetIndex]);

    int iSize = 0;
    iSize += sizeof(unsigned int);
    iSize += (iMaxProbeSetNameLength + 1);
    iSize += sizeof(unsigned char);
    iSize += sizeof(unsigned int);
    iSize += sizeof(float);
    iSize += sizeof(float);
    iSize += sizeof(float);
    int iStructCount = 100000;
    char* pBuffer = new char[iSize * iStructCount];
    int iIndex = 0;
    int iCount = 0;
    AffxString str;
    for (int iProbeSetIndex = 0; (iProbeSetIndex < arProbeSets.getCount()); iProbeSetIndex++)
    {
        CNProbeSet* p = arProbeSets.getAt(iProbeSetIndex);
        if (!p->isProcess()) {continue;}
        if( !p->processAsCN())
        {
                str = prepareAscii(p->getProbeSetName(), iMaxProbeSetNameLength + 1);
                loadBuffer(pBuffer, iIndex, str, (iMaxProbeSetNameLength + 1));
                loadBuffer(pBuffer, iIndex, (unsigned char)p->getChromosome());
                loadBuffer(pBuffer, iIndex, (unsigned int)p->getPosition());
                loadBuffer(pBuffer, iIndex, numeric_limits<float>::quiet_NaN());
                loadBuffer(pBuffer, iIndex, numeric_limits<float>::quiet_NaN());
                loadBuffer(pBuffer, iIndex, numeric_limits<float>::quiet_NaN());
                iCount++;
        }
        else
        {
                str = prepareAscii(p->getProbeSetName(), iMaxProbeSetNameLength + 1);
                loadBuffer(pBuffer, iIndex, str, (iMaxProbeSetNameLength + 1));
                loadBuffer(pBuffer, iIndex, (unsigned char)p->getChromosome());
                loadBuffer(pBuffer, iIndex, (unsigned int)p->getPosition());
                loadBuffer(pBuffer, iIndex, p->getLog2Ratio());
                loadBuffer(pBuffer, iIndex, p->getLog2RatioMedianSmooth());
                if (p->getChromosome() < iXChromosome) {
                    loadBuffer(pBuffer, iIndex, p->getCalibratedSmoothedLog2Ratio(m_pEngine->getOptDouble("alpha-cn-calibrate"), m_pEngine->getOptDouble("beta-cn-calibrate")));
                }
                else if (p->getChromosome() == iXChromosome) {
                    loadBuffer(pBuffer, iIndex, p->getCalibratedSmoothedLog2Ratio(m_pEngine->getOptDouble("alpha-X-cn-calibrate"), m_pEngine->getOptDouble("beta-X-cn-calibrate")));
                }
                else if (p->getChromosome() == iYChromosome) {
                    loadBuffer(pBuffer, iIndex, p->getCalibratedSmoothedLog2Ratio(m_pEngine->getOptDouble("alpha-Y-cn-calibrate"), m_pEngine->getOptDouble("beta-Y-cn-calibrate")));
                }
                else {
                    loadBuffer(pBuffer, iIndex, (float)0.0);
                }

                iCount++;
        }
        if (iCount >= iStructCount)
        {
            set.WriteBuffer(pBuffer, iIndex);
            iCount = 0;
            iIndex = 0;
        }
    }
    if (iCount > 0) {set.WriteBuffer(pBuffer, iIndex);}
    delete[] pBuffer;
}

void CNReporterCychp::writeAlgorthmDataMarkerABSignalDataSet(  affymetrix_calvin_io::GenericFileWriter* pWriter,
                                                                int iGroupIndex,
                                                                int iSetIndex,
                                                                std::vector<int>& vOffsets,
                                                                CNProbeSetArray& arProbeSets)
{
    const bool cyto2flag = m_pEngine->getOptBool("cyto2");

    affymetrix_calvin_io::DataSetWriter& set = pWriter->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    pWriter->SeekFromBeginPos(vOffsets[iSetIndex]);

    int iSize = 0;
    iSize += sizeof(unsigned int);
    if (cyto2flag) {
        iSize += sizeof(float);
        iSize += sizeof(float);
    }
    iSize += sizeof(float);
    int iStructCount = 100000;
    char* pBuffer = new char[iSize * iStructCount];
    int iIndex = 0;
    int iCount = 0;
    unsigned int uiIndex = 0;
    for (int iProbeSetIndex = 0; (iProbeSetIndex < arProbeSets.getCount()); iProbeSetIndex++)
    {
        CNProbeSet* p = arProbeSets.getAt(iProbeSetIndex);
        if (!p->isProcess()) {continue;}
        if (p->processAsVisualization())
        {
            loadBuffer(pBuffer, iIndex, uiIndex);
            if (cyto2flag) {
                loadBuffer(pBuffer, iIndex, p->getAAlleleSignal());
                loadBuffer(pBuffer, iIndex, p->getBAlleleSignal());
                loadBuffer(pBuffer, iIndex, p->getSCAR());
            } else {
                loadBuffer(pBuffer, iIndex, p->getAllelicDifference());
            }

            iCount++;
            if (iCount >= iStructCount)
            {
                set.WriteBuffer(pBuffer, iIndex);
                iCount = 0;
                iIndex = 0;
            }
        }
        uiIndex++;
    }
    if (iCount > 0) {set.WriteBuffer(pBuffer, iIndex);}
    delete[] pBuffer;
}

void CNReporterCychp::writeProbeSetsAllelePeaksDataSet(    affymetrix_calvin_io::GenericFileWriter* pWriter,
                                                                int iGroupIndex,
                                                                int iSetIndex,
                                                                std::vector<int>& vOffsets,
                                                                CNProbeSetArray& arProbeSets)
{
    affymetrix_calvin_io::DataSetWriter& set = pWriter->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    pWriter->SeekFromBeginPos(vOffsets[iSetIndex]);

    int iMaxProbeSetNameLength = getMaxProbeSetNameLength();

    int iSize = 0;
    iSize += sizeof(unsigned int);
    iSize += (iMaxProbeSetNameLength + 1);
    iSize += sizeof(unsigned char);
    iSize += sizeof(unsigned int);
    iSize += sizeof(unsigned int);
    iSize += sizeof(unsigned int);
    int iStructCount = 100000;
    char* pBuffer = new char[iSize * iStructCount];
    int iIndex = 0;
    int iCount = 0;
    AffxString str;
    for (int i = 0; (i < m_pvProbeSets->getCount()); i++)
    {
        CNProbeSet* p = m_pvProbeSets->getAt(i);
        if ((p->getAllelePeaks1() > 0) || (p->getAllelePeaks2() > 0))
        {
            str = CNReporter::prepareAscii(p->getProbeSetName(), iMaxProbeSetNameLength + 1);
            CNReporter::loadBuffer(pBuffer, iIndex, str, (iMaxProbeSetNameLength + 1));
            CNReporter::loadBuffer(pBuffer, iIndex, (unsigned char)p->getChromosome());
            CNReporter::loadBuffer(pBuffer, iIndex, (unsigned int)p->getPosition());
            CNReporter::loadBuffer(pBuffer, iIndex, p->getAllelePeaks1());
            CNReporter::loadBuffer(pBuffer, iIndex, p->getAllelePeaks2());

            iCount++;
            if (iCount >= iStructCount)
            {
                set.WriteBuffer(pBuffer, iIndex);
                iCount = 0;
                iIndex = 0;
            }
        }
    }
    if (iCount > 0) {set.WriteBuffer(pBuffer, iIndex);}
    delete[] pBuffer;
}

void CNReporterCychp::writeSegmentsDataSet(     affymetrix_calvin_io::GenericFileWriter* pWriter,
                                                int iGroupIndex,
                                                int iSetIndex,
                                                std::vector<int>& vOffsets,
                                                CNSegmentArray& arSegments,
                                                short nSegmentType)
{
    AffxString strName = CNSegment::getSegmentTypeString(nSegmentType);
    bool bMosaicism = (strName == "Mosaicism");
    bool bCN = (strName == "CN");
    affymetrix_calvin_io::DataSetWriter& set = pWriter->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    pWriter->SeekFromBeginPos(vOffsets[iSetIndex]);

    int iSize = 0;
    iSize += sizeof(unsigned int);
    iSize += sizeof(unsigned char);
    iSize += sizeof(unsigned int);
    iSize += sizeof(unsigned int);
    iSize += sizeof(unsigned int);
    iSize += sizeof(unsigned int);
    if (!bCN)
    {
        iSize += sizeof(unsigned char);
        iSize += sizeof(float);
    }
    if (bCN)
    {
        iSize += sizeof(float);
        iSize += sizeof(float);
    }
    if (bMosaicism)
    {
        iSize += sizeof(float);
        iSize += sizeof(float);
    }
    int iStructCount = 100000;
    char* pBuffer = new char[iSize * iStructCount];
    int iIndex = 0;
    int iCount = 0;
    for (int iSegmentIndex = 0; (iSegmentIndex < arSegments.getCount()); iSegmentIndex++)
    {
        CNSegment* p = arSegments.getAt(iSegmentIndex);
        if (p->getSegmentType() == nSegmentType)
        {
            m_uiSegmentID++;
            loadBuffer(pBuffer, iIndex, m_uiSegmentID);
            loadBuffer(pBuffer, iIndex, (unsigned char)p->getChromosome());
            loadBuffer(pBuffer, iIndex, (unsigned int)p->getStartPosition());
            loadBuffer(pBuffer, iIndex, (unsigned int)p->getEndPosition());
            loadBuffer(pBuffer, iIndex, (unsigned int)p->getMarkerCount());
            loadBuffer(pBuffer, iIndex, (unsigned int)p->getMeanMarkerDistance()); // In base pairs
            if (!bCN)
            {
                loadBuffer(pBuffer, iIndex, (unsigned char)p->getCall());
                loadBuffer(pBuffer, iIndex, p->getConfidence());
            }
            if (bCN)
            {
                loadBuffer(pBuffer, iIndex, (float)p->getCall());
                loadBuffer(pBuffer, iIndex, p->getConfidence());
            }
            if (bMosaicism)
            {
                loadBuffer(pBuffer, iIndex, p->getCalibratedCN());
                // convert our enum.
                float tmp=MosaicClassToFloat(p->getMixture());
                loadBuffer(pBuffer, iIndex, tmp);
            }

            iCount++;
            if (iCount >= iStructCount)
            {
                set.WriteBuffer(pBuffer, iIndex);
                iCount = 0;
                iIndex = 0;
            }
        }
    }
    if (iCount > 0) {set.WriteBuffer(pBuffer, iIndex);}
    delete[] pBuffer;
}

void CNReporterCychp::writeGenInfoDataSet(    affymetrix_calvin_io::GenericFileWriter* pWriter,
                                               int iGroupIndex,
                                               int iSetIndex,
                                               std::vector<int>& vOffsets,
                                               CNProbeSetArray& arProbeSets)
{
    affymetrix_calvin_io::DataSetWriter& set = pWriter->GetDataGroupWriter(iGroupIndex).GetDataSetWriter(iSetIndex);
    pWriter->SeekFromBeginPos(vOffsets[iSetIndex]);

    int iSize = 0;
    iSize += sizeof(unsigned int);
    iSize += sizeof(unsigned char);
    iSize += sizeof(float);
    iSize += sizeof(unsigned char);
    iSize += sizeof(float);
    iSize += sizeof(float);
    iSize += sizeof(float);
    iSize += sizeof(float);
    int iStructCount = 100000;
    char* pBuffer = new char[iSize * iStructCount];
    int iIndex = 0;
    int iCount = 0;
    unsigned int uiIndex = 0;
    float fConfidenceThreshold = CNAnalysisMethod::getConfidenceThreshold(m_pEngine->getOpt("brlmmp-parameters"));

    for (int iProbeSetIndex = 0; (iProbeSetIndex < arProbeSets.getCount()); iProbeSetIndex++)
    {
        CNProbeSet* p = arProbeSets.getAt(iProbeSetIndex);
        // the three if clauses here (isProcess, SNP, confidence)
        // should sync up with if clauses in CNAnalysisMethodGenotype::calculateCallRate()

        if (!p->isProcess()) {continue;}
        if (p->processAsSNP())
        {
            loadBuffer(pBuffer, iIndex, uiIndex);
            if (p->getGenotypeConfidence() >= fConfidenceThreshold) {
                loadBuffer(pBuffer, iIndex, (unsigned char)ALLELE_NO_CALL);
            } else {
                loadBuffer(pBuffer, iIndex, (unsigned char)p->getGenotypeCallCode());
            }
            loadBuffer(pBuffer, iIndex, p->getGenotypeConfidence());
            loadBuffer(pBuffer, iIndex, (unsigned char)p->getGenotypeCallCode());
            loadBuffer(pBuffer, iIndex, p->getAAlleleSignal());
            loadBuffer(pBuffer, iIndex, p->getBAlleleSignal());
            loadBuffer(pBuffer, iIndex, p->getSignalStrengthMvA());
            loadBuffer(pBuffer, iIndex, p->getSignalContrastMvA());

            iCount++;
            if (iCount >= iStructCount)
            {
                set.WriteBuffer(pBuffer, iIndex);
                iCount = 0;
                iIndex = 0;
            }
        }
        uiIndex++;
    }
    if (iCount > 0) {set.WriteBuffer(pBuffer, iIndex);}
    delete[] pBuffer;
}
