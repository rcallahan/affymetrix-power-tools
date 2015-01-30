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
 * @file CNReporterCnchp.cpp
 * @brief This file contains the CNReporterCnchp class members.
 */

//
#include "copynumber/CNReporterCnchp.h"
//
#include "copynumber/Annotation.h"
#include "copynumber/CNAnalysisMethod.h"
//
#include "calvin_files/data/src/CHPMultiDataData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileBufferWriter.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
#include "chipstream/BioTypes.h"
#include "util/AffxStatistics.h"
#include "util/BaseEngine.h"
#include "util/Fs.h"

CNReporterCnchp::CNReporterCnchp()
{
    excludedOptionNames.insert("file-name-prefix");
    excludedOptionNames.insert("file-name-suffix");
    excludedOptionNames.insert("file-name-ext");
    excludedOptionNames.insert("waviness-seg-count-loss");
    excludedOptionNames.insert("waviness-seg-count-gain");
    excludedOptionNames.insert("waviness-seg-count");
    excludedOptionNames.insert("wave-count");
    excludedOptionNames.insert("reference-wave-count");
    excludedOptionNames.insert("cytoscan-hd");
}

/**
 * @brief Define the options used by this reporter
 * @param PgOptions& - The options to define into
 */
void CNReporterCnchp::defineOptions(BaseEngine &e)
{
    e.defineOption("","cnchp-output", PgOpt::BOOL_OPT,
                    "Report CNCHP files",
                    "true");
}

/**
 * @brief Validate the options
 * @param PgOPtions& - The options to validate
 */
void CNReporterCnchp::checkOptions(BaseEngine &e)
{
}

/**
 * @brief Run the reporter
 */
void CNReporterCnchp::run()
{
    isSetup();
    if (!m_pEngine->getOptBool("cnchp-output")) {return;}
    Verbose::out(1, "CNReporterCnchp::report(...) start");

    int iMaxChromosome = 0;
    int iMaxProbeSetNameLength = 0;
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
        iMaxChromosome = Max(iMaxChromosome, (int)pobjProbeSet->getChromosome());
    }
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
                if (Fs::basename(str).find(strCelFileName) == 0)
                {
                    strCelFileName = str;
                    break;
                }
            }
        }
    }

    AffxString strExperimentName = getExperiment()->getExperimentName().toLowerCase();
    if (strExperimentName.endsWith(".cel"))
    {
        strExperimentName = getExperiment()->getExperimentName().substring(0, (unsigned int)strExperimentName.length() - 4);
    }
    else {strExperimentName = getExperiment()->getExperimentName();}

    AffxString strFileName;
    AffxString strAnalysisName = m_pEngine->getOpt("set-analysis-name");
    if (strAnalysisName == "CN5")
    {
        strFileName = Fs::join(m_pEngine->getOpt("out-dir"),strExperimentName + ".CN5.CNCHP");
    }
    else
    {
        strFileName = Fs::join(m_pEngine->getOpt("out-dir"),strExperimentName + "." + strAnalysisName + ".CN5.CNCHP");
    }

    if ((m_pEngine->isOptDefined("result-files")) && (m_pEngine->getOptVector("result-files").size() > 0))
    {
        AffxString strPotentialName = m_pEngine->getOptVector("result-files")[getExperiment()->getIndex()];
        if (strPotentialName != "") { strFileName = strPotentialName; }
    }

    Verbose::out(1, "Writing file: " + strFileName);
    try
    {
        affymetrix_calvin_parameter::ParameterNameValueTypeList params;
        affymetrix_calvin_parameter::ParameterNameValueType param;

        affymetrix_calvin_io::CHPMultiDataData objChpData(strFileName);

        // Copy over headers from the cel file.
        try
        {
        affymetrix_fusion_io::FusionCELData objCelFileData;
        objCelFileData.SetFileName(strCelFileName.c_str());
        objCelFileData.ReadHeader();
        if (objCelFileData.GetGenericData() != NULL)
        {
            objChpData.GetFileHeader()->GetGenericDataHdr()->AddParent(*objCelFileData.GetGenericData()->Header().GetGenericDataHdr());
        }
        objCelFileData.Close();
        } catch(...) {}

        objChpData.SetAlgName(L"CN5");
        objChpData.SetAlgVersion(L"5.0.0");
        objChpData.SetArrayType(StringUtils::ConvertMBSToWCS(getArrayName()));

        affymetrix_calvin_io::GenericDataHeader* pobjGenericDataHeader = objChpData.GetFileHeader()->GetGenericDataHdr();
        pobjGenericDataHeader->SetFileId(affxutil::Guid::GenerateNewGuid());

        pobjGenericDataHeader->SetFileCreationTime(DateTime::GetCurrentDateTime().ToString());
        param.SetName(L"affymetrix-array-type");
        param.SetValueText(StringUtils::ConvertMBSToWCS(getArrayName()));
        pobjGenericDataHeader->AddNameValParam(param);

        // This one is the one GTC uses!
        pobjGenericDataHeader->SetFileCreationTime(DateTime::GetCurrentDateTime().ToString());
        param.SetName(L"affymetrix-algorithm-param-ArraySet");
        param.SetValueText(StringUtils::ConvertMBSToWCS(getArrayName()));
        pobjGenericDataHeader->AddNameValParam(param);

        param.SetName(L"program-name");
        param.SetValueText(StringUtils::ConvertMBSToWCS(m_pEngine->getOpt("program-name")));
        pobjGenericDataHeader->AddNameValParam(param);

        param.SetName(L"program-version");
        param.SetValueText(StringUtils::ConvertMBSToWCS(m_pEngine->getOpt("program-version")));
        pobjGenericDataHeader->AddNameValParam(param);

        param.SetName(L"program-company");
        param.SetValueText(StringUtils::ConvertMBSToWCS("Affymetrix"));
        pobjGenericDataHeader->AddNameValParam(param);

        param.SetName(L"create_date");
        param.SetValueText(StringUtils::ConvertMBSToWCS(Util::getTimeStamp()));
        pobjGenericDataHeader->AddNameValParam(param);

        vector<affymetrix_calvin_io::ColumnInfo> cols;
        cols.push_back(affymetrix_calvin_io::FloatColumn(L"CNState"));
        cols.push_back(affymetrix_calvin_io::FloatColumn(L"Log2Ratio"));
        if (m_pEngine->getOptInt("gaussian-smooth-exp") >= 0)
        {
            cols.push_back(affymetrix_calvin_io::FloatColumn(L"SmoothSignal"));
        }
        cols.push_back(affymetrix_calvin_io::FloatColumn(L"LOH"));
        cols.push_back(affymetrix_calvin_io::FloatColumn(L"Allele Difference"));
        objChpData.SetEntryCount(    affymetrix_calvin_io::CopyNumberMultiDataType,
                                                (int32_t)getProbeSets()->getCNProcessCount(),
                                                iMaxProbeSetNameLength + 1,
                                                cols);

        param.SetName(L"cnchp-version");
        param.SetValueInt32(3);
        params.push_back(param);

        // dump options
        vector<string> optionNames;
        m_pEngine->getOptionNames(optionNames,1);
        std::vector<PgOpt::PgOptType> optionTypes;
        m_pEngine->getOptionTypes(optionTypes,1);
        for(int i=0; i< optionNames.size(); i++) {
            std::string name = optionNames[i];
            if (isExcluded(name)) {
                continue;
            }

            std::vector<std::string> vals = m_pEngine->getOptVector(name,1);
            if (vals.size() > 1)
            {
                for(int j=0; j< vals.size(); j++)
                {
                    std::string val = vals[j];
                    loadParam("option-" + name + "-" + ::getInt(j+1), optionTypes[i], val, param);
                    params.push_back(param);
                }
            }
            else
            {
                std::string val = m_pEngine->getOpt(name,1);
                loadParam("option-" + name, optionTypes[i], val, param);
                params.push_back(param);
            }
        }
        // dump state.
        std::vector<std::string> stateNames;
        m_pEngine->getOptionNames(stateNames);
        std::vector<PgOpt::PgOptType> stateTypes;
        m_pEngine->getOptionTypes(stateTypes);
        for(int i=0; i< stateNames.size(); i++) {
            std::string name = stateNames[i];
            if (isExcluded(name)) {
                continue;
            }

            std::vector<std::string> vals = m_pEngine->getOptVector(name);
            if (vals.size() > 1)
            {
                for(int j=0; j< vals.size(); j++)
                {
                    std::string val = vals[j];
                    loadParam("state-" + name + "-" + ::getInt(j+1), stateTypes[i], val, param);
                    params.push_back(param);
                }
            }
            else
            {
                std::string val = m_pEngine->getOpt(name);
                loadParam("state-" + name, stateTypes[i], val, param);
                params.push_back(param);
            }
        }

        objChpData.AddAlgParams(params);

        vector<pair<string, string> > metaData = m_pEngine->getMetaDataDescription();
        ParameterNameValueTypeList appMetaInfoList;
        for (int i = 0; i < metaData.size(); i++) 
        {
                pair<string,string> p = metaData[i];
                param.SetName(StringUtils::ConvertMBSToWCS(p.first));
                param.SetValueText(StringUtils::ConvertMBSToWCS(p.second));
                appMetaInfoList.push_back(param);
        }
        objChpData.AddAppMetaInfo(appMetaInfoList);

        for (int iIndex = 0; (iIndex < (int)CNAnalysisMethod::getParams()->size()); iIndex++)
        {
            ParameterNameValueType param = CNAnalysisMethod::getParams()->at(iIndex);
            wstring name = param.GetName();
            if (name == L"affymetrix-algorithm-param-wave-count") {
                continue;
            }
            pobjGenericDataHeader->AddNameValParam(param);
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
                        pobjGenericDataHeader->AddNameValParam(CNAnalysisMethod::getCelFileParams()->at(uiParamIndex));
                    }
                }
            }
        }

        // Copy over headers from the snp annotation file.
        if (Annotation::getParams()->size() > 0)
        {
            for (unsigned int uiParamIndex = 0; (uiParamIndex < Annotation::getParams()->size()); uiParamIndex++)
            {
                pobjGenericDataHeader->AddNameValParam(Annotation::getParams()->at(uiParamIndex));
            }
        }

        params.clear();
        /*
        param.SetName(L"CN2Gender-ChrX-mean");
        param.SetValueFloat(m_pobjExperiment->getChrXMean());
        params.push_back(param);

        param.SetName(L"CN2Gender-ChrY-mean");
        param.SetValueFloat(m_pobjExperiment->getChrYMean());
        params.push_back(param);

        param.SetName(L"CN2Gender");
        param.SetValueFloat(m_pobjExperiment->getGenderAsInt());
        params.push_back(param);
        */

        param.SetName(L"Gender");
        param.SetValueFloat(m_pobjExperiment->getCNCallGenderAsInt());
        params.push_back(param);

        param.SetName(L"MedianAutosomeMedian");
        param.SetValueFloat(m_pobjExperiment->getMedianAutosomeMedian());
        params.push_back(param);
        param.SetName(L"MAPD");
        param.SetValueFloat(m_pobjExperiment->getMadDiffCN());
        params.push_back(param);
        param.SetName(L"iqr");
        param.SetValueFloat(m_pobjExperiment->getIqr());
        params.push_back(param);
        param.SetName(L"all_probeset_rle_mean");
        param.SetValueFloat(m_pobjExperiment->getMeanAbsRle());
        params.push_back(param);
        param.SetName(L"gc-correction-size");
        param.SetValueFloat(m_pobjExperiment->getGCCorrectionMetric());
        params.push_back(param);

        param.SetName(L"sample-median-cn-state");
        param.SetValueFloat(m_pobjExperiment->getMedianCnState());
        params.push_back(param);

        param.SetName(L"sample-hom-frequency");
        param.SetValueFloat(m_pobjExperiment->getHomFrequency());
        params.push_back(param);

        param.SetName(L"sample-het-frequency");
        param.SetValueFloat(m_pobjExperiment->getHetFrequency());
        params.push_back(param);

//        param.SetName(L"median-raw-intensity");
//        param.SetValueFloat(m_pobjExperiment->getMedianRawIntensity());
//        params.push_back(param);

        if (CNExperiment::getQCMetricColumnNames()->getCount() != m_pobjExperiment->getQCMetricColumnValues()->getCount())
        {
            Err::errAbort("Mismatch occurred between QC Metric column names count and QC Metric column values count.");
        }
        for (int iMetricIndex = 0; (iMetricIndex < CNExperiment::getQCMetricColumnNames()->getCount()); iMetricIndex++)
        {
            param.SetName(StringUtils::ConvertMBSToWCS(*CNExperiment::getQCMetricColumnNames()->getAt(iMetricIndex)));
            if ((*CNExperiment::getQCMetricColumnNames()->getAt(iMetricIndex)).endsWith("_gender"))
            {
                float f = affx::UnknownGender;
                AffxString str = *m_pobjExperiment->getQCMetricColumnValues()->getAt(iMetricIndex);
                str = str.toLowerCase();
                if (str == "male" ) {f = affx::Male;} else if (str == "female") {f = affx::Female;}
                param.SetName(StringUtils::ConvertMBSToWCS(*CNExperiment::getQCMetricColumnNames()->getAt(iMetricIndex)));
                param.SetValueFloat(f);
                params.push_back(param);
            }
            else
            {
                param.SetName(StringUtils::ConvertMBSToWCS(*CNExperiment::getQCMetricColumnNames()->getAt(iMetricIndex)));
                param.SetValueFloat((float)::getDouble(*m_pobjExperiment->getQCMetricColumnValues()->getAt(iMetricIndex)));
                params.push_back(param);
            }
        }
        for (int iChromosome = 0; iChromosome <= iMaxChromosome; iChromosome++)
        {
            float fMinLog2Ratio = 1000000;
            float fMaxLog2Ratio = -1000000;
            int iProbeSetCount = 0;
            for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
                if (pobjProbeSet->getChromosome() == iChromosome)
                {
                    fMinLog2Ratio = Min(fMinLog2Ratio, pobjProbeSet->getLog2Ratio());
                    fMaxLog2Ratio = Max(fMaxLog2Ratio, pobjProbeSet->getLog2Ratio());
                    iProbeSetCount++;
                }
            }
            if (iProbeSetCount > 0)
            {
                AffxString str = "chrom_" + ::getInt(iChromosome) + "_MinSignal";
                param.SetName(StringUtils::ConvertMBSToWCS(str));
                param.SetValueFloat(fMinLog2Ratio);
                params.push_back(param);
                str = "chrom_" + ::getInt(iChromosome) + "_MaxSignal";
                param.SetName(StringUtils::ConvertMBSToWCS(str));
                param.SetValueFloat(fMaxLog2Ratio);
                params.push_back(param);
            }
        }

        // CN waviness parameters
        int segCountLoss, segCountGain;
        float sd;
        wavinessSegCounts(segCountLoss, segCountGain, sd);

        param.SetName(L"waviness-sd");
        param.SetValueFloat(sd);
        params.push_back(param);

        objChpData.AddSummaryParams(params);
        objChpData.GetGenericData().Header().GetGenericDataHdr()->SetFileTypeId(CHP_MULTI_DATA_TYPE);

        affymetrix_calvin_io::DataSetHeader* pobjDataSetHeader = objChpData.GetDataSetHeader(affymetrix_calvin_io::CopyNumberMultiDataType);
        for (int iChromosome = 0; iChromosome <= iMaxChromosome; iChromosome++)
        {
            int iStartIndex = -1;
            int iProbeSetCount = 0;
            for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
                if (pobjProbeSet->getChromosome() == iChromosome)
                {
                    if (iStartIndex == -1) {iStartIndex = iIndex;}
                    iProbeSetCount++;
                }
            }
            if (iProbeSetCount > 0)
            {
                std::wstring wstr = StringUtils::ConvertMBSToWCS(::getInt(iChromosome));
                param.SetName(wstr + L":start");
                param.SetValueInt32(iStartIndex);
                pobjDataSetHeader->AddNameValParam(param);
                param.SetName(wstr + L":count");
                param.SetValueInt32(iProbeSetCount);
                pobjDataSetHeader->AddNameValParam(param);
                param.SetName(wstr + L":display");
                if (iChromosome == m_pEngine->getOptInt("xChromosome")) {param.SetValueAscii("X");}
                else if (iChromosome == m_pEngine->getOptInt("yChromosome")) {param.SetValueAscii("Y");}
                else {param.SetValueAscii(::getInt(iChromosome));}
                pobjDataSetHeader->AddNameValParam(param);
            }
        }

        affymetrix_calvin_io::CHPMultiDataFileWriter *objChpWriter = new affymetrix_calvin_io::CHPMultiDataFileWriter(objChpData);
        delete objChpWriter;

        affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData objCopyNumberData;

        param.SetName(L"CNState");
        param.SetValueFloat(0);
        objCopyNumberData.metrics.push_back(param);
        param.SetName(L"Log2Ratio");
        param.SetValueFloat(0);
        objCopyNumberData.metrics.push_back(param);
        if (m_pEngine->getOptInt("gaussian-smooth-exp") >= 0)
        {
            param.SetName(L"SmoothSignal");
            param.SetValueFloat(0);
            objCopyNumberData.metrics.push_back(param);
        }
        param.SetName(L"LOH");
        param.SetValueFloat(0);
        objCopyNumberData.metrics.push_back(param);
        param.SetName(L"Allele Difference");
        param.SetValueFloat(0);
        objCopyNumberData.metrics.push_back(param);

        affymetrix_calvin_io::CHPMultiDataFileBufferWriter bufferWriter;
        std::vector<affymetrix_calvin_io::MultiDataType> dataTypes;
        std::map<affymetrix_calvin_io::MultiDataType, int> maxLen;
        dataTypes.push_back(affymetrix_calvin_io::CopyNumberMultiDataType);
        maxLen[affymetrix_calvin_io::CopyNumberMultiDataType] = iMaxProbeSetNameLength + 1;
        std::vector<std::string> chpFiles;
        chpFiles.push_back(strFileName);
        bufferWriter.Initialize(&chpFiles, dataTypes, maxLen);

        //objChpWriter.SeekToDataSet(affymetrix_calvin_io::CopyNumberMultiDataType);
        for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
        {
            CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
            if (!pobjProbeSet->isProcess()) {continue;}
            AffxString strProbeSetName = pobjProbeSet->getProbeSetName();
            //int iLength = (int)strProbeSetName.length();
            objCopyNumberData.name = strProbeSetName;
            objCopyNumberData.chr = pobjProbeSet->getChromosome();
            objCopyNumberData.position = pobjProbeSet->getPosition();
            objCopyNumberData.metrics[0].SetValueFloat((float)pobjProbeSet->getCNState());
            objCopyNumberData.metrics[1].SetValueFloat(pobjProbeSet->getLog2Ratio());
            if (m_pEngine->getOptInt("gaussian-smooth-exp") >= 0)
            {
                if (m_pEngine->getOptInt("gaussian-smooth-exp") == 1)
                {
                    int iXChromosome = m_pEngine->getOptInt("xChromosome");
                    int iYChromosome = m_pEngine->getOptInt("yChromosome");
                    if (pobjProbeSet->getChromosome() < iXChromosome) {
                        objCopyNumberData.metrics[2].SetValueFloat(pobjProbeSet->getCalibratedSmoothedLog2Ratio(m_pEngine->getOptDouble("alpha-cn-calibrate"), m_pEngine->getOptDouble("beta-cn-calibrate")));
                    }
                    else if (pobjProbeSet->getChromosome() == iXChromosome) {
                        objCopyNumberData.metrics[2].SetValueFloat(pobjProbeSet->getCalibratedSmoothedLog2Ratio(m_pEngine->getOptDouble("alpha-X-cn-calibrate"), m_pEngine->getOptDouble("beta-X-cn-calibrate")));
                    }
                    else if (pobjProbeSet->getChromosome() == iYChromosome) {
                        objCopyNumberData.metrics[2].SetValueFloat(pobjProbeSet->getCalibratedSmoothedLog2Ratio(m_pEngine->getOptDouble("alpha-Y-cn-calibrate"), m_pEngine->getOptDouble("beta-Y-cn-calibrate")));
                    }
                    else {
                        objCopyNumberData.metrics[2].SetValueFloat(0.0);
                    }
                }
                else
                {
                    objCopyNumberData.metrics[2].SetValueFloat(pobjProbeSet->getSmoothedLog2Ratio());
                }
                objCopyNumberData.metrics[3].SetValueFloat(pobjProbeSet->getLoh());
                objCopyNumberData.metrics[4].SetValueFloat(pobjProbeSet->getAllelicDifference());
            }
            else
            {
                objCopyNumberData.metrics[2].SetValueFloat(pobjProbeSet->getLoh());
                objCopyNumberData.metrics[3].SetValueFloat(pobjProbeSet->getAllelicDifference());
            }
            //objChpWriter.WriteEntry(objCopyNumberData);
            bufferWriter.WriteMultiDataCopyNumberEntry(affymetrix_calvin_io::CopyNumberMultiDataType, 0, objCopyNumberData);
        }
        bufferWriter.FlushBuffer();

        // Convert from calvin to text format?
        if (m_pEngine->getOptBool("text-output"))
        {
            CalvinToText converter;
            converter.run(strFileName, strFileName + ".txt", true, true, true);
        }
    } catch(...) {throw(Except("CNCHP write failed for file: " + strFileName));}
    Verbose::out(1, "CNReporterCnchp::report(...) end");
}

bool CNReporterCnchp::isExcluded(const string& name) {
    return excludedOptionNames.count(name);
}
