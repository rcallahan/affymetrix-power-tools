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
 * @file   QuantMethodExprCCCHPReport.cpp
 * @author David Le
 * @date   Mon May 15 12:09:42 2006
 * 
 * @brief  Class for reporting results of quantification methods.
 */

//
#include "chipstream/QuantMethodExprCCCHPReport.h"
//
#include "calvin_files/data/src/CHPQuantificationData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/parsers/src/CHPQuantificationFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/GenericDataHeaderUpdater.h"
#include "portability/affy-base-types.h"
#include "util/Fs.h"
//
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

using namespace std;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_io;
using namespace affymetrix_fusion_io;

/** Constructor. */
QuantMethodExprCCCHPReport::QuantMethodExprCCCHPReport(AnalysisInfo& chpInfo,
                                                       const std::string& prefix,
                                                       const std::string algName)
{
    m_Info = chpInfo;
    m_AlgName = algName;
    m_Prefix = prefix;
    m_CurrentProbeSetCount = 0;
}

/** Destructor. */
QuantMethodExprCCCHPReport::~QuantMethodExprCCCHPReport() {
}

/**
 * Register a Chip Summary interface instance to pull metrics from
 */
void QuantMethodExprCCCHPReport::registerChipSummary(ChipSummary *summary) {
    m_ChipSummaries.push_back(summary);
}

/** Do a check that the original id we wrote out and the one that
    actually came with the probeset group match. */
void QuantMethodExprCCCHPReport::checkCurrentId(ProbeSetGroup &psGroup) {
    if (psGroup.name != NULL) {
        if (!Util::sameString(psGroup.name, m_Info.m_ProbesetNames[m_CurrentProbeSetCount])) {
            Err::errAbort("QuantMethodExprCCCHPReport::checkCurrentId() - Expecting to get name: '" + 
                          ToStr(m_Info.m_ProbesetNames[m_CurrentProbeSetCount]) + "' but got name: '" + ToStr(psGroup.name) +
                          "' instead.");
        }
    }
    else {
        Err::errAbort("QuantMethodExprCCCHPReport::checkCurrentId() - probeset with no name?'");
    }
}


/**
 * Set the CHP filenames to use for output. This overrides any prefix and allows
 * the caller to place the chp files anyplace they are desired.
 * @param fileNames - Vector of fileNames with full path to output
 * file, one for each cel file.
 */
void QuantMethodExprCCCHPReport::setChpFileNames(std::vector<std::string> &fileNames) {
    m_CHPFileNames = fileNames;
}

/** Swap out the .cel for .chp */
void QuantMethodExprCCCHPReport::setupFileNames(const IntensityMart &iMart)
{
    m_CELFileNames = iMart.getCelFileNames();
    if (m_CHPFileNames.size() == 0) {
        m_CHPFileNames = Fs::changeDirAndExt(m_Prefix,m_CELFileNames,"."+m_Info.m_AlgName+".chp");
    }
    if (m_CHPFileNames.size() != iMart.getCelFileNames().size()) 
        Err::errAbort("Must be same number of output CHP files (" + 
                      ToStr(m_CHPFileNames.size()) + ") as input CEL files (" + ToStr(iMart.getCelFileNames().size()));
}

void QuantMethodExprCCCHPReport::removeTmpChps() {
    // Remove any existing TMP CHP file
    std::vector<std::string> filesToRemove;
    for (int chip=0; chip<m_CHPFileNames.size(); chip++) {
        std::string tmpChpName = m_CHPFileNames[chip] + ".tmp";
        filesToRemove.push_back(tmpChpName);
    }
    for ( int i=0; i < filesToRemove.size();  i++ ) { 
        Fs::rm(filesToRemove[i], false);
    } 
}

void QuantMethodExprCCCHPReport::removeAllChps() {
    // Remove any existing CHP file
    std::vector<std::string> filesToRemove;
    for (int chip=0; chip<m_CHPFileNames.size(); chip++) {
        std::string chpName = m_CHPFileNames[chip];
        filesToRemove.push_back(chpName);
    }
    for ( int i=0; i < filesToRemove.size();  i++ ) { 
        Fs::rm(filesToRemove[i], false);
    }
    removeTmpChps();
}

/** 
 * Get set up for a run of reporting probesets. Often used to open file
 * streams and print headers to files etc.
 * 
 * @param qMethod - Quantification method to be used.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodExprCCCHPReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart) 
{
    QuantExprMethod *eMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
    if (eMethod == NULL) { Err::errAbort("Can only use a QuantMethodExprReport with a QuantExprMethod."); }
    
    setupFileNames(iMart);
    int nfiles = m_CHPFileNames.size();
    
    // Make sure our output directory exists.
    if (!Fs::isWriteableDir(m_Prefix.c_str()) &&
        (Fs::mkdirPath(m_Prefix, false) != APT_OK)) {
        APT_ERR_ABORT("Can't make or write to directory: " + m_Prefix);
    }
    
    removeAllChps();
    
    // Get CEL file GUIDs
    ///@todo This be computed by the engine and passed in via AnalysisInfo
    m_celGuids.resize(nfiles);
    std::string tmp_unc_name;
    for (int chip=0; chip<nfiles; chip++) {
        FusionCELData cel;
        try {
            tmp_unc_name=Fs::convertToUncPath(m_CELFileNames[chip]);
            cel.SetFileName(tmp_unc_name.c_str());
            if (!cel.ReadHeader()) {
                Err::errAbort("Unable to read CEL file: "+FS_QUOTE_PATH(tmp_unc_name));
            }
            GenericData *gdata = cel.GetGenericData();
            if (gdata != NULL) {
                m_celGuids[chip] = gdata->Header().GetGenericDataHdr()->GetFileId();
            }
            cel.Close();
        }
        catch (...) {
            Err::errAbort("Unable to read CEL file " + tmp_unc_name);
        }
    }

    int maxProbeSetNameLength = 0;
    for (int i=0; i<m_Info.m_ProbesetNames.size(); i++) {
        int len = (int)strlen(m_Info.m_ProbesetNames.at(i));
        if (m_Info.m_ProbesetDisplayNames.size() > 0 && m_Info.m_ProbesetDisplayNames.at(i) != NULL)
            len = (int)strlen(m_Info.m_ProbesetDisplayNames.at(i));
        maxProbeSetNameLength = Max(maxProbeSetNameLength, len);
    }

    // Prepare headers for all CHP files.
    wstring algName = StringUtils::ConvertMBSToWCS(m_Info.m_AlgName);
    wstring algVersion = StringUtils::ConvertMBSToWCS(m_Info.m_AlgVersion);
    wstring chipType = StringUtils::ConvertMBSToWCS(m_Info.m_ChipType);

    // For each chip, precreate all probeset signal entries (default to 0.0).
    Verbose::out(1,"QuantMethodExprCCCHPReport: Creating temporary files for CHP output");
    for (int chip=0; chip<nfiles; chip++) {
        try {
            ParameterNameValueType param;

            // Create tmp chp file
            std::string tmp_chp_name=m_CHPFileNames[chip] + ".tmp";
            CHPQuantificationData *data = new CHPQuantificationData(tmp_chp_name);
            m_TmpChpFiles.push_back(tmp_chp_name);

            // set parent header
            FusionCELData cel;
            try {
                tmp_unc_name=Fs::convertToUncPath(m_CELFileNames[chip]);
                cel.SetFileName(tmp_unc_name.c_str());
                if (!cel.ReadHeader()) {
                  Err::errAbort("Unable to read CEL file: "+FS_QUOTE_PATH(tmp_unc_name));
                }
                GenericData *gdata = cel.GetGenericData();
                if (gdata != NULL) {
                    data->GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr());
                }
                cel.Close();
            }
            catch (...) {
              Err::errAbort("Unable to read CEL file: "+FS_QUOTE_PATH(tmp_unc_name));
            }

            data->SetEntryCount(m_Info.m_NumProbeSets, maxProbeSetNameLength); 
            data->SetAlgName(algName);
            data->SetAlgVersion(algVersion);
            data->SetArrayType(chipType);

            param.SetName(L"program-name");
            param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ProgramName));
            data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
            param.SetName(L"program-version");
            param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ProgramVersion));
            data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
            param.SetName(L"program-company");
            param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ProgramCompany));
            data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);

            // Add algorithm parameters to list.
            ParameterNameValueTypeList paramList;
            assert(m_Info.m_ParamNames.size() == m_Info.m_ParamValues.size());
            for (int i=0; i<m_Info.m_ParamNames.size(); i++) {
                if (m_Info.m_ParamValues[i].length() > 0) {
                    param.SetName(StringUtils::ConvertMBSToWCS(m_Info.m_ParamNames[i]));
                    param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ParamValues[i]));
                    paramList.push_back(param);
                }
            }

            // Add list of all CEL GUIDs in batch
            ///@todo should this be computed by the engine and passed in via AnalysisInfo?
            string prefix = "apt-opt-";
            for (int chip=0; chip<m_CHPFileNames.size(); chip++) {
                if (m_celGuids[chip].empty() == false) {
                    string paramName = prefix + "cel-guid-" + ToStr(chip+1);
                    param.SetName(StringUtils::ConvertMBSToWCS(paramName));
                    param.SetValueText(StringUtils::ConvertMBSToWCS(m_celGuids[chip]));
                    paramList.push_back(param);
                }
            }
            data->AddAlgParams(paramList);

            // Add the run report parameters to the list
            ParameterNameValueTypeList summaryParamList;
            std::string blankStr(256, ' ');
            for (int source=0; source<m_ChipSummaries.size(); source++) {
                ChipSummary::metricDefVec_t metricDefs = m_ChipSummaries[source]->getMetricDefs();
                for (int i = 0; i < metricDefs.size(); i++) {
                    param.SetName(StringUtils::ConvertMBSToWCS(metricDefs[i].m_name));
                    if (metricDefs[i].m_type == ChipSummary::Metric::Double) {
                        param.SetValueFloat(-1.0);
                    } 
                    else if (metricDefs[i].m_type == ChipSummary::Metric::Integer) {
                        param.SetValueInt32(-1);
                    } 
                    else if (metricDefs[i].m_type == ChipSummary::Metric::String) {
                        param.SetValueAscii(blankStr);
                    } 
                    else {
                        Err::errAbort("QuantMethodExprCCCHPReport: Unable to handle unknown type: " + 
                                      ToStr(metricDefs[i].m_type) );
                    }
                    summaryParamList.push_back(param);
                }
            }
            data->AddSummaryParams(summaryParamList);
			
            ProbeSetQuantificationData entry;
            CHPQuantificationFileWriter writer(*data);
            writer.SeekToDataSet();        // seek to data table location
            for (int index=0; index<m_Info.m_ProbesetNames.size(); index++) {
                if (m_Info.m_ProbesetDisplayNames.size() > 0 && m_Info.m_ProbesetDisplayNames[index] != NULL)
                    entry.name = m_Info.m_ProbesetDisplayNames[index];
                else
                    entry.name = m_Info.m_ProbesetNames[index];
                entry.quantification = 0.0f;
                writer.WriteEntry(entry);
            }
            
            delete data;
        }
        catch (...) {
            Err::errAbort("QuantMethodExprCHPReport::prepare() - Unable to write header and/or precreate signal entries to file: " + m_CHPFileNames[chip] + ".tmp");
        }
    }
    
    // initialize expression signal buffer writer
    m_ExpressionQuantificationBufferWriter.Initialize(&m_TmpChpFiles);

    return true;
}

/** 
 * After every probeset computation this function is called an is an opportunity
 * to query the quantification method for results, residuals, etc.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodExprCCCHPReport::report(ProbeSetGroup &psGroup, 
                                        QuantMethod &qMethod, 
                                        const IntensityMart &iMart, 
                                        std::vector<ChipStream *> &iTrans, 
                                        PmAdjuster &pmAdjust) 
{
    QuantExprMethod *eMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
    if (eMethod == NULL) { Err::errAbort("Can only use a QuantMethodExprReport with a QuantExprMethod."); }
    assert(!psGroup.probeSets.empty());
    checkCurrentId(psGroup);
    // We only report expression probe sets.
    if (psGroup.probeSets[0]->psType != ProbeSet::Expression) { return false; }

    // Write signal entry to buffer writer.
    for (int chip=0; chip<eMethod->getNumTargets(); chip++) {
        m_ExpressionQuantificationBufferWriter.WriteQuantificationEntry(chip, eMethod->getSignalEstimate(chip));
    }

    m_GoodProbesets.push_back(m_CurrentProbeSetCount);
    m_CurrentProbeSetCount++;
    return true;
}

/** 
 * If a probeset fails to compute for whatever reason then this method is
 * called rather than the normal report call above. By default does nothing.
 * 
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 * @param layout - Where the probesets, probes, etc are on the chip.
 * 
 * @return true if success, false otherwise.
 */
bool QuantMethodExprCCCHPReport::reportFailure(ProbeSetGroup &psGroup,
                                               QuantMethod &qMethod, 
                                               const IntensityMart &iMart, 
                                               std::vector<ChipStream *> &iTrans, 
                                               PmAdjuster &pmAdjust)
{
    QuantExprMethod *eMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
    if (eMethod == NULL) { Err::errAbort("Can only use a QuantMethodExprReport with QuantExprMethod."); }
    assert(!psGroup.probeSets.empty());
    checkCurrentId(psGroup);
    m_CurrentProbeSetCount++;

    // We only report expression probe sets.
    if (psGroup.probeSets[0]->psType != ProbeSet::Expression) { return false; }

    // Write a blank entry in order to be consistent with psGroups list.
    for (int chip=0; chip<m_CHPFileNames.size(); chip++) {
        m_ExpressionQuantificationBufferWriter.WriteQuantificationEntry(chip, 0.0f);
    }
    return true;
}

/** 
 * No more probesets will be processed, this is a chance to finish outputting
 * results and clean up.
 * @param qMethod - Quantification method that was used.
 * @return true if success, false otherwise.
 */
bool QuantMethodExprCCCHPReport::finish(QuantMethod &qMethod) 
{
    // Sanity to check we saw all the probe sets we were expecting.
    if (m_CurrentProbeSetCount != m_Info.m_NumProbeSets) {
        Err::errAbort("QuantMethodExprCCCHPReport::finish() - Expecting: " + ToStr(m_Info.m_NumProbeSets) +
            " but got: " + ToStr(m_CurrentProbeSetCount) + ". Command Console CHP file will be corrupt.");
    }

    // Flush remaining signal entries in the buffer.
    m_ExpressionQuantificationBufferWriter.FlushBuffer();

    // Rewrite CHP files to get chip summary entires
    Verbose::out(1,"Creating final files for CHP output");
    Verbose::progressBegin(1, ToStr("Finalizing Expression CHP Files"), 
                           m_CHPFileNames.size(), 1, m_CHPFileNames.size());
    try {
        for (unsigned int chip = 0; chip < m_CHPFileNames.size(); chip++) {
            // open up tmp chp file to pull results from
            GenericData data;
            GenericFileReader reader;
            std::string filename = m_CHPFileNames[chip]+".tmp";
            reader.SetFilename(filename);
            reader.ReadHeader(data);

            GenericDataHeader* hdr = data.Header().GetGenericDataHdr();
            GenericDataHeader updateHdr;
            for (int source = 0; source < m_ChipSummaries.size(); source++) {
                ChipSummary::metricDefVec_t metricDefs = m_ChipSummaries[source]->getMetricDefs();
                for (int i = 0; i < metricDefs.size(); i++) {
                    ChipSummary::Metric metric;
                    if (!m_ChipSummaries[source]->getMetric(chip, metricDefs[i].m_name, metric)) {
                        Err::errAbort("QuantMethodExprCCCHPReport: metric '" + metricDefs[i].m_name +
                                      "' was not found");
                    }
                    std::wstring mName(CHIP_SUMMARY_PARAMETER_NAME_PREFIX);
                    mName += StringUtils::ConvertMBSToWCS(metric.m_Name);
                    ParameterNameValueType param;
                    if (hdr->FindNameValParam(mName, param) == false) {
                        Err::errAbort("QuantMethodExprCCCHPReport: metric name '" + StringUtils::ConvertWCSToMBS(mName) +
                                      "' could not be found in the header of " + filename);
                    }

                    switch (param.GetParameterType()) {
                    case ParameterNameValueType::Int8Type:
                        param.SetValueInt8((int8_t)metric.m_Integer);
                        break;
                    
                    case ParameterNameValueType::UInt8Type:
                        param.SetValueUInt8((u_int8_t)metric.m_Integer);
                        break;
                    
                    case ParameterNameValueType::Int16Type:
                        param.SetValueInt16((int16_t)metric.m_Integer);
                        break;
                    
                    case ParameterNameValueType::UInt16Type:
                        param.SetValueUInt16((u_int16_t)metric.m_Integer);
                        break;
                    
                    case ParameterNameValueType::Int32Type:
                        param.SetValueInt32((int32_t)metric.m_Integer);
                        break;
                    
                    case ParameterNameValueType::UInt32Type:
                        param.SetValueUInt32((u_int32_t)metric.m_Integer);
                        break;
                
                    case ParameterNameValueType::FloatType:
                        param.SetValueFloat((float)metric.m_Double);
                        break;
                
                    case ParameterNameValueType::TextType:
                        param.SetValueText(StringUtils::ConvertMBSToWCS(metric.m_String), (int) metric.m_String.length());
                        break;
                
                    case ParameterNameValueType::AsciiType:
                        if (metric.m_String.size() > 256) {
                            Err::errAbort("QuantMethodExprCCCHPReport: string header parameter too long, name = '" +
                                          metric.m_Name + "', value = '" + metric.m_String + "'");
                        }
                        param.SetValueAscii(metric.m_String, (int) metric.m_String.length());
                        break;

                    default:
                        Err::errAbort("QuantMethodExprCCCHPReport: unknown header parameter type found in file " +
                                      filename);
                    }
                    updateHdr.AddNameValParam(param);
                }
            }
            std::ofstream os;
            Fs::aptOpen(os, filename, std::ios::out|std::ios::binary|std::ios::in);
            if (!os) {
                Err::errAbort("QuantMethodExprCCCHPReport: file " + filename +
                              " could not be opened for writing");
            }
            GenericDataHeaderUpdater updater;
            updater.Update(os, updateHdr, *hdr);
            os.close();

            Verbose::progressStep(1);
        }
    } catch (...) {
        removeAllChps();
        Err::errAbort("Error in creating final CHP output.");
    }
    Verbose::progressEnd(1, ToStr("Done."));

    // Remove .tmp extension
    for (unsigned int i = 0; i < m_CHPFileNames.size(); i++) {
        std::string from = m_CHPFileNames[i] + ".tmp";
        std::string to = m_CHPFileNames[i];
        if (!Fs::fileRename(from.c_str(),to.c_str())) {
            removeAllChps();
            Err::errAbort("Unable to rename '" + from + "' to '" + to + "'");
        }
    }
    removeTmpChps();

    return true;
}
