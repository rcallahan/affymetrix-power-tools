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
 * @file   QuantMethodMultiDataCCCHPReport.cpp
 * @author David Le
 * @date   Wed Mar 15 23:01:47 2006
 *
 * @brief Reporter for genotyping probe sets that outputs CC CHP files. There is
 * one CC CHP file produced for each cel file and the order of the CC CHP file is
 * guaranteed to be the same as the cdf file. It is very important that the
 * number and order of the probeset results is the exact same as the CDF
 * file. There have been a number of high impact bugs on this front and it is
 * important to be very careful.
 */

//
#include "chipstream/QuantMethodMultiDataCCCHPReport.h"
//
#include "chipstream/BioTypes.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/QuantBRLMM.h"
//
#include "calvin_files/data/src/CHPMultiDataData.h"
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/parsers/src/CHPMultiDataFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/GenericDataHeaderUpdater.h"
#include "file/CHPFileData.h"
#include "portability/affy-base-types.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cstring>
#include <ctime>
#include <string>


using namespace std;
using namespace affymetrix_calvin_data;
using namespace affymetrix_calvin_io;
using namespace affymetrix_fusion_io;
using namespace affx;

/** Swap out the .cel for .chp */
void QuantMethodMultiDataCCCHPReport::setupFileNames(const IntensityMart &iMart)
{
    m_CELFileNames = iMart.getCelFileNames();
    if (m_CHPFileNames.size() == 0) {
        m_CHPFileNames = Fs::changeDirAndExt(m_Prefix,m_CELFileNames,"."+m_Info.m_AlgName+".chp");
    }
    if (m_CHPFileNames.size() != iMart.getCelFileNames().size()) 
        Err::errAbort("Must be same number of output CHP files (" + 
                      ToStr(m_CHPFileNames.size()) + ") as input CEL files (" + ToStr(iMart.getCelFileNames().size()));
}

void QuantMethodMultiDataCCCHPReport::removeTmpChps() {
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

void QuantMethodMultiDataCCCHPReport::removeAllChps() {
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
 *
 * @return true if success, false otherwise.
 */
bool QuantMethodMultiDataCCCHPReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart)
{
    setupFileNames(iMart);
    size_t nfiles = m_CHPFileNames.size();

    // Make sure our output directory exists.
    if (!Fs::isWriteableDir(m_Prefix.c_str()) &&
        (Fs::mkdirPath(m_Prefix, false) != APT_OK))  {
        Err::errAbort("Can't make or write to directory: " + m_Prefix);
    }

    removeAllChps();

    // Get CEL file GUIDs
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
            if (gdata != NULL)
            {
                m_celGuids[chip] = gdata->Header().GetGenericDataHdr()->GetFileId();
            }
            cel.Close();
        }
        catch (...)
        {
          Err::errAbort("Unable to read CEL file: "+FS_QUOTE_PATH(tmp_unc_name));
        }
    }

    int maxProbeSetNameLength = m_Info.m_MaxPsNameLength;
    ProbeListPacked plp;
    if (genoTypeOnly == true) {
        nExpression = 0;
    }
    for (int i=0; i<m_Info.m_ProbesetNames.size(); i++) {
        maxProbeSetNameLength = Max(maxProbeSetNameLength,(int)strlen(m_Info.m_ProbesetNames.at(i)));
    }

        //If cdf says expression OR cdf says genotyping but probeset has only one block
        //>             Analyze as expression
        //> Else if cdf says one of many genotyping types
        //>             Analyze as genotyping
        //> Else ignore

        ///@todo find a better way to get probeset type and allele
        ///count than to directly access ProbeListPacked object in
        ///ChipLayout
        // plp = layout.m_PlFactory.getProbeListByName(m_Info.m_ProbesetNames.at(i));
//         if (plp.isNull()) {
//             APT_ERR_ABORT("psIndex: '"+ ToStr(i) +"' into ProbeListPacked list is out of range.");
//         }
//         int block_count = plp.block_cnt();        
//         ProbeSet::ProbeSetType ps_type=(ProbeSet::ProbeSetType)plp.get_type();
//         if (genoTypeOnly == false && 
//             (ps_type == ProbeSet::Expression ||
//              (ps_type == ProbeSet::GenoType && block_count != 2 && block_count != 4)
//              )
//             ) {
//             ++nExpression;
//         }
//         else if (block_count > 1) {
//             int first_allele = plp.get_blockAllele(0);
//             bool more_than_one_allele = false;
//             for (int i = 1; i < block_count; i++) {
//                 if (plp.get_blockAllele(i) != first_allele) {
//                     more_than_one_allele = true;
//                     break;
//                 }
//             }
//             if (more_than_one_allele && 
//                 (ps_type == ProbeSet::GenoType || 
//                  ps_type == ProbeSet::Marker || 
//                  ps_type == ProbeSet::MultichannelMarker)) {
//               ++nGenotype;
//             }
//         }
//    }

    // Prepare headers for all CHP files.
    wstring algName = StringUtils::ConvertMBSToWCS(m_Info.m_AlgName);
    wstring algVersion = StringUtils::ConvertMBSToWCS(m_Info.m_AlgVersion);
    wstring chipType = StringUtils::ConvertMBSToWCS(m_Info.m_ChipType);

    // Create columns for the force, signals a and b.
    vector<ColumnInfo> extraColumns;
    if (nGenotype > 0) {
        QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(&qMethod);
        if (gMethod == NULL) { Err::errAbort("Can only use a QuantMethodGTypeReport with a QuantGTypeMethod."); }

        ParameterNameValueType nv;
        string aName;
        string bName;
        gMethod->getAlleleValueNames(aName, bName);
        FloatColumn contrastCol(StringUtils::ConvertMBSToWCS(aName));
        extraColumns.push_back(contrastCol);
        FloatColumn strengthCol(StringUtils::ConvertMBSToWCS(bName));
        extraColumns.push_back(strengthCol);
        UByteColumn forcedCol(L"Forced Call");
        extraColumns.push_back(forcedCol);
    }

    // For each chip, precreate all probeset signal entries (default to 0.0).
    Verbose::out(1,"QuantMethodMultiDataCCCHPReport Creating temporary files for CHP output");
    std::string tmp_chp_name;
    for (unsigned int i=0; i<nfiles; i++)
    {
        try
        {
            ParameterNameValueType param;

            // Create tmp chp file
            std::string tmp_chp_name=m_CHPFileNames[i] + ".tmp";
            CHPMultiDataData *data = new CHPMultiDataData(tmp_chp_name);
            m_TmpChpFiles.push_back(tmp_chp_name);

            // set parent header
            FusionCELData cel;
            try {
                std::string tmp_unc_name=Fs::convertToUncPath(m_CELFileNames[i]);
                cel.SetFileName(tmp_unc_name.c_str());
                if (!cel.ReadHeader()) {
                    Err::errAbort("Unable to read CEL file: "+FS_QUOTE_PATH(tmp_unc_name));
                }
                GenericData *gdata = cel.GetGenericData();
                if (gdata != NULL)
                {
                    data->GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr());
                }
                cel.Close();
            }
            catch (...)
            {
              Err::errAbort("Unable to read CEL file: "+FS_QUOTE_PATH(tmp_unc_name));
            }
            if (nGenotype > 0) {
              int entryCount = nGenotype;
              if (nReporting > 0) {
                entryCount = nReporting;
              }
              data->SetEntryCount(GenotypeMultiDataType, entryCount, maxProbeSetNameLength, extraColumns);
            }
            if (nExpression > 0)
                data->SetEntryCount(ExpressionMultiDataType, nExpression, maxProbeSetNameLength);
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
            for (int i=0; i<m_Info.m_ParamNames.size(); i++) 
            {
                // reasons not to add this param...
                // ...no value for this name.
                if (m_Info.m_ParamValues[i].length()==0) {
                    continue;
                }
                // ...is it in the set of CHP supressed headers?
                if ((m_output_chip_view_md5sum==true) &&
                    (m_Info.m_ParamNames[i].find("chip-view-md5sum")!=std::string::npos)) {
                    continue;
                }
                // add it.
                param.SetName(StringUtils::ConvertMBSToWCS(m_Info.m_ParamNames[i]));
                param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ParamValues[i]));
                paramList.push_back(param);
            }
            
            // Add list of all CEL GUIDs in batch
            ///@todo This be computed by the engine and passed in via AnalysisInfo
            string prefix = "apt-opt-";
            for (int chip=0; chip<m_CHPFileNames.size(); chip++)
            {
                if (m_celGuids[chip].empty() == false)
                {
                    string paramName = prefix + "cel-guid-" + ToStr(chip+1);
                    param.SetName(StringUtils::ConvertMBSToWCS(paramName));
                    param.SetValueText(StringUtils::ConvertMBSToWCS(m_celGuids[chip]));
                    paramList.push_back(param);
                }
            }

            data->AddAlgParams(paramList);

            ///@todo This be computed by the engine and passed in via AnalysisInfo                                                          

            // Add client meta data information to list                                                                                                 
            ParameterNameValueTypeList appMetaInfoList;
            assert(m_Info.m_ClientInfoNames.size() == m_Info.m_ClientInfoValues.size());
            for (int i=0; i<m_Info.m_ClientInfoNames.size(); i++)
            {
                if (m_Info.m_ClientInfoValues[i].length() > 0)
                {
                    param.SetName(StringUtils::ConvertMBSToWCS(m_Info.m_ClientInfoNames[i]));
                    param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ClientInfoValues[i]));
                    appMetaInfoList.push_back(param);
                    //                    appMetaInfoList.push_back(param);
                }
            }
            data->AddAppMetaInfo(appMetaInfoList);
            
            ParameterNameValueTypeList summaryParamList;
            std::string blankStr(256, ' ');
            for (int source = 0; source < m_ChipSummaries.size(); source++) {
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
                        Err::errAbort("QuantMethodMultiDataCCCHPReport: Unable to handle unknown type: " +
                                      ToStr(metricDefs[i].m_type) );
                    }
                    summaryParamList.push_back(param);
                }
            }
            data->AddSummaryParams(summaryParamList);
            
            CHPMultiDataFileWriter writer(*data);
            delete data;
        }
        catch (...)
        {
            Err::errAbort("QuantMethodMultiDataCHPReport::prepare() - Unable to write header and/or precreate signal entries to file: " + m_CHPFileNames[i] + ".tmp");
        }
    }

    // Initialize genotype entry buffer writer.
    std::vector<MultiDataType> dataTypes;
    std::map<MultiDataType, int> maxLen;
    if (nGenotype > 0)
    {
        dataTypes.push_back(GenotypeMultiDataType);
        maxLen[GenotypeMultiDataType] = maxProbeSetNameLength;
    }
    if (nExpression > 0)
    {
        dataTypes.push_back(ExpressionMultiDataType);
        maxLen[ExpressionMultiDataType] = maxProbeSetNameLength;
    }
    m_GenotypeEntryBufferWriter.Initialize(&m_TmpChpFiles, dataTypes, maxLen);
    int buffer_flush_limit = 2*MAX_BUFFER_SIZE;
    // 2**17 == 131072 -- a power of 2 close to 100k
    if (131072 * nfiles > buffer_flush_limit) {
        buffer_flush_limit = 131072 * nfiles;
    }
    m_GenotypeEntryBufferWriter.SetMaxBufferSize(buffer_flush_limit);
    return true;
}

/**
 * After every probeset computation this function is called an is an opportunity
 * to query the quantification method for results, residuals, etc.
 *
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 *
 * @return true if success, false otherwise.
 */
bool QuantMethodMultiDataCCCHPReport::report(ProbeSetGroup &psGroup, 
                                             QuantMethod &qMethod,
                                             const IntensityMart &iMart,
                                             std::vector<ChipStream *> &iTrans,
                                             PmAdjuster &pmAdjust)
{
    QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(&qMethod);
    QuantExprMethod *eMethod = dynamic_cast<QuantExprMethod *>(&qMethod);

    assert(psGroup.probeSets.size() > 0 && psGroup.probeSets[0]);

    if (m_ProbeSetsToReport != NULL) {
        if (m_ProbeSetsToReport->find(psGroup.name)==m_ProbeSetsToReport->end()) {
            // @todo: should this return true or false?
            m_ProcessedProbeSetCount++;
            return false;
        }
    }

    //checkCurrentId(psGroup); we are rewriting things, so skip overhead
    // Create an expression CHP entry from analysis results.

    const ProbeSet::ProbeSetType psType = psGroup.probeSets[0]->psType;
    if (psType == ProbeSet::Expression && nExpression > 0)
    {
        int nChips = 0;

        if (gMethod != NULL) {
            nChips = gMethod->getNumCalls();
        } 
        else if (eMethod != NULL) {
            nChips = eMethod->getNumTargets();
        } 
        else {
            Err::errAbort("Unknown quantification type " + qMethod.getType());
        }

        ProbeSetMultiDataExpressionData entry;
        for (int target=0; target<nChips; target++)
        {
            entry.name = psGroup.name;

            if (gMethod != NULL) {
                entry.quantification = gMethod->getConfidence(target);
            } 
            else { // eMethod != NULL because of above condition
                entry.quantification = eMethod->getSignalEstimate(target);
            } 

            m_GenotypeEntryBufferWriter.WriteMultiDataExpressionEntry(ExpressionMultiDataType, target, entry);
        }
        m_GoodProbesetsExpr.push_back(m_ProcessedProbeSetCount);
        m_GoodProbesetsExprIndex.push_back(m_ProcessedExprProbeSetCount);

        m_ProcessedExprProbeSetCount++;
    }
    // Create a genotype CHP entry from analysis results.
    else if ((psType == ProbeSet::GenoType || psType == ProbeSet::Marker || psType == ProbeSet::MultichannelMarker) && nGenotype > 0)
    {
        if (gMethod == NULL) { Err::errAbort("Can only generate genotype output in MutiData CHP file when using a QuantGTypeMethod."); }

        ProbeSetMultiDataGenotypeData entry;
        entry.metrics.resize(3);
        double avalue;
        double bvalue;
        for (int target=0; target<gMethod->getNumCalls(); target++)
        {
            ///@todo should check that the call encoding is the expected GType encoding
            int call = gMethod->getCall(target);
            u_int8_t forcedCall = (u_int8_t)GTypeCallToChpValue(GType_from_int(gMethod->getForcedCall(target)));
            double confidence = gMethod->getConfidence(target);
            entry.name = psGroup.name;

            // Translate our representation to the CHP call representation.
            entry.call = (u_int8_t)GTypeCallToChpValue(GType_from_int(call));
            entry.confidence = confidence;
                
            gMethod->getAlleleValues(target, avalue, bvalue);
            entry.metrics[0].SetValueFloat(avalue);
            entry.metrics[1].SetValueFloat(bvalue);
            entry.metrics[2].SetValueUInt8(forcedCall);
            m_GenotypeEntryBufferWriter.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, target, entry);
        }
        m_GoodProbesetsGType.push_back(m_ProcessedProbeSetCount);
        m_GoodProbesetsGTypeIndex.push_back(m_ProcessedGTypeProbeSetCount);

        m_ProcessedGTypeProbeSetCount++;
    }

    m_ProcessedProbeSetCount++;

    return true;
}

/**
 * For genotyping CHP files GTYPE likes to have the median of the raw intensity
 * probes stuffed into the confidence field of the CHP files. I'm not sure why
 * but Richard C asked for it so here it is.
 *
 * @param psGroup - Groupt of probes.
 * @param iMart - Raw intensity data for chips.
 * @param chipIx - Which chip we are taking the median for.
 *
 * @return median of all raw intensities for the PM probes in psGroup.
 */
float QuantMethodMultiDataCCCHPReport::medianOfPmProbes(ProbeSetGroup &psGroup, 
const IntensityMart &iMart, int chipIx) {
    vector<float> pm;
    float med = 0;
    for (unsigned int psIx = 0; psIx < psGroup.probeSets.size(); psIx++) {
        const ProbeSet *ps = psGroup.probeSets[psIx];
        for (unsigned int atomIx = 0; atomIx < ps->atoms.size(); atomIx++) {
            Atom *a = ps->atoms[atomIx];
            for (unsigned int probeIx = 0; probeIx < a->probes.size(); probeIx++) {
                Probe *p = a->probes[probeIx];
                if (Probe::isPm(*p)) {
                    pm.push_back(iMart.getProbeIntensity(p->id, chipIx));
                }
            }
        }
    }
    if (pm.empty()) {
        Err::errAbort("QuantMethodMultiDataCHPReport::medianOfPmProbes() - No Pm Probes.");
    }
    med = median(pm.begin(), pm.end());
    return med;
}

/**
 * If a probeset fails to compute for whatever reason then this method is
 * called rather than the normal report call above. By default does nothing.
 *
 * @param psGroup - List of probesets from which probes were used.
 * @param qMethod - Quantification method with compute method called.
 *
 * @return true if success, false otherwise.
 */
bool QuantMethodMultiDataCCCHPReport::reportFailure(ProbeSetGroup &psGroup, 
                                                    QuantMethod &qMethod,
                                                    const IntensityMart &iMart, 
                                                    std::vector<ChipStream *> &iTrans,
                                                    PmAdjuster &pmAdjust)
{
    /// @todo CLARIFY WHAT SHOULD HAPPEN TO CONTROL PROBESETS (ie AFFX-5Q-123)
    /// currently none of them make it out because:
    ///   (1) these probesets get flagged as Expression
    ///   (2) in the reporter construction in Genotype Engine, expression output is disabled
    ///   (3) when we rewrite the chp files in finish() these are dropped because they are not flagged as "good"
    assert(psGroup.probeSets.size() > 0);
    //checkCurrentId(psGroup); we are rewriting things, so skip overhead

    // Create a blank entry in order to be consistent with psGroups list.
    const ProbeSet::ProbeSetType psType = ChipLayout::getProbeSetType(psGroup.probeSets[0]);
    if (psType == ProbeSet::Expression && nExpression > 0)
    {
        ProbeSetMultiDataExpressionData entry;
        for (int target=0; target<m_CHPFileNames.size(); target++)
        {
            entry.name = psGroup.name;

            entry.quantification = medianOfPmProbes(psGroup, iMart, target);
            m_GenotypeEntryBufferWriter.WriteMultiDataExpressionEntry(ExpressionMultiDataType, target, entry);
        }
    }
    else if (psType == ProbeSet::GenoType && nGenotype > 0)
    {
        ProbeSetMultiDataGenotypeData entry;
        entry.metrics.resize(3);
        for (int target=0; target<m_CHPFileNames.size(); target++)
        {
            entry.name = psGroup.name;

            entry.call = ALLELE_NO_CALL;
            entry.confidence = medianOfPmProbes(psGroup, iMart, target);
            entry.metrics[0].SetValueFloat(0.0f);
            entry.metrics[1].SetValueFloat(0.0f);
            entry.metrics[2].SetValueUInt8(ALLELE_NO_CALL);
            m_GenotypeEntryBufferWriter.WriteMultiDataGenotypeEntry(GenotypeMultiDataType, target, entry);
        }
    }
    // update counter to reflect that this probeset has been accounted for
    // this is also incremented in ::report() method
    m_ProcessedProbeSetCount++;

    return true;
}

/**
 * No more probesets will be processed, this is a chance to finish outputting
 * results and clean up.
 * @param qMethod - Quantification method that was used.
 * @return true if success, false otherwise.
 */
bool QuantMethodMultiDataCCCHPReport::finish(QuantMethod &qMethod)
{
    // Sanity to check we saw all the probe sets we were expecting.
    if (m_ProcessedProbeSetCount != m_Info.m_NumProbeSets)
    {
        Err::errAbort("QuantMethodMultiDataCCCHPReport::finish() - Expecting: " + ToStr(m_Info.m_NumProbeSets) +
                      " but got: " + ToStr(m_ProcessedProbeSetCount) + 
                      ". Command Console CHP file will be corrupt.");
    }

    // Flush remaining signal entries in the buffer.
    m_GenotypeEntryBufferWriter.FlushBuffer();

    // Rewrite CHP files to get chip summary entires
    Verbose::out(1,"Creating final files for CHP output");
    Verbose::progressBegin(1, ToStr("Finalizing Multi Data CHP Files"), 
                           m_CHPFileNames.size(), 1, m_CHPFileNames.size());
    try {
        for (unsigned int chip = 0; chip < m_CHPFileNames.size(); chip++) 
        {
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
                        Err::errAbort("QuantMethodMultiDataCCCHPReport: metric '" + metricDefs[i].m_name + "' was not found");
                    }
                    std::wstring mName(CHIP_SUMMARY_PARAMETER_NAME_PREFIX);
                    mName += StringUtils::ConvertMBSToWCS(metric.m_Name);
                    ParameterNameValueType param;
                    if (hdr->FindNameValParam(mName, param) == false) {
                        Err::errAbort("QuantMethodMultiDataCCCHPReport: metric name '" + StringUtils::ConvertWCSToMBS(mName) + 
                                      "' could not be found in the header of " + filename);

                    }

                    switch (param.GetParameterType())
                    {
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
                            Err::errAbort(std::string("QuantMethodMultiDataCCCHPReport: string header parameter too long, name = '") +
                                          metric.m_Name + 
                                          "', value = '" + 
                                          metric.m_String + "'");
                        }
                        param.SetValueAscii(metric.m_String, (int) metric.m_String.length());
                        break;

                    default:
                        Err::errAbort("QuantMethodMultiDataCCCHPReport: unknown header parameter type found in file " + filename);
                    }
                    updateHdr.AddNameValParam(param);
                }
            }
            std::ofstream os;
            os.open(filename.c_str(), std::ios::out|std::ios::binary|std::ios::in);
            if (!os) {
                Err::errAbort("QuantMethodMultiDataCCCHPReport: file " + 
                              filename + " could not be opened for writing");
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
    for (unsigned int i = 0; i < m_CHPFileNames.size(); i++) 
    {
        std::string from = m_CHPFileNames[i] + ".tmp";
        std::string to = m_CHPFileNames[i];
        if (!Fs::fileRename(from.c_str(),to.c_str())) {
            removeAllChps();
            Err::errAbort("Unable to rename '" + from + "' to '" + to +"'");
        }
    }
    removeTmpChps();

    return true;
}

void QuantMethodMultiDataCCCHPReport::addStdHeaders(QuantMethodReport *qReport,
                   const std::string& execGuid, 
                   const std::string& reportGuid,
                   const std::string& timeStr,
                   const std::string& commandLine,
                   const std::string& execVersion,
                   const AnalysisInfo& info) {
  nExpression = info.m_NumExpression;
  nGenotype = info.m_NumGenotyping;
  nReporting = info.m_NumReporting;
  m_Info = info;
  m_AlgName = info.m_AlgName;
  /* Most of the headers get written in the prepare statement... */
}
