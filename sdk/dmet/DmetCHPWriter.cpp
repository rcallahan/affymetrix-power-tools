////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
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

//
#include "dmet/DmetCHPWriter.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPMultiDataFileWriter.h"
#include "chipstream/ProbeListFactory.h"
#include "chipstream/TsvReport.h"
#include "file/TsvFile/PgfFile.h"
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "stats/stats.h"
#include "util/AffxConv.h" 
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Fs.h"
#include "util/RowFile.h"
#include "util/Verbose.h"
//
#include <cerrno>
#include <cstring>
#include <set>
#include <string>
#include <vector>

using namespace std;
using namespace affx;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_utilities;

DmetCHPWriter::Reg DmetCHPWriter::reg;

DmetCHPWriter * DmetCHPWriter::FromBase(BaseEngine *engine)
{
	if (engine != NULL && engine->getEngineName() == DmetCHPWriter::EngineName())
		return (DmetCHPWriter *)engine;
	return NULL;
}

// Constructor

DmetCHPWriter::DmetCHPWriter() {
  defineOptions();
  defineStates();
  m_CallCoder = NULL;
}

void DmetCHPWriter::setHeader(CHPMultiDataData  *data, ParameterNameValueTypeList &params, ParameterNameValueTypeList &metrics){

  ParameterNameValueType param;

  // Prepare headers for all CHP files.
  data->SetAlgName(StringUtils::ConvertMBSToWCS(m_Info.m_AlgName));
  data->SetAlgVersion(StringUtils::ConvertMBSToWCS(m_Info.m_AlgVersion));
  data->SetArrayType(StringUtils::ConvertMBSToWCS(m_Info.m_ChipType));

  param.SetName(StringUtils::ConvertMBSToWCS("program-name"));
  param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ProgramName));
  data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
  param.SetName(StringUtils::ConvertMBSToWCS("program-version"));
  param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ProgramVersion));
  data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
  param.SetName(StringUtils::ConvertMBSToWCS("program-company"));
  param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ProgramCompany));
  data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);

  // Encoding parameters
  ///@todo we should check these against what was used for the APG output
  param.SetName(StringUtils::ConvertMBSToWCS("call-encoding-version"));
  param.SetValueAscii(getOpt("call-coder-version"));
  data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
  param.SetName(StringUtils::ConvertMBSToWCS("call-encoding-max-alleles"));
  param.SetValueInt32(getOptInt("call-coder-max-alleles"));
  data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
  param.SetName(StringUtils::ConvertMBSToWCS("call-encoding-data-size"));
  param.SetValueAscii(getOpt("call-coder-type"));
  data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);

  // Transformation
  ///@todo do not hard code
  param.SetName(StringUtils::ConvertMBSToWCS("apt-state-transform"));
  param.SetValueText(StringUtils::ConvertMBSToWCS("MvA"));
  params.push_back(param);

  // Other parameters from m_Info
  for(int i=0; i<m_Info.m_ParamNames.size(); i++) 
  {
      if (m_Info.m_ParamValues[i].length() > 0)
      {
          param.SetName(StringUtils::ConvertMBSToWCS(m_Info.m_ParamNames[i]));
          param.SetValueText(StringUtils::ConvertMBSToWCS(m_Info.m_ParamValues[i]));
          params.push_back(param);
      }
  }

  data->AddAlgParams(params);
  data->AddSummaryParams(metrics);

  ///@todo propigate meta info from a5 input files for each sub-engine
}  

void DmetCHPWriter::setCopyNumberData(
        CHPMultiDataFileWriter &writer, 
        vector<string> copyNumberRegionNames,
        string celFileName ){

  map<string,DmetCopyNumberData> cnData;

  string copyNumberFileName = getOpt("a5-copynumber");
  bool haveCnData = true;
  if(copyNumberFileName == "") {
      haveCnData = false;
  }

  if(haveCnData) {

    affx::File5_File copyNumberFile5;
    affx::File5_Tsv* copyNumberTsv5 = NULL;

    try {

        copyNumberFile5.open(copyNumberFileName, affx::FILE5_OPEN_RO);
        string copyNumberPrefix = getOpt("set-analysis-name") + ".copynumber";
        copyNumberTsv5= copyNumberFile5.openTsv(copyNumberPrefix);
        if(copyNumberTsv5==NULL)
            APT_ERR_ABORT("Unable to get copy number data from a5 file while writing CHP file.");

        string justReadCelFileName;

        while(copyNumberTsv5->nextLevel(0)==affx::FILE5_OK){
            DmetCopyNumberData data;
            int dataInt;
            double dataDouble;
    
            if(copyNumberTsv5->get(0, "cel", &justReadCelFileName)!=affx::FILE5_OK)
                APT_ERR_ABORT("Unable to get column 'cel' from CN file '"+copyNumberFileName+"'.");

            if(copyNumberTsv5->get(0, "region", &data.name)!=affx::FILE5_OK)
                APT_ERR_ABORT("Unable to get column 'cel' from CN file '"+copyNumberFileName+"'.");

            dataInt = m_CallCoder->abstractAlleleToGenotypeCallNum("NotAvailable");   
            if(copyNumberTsv5->get(0, "cn_call", &dataInt) != FILE5_OK)
                APT_ERR_ABORT("Unable to get column 'cn_call' from CN file '"+copyNumberFileName+"'.");
            data.call = dataInt;

            dataDouble = -2.0;
            if(copyNumberTsv5->get(0, "pvalue", &dataDouble) != FILE5_OK)
                APT_ERR_ABORT("Unable to get column 'pvalue' from CN file '"+copyNumberFileName+"'.");
            data.confidence = dataDouble;

            dataInt = m_CallCoder->abstractAlleleToGenotypeCallNum("NotAvailable");   
            if(copyNumberTsv5->get(0, "cn_force", &dataInt) != FILE5_OK)
                APT_ERR_ABORT("Unable to get column 'cn_force' from CN file '"+copyNumberFileName+"'.");
            data.force = dataInt;

            dataDouble = -2.0;
            if(copyNumberTsv5->get(0, "cn_estimate", &dataDouble) != FILE5_OK)
                APT_ERR_ABORT("Unable to get column 'cn_estimate' from CN file '"+copyNumberFileName+"'.");
            data.estimate=dataDouble;

            dataInt = -2;
            if(copyNumberTsv5->get(0, "cn_lower", &dataInt) != FILE5_OK)
                APT_ERR_ABORT("Unable to get column 'cn_lower' from CN file '"+copyNumberFileName+"'.");
            data.lower = dataInt;

            dataInt = -2;
            if(copyNumberTsv5->get(0, "cn_upper", &dataInt) != FILE5_OK)
                APT_ERR_ABORT("Unable to get column 'cn_upper' from CN file '"+copyNumberFileName+"'.");
            data.upper = dataInt;

            if(justReadCelFileName == celFileName){
                if(cnData.find(data.name) != cnData.end())
                    APT_ERR_ABORT("Duplicate entries found for CN region '"+data.name+"' in sample '"+celFileName+"'");
                cnData[data.name] = data;
            }
        }

        copyNumberTsv5->close();
        Freez(copyNumberTsv5);
        copyNumberFile5.close(); 
    }
    catch(...)
    {
        Freez(copyNumberTsv5);
        throw;
    }
  }

  string regionName;
  vector<string>::const_iterator beginIter=copyNumberRegionNames.begin();
  vector<string>::const_iterator endIter=copyNumberRegionNames.end();
  writer.SeekToDataSet(DmetCopyNumberMultiDataType);
  for(; beginIter!=endIter; beginIter++){
    regionName = *beginIter;

    if(cnData.find(regionName) != cnData.end()) { // We have valid data and write to CHP file.
      writer.WriteEntry(cnData[regionName]);
    } else {  // We write dummy data  as the probeset has been hidden.
      DmetCopyNumberData entry;
      entry.name =regionName ;
      entry.call = DMET_INVALID_INTEGER;
      entry.confidence = DMET_INVALID_DOUBLE;
      entry.force = DMET_INVALID_INTEGER;
      entry.estimate = DMET_INVALID_DOUBLE;
      entry.lower = DMET_INVALID_INTEGER;
      entry.upper = DMET_INVALID_INTEGER;

      writer.WriteEntry(entry);
    }// end of else of copyNumberRegionIsInCNFile check     

  }// end of the loop over the complete list of probesets.
}

void parseContext(const string &context, string &aAllele, string &bAllele, string &aContext, string &bContext) {
  bContext = context;
  bAllele = bContext.substr(0,bContext.rfind(","));
  aContext = bAllele.substr(0,bAllele.rfind(","));
  aAllele = aContext.substr(0,aContext.rfind(","));
  bContext = bContext.substr(bContext.rfind(",")+1);
  bAllele = bAllele.substr(bAllele.rfind(",")+1);
  aContext = aContext.substr(aContext.rfind(",")+1);
}

void computeMaxSignal(const string &psName, const int maxAlleles, const int maxContext, map<string, double> &summaries, vector<double> &maxSignals, vector<int> &maxContexts) {
    maxSignals.resize(maxAlleles);
    maxContexts.resize(maxAlleles);
    for(int aId = 0; aId < maxAlleles; aId++) {
        int context = DMET_INVALID_UINT;
        double signal = DMET_INVALID_DOUBLE;
        for(int cId = 0; cId <= maxContext; cId++) {
            string name = psName + "-" + ToStr(aId) + "-" + ToStr(cId);
            if(summaries.find(name) != summaries.end()) {
                double newSignal = summaries[name];
                if(newSignal > signal) {
                    signal = newSignal;
                    context = cId;
                }
            }
        }
        maxContexts[aId] = context;
        maxSignals[aId] = signal;
    }
}

void DmetCHPWriter::setGenotypeData(
        CHPMultiDataFileWriter &writer, 
        vector<string> const &probesetNames,
        vector<int> const &probesetAlleles,
        map<string,int> &alleleMap,
        string celFileName,
        ChipLayout &layout,
        map<string, int> &calls,
        map<string, int> &forcedCalls,
        map<string, double> &confidences,
        map<string, double> &summaries,
        map<string, string> &contexts
        ){

  // Write out bi-alleleic markers on first pass
  writer.SeekToDataSet(DmetBiAllelicMultiDataType);
  for(int pIx = 0; pIx < probesetNames.size(); pIx++) {
      string psName = probesetNames[pIx];
      int allele = probesetAlleles[alleleMap[psName]];
      if(allele==2){
        if(calls.find(psName) == calls.end()) {
            DmetBiAllelicData entry;
            entry.name = psName;
            entry.call = m_CallCoder->abstractAlleleToGenotypeCallNum("NotAvailable");   
            entry.confidence = DMET_INVALID_DOUBLE;
            entry.force = entry.call;
            entry.signalA = DMET_INVALID_DOUBLE;
            entry.signalB = DMET_INVALID_DOUBLE;
            entry.contextA = DMET_INVALID_UINT;
            entry.contextB = DMET_INVALID_UINT;
            writer.WriteEntry(entry);
        } else {
            DmetBiAllelicData entry;
            entry.name = psName;
            entry.call = calls[psName];
            entry.force = forcedCalls[psName];
            entry.confidence = confidences[psName];
    
            // Parse out allele and context
            string context = contexts[psName];
            string aAllele, bAllele, aContext, bContext;
            parseContext(context, aAllele, bAllele, aContext, bContext);

            // Set allele signal and context
            entry.signalA = summaries[psName + "-" + aAllele + "-" + aContext];
            entry.signalB = summaries[psName + "-" + bAllele + "-" + bContext];
            if(getOptBool("null-context")) {
                entry.contextA = DMET_INVALID_UINT;
                entry.contextB = DMET_INVALID_UINT;
            } else {
                entry.contextA = Convert::toInt(aContext);
                entry.contextB = Convert::toInt(bContext);
            }
            writer.WriteEntry(entry);
        }
      }
  }
  // Write out multi-alleleic markers on second pass
  writer.SeekToDataSet(DmetMultiAllelicMultiDataType);
  for(int pIx = 0; pIx < probesetNames.size(); pIx++) {
      string psName = probesetNames[pIx];
      int allele = probesetAlleles[alleleMap[psName]];
      if(allele>2){
        if(calls.find(psName) == calls.end()) {
            DmetMultiAllelicData entry;
            entry.name = psName;
            entry.call = m_CallCoder->abstractAlleleToGenotypeCallNum("NotAvailable");   
            entry.confidence = DMET_INVALID_DOUBLE;
            entry.force = entry.call;
            entry.alleleCount = allele;
            ///@todo cannot assume no gaps in allele designation
            entry.signalA = (allele >= 1 ? DMET_INVALID_DOUBLE : 0.0);
            entry.signalB = (allele >= 2 ? DMET_INVALID_DOUBLE : 0.0);
            entry.signalC = (allele >= 3 ? DMET_INVALID_DOUBLE : 0.0);
            entry.signalD = (allele >= 4 ? DMET_INVALID_DOUBLE : 0.0);
            entry.signalE = (allele >= 5 ? DMET_INVALID_DOUBLE : 0.0);
            entry.signalF = (allele >= 6 ? DMET_INVALID_DOUBLE : 0.0);
            entry.contextA = DMET_INVALID_UINT;
            entry.contextB = DMET_INVALID_UINT;
            entry.contextC = DMET_INVALID_UINT;
            entry.contextD = DMET_INVALID_UINT;
            entry.contextE = DMET_INVALID_UINT;
            entry.contextF = DMET_INVALID_UINT;
            writer.WriteEntry(entry);
        } else {
            DmetMultiAllelicData entry;
            entry.name = psName;
            entry.call = calls[psName];
            entry.force = forcedCalls[psName];
            entry.confidence = confidences[psName];
            entry.alleleCount = allele;
    
            // Parse out allele and context
            string context = contexts[psName];
            string aAllele, bAllele, aContext, bContext;
            parseContext(context, aAllele, bAllele, aContext, bContext);
    
            // Set allele signal and context
            entry.signalA = summaries[psName + "-" + aAllele + "-" + aContext];
            entry.signalB = summaries[psName + "-" + bAllele + "-" + bContext];
            if(getOptBool("null-context")) {
                entry.contextA = DMET_INVALID_UINT;
                entry.contextB = DMET_INVALID_UINT;
            } else {
                entry.contextA = Convert::toInt(aContext);
                entry.contextB = Convert::toInt(bContext);
            }

            vector<double> maxSignals;
            vector<int>    maxContexts;

            // There could be gaps in allele and context IDs
            ProbeListPacked pList = layout.getProbeListByName(psName);
            if (pList.isNull()) {
              APT_ERR_ABORT("Unable to find probeset '"+psName+"' in layout");
            }

            int maxContext = 0;
            for(int bIx=0; bIx < pList.block_cnt(); bIx++)
                if(pList.get_blockContext(bIx) > maxContext)
                    maxContext = pList.get_blockContext(bIx);
            computeMaxSignal(psName, 6, maxContext, summaries, maxSignals, maxContexts);
    
            // set max signal for each allele over all context
            ///@todo cannot assume no gaps in allele designation
            entry.signalA = (allele >= 1 ? maxSignals[0] : 0.0);
            entry.signalB = (allele >= 2 ? maxSignals[1] : 0.0);
            entry.signalC = (allele >= 3 ? maxSignals[2] : 0.0);
            entry.signalD = (allele >= 4 ? maxSignals[3] : 0.0);
            entry.signalE = (allele >= 5 ? maxSignals[4] : 0.0);
            entry.signalF = (allele >= 6 ? maxSignals[5] : 0.0);
            // set context with max signal for each allele
            ///@todo cannot assume no gaps in allele designation
            if(getOptBool("null-context")) {
                entry.contextA = DMET_INVALID_UINT;
                entry.contextB = DMET_INVALID_UINT;
                entry.contextC = DMET_INVALID_UINT;
                entry.contextD = DMET_INVALID_UINT;
                entry.contextE = DMET_INVALID_UINT;
                entry.contextF = DMET_INVALID_UINT;
            } else {
                entry.contextA = (allele >= 1 ? maxContexts[0] : DMET_INVALID_UINT);
                entry.contextB = (allele >= 2 ? maxContexts[1] : DMET_INVALID_UINT);
                entry.contextC = (allele >= 3 ? maxContexts[2] : DMET_INVALID_UINT);
                entry.contextD = (allele >= 4 ? maxContexts[3] : DMET_INVALID_UINT);
                entry.contextE = (allele >= 5 ? maxContexts[4] : DMET_INVALID_UINT);
                entry.contextF = (allele >= 6 ? maxContexts[5] : DMET_INVALID_UINT);
            }
            writer.WriteEntry(entry);
        }
      } 
  }// end of the loop over the complete list of probesets.
}

void DmetCHPWriter::fillInReportMetrics(const string &reportFile, const string &celFile, ParameterNameValueTypeList &metrics) {
    // a bit of a waste in that we keep re-reading this file for each cel file
    affx::TsvFile tsv;

    /* Some values one might see:
        cel_files
        computed_gender
        call_rate
        het_rate
        hom_rate
        cluster_distance_mean
        cluster_distance_stdev
        raw_intensity_mean
        raw_intensity_stdev
        allele_summarization_mean
        allele_summarization_stdev
        allele_deviation_mean
        allele_deviation_stdev
        allele_mad_residuals_mean
        allele_mad_residuals_stdev
        cn-probe-chrXY-ratio_gender_meanX
        cn-probe-chrXY-ratio_gender_meanY
        cn-probe-chrXY-ratio_gender_ratio
        cn-probe-chrXY-ratio_gender
    */

    //WARNING: We do not automatically load everything. One reason would be that if
    //         some one consented everything but one marker, one could probably use
    //         certain fields to infer the genotype of the missing marker

    string foundCelFile;
    double callRate = 0.0;
    string computedGender = "unknown";
    double cnMeanX = 0.0, cnMeanY = 0.0, cnRatio = 0.0;
    string cnGender = "unknown";

    tsv.bind(0,"cel_files",      &foundCelFile,         affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"total_call_rate",      &callRate,        affx::TSV_BIND_REQUIRED);
    tsv.bind(0,"computed_gender",&computedGender,  affx::TSV_BIND_REQUIRED);

    tsv.bind(0,"cn-probe-chrXY-ratio_gender_meanX",&cnMeanX,   affx::TSV_BIND_OPTIONAL);
    tsv.bind(0,"cn-probe-chrXY-ratio_gender_meanY",&cnMeanY,   affx::TSV_BIND_OPTIONAL);
    tsv.bind(0,"cn-probe-chrXY-ratio_gender_ratio",&cnRatio,   affx::TSV_BIND_OPTIONAL);
    tsv.bind(0,"cn-probe-chrXY-ratio_gender",      &cnGender,  affx::TSV_BIND_OPTIONAL);

    if(tsv.open(reportFile) != affx::TSV_OK)
        APT_ERR_ABORT("Couldn't open file: '" + reportFile + "' to read.");
    bool found = false;
    while((tsv.nextLevel(0) == affx::TSV_OK) && !found) {
        if(foundCelFile == celFile) {

            ParameterNameValueType param;

            param.SetName(StringUtils::ConvertMBSToWCS("qc_call_rate"));
            callRate = roundDouble(callRate, 2);
            param.SetValueFloat(callRate);
            metrics.push_back(param);

            param.SetName(StringUtils::ConvertMBSToWCS("computed_gender"));
            param.SetValueAscii(computedGender);
            metrics.push_back(param);

            if(tsv.cname2cidx(0, "cn-probe-chrXY-ratio_gender_meanX") != TSV_ERR_NOTFOUND) {
                param.SetName(StringUtils::ConvertMBSToWCS("cn-probe-chrXY-ratio_gender_meanX"));
                param.SetValueFloat(cnMeanX);
                metrics.push_back(param);
            }

            if(tsv.cname2cidx(0, "cn-probe-chrXY-ratio_gender_meanY") != TSV_ERR_NOTFOUND) {
                param.SetName(StringUtils::ConvertMBSToWCS("cn-probe-chrXY-ratio_gender_meanY"));
                param.SetValueFloat(cnMeanY);
                metrics.push_back(param);
            }

            if(tsv.cname2cidx(0, "cn-probe-chrXY-ratio_gender_ratio") != TSV_ERR_NOTFOUND) {
                param.SetName(StringUtils::ConvertMBSToWCS("cn-probe-chrXY-ratio_gender_ratio"));
                param.SetValueFloat(cnRatio);
                metrics.push_back(param);
            }

            if(tsv.cname2cidx(0, "cn-probe-chrXY-ratio_gender") != TSV_ERR_NOTFOUND) {
                param.SetName(StringUtils::ConvertMBSToWCS("cn-probe-chrXY-ratio_gender"));
                param.SetValueAscii(cnGender);
                metrics.push_back(param);
            }

            found = true;
        }
    }

    tsv.close();
}

void DmetCHPWriter::parseRegionsFile(vector<string> &copyNumberRegionNames, int &regionCount) {
  // We first read the complete list of CN regions from the region-model.txt file.
  affx::TsvFile tsv;
  string regionFileName = getOpt("region-model");
  if (tsv.open(regionFileName) != TSV_OK)
    APT_ERR_ABORT("Couldn't open " + regionFileName + " to read while writing CHP files.");

  string regionName;
  while (tsv.nextLevel(0)==TSV_OK){
    if ((tsv.get(0,"region",regionName)!=TSV_OK) || (regionName=="")) {
      APT_ERR_ABORT("No region name or empty region name was found while reading the region-model file while writing CHP files.");
    }
    copyNumberRegionNames.push_back(regionName);
    regionCount++;
  }
} 

void DmetCHPWriter::loadLayout(ChipLayout &layout, std::vector<string> &probesetNames, std::vector<int> &alleles) {

  string cdfFile = getOpt("cdf-file");
  string spfFile = getOpt("spf-file");

  std::vector<bool> probeSubset;
  probeidmap_t killList;
  std::set<const char *, Util::ltstr> probeSetsToLoad;
  std::set<affxcdf::GeneChipProbeSetType> psTypesToLoad;
  psTypesToLoad.insert(affxcdf::MarkerProbeSetType);
  vector<const char *> psNames;

  try {
    if(cdfFile != "") {
      Verbose::out(1, ToStr("Opening layout file: ") + cdfFile);
      if(!layout.openCdf(cdfFile, probeSetsToLoad, &psNames,
                         probeSubset, "", killList, false, psTypesToLoad)) {
        APT_ERR_ABORT("Couldn't open layout file: " + cdfFile);
      }
    }
    else if(spfFile != "") {
      Verbose::out(1, ToStr("Opening layout file: ") + spfFile);
      layout.openSpf(spfFile, probeSetsToLoad, &psNames, probeSubset,
                           "", false, psTypesToLoad);
    }
    else {
      APT_ERR_ABORT("Must have either a cdf file or spf file.");
    }
  }
  catch(const Except &e) {
    for(uint32_t i = 0; i < psNames.size(); i++) {
      Freez(psNames[i]);
    }
    psNames.clear();
    APT_ERR_ABORT(e.what());
  }
  catch(...) {
    for(uint32_t i = 0; i < psNames.size(); i++) {
      Freez(psNames[i]);
    }
    psNames.clear();
    APT_ERR_ABORT("Unknown exception caught - terminating.");
  }
  Verbose::out(2, "Loaded " + ToStr(psNames.size()) + " probesets.");

  for(int i=0; i<psNames.size(); i++){
      probesetNames.push_back(string(psNames[i]));
  }

  for(uint32_t i = 0; i < psNames.size(); i++) {
      Freez(psNames[i]);
  }
  psNames.clear();

  // Figure out how many alleles for each probeset
  for(int pIx=0; pIx < layout.getProbeSetCount(); pIx++) {
    ProbeListPacked pList = layout.getProbeListAtIndex(pIx);
    map<int,int> alleleMap;
    for(int bIx=0; bIx < pList.block_cnt(); bIx++) {
        alleleMap[pList.get_blockAllele(bIx)] = 1;
    }
    alleles.push_back(alleleMap.size());
  }
}

void DmetCHPWriter::fillInCallMetrics(map<string,int> &calls, ParameterNameValueTypeList &metrics) {
    map<string,int>::iterator iter = calls.begin();
    int count=0, het=0, hom=0, pra=0, nc=0, zcn=0;
    for(;iter != calls.end(); iter++) {
        count++;
        if(m_CallCoder->abstractAlleleToGenotypeCallNum("PossibleRareAllele") == iter->second){
            pra++;
        } else if(m_CallCoder->abstractAlleleToGenotypeCallNum("ZeroCopyNumber") == iter->second){
            zcn++;
        } else if(m_CallCoder->abstractAlleleToGenotypeCallNum("NoCall") == iter->second){
            nc++;
        } else if(m_CallCoder->isHom(iter->second)){
            hom++;
        } else if(m_CallCoder->isHet(iter->second)){
            het++;
        } else {
            Verbose::out(1,"Unaccounted for call (" + ToStr(iter->second) + "). Treating as no-call for summary stats.");
            nc++;
        }
    }

    double callRate = (count - nc) / double(count) * 100;
    double homRate = (hom) / double(count) * 100;
    double hetRate = (het) / double(count) * 100;
    //double praRate = (pra) / double(count) * 100;

    callRate = roundDouble(callRate, 2);
    homRate = roundDouble(homRate, 2);
    hetRate = roundDouble(hetRate, 2);

    // Add Call Rate
    ParameterNameValueType param;

    param.SetName(StringUtils::ConvertMBSToWCS("call_rate"));
    param.SetValueFloat(callRate);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("het_rate"));
    param.SetValueFloat(hetRate);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("hom_rate"));
    param.SetValueFloat(homRate);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("call_count"));
    param.SetValueInt32(count - nc);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("hom_count"));
    param.SetValueInt32(hom);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("het_count"));
    param.SetValueInt32(het);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("pra_count"));
    param.SetValueInt32(pra);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("cn0_count"));
    param.SetValueInt32(zcn);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("no_call_count"));
    param.SetValueInt32(nc);
    metrics.push_back(param);
    param.SetName(StringUtils::ConvertMBSToWCS("probeset_count"));
    param.SetValueInt32(count);
    metrics.push_back(param);

    Verbose::out(3,"Call Metrics: " +
            ToStr(count) + " markers, " +
            ToStr(het) + " het, " +
            ToStr(hom) + " hom, " +
            ToStr(pra) + " pra, " +
            ToStr(zcn) + " cn0, " +
            ToStr(nc) + " nc");
}

void DmetCHPWriter::runImp(){
 
  vector<string> celFiles = getOptVector("cels");
  setOpt("cels",celFiles);

  // Load up chip layout info
  int biAlleleCount=0;
  int multiAlleleCount=0;
  ChipLayout layout;
  vector<string> probesetNames; /// The names of the probesets (lib file order)
  map<string,int> alleleMap; 
  vector<int> probesetAlleles; /// The number of alleles in the probesets (lib file order)
  CHPMultiDataData *data = NULL;

  try {

    m_CallCoder = new GenoCallCoder(getOptInt("call-coder-max-alleles"),getOpt("call-coder-type"),getOpt("call-coder-version"), '\0');  

    loadLayout(layout, probesetNames, probesetAlleles);
    for(int i=0; i<probesetAlleles.size(); i++) {
        if(probesetAlleles[i] > 2)
            multiAlleleCount++;
        else if(probesetAlleles[i] == 2)
            biAlleleCount++;
    }
    
    // Figure out max probeset name length
    int maxProbeSetNameLength = 0;
    for(int i=0; i<probesetNames.size(); i++) {
        if(probesetNames[i].size() > maxProbeSetNameLength)
            maxProbeSetNameLength = probesetNames[i].size();
        alleleMap[probesetNames[i]] = i;
    }

    sort(probesetNames.begin(),probesetNames.end());
    
    //  These are the copy number regions as defined in the regions file. 
    //  Again, when loaded into the vector they maintain the order as seen in the regions file.
    int regionCount=0;
    vector<string> copyNumberRegionNames;
    parseRegionsFile(copyNumberRegionNames, regionCount);
    
    // Figure out max region name
    int maxRegionNameLength = 0;
    for(int i=0; i<copyNumberRegionNames.size(); i++)
        if(copyNumberRegionNames[i].length() > maxRegionNameLength)
            maxRegionNameLength = copyNumberRegionNames[i].size();

   sort(copyNumberRegionNames.begin(),copyNumberRegionNames.end()); 
    
    // We will first read the summary file to get a list of the CEL files we are dealing with.
    vector<string> celFileNames = getOptVector("cels");
    
    // Setup output folder
    if(!Fs::isWriteableDir(getOpt("out-dir")))
      if(Fs::mkdirPath(getOpt("out-dir"), false) != APT_OK) {
        APT_ERR_ABORT("Can't make or write to directory: " + getOpt("out-dir"));
      }
    
    // Iterate over each sample and generate a CHP file
    vector<string>::iterator begin = celFileNames.begin();
    vector<string>::iterator end = celFileNames.end();
    Verbose::out(1,"Generating CHP files with " + ToStr(biAlleleCount) + " bi-allele, " + ToStr(multiAlleleCount) + " multi-allele, and " + ToStr(regionCount) + " region entries.");
    Verbose::progressBegin(1, ToStr("Generating CHP files"), celFileNames.size(), 1, celFileNames.size());
    for(; begin!=end; begin++){
        Verbose::progressStep(1);
        string celFileName = Fs::basename(*begin);
        string celFilePrefix = Fs::noextname1(celFileName);
        string outFileName = Fs::join(getOpt("out-dir"), celFilePrefix + "." + getOpt("set-analysis-name") + ".chp.tmp");
    
        { // put CHP writer/data in nested scope so that they are destroyed
          // and the file closed before we try and rename.

            data = new CHPMultiDataData(outFileName);

            // set parent header
            FusionCELData cel;
            try {
                cel.SetFileName(begin->c_str());
                if(!cel.ReadHeader()) {
                    APT_ERR_ABORT("Unable to read CEL file " + celFileName);
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
                APT_ERR_ABORT("Unable to read CEL file " + celFileName);
            }
     
            // Initialize the size of the 3 CHP sections and the max length of the probeset
            // name.  These are Calvin format required values. 
            data->SetEntryCount(DmetBiAllelicMultiDataType, biAlleleCount, maxProbeSetNameLength);
            data->SetEntryCount(DmetMultiAllelicMultiDataType, multiAlleleCount, maxProbeSetNameLength);
            data->SetEntryCount(DmetCopyNumberMultiDataType, regionCount, maxRegionNameLength);
    
            // Load up data
            // perhaps we should really keep the File5/Tsv5 objects around rather than open/close for each sample
            // perhaps we should really only load up the probeset names once, not for each sample
            map<string, int> calls;
            fillVectorFromA5<int>(getOpt("a5-calls"), getOpt("set-analysis-name") + ".calls", celFileName, calls);
            map<string, int> forcedCalls;
            fillVectorFromA5<int>(getOpt("a5-forced-calls"), getOpt("set-analysis-name") + ".forced-calls", celFileName, forcedCalls);
            map<string, double> confidences;
            fillVectorFromA5<double>(getOpt("a5-confidences"), getOpt("set-analysis-name") + ".confidences", celFileName, confidences);
            map<string, double> summaries;
            fillVectorFromA5<double>(getOpt("a5-summaries"), getOpt("set-analysis-name") + ".summary", celFileName, summaries);
            map<string, string> contexts;
            fillVectorFromA5<string>(getOpt("a5-context"), getOpt("set-analysis-name") + ".context", celFileName, contexts);
    
            ParameterNameValueTypeList params;
            ParameterNameValueTypeList metrics;

            // Filter Calls
            map<string, int>::iterator iter = calls.begin();
            float cutoff = getOptDouble("geno-call-thresh");
            for(;iter != calls.end(); iter++) {
                if(m_CallCoder->abstractAlleleToGenotypeCallNum("PossibleRareAllele") != iter->second && 
                       m_CallCoder->abstractAlleleToGenotypeCallNum("ZeroCopyNumber")  != iter->second) {
                    if(confidences[iter->first] > cutoff) {
                        ///@todo should we handle cn aware no calls? eg N NN NNN ...
                        iter->second = m_CallCoder->abstractAlleleToGenotypeCallNum("NoCall");
                    } else {
                        iter->second = forcedCalls[iter->first];
                    }
                }
            }
    
            // Write Header
            fillInCallMetrics(calls,metrics);
            fillInReportMetrics(getOpt("report-file"),celFileName,metrics);
    
            setHeader(data, params, metrics);
     
            CHPMultiDataFileWriter writer(*data);
    
            // Populate Data
            setGenotypeData(writer, probesetNames, probesetAlleles, alleleMap, celFileName, layout, calls, forcedCalls, confidences, summaries, contexts); 
            setCopyNumberData(writer, copyNumberRegionNames, celFileName);
  
            // close chp file
            Freez(data);

        }
        std::string from = outFileName;
        std::string to = outFileName.substr(0,outFileName.rfind(".tmp"));
        if (!Fs::fileRename(from,to,false)) {
            APT_ERR_ABORT("Unable to rename '" + from + "' to '" + to +"' ("+strerror(errno)+")");
        }
    } // end for loop over cel files
    Freez(m_CallCoder);
  } // end try
  catch(...)
  {
        Freez(data);
        Freez(m_CallCoder);
        throw;
  }
  Verbose::progressEnd(1, "Done.");
}

void DmetCHPWriter::defineOptions() {

  defineOptionSection("Input Options");
    defineOption("", "a5-copynumber", PgOpt::STRING_OPT,
                     "CopyNumber Engine output in a5 format. ",
                     "");
    defineOption("", "a5-summaries", PgOpt::STRING_OPT,
                     "Probeset summary values in a5 format. ",
                     "");
    defineOption("", "a5-calls", PgOpt::STRING_OPT,
                     "Probeset call values in a5 format. ",
                     "");
    defineOption("", "a5-forced-calls", PgOpt::STRING_OPT,
                     "Probeset call values in a5 format. ",
                     "");
    defineOption("", "a5-confidences", PgOpt::STRING_OPT,
                     "Probeset confidence values in a5 format. ",
                     "");
    defineOption("", "a5-context", PgOpt::STRING_OPT,
                     "Probeset context values in a5 format. ",
                     "");
    defineOption("", "region-model", PgOpt::STRING_OPT,
                     "Regions model parameters. ",
                     "");
    defineOption("c", "cdf-file", PgOpt::STRING_OPT,
                     "File defining probe sets. Use either --cdf-file or --spf-file. ",
                     "");
    defineOption("", "spf-file", PgOpt::STRING_OPT,
                     "File defining probe sets in spf (simple probe format) which is like a text cdf file.",
                     "");
    defineOption("", "report-file", PgOpt::STRING_OPT,
                     "Report file from ProbesetGenotype engine.",
                     "");


  defineOptionSection("Analysis Options");
    defineOption("", "set-analysis-name", PgOpt::STRING_OPT,
                     "Explicitly set the analysis name. "
                     "This affects output file names (ie prefix) and various meta info.",
                     "dmet");
    defineOption("", "geno-call-thresh", PgOpt::DOUBLE_OPT,
                  "The confidence threshold for reporting calls in the CHP file.",
                  "0.1");
    defineOption("", "null-context", PgOpt::BOOL_OPT,
                    "Indicates whether or not context info should be populated in the CHP files.",
                    "true");

  defineOptionSection("Advanced Options");
    defineOption("", "call-coder-max-alleles", PgOpt::INT_OPT,
                    "For encoding/decoding calls, the max number of alleles per marker to allow.",
                    "6");
    defineOption("", "call-coder-type", PgOpt::STRING_OPT,
                    "The data size used to encode the call.",
                    "UCHAR");
    defineOption("", "call-coder-version", PgOpt::STRING_OPT,
                    "The version of the encoder/decoder to use",
                    "1.0");

  defineOptionSection("Engine Options (Not used on command line)");
    defOptMult("", "cels", PgOpt::STRING_OPT,
                     "Cel files to process.",
                     "");
    defineOption("", "batch-name", PgOpt::STRING_OPT,
					 "The name of the batch for the dynamic cluster analysis. ",
					 "");
}

void DmetCHPWriter::defineStates() { }

void DmetCHPWriter::checkOptionsImp() {

    defineStates();

    setLibFileOpt("region-model");
    setLibFileOpt("cdf-file");
    setLibFileOpt("spf-file");
}


