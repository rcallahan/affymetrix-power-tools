////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
 * @file CNLog2RatioData.cpp
 *
 * @brief This file contains the CNLog2RatioData class members.
 */

//
/*
#include "chipstream/TsvReport.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
*/
#include "copynumber/CNLog2RatioData.h"
//
#include "copynumber/Annotation.h"
#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNLog2RatioEngine.h"
#include "copynumber/CNAnalysisMethodCovariateParams.h"
//
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "file5/File5_File.h"
#include "file5/File5_Tsv.h"
#include "portability/affy-base-types.h"
#include "stats/stats.h"
#include "util/AffxArray.h"
#include "util/AffxBinaryFile.h"
#include "util/AffxStatistics.h"
#include "util/AffxString.h"
#include "util/Fs.h"
#include "util/Guid.h"
#include "util/TmpFileFactory.h"
//
#include <cmath>
#include <limits>
//

using namespace std;
using namespace affx;
/**
 * @brief Constructor
 */
CNLog2RatioData::CNLog2RatioData()
{
    m_bLog2RatioEngine = false;
    m_pEngine = NULL;
    m_vProbeSets.clear();

        m_dMuPrimeAA=0;
        m_dMuPrimeAB=0;
        m_dMuPrimeBB=0;

        m_dSigmaPrimeAA=0;
        m_dSigmaPrimeAB=0;
        m_dSigmaPrimeBB=0;

}

/**
 * @brief Destructor
 */
CNLog2RatioData::~CNLog2RatioData()
{
    clear();
}

/**
 * @brief Release the memory associated with the object
 */
void CNLog2RatioData::clear()
{
    m_arExperiments.deleteAll();
    clearProbeSets();
    m_mxAAlleleEstimates.clear();
    m_mxBAlleleEstimates.clear();
    m_mxGenotypeCalls.clear();
    m_mxGenotypeConfidences.clear();

    m_v1.clear();
}

void CNLog2RatioData::clearProbeSets()
{
    m_arProbeSetsAlternateSort.nullAll();
    m_arProbeSets.nullAll();
    m_vProbeSets.clear();
}

/**
 * @brief Load the probe set names from the restrict list file
 * @param const AffxString& - The signal summary file name
 * @param AffxArray<AffxString>& - The vector to load the string values into
 */
void CNLog2RatioData::loadProbeSetNamesFromRestrictList(AffxArray<AffxString>& ar)
{
    Verbose::out(1, "CNLog2RatioData::loadProbeSetNamesFromRestrictList");
    AffxString strProbeListFileName = m_pEngine->getOpt("probeset-ids");
    if (strProbeListFileName != "") {
      TsvFile tsv;
      tsv.m_optAutoTrim = true;
      std::string strProbeSetName;
      tsv.open(strProbeListFileName);
      tsv.bind(0,0, &strProbeSetName,  affx::TSV_BIND_REQUIRED);
      while( tsv.nextLevel(0) == affx::TSV_OK) {
        ar.add(new AffxString(strProbeSetName));
      }
      tsv.clear();
      ar.quickSort(0);
    }
}

/**
 * @brief Load the model file and set up the Probe Set array with the pre-calculated median values from the log2 ratio model.
 * @param const AffxString& - The name of the file to load.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadModelFile(const AffxString& strFileName)
{
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Loading Copy Number Reference File.");
    Verbose::out(3, "CNLog2RatioData::loadModelFile");
    affymetrix_calvin_parameter::ParameterNameValueType param;
    bool bSuccessful = false;
    CNProbeSet objSearch;
    int iSearchIndex = -1;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        try
        {
            AffxString str;
            double d = 0;
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_OPEN);
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                AffxString str;
                tsv5->get(0, 0, &str);
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if (str.startsWith("#%affymetrix-algorithm-param-apt-opt-cel-file-"))
                {
                    int iFindIndex = str.indexOf("=");
                    AffxString strPrompt = str.substring(2, iFindIndex);
                    AffxString strValue = str.substring(iFindIndex + 1);
                    param.SetName(StringUtils::ConvertMBSToWCS(strPrompt));
                    param.SetValueAscii(Fs::basename(strValue));
                    CNAnalysisMethod::getCelFileParams()->push_back(param);
                }
                else if (str.startsWith(CovariateParams::m_CNReferencePrefix))
                {
                    CovariateParams::restoreCovariateParams(str);
                    if (str.find("signal-command") != std::string::npos)
                    {
                        m_pEngine->setOpt("signal-adjustment-covariates", CovariateParams::m_signalCommand);
                    }
                    else if (str.find("log2ratio-command") != std::string::npos)
                    {
                        m_pEngine->setOpt("lr-adjustment-covariates", CovariateParams::m_LRCommand);
                    }
                    else if (str.find("allele-peaks-command") != std::string::npos)
                    {
                        m_pEngine->setOpt("allele-peaks-adjustment-covariates", CovariateParams::m_APCommand);
                    }
                }
            }
            tsv5->close();
            delete tsv5;

            bool useCovariateTsv = m_pEngine->getEngineName() == "CNCytoEngine"  &&
                                   group5->name_exists("Covariates");
            tsv5 = group5->openTsv("Reference", affx::FILE5_OPEN);
            affx::File5_Tsv* tsv5_C = NULL;
            int numCovariates;
            if (useCovariateTsv)
            {
                tsv5_C = group5->openTsv("Covariates", affx::FILE5_OPEN);
                numCovariates = tsv5_C->getColumnCount(0);
            }
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                AffxString strProbeSetName;
                tsv5->get(0, 0, &strProbeSetName);
                objSearch.setProbeSetName(strProbeSetName);
                iSearchIndex = m_arProbeSets.binarySearch(objSearch, 0);
                if (iSearchIndex != -1)
                {
                    CNProbeSet* p = m_arProbeSets.getAt(iSearchIndex);
                    tsv5->get(0, 1, &d); p->setMedianSignal((float)d);
                    tsv5->get(0, 2, &d); p->setXXMedianSignal((float)d);
                    tsv5->get(0, 3, &d); p->setYMedianSignal((float)d);
                    tsv5->get(0, 4, &d); p->setAAMedianSignal((float)d);
                    tsv5->get(0, 5, &d); p->setABMedianSignal((float)d);
                    tsv5->get(0, 6, &d); p->setBBMedianSignal((float)d);
                    if (useCovariateTsv)
                    {
                        if (tsv5_C->nextLine() == affx::FILE5_OK)
                        {
                            for (int i = 0; i < numCovariates; i++)
                            {
                                tsv5_C->get(0, i, &d);
                                p->setCovariateValue(i, (float)d);
                            }
                        }
                    }
                }
            }
            tsv5->close();
            delete tsv5;
            if (useCovariateTsv)
            {
                tsv5_C->close();
                delete tsv5_C;
            }

            group5->close();
            delete group5;
            file5.close();
            bSuccessful = true;
        } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    }
    else
    {
        unsigned int uiRowCount = affx::TsvFile::getLineCountInFile(strFileName);
        if (uiRowCount > 1)
        {
            uiRowCount--;
            affx::TsvFile tsv;
            AffxString key;
            AffxString value;
            tsv.m_optAutoTrim = true;
            if (tsv.open(strFileName) == affx::TSV_OK) {
              for ( tsv.headersBegin(); tsv.headersNext(key,value) == affx::TSV_OK; ) {
                ///@todo AW: should this be "state" not "option"? Do we need to handle both for compatability reasons?
                if (key.startsWith("affymetrix-algorithm-param-apt-opt-cel-file-") ){
                  param.SetName(StringUtils::ConvertMBSToWCS(key));
                  param.SetValueAscii(Fs::basename(value));
                  CNAnalysisMethod::getCelFileParams()->push_back(param);
                }
              }
              const std::string arModelHeaders[] = {
                "probeset_id", "mediansignal","xxmediansignal","ymediansignal",
                "aamediansignal",    "abmediansignal", "bbmediansignal",                
              };
              int iModelHeadersCount = (sizeof( arModelHeaders ) / sizeof(arModelHeaders[0] ) );
              bool okColumns = false;
              if ( tsv.getColumnCount(0) == iModelHeadersCount) {
                AffxString cname;
                okColumns = true;
                for ( int cidx = 0; okColumns && (cidx < iModelHeadersCount) &&
                        (okColumns = (tsv.cidx2cname(0,cidx, cname) == affx::TSV_OK)); cidx++ ) {
                  if ( cname.toLowerCase() != arModelHeaders[cidx] ) {
                    okColumns = false;
                  }
                }
              } 
              if ( !okColumns ) {
                throw(Except("The file specified does not appear to be a model file: " + strFileName));
              }
              else {
                bSuccessful = true;
                int iDotCount = 10;
                Verbose::progressBegin(1, "Loading Model file ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
                std::string col;
                while (tsv.nextLevel(0) == affx::TSV_OK) {
                  Verbose::progressStep(1);
                  int cidx = 0;
                  std::string col;
                  AffxString strProbeSetName;
                  tsv.get(0,cidx++, strProbeSetName);
                  objSearch.setProbeSetName(strProbeSetName);
                  iSearchIndex = m_arProbeSets.binarySearch(objSearch, 0);
                  if (iSearchIndex != -1) {
                    CNProbeSet* p = m_arProbeSets.getAt(iSearchIndex);
                    tsv.get(0,cidx++, col);
                    p->setMedianSignal((float)AffxByteArray(col).parseDouble());
                    tsv.get(0,cidx++, col);
                    p->setXXMedianSignal((float)AffxByteArray(col).parseDouble());
                    tsv.get(0,cidx++, col);
                    p->setYMedianSignal((float)AffxByteArray(col).parseDouble());
                    tsv.get(0,cidx++, col);
                    p->setAAMedianSignal((float)AffxByteArray(col).parseDouble());
                    tsv.get(0,cidx++, col);
                    p->setABMedianSignal((float)AffxByteArray(col).parseDouble());
                    tsv.get(0,cidx++, col);
                    p->setBBMedianSignal((float)AffxByteArray(col).parseDouble());
                  }
                }
                Verbose::progressEnd(1, "Done.");
              }
              tsv.clear();
            } else {throw(Except("Cannot open file: " + strFileName));}
        } else {throw(Except("Model file contains no data: " + strFileName));}
    }
    return bSuccessful;
}

/**
 * @brief Load the model file and set up the Probe Set array with the pre-calculated median values from the log2 ratio model.
 * @param const AffxString& - The name of the file to load.
 * @return bool - true if successful
 */
void CNLog2RatioData::loadCyto2ModelFile(const AffxString& strFileName)
{
    Verbose::out(3, "CNLog2RatioData::loadCyto2ModelFile");
    CNProbeSet objSearch;
    int iSearchIndex = -1;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        try
        {
            double d = 0;
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
            affx::File5_Tsv* tsv5 = group5->openTsv("MedianSignals", affx::FILE5_OPEN);
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                AffxString strProbeSetName;
                tsv5->get(0, 0, &strProbeSetName);
                objSearch.setProbeSetName(strProbeSetName);
                iSearchIndex = m_arProbeSets.binarySearch(objSearch, 0);
                if (iSearchIndex != -1)
                {
                    CNProbeSet* p = m_arProbeSets.getAt(iSearchIndex);
                    tsv5->get(0, 1, &d); p->setMedianSignal((float)d);
                    tsv5->get(0, 2, &d); p->setXXMedianSignal((float)d);
                    tsv5->get(0, 3, &d); p->setYMedianSignal((float)d);
                    tsv5->get(0, 4, &d); p->setAAMedianSignal((float)d);
                    tsv5->get(0, 5, &d); p->setABMedianSignal((float)d);
                    tsv5->get(0, 6, &d); p->setBBMedianSignal((float)d);
                }
            }
            tsv5->close();
            delete tsv5;
            group5->close();
            delete group5;
            file5.close();
        } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    } else {throw(Except("Reference file is not in HDF5 format: " + strFileName));}
}

/**
 * @brief Load the SNP reference file and set up the Probe Set array with the pre-calculated median values from the log2 ratio model.
 * @param const AffxString& - The name of the file to load.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadSnpReferenceFile(const AffxString& strFileName)
{

    Verbose::out(3, "CNLog2RatioData::loadSnpReferenceFile");
    bool bSuccessful = false;
    CNProbeSet objSearch;
    int iSearchIndex = -1;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        try
        {
            AffxString str;
            double d = 0;
            int iInputInteger=0;
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("SNPReference", affx::FILE5_OPEN);
            affx::File5_Tsv* tsv5 = group5->openTsv("SNPReference", affx::FILE5_OPEN);
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                AffxString strProbeSetName;
                tsv5->get(0, 0, &strProbeSetName);
                objSearch.setProbeSetName(strProbeSetName);
                iSearchIndex = m_arProbeSets.binarySearch(objSearch, 0);
                if (iSearchIndex != -1)
                {
                    CNProbeSet* p = m_arProbeSets.getAt(iSearchIndex);
                    tsv5->get(0, 1, &d); p->setMuAB((float)d);
                    tsv5->get(0, 2, &d); p->setMuAA((float)d);
                    tsv5->get(0, 3, &d); p->setMuBB((float)d);
                    tsv5->get(0, 4, &d); p->setInformation((float)d);
                    tsv5->get(0, 5, &iInputInteger); p->setUseInEMAlgorithm(iInputInteger);
                }
            }
            tsv5->close();
            delete tsv5;

            tsv5 = group5->openTsv("FLDValues", affx::FILE5_OPEN);
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                AffxString strProbeSetName;
                tsv5->get(0, 0, &strProbeSetName);
                objSearch.setProbeSetName(strProbeSetName);
                iSearchIndex = m_arProbeSets.binarySearch(objSearch, 0);
                if (iSearchIndex != -1)
                {
                    CNProbeSet* p = m_arProbeSets.getAt(iSearchIndex);
                    tsv5->get(0, 1, &d); p->setFLD((float)d);
                    p->setValidFLDExists(true);
                }
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

/**
 * @brief Load the allele summary table file and set up the experiment array with experiment name.
 * @param const AffxString& - The name of the file to load.
 * @return int - The number of experiments loaded
 */
int CNLog2RatioData::loadExperiments(const AffxString& strFileName)
{
    getExperiments()->deleteAll();
    unsigned int uiColCount = 0;
    if (affx::File5_File::isHdf5file(strFileName))
    {
        try
        {
            AffxString str;
            affx::File5_File file5;
            file5.open(strFileName, affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = NULL;
            affx::File5_Tsv* tsv5 = NULL;
            AffxString strAnalysisName = m_pEngine->getOpt("set-analysis-name");
            if (strFileName == m_pEngine->getOpt("reference-file"))
            {
                group5 = file5.openGroup(strAnalysisName, affx::FILE5_OPEN);
                tsv5 = group5->openTsv(strAnalysisName + ".plier.summary", affx::FILE5_OPEN);
            }
            else
            {
                tsv5 = file5.openTsv(strAnalysisName + ".plier.summary", affx::FILE5_OPEN);
            }
            uiColCount = tsv5->getColumnCount(0);
            for (unsigned int uiColIndex = 1; (uiColIndex < uiColCount); uiColIndex++)
            {
                CNExperiment* pobjExperiment = new CNExperiment;
                tsv5->getColumnName(0, uiColIndex, &str);
                pobjExperiment->setIndex(uiColIndex - 1);
                pobjExperiment->setExperimentName(str);
                getExperiments()->add(pobjExperiment);
            }
            tsv5->close();
            delete tsv5;
            if (group5 != NULL)
            {
                group5->close();
                delete group5;
            }
            file5.close();
            uiColCount--;
        } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    }
    else
    {
      affx::TsvFile tsv;
      tsv.m_optAutoTrim = true;
      if (tsv.openTable(strFileName) == affx::TSV_OK) {
        std::string col;
        std::string columnHeader("probeset_id");
        while ( (tsv.nextLevel(0) == affx::TSV_OK) &&
                (tsv.get(0,0,col) == affx::TSV_OK) &&
                ( col != columnHeader) ) { }
                
        if ( col == columnHeader ){
          uiColCount = tsv.getColumnCount(0) - 1;
          int cidx = 1;
          for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++) {
            CNExperiment* pobjExperiment = new CNExperiment;
            pobjExperiment->setIndex(uiColIndex);
            tsv.get(0,cidx++, col);
            pobjExperiment->setExperimentName(col);
            getExperiments()->add(pobjExperiment);
          }
        }
        tsv.clear();
      } else {throw(Except("Cannot open file: " + strFileName));}
    }
    return uiColCount;
}


/**
 * @brief Load externally computed experiment genders.
 * @param const AffxString& - The name of the file to load from
 */

bool CNLog2RatioData::loadGenderOverrideFile(const AffxString& strFileName)
{
        //std::vector< std::pair< std:string, std::string> > vGenders;
        std::vector< pair<string,string>  > vGenders;
        std::pair<std::string, std::string> genderCall;
        affx::TsvFile tsv;
        tsv.m_optAutoTrim = true;
        if(tsv.open(strFileName) != affx::TSV_OK)
                Err::errAbort("Unable to open gender-override-file: "  + m_pEngine->getOpt("gender-override-file"));


        // We first load the gender information in to a vector from the gender-override-file.
        while(tsv.nextLevel(0)==TSV_OK){
        string inputString1;
        string inputString2;

                if(tsv.get(0,"experiment",inputString1)==TSV_OK)
                {
                        if(tsv.get(0,"gender",inputString2)==TSV_OK)
                        {

                               std::pair<string,string> pairGender(inputString1, inputString2);
                               vGenders.push_back(pairGender);
                        } else
                        {
                                Err::errAbort("Unable to parse gender-override-file: missing column: gender.");
                        }

                } else
                {
                        Err::errAbort("Unable to parse gender-override-file: missing column: gender.");
                }

        }
        //  Now iterate over the experiments that have been loaded and see if we have gender information from the gender-override-file.
        //  If not we issue a warning that no information has been found.
        int iExperimentCount = getExperiments()->getCount();
        for (int iExperimentIndex = 0; (iExperimentIndex < iExperimentCount); iExperimentIndex++)
        {
               CNExperiment* pobjExperiment = getExperiments()->getAt(iExperimentIndex);
               AffxString strExperimentName = pobjExperiment->getExperimentName();

               vector< pair <string, string> >::iterator beginIter= vGenders.begin();
               vector< pair <string, string> >::iterator endIter= vGenders.end();
               bool foundFlag=false;
               for( ; beginIter!=endIter; beginIter++)
               {
                       if(strExperimentName == (*beginIter).first)
                       {
                               if((*beginIter).second == "M")
                               {
                                   pobjExperiment->setRawIntensityRatioGender(affx::Male);
                                   Verbose::out(1, "CNLog2RatioData::loadGenderOverrideFile(...)\t" + pobjExperiment->getExperimentName() + "\t" + ((pobjExperiment->getRawIntensityRatioGender() == "male") ? "XY" : "XX"));

                               }
                               if((*beginIter).second == "F")
                               {
                                       pobjExperiment->setRawIntensityRatioGender(affx::Female);
                                       Verbose::out(1, "CNLog2RatioData::loadGenderOverrideFile(...)\t" + pobjExperiment->getExperimentName() + "\t" + ((pobjExperiment->getRawIntensityRatioGender() == "male") ? "XY" : "XX"));
                               }
                               if((*beginIter).second != "M" && (*beginIter).second !="F")
                               {
                                       Err::errAbort("Invalid gender found in gender-override-file for sample: " + strExperimentName);
                               }
                               foundFlag=true;
                        }
               }

               if(!foundFlag)
               {

                       Err::errAbort("No gender information found in gender-override-file for sample: " + strExperimentName);
               }


        }

        return true;
}


/**
 * @brief Load the QC metrics from the report file
 * @param const AffxString& - The name of the file to load from
 */
void CNLog2RatioData::loadQCMetrics(const AffxString& strFileName)
{
    if (strFileName == "") {return;}
    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
    
    if (tsv.open(strFileName) == affx::TSV_OK) {
      std::string strColumnName;
      for (int iColumn = 1; (iColumn < tsv.getColumnCount(0)); iColumn++){
        strColumnName = tsv.getColumnName(0, iColumn);
        CNExperiment::getQCMetricColumnNames()->add(new AffxString(strColumnName));
      }
      while (tsv.nextLevel(0) == affx::TSV_OK)  {
        int cidx = 0;
        AffxString strExperimentName;
        std::string strColumnValue;
        tsv.get(0,cidx++, strExperimentName );
        for (int iExperimentIndex = 0; (iExperimentIndex < getExperiments()->getCount()); iExperimentIndex++)  {
          CNExperiment* pobjExperiment = getExperiments()->getAt(iExperimentIndex);
          if (Fs::basename(pobjExperiment->getExperimentName()) == Fs::basename(strExperimentName))  {
            for (; cidx < tsv.getColumnCount(0); cidx++) {
              tsv.get(0,cidx, strColumnValue);
              pobjExperiment->getQCMetricColumnValues()->add(new AffxString(strColumnValue));
            }
          }
        }
      }
      tsv.clear();
    } else {throw(Except("Cannot open file: " + strFileName));}
}

/**
 * @brief Load the QC metric from the genotype report file
 * @param const AffxString& - The name of the file to load from
 */
void CNLog2RatioData::loadGenotypeReport(const AffxString& strFileName)
{
    if (strFileName == "") {throw(Except("Genotype report file not specied."));}
    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
    
    if (tsv.open(strFileName) == affx::TSV_OK)  {
      int iColumnNumber = 0;
      int ret = tsv.cname2cidx(0,std::string("computed_gender"));
      
      if ( ret > TSV_OK ) {
        iColumnNumber = ret;
      }
      else {
        throw(Except("Genotype report file must contain the computed_gender column."));
      }
      int iMaleCount = 0;
      int iFemaleCount = 0;
      int iUnknownCount = 0;
      while (tsv.nextLevel(0) == affx::TSV_OK){
        int cidx = 0;
        std::string col;
        tsv.get(0, cidx++, col);
        AffxString strExperimentName = col;
        for (int iExperimentIndex = 0; (iExperimentIndex < getExperiments()->getCount()); iExperimentIndex++) {
          CNExperiment* pobjExperiment = getExperiments()->getAt(iExperimentIndex);
          AffxString str1 = Fs::basename(pobjExperiment->getExperimentName());
          str1 = str1.toLowerCase();
          AffxString str2 = Fs::basename(strExperimentName);
          str2 = str2.toLowerCase();
          if (str1 == str2) {
            AffxString strColumnValue;
            tsv.get(0,iColumnNumber, strColumnValue);
            pobjExperiment->setRawIntensityRatioGenderFromString(strColumnValue);
            strColumnValue = strColumnValue.toLowerCase();
            if (strColumnValue == "male") 
            {
                pobjExperiment->setY(true); 
                pobjExperiment->setCNCallGenderFromString("male"); 
                iMaleCount++;
                }
            else if (strColumnValue == "female") 
            {
                pobjExperiment->setXX(true); 
                pobjExperiment->setCNCallGenderFromString("female"); 
                iFemaleCount++;
            }
            else  {
              Verbose::out(1, "WARNING: Cannot determine sex for experiment: " + pobjExperiment->getExperimentName());
              iUnknownCount++;
            }
          }
        }
      }
      tsv.clear();
      if (m_pEngine->getOptBool("create-reference")) {
        Verbose::out(1, "Copynumber Reference contains " + ::getInt(iFemaleCount) + " female, " + ::getInt(iMaleCount) + " male, and " + ::getInt(iUnknownCount) + " unknown genders.");
      }
    } else {throw(Except("Cannot open file: " + strFileName));}
}

/**
 * @brief Dump the QC metrics into the copynumber report file
 * @param const AffxString& - The name of the file to dump to
 */
void CNLog2RatioData::dumpGenotypeReport(const AffxString& strFileName)
{
  AffxString strOutputFileName = Fs::join(m_pEngine->getOpt("out-dir"),
                                          m_pEngine->getOpt("set-analysis-name") + ".CopyNumber.Report.txt");
  affx::TsvFile tsvOut;

  if (strFileName != "")  {
    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
      
    if (tsv.open(strFileName) == affx::TSV_OK) {
      tsvOut.addHeader("version-to-report", m_pEngine->getOpt("version-to-report"));
      TsvFileHeaderLine * tsvHeader = NULL;
      for ( tsv.headersBegin(); (tsvHeader = tsv.nextHeaderPtr()) != NULL; ) {
        if ( tsvHeader->m_key.empty() ) {
          tsvOut.addHeaderComment(tsvHeader->m_value);
        }
        else {
          if ( tsvHeader->m_key == "guid") {
            tsvOut.addHeader("guid" , affxutil::Guid::GenerateNewGuid());
          }
          else {
            tsvOut.addHeader(tsvHeader->m_key,tsvHeader->m_value);
          }
        }
      }
      const std::string additionalHeaders[] = {
        "CN2Gender-ChrX-mean", "CN2Gender-ChrY-mean",
        "CN2Gender", "MedianAutosomeMedian",
        "MAPD", "iqr", "all_probeset_rle_mean",
        "gc-correction-size", "sample-median-cn-state",
        "sample-hom-frequency", "sample-het-frequency",
        "median-raw-intensity", "SNPQC", "RawSNPQC", "XX",
        "Y", "AntigenomicRatio", "waviness-seg-count" , "waviness-sd"
      };
      int additionalHeadersCount = sizeof( additionalHeaders) / sizeof(additionalHeaders[0]);
      std::string cname;
      int cidx = 0;
      for (; cidx < tsv.getColumnCount(0); cidx++ ) {
        tsv.cidx2cname( 0, cidx, cname );
        tsvOut.defineColumn(0, cidx, cname);
      }
      int totalColumns = tsv.getColumnCount(0) + additionalHeadersCount;
      for ( ;cidx < totalColumns; cidx++ ) {
        cname = additionalHeaders[cidx - tsv.getColumnCount(0)];
        tsvOut.defineColumn(0, cidx, cname);
      }
      tsvOut.writeTsv_v1(strOutputFileName);
        
      while (tsv.nextLevel(0) == affx::TSV_OK) {
        std::string col;
        for ( cidx = 0; cidx < tsv.getColumnCount(0); cidx++ ) {
          tsv.get(0,cidx, col);
          tsvOut.set(0,cidx, col);
        }

        tsv.get(0,0,col);
        AffxString strExperimentName = col;
        for (int iExperimentIndex = 0; (iExperimentIndex < getExperiments()->getCount()); iExperimentIndex++) {
          CNExperiment* pobjExperiment = getExperiments()->getAt(iExperimentIndex);
          AffxString str1 = Fs::basename(pobjExperiment->getExperimentName());
          str1 = str1.toLowerCase();
          AffxString str2 = Fs::basename(strExperimentName);
          str2 = str2.toLowerCase();
          if (str1 == str2) {
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getChrXMean(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getChrYMean(), 6));
            tsvOut.set(0,cidx++,pobjExperiment->getCNCallGender());
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianAutosomeMedian(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMadDiffCN(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getIqr(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMeanAbsRle(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getGCCorrectionMetric(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianCnState(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getHomFrequency(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getHetFrequency(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianRawIntensity(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getSNPQC(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getRawSNPQC(), 6));
            tsvOut.set(0,cidx++,::getInt(pobjExperiment->hasXX() ? 1 : 0));
            tsvOut.set(0,cidx++,::getInt(pobjExperiment->hasY() ? 1 : 0));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getAntigenomicRatio(), 6));
            tsvOut.set(0,cidx++,::getInt(pobjExperiment->getWavinessSegCountTotal()));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getWavinessSd(), 6));
          }
        }
        tsvOut.writeLevel(0);
      }
      tsv.clear();
      tsvOut.close();
      tsvOut.clear();
    } else {throw(Except("Cannot open file: " + strFileName));}
  }
  else {
    tsvOut.addHeader("version-to-report" , m_pEngine->getOpt("version-to-report"));
    tsvOut.addHeader("guid" , affxutil::Guid::GenerateNewGuid());
    const std::string headers[] = {
      "cel_files","CN2Gender-ChrX-mean","CN2Gender-ChrY-mean",
      "CN2Gender","MedianAutosomeMedian","MAPD","iqr",
      "all_probeset_rle_mean","gc-correction-size",
      "sample-median-cn-state","sample-hom-frequency",
      "sample-het-frequency","median-raw-intensity","SNPQC",
      "RawSNPQC","XX","Y","AntigenomicRatio","waviness-seg-count",
      "waviness-sd",
    };
    int headerCount = sizeof(headers) / sizeof(headers[0]);
    for ( int i = 0; i < headerCount; i++ ) {
      tsvOut.defineColumn(0,i,headers[i]);
    }
    tsvOut.writeTsv_v1(strOutputFileName);
    
    for (int iExperimentIndex = 0; (iExperimentIndex < getExperiments()->getCount()); iExperimentIndex++){
      int cidx =0;
      CNExperiment* pobjExperiment = getExperiments()->getAt(iExperimentIndex);
      tsvOut.set(0,cidx++, Fs::basename(pobjExperiment->getExperimentName()));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getChrXMean(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getChrYMean(), 6));
      tsvOut.set(0,cidx++,pobjExperiment->getCNCallGender());
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianAutosomeMedian(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMadDiffCN(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getIqr(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMeanAbsRle(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getGCCorrectionMetric(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianCnState(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getHomFrequency(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getHetFrequency(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianRawIntensity(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getSNPQC(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getRawSNPQC(), 6));
      tsvOut.set(0,cidx++,::getInt(pobjExperiment->hasXX() ? 1 : 0));
      tsvOut.set(0,cidx++,::getInt(pobjExperiment->hasY() ? 1 : 0));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getAntigenomicRatio(), 6));
      tsvOut.set(0,cidx++,::getInt(pobjExperiment->getWavinessSegCountTotal()));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getWavinessSd(), 6));
      tsvOut.writeLevel(0);
    }
    tsvOut.close();
    tsvOut.clear();
  }
}

/**
 * @brief Dump the QC metrics into the copynumber report file
 * @param const AffxString& - The name of the file to dump to
 */
void CNLog2RatioData::dumpCyto2Report(const AffxString& strFileName)
{
  wstring wDateTime = DateTime::GetCurrentDateTime().ToString();
  AffxString dateTime = StringUtils::ConvertWCSToMBS(wDateTime);
  dateTime.replace(':', '-');
  AffxString strOutputFileName = Fs::join(m_pEngine->getOpt("out-dir"),
                                          dateTime + "." + m_pEngine->getOpt("set-analysis-name") + ".CopyNumber.Report.txt");
  affx::TsvFile tsvOut;
  if (strFileName != "")  {
    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
      
    if (tsv.open(strFileName) == affx::TSV_OK) {
      tsvOut.addHeader("version-to-report", m_pEngine->getOpt("version-to-report"));
      TsvFileHeaderLine * tsvHeader = NULL;
      for ( tsv.headersBegin(); (tsvHeader = tsv.nextHeaderPtr()) != NULL; ) {
        if ( tsvHeader->m_key.empty() ) {
          tsvOut.addHeaderComment(tsvHeader->m_value);
        }
        else {
          if ( tsvHeader->m_key == "guid") {
            tsvOut.addHeader("guid" , affxutil::Guid::GenerateNewGuid());
          }
          else {
            tsvOut.addHeader(tsvHeader->m_key,tsvHeader->m_value);
          }
        }
      }
      const std::string additionalHeaders[] = {
        "Y-Gender" , "MedianAutosomeMedian" , "MAPD" , "iqr" ,
        "all_probeset_rle_mean" , "gc-correction-size" ,
        "sample-median-cn-state" , "median-raw-intensity" ,
        "call-rate", "SNPQC" , "RawSNPQC" , "AntigenomicRatio" ,
        "waviness-seg-count" , "waviness-sd"
      };
      int additionalHeadersCount = sizeof( additionalHeaders) / sizeof(additionalHeaders[0]);
      std::string cname;
      int cidx = 0;
      for (; cidx < tsv.getColumnCount(0); cidx++ ) {
        tsv.cidx2cname( 0, cidx, cname );
        tsvOut.defineColumn(0, cidx, cname);
      }
      int totalColumns = tsv.getColumnCount(0) + additionalHeadersCount;
      for ( ;cidx < totalColumns; cidx++ ) {
        cname = additionalHeaders[cidx - tsv.getColumnCount(0)];
        tsvOut.defineColumn(0, cidx, cname);
      }
      tsvOut.writeTsv_v1(strOutputFileName);

      while (tsv.nextLevel(0) == affx::TSV_OK) {

        std::string col;
        for ( cidx = 0; cidx < tsv.getColumnCount(0); cidx++ ) {
          tsv.get(0,cidx, col);
          tsvOut.set(0,cidx, col);
        }
        tsv.get(0,0, col);
        AffxString strExperimentName = col;

        for (int iExperimentIndex = 0; (iExperimentIndex < getExperiments()->getCount()); iExperimentIndex++)  {
          CNExperiment* pobjExperiment = getExperiments()->getAt(iExperimentIndex);
          AffxString str1 = Fs::basename(pobjExperiment->getExperimentName());
          str1 = str1.toLowerCase();
          AffxString str2 = Fs::basename(strExperimentName);
          str2 = str2.toLowerCase();
          if (str1 == str2) {
            tsvOut.set(0,cidx++,pobjExperiment->getCNCallGender());
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianAutosomeMedian(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMadDiffCN(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getIqr(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMeanAbsRle(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getGCCorrectionMetric(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianCnState(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianRawIntensity(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getCallRate(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getSNPQC(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getRawSNPQC(), 6));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getAntigenomicRatio(), 6));
            tsvOut.set(0,cidx++,::getInt(pobjExperiment->getWavinessSegCountTotal()));
            tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getWavinessSd(), 6));
          }
        }
        tsvOut.writeLevel(0);
      }
      tsv.clear();
      tsvOut.close();
      tsvOut.clear();
    } else {throw(Except("Cannot open file: " + strFileName));}
  }
  else  {
    tsvOut.addHeader("version-to-report" , m_pEngine->getOpt("version-to-report"));
    tsvOut.addHeader("guid" , affxutil::Guid::GenerateNewGuid());
    const std::string headers[] = {
      "cel_files" , "Y-Gender" , "MedianAutosomeMedian" ,
      "MAPD" , "iqr" , "all_probeset_rle_mean" ,
      "gc-correction-size" , "sample-median-cn-state" , "median-raw-intensity" ,
      "call-rate", "SNPQC" , "RawSNPQC" , "AntigenomicRatio" , "waviness-seg-count" , "waviness-sd"
    };
    int headerCount = sizeof(headers) / sizeof(headers[0]);
    for ( int i = 0; i < headerCount; i++ ) {
      tsvOut.defineColumn(0,i,headers[i]);
    }
    tsvOut.writeTsv_v1(strOutputFileName);

    for (int iExperimentIndex = 0; (iExperimentIndex < getExperiments()->getCount()); iExperimentIndex++)  {
      int cidx = 0;
      CNExperiment* pobjExperiment = getExperiments()->getAt(iExperimentIndex);
      tsvOut.set(0,cidx++,Fs::basename(pobjExperiment->getExperimentName()));
      tsvOut.set(0,cidx++,pobjExperiment->getCNCallGender());
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianAutosomeMedian(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMadDiffCN(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getIqr(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMeanAbsRle(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getGCCorrectionMetric(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianCnState(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getMedianRawIntensity(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getCallRate(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getSNPQC(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getRawSNPQC(), 6));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getAntigenomicRatio(), 6));
      tsvOut.set(0,cidx++,::getInt(pobjExperiment->getWavinessSegCountTotal()));
      tsvOut.set(0,cidx++,::getDouble(pobjExperiment->getWavinessSd(), 6));
      tsvOut.writeLevel(0);
    }
    tsvOut.close();
    tsvOut.clear();
  }
}

/**
 * Load the allele summary table file and set up the A Allele signal matrix and the B Allele signal matrix.
 * The data is load by Probe set meaning the all the experiment data is loaded for a subset of the probe sets.
 * @param const AffxString& - The name of the file to load.
 * @param int - The probe set index to start from.
 * @param int - The number of probe sets to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadAlleleEstimatesByProbeSet(const AffxString& strFileName, int iProbeSetIndex, int iProbeSetsToProcess)
{
    if (affx::File5_File::isHdf5file(strFileName)) {return loadAlleleEstimatesByProbeSetHdf5(strFileName, iProbeSetIndex, iProbeSetsToProcess);}
    bool bSuccessful = false;
    affx::TsvFile tsv;
    unsigned int uiRowCount = m_arProbeSets.getCount();
    CNProbeSet objSearch;
    tsv.m_optAutoTrim = true;
    
    if (tsv.openTable(strFileName) == affx::TSV_OK) {
      bSuccessful = true;
      std::string col;
      std::string columnHeader("probeset_id");
      while ( (tsv.nextLevel(0) == affx::TSV_OK) &&
              (tsv.get(0,0,col) == affx::TSV_OK) &&
              ( col != columnHeader) ) { }
                
      if ( col == columnHeader ){
        unsigned int uiColCount = tsv.getColumnCount(0) - 1;
        int iDotCount = 10;
        Verbose::progressBegin(1, "Loading Allele Estimates file ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
        int iRowIndex = 0;
        int iProbeSetsProcessed = 0;
        while (tsv.nextLevel(0) == affx::TSV_OK)    {
          int cidx = 0;
          Verbose::progressStep(1);
          std::string col;
          tsv.get(0, cidx++, col);
          AffxString strProbeSetName = col;
          if ((strProbeSetName.startsWith("SNP_A-")) || (strProbeSetName.startsWith("AFFX-SNP_A-"))) {
            strProbeSetName = strProbeSetName.substring(0, (int)strProbeSetName.length() - 2);
          }
          bool bCn = ((strProbeSetName.startsWith("CN")) || (strProbeSetName.startsWith("chr")));
          objSearch.setProbeSetName(strProbeSetName);
          iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
          if ((iRowIndex == -1) || ((iRowIndex >= iProbeSetIndex) && (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), (int)uiRowCount)))){
            if (iRowIndex != -1) {iRowIndex -= iProbeSetIndex;}
            if (iRowIndex != -1) {
              for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++)  {
                if (m_pEngine->getOptBool("log2-input"))   {
                  tsv.get(0,cidx++, col);
                  m_mxAAlleleEstimates.set(uiColIndex, iRowIndex, (float)(exp(log(2.0) * AffxByteArray(col).parseDouble())));
                }
                else {
                  m_mxAAlleleEstimates.set(uiColIndex, iRowIndex, (float)AffxByteArray(col).parseDouble());
                }
                if (bCn) {m_mxAAlleleEstimates.set(uiColIndex, iRowIndex, (float)(m_mxAAlleleEstimates.get(uiColIndex, iRowIndex) / 2.0));}
              }
            }
            if ((strProbeSetName.startsWith("SNP_A-")) || (strProbeSetName.startsWith("AFFX-SNP_A-"))) {
              if ( tsv.nextLevel(0) == affx::TSV_OK ) {
                if (iRowIndex != -1)  {
                  int cidx2 = 1;
                  for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++) {
                    tsv.get(0,cidx2++, col);
                    if (m_pEngine->getOptBool("log2-input")) {
                      m_mxBAlleleEstimates.set(uiColIndex, iRowIndex, (float)(exp(log(2.0) * AffxByteArray(col).parseDouble())));
                    }
                    else {
                      m_mxBAlleleEstimates.set(uiColIndex, iRowIndex, (float)AffxByteArray(col).parseDouble());
                    }
                  }
                }
              }
            }
            else {
              if (iRowIndex != -1) {
                for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++)  {
                  m_mxBAlleleEstimates.set(uiColIndex, iRowIndex, m_mxAAlleleEstimates.get(uiColIndex, iRowIndex));
                }
              }
            }
            if (iRowIndex != -1) {iProbeSetsProcessed++; if (iProbeSetsProcessed == iProbeSetsToProcess) {break;}}
          }
        }
        Verbose::progressEnd(1, "Done.");
      }
      tsv.clear();
    } else {throw(Except("Cannot open file: " + strFileName));}
    return bSuccessful;
}

/**
 * Load the allele summary table file and set up the A Allele signal matrix and the B Allele signal matrix.
 * The data is load by Probe set meaning the all the experiment data is loaded for a subset of the probe sets.
 * @param const AffxString& - The name of the file to load.
 * @param int - The probe set index to start from.
 * @param int - The number of probe sets to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadAlleleEstimatesByProbeSetHdf5(const AffxString& strFileName, int iProbeSetIndex, int iProbeSetsToProcess)
{
    bool bSuccessful = false;
    unsigned int uiRowCount = m_arProbeSets.getCount();
    CNProbeSet objSearch;
    try
    {
        AffxString str;
        double d = 0;
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_OPEN_RO);
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;
        AffxString strAnalysisName = m_pEngine->getOpt("set-analysis-name");
        if (strFileName == m_pEngine->getOpt("reference-file"))
        {
            group5 = file5.openGroup(strAnalysisName, affx::FILE5_OPEN);
            tsv5 = group5->openTsv(strAnalysisName + ".plier.summary", affx::FILE5_OPEN);
        }
        else
        {
            tsv5 = file5.openTsv(strAnalysisName + ".plier.summary", affx::FILE5_OPEN);
        }
        bSuccessful = true;
        unsigned int uiColCount = tsv5->getColumnCount(0) - 1;
        int iDotCount = 10;
        Verbose::progressBegin(1, "Loading Allele Estimates file ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
        int iRowIndex = 0;
        int iProbeSetsProcessed = 0;
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            Verbose::progressStep(1);
            AffxString strProbeSetName; tsv5->get(0, 0, &strProbeSetName);
            if ((strProbeSetName.startsWith("SNP_A-")) || (strProbeSetName.startsWith("AFFX-SNP_A-")))
            {
                strProbeSetName = strProbeSetName.substring(0, (int)strProbeSetName.length() - 2);
            }
            bool bCn = ((strProbeSetName.startsWith("CN")) || (strProbeSetName.startsWith("chr")));
            objSearch.setProbeSetName(strProbeSetName);
            iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
            if ((iRowIndex == -1) || ((iRowIndex >= iProbeSetIndex) && (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), (int)uiRowCount))))
            {
                if (iRowIndex != -1) {iRowIndex -= iProbeSetIndex;}
                if (iRowIndex != -1)
                {
                    for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++)
                    {
                        tsv5->get(0, (uiColIndex + 1), &d);
                        if (m_pEngine->getOptBool("log2-input"))
                        {
                            m_mxAAlleleEstimates.set(uiColIndex, iRowIndex, (float)(exp(log(2.0) * d)));
                        }
                        else
                        {
                            m_mxAAlleleEstimates.set(uiColIndex, iRowIndex, (float)d);
                        }
                        if (bCn) {m_mxAAlleleEstimates.set(uiColIndex, iRowIndex, (float)(m_mxAAlleleEstimates.get(uiColIndex, iRowIndex) / 2.0));}
                    }
                }
                if ((strProbeSetName.startsWith("SNP_A-")) || (strProbeSetName.startsWith("AFFX-SNP_A-")))
                {
                    tsv5->nextLine();
                    if (iRowIndex != -1)
                    {
                        for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++)
                        {
                            tsv5->get(0, (uiColIndex + 1), &d);
                            if (m_pEngine->getOptBool("log2-input"))
                            {
                                m_mxBAlleleEstimates.set(uiColIndex, iRowIndex, (float)(exp(log(2.0) * d)));
                            }
                            else
                            {
                                m_mxBAlleleEstimates.set(uiColIndex, iRowIndex, (float)d);
                            }
                        }
                    }
                }
                else
                {
                    if (iRowIndex != -1)
                    {
                        for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++)
                        {
                            m_mxBAlleleEstimates.set(uiColIndex, iRowIndex, m_mxAAlleleEstimates.get(uiColIndex, iRowIndex));
                        }
                    }
                }
                if (iRowIndex != -1) {iProbeSetsProcessed++; if (iProbeSetsProcessed == iProbeSetsToProcess) {break;}}
            }
        }
        Verbose::progressEnd(1, "Done.");
        tsv5->close();
        delete tsv5;
        if (group5 != NULL)
        {
            group5->close();
            delete group5;
        }
        file5.close();
        bSuccessful = true;
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    return bSuccessful;
}

/**
 * Load the allele summary table file and set up the A Allele signal matrix and the B Allele signal matrix.
 * The data is load by experiment meaning the all the probe set data is loaded for a subset of the experiments.
 * @param const AffxString& - The name of the file to load.
 * @param int - The experiment index to start from.
 * @param int - The number of exeriments to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadAlleleEstimatesByExperiment(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess)
{
    if (affx::File5_File::isHdf5file(strFileName)) {return loadAlleleEstimatesByExperimentHdf5(strFileName, iExperimentIndex, iExperimentsToProcess);}
    Verbose::out(3, "CNLog2RatioData::loadAlleleEstimatesByExperiment");
    bool bSuccessful = false;

    affx::TsvFile tsv;
    unsigned int uiRowCount = m_arProbeSets.getCount();
    int iRowIndex = -1;
    CNProbeSet objSearch;
    tsv.m_optAutoTrim = true;
    
    if (tsv.openTable(strFileName) == affx::TSV_OK) {
      bSuccessful = true;
      std::string col;
      std::string columnHeader("probeset_id");
      while ( (tsv.nextLevel(0) == affx::TSV_OK) &&
              (tsv.get(0,0,col) == affx::TSV_OK) &&
              ( col != columnHeader) ) { }
                
      if ( col == columnHeader ){
        bSuccessful = true;
        unsigned int uiColCount = tsv.getColumnCount(0) - 1;
        getExperiments()->reserve(uiColCount);
        int iDotCount = 10;
        Verbose::progressBegin(1, "Loading Allele Estimates file ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
        while (tsv.nextLevel(0) == affx::TSV_OK) {
          Verbose::progressStep(1);
          tsv.get(0,0, col);
          AffxString strProbeSetName = col;
          if ((strProbeSetName.startsWith("SNP_A-")) || (strProbeSetName.startsWith("AFFX-SNP_A-"))) {
            strProbeSetName = strProbeSetName.substring(0, (int)strProbeSetName.length() - 2);
          }
          bool bCn = ((strProbeSetName.startsWith("CN")) || (strProbeSetName.startsWith("chr")));
          objSearch.setProbeSetName(strProbeSetName);
          iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
          if (iRowIndex != -1) {
            for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++){
              tsv.get(0, iColIndex -1, col);
              if (m_pEngine->getOptBool("log2-input"))  {
                m_mxAAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)(exp(log(2.0) * AffxByteArray(col).parseDouble())));
              }
              else {
                m_mxAAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float) AffxByteArray(col).parseDouble());
              }
              if (bCn) {m_mxAAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)(m_mxAAlleleEstimates.get((iColIndex - iExperimentIndex), iRowIndex) / 2.0));}
            }
          }
          if ((strProbeSetName.startsWith("SNP_A-")) || (strProbeSetName.startsWith("AFFX-SNP_A-")))  {
            if ( tsv.nextLevel(0) == affx::TSV_OK ) {
              if (iRowIndex != -1) {
                for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++){
                  tsv.get(0,iColIndex - 1, col);
                  if (m_pEngine->getOptBool("log2-input")) {
                    m_mxBAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)(exp(log(2.0) * AffxByteArray(col).parseDouble())));
                  }
                  else {
                    m_mxBAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)AffxByteArray(col).parseDouble());
                  }
                }
              }
            }
          }
          else {
            if (iRowIndex != -1)  {
              for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++){
                m_mxBAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, m_mxAAlleleEstimates.get((iColIndex - iExperimentIndex), iRowIndex));
              }
            }
          }
        }
      }
      Verbose::progressEnd(1, "Done.");
    } else {throw(Except("Cannot open file: " + strFileName));}
    tsv.clear();
    return bSuccessful;
}

/* Load the allele summary table file and set up the A Allele signal matrix and the B Allele signal matrix.
 * The data is load by experiment meaning the all the probe set data is loaded for a subset of the experiments.
 * @param const AffxString& - The name of the file to load.
 * @param int - The experiment index to start from.
 * @param int - The number of exeriments to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadAlleleEstimatesByExperimentHdf5(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess)
{
    Verbose::out(3, "CNLog2RatioData::loadAlleleEstimatesByExperimentHdf5");
    bool bSuccessful = false;
    unsigned int uiRowCount = m_arProbeSets.getCount();
    int iRowIndex = -1;
    CNProbeSet objSearch;
    try
    {
        AffxString str;
        double d = 0;
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_OPEN_RO);
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;
        AffxString strAnalysisName = m_pEngine->getOpt("set-analysis-name");
        if (strFileName == m_pEngine->getOpt("reference-file"))
        {
            group5 = file5.openGroup(strAnalysisName, affx::FILE5_OPEN);
            tsv5 = group5->openTsv(strAnalysisName + ".plier.summary", affx::FILE5_OPEN);
        }
        else
        {
            tsv5 = file5.openTsv(strAnalysisName + ".plier.summary", affx::FILE5_OPEN);
        }
        bSuccessful = true;
        unsigned int uiColCount = tsv5->getColumnCount(0) - 1;
        getExperiments()->reserve(uiColCount);
        int iDotCount = 10;
        Verbose::progressBegin(1, "Loading Allele Estimates file ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            Verbose::progressStep(1);
            AffxString strProbeSetName; tsv5->get(0, 0, &strProbeSetName);
            if ((strProbeSetName.startsWith("SNP_A-")) || (strProbeSetName.startsWith("AFFX-SNP_A-")))
            {
                strProbeSetName = strProbeSetName.substring(0, (int)strProbeSetName.length() - 2);
            }
            bool bCn = ((strProbeSetName.startsWith("CN")) || (strProbeSetName.startsWith("chr")));
            objSearch.setProbeSetName(strProbeSetName);
            iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
            if (iRowIndex != -1)
            {
                for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)
                {
                    tsv5->get(0, (iColIndex + 1), &d);
                    if (m_pEngine->getOptBool("log2-input"))
                    {
                        m_mxAAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)(exp(log(2.0) * d)));

                    }
                    else
                    {
                        m_mxAAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)d);
                    }
                    if (bCn) {m_mxAAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)(m_mxAAlleleEstimates.get((iColIndex - iExperimentIndex), iRowIndex) / 2.0));}
                }
            }
            if ((strProbeSetName.startsWith("SNP_A-")) || (strProbeSetName.startsWith("AFFX-SNP_A-")))
            {
                tsv5->nextLine();
                if (iRowIndex != -1)
                {
                    for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)
                    {
                        tsv5->get(0, (iColIndex + 1), &d);
                        if (m_pEngine->getOptBool("log2-input"))
                        {
                            m_mxBAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)(exp(log(2.0) * d)));
                        }
                        else
                        {
                            m_mxBAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, (float)d);
                        }
                    }
                }
            }
            else
            {
                if (iRowIndex != -1)
                {
                    for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)
                    {
                        m_mxBAlleleEstimates.set((iColIndex - iExperimentIndex), iRowIndex, m_mxAAlleleEstimates.get((iColIndex - iExperimentIndex), iRowIndex));
                    }
                }
            }
        }
        Verbose::progressEnd(1, "Done.");
        tsv5->close();
        delete tsv5;
        if (group5 != NULL)
        {
            group5->close();
            delete group5;
        }
        file5.close();
        bSuccessful = true;
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    return bSuccessful;
}

/**
 * Load the genotype call table file and set up the genotype call matrix.
 * The data is load by Probe set meaning the all the experiment data is loaded for a subset of the probe sets.
 * @param const AffxString& - The name of the file to load.
 * @param int - The probe set index to start from.
 * @param int - The number of probe sets to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadGenotypeCallsByProbeSet(const AffxString& strFileName, int iProbeSetIndex, int iProbeSetsToProcess)
{
    if (affx::File5_File::isHdf5file(strFileName)) {return loadGenotypeCallsByProbeSetHdf5(strFileName, iProbeSetIndex, iProbeSetsToProcess);}
    bool bSuccessful = false;
    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
    
    unsigned int uiRowCount = m_arProbeSets.getCount();
    CNProbeSet objSearch;
    tsv.openTable(strFileName);
    std::string col;
    std::string columnHeader("probeset_id");
    while ( (tsv.nextLevel(0) == affx::TSV_OK) &&
            (tsv.get(0,0,col) == affx::TSV_OK) &&
            ( col != columnHeader) ) { }
    
    if ( col == columnHeader ){
      bSuccessful = true;
      unsigned int uiColCount = tsv.getColumnCount(0) - 1;
      int cidx = 1;
      if ((int)uiColCount != getExperiments()->getCount()) {throw(Except("Genotype Calls table does not have the correct number of experiments."));}
      for (int iColIndex = 0; (iColIndex < (int)uiColCount); iColIndex++)  {
        AffxString strExperimentName = getExperiments()->getAt(iColIndex)->getExperimentName();
        int iFindIndex = strExperimentName.indexOf(".");
        if (iFindIndex != -1) {strExperimentName = strExperimentName.substring(0, iFindIndex);}
        tsv.get(0,cidx++, col);
        AffxString strColumnName = Fs::basename(col);
        iFindIndex = strColumnName.indexOf(".");
        if (iFindIndex != -1) {strColumnName = strColumnName.substring(0, iFindIndex);}
        if (strColumnName != strExperimentName) {throw(Except("Experiment name mismatch between Allele Summary table and Genotype Calls table: " + strExperimentName + " vs. " + strColumnName));}
      }
      for (int iRowIndex = iProbeSetIndex; (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), (int)uiRowCount)); iRowIndex++)  {
        for (int iColIndex = 0; (iColIndex < (int)uiColCount); iColIndex++) {
          m_mxGenotypeCalls.set(iColIndex, (iRowIndex - iProbeSetIndex), (char)-1);
        }
      }
      int iDotCount = 10;
      Verbose::progressBegin(1, "Loading Genotype Calls ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
      while (tsv.nextLevel(0) == affx::TSV_OK) {
        Verbose::progressStep(1);
        cidx = 0;
        tsv.get(0,cidx++, col);
        AffxString strProbeSetName = col;
        objSearch.setProbeSetName(strProbeSetName);
        int iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
        if ((iRowIndex == -1) || ((iRowIndex >= iProbeSetIndex) && (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), (int)uiRowCount)))){
          if (iRowIndex != -1) {
            for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++) {
              tsv.get(0,cidx++, col);
              m_mxGenotypeCalls.set(uiColIndex, (iRowIndex - iProbeSetIndex), (char)AffxByteArray(col).parseInt());
            }
          }
        }
      }
      Verbose::progressEnd(1, "Done.");
    }
    tsv.clear();
    return bSuccessful;
}

/**
 * Load the genotype call table file and set up the genotype call matrix.
 * The data is load by Probe set meaning the all the experiment data is loaded for a subset of the probe sets.
 * @param const AffxString& - The name of the file to load.
 * @param int - The probe set index to start from.
 * @param int - The number of probe sets to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadGenotypeCallsByProbeSetHdf5(const AffxString& strFileName, int iProbeSetIndex, int iProbeSetsToProcess)
{
    bool bSuccessful = false;
    unsigned int uiRowCount = m_arProbeSets.getCount();
    CNProbeSet objSearch;
    try
    {
        AffxString str;
        int i = 0;
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_OPEN_RO);
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;
        AffxString strAnalysisName = m_pEngine->getOpt("set-analysis-name");
        if (strFileName == m_pEngine->getOpt("reference-file"))
        {
            group5 = file5.openGroup(strAnalysisName, affx::FILE5_OPEN);
            tsv5 = group5->openTsv(strAnalysisName + ".calls", affx::FILE5_OPEN);
        }
        else
        {
            tsv5 = file5.openTsv(strAnalysisName + ".calls", affx::FILE5_OPEN);
        }
        bSuccessful = true;
        unsigned int uiColCount = tsv5->getColumnCount(0) - 1;
        if ((int)uiColCount != getExperiments()->getCount()) {throw(Except("Genotype Calls table does not have the correct number of experiments."));}
        for (int iColIndex = 0; (iColIndex < (int)uiColCount); iColIndex++)
        {
            AffxString strExperimentName = getExperiments()->getAt(iColIndex)->getExperimentName();
            int iFindIndex = strExperimentName.indexOf(".");
            if (iFindIndex != -1) {strExperimentName = strExperimentName.substring(0, iFindIndex);}
            AffxString strColumnName = tsv5->getColumnPtr(0, iColIndex + 1)->getColumnName();
            iFindIndex = strColumnName.indexOf(".");
            if (iFindIndex != -1) {strColumnName = strColumnName.substring(0, iFindIndex);}
            if (strColumnName != strExperimentName) {throw(Except("Experiment name mismatch between Allele Summary table and Genotype Calls table: " + strExperimentName + " vs. " + strColumnName));}
        }
        for (int iRowIndex = iProbeSetIndex; (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), (int)uiRowCount)); iRowIndex++)
        {
            for (int iColIndex = 0; (iColIndex < (int)uiColCount); iColIndex++)
            {
                m_mxGenotypeCalls.set(iColIndex, (iRowIndex - iProbeSetIndex), (char)-1);
            }
        }
        int iDotCount = 10;
        Verbose::progressBegin(1, "Loading Genotype Calls ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            Verbose::progressStep(1);
            AffxString strProbeSetName; tsv5->get(0, 0, &strProbeSetName);
            objSearch.setProbeSetName(strProbeSetName);
            int iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
            if ((iRowIndex == -1) || ((iRowIndex >= iProbeSetIndex) && (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), (int)uiRowCount))))
            {
                if (iRowIndex != -1)
                {
                    for (unsigned int uiColIndex = 0; (uiColIndex < uiColCount); uiColIndex++)
                    {
                        tsv5->get(0, (uiColIndex + 1), &i);
                        m_mxGenotypeCalls.set(uiColIndex, (iRowIndex - iProbeSetIndex), (char)i);
                    }
                }
            }
        }
        Verbose::progressEnd(1, "Done.");
        tsv5->close();
        delete tsv5;
        if (group5 != NULL)
        {
            group5->close();
            delete group5;
        }
        file5.close();
        bSuccessful = true;
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    return bSuccessful;
}

/**
 * Load the genotype call table file and set up the genotype call matrix.
 * The data is load by experiment meaning the all the probe set data is loaded for a subset of the experiments.
 * @param const AffxString& - The name of the file to load.
 * @param int - The experiment index to start from.
 * @param int - The number of exeriments to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadGenotypeCallsByExperiment(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess)
{
    if (affx::File5_File::isHdf5file(strFileName)) {return loadGenotypeCallsByExperimentHdf5(strFileName, iExperimentIndex, iExperimentsToProcess);}
    bool bSuccessful = false;
    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
    
    unsigned int uiRowCount = m_arProbeSets.getCount();
    CNProbeSet objSearch;
    std::string col;
    
    tsv.openTable(strFileName);
    bSuccessful = true;
    std::string columnHeader("probeset_id");
    while ( (tsv.nextLevel(0) == affx::TSV_OK) &&
            (tsv.get(0,0,col) == affx::TSV_OK) &&
            ( col != columnHeader) ) { }
    
    if ( col == columnHeader ){
      unsigned int uiColCount = tsv.getColumnCount(0) - 1;
      int cidx = 1;
      if ((int)uiColCount != getExperiments()->getCount()) {throw(Except("Genotype Calls table does not have the correct number of experiments."));}
      for (int iColIndex = 0; (iColIndex < (int)uiColCount); iColIndex++)  {
        AffxString strExperimentName = getExperiments()->getAt(iColIndex)->getExperimentName();
        int iFindIndex = strExperimentName.indexOf(".");
        if (iFindIndex != -1) {strExperimentName = strExperimentName.substring(0, iFindIndex);}
        tsv.get(0,cidx++, col);
        AffxString strColumnName = Fs::basename(col);
        iFindIndex = strColumnName.indexOf(".");
        if (iFindIndex != -1) {strColumnName = strColumnName.substring(0, iFindIndex);}
        if (strColumnName != strExperimentName) {throw(Except("Experiment name mismatch between Allele Summary table and Genotype Calls table: " + strExperimentName + " vs. " + strColumnName));}
      }
      for (int iRowIndex = 0; (iRowIndex < (int)uiRowCount); iRowIndex++) {
        for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++) {
          m_mxGenotypeCalls.set((iColIndex - iExperimentIndex), iRowIndex, (char)-1);
        }
      }
      int iDotCount = 10;
      Verbose::progressBegin(1, "Loading Genotype Calls ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
      while (tsv.nextLevel(0) == affx::TSV_OK) {
        Verbose::progressStep(1);
        int cidx = 0;
        tsv.get(0, cidx++, col);
        AffxString strProbeSetName = col;
        objSearch.setProbeSetName(strProbeSetName);
        int iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
        if (iRowIndex != -1){
          for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)  {
            tsv.get(0,iColIndex - 1, col);
            m_mxGenotypeCalls.set((iColIndex - iExperimentIndex), iRowIndex, (char)AffxByteArray(col).parseInt());
          }
        }
      }
      Verbose::progressEnd(1, "Done.");
    }
    tsv.clear();
    return bSuccessful;
}

/**
 * Load the genotype call table file and set up the genotype call matrix.
 * The data is load by experiment meaning the all the probe set data is loaded for a subset of the experiments.
 * @param const AffxString& - The name of the file to load.
 * @param int - The experiment index to start from.
 * @param int - The number of exeriments to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadGenotypeCallsByExperimentHdf5(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess)
{
    bool bSuccessful = false;
    unsigned int uiRowCount = m_arProbeSets.getCount();
    CNProbeSet objSearch;
    try
    {
        AffxString str;
        int i = 0;
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_OPEN_RO);
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;
        AffxString strAnalysisName = m_pEngine->getOpt("set-analysis-name");
        if (strFileName == m_pEngine->getOpt("reference-file"))
        {
            group5 = file5.openGroup(strAnalysisName, affx::FILE5_OPEN);
            tsv5 = group5->openTsv(strAnalysisName + ".calls", affx::FILE5_OPEN);
        }
        else
        {
            tsv5 = file5.openTsv(strAnalysisName + ".calls", affx::FILE5_OPEN);
        }
        bSuccessful = true;
        unsigned int uiColCount = tsv5->getColumnCount(0) - 1;
        if ((int)uiColCount != getExperiments()->getCount()) {throw(Except("Genotype Calls table does not have the correct number of experiments."));}
        for (int iColIndex = 0; (iColIndex < (int)uiColCount); iColIndex++)
        {
            AffxString strExperimentName = getExperiments()->getAt(iColIndex)->getExperimentName();
            int iFindIndex = strExperimentName.indexOf(".");
            if (iFindIndex != -1) {strExperimentName = strExperimentName.substring(0, iFindIndex);}
            AffxString strColumnName = tsv5->getColumnPtr(0, iColIndex + 1)->getColumnName();
            iFindIndex = strColumnName.indexOf(".");
            if (iFindIndex != -1) {strColumnName = strColumnName.substring(0, iFindIndex);}
            if (strColumnName != strExperimentName) {throw(Except("Experiment name mismatch between Allele Summary table and Genotype Calls table: " + strExperimentName + " vs. " + strColumnName));}
        }
        for (int iRowIndex = 0; (iRowIndex < (int)uiRowCount); iRowIndex++)
        {
            for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)
            {
                m_mxGenotypeCalls.set((iColIndex - iExperimentIndex), iRowIndex, (char)-1);
            }
        }
        int iDotCount = 10;
        Verbose::progressBegin(1, "Loading Genotype Calls ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            Verbose::progressStep(1);
            AffxString strProbeSetName; tsv5->get(0, 0, &strProbeSetName);
            objSearch.setProbeSetName(strProbeSetName);
            int iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
            if (iRowIndex != -1)
            {
                for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)
                {
                    tsv5->get(0, (iColIndex + 1), &i);
                    m_mxGenotypeCalls.set((iColIndex - iExperimentIndex), iRowIndex, (char)i);
                }
            }
        }
        Verbose::progressEnd(1, "Done.");
        tsv5->close();
        delete tsv5;
        if (group5 != NULL)
        {
            group5->close();
            delete group5;
        }
        file5.close();
        bSuccessful = true;
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    return bSuccessful;
}

/**
 * Load the genotype confidence table file and set up the genotype confidence matrix.
 * The data is load by experiment meaning the all the probe set data is loaded for a subset of the experiments.
 * @param const AffxString& - The name of the file to load.
 * @param int - The experiment index to start from.
 * @param int - The number of exeriments to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadGenotypeConfidencesByExperiment(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess)
{
  if (affx::File5_File::isHdf5file(strFileName)) {return loadGenotypeConfidencesByExperimentHdf5(strFileName, iExperimentIndex, iExperimentsToProcess);}
  bool bSuccessful = false;
  affx::TsvFile tsv;
  tsv.m_optAutoTrim = true;
  unsigned int uiRowCount = m_arProbeSets.getCount();
  int iRowIndex = -1;
  CNProbeSet objSearch;
  tsv.openTable(strFileName);
  bSuccessful = true;
  std::string col;
  std::string columnHeader("probeset_id");
  while ( (tsv.nextLevel(0) == affx::TSV_OK) &&
          (tsv.get(0,0,col) == affx::TSV_OK) &&
          ( col != columnHeader) ) { }
    
  if ( col == columnHeader ){
    unsigned int uiColCount = tsv.getColumnCount(0) - 1;
    int cidx = 1;
    if ((int)uiColCount != getExperiments()->getCount()) {throw(Except("Genotype Calls table does not have the correct number of experiments."));}
    for (int iColIndex = 0; (iColIndex < (int)uiColCount); iColIndex++){
      tsv.get(0, cidx++, col);
      AffxString strExperimentName = getExperiments()->getAt(iColIndex)->getExperimentName();
      int iFindIndex = strExperimentName.indexOf(".");
      if (iFindIndex != -1) {strExperimentName = strExperimentName.substring(0, iFindIndex);}
      AffxString strColumnName = Fs::basename(col);
      iFindIndex = strColumnName.indexOf(".");
      if (iFindIndex != -1) {strColumnName = strColumnName.substring(0, iFindIndex);}
      if (strColumnName != strExperimentName) {throw(Except("Experiment name mismatch between Allele Summary table and Genotype Calls table: " + strExperimentName + " vs. " + strColumnName));}
    }
    int iDotCount = 10;
    Verbose::progressBegin(1, "Loading Genotype Confidences ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
    while (tsv.nextLevel(0) == affx::TSV_OK) {
      Verbose::progressStep(1);
      int cidx = 0;
      tsv.get(0,cidx++, col);
      AffxString strProbeSetName = col;
      objSearch.setProbeSetName(strProbeSetName);
      iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
      if (iRowIndex != -1) {
        for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)  {
          tsv.get(0, iColIndex - 1, col);
          m_mxGenotypeConfidences.set((iColIndex - iExperimentIndex), iRowIndex, (float)AffxByteArray(col).parseDouble());
        }
      }
    }
    Verbose::progressEnd(1, "Done.");
  }
  tsv.clear();
  return bSuccessful;
}

/**
 * Load the genotype confidence table file and set up the genotype confidence matrix.
 * The data is load by experiment meaning the all the probe set data is loaded for a subset of the experiments.
 * @param const AffxString& - The name of the file to load.
 * @param int - The experiment index to start from.
 * @param int - The number of exeriments to process each time this function is called.
 * @return bool - true if successful
 */
bool CNLog2RatioData::loadGenotypeConfidencesByExperimentHdf5(const AffxString& strFileName, int iExperimentIndex, int iExperimentsToProcess)
{
    bool bSuccessful = false;
    unsigned int uiRowCount = m_arProbeSets.getCount();
    CNProbeSet objSearch;
    try
    {
        AffxString str;
        double d = 0;
        affx::File5_File file5;
        file5.open(strFileName, affx::FILE5_OPEN_RO);
        affx::File5_Group* group5 = NULL;
        affx::File5_Tsv* tsv5 = NULL;
        AffxString strAnalysisName = m_pEngine->getOpt("set-analysis-name");
        if (strFileName == m_pEngine->getOpt("reference-file"))
        {
            group5 = file5.openGroup(strAnalysisName, affx::FILE5_OPEN);
            tsv5 = group5->openTsv(strAnalysisName + ".confidences", affx::FILE5_OPEN);
        }
        else
        {
            tsv5 = file5.openTsv(strAnalysisName + ".confidences", affx::FILE5_OPEN);
        }
        bSuccessful = true;
        unsigned int uiColCount = tsv5->getColumnCount(0) - 1;
        if ((int)uiColCount != getExperiments()->getCount()) {throw(Except("Genotype Confidences table does not have the correct number of experiments."));}
        for (int iColIndex = 0; (iColIndex < (int)uiColCount); iColIndex++)
        {
            AffxString strExperimentName = getExperiments()->getAt(iColIndex)->getExperimentName();
            int iFindIndex = strExperimentName.indexOf(".");
            if (iFindIndex != -1) {strExperimentName = strExperimentName.substring(0, iFindIndex);}
            AffxString strColumnName = tsv5->getColumnPtr(0, iColIndex + 1)->getColumnName();
            iFindIndex = strColumnName.indexOf(".");
            if (iFindIndex != -1) {strColumnName = strColumnName.substring(0, iFindIndex);}
            if (strColumnName != strExperimentName) {throw(Except("Experiment name mismatch between Allele Summary table and Genotype Confidences table: " + strExperimentName + " vs. " + strColumnName));}
        }
        for (int iRowIndex = 0; (iRowIndex < (int)uiRowCount); iRowIndex++)
        {
            for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)
            {
                m_mxGenotypeConfidences.set((iColIndex - iExperimentIndex), iRowIndex, 0);
            }
        }
        int iDotCount = 10;
        Verbose::progressBegin(1, "Loading Genotype Confidences ", iDotCount, Max((int)(uiRowCount / iDotCount), 1), uiRowCount);
        while (tsv5->nextLine() == affx::FILE5_OK)
        {
            Verbose::progressStep(1);
            AffxString strProbeSetName; tsv5->get(0, 0, &strProbeSetName);
            objSearch.setProbeSetName(strProbeSetName);
            int iRowIndex = m_arProbeSets.binarySearch(objSearch, 0);
            if (iRowIndex != -1)
            {
                for (int iColIndex = iExperimentIndex; (iColIndex < Min((iExperimentIndex + iExperimentsToProcess), (int)uiColCount)); iColIndex++)
                {
                    tsv5->get(0, (iColIndex + 1), &d);
                    m_mxGenotypeConfidences.set((iColIndex - iExperimentIndex), iRowIndex, (float)d);
                }
            }
        }
        Verbose::progressEnd(1, "Done.");
        tsv5->close();
        delete tsv5;
        if (group5 != NULL)
        {
            group5->close();
            delete group5;
        }
        file5.close();
        bSuccessful = true;
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
    return bSuccessful;
}

/**
 * Calculate the median singal for the log2 ratio model.
 * The data is processed by probes set meaning the all the experiment data is processed for a subset of the probe sets.
 * Note that the median signal is A + B.
 * @param int - The probe set index to start from.
 * @param int - The number of probe sets to process each time this function is called.
 */
void CNLog2RatioData::calculateMedianSignal(int iProbeSetIndex, int iProbeSetsToProcess)
{
    Verbose::out(3, "CNLog2RatioData::calculateMedianSignal");
    AffxMultiDimensionalArray<float> vSignals(getExperiments()->getCount());
    for (int iRowIndex = iProbeSetIndex; (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), m_arProbeSets.getCount())); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = m_arProbeSets.getAt(iRowIndex);
        for (int iColIndex = 0; (iColIndex < getExperiments()->getCount()); iColIndex++)
        {
            vSignals.set(iColIndex, (float)getSnpSignalEstimate(iColIndex, (iRowIndex - iProbeSetIndex), pobjProbeSet));
        }
        pobjProbeSet->setMedianSignal((float)::roundDouble(vSignals.median(), 10));
    }
}

/**
 * Calculate the median singals for the log2 ratio model.
 * The data is processed by probes set meaning the all the experiment data is processed for a subset of the probe sets.
 * Note that the median signals broken out by genotype is A - B.
 * @param int - The probe set index to start from.
 * @param int - The number of probe sets to process each time this function is called.
 */
void CNLog2RatioData::calculateMedianSignalsPerGenotype(int iProbeSetIndex, int iProbeSetsToProcess)
{
    Verbose::out(3, "CNLog2RatioData::calculateMedianSignalsPerGenotype");
    int iExperimentCount = getExperiments()->getCount();
    AffxMultiDimensionalArray<float> vAASignals(iExperimentCount);
    AffxMultiDimensionalArray<float> vABSignals(iExperimentCount);
    AffxMultiDimensionalArray<float> vBBSignals(iExperimentCount);
    for (int iRowIndex = iProbeSetIndex; (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), m_arProbeSets.getCount())); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = m_arProbeSets.getAt(iRowIndex);
        int iAAIndex = 0;
        int iABIndex = 0;
        int iBBIndex = 0;
        for (int iColIndex = 0; (iColIndex < getExperiments()->getCount()); iColIndex++)
        {
            CNExperiment* pobjExperiment = getExperiments()->getAt(iColIndex);
            if ((pobjProbeSet->getChromosome() != m_pEngine->getOptInt("xChromosome")) || (pobjExperiment->hasXX()))
            {
                float fAllelicDifference = (float)getAllelicDifference(iColIndex, 0, (iRowIndex - iProbeSetIndex));
                switch (m_mxGenotypeCalls.get(iColIndex, (iRowIndex - iProbeSetIndex)))
                {
                case 0: vAASignals.set(iAAIndex, fAllelicDifference); iAAIndex++; break;
                case 1: vABSignals.set(iABIndex, fAllelicDifference); iABIndex++; break;
                case 2: vBBSignals.set(iBBIndex, fAllelicDifference); iBBIndex++; break;
                }
            }
        }
        pobjProbeSet->setAAMedianSignal((float)::roundDouble(vAASignals.median(iAAIndex), 10));
        pobjProbeSet->setABMedianSignal((float)::roundDouble(vABSignals.median(iABIndex), 10));
        pobjProbeSet->setBBMedianSignal((float)::roundDouble(vBBSignals.median(iBBIndex), 10));
    }
}

/**
 * Calculate the median singals for the log2 ratio model.
 * The data is processed by probes set meaning the all the experiment data is processed for a subset of the probe sets.
 * Note that the median signal is A + B. The hasXX and hasY is set in the calculateXY function.
 * @param int - The probe set index to start from.
 * @param int - The number of probe sets to process each time this function is called.
 */
void CNLog2RatioData::calculateXYMedianSignal(int iProbeSetIndex, int iProbeSetsToProcess)
{
    Verbose::out(3, "CNLog2RatioData::calculateXYMedianSignal");
    int iXXCount = 0;
    int iYCount = 0;
    for (int iColIndex = 0; (iColIndex < getExperiments()->getCount()); iColIndex++)
    {
        CNExperiment* pobjExperiment = getExperiments()->getAt(iColIndex);
        if (pobjExperiment->hasXX()) {iXXCount++;}
        if (pobjExperiment->hasY()) {iYCount++;}
    }
    AffxMultiDimensionalArray<float> vXX(iXXCount);
    AffxMultiDimensionalArray<float> vY(iYCount);
    for (int iRowIndex = iProbeSetIndex; (iRowIndex < Min((iProbeSetIndex + iProbeSetsToProcess), m_arProbeSets.getCount())); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = m_arProbeSets.getAt(iRowIndex);
        int iXXIndex = 0;
        int iYIndex = 0;
        for (int iColIndex = 0; (iColIndex < getExperiments()->getCount()); iColIndex++)
        {
            CNExperiment* pobjExperiment = getExperiments()->getAt(iColIndex);
            if (pobjExperiment->hasXX())
            {
                vXX.set(iXXIndex, (float)getSnpSignalEstimate(iColIndex, (iRowIndex - iProbeSetIndex), pobjProbeSet));
                iXXIndex++;
            }
            if (pobjExperiment->hasY())
            {
                vY.set(iYIndex, (float)getSnpSignalEstimate(iColIndex, (iRowIndex - iProbeSetIndex), pobjProbeSet));
                iYIndex++;
            }
        }
        // The values are rounded here to insure that the values stored in the model will match the values as caluclated using the same inputs.
        pobjProbeSet->setXXMedianSignal((float)::roundDouble(vXX.median(), 10));
        pobjProbeSet->setYMedianSignal((float)::roundDouble(vY.median(), 10));
    }
}

/**
 * Write the model file in text format. All the median signals asssociated with the probe sets are written out.
 * Values are rounded to insure that the stored values will match the calculated values when using the same input.
 * @param const AffxString& - The file name to write to.
 * @return bool - true if successful
 */
bool CNLog2RatioData::writeModelFile(const AffxString& strFileName)
{
    Verbose::out(3, "CNLog2RatioData::writeModelFile");
    if ( !Fs::dirExists(Fs::dirname(strFileName)) ) {
      Fs::mkdirPath(Fs::dirname(strFileName), false);
    }
    bool bSuccessful = false;
    affx::TsvFile tsvOut;
    
    AffxByteArray baLine;
    int iDotCount = 10;
    if (m_pEngine->isOptDefined("reference-text-output") && m_pEngine->getOptBool("reference-text-output"))  {
      // write out the options used.
      tsvOut.addHeader("guid" ,  affxutil::Guid::GenerateNewGuid());

      // dump options.
      tsvOut.addHeader("affymetrix-algorithm-param-apt-opt-celFileCount", m_pEngine->getOpt("sample-count"));
      tsvOut.addHeader("affymetrix-algorithm-param-apt-analysis-name" , m_pEngine->getOpt("set-analysis-name"));
      tsvOut.addHeader( "affymetrix-algorithm-param-apt-opt-adapter-normalization" ,
                        (m_pEngine->getOptBool("adapter-type-normalization") ? "true" : "false"));

      // dump options
      vector<string> optionNames;
      m_pEngine->getOptionNames(optionNames,1);
      for(int i=0; i< optionNames.size(); i++) {
        std::string name = optionNames[i];
        std::vector<std::string> vals = m_pEngine->getOptVector(name,1);
        if (vals.size() > 1)   {
          for(int j=0; j< vals.size(); j++) {
            std::string val = vals[j];
            tsvOut.addHeader("affymetrix-algorithm-param-apt-opt-" + name,  val);
          }
        }
        else {
          std::string val = m_pEngine->getOpt(name,1);
          //aw: hack because old versions of the code use "opt-set-analysis-name" instead of "state-analysis-name" to load the global prefix
          if(name == "set-analysis-name") { val = m_pEngine->getOpt("set-analysis-name"); }
          tsvOut.addHeader("affymetrix-algorithm-param-apt-opt-" + name , val);
        }
      }
      // dump state.
      std::vector<std::string> stateNames;
      m_pEngine->getOptionNames(stateNames);
      for(int i=0; i< stateNames.size(); i++) {
        std::string name = stateNames[i];
        std::vector<std::string> vals = m_pEngine->getOptVector(name);
        if (vals.size() > 1) {
          for(int j=0; j< vals.size(); j++) {
            std::string val = vals[j];
            tsvOut.addHeader(std::string("affymetrix-algorithm-param-apt-state-") + name ,  val);
          }
        }
        else {
          std::string val = m_pEngine->getOpt(name);
          tsvOut.addHeader(std::string("affymetrix-algorithm-param-apt-state-") + name , val);
        }
      }

      // write out the reference data.
      const std::string headers[] = {
        "probeset_id" , "MedianSignal" , "XXMedianSignal" , "YMedianSignal" ,
        "AAMedianSignal" , "ABMedianSignal" , "BBMedianSignal"
      };
      int headersCount = sizeof( headers ) / sizeof( headers[0] );
      for ( int i = 0; i < headersCount ; i++ ) {
        tsvOut.defineColumn(0, i, headers[i]);
      }
      tsvOut.writeTsv_v1(strFileName);
      
      Verbose::progressBegin(1, "Writing Model file ", iDotCount, Max((int)(m_arProbeSets.getCount() / iDotCount), 1), m_arProbeSets.getCount());
      for (int iRowIndex = 0; (iRowIndex < m_arProbeSets.getCount()); iRowIndex++)  {
        Verbose::progressStep(1);
        int cidx = 0;
        CNProbeSet* pobjProbeSet = m_arProbeSets.getAt(iRowIndex);
        tsvOut.set(0, cidx++, pobjProbeSet->getProbeSetName());
        tsvOut.set(0, cidx++, ::getDouble(pobjProbeSet->getMedianSignal(), 10));
        tsvOut.set(0, cidx++, ::getDouble(pobjProbeSet->getXXMedianSignal(), 10));
        tsvOut.set(0, cidx++, ::getDouble(pobjProbeSet->getYMedianSignal(), 10));
        tsvOut.set(0, cidx++, ::getDouble(pobjProbeSet->getAAMedianSignal(), 10));
        tsvOut.set(0, cidx++, ::getDouble(pobjProbeSet->getABMedianSignal(), 10));
        tsvOut.set(0, cidx++, ::getDouble(pobjProbeSet->getBBMedianSignal(), 10));
        tsvOut.writeLevel(0);
      }
      Verbose::progressEnd(1, "Done.");
      tsvOut.close();
      tsvOut.clear();
      bSuccessful = true;
    }
    else
    {
        try
        {
            AffxString str;
            affx::File5_File file5;
            bool bFileExists = Fs::fileExists(strFileName);
            
      // @todo I think this could be "affx::FILE5_OPEN||affx::FILE5_CREATE" -jhg
            file5.open(strFileName, (bFileExists ? affx::FILE5_OPEN : affx::FILE5_CREATE));
            affx::File5_Group* group5 = file5.openGroup("CopyNumber", affx::FILE5_REPLACE);
            affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_REPLACE);
            tsv5->defineColumn(0, 0, "Parameter", affx::FILE5_DTYPE_STRING, 1024);
            tsv5->set_string(0, 0, "#%guid=" + affxutil::Guid::GenerateNewGuid()); tsv5->writeLevel(0);

            // dump options.
            tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-opt-celFileCount=" + m_pEngine->getOpt("sample-count")); tsv5->writeLevel(0);
            tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-analysis-name=" + m_pEngine->getOpt("set-analysis-name")); tsv5->writeLevel(0);
            tsv5->set_string(0, 0, std::string() + "#%affymetrix-algorithm-param-apt-opt-adapter-normalization=" + (m_pEngine->getOptBool("adapter-type-normalization") ? "true" : "false")); tsv5->writeLevel(0);

            // dump options
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
                        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-opt-" + name + "=" + val); tsv5->writeLevel(0);
                    }
                }
                else
                {
                    std::string val = m_pEngine->getOpt(name,1);
                    //aw: hack because old versions of the code use "opt-set-analysis-name" instead of "state-analysis-name" to load the global prefix
                    if(name == "set-analysis-name") { val = m_pEngine->getOpt("set-analysis-name"); }
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-opt-" + name + "=" + val); tsv5->writeLevel(0);
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
                        tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-state-" + name + "=" + val); tsv5->writeLevel(0);
                    }
                }
                else
                {
                    std::string val = m_pEngine->getOpt(name);
                    tsv5->set_string(0, 0, "#%affymetrix-algorithm-param-apt-state-" + name + "=" + val); tsv5->writeLevel(0);
                }
            }

            tsv5->close();
            delete tsv5;

            int iMaxProbeSetNameLength = 0;
            for (int iRowIndex = 0; (iRowIndex < m_arProbeSets.getCount()); iRowIndex++)
            {
                CNProbeSet* pobjProbeSet = m_arProbeSets.getAt(iRowIndex);
                iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
            }
            Verbose::progressBegin(1, "Writing Model file ", iDotCount, Max((int)(m_arProbeSets.getCount() / iDotCount), 1), m_arProbeSets.getCount());
            tsv5 = group5->openTsv("Reference", affx::FILE5_REPLACE);
            tsv5->defineColumn(0, 0, "probeset_id", affx::FILE5_DTYPE_STRING, iMaxProbeSetNameLength);
            tsv5->defineColumn(0, 1, "MedianSignal", affx::FILE5_DTYPE_DOUBLE, 0);
            tsv5->defineColumn(0, 2, "XXMedianSignal", affx::FILE5_DTYPE_DOUBLE, 0);
            tsv5->defineColumn(0, 3, "YMedianSignal", affx::FILE5_DTYPE_DOUBLE, 0);
            tsv5->defineColumn(0, 4, "AAMedianSignal", affx::FILE5_DTYPE_DOUBLE, 0);
            tsv5->defineColumn(0, 5, "ABMedianSignal", affx::FILE5_DTYPE_DOUBLE, 0);
            tsv5->defineColumn(0, 6, "BBMedianSignal", affx::FILE5_DTYPE_DOUBLE, 0);
            tsv5->defineColumn(0, 7, "Chromosome", affx::FILE5_DTYPE_INT, 0);
            tsv5->defineColumn(0, 8, "Position", affx::FILE5_DTYPE_INT, 0);
            for (int iRowIndex = 0; (iRowIndex < m_arProbeSets.getCount()); iRowIndex++)
            {
                Verbose::progressStep(1);
                CNProbeSet* pobjProbeSet = m_arProbeSets.getAt(iRowIndex);
                tsv5->set_string(0, 0, pobjProbeSet->getProbeSetName());
                tsv5->set_d(0, 1, pobjProbeSet->getMedianSignal());
                tsv5->set_d(0, 2, pobjProbeSet->getXXMedianSignal());
                tsv5->set_d(0, 3, pobjProbeSet->getYMedianSignal());
                tsv5->set_d(0, 4, pobjProbeSet->getAAMedianSignal());
                tsv5->set_d(0, 5, pobjProbeSet->getABMedianSignal());
                tsv5->set_d(0, 6, pobjProbeSet->getBBMedianSignal());
                tsv5->set_i(0, 7, pobjProbeSet->getChromosome());
                tsv5->set_i(0, 8, pobjProbeSet->getPosition());
                tsv5->writeLevel(0);
            }
            Verbose::progressEnd(1, "Done.");
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

/**
 * Write the model file in text format. All the median signals asssociated with the probe sets are written out.
 * Values are rounded to insure that the stored values will match the calculated values when using the same input.
 * @param const AffxString& - The file name to write to.
 * @return bool - true if successful
 */
void CNLog2RatioData::writeCyto2ModelFile(const AffxString& strFileName)
{
    Verbose::out(3, "CNLog2RatioData::writeCyto2ModelFile");
    if ( !Fs::dirExists(Fs::dirname(strFileName)) ) {
      Fs::mkdirPath(Fs::dirname(strFileName), false);
    }
    AffxByteArray baLine;
    int iDotCount = 10;
    try
    {
        AffxString str;
        affx::File5_File file5;
        bool bFileExists = Fs::fileExists(strFileName);

        // @todo I think this could be "affx::FILE5_OPEN||affx::FILE5_CREATE" -jhg
        file5.open(strFileName, (bFileExists ? affx::FILE5_OPEN : affx::FILE5_CREATE));
        affx::File5_Group* group5 = file5.openGroup("Cyto2", affx::FILE5_REPLACE);
        affx::File5_Tsv* tsv5 = group5->openTsv("Parameters", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "Parameter", affx::FILE5_DTYPE_STRING, 1024);
        tsv5->set_string(0, 0, "#%guid=" + affxutil::Guid::GenerateNewGuid()); tsv5->writeLevel(0);

        // dump options
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
        for (int iRowIndex = 0; (iRowIndex < m_arProbeSets.getCount()); iRowIndex++)
        {
            CNProbeSet* pobjProbeSet = m_arProbeSets.getAt(iRowIndex);
            iMaxProbeSetNameLength = Max(iMaxProbeSetNameLength, (int)pobjProbeSet->getProbeSetName().length());
        }
        Verbose::progressBegin(1, "Writing Model file ", iDotCount, Max((int)(m_arProbeSets.getCount() / iDotCount), 1), m_arProbeSets.getCount());
        tsv5 = group5->openTsv("MedianSignals", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "probeset_id", affx::FILE5_DTYPE_STRING, iMaxProbeSetNameLength);
        tsv5->defineColumn(0, 1, "MedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 2, "XXMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 3, "YMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 4, "AAMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 5, "ABMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 6, "BBMedianSignal", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 7, "Chromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 8, "Position", affx::FILE5_DTYPE_INT);
        for (int iRowIndex = 0; (iRowIndex < m_arProbeSets.getCount()); iRowIndex++)
        {
            Verbose::progressStep(1);
            CNProbeSet* pobjProbeSet = m_arProbeSets.getAt(iRowIndex);
            tsv5->set_string(0, 0, pobjProbeSet->getProbeSetName());
            tsv5->set_d(0, 1, pobjProbeSet->getMedianSignal());
            tsv5->set_d(0, 2, pobjProbeSet->getXXMedianSignal());
            tsv5->set_d(0, 3, pobjProbeSet->getYMedianSignal());
            tsv5->set_d(0, 4, pobjProbeSet->getAAMedianSignal());
            tsv5->set_d(0, 5, pobjProbeSet->getABMedianSignal());
            tsv5->set_d(0, 6, pobjProbeSet->getBBMedianSignal());
            tsv5->set_i(0, 7, pobjProbeSet->getChromosome());
            tsv5->set_i(0, 8, pobjProbeSet->getPosition());
            tsv5->writeLevel(0);
        }
        Verbose::progressEnd(1, "Done.");
        tsv5->close();
        delete tsv5;
        group5->close();
        delete group5;
        file5.close();
    } catch(...) {throw(Except("Cannot open file: " + strFileName));}
}

/**
 * @brief Return the signal estimate for the SNP
 * @param int - The column index
 * @param int - The row idex
 * @return double - The signal estimate for the SNP
 */
double CNLog2RatioData::getSnpSignalEstimate(int iColIndex, int iRowIndex, CNProbeSet* p)
{
    double d = (m_mxAAlleleEstimates.get(iColIndex, iRowIndex) + m_mxBAlleleEstimates.get(iColIndex, iRowIndex));
    if (d == 0) {Err::errAbort("Zero SNP signal estimate found. CEL file may be truncated.\t" + p->getProbeSetName());}
    if (d != d) {Err::errAbort("NaN SNP signal estimate found. CEL file may be truncated.");}
    if (!Util::isFinite(d)) {Err::errAbort("Non finite SNP signal estimate found. CEL file may be truncated.");}
    return d;
}

/**
 * @brief Return the allelic difference for the SNP
 * @param int - The column index
 * @param int - The experiment index
 * @param int - The row idex
 * @return double - The allelic difference for the SNP
 */
double CNLog2RatioData::getAllelicDifference(int iColIndex, int iExperimentIndex, int iRowIndex)
{
    double d = (m_mxAAlleleEstimates.get((iColIndex - iExperimentIndex), iRowIndex) - m_mxBAlleleEstimates.get((iColIndex - iExperimentIndex), iRowIndex));
    if (!Util::isFinite(d)) {d = 0.0;}
    return d;
}

void CNLog2RatioData::loadProbeSets()
{
    AffxArray<AffxString> arRestrictList;
    bool bRestrictList = (m_pEngine->getOpt("probeset-ids") != "");
    if (bRestrictList)
    {
        loadProbeSetNamesFromRestrictList(arRestrictList);
    }
    AffxString strReferenceFileName = m_pEngine->getOpt("reference-file");
    affx::File5_File file5;
    affx::File5_Group* group5 = NULL;
    affx::File5_Tsv* tsv5 = NULL;

    file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
    if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
    {
        group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
        tsv5 = group5->openTsv("MedianSignals", affx::FILE5_OPEN);
    }
    else
    {
        group5 = file5.openGroup("CopyNumber", affx::FILE5_OPEN);
        tsv5 = group5->openTsv("Reference", affx::FILE5_OPEN);
    }
    unsigned int uiCount = 0;
    AffxString str;
    int iSearchIndex = -1;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        uiCount++;
    }
    tsv5->close();
    delete tsv5;

    m_arProbeSets.nullAll();
    m_arProbeSetsAlternateSort.nullAll();
    m_vProbeSets.clear();
    m_vProbeSets.allocate(uiCount);
    m_arProbeSets.reserve(uiCount);

    int iInputValue=0;
    double d = 0;
    if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
    {
        tsv5 = group5->openTsv("MedianSignals", affx::FILE5_OPEN);
    }
    else
    {
        tsv5 = group5->openTsv("Reference", affx::FILE5_OPEN);
    }
    unsigned int uiIndex = 0;
    while (tsv5->nextLine() == affx::FILE5_OK)
    {
        CNProbeSet *pset = m_vProbeSets.getAt(uiIndex);
        tsv5->get(0, 0, &str);
        if (bRestrictList)
        {
            iSearchIndex = arRestrictList.binarySearch(str, 0);
            if (iSearchIndex == -1) {pset->setProcess(false);}
        }
        pset->setProbeSetName(str);
        tsv5->get(0, 1, &d); pset->setMedianSignal(d);
        tsv5->get(0, 2, &d); pset->setXXMedianSignal(d);
        tsv5->get(0, 3, &d); pset->setYMedianSignal(d);
        tsv5->get(0, 4, &d); pset->setAAMedianSignal(d);
        tsv5->get(0, 5, &d); pset->setABMedianSignal(d);
        tsv5->get(0, 6, &d); pset->setBBMedianSignal(d);
        tsv5->get(0, 7, &iInputValue); pset->setChromosome(iInputValue);
        m_arProbeSets.add(pset);
        m_arProbeSetsAlternateSort.add(pset);
        uiIndex++;
    }
    tsv5->close();
    delete tsv5;

    group5->close();
    delete group5;
    file5.close();

    m_arProbeSets.quickSort(0); // By ProbeSetName
    m_arProbeSetsAlternateSort.quickSort(1); // By Chromosome, Position

    // Warn if requested probe set not found.
    CNProbeSet objSearch;
    if (arRestrictList.getCount() > 0)
    {
        for (int iIndex = 0; (iIndex < arRestrictList.getCount()); iIndex++)
        {
            AffxString str = *arRestrictList.getAt(iIndex);
            objSearch.setProbeSetName(str);
            int iSearchIndex = m_arProbeSets.binarySearch(objSearch, 0);
            if (iSearchIndex  == -1)
            {
                Verbose::out(1, "WARNING: ProbeSet specified in probeset-ids file will not appear in output data. ProbeSet: " + str);
            }
        }
    }
    arRestrictList.deleteAll();
}

void CNLog2RatioData::loadAnnotation(bool bNewProbeSets)
{
    if (bNewProbeSets)
    {
        m_arProbeSets.nullAll();
        m_arProbeSetsAlternateSort.nullAll();
    }
    else {loadProbeSets();}

    Annotation().loadAnnotationForCopynumber(*m_pEngine, bNewProbeSets, m_vProbeSets, m_arProbeSets);

    if (bNewProbeSets)
    {
        m_arProbeSets.reserve(m_vProbeSets.size());
        m_arProbeSetsAlternateSort.reserve(m_vProbeSets.size());
        for (unsigned int iIndex = 0; (iIndex < m_vProbeSets.size()); iIndex++)
        {
            CNProbeSet *pset = m_vProbeSets.getAt(iIndex);
            if (pset->getChromosome() != 0)
            {
                m_arProbeSets.add(pset);
                m_arProbeSetsAlternateSort.add(pset);
            }
        }
        m_arProbeSets.quickSort(0);
        m_arProbeSetsAlternateSort.quickSort(0);

        readCovariatesFile();

        for (int i = 1; (i < m_arProbeSets.getCount()); i++)
        {
            CNProbeSet* pPrev = m_arProbeSets.getAt(i - 1);
            CNProbeSet* pNext = m_arProbeSets.getAt(i);
            if (pPrev->getProbeSetName() == pNext->getProbeSetName())
            {
                Verbose::out(1, "WARNING: Duplicate probe set exists in annotation-file. The duplicate is being removed: " + pPrev->getProbeSetName());
                m_arProbeSets.removeAt(i);
                m_arProbeSetsAlternateSort.removeAt(i);
                i--;
            }
        }
        m_arProbeSetsAlternateSort.quickSort(1);
    }
    for (int i = 0; (i < m_arProbeSets.getCount()); i++)
    {
        CNProbeSet* p = m_arProbeSets.getAt(i);
        if ((p->getChromosome() == 0) || (p->getPosition() == 0))
        {
            p->setChromosome(255);
            p->setPosition(0);
            p->setGCContent(0);
            p->setReplicateCount(1);
            p->setPseudoAutosomalRegion(false);
            p->setStyAdapterCode(false);
            p->setNspAdapterCode(false);
        }
    }
    overrideGcContent();

   // NOTE: The set up of the GC correction bins was moved to the CNAnalysisMethodLog2Ratio module to allow
   // a subset of probe sets to be binned rather than the complete set as was original done at this location.

    Verbose::out(1, ::getInt(m_arProbeSets.getCount()) + "\tProbeSets to process.");
}

void CNLog2RatioData::overrideGcContent()
{
  if (!m_pEngine->isOptDefined("gc-content-override-file")) {return;}
  AffxString strGcContentOverrideFileName = m_pEngine->getOpt("gc-content-override-file");
  if (strGcContentOverrideFileName != "")  {
    CNProbeSet objSearch;
    affx::TsvFile tsv;
    tsv.m_optAutoTrim = true;
    std::string col;

    if (tsv.open(strGcContentOverrideFileName) == affx::TSV_OK)  {
      while (tsv.nextLevel(0) == affx::TSV_OK)  {
        int cidx = 0;
        tsv.get(0, cidx++, col);
        AffxString strProbeSetName = col;
        objSearch.setProbeSetName(strProbeSetName);
        int iSearchIndex = m_arProbeSets.binarySearch(objSearch, 0);
        if (iSearchIndex != -1) {
          tsv.get(0, cidx++, col);
          CNProbeSet* p = m_arProbeSets.getAt(iSearchIndex);
          p->setGCContent((float)AffxByteArray(col).parseDouble());
        }
      }
      tsv.clear();
    }
    else {throw(Except("Cannot open file: " + strGcContentOverrideFileName));}
  }
}

// This method is located in CNLog2RatioData because at this point we have a sorted
// CNProbeSetArray (m_arProbeSets) which is needed in order to do binarySearches based on probe set name
void CNLog2RatioData::readCovariatesFile()
{
    // Both annotation-file and covariate-file are required for reading the covariate-file to make any sense.
    if ((!m_pEngine->isOptDefined("annotation-file")) || (m_pEngine->getOpt("annotation-file") == "")) {return ;}
    if ((m_pEngine->isOptDefined("covariates-file")) && (m_pEngine->getOpt("covariates-file") != ""))
    {
        // Ensure that space has been allocated for the covariates.
        if (m_arProbeSets.size() > 0 && m_arProbeSets[0]->getNumCovariates() == 0)
        {
            Verbose::warn(1, "covariates-file option is ignored. annotation-file does not support covariates.");
            return;
        }

        try
        {
            // Open and read the file.
            affx::File5_File file5;
            file5.open(m_pEngine->getOpt("covariates-file"), affx::FILE5_OPEN_RO);
            affx::File5_Group* group5 = file5.openGroup("Covariates", affx::FILE5_OPEN);
            affx::File5_Tsv* tsv5 = group5->openTsv("Covariate-Values", affx::FILE5_OPEN);

            // Filter the names of covariates to those that we would expect to find in the covariates file (as opposed to the annotation db file).
            vector<pair<int,int> > columnIdxToCovariateIdx;
            for (map<string, int>::iterator ii = CovariateParams::m_allCovariateMap.begin(); ii != CovariateParams::m_allCovariateMap.end(); ++ii)
            {
                if (CovariateParams::isReservedCovariateName(ii->first) == false)
                {
                    // Build a vector of indexes (column index in the file, covariate index in the array).
                    int colIdx = tsv5->getColumnIdx(0, ii->first);
                    if (colIdx < 0) Err::errAbort("Covariate [" + ii->first + "] missing from file - " + m_pEngine->getOpt("covariates-file"));

                    columnIdxToCovariateIdx.push_back(make_pair(colIdx, CovariateParams::mapCovariateNameToIndex(ii->first)));
                }
            }

            // For each row in the file.
                // Find matching probe set.
                // Loop through all unreserved column and assign covariates.
            int iSearchIndex = -1;
            CNProbeSet objSearch;
            while (tsv5->nextLine() == affx::FILE5_OK)
            {
                AffxString strProbeSetName;
                tsv5->get(0, 0, &strProbeSetName);
                objSearch.setProbeSetName(strProbeSetName);
                iSearchIndex = m_arProbeSets.binarySearch(objSearch, 0);
                if (iSearchIndex != -1)
                {
                    CNProbeSet* p = m_arProbeSets.getAt(iSearchIndex);
                    for (vector<pair<int,int> >::iterator ii = columnIdxToCovariateIdx.begin(); ii != columnIdxToCovariateIdx.end(); ++ii)
                    {
                        float d;
                        tsv5->get(0, ii->first, &d);
                        p->setCovariateValue(ii->second, (float)d);
                    }
                }
            }
            tsv5->close();
            delete tsv5;
            group5->close();
            delete group5;
            file5.close();
        } catch(...) {throw(Except("Cannot open file: " + m_pEngine->getOpt("covariates-file")));}
    }
}
