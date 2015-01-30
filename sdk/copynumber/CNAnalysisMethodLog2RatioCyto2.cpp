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
 * @file CNAnalysisMethodLog2RatioCyto2.cpp
 *
 * @brief This file contains the CNAnalysisMethodLog2RatioCyto2 class members.
 */
#include "copynumber/CNAnalysisMethodLog2RatioCyto2.h"
//
#include "copynumber/CNAnalysisMethodFactory.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/AffxStatistics.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodLog2RatioCyto2::explainSelf()
{
    CNAnalysisMethodLog2RatioCyto2 obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodLog2RatioCyto2::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)
  SelfDoc::Opt opt2 = {"gc-correction", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA", "Log2RatioCyto2 GC Correction"};
  opts.push_back(opt2);

  SelfDoc::Opt opt3 = {"median-autosome-median-normalization", SelfDoc::Opt::Boolean, "true", "true", "NA", "NA", "Log2RatioCyto2 Median Autosmome Median Normalization"};
  opts.push_back(opt3);

  SelfDoc::Opt opt4 = {"median-smooth-marker-count", SelfDoc::Opt::Integer, "5", "5", "3", "NA", "Log2Ratio Median Smooth marker count."};
  opts.push_back(opt4);

  SelfDoc::Opt opt5 = {"trim-high", SelfDoc::Opt::Double, "2.0", "2.0", "0", "NA", "High trim value for Log2Ratios."};
  opts.push_back(opt5);

  SelfDoc::Opt opt6 = {"trim-low", SelfDoc::Opt::Double, "-2.5", "-2.5", "NA", "0", "Low trim value for Log2Ratios."};
  opts.push_back(opt6);
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
SelfCreate* CNAnalysisMethodLog2RatioCyto2::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodLog2RatioCyto2* pMethod = new CNAnalysisMethodLog2RatioCyto2();
    std::string strPrefix = getPrefix();

    pMethod->m_bGCCorrection = setupBoolParameter("gc-correction", strPrefix, params, doc, opts);
    pMethod->m_bMedianAutosomeMedianNormalization = setupBoolParameter("median-autosome-median-normalization", strPrefix, params, doc, opts);
    pMethod->m_iMedianSmoothWindowSize = setupIntParameter("median-smooth-marker-count", strPrefix, params, doc, opts);
    pMethod->m_dTrimHigh = setupDoubleParameter("trim-high", strPrefix, params, doc, opts);
    pMethod->m_dTrimLow = setupDoubleParameter("trim-low", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodLog2RatioCyto2::run()
{
//    Verbose::out(1, "CNAnalysisMethodLog2RatioCyto2::run(...) start");
    isSetup();
    m_dYTarget = m_pEngine->getOptDouble("yTarget");
    Verbose::progressBegin(1, "CNAnalysisMethodLog2RatioCyto2::run(...) ", 5, 1, 5);

    determineLocalProbeSets();

    calculateLog2Ratios(true);
    trim();

        if (m_pEngine->getOptBool("keep-intermediate-data"))
        {
            writeLog2Ratios("original");
            writeLog2RatiosTsv("original", getProbeSets());
        }


    applyWaveCorrection();
    Verbose::progressStep(1);

    if (m_bGCCorrection)
    {
        GCCorrection(getExperiment());
                if(m_pEngine->getOptBool("keep-intermediate-data"))
                {
                        writeLog2Ratios("afterGCCorrection");
                        writeLog2RatiosTsv("afterGCCorrection", getProbeSets());

                }
    }

    Verbose::progressStep(1);
    if (m_bMedianAutosomeMedianNormalization)
    {
        calculateMedianAutosomeMedian(getExperiment());
        Verbose::progressStep(1);
        normalizeLog2Ratios(getExperiment());
    }
    Verbose::progressStep(1);
    calculateQCMetrics(getExperiment());
    getProbeSets()->calculateLog2RatioMedianSmooth(m_iMedianSmoothWindowSize);
    Verbose::progressEnd(1, "Done");
//    Verbose::out(1, "CNAnalysisMethodLog2RatioCyto2::run(...) end");
}

void CNAnalysisMethodLog2RatioCyto2::applyWaveCorrection()
{
    if (!m_pEngine->isOptDefined("cyto2")) {return;}
    if (!m_pEngine->getOptBool("cyto2")) {return;}
    if (!m_pEngine->isOptDefined("wave-correction-log2ratio-adjustment-method")) {return;}
    if (m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") == "none") {return;}
    if (m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method") == "") {return;}

    CNAnalysisMethodFactory amFactory;
    CNAnalysisMethod* am = NULL;
    am = amFactory.CNAnalysisMethodForString(m_pEngine->getOpt("wave-correction-log2ratio-adjustment-method"));
    am->setEngine(m_pEngine);
    am->setup(*getExperiment(), *(CNAnalysisMethod::getProbeSets()));
    am->run();
    delete am;
}
