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
 * @file CNAnalysisMethodCNGender.cpp
 *
 * @brief This file contains the CNAnalysisMethodCNGender class members.
 */
#include "calvin_files/utils/src/StringUtils.h"
#include "copynumber/CNAnalysisMethodCNGender.h"
#include "util/AffxStatistics.h"
#include "util/Fs.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodCNGender::explainSelf()
{
    CNAnalysisMethodCNGender obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodCNGender::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)
  SelfDoc::Opt opt1 = {"male-chrX-lower-threshold", SelfDoc::Opt::Double, "0.8", "0.8", "0", "NA", "Male ChrX Lower Threshold"};
  SelfDoc::Opt opt2 = {"male-chrX-upper-threshold", SelfDoc::Opt::Double, "1.3", "1.3", "0", "NA", "Male ChrX Upper Threshold"};
  SelfDoc::Opt opt3 = {"male-chrY-lower-threshold", SelfDoc::Opt::Double, "0.8", "0.8", "0", "NA", "Male ChrY Lower Threshold"};
  SelfDoc::Opt opt4 = {"male-chrY-upper-threshold", SelfDoc::Opt::Double, "1.2", "1.2", "0", "NA", "Male ChrY Upper Threshold"};

  SelfDoc::Opt opt5 = {"female-chrX-lower-threshold", SelfDoc::Opt::Double, "1.9", "1.9", "0", "NA", "Female ChrX Lower Threshold"};
  SelfDoc::Opt opt6 = {"female-chrX-upper-threshold", SelfDoc::Opt::Double, "2.1", "2.1", "0", "NA", "Female ChrX Upper Threshold"};
  SelfDoc::Opt opt7 = {"female-chrY-lower-threshold", SelfDoc::Opt::Double, "0", "0", "0", "NA", "Female ChrY Lower Threshold"};
  SelfDoc::Opt opt8 = {"female-chrY-upper-threshold", SelfDoc::Opt::Double, "0.4", "0.4", "0", "NA", "Female ChrY Upper Threshold"};

  SelfDoc::Opt opt9 = {"mapd-threshold", SelfDoc::Opt::Double, "0.5", "0.5", "0", "1", "MAPD Threshold"};

  opts.push_back(opt1);
  opts.push_back(opt2);
  opts.push_back(opt3);
  opts.push_back(opt4);
  opts.push_back(opt5);
  opts.push_back(opt6);
  opts.push_back(opt7);
  opts.push_back(opt8);
  opts.push_back(opt9);

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
SelfCreate* CNAnalysisMethodCNGender::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodCNGender* pMethod = new CNAnalysisMethodCNGender();
    std::string strPrefix = getPrefix();

    pMethod->m_dMaleChrXLowerThreshold = setupDoubleParameter("male-chrX-lower-threshold", strPrefix, params, doc, opts);
    pMethod->m_dMaleChrXUpperThreshold = setupDoubleParameter("male-chrX-upper-threshold", strPrefix, params, doc, opts);
    pMethod->m_dMaleChrYLowerThreshold = setupDoubleParameter("male-chrY-lower-threshold", strPrefix, params, doc, opts);
    pMethod->m_dMaleChrYUpperThreshold = setupDoubleParameter("male-chrY-upper-threshold", strPrefix, params, doc, opts);

    pMethod->m_dFemaleChrXLowerThreshold = setupDoubleParameter("female-chrX-lower-threshold", strPrefix, params, doc, opts);
    pMethod->m_dFemaleChrXUpperThreshold = setupDoubleParameter("female-chrX-upper-threshold", strPrefix, params, doc, opts);
    pMethod->m_dFemaleChrYLowerThreshold = setupDoubleParameter("female-chrY-lower-threshold", strPrefix, params, doc, opts);
    pMethod->m_dFemaleChrYUpperThreshold = setupDoubleParameter("female-chrY-upper-threshold", strPrefix, params, doc, opts);

    pMethod->m_dMAPDThreshold = setupDoubleParameter("mapd-threshold", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNAnalysisMethodCNGender::CNAnalysisMethodCNGender()
{
    m_dMaleChrXLowerThreshold = 0;
    m_dMaleChrXUpperThreshold = 0;
    m_dMaleChrYLowerThreshold = 0;
    m_dMaleChrYUpperThreshold = 0;

    m_dFemaleChrXLowerThreshold = 0;
    m_dFemaleChrXUpperThreshold = 0;
    m_dFemaleChrYLowerThreshold = 0;
    m_dFemaleChrYUpperThreshold = 0;

    m_dMAPDThreshold = 0;
}

/**
 * @brief Run the analysis.
 */
void CNAnalysisMethodCNGender::run()
{
    Verbose::out(1, "CNAnalysisMethodCNGender::run(...) start");
    isSetup();

    // This method is revised to use new methods available in the base class
    // CNAnalysisMethod for finding positions and other meta data on
    // probes in the probe vector got from getProbeSets.

    double dChrXMean = 0.0;  // Accumulate then divide to get the mean.

    // .first gives the start and .second gives one past the end.
    pair<int,int>chrX = getChrBounds(m_iXChromosome, getProbeSets());
    for (int i=chrX.first; i<chrX.second; i++)
    {
        CNProbeSet * pobjProbeSet = getProbeSets()->at(i);
        if (!pobjProbeSet->isPseudoAutosomalRegion())
        {
            dChrXMean += pobjProbeSet->getCNState();  // Accumulate
        }
    }
    dChrXMean /= getProbeSetCount(m_iXChromosome, getProbeSets());  // Divide the sum
    getExperiment()->setChrXMean((float)dChrXMean); // Save it


    // Chromosome Y is done the same as chromosome X.
    double dChrYMean = 0.0;
    pair<int,int>chrY = getChrBounds(m_iYChromosome, getProbeSets());
    for (int i=chrY.first; i<chrY.second; i++)
    {
        CNProbeSet * pobjProbeSet = getProbeSets()->at(i);
        if (!pobjProbeSet->isPseudoAutosomalRegion())
        {
            dChrYMean += pobjProbeSet->getCNState();
        }
    }
    dChrYMean /= getProbeSetCount(m_iYChromosome, getProbeSets());
    getExperiment()->setChrYMean((float)dChrYMean);

    // The rest is the same as the old version.
    // Follow the rules for gender calling the same as the old version.
    if ((dChrXMean >= m_dMaleChrXLowerThreshold) && (dChrXMean <= m_dMaleChrXUpperThreshold) &&
        (dChrYMean >= m_dMaleChrYLowerThreshold) && (dChrYMean <= m_dMaleChrYUpperThreshold))
    {
        getExperiment()->setCNCallGender(affx::Male);
    }
    else if ((dChrXMean >= m_dFemaleChrXLowerThreshold) && (dChrXMean <= m_dFemaleChrXUpperThreshold) &&
        (dChrYMean >= m_dFemaleChrYLowerThreshold) && (dChrYMean <= m_dFemaleChrYUpperThreshold))
    {
        getExperiment()->setCNCallGender(affx::Female);
    }
    else
    {
        getExperiment()->setCNCallGender(affx::UnknownGender);
    }
    if (getExperiment()->getMadDiffCN() >= m_dMAPDThreshold)
    {
        getExperiment()->setCNCallGender(affx::UnknownGender);
        Verbose::out(1, "WARNING: MAPD over threshold for " + Fs::basename(getExperiment()->getExperimentName()) + ", using CN2Gender call of 'unknown'.");
    }
    Verbose::out(1, "CNAnalysisMethodCNGender::run(...) end");
}

/**
 * @brief Run the analysis. (The old version of run commented out)
 */
/*  Commented out on Dec 18, 2008 but left in for comparison.
    Feel free to delete if the above version is satisfactory.
void CNAnalysisMethodCNGender::run()
{
    Verbose::out(1, "CNAnalysisMethodCNGender::run(...) start");
    isSetup();

    double dChrXSum = 0;
    double dChrYSum = 0;
    int iChrXCount = 0;
    int iChrYCount = 0;
    for (int iRowIndex = 0; (iRowIndex < (int)getProbeSets()->size()); iRowIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iRowIndex);
        if (pobjProbeSet->getChromosome() == m_iXChromosome)
        {
            if ((!pobjProbeSet->isChrXPar1()) && (!pobjProbeSet->isChrXPar2()))
            {
                dChrXSum += pobjProbeSet->getCNState();
                iChrXCount++;
            }
        }
        else if (pobjProbeSet->getChromosome() == m_iYChromosome)
        {
            if ((!pobjProbeSet->isChrXPar1()) && (!pobjProbeSet->isChrXPar2()))
            {
                dChrYSum += pobjProbeSet->getCNState();
                iChrYCount++;
            }
        }
    }
    double dChrXMean = (dChrXSum / (double)iChrXCount);
    double dChrYMean = (dChrYSum / (double)iChrYCount);
    getExperiment()->setChrXMean((float)dChrXMean);
    getExperiment()->setChrYMean((float)dChrYMean);

    if ((dChrXMean >= m_dMaleChrXLowerThreshold) && (dChrXMean <= m_dMaleChrXUpperThreshold) &&
        (dChrYMean >= m_dMaleChrYLowerThreshold) && (dChrYMean <= m_dMaleChrYUpperThreshold))
    {
        getExperiment()->setGender(affx::Male);
    }
    else if ((dChrXMean >= m_dFemaleChrXLowerThreshold) && (dChrXMean <= m_dFemaleChrXUpperThreshold) &&
        (dChrYMean >= m_dFemaleChrYLowerThreshold) && (dChrYMean <= m_dFemaleChrYUpperThreshold))
    {
        getExperiment()->setGender(affx::Female);
    }
    else
    {
        getExperiment()->setGender(affx::UnknownGender);
    }
    if (getExperiment()->getMadDiffCN() >= m_dMAPDThreshold)
    {
        getExperiment()->setGender(affx::UnknownGender);
        Verbose::out(1, "WARNING: MAPD over threshold for " + Fs::basename(getExperiment()->getExperimentName()) + ", using CN2Gender call of 'unknown'.");
    }
    Verbose::out(1, "CNAnalysisMethodCNGender::run(...) end");
}
*/
