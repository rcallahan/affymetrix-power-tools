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
 * @file CNFamilialAnalysisMethod.cpp
 *
 * @brief This file contains the CNFamilialAnalysisMethod class members.
 */

#include "copynumber/CNFamilialAnalysisMethod.h"
//
#include "copynumber/CNAnalysisMethodLOH.h"
//
#include "portability/affy-base-types.h"
#include "util/AffxStatistics.h"
//
#include <limits>
//

std::vector<affymetrix_calvin_parameter::ParameterNameValueType> CNFamilialAnalysisMethod::m_vParams;

/**
 * @brief Constructor
 */
CNFamilialAnalysisMethod::CNFamilialAnalysisMethod()
{
    m_vDataSetOffset.clear();
    m_pcychpIndex = NULL;
    m_pcychpMother = NULL;
    m_pcychpFather = NULL;
}

/**
 * @brief Destructor
 */
CNFamilialAnalysisMethod::~CNFamilialAnalysisMethod()
{
    getSegmentOverlaps().deleteAll();

    getIndexSegments().deleteAll();
    getMotherSegments().deleteAll();
    getFatherSegments().deleteAll();
}

/**
 * Return SelfDoc option associated with a specified name.
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc options to search.
 * @param const std::string& - The specified name to search for.
 * @return SelfDoc::Opt* - The SelfDoc option pointer or NULL if not found.
 */
SelfDoc::Opt* CNFamilialAnalysisMethod::getSelfDocOpt(std::vector<SelfDoc::Opt>& opts, const std::string& strName)
{
    for (int iIndex = 0; (iIndex < opts.size()); iIndex++)
    {
        if (opts[iIndex].name == strName) {return &opts[iIndex];}
    }
    return NULL;
}

/**
 * Setup a bool parameter for this analysis method.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return bool - The value of the parameter as setup by this function.
 */
bool CNFamilialAnalysisMethod::setupBoolParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string,std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
    SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
    if (popt == NULL) {throw(Except("SelfDoc::Opt not found: " + strName));}
    bool b = AffxByteArray(popt->defaultVal).parsebool();
    fillInValue(b, std::string(strName), params, doc);
    std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
    affymetrix_calvin_parameter::ParameterNameValueType param;
    param.SetName(wstr);
    param.SetValueInt8(b);
    m_vParams.push_back(param);
    return b;
}

/**
 * Setup an int parameter for this analysis method, including min amd max validation.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return int - The value of the parameter as setup by this function.
 */
int CNFamilialAnalysisMethod::setupIntParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string,std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
    SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
    if (popt == NULL) {throw(Except("SelfDoc::Opt not found: " + strName));}
    int i = ::getInt(popt->defaultVal);
    fillInValue(i, std::string(strName), params, doc);
    if ((popt->minVal != "NA") && (popt->minVal != "") && (i < ::getInt(popt->minVal))) {throw(Except("SelfDoc::Opt " + strName + " below minimum value of " + popt->minVal));}
    if ((popt->maxVal != "NA") && (popt->maxVal != "") && (i > ::getInt(popt->maxVal))) {throw(Except("SelfDoc::Opt " + strName + " above maximum value of " + popt->maxVal));}
    std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
    affymetrix_calvin_parameter::ParameterNameValueType param;
    param.SetName(wstr);
    param.SetValueInt32(i);
    m_vParams.push_back(param);
    return i;
}

/**
 * Setup a float parameter for this analysis method, including min amd max validation.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return float - The value of the parameter as setup by this function.
 */
float CNFamilialAnalysisMethod::setupFloatParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string,std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
    SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
    if (popt == NULL) {throw(Except("SelfDoc::Opt not found: " + strName));}
    float f = (float)::getDouble(popt->defaultVal);
    fillInValue(f, std::string(strName), params, doc);
    if ((popt->minVal != "NA") && (popt->minVal != "") && (f < (float)::getDouble(popt->minVal))) {throw(Except("SelfDoc::Opt " + strName + " below minimum value of " + popt->minVal));}
    if ((popt->maxVal != "NA") && (popt->maxVal != "") && (f > (float)::getDouble(popt->maxVal))) {throw(Except("SelfDoc::Opt " + strName + " above maximum value of " + popt->maxVal));}
    std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
    affymetrix_calvin_parameter::ParameterNameValueType param;
    param.SetName(wstr);
    param.SetValueFloat(f);
    m_vParams.push_back(param);
    return f;
}

/**
 * Setup a string parameter for this analysis method, including min amd max validation.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return std::string - The value of the parameter as setup by this function.
 */
std::string CNFamilialAnalysisMethod::setupStringParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string,std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
    SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
    if (popt == NULL) {throw(Except("SelfDoc::Opt not found: " + strName));}
    std::string str = popt->defaultVal;
    fillInValue(str, std::string(strName), params, doc);
    if ((popt->minVal != "") && (str < popt->minVal)) {throw(Except("SelfDoc::Opt " + strName + " below minimum value of " + popt->minVal));}
    if ((popt->maxVal != "") && (str > popt->maxVal)) {throw(Except("SelfDoc::Opt " + strName + " above maximum value of " + popt->maxVal));}
    std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
    affymetrix_calvin_parameter::ParameterNameValueType param;
    param.SetName(wstr);
    param.SetValueAscii(str);
    m_vParams.push_back(param);
    return str;
}

/**
 * @brief Create new segments for the familial file
 * @param CNCychp& - Index's cychp data
 * @param CNCychp& - Parent's cychp data
 * @param int - The sample key for the parent
 * @param CNSegmentArray& - The segment array to load
 * @param int - The segment type to use
 * @param int - The desired call to load
 * @param int - The marker count cutoff
 * @param int - The call cutoff
 * @param int - The minimum genomic span for this type of segment
 */
void CNFamilialAnalysisMethod::newSegments(CNCychp& cychpReference, CNCychp& cychpFamilial, int iFamilialSampleKey, CNSegmentArray& vSegments, int iSegmentType, int iDesiredCall, int iMarkerCountCutoff, int iHetCutoff, int iMinGenomicSpan)
{
    vSegments.deleteAll();
    int iSnpCount = cychpReference.getCychpProbeSetsCopyNumbers().getCount();
    for (int iIndex = 0; (iIndex < iSnpCount); iIndex++)
    {
        CNCychpProbeSetsCopyNumber* pobjSnp = cychpFamilial.getCychpProbeSetsCopyNumbers().at(iIndex);
        pobjSnp->setSegmentOutput(-1);
    }
    CNAnalysisMethodLOH objSegmentAnalysis;
    int iLastChromosome = cychpReference.getCychpProbeSetsCopyNumbers().at(iSnpCount - 1)->getChromosome();
    for (int iChromosome = 1; (iChromosome <= iLastChromosome); iChromosome++)
    {
        int iChromosomeSnpCount = 0;
        for (int iIndex = 0; (iIndex < iSnpCount); iIndex++)
        {
            CNCychpProbeSetsCopyNumber* pobjSnp = cychpFamilial.getCychpProbeSetsCopyNumbers().at(iIndex);
            if ((pobjSnp->getSegmentInput() >= 0) && (pobjSnp->getChromosome() == iChromosome))
            {
                iChromosomeSnpCount++;
            }
        }
        std::vector<char> vChromosomeGenotypeCalls(iChromosomeSnpCount);
        std::vector<int> vChromosomePositions(iChromosomeSnpCount);
        std::vector<int> vChromosomeSnpIndexes(iChromosomeSnpCount);
        int iSnpIndex = 0;
        for (int iIndex = 0; (iIndex < iSnpCount); iIndex++)
        {
            CNCychpProbeSetsCopyNumber* pobjSnp = cychpFamilial.getCychpProbeSetsCopyNumbers().at(iIndex);
            if ((pobjSnp->getSegmentInput() >= 0) && (pobjSnp->getChromosome() == iChromosome))
            {
                vChromosomeGenotypeCalls[iSnpIndex] = pobjSnp->getSegmentInput();
                vChromosomePositions[iSnpIndex] = pobjSnp->getPosition();
                vChromosomeSnpIndexes[iSnpIndex] = iIndex;
                iSnpIndex++;
            }
        }
        int iStartIndex = 0;
        for (int iIndex = 1; (iIndex <= iChromosomeSnpCount); iIndex++)
        {
            if ((iIndex == iChromosomeSnpCount) || ((cychpFamilial.getCychpProbeSetsCopyNumbers().at(vChromosomeSnpIndexes[iIndex])->getPosition() - cychpFamilial.getCychpProbeSetsCopyNumbers().at(vChromosomeSnpIndexes[iIndex - 1])->getPosition()) > 1000000))
            {
                int iSegmentSnpCount = ((iIndex - 1) - iStartIndex + 1);
                if (iSegmentSnpCount >= iMarkerCountCutoff)
                {
                    std::vector<char> vSegmentGenotypeCalls(iSegmentSnpCount);
                    std::vector<int> vSegmentPositions(iSegmentSnpCount);
                    std::vector<char> vLoh(iSegmentSnpCount);
                    for (int iSegmentIndex = 0; (iSegmentIndex < iSegmentSnpCount); iSegmentIndex++)
                    {
                        vSegmentGenotypeCalls[iSegmentIndex] = vChromosomeGenotypeCalls[iStartIndex + iSegmentIndex];
                        vSegmentPositions[iSegmentIndex] = vChromosomePositions[iStartIndex + iSegmentIndex];
                    }
                    objSegmentAnalysis.lohFind(vSegmentGenotypeCalls, vSegmentPositions, iMarkerCountCutoff, iHetCutoff, iMinGenomicSpan, vLoh);
                    for (int iSegmentIndex = 0; (iSegmentIndex < iSegmentSnpCount); iSegmentIndex++)
                    {
                        cychpFamilial.getCychpProbeSetsCopyNumbers().at(vChromosomeSnpIndexes[iStartIndex + iSegmentIndex])->setSegmentOutput(vLoh[iSegmentIndex]);
                    }
                }
                iStartIndex = iIndex;
            }
        }
    }

    int iMarkerCount = 0;
    CNSegment* pobjSegment = NULL;
    CNCychpProbeSetsCopyNumber* pobjSnp = NULL;
    CNCychpProbeSetsCopyNumber* pIndexSnp = NULL;
    int iStartIndex = 0;
    for (int iIndex = 0; (iIndex < cychpFamilial.getCychpProbeSetsCopyNumbers().getCount()); iIndex++)
    {
        pobjSnp = cychpFamilial.getCychpProbeSetsCopyNumbers().getAt(iIndex);
        pIndexSnp = cychpReference.getCychpProbeSetsCopyNumbers().getAt(iIndex);
        if (pobjSnp->getSegmentOutput() != -1)
        {
            iStartIndex = iIndex + 1;
            break;
        }
        pobjSnp = NULL;
        pIndexSnp = NULL;
    }
    int iHetCount = 0;
    int iHomCount = 0;
    int iAgreeMarkerCount = 0;
    if (pobjSnp != NULL)
    {
        iHetCount = 0;
        iHomCount = 0;
        pobjSegment = new CNSegment;
        pobjSegment->setSegmentType(iSegmentType);
        pobjSegment->setFamilialSampleKey(iFamilialSampleKey);
        pobjSegment->setChromosome(pobjSnp->getChromosome());
        pobjSegment->setStartPosition(pobjSnp->getPosition());
        pobjSegment->setEndPosition(pobjSnp->getPosition());
        pobjSegment->setCall((char)pobjSnp->getSegmentOutput());
        if ((pobjSnp->getSegmentInput() == 0) && (pobjSnp->getSegmentOutput() == 1)) {iAgreeMarkerCount++;}
        if ((pobjSnp->getSegmentInput() == 1) && (pobjSnp->getSegmentOutput() == 0)) {iAgreeMarkerCount++;}
        if (pIndexSnp->getGenotypeCall() == 1) {iHetCount++;} else if (pIndexSnp->getGenotypeCall() != -1) {iHomCount++;}
        iMarkerCount = 1;
        for (int iIndex = iStartIndex; (iIndex < cychpFamilial.getCychpProbeSetsCopyNumbers().getCount()); iIndex++)
        {
            pobjSnp = cychpFamilial.getCychpProbeSetsCopyNumbers().getAt(iIndex);
            pIndexSnp = cychpReference.getCychpProbeSetsCopyNumbers().getAt(iIndex);
            if (pobjSnp->getSegmentOutput() != -1)
            {
                if ((pobjSnp->getChromosome() != pobjSegment->getChromosome()) || (pobjSnp->getSegmentOutput() != pobjSegment->getCall()))
                {
                    pobjSegment->setSegmentName(affxutil::Guid::GenerateNewGuid());
                    pobjSegment->setMarkerCount(iMarkerCount);
                    pobjSegment->setHeterozygosity((double)iHetCount / (double)iMarkerCount);
                    pobjSegment->setHomozygosity((double)iHomCount / (double)iMarkerCount);
                    pobjSegment->setConfidence((double)iAgreeMarkerCount / (double)iMarkerCount);
                    if ((iDesiredCall == -1) || (pobjSegment->getCall() == iDesiredCall))
                    {
                        pobjSegment->setCall(1);
                        vSegments.add(pobjSegment);
                    } else {delete pobjSegment;}
                    iHetCount = 0;
                    iHomCount = 0;
                    pobjSegment = new CNSegment;
                    pobjSegment->setSegmentType(iSegmentType);
                    pobjSegment->setFamilialSampleKey(iFamilialSampleKey);
                    pobjSegment->setChromosome(pobjSnp->getChromosome());
                    pobjSegment->setStartPosition(pobjSnp->getPosition());
                    pobjSegment->setCall((char)pobjSnp->getSegmentOutput());
                    iAgreeMarkerCount = 0;
                    iMarkerCount = 0;
                }
                if ((pobjSnp->getSegmentInput() == 0) && (pobjSnp->getSegmentOutput() == 1)) {iAgreeMarkerCount++;}
                if ((pobjSnp->getSegmentInput() == 1) && (pobjSnp->getSegmentOutput() == 0)) {iAgreeMarkerCount++;}
                if (pIndexSnp->getGenotypeCall() == 1) {iHetCount++;} else if (pIndexSnp->getGenotypeCall() != -1) {iHomCount++;}
                iMarkerCount++;
                pobjSegment->setEndPosition(pobjSnp->getPosition());
            }
        }
    }
    if (pobjSegment != NULL)
    {
        pobjSegment->setMarkerCount(iMarkerCount);
        pobjSegment->setSegmentName(affxutil::Guid::GenerateNewGuid());
        pobjSegment->setHeterozygosity((double)iHetCount / (double)iMarkerCount);
        pobjSegment->setHomozygosity((double)iHomCount / (double)iMarkerCount);
        pobjSegment->setConfidence((double)iAgreeMarkerCount / (double)iMarkerCount);
        if (pobjSegment->getCall() == iDesiredCall)
        {
            pobjSegment->setCall(1);
            vSegments.add(pobjSegment);
        } else {delete pobjSegment;}
    }
}

/**
 * @brief Anneal segment input field in the specified cychp data
 * @param CNCychp& - The cychp data to anneal
 */
void CNFamilialAnalysisMethod::anneal(CNCychp& cychp)
{
    for (int i = 1; (i < (cychp.getCychpProbeSetsCopyNumbers().getCount() - 1)); i++)
    {
        CNCychpProbeSetsCopyNumber* p0 = cychp.getCychpProbeSetsCopyNumbers().getAt(i - 1);
        CNCychpProbeSetsCopyNumber* p1 = cychp.getCychpProbeSetsCopyNumbers().getAt(i);
        CNCychpProbeSetsCopyNumber* p2 = cychp.getCychpProbeSetsCopyNumbers().getAt(i + 1);
        if (p1->getSegmentInput() == -1)
        {
            if (p0->getSegmentInput() == p2->getSegmentInput())
            {
                p1->setSegmentInput(p0->getSegmentInput());
            }
        }
    }
}


void CNFamilialAnalysisMethod::addGroupDataset(const string& group, const string& dataset)
{
    m_groupDatasets[group].insert(dataset);
}
