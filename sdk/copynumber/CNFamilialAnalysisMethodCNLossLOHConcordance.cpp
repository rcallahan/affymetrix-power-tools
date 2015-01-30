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
 * @file CNFamilialAnalysisMethodCNLossLOHConcordance.cpp
 *
 * @brief This file contains the CNFamilialAnalysisMethodCNLossLOHConcordance class members.
 */
#include "copynumber/CNFamilialAnalysisMethodCNLossLOHConcordance.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/AffxStatistics.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNFamilialAnalysisMethodCNLossLOHConcordance::explainSelf()
{
    CNFamilialAnalysisMethodCNLossLOHConcordance obj;
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
std::vector<SelfDoc::Opt> CNFamilialAnalysisMethodCNLossLOHConcordance::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  SelfDoc::Opt opt1 = {getType() + "-marker-count-cutoff", SelfDoc::Opt::Integer, "21", "21", "1", "NA", "CNLossLOHConcordance Marker Count Cutoff"};
  opts.push_back(opt1);

  SelfDoc::Opt opt2 = {getType() + "-call-cutoff", SelfDoc::Opt::Integer, "4", "4", "0", "NA",  "CNLossLOHConcordance Call Cutoff"};
  opts.push_back(opt2);

  SelfDoc::Opt opt3 = {getType() + "-min-genomic-span", SelfDoc::Opt::Integer, "1000000", "1000000", "1", "NA", "CNLossLOHConcordance Minimum Genomic Span"};
  opts.push_back(opt3);

  SelfDoc::Opt opt4 = {getType() + "-xChromosome", SelfDoc::Opt::Integer, "24", "24", "1", "NA", "CNLossLOHConcordance X Chromosome"};
  opts.push_back(opt4);

  SelfDoc::Opt opt5 = {getType() + "-yChromosome", SelfDoc::Opt::Integer, "25", "25", "1", "NA", "CNLossLOHConcordance Y Chromosome"};
  opts.push_back(opt5);

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
SelfCreate* CNFamilialAnalysisMethodCNLossLOHConcordance::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNFamilialAnalysisMethodCNLossLOHConcordance* pMethod = new CNFamilialAnalysisMethodCNLossLOHConcordance();
    std::string strPrefix = getPrefix();

    pMethod->m_iMarkerCountCutoff = setupIntParameter(getType() + "-marker-count-cutoff", strPrefix, params, doc, opts);
    pMethod->m_iCallCutoff = setupIntParameter(getType() + "-call-cutoff", strPrefix, params, doc, opts);
    pMethod->m_iMinGenomicSpan = setupIntParameter(getType() + "-min-genomic-span", strPrefix, params, doc, opts);
    pMethod->m_iXChromosome = setupIntParameter(getType() + "-xChromosome", strPrefix, params, doc, opts);
    pMethod->m_iYChromosome = setupIntParameter(getType() + "-yChromosome", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNFamilialAnalysisMethodCNLossLOHConcordance::CNFamilialAnalysisMethodCNLossLOHConcordance()
{
    m_iMarkerCountCutoff = 0;
    m_iCallCutoff = 0;
    m_iMinGenomicSpan = 0;
    m_iXChromosome = 0;
    m_iYChromosome = 0;
}

/**
 * @brief Run the analysis.
 */
void CNFamilialAnalysisMethodCNLossLOHConcordance::run()
{
    Verbose::out(1, "CNFamilialAnalysisMethodCNLossLOHConcordance::run(...) start");
    CNCychp& cychpIndex = getCychpIndex();
    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();

    getIndexSegments().deleteAll();
    getMotherSegments().deleteAll();
    getFatherSegments().deleteAll();

    int iProbeSetCount = cychpIndex.getCychpProbeSetsCopyNumbers().getCount();
    if (cychpMother.getFamilialCall() == 1)
    {
        Verbose::progressBegin(1, "CNFamilialAnalysisMethodCNLossLOHConcordance::run(...) Mother's segments ", iProbeSetCount, 100000, (iProbeSetCount / 100000));
        for (int i = 0; (i < iProbeSetCount); i++)
        {
            Verbose::progressStep(1);
            CNCychpProbeSetsCopyNumber* pIndex = cychpIndex.getCychpProbeSetsCopyNumbers().getAt(i);
            CNCychpProbeSetsCopyNumber* pMother = cychpMother.getCychpProbeSetsCopyNumbers().getAt(i);
            pMother->setSegmentInput(1);
            if ((pMother->getGenotypeCall() == -1) || (pIndex->getGenotypeCall() == -1))
            {
                pMother->setSegmentInput(-1);
            }
            else
            {
                pMother->setSegmentInput(1);
                if (pIndex->getLohCall() == 1)
                {
                    if (pIndex->getCNLoss(cychpIndex.getCychpHeader().getGender(), m_iXChromosome, m_iYChromosome) == 1)
                    {
                        if (isGenotypeConcordance(pIndex, pMother)) {pMother->setSegmentInput(0);}
                    }
                }
            }
        }
        anneal(cychpMother);
        newSegments(cychpIndex, cychpMother, 1, getMotherSegments(), getSegmentType(), 1, m_iMarkerCountCutoff, m_iCallCutoff, m_iMinGenomicSpan);
        Verbose::progressEnd(1, "Done");
    }
    if (cychpFather.getFamilialCall() == 1)
    {
        Verbose::progressBegin(1, "CNFamilialAnalysisMethodCNLossLOHConcordance::run(...) Father's segments ", iProbeSetCount, 100000, (iProbeSetCount / 100000));
        for (int i = 0; (i < iProbeSetCount); i++)
        {
            Verbose::progressStep(1);
            CNCychpProbeSetsCopyNumber* pIndex = cychpIndex.getCychpProbeSetsCopyNumbers().getAt(i);
            CNCychpProbeSetsCopyNumber* pFather = cychpFather.getCychpProbeSetsCopyNumbers().getAt(i);
            if ((pFather->getGenotypeCall() == -1) || (pIndex->getGenotypeCall() == -1))
            {
                pFather->setSegmentInput(-1);
            }
            else
            {
                pFather->setSegmentInput(1);
                if (pIndex->getLohCall() == 1)
                {
                    if (pIndex->getCNLoss(cychpIndex.getCychpHeader().getGender(), m_iXChromosome, m_iYChromosome) == 1)
                    {
                        if (isGenotypeConcordance(pIndex, pFather)) {pFather->setSegmentInput(0);}
                    }
                }
            }
        }
        anneal(cychpFather);
        newSegments(cychpIndex, cychpFather, 2, getFatherSegments(), getSegmentType(), 1, m_iMarkerCountCutoff, m_iCallCutoff, m_iMinGenomicSpan);
        Verbose::progressEnd(1, "Done");
    }

    Verbose::out(1, "CNFamilialAnalysisMethodCNLossLOHConcordance::run(...) end");
}
