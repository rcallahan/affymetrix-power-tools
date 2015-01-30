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
 * @file CNFamilialAnalysisMethodSegmentOverlap.cpp
 *
 * @brief This file contains the CNFamilialAnalysisMethodSegmentOverlap class members.
 */
#include "copynumber/CNFamilialAnalysisMethodSegmentOverlap.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "util/AffxStatistics.h"
//

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNFamilialAnalysisMethodSegmentOverlap::explainSelf()
{
    CNFamilialAnalysisMethodSegmentOverlap obj;
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
std::vector<SelfDoc::Opt> CNFamilialAnalysisMethodSegmentOverlap::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  SelfDoc::Opt opt1 = {getType() + "-threshold", SelfDoc::Opt::Float, "0.75", "0.75", "0", "1", "SegmentOverlap threshold"};
  opts.push_back(opt1);

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
SelfCreate* CNFamilialAnalysisMethodSegmentOverlap::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNFamilialAnalysisMethodSegmentOverlap* pMethod = new CNFamilialAnalysisMethodSegmentOverlap();
    std::string strPrefix = getPrefix();

    pMethod->m_fSegmentOverlapThreshold = setupFloatParameter(getType() + "-threshold", strPrefix, params, doc, opts);

    return pMethod;
}

/**
 * @brief Constructor
 */
CNFamilialAnalysisMethodSegmentOverlap::CNFamilialAnalysisMethodSegmentOverlap()
{
    m_fSegmentOverlapThreshold = 0;
}

/**
 * @brief Run the analysis.
 */
void CNFamilialAnalysisMethodSegmentOverlap::run()
{
    Verbose::out(1, "CNFamilialAnalysisMethodSegmentOverlap::run(...) start");
    CNCychp& cychpIndex = getCychpIndex();
    CNCychp& cychpMother = getCychpMother();
    CNCychp& cychpFather = getCychpFather();

    getIndexSegments().deleteAll();
    getMotherSegments().deleteAll();
    getFatherSegments().deleteAll();

    CNCychpSegments& vIndexSegments = cychpIndex.getCychpSegments();
    CNCychpSegments& vMotherSegments = cychpMother.getCychpSegments();
    CNCychpSegments& vFatherSegments = cychpFather.getCychpSegments();

    if (cychpMother.getFamilialCall() == 1)
    {
        Verbose::progressBegin(1, "CNFamilialAnalysisMethodSegmentOverlap::run(...) Mother's segments ", vIndexSegments.getCount(), 100, (vIndexSegments.getCount() / 100));
        for (int i0 = 0; (i0 < vIndexSegments.getCount()); i0++)
        {
            Verbose::progressStep(1);
            CNCychpSegment* p0 = vIndexSegments.getAt(i0);
            for (int i1 = 0; (i1 < vMotherSegments.getCount()); i1++)
            {
                CNCychpSegment* p1 = vMotherSegments.getAt(i1);
                if (p0->getDataSetIndex() == p1->getDataSetIndex())
                {
                    if (p0->getChromosome() == p1->getChromosome())
                    {
                        double dSegmentOverlap = 0;
                        if ((p0->getStartPosition() <= p1->getStartPosition()) && (p0->getStopPosition() >= p1->getStartPosition()))
                        {
                            double d1 = Min(p0->getStopPosition(), p1->getStopPosition()) - p1->getStartPosition() + 1;
                            double d2 = p0->getStopPosition() - p0->getStartPosition() + 1;
                            dSegmentOverlap = (d1 / d2);
                        }
                        else if ((p1->getStartPosition() <= p0->getStartPosition()) && (p1->getStopPosition() >= p0->getStartPosition()))
                        {
                            double d1 = Min(p0->getStopPosition(), p1->getStopPosition()) - p0->getStartPosition() + 1;
                            double d2 = p0->getStopPosition() - p0->getStartPosition() + 1;
                            dSegmentOverlap = (d1 / d2);
                        }
                        if ((dSegmentOverlap >= m_fSegmentOverlapThreshold) && (p0->getSegmentType() != "Unknown"))
                        {
                            CNSegmentOverlap* p = new CNSegmentOverlap;
                            p->ReferenceSampleKey = 0;
                            p->FamilialSampleKey = 1;
                            p->ReferenceSegmentID = p0->getSegmentID();
                            p->FamilialSegmentID = p1->getSegmentID();
                            p->SegmentType = p0->getSegmentType();
                            m_vSegmentOverlaps.add(p);
                        }
                    }
                }
            }
        }
        Verbose::progressEnd(1, "Done");
    }
    if (cychpFather.getFamilialCall() == 1)
    {
        Verbose::progressBegin(1, "CNFamilialAnalysisMethodSegmentOverlap::run(...) Father's segments ", vIndexSegments.getCount(), 100, (vIndexSegments.getCount() / 100));
        for (int i0 = 0; (i0 < vIndexSegments.getCount()); i0++)
        {
            Verbose::progressStep(1);
            CNCychpSegment* p0 = vIndexSegments.getAt(i0);
            for (int i1 = 0; (i1 < vFatherSegments.getCount()); i1++)
            {
                CNCychpSegment* p1 = vFatherSegments.getAt(i1);
                if (p0->getDataSetIndex() == p1->getDataSetIndex())
                {
                    if (p0->getChromosome() == p1->getChromosome())
                    {
                        double dSegmentOverlap = 0;
                        if ((p0->getStartPosition() <= p1->getStartPosition()) && (p0->getStopPosition() >= p1->getStartPosition()))
                        {
                            double d1 = Min(p0->getStopPosition(), p1->getStopPosition()) - p1->getStartPosition() + 1;
                            double d2 = p0->getStopPosition() - p0->getStartPosition() + 1;
                            dSegmentOverlap = (d1 / d2);
                        }
                        else if ((p1->getStartPosition() <= p0->getStartPosition()) && (p1->getStopPosition() >= p0->getStartPosition()))
                        {
                            double d1 = Min(p0->getStopPosition(), p1->getStopPosition()) - p0->getStartPosition() + 1;
                            double d2 = p0->getStopPosition() - p0->getStartPosition() + 1;
                            dSegmentOverlap = (d1 / d2);
                        }
                        if ((dSegmentOverlap >= m_fSegmentOverlapThreshold) && (p0->getSegmentType() != "Unknown"))
                        {
                            CNSegmentOverlap* p = new CNSegmentOverlap;
                            p->ReferenceSampleKey = 0;
                            p->FamilialSampleKey = 2;
                            p->ReferenceSegmentID = p0->getSegmentID();
                            p->FamilialSegmentID = p1->getSegmentID();
                            p->SegmentType = p0->getSegmentType();
                            m_vSegmentOverlaps.add(p);
                        }
                    }
                }
            }
        }
        Verbose::progressEnd(1, "Done");
    }
    Verbose::out(1, "CNFamilialAnalysisMethodSegmentOverlap::run(...) end");
}
