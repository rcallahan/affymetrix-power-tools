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

#ifndef _CNFamilialAnalysisMethod_H_
#define _CNFamilialAnalysisMethod_H_
/**
 * @file CNFamilialAnalysisMethod.h
 *
 * @brief This header contains the CNFamilialAnalysisMethod class definition.
 */

#include "copynumber/CNCychp.h"
#include "copynumber/CNSegment.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
#include "util/Err.h"
#include "util/Guid.h"
#include "util/Util.h"
#include "util/Verbose.h"
#include <vector>
#include <string>
#include <map>
#include <set>
//

/**
 * @brief  A base class for copy number quant methods.
 *
 */
class CNFamilialAnalysisMethod : public SelfDoc, public SelfCreate
{
protected:
    static std::vector<affymetrix_calvin_parameter::ParameterNameValueType> m_vParams;

protected:
    static AffxString getPrefix() {return "affymetrix-algorithm-param-";}
    static SelfDoc::Opt* getSelfDocOpt(std::vector<SelfDoc::Opt>& opts, const std::string& strName);
    static bool setupBoolParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string,std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
    static int setupIntParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string,std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
    static float setupFloatParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string,std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
    static std::string setupStringParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string,std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);

public:
    static std::vector<affymetrix_calvin_parameter::ParameterNameValueType>& getParams() {return m_vParams;}

protected:
    std::vector<unsigned int> m_vDataSetOffset;
    CNCychp* m_pcychpIndex;
    CNCychp* m_pcychpMother;
    CNCychp* m_pcychpFather;
    AffxString m_strAlleleFrequencyFileName;

    CNSegmentOverlapArray m_vSegmentOverlaps;

    CNSegmentArray m_vIndexSegments;
    CNSegmentArray m_vMotherSegments;
    CNSegmentArray m_vFatherSegments;

    std::map<std::string, std::set<std::string> > m_groupDatasets;

public:
    CNFamilialAnalysisMethod();
    virtual ~CNFamilialAnalysisMethod();

    virtual bool isFamilialAnalysis() = 0;
    virtual bool isSegmentOverlapAnalysis() = 0;
    virtual bool isSegmentTypeAnalysis() = 0;

    virtual AffxString getName() {return "";}
    int getSegmentType() {return CNSegment::getSegmentType(getName());}
    AffxString getSegmentTypeString() {return CNSegment::getSegmentTypeString(getSegmentType());}

    virtual void setup(CNCychp& cychpIndex, CNCychp& cychpMother, CNCychp& cychpFather, const AffxString& fileName)
    {
        m_pcychpIndex = &cychpIndex;
        m_pcychpMother = &cychpMother;
        m_pcychpFather = &cychpFather;
        m_strAlleleFrequencyFileName = fileName;
    }

    CNCychp& getCychpIndex() {return *m_pcychpIndex;}
    CNCychp& getCychpMother() {return *m_pcychpMother;}
    CNCychp& getCychpFather() {return *m_pcychpFather;}

    virtual void run() = 0;

    CNSegmentOverlapArray& getSegmentOverlaps() {return m_vSegmentOverlaps;}
    int getSegmentOverlapCount() {return m_vSegmentOverlaps.getCount();}

    CNSegmentArray& getIndexSegments() {return m_vIndexSegments;}
    CNSegmentArray& getMotherSegments() {return m_vMotherSegments;}
    CNSegmentArray& getFatherSegments() {return m_vFatherSegments;}

    int getSegmentCount() {return m_vIndexSegments.getCount() + m_vMotherSegments.getCount() + m_vFatherSegments.getCount();}

    void setDataSetOffsetVec(unsigned int ui) {m_vDataSetOffset.push_back(ui);}
    std::vector<unsigned int> getDataSetOffsetVec() {return m_vDataSetOffset;}

    void addGroupDataset(const std::string& group, const std::string& dataset);
    std::map<std::string, std::set<std::string> >& getGroupDatasets() {return m_groupDatasets;}

protected:
    void newSegments(CNCychp& cychpReference, CNCychp& cychpFamilial, int iFamilialSampleKey, CNSegmentArray& vSegments, int iSegmentType, int iDesiredCall, int iMarkerCountCutoff, int iHetCutoff, int iMinGenomicSpan);
    void anneal(CNCychp& cychp);

    bool isGenotypeConcordance(CNCychpProbeSetsCopyNumber* pProbeSetIndex, CNCychpProbeSetsCopyNumber* pProbeSetParent)
    {
        if ((pProbeSetParent->getGenotypeCall() == 0) && (pProbeSetIndex->getGenotypeCall() == 0)) {return true;}
        if ((pProbeSetParent->getGenotypeCall() == 1) && (pProbeSetIndex->getGenotypeCall() == 0)) {return true;}
        if ((pProbeSetParent->getGenotypeCall() == 1) && (pProbeSetIndex->getGenotypeCall() == 1)) {return true;}
        if ((pProbeSetParent->getGenotypeCall() == 1) && (pProbeSetIndex->getGenotypeCall() == 2)) {return true;}
        if ((pProbeSetParent->getGenotypeCall() == 2) && (pProbeSetIndex->getGenotypeCall() == 2)) {return true;}
        return false;
    }
};

#endif


