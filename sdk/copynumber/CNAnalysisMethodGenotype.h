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

#ifndef _CNAnalysisMethodGenotype_H_
#define _CNAnalysisMethodGenotyhpe_H_
/**
 * @file CNAnalysisMethodGenotype.h
 *
 * @brief This header contains the CNAnalysisMethodGenotype class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

#include "copynumber/CNAnalysisMethodChipstream.h"
//
#include "copynumber/CNAnalysisMethodCovariateSignalAdjuster.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "copynumber/CNAnalysisMethodReference.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/BioTypes.h"
#include "chipstream/CelReader.h"
#include "chipstream/ChipLayout.h"
#include "chipstream/ChipStreamFactory.h"
#include "chipstream/HomHiLoCelListener.h"
#include "chipstream/QuantLabelZIO.h"
#include "chipstream/SketchQuantNormTran.h"
#include "chipstream/SparseMart.h"
#include "file/TsvFile/TsvFile.h"
#include "normalization/normalization.h"
#include "util/AffxStatistics.h"
#include "util/RowFile.h"




/**
 * @brief  The Genotype analysis method.
 *
 */
class CNAnalysisMethodGenotype : public CNAnalysisMethod
{
public:

    CNProbeSetArray  m_vLocalProbeSets;
    bool m_bLocalProbeSetsDetermined;

public:
    CNAnalysisMethodGenotype();
    virtual ~CNAnalysisMethodGenotype() {}

    static std::string getType() {return "genotype";}
    static std::string getDescription() {return "Genotype";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    virtual AffxString getName() {return getType();}
    virtual void run();
    virtual bool isSegmentTypeAnalysis() {return false;}

private:
    int m_iPrevGenderCode;
    int m_iGenderCode;
    int m_iSeparation;


public:

private:
    void calculateGenotypes();
    float calculateCallRate();
    void initializeSnpParameters(snp_param& sp, int iCallMethod);
    void setBrlmmpParameters(snp_param& sp);
    void loadSnpDistributions();
    void imputeUncalledCNstateForSNPs(CNProbeSetArray* pvProbeSets);
    void impute(int iBeginIndex, int iEndIndex, CNProbeSetArray* pvProbeSets );

};

#endif


