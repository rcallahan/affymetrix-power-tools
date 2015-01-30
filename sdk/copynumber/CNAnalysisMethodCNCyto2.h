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

#ifndef _CNAnalysisMethodCNCyto2_H_
#define _CNAnalysisMethodCNCyto2_H_
/**
 * @file CNAnalysisMethodCNCyto2.h
 *
 * @brief This header contains the CNAnalysisMethodCNCyto2 class definition.
 */

#include "copynumber/CNAnalysisMethod.h"
//
#include <valarray>
//

/**
 * @brief  The CN analysis method.
 */

class CNAnalysisMethodCNCyto2 : public CNAnalysisMethod
{
private:
    // Copy Number State Parameters
    std::vector<int> m_vCNState;
    std::valarray<double> m_Mu;
    std::valarray<double> m_Sigma;
    double m_diagWeight;
    double m_mapdWeight;
    int m_minSegmentSize;
    double m_hmm_confidence_weight;

    std::vector<int> m_vCNState_X;
    std::valarray<double> m_Mu_X;
    std::valarray<double> m_Sigma_X;
    double m_diagWeight_X;
    double m_mapdWeight_X;
    int m_minSegmentSize_X;
    double m_hmm_confidence_weight_X;

    std::vector<int> m_vCNState_Y;
    std::valarray<double> m_Mu_Y;
    std::valarray<double> m_Sigma_Y;
    double m_diagWeight_Y;
    double m_mapdWeight_Y;
    int m_minSegmentSize_Y;
    double m_hmm_confidence_weight_Y;

    bool m_shrink;
    std::vector<double> m_shrink_lprec;
    double m_shrink_converge;
    bool m_shrink_downweight_outlier;
    double m_shrink_downweight_df;
    int m_shrink_downweight_maxiter;

    CNProbeSetArray  m_vLocalProbeSets;
    bool m_bLocalProbeSetsDetermined;

public:
    CNAnalysisMethodCNCyto2();
    virtual ~CNAnalysisMethodCNCyto2() {}

    static std::string getType() {return "cn-cyto2";}
    static std::string getDescription() {return "CopyNumber CNCyto2";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate * newObject(std::map<std::string, std::string>& param);
    void setCNCallsToInvalidState();
    virtual AffxString getName() { return getType(); }

    virtual void setEngine(BaseEngine* p)
    {
        AffxString str;
        m_pEngine = p;
        for (unsigned int i = 0; (i < m_vCNState.size()); i++)
        {
            if (i > 0) {str += ",";}
            str += ::getInt(m_vCNState[i]);
        }
        p->setOpt("hmmCN_state", str);
        str = "";
        for (unsigned int i = 0; (i < m_Mu.size()); i++)
        {
            if (i > 0) {str += ",";}
            str += ::getDouble(m_Mu[i], 6);
        }
        p->setOpt("hmmCN_mu",  str);

        // Set yTarget if m_Mu_Y has the second element
        if (m_Mu_Y.size() > 1) {m_pEngine->setOpt("yTarget", ToStr(m_Mu_Y[1]));}
    }

    virtual void run();
    virtual bool isSegmentTypeAnalysis() { return false; }

    void determineLocalProbeSets();

    CNProbeSetArray* getProbeSets()
    {
            return &m_vLocalProbeSets;
    }

protected:
    int getLastAutosomeChromosome();

    void copyNumberStatePostProcessing(
        int iChromosome,
        std::valarray<int> & vCNStates,
        std::valarray<float> & vCNConfidences,
        std::valarray<double> & shrunk_log2_ratios);

    void copyNumberStatePostProcessing(
        int iChromosome,
        std::valarray<int> & vCNStates,
        std::valarray<float> & vCNConfidences);

    double computeMAPD();

    int hmmPriorsVerify(const valarray<double>& mu_arr, int nRows, int nCols);

    void compute_confidence_mu_fixed(
        const int iExpectedCNState,
        const std::valarray<double> & samplr,
        std::valarray<int> & state,
        std::valarray<float> & confidence,
        const std::valarray<double> mu,
        const std::valarray<double> sigma,
        const double hmm_weight
    );


    void compute_confidence_mu_variable(
        const int iExpectedCNState,
        const std::valarray<double> & samplr,
        std::valarray<int> & state,
        std::valarray<float> & confidence,
        const std::valarray<double> mu,
        const std::valarray<double> sigma,
        const double hmm_weight
    );


    float compute_confidence(
        const int iExpectedCNState,
        const double logratio,
        const int cnstate,
        const std::valarray<double> mu,
        const std::valarray<double> sigma,
        const double hmm_weight
    );

};

#endif
