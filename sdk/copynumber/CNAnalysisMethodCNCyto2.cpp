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
 * @file CNAnalysisMethodCNCyto2.cpp
 *
 * @brief This file contains the CNAnalysisMethodCNCyto2 class members.
 */

// Microsoft seems to defeat some of the purpose of using <cmath> leaving
// certain things undefined such as "M_PI".
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include "copynumber/CNAnalysisMethodCNCyto2.h"

#include "copynumber/CNSmoother.h"
#include "copynumber/TransitionMatrix.h"
#include "copynumber/CopyNumberHMM.h"
#include "copynumber/CNAnalysisMethodFactory.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"

#include "copynumber/WaveletShrink.h"

//
#include <valarray>
//

//  All the includes should be done with before this line
using namespace std;


CNAnalysisMethodCNCyto2::CNAnalysisMethodCNCyto2()
{
        m_bLocalProbeSetsDetermined=false;
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodCNCyto2::explainSelf()
    {
    SelfDoc doc;
    doc.setDocName(getType());
    doc.setDocDescription(getDescription());
    doc.setDocOptions(getDefaultDocOptions());
    return doc;
    }

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
vector<SelfDoc::Opt> CNAnalysisMethodCNCyto2::getDefaultDocOptions()
    {
    vector<SelfDoc::Opt> opts;

    SelfDoc::Opt hmmCN_state = {"hmmCN_state", SelfDoc::Opt::String,
        "0,1,2,3,4,5", "0,1,2,3,4,5", "", "", "CN HMM State"};
    opts.push_back(hmmCN_state);

    SelfDoc::Opt hmmCN_mu = {"hmmCN_mu", SelfDoc::Opt::String,
        "-1.63,-0.58,0,0.45,0.72,0.93", "-1.63,-0.58,0,0.45,0.72,0.93", "", "",
        "CN HMM Mu"};
    opts.push_back(hmmCN_mu);

    SelfDoc::Opt hmmCN_sigma = {"hmmCN_sigma", SelfDoc::Opt::String,
        "0.3,0.3,0.3,0.3,0.3,0.3", "0.3,0.3,0.3,0.3,0.3,0.3", "", "",
        "CN HMM Sigma"};
    opts.push_back(hmmCN_sigma);

    SelfDoc::Opt diagonal_weight = {"diagonal-weight", SelfDoc::Opt::Double,
        "0.995", "0.995", "0.0", "1.0",
        "CN HMM Diagonal Weight"};
    opts.push_back(diagonal_weight);

    SelfDoc::Opt mapd_weight = {"mapd-weight", SelfDoc::Opt::Double,
        "0.22", "0.22", "0.0", "10000.0",
        "CN MAPD Weight"};
    opts.push_back(mapd_weight);

    SelfDoc::Opt min_segment_size = {"min-segment-size", SelfDoc::Opt::Integer,
        "5", "5", "1", "NA",
        "CN Min Segment Size"};
    opts.push_back(min_segment_size);

    SelfDoc::Opt hmm_confidence_weight = {"hmm-confidence-weight",
        SelfDoc::Opt::Double, "0.6", "0.6", "0.0", "1.0",
        "HMM influence on confidence"};
    opts.push_back(hmm_confidence_weight);

    // For X chromosome
    SelfDoc::Opt hmmCN_state_X = {"hmmCN_state-X", SelfDoc::Opt::String,
        "0,1,2,3,4,5", "0,1,2,3,4,5", "", "", "CN HMM State for X"};
    opts.push_back(hmmCN_state_X);

    SelfDoc::Opt hmmCN_mu_X = {"hmmCN_mu-X", SelfDoc::Opt::String,
        "-1.63,-0.58,0,0.45,0.72,0.93", "-1.63,-0.58,0,0.45,0.72,0.93", "", "",
        "CN HMM Mu for X"};
    opts.push_back(hmmCN_mu_X);

    SelfDoc::Opt hmmCN_sigma_X = {"hmmCN_sigma-X", SelfDoc::Opt::String,
        "0.3,0.3,0.3,0.3,0.3,0.3", "0.3,0.3,0.3,0.3,0.3,0.3", "", "",
        "CN HMM Sigma for X"};
    opts.push_back(hmmCN_sigma_X);

    SelfDoc::Opt diagonal_weight_X = {"diagonal-weight-X", SelfDoc::Opt::Double,
        "0.995", "0.995", "0.0", "1.0",
        "CN HMM Diagonal Weight for X"};
    opts.push_back(diagonal_weight_X);

    SelfDoc::Opt mapd_weight_X = {"mapd-weight-X", SelfDoc::Opt::Double,
        "0.22", "0.22", "0.0", "10000.0",
        "CN MAPD Weight for X"};
    opts.push_back(mapd_weight_X);

    SelfDoc::Opt min_segment_size_X = {"min-segment-size-X", SelfDoc::Opt::Integer,
        "5", "5", "1", "NA",
        "CN Min Segment Size for X"};
    opts.push_back(min_segment_size_X);

    SelfDoc::Opt hmm_confidence_weight_X = {"hmm-confidence-weight-X",
        SelfDoc::Opt::Double, "0.6", "0.6", "0.0", "1.0",
        "HMM influence on confidence for X"};
    opts.push_back(hmm_confidence_weight_X);

    // For Y chromosome
    SelfDoc::Opt hmmCN_state_Y = {"hmmCN_state-Y", SelfDoc::Opt::String,
        "0,1,2,3,4,5", "0,1,2,3,4,5", "", "", "CN HMM State for Y"};
    opts.push_back(hmmCN_state_Y);

    SelfDoc::Opt hmmCN_mu_Y = {"hmmCN_mu-Y", SelfDoc::Opt::String,
        "-1.63,-0.58,0,0.45,0.72,0.93", "-1.63,-0.58,0,0.45,0.72,0.93", "", "",
        "CN HMM Mu for Y"};
    opts.push_back(hmmCN_mu_Y);

    SelfDoc::Opt hmmCN_sigma_Y = {"hmmCN_sigma-Y", SelfDoc::Opt::String,
        "0.3,0.3,0.3,0.3,0.3,0.3", "0.3,0.3,0.3,0.3,0.3,0.3", "", "",
        "CN HMM Sigma for Y"};
    opts.push_back(hmmCN_sigma_Y);

    SelfDoc::Opt diagonal_weight_Y = {"diagonal-weight-Y", SelfDoc::Opt::Double,
        "0.995", "0.995", "0.0", "1.0",
        "CN HMM Diagonal Weight for Y"};
    opts.push_back(diagonal_weight_Y);

    SelfDoc::Opt mapd_weight_Y = {"mapd-weight-Y", SelfDoc::Opt::Double,
        "0.22", "0.22", "0.0", "10000.0",
        "CN MAPD Weight for Y"};
    opts.push_back(mapd_weight_Y);

    SelfDoc::Opt min_segment_size_Y = {"min-segment-size-Y", SelfDoc::Opt::Integer,
        "5", "5", "1", "NA",
        "CN Min Segment Size for Y"};
    opts.push_back(min_segment_size_Y);

    SelfDoc::Opt hmm_confidence_weight_Y = {"hmm-confidence-weight-Y",
        SelfDoc::Opt::Double, "0.6", "0.6", "0.0", "1.0",
        "HMM influence on confidence for Y"};
    opts.push_back(hmm_confidence_weight_Y);

//Shrinkage
    SelfDoc::Opt shrink = {"shrink",
        SelfDoc::Opt::Boolean, "false", "false", "NA", "NA",
        "Apply wavelet shrinkage to the log2 ratios."};
    opts.push_back(shrink);

    SelfDoc::Opt shrink_lprec = {"shrink-lprec", SelfDoc::Opt::String,
        "0.5,4.0", "0.5,4.0", "", "",
        "Precision of the likelihood for wavelet shrinkage"};
    opts.push_back(shrink_lprec);

    SelfDoc::Opt shrink_converge = {"shrink-converge", SelfDoc::Opt::Double,
        "0.000001", "0.000001", "0.00000001", "10000.0",
        "Convergence criteria for shrinkage"};
    opts.push_back(shrink_converge);

    SelfDoc::Opt shrink_downweight_outlier = {"shrink-downweight-outlier",
        SelfDoc::Opt::Boolean, "true", "true", "NA", "NA",
        "If true, then outlier values are downweighted when shrunk."};
    opts.push_back(shrink_downweight_outlier);

    SelfDoc::Opt shrink_downweight_df = {"shrink-downweight-df",
        SelfDoc::Opt::Double, "6.5", "6.5", "1.0", "100.0",
        "Degrees of freedom for downweighting outliers during shrinkage"};
    opts.push_back(shrink_downweight_df);

    SelfDoc::Opt shrink_downweight_maxiter = {"shrink-downweight-maxiter",
        SelfDoc::Opt::Integer, "15", "15", "1", "100",
        "Maximum iterations for downweighting outliers during shrinkage"};
    opts.push_back(shrink_downweight_maxiter);

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
SelfCreate * CNAnalysisMethodCNCyto2::
    newObject(map<string,string>& params)
    {
    SelfDoc doc = explainSelf();
    vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodCNCyto2 * pMethod = new CNAnalysisMethodCNCyto2();
    string strPrefix = getPrefix();

    AffxString str;
    AffxByteArray ba;

    // Diagonal
    pMethod->m_diagWeight = setupDoubleParameter("diagonal-weight", strPrefix,
        params, doc, opts);

    // MAPD weight used in HMM dispersion estimate
    pMethod->m_mapdWeight = setupDoubleParameter("mapd-weight", strPrefix,
        params, doc, opts);


    // Minimum size for a segment
    pMethod->m_minSegmentSize = setupIntParameter("min-segment-size",
        strPrefix, params, doc, opts);

    // m_vCNState
    str = setupStringParameter("hmmCN_state", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_vCNState.resize(ba.parameterCount());
    for (int k=0; k<ba.parameterCount(); k++)
        {
        pMethod->m_vCNState[k] = ba.getParameter(k+1).parseInt();
        }

    // hmmCN_mu
    str = setupStringParameter("hmmCN_mu", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_Mu.resize(ba.parameterCount());
    for (int k=0; k<ba.parameterCount(); k++)
        {
        pMethod->m_Mu[k] = ba.getParameter(k+1).parseDouble();
        }

    // hmmCN_sigma
    str = setupStringParameter("hmmCN_sigma", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_Sigma.resize(ba.parameterCount());
    for (int k=0; k<ba.parameterCount(); k++)
        {
        pMethod->m_Sigma[k] = ba.getParameter(k+1).parseDouble();
        }

    if (pMethod->m_vCNState.size() != pMethod->m_Mu.size())
        {
        throw(Except("Mu must conform to CNState."));
        }

    if (pMethod->m_vCNState.size() != pMethod->m_Sigma.size())
        {
        throw(Except("Sigma must conform to CNState."));
        }


    for (int k=0; k<pMethod->m_vCNState.size(); k++)
        {
        if (pMethod->m_Sigma[k] <= 0.0)
            throw(Except("Sigma parameters must be > 0."));

        if ((k > 0) && (pMethod->m_Mu[k-1] > pMethod->m_Mu[k]))
            throw(Except("Mu parameters must be ordered increasing."));
        }

    // MAPD weight used in HMM dispersion estimate
    pMethod->m_hmm_confidence_weight = setupDoubleParameter(
            "hmm-confidence-weight", strPrefix, params, doc, opts);

    if (pMethod->m_hmm_confidence_weight <= 0.0 ||
            pMethod->m_hmm_confidence_weight >= 1.0)
        {
        throw(Except("HMM Confidence weight must be in [0,1]"));
        }

    // Same parameters for X chromosome (if separate set of parameters is required for X)

    // NB: The X- and Y-specific values of diagonal-weight, mapd-weight, min-segment-size,
    // hmm-confidence-weight, m_vCNState, m_Mu, m_Sigma, must be set to the corresponding
    // autosomal values (NOT the defaults) if they haven't been specified on the cn-cyto2 analysis line.
    std::map<string, string>::iterator iter;

    // Diagonal
    pMethod->m_diagWeight_X = setupDoubleParameter("diagonal-weight-X", strPrefix,
        params, doc, opts);
    iter = params.find("diagonal-weight-X");
    if (iter == params.end())
    {
        pMethod->m_diagWeight_X = pMethod->m_diagWeight;
    }

    // MAPD weight used in HMM dispersion estimate
    pMethod->m_mapdWeight_X = setupDoubleParameter("mapd-weight-X", strPrefix,
        params, doc, opts);
    iter = params.find("mapd-weight-X");
    if (iter == params.end())
    {
        pMethod->m_mapdWeight_X = pMethod->m_mapdWeight;
    }

    // Minimum size for a segment
    pMethod->m_minSegmentSize_X = setupIntParameter("min-segment-size-X",
        strPrefix, params, doc, opts);
    iter = params.find("min-segment-size-X");
    if (iter == params.end())
    {
        pMethod->m_minSegmentSize_X = pMethod->m_minSegmentSize;
    }

    // hmmCN_state-X
    str = setupStringParameter("hmmCN_state-X", strPrefix, params, doc, opts);
    iter = params.find("hmmCN_state-X");
    if (iter != params.end())
    {
        ba.assign(str);
        ba.replace(",", " ");
        pMethod->m_vCNState_X.resize(ba.parameterCount());
        for (int k=0; k<ba.parameterCount(); k++)
            {
            pMethod->m_vCNState_X[k] = ba.getParameter(k+1).parseInt();
            }
    }
    else {
        pMethod->m_vCNState_X.resize(pMethod->m_vCNState.size());
        pMethod->m_vCNState_X = pMethod->m_vCNState;
    }

    // hmmCN_mu-X
    str = setupStringParameter("hmmCN_mu-X", strPrefix, params, doc, opts);
    iter = params.find("hmmCN_mu-X");
    if (iter != params.end())
    {
        ba.assign(str);
        ba.replace(",", " ");
        pMethod->m_Mu_X.resize(ba.parameterCount());
        for (int k=0; k<ba.parameterCount(); k++)
            {
            pMethod->m_Mu_X[k] = ba.getParameter(k+1).parseDouble();
            }
    }
    else {
        pMethod->m_Mu_X.resize(pMethod->m_Mu.size());
        pMethod->m_Mu_X = pMethod->m_Mu;
    }

    // hmmCN_sigma-X
    str = setupStringParameter("hmmCN_sigma-X", strPrefix, params, doc, opts);
    iter = params.find("hmmCN_sigma-X");
    if (iter != params.end())
    {
        ba.assign(str);
        ba.replace(",", " ");
        pMethod->m_Sigma_X.resize(ba.parameterCount());
        for (int k=0; k<ba.parameterCount(); k++)
            {
            pMethod->m_Sigma_X[k] = ba.getParameter(k+1).parseDouble();
            }
    }
    else {
        pMethod->m_Sigma_X.resize(pMethod->m_Sigma.size());
        pMethod->m_Sigma_X = pMethod->m_Sigma;
    }

    if (pMethod->m_vCNState_X.size() != pMethod->m_Mu_X.size())
        {
        throw(Except("Mu-X must conform to CNState-X."));
        }

    if (pMethod->m_vCNState_X.size() != pMethod->m_Sigma_X.size())
        {
        throw(Except("Sigma-X must conform to CNState-X."));
        }


    for (int k=0; k<pMethod->m_vCNState_X.size(); k++)
        {
        if (pMethod->m_Sigma_X[k] <= 0.0)
            throw(Except("Sigma-X parameters must be > 0."));

        if ((k > 0) && (pMethod->m_Mu_X[k-1] > pMethod->m_Mu_X[k]))
            throw(Except("Mu-X parameters must be ordered increasing."));
        }

    // MAPD weight used in HMM dispersion estimate
    pMethod->m_hmm_confidence_weight_X = setupDoubleParameter(
            "hmm-confidence-weight-X", strPrefix, params, doc, opts);
    iter = params.find("hmm-confidence-weight-X");
    if (iter == params.end())
    {
        pMethod->m_hmm_confidence_weight_X = pMethod->m_hmm_confidence_weight;
    }
    else if (pMethod->m_hmm_confidence_weight_X <= 0.0 ||
             pMethod->m_hmm_confidence_weight_X >= 1.0)
    {
        throw(Except("HMM Confidence-X weight must be in [0,1]"));
    }

    // And more parameters for Y chromosome (if separate set of parameters is required for Y)
    // Diagonal
    pMethod->m_diagWeight_Y = setupDoubleParameter("diagonal-weight-Y", strPrefix,
        params, doc, opts);
    iter = params.find("diagonal-weight-Y");
    if (iter == params.end())
    {
        pMethod->m_diagWeight_Y = pMethod->m_diagWeight;
    }

    // MAPD weight used in HMM dispersion estimate
    pMethod->m_mapdWeight_Y = setupDoubleParameter("mapd-weight-Y", strPrefix,
        params, doc, opts);
    iter = params.find("mapd-weight-Y");
    if (iter == params.end())
    {
        pMethod->m_mapdWeight_Y = pMethod->m_mapdWeight;
    }


    // Minimum size for a segment
    pMethod->m_minSegmentSize_Y = setupIntParameter("min-segment-size-Y",
        strPrefix, params, doc, opts);
    iter = params.find("min-segment-size-Y");
    if (iter == params.end())
    {
        pMethod->m_minSegmentSize_Y = pMethod->m_minSegmentSize;
    }

    // hmmCN_state-Y
    str = setupStringParameter("hmmCN_state-Y", strPrefix, params, doc, opts);
    iter = params.find("hmmCN_state-Y");
    if (iter != params.end())
    {
        ba.assign(str);
        ba.replace(",", " ");
        pMethod->m_vCNState_Y.resize(ba.parameterCount());
        for (int k=0; k<ba.parameterCount(); k++)
            {
            pMethod->m_vCNState_Y[k] = ba.getParameter(k+1).parseInt();
            }
    }
    else {
        pMethod->m_vCNState_Y.resize(pMethod->m_vCNState.size());
        pMethod->m_vCNState_Y = pMethod->m_vCNState;
    }

    // hmmCN_mu-Y
    str = setupStringParameter("hmmCN_mu-Y", strPrefix, params, doc, opts);
    iter = params.find("hmmCN_mu-Y");
    if (iter != params.end())
    {
        ba.assign(str);
        ba.replace(",", " ");
        pMethod->m_Mu_Y.resize(ba.parameterCount());
        for (int k=0; k<ba.parameterCount(); k++)
            {
            pMethod->m_Mu_Y[k] = ba.getParameter(k+1).parseDouble();
            }
    }
    else {
        pMethod->m_Mu_Y.resize(pMethod->m_Mu.size());
        pMethod->m_Mu_Y = pMethod->m_Mu;
    }

    // hmmCN_sigma-Y
    str = setupStringParameter("hmmCN_sigma-Y", strPrefix, params, doc, opts);
    iter = params.find("hmmCN_sigma-Y");
    if (iter != params.end())
    {
        ba.assign(str);
        ba.replace(",", " ");
        pMethod->m_Sigma_Y.resize(ba.parameterCount());
        for (int k=0; k<ba.parameterCount(); k++)
            {
            pMethod->m_Sigma_Y[k] = ba.getParameter(k+1).parseDouble();
            }
    }
    else {
        pMethod->m_Sigma_Y.resize(pMethod->m_Sigma.size());
        pMethod->m_Sigma_Y = pMethod->m_Sigma;
    }

    if (pMethod->m_vCNState_Y.size() != pMethod->m_Mu_Y.size())
        {
        throw(Except("Mu-Y must conform to CNState-Y."));
        }

    if (pMethod->m_vCNState_Y.size() != pMethod->m_Sigma_Y.size())
        {
        throw(Except("Sigma-Y must conform to CNState-Y."));
        }


    for (int k=0; k<pMethod->m_vCNState_Y.size(); k++)
        {
        if (pMethod->m_Sigma_Y[k] <= 0.0)
            throw(Except("Sigma-Y parameters must be > 0."));

        if ((k > 0) && (pMethod->m_Mu_Y[k-1] > pMethod->m_Mu_Y[k]))
            throw(Except("Mu-Y parameters must be ordered increasing."));
        }

    // MAPD weight used in HMM dispersion estimate
    pMethod->m_hmm_confidence_weight_Y = setupDoubleParameter(
            "hmm-confidence-weight-Y", strPrefix, params, doc, opts);
    iter = params.find("hmm-confidence-weight-Y");
    if (iter == params.end())
    {
        pMethod->m_hmm_confidence_weight_Y = pMethod->m_hmm_confidence_weight;
    }
    else if (pMethod->m_hmm_confidence_weight_Y <= 0.0 ||
             pMethod->m_hmm_confidence_weight_Y >= 1.0)
    {
        throw(Except("HMM Confidence-Y weight must be in [0,1]"));
    }


//Shrinkage
	// true or false
    pMethod->m_shrink = setupBoolParameter("shrink",
        strPrefix, params, doc, opts);

    // shrink-lprec
    str = setupStringParameter("shrink-lprec", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_shrink_lprec.resize(ba.parameterCount());
    for (int k=0; k<ba.parameterCount(); k++)
        {
        pMethod->m_shrink_lprec[k] = ba.getParameter(k+1).parseDouble();
        }

    // Convergence criterion for shrinkage
    pMethod->m_shrink_converge = setupDoubleParameter("shrink-converge",
        strPrefix, params, doc, opts);

    pMethod->m_shrink_downweight_outlier =
        setupBoolParameter("shrink-downweight-outlier",
        strPrefix, params, doc, opts);

    pMethod->m_shrink_downweight_df =
        setupDoubleParameter("shrink-downweight-df",
        strPrefix, params, doc, opts);

    pMethod->m_shrink_downweight_maxiter =
        setupIntParameter("shrink-downweight-maxiter",
        strPrefix, params, doc, opts);

    return pMethod;
}



/**
 * @brief Run the analysis
 */
void CNAnalysisMethodCNCyto2::run()
{
    Verbose::out(1, "MAJOR PROGRESS UPDATE: Running Copy Number Analysis.");
    // An exception is thrown if there is a Setup problem.
    isSetup();

    if (m_pEngine->getOptBool("keep-intermediate-data-local"))
    {
        writeLog2Ratios("before");
    }

    if (m_pEngine->getOptBool("keep-intermediate-data-local"))
    {
        writeLog2Ratios("after");
    }


    setCNCallsToInvalidState();

    determineLocalProbeSets();

    double mapd = computeMAPD();

    vector<int>chromosomes = getChromosomes(getProbeSets());
    int last_chr = chromosomes[chromosomes.size() - 1];

    Verbose::progressBegin(1,"CNAnalysisMethodCNCyto2::run(...) ",
         last_chr, 1, last_chr);

    CNSmoother CopyNumberSmoother(m_minSegmentSize);

    vector<string> vProbeSetNamesRef;
    string referenceFileName = m_pEngine->getOpt("reference-file");

    for (int i=0; i<chromosomes.size(); i++)
    {
        int chr = chromosomes[i];
        Verbose::progressStep(1);

        valarray<double> *mu_ptr;
        valarray<double> *sigma_ptr;
        double diagWeight;
        double mapdWeight;
        double hmm_confidence_weight;
        int numMu;

        if (chr < m_iXChromosome)
        {
            mu_ptr = &m_Mu;
            sigma_ptr = &m_Sigma;
            diagWeight = m_diagWeight;
            mapdWeight = m_mapdWeight;
            hmm_confidence_weight = m_hmm_confidence_weight;
            numMu = m_Mu.size();
        }
        else if (chr == m_iXChromosome)
        {
            mu_ptr = &m_Mu_X;
            sigma_ptr = &m_Sigma_X;
            diagWeight = m_diagWeight_X;
            mapdWeight = m_mapdWeight_X;
            hmm_confidence_weight = m_hmm_confidence_weight_X;
            numMu = m_Mu_X.size();
        }
        else if (chr == m_iYChromosome)
        {
            mu_ptr = &m_Mu_Y;
            sigma_ptr = &m_Sigma_Y;
            diagWeight = m_diagWeight_Y;
            mapdWeight = m_mapdWeight_Y;
            hmm_confidence_weight = m_hmm_confidence_weight_Y;
            numMu = m_Mu_Y.size();
        }

        valarray<double>sigma(numMu);
        for (int i=0; i<numMu; i++)
            {
            sigma[i] = (*sigma_ptr)[i] + mapd * mapdWeight;
            }
        valarray<double> sigma_sq(numMu);
        sigma_sq = sigma * sigma;

        int iProbeSetCount = getProbeSetCount(chr, getProbeSets());
        if (iProbeSetCount == 0) continue; // probably 23 or Y

        int iCNCount = 0;
        int chr_start = getChrBounds(chr, getProbeSets()).first;
        CNProbeSetArray * psets = getProbeSets();
        for (int i=0; i<iProbeSetCount; i++)
        {
            if (psets->at(chr_start + i)->processAsCN()) {iCNCount++;}
        }

        valarray<double> vLog2Ratios(iCNCount);
        int iExpectedCNState = 2;
        if (chr==m_iXChromosome && getExperiment()->hasY()) iExpectedCNState = 1;
        if (chr==m_iYChromosome && getExperiment()->hasY()) iExpectedCNState = 1;
        if (chr==m_iYChromosome && !getExperiment()->hasY()) {iExpectedCNState = 0;}

        // For males, need to find the par boundaries on X
        // because par regions look female due to duplication on Y.
        int nopar_start=0, nopar_end=0;  // initialize to kill the warnings
        if (chr==m_iXChromosome && getExperiment()->hasY())
        {
            bool found_start=false, found_end=false;
            int idx = 0;
            for (int i=0; i<iProbeSetCount; i++)
            {
                if (psets->at(chr_start + i)->processAsCN()) {
                    if (!found_start &&
                        !psets->at(chr_start + i)->isPseudoAutosomalRegion()) {
                        found_start = true;
                        nopar_start = idx;
                    }

                    if (found_start &&
                        psets->at(chr_start + i)->isPseudoAutosomalRegion()) {
                        found_end = true;
                        nopar_end = idx;
                    }

                    idx++;
                }
                if (found_end)
                    break;
            }
        }

        int iCNIndex = 0;
        for (int i=0; i<iProbeSetCount; i++) {
            if (psets->at(chr_start + i)->processAsCN()) {
                vLog2Ratios[iCNIndex] = psets->at(chr_start+i)->getLog2Ratio();
                iCNIndex++;
            }
        }

        if (iCNIndex == 0) {continue;}

        TransitionMatrix transMatrix(numMu, diagWeight);
        CopyNumberHMM cnHMM(transMatrix);

        valarray<int> vCNStates(0,iCNCount);
        valarray<float> vCNConfidences((float)0, iCNCount);

        valarray<double> log2r_shrunk(vLog2Ratios);

        if (m_shrink) {
            Wavelet * W = wavelet('h',1);
            MODE mode=MODE_ZEROPAD;

            wavelet_sandbox(log2r_shrunk, W, m_shrink_lprec, 1.0,
                mode, m_shrink_converge, m_shrink_downweight_outlier,
                m_shrink_downweight_df, m_shrink_downweight_maxiter);

            free_wavelet(W);
            }

        cnHMM.add_data(log2r_shrunk, *mu_ptr, sigma_sq);
        vCNStates = cnHMM.solve_map();

        CopyNumberSmoother.smooth(vCNStates);

        compute_confidence_mu_fixed(iExpectedCNState, vLog2Ratios,
            vCNStates, vCNConfidences, *mu_ptr, sigma,
            hmm_confidence_weight);

        // ParAlagous Regions (PAR) in females are duplicated on males.
        // The duplication of the markers on Y will make a male X look
        // like a female XX.
        if (chr==m_iXChromosome && getExperiment()->hasY())
            {
            valarray<float> tmp_conf((float)0, iCNCount);

            compute_confidence_mu_fixed(2, vLog2Ratios, vCNStates,
                tmp_conf, *mu_ptr, sigma, hmm_confidence_weight);

            for (int ipos=0; ipos<nopar_start; ipos++) {
                vCNConfidences[ipos] = tmp_conf[ipos];
                }

            for (int ipos=nopar_end; ipos<iCNCount; ipos++) {
                vCNConfidences[ipos] = tmp_conf[ipos];
                }
            }
        if(m_pEngine->getOptBool("cyto2"))
        {
            copyNumberStatePostProcessing(chr, vCNStates, vCNConfidences);
        }
        else
        {
            copyNumberStatePostProcessing(chr, vCNStates, vCNConfidences,
            log2r_shrunk);
        }
    }

    Verbose::progressEnd(1, "Done");
}

// Check the means for correctness, return offending row index if error,
// return -1 if OK.
int CNAnalysisMethodCNCyto2::hmmPriorsVerify(
                                    const valarray<double>& mu_arr,
                                    int nRows,
                                    int nCols
                                    )
{
    for (int i = 0; i < nRows; i++) {
        for (int j = 1; j < nCols; j++) {
            if (mu_arr[i*nCols + j] <= mu_arr[i*nCols + j - 1]) {
                return i;    // error
            }
        }
    }
    return -1;   // OK
}

/**
 * @brief Anneal the CN calls and load them into the CNProbeSet vector.
 * @param int - The chromosome to process.
 * @param vector<float>& - The log2 ration vector
 * @param vector<int>& - The CN calls vector
 */

// This seems to assigne the copy number state to a probe set.
// The trick with assuming a sort makes this algorithm fragile.
void CNAnalysisMethodCNCyto2::copyNumberStatePostProcessing(
        int iChromosome,
        valarray<int> & vCNStates,
        valarray<float> & vCNConfidences,
        valarray<double> & shrunk_log2r)
    {
    CNProbeSetArray * psets = getProbeSets();
    pair<int,int>chr_bounds = getChrBounds(iChromosome, psets);


    int pos = 0;
    for (int i=chr_bounds.first; i<chr_bounds.second; i++) {
        CNProbeSet * pobjProbeSet = psets->at(i);
        if (!pobjProbeSet->processAsCN()) continue;
        pobjProbeSet->setCNState(vCNStates[pos]);
        pobjProbeSet->setCNConfidence(vCNConfidences[pos]);
        pobjProbeSet->setLog2RatioMedianSmooth((float)shrunk_log2r[pos]);
        pos++;
        }


/* ye olde loope
    int iIndex = 0;
    for (int i=0; i<getProbeSets()->size(); i++) {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(i);
        if (pobjProbeSet->getChromosome() == iChromosome) {
            if (!pobjProbeSet->processAsCN()) continue;
            pobjProbeSet->setCNState(vCNStates[iIndex]);
            pobjProbeSet->setCNConfidence(vCNConfidences[iIndex]);
//          pobjProbeSet->setLog2RatioMedianSmooth((float)shrunk_log2r[iIndex]);
            iIndex++;
            }
        }
*/
    }


void CNAnalysisMethodCNCyto2::copyNumberStatePostProcessing(
        int iChromosome,
        valarray<int> & vCNStates,
        valarray<float> & vCNConfidences)
        {
        CNProbeSetArray * psets = getProbeSets();
        pair<int,int>chr_bounds = getChrBounds(iChromosome, psets);


    int pos = 0;
    for (int i=chr_bounds.first; i<chr_bounds.second; i++) {
        CNProbeSet * pobjProbeSet = psets->at(i);
        if (!pobjProbeSet->processAsCN()) continue;
        pobjProbeSet->setCNState(vCNStates[pos]);
        pobjProbeSet->setCNConfidence(vCNConfidences[pos]);
        pos++;
        }

}


/**
 * @brief Performs the Confidence algorithm on log2 ratios of a sample on
 * a per chromosome basis.
 *
 * @param int.................-.Expected CN State.
 * @param samplr              - Sample log2 ratios.
 * @param confidence  .       - in/out: CN state confidences.
 * @param mu                  - in/out: Mu values.
 * @param sigma               - in/out: Sigma values.
 */


void CNAnalysisMethodCNCyto2::compute_confidence_mu_fixed(
        const int iExpectedCNState,
        const valarray<double> & samplr,
        valarray<int> & state,
        valarray<float> & confidence,
        const valarray<double> mu,
        const valarray<double> sigma,
        const double hmm_weight)
    {
    int N=samplr.size();
    for (int i=0; i<N; i++) {
        confidence[i] = compute_confidence(iExpectedCNState, samplr[i],
            state[i], mu, sigma, hmm_weight);
        }
    }




/**
 * @brief Performs the Confidence algorithm on log2 ratios of a sample on
 * a per chromosome basis.
 *
 * @param int.................-.Expected CN State.
 * @param samplr              - Sample log2 ratios.
 * @param confidence  .       - in/out: CN state confidences.
 * @param mu                  - in/out: Mu values.
 * @param sigma               - in/out: Sigma values.
 */


void CNAnalysisMethodCNCyto2::compute_confidence_mu_variable(
        const int iExpectedCNState,
        const valarray<double> & samplr,
        valarray<int> & state,
        valarray<float> & confidence,
        const valarray<double> mu,
        const valarray<double> sigma,
        const double hmm_weight)
    {
    int N=samplr.size();
    int P = sigma.size();

    valarray<double>mu_marker(P);

    for (int i=0; i<N; i++) {
        mu_marker = mu[slice(i*P,P,1)];
        confidence[i] = compute_confidence(iExpectedCNState, samplr[i],
            state[i], mu_marker, sigma, hmm_weight);
        }
    }



// for computing confidence one marker at a time
float CNAnalysisMethodCNCyto2::compute_confidence(
    const int iExpectedCNState,
    const double logratio,
    const int cnstate,
    const valarray<double> mu,
    const valarray<double> sigma,
    const double hmm_weight)
    {
    int K=mu.size();
    float confidence;

    // Anything outside this range is an outlier.  A reasonable
    // range to leave hard coded inside this specific purpose function.
    // If it is an outlier, return that all possibilities are equal.
    if ( (logratio < mu[0] - 5.0*sigma[0]) ||
         (logratio > mu[K-1] + 5.0*sigma[K-1]) ) {
        confidence = 1.0/K;
        double hmm_adj = (cnstate == iExpectedCNState) ? 1.0 : 0.0;
        confidence = hmm_adj*hmm_weight + confidence*(1.0 - hmm_weight);
        return 1.0 - confidence;
        }


    // Probabilities for each possible state not summing to 1.
    valarray<double> probs(K);


    for (int k=0; k<K; k++) {
        double prec = 1.0/(sigma[k]*sigma[k]);
        double diff = logratio - mu[k];
        diff = prec*diff*diff;

        if (diff > 100.0)
            {
            probs[k] = 0.0;
            continue;
            }

        probs[k] = sqrt(prec)*exp(-0.5*diff);
        }

    confidence = (float)(probs[iExpectedCNState]/probs.sum());

    // The confidence is a weighted confidence of the hmm, which is always
    // certain it is right , and the multinomial derived from the likelihood.
    double hmm_adj = (cnstate == iExpectedCNState) ? 1.0 : 0.0;
    confidence = hmm_adj*hmm_weight + confidence*(1.0 - hmm_weight);

    // a quick fix to flip interpretation of confidence
    return 1.0 - confidence;
    }


/**
 * @brief Compute MAPD
 */
double CNAnalysisMethodCNCyto2::computeMAPD()
    {
    Verbose::out(1, "CNAnalysisMethodCNCyto2::MAPD(...) start");

    // Load up the absolute differences in an array.
    int N = getProbeSets()->size();
    vector<double> log2ratio(N-1);

    double previous = getProbeSets()->at(0)->getLog2Ratio();
    for (int i=0; i<N-1; i++)
        {
        double current = getProbeSets()->at(i+1)->getLog2Ratio();
        log2ratio[i] = abs(current-previous);
        previous = current;
        }

    // don't be picky about computing the median with these many numbers
    nth_element(log2ratio.begin(), log2ratio.begin()+(N+1)/2, log2ratio.end());

    return log2ratio[(N+1)/2];
    }

void CNAnalysisMethodCNCyto2::determineLocalProbeSets()
{
        if(m_bLocalProbeSetsDetermined)
        {
                return;
        }
        const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
        int iNumberOfProbeSets=CNAnalysisMethod::getProbeSets()->getCount();
        for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
        {
                if( CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAll() )
                {                        
                    if (isCytoScanHD && (!CNAnalysisMethod::getProbeSets()->getAt(iIndex)->processAsCN())) {continue;}
                    getProbeSets()->add( CNAnalysisMethod::getProbeSets()->getAt(iIndex));
                }
        }
        m_bLocalProbeSetsDetermined=true;
}


void CNAnalysisMethodCNCyto2::setCNCallsToInvalidState()
{

        int iNumberOfProbeSets=CNAnalysisMethod::getProbeSets()->getCount();
        for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
        {

            CNAnalysisMethod::getProbeSets()->getAt(iIndex)->setCNState(-1);
        }
}
