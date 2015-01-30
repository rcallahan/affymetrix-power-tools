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
 * @file CNAnalysisMethodCN.cpp
 *
 * @brief This file contains the CNAnalysisMethodCN class members.
 */

// Microsoft seems to defeat some of the purpose of using <cmath> leaving
// certain things undefined such as "M_PI".
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include "copynumber/CNAnalysisMethodCN.h"
//
#include "stats/stats.h"
#include "util/AffxMultiDimensionalArray.h"
//

/**
 * Constructor.
 */
CNAnalysisMethodCN::CNAnalysisMethodCN()
{
    m_dTransMatDecay = 0.0;
    m_uiEMIterations = 0;
    m_dEMConvergentThreshold = 0.0;
    m_iNormalState = 0;
    m_iPostCNFitMaxOutlierRemoveRunSize = 0;
    m_fwdOnly = 0;
    m_inormalStateMinObservations = 0;
    m_iSmoothOutliers = 0;
    m_transTypeStat = 0;
    m_updatePrior = false;
    m_updateMu = false;
    m_updateSigma = false;
    m_updateTransMat = false;
    m_scaled = false;
}

CNAnalysisMethodCN::~CNAnalysisMethodCN() {}

std::string CNAnalysisMethodCN::getType() {
  return "cn-state";
}
std::string CNAnalysisMethodCN::getDescription() {
  return "CopyNumber CNState";
}

std::string CNAnalysisMethodCN::getVersion() {
  return CNAnalysisMethodCN_Version;
}


AffxString CNAnalysisMethodCN::getName() {
  return getType();
}

void CNAnalysisMethodCN::setEngine(BaseEngine* p)
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
  for (unsigned int i = 0; (i < m_vMu.size()); i++)
    {
      if (i > 0) {str += ",";}
      if (i == 1)
      {
          if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
          {
            if (i == 1) {m_pEngine->setOpt("yTarget", ToStr(m_vMu[i]));}
          }
          else
          {
            m_pEngine->setOpt("yTarget", "-0.56747");
          }
      }
      str += ::getDouble(m_vMu[i], 6);
    }
  p->setOpt("hmmCN_mu",  str);
}

bool CNAnalysisMethodCN::isSegmentTypeAnalysis() {
  return false;
}

void CNAnalysisMethodCN::setup(
                               bool updatePrior,
                               bool updateMu,
                               bool updateSigma,
                               bool updateTransMat,
                               bool fwdOnly,
                               bool scaled,
                               bool transTypeStat
                               )
{
  m_updatePrior    = updatePrior;
  m_updateMu       = updateMu;
  m_updateSigma    = updateSigma;
  m_updateTransMat = updateTransMat;
  m_fwdOnly        = fwdOnly;
  m_scaled         = scaled;
  m_transTypeStat  = transTypeStat;
}

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNAnalysisMethodCN::explainSelf()
{
    CNAnalysisMethodCN obj;
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
std::vector<SelfDoc::Opt> CNAnalysisMethodCN::getDefaultDocOptions()
{
  std::vector<SelfDoc::Opt> opts;

  // SelfDoc::Opt(name, type, value, default, min, max, description)

  SelfDoc::Opt hmmCN_state = {"hmmCN_state", SelfDoc::Opt::String, "0,1,2,3,4", "0,1,2,3,4", "", "", "CN HMM State"};
  opts.push_back(hmmCN_state);

  SelfDoc::Opt hmmCN_prior_prob = {"hmmCN_prior_prob", SelfDoc::Opt::String, "0.2,0.2,0.2,0.2,0.2", "0.2,0.2,0.2,0.2,0.2", "", "", "CN HMM Prior Probability"};
  opts.push_back(hmmCN_prior_prob);

  // gc-correction false = -2,-0.552,0,0.339,0.543
  // gc-correction true  = -2,-0.533,0,0.363,0.567
  SelfDoc::Opt hmmCN_mu = {"hmmCN_mu", SelfDoc::Opt::String, "-2,-0.533,0,0.363,0.567", "-2,-0.533,0,0.363,0.567", "", "", "CN HMM Mu"};
  opts.push_back(hmmCN_mu);

  SelfDoc::Opt hmmCN_sigma = {"hmmCN_sigma", SelfDoc::Opt::String, "0.2,0.2,0.2,0.2,0.2", "0.2,0.2,0.2,0.2,0.2", "", "", "CN HMM Sigma"};
  opts.push_back(hmmCN_sigma);

  SelfDoc::Opt hmmCN_TransitionDecay = {"hmmCN_TransitionDecay", SelfDoc::Opt::Double, "1e+9", "1e+9", "NA", "NA", "CN HMM Transition Decay"};
  opts.push_back(hmmCN_TransitionDecay);

  SelfDoc::Opt hmmCN_StateEstimationMethod = {"hmmCN_StateEstimationMethod", SelfDoc::Opt::String, "EM", "EM", "", "", "CN HMM StateEstimation Method"};
  opts.push_back(hmmCN_StateEstimationMethod);

  SelfDoc::Opt hmmCN_EMIterations = {"hmmCN_EMIterations", SelfDoc::Opt::Integer, "1", "1", "NA", "NA", "CN HMM Estimation Method Iterations"};
  opts.push_back(hmmCN_EMIterations);

  SelfDoc::Opt hmmCN_EMConvergenceThreshold = {"hmmCN_EMConvergenceThreshold", SelfDoc::Opt::Double, "0.0001", "0.0001", "NA", "NA", "CN HMM Estimation Method Convergence Threshold"};
  opts.push_back(hmmCN_EMConvergenceThreshold);

  SelfDoc::Opt hmmCN_NormalState = {"hmmCN_NormalState", SelfDoc::Opt::Integer, "2", "2", "NA", "NA", "CN HMM Normal State"};
  opts.push_back(hmmCN_NormalState);

  SelfDoc::Opt hmmCN_ForwardOnly = {"hmmCN_ForwardOnly", SelfDoc::Opt::Integer, "0", "0", "NA", "NA", "CN HMM Forward Only"};
  opts.push_back(hmmCN_ForwardOnly);

  SelfDoc::Opt hmmCN_NormalStateMinObservations = {"hmmCN_NormalStateMinObservations", SelfDoc::Opt::Integer, "2", "2", "NA", "NA", "CN HMM Normal State Min Observations"};
  opts.push_back(hmmCN_NormalStateMinObservations);

  SelfDoc::Opt hmmCN_SmoothOutliers = {"hmmCN_SmoothOutliers", SelfDoc::Opt::Integer, "1", "1", "NA", "NA", "CN HMM Smooth Outliers"};
  opts.push_back(hmmCN_SmoothOutliers);

  SelfDoc::Opt hmmCN_TransTypeStat = {"hmmCN_TransTypeStat", SelfDoc::Opt::Integer, "0", "0", "NA", "NA", "CN HMM Trans Type Stat"};
  opts.push_back(hmmCN_TransTypeStat);

  SelfDoc::Opt PostCNFitMaxOutlierRemoveRunSize = {"PostCNFitMaxOutlierRemoveRunSize", SelfDoc::Opt::Integer, "1", "1", "NA", "NA", "Post CN Fit Maximum Outlier Remove Run Size"};
  opts.push_back(PostCNFitMaxOutlierRemoveRunSize);

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
SelfCreate* CNAnalysisMethodCN::newObject(std::map<std::string,std::string>& params)
{
    SelfDoc doc = explainSelf();
    std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
    CNAnalysisMethodCN* pMethod = new CNAnalysisMethodCN();
    std::string strPrefix = getPrefix();

    AffxString str;
    AffxByteArray ba;

    str = setupStringParameter("hmmCN_state", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_vCNState.resize(ba.parameterCount());
    for (int iIndex = 1; (iIndex <= ba.parameterCount()); iIndex++)
    {
        pMethod->m_vCNState[iIndex - 1] = ba.getParameter(iIndex).parseInt();
    }

    str = setupStringParameter("hmmCN_prior_prob", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_vPriorProb.resize(ba.parameterCount());
    for (int iIndex = 1; (iIndex <= ba.parameterCount()); iIndex++)
    {
        pMethod->m_vPriorProb[iIndex - 1] = ba.getParameter(iIndex).parseDouble();
    }

    str = setupStringParameter("hmmCN_mu", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_vMu.resize(ba.parameterCount());
    for (int iIndex = 1; (iIndex <= ba.parameterCount()); iIndex++)
    {
        pMethod->m_vMu[iIndex - 1] = ba.getParameter(iIndex).parseDouble();
    }

    str = setupStringParameter("hmmCN_sigma", strPrefix, params, doc, opts);
    ba.assign(str);
    ba.replace(",", " ");
    pMethod->m_vSigma.resize(ba.parameterCount());
    for (int iIndex = 1; (iIndex <= ba.parameterCount()); iIndex++)
    {
        pMethod->m_vSigma[iIndex - 1] = ba.getParameter(iIndex).parseDouble();
    }

    pMethod->m_dTransMatDecay = setupDoubleParameter("hmmCN_TransitionDecay", strPrefix, params, doc, opts);
    pMethod->m_strStateEstMethod = setupStringParameter("hmmCN_StateEstimationMethod", strPrefix, params, doc, opts);
    pMethod->m_fwdOnly = setupIntParameter("hmmCN_ForwardOnly", strPrefix, params, doc, opts);
    pMethod->m_inormalStateMinObservations = setupIntParameter("hmmCN_NormalStateMinObservations", strPrefix, params, doc, opts);
    pMethod->m_iSmoothOutliers = setupIntParameter("hmmCN_SmoothOutliers", strPrefix, params, doc, opts);
    pMethod->m_transTypeStat = setupIntParameter("hmmCN_TransTypeStat", strPrefix, params, doc, opts);
    pMethod->m_uiEMIterations = setupIntParameter("hmmCN_EMIterations", strPrefix, params, doc, opts);
    pMethod->m_dEMConvergentThreshold = setupDoubleParameter("hmmCN_EMConvergenceThreshold", strPrefix, params, doc, opts);
    pMethod->m_iNormalState = setupIntParameter("hmmCN_NormalState", strPrefix, params, doc, opts);
    pMethod->m_iPostCNFitMaxOutlierRemoveRunSize = setupIntParameter("PostCNFitMaxOutlierRemoveRunSize", strPrefix, params, doc, opts);

    if (pMethod->m_vCNState.size() != pMethod->m_vPriorProb.size()) {throw(Except("PriorProb parameter must be the same size as the CNState parameter."));}
    if (pMethod->m_vCNState.size() != pMethod->m_vMu.size()) {throw(Except("Mu parameter must be the same size as the CNState parameter."));}
    if (pMethod->m_vCNState.size() != pMethod->m_vSigma.size()) {throw(Except("Sigma parameter must be the same size as the CNState parameter."));}
    double dSumPriorProb = 0;
    for (int iIndex = 0; (iIndex < (int)pMethod->m_vCNState.size()); iIndex++)
    {
        dSumPriorProb += pMethod->m_vPriorProb[iIndex];
        if (pMethod->m_vPriorProb[iIndex] <= 0) {throw(Except("PriorProb parameters must be greater than or equal to zero."));}
        if (pMethod->m_vSigma[iIndex] <= 0) {throw(Except("Sigma parameters must be greater than or equal to zero."));}
        if ((iIndex > 0) && (pMethod->m_vMu[iIndex - 1] > pMethod->m_vMu[iIndex])) {throw(Except("Mu parameters must be sorted and increasing in order."));}
    }
    for (int iIndex = 0; (iIndex < (int)pMethod->m_vPriorProb.size()); iIndex++)
    {
        pMethod->m_vPriorProb[iIndex] /= dSumPriorProb;
    }
    if (pMethod->m_iNormalState > (int)(pMethod->m_vCNState.size() - 1)) {throw(Except("Normal State parameter must be less than the number of CN States specified."));}

    return pMethod;
}

/**
 * @brief Run the analysis
 */
void CNAnalysisMethodCN::run()
{
    isSetup();

    int iLastChromosome = getProbeSets()->at(getProbeSets()->size() - 1)->getChromosome();
    int iLastAutosomeChromosome = getLastAutosomeChromosome();
    AffxMultiDimensionalArray<float> arAutosomeMedians(iLastAutosomeChromosome);
//    Verbose::progressBegin(1, "CNAnalysisMethodCN::run(...) ", iLastChromosome, 1, iLastChromosome);
    for (int iChromosome = 1; (iChromosome <= iLastChromosome); iChromosome++)
    {
//        Verbose::progressStep(1);
        int iProbeSetCount = getProbeSetCount(iChromosome, getProbeSets());
        if (iProbeSetCount > 0)
        {
            std::vector<int> vPositions(iProbeSetCount);
            std::vector<float> vLog2Ratios(iProbeSetCount);
            std::vector<int> vCNStates(iProbeSetCount);
            int iIndex = 0;
            for (int iProbeSetIndex = 0; (iProbeSetIndex < (int)getProbeSets()->size()); iProbeSetIndex++)
            {
                CNProbeSet* pobjProbeSet = getProbeSets()->at(iProbeSetIndex);
                if (pobjProbeSet->getChromosome() == iChromosome)
                {
                    vPositions[iIndex] = pobjProbeSet->getPosition();
                    vLog2Ratios[iIndex] = pobjProbeSet->getLog2Ratio();
                    vCNStates[iIndex] = 0;
                    iIndex++;
                }
            }
            std::vector<double> vDiagProb(5, 0.94);
            setup(true, true, true, true, (m_fwdOnly == 1), true, (m_transTypeStat == 1));
            diagProb(vDiagProb);
            std::vector<double> vStateMargProb; // unused output
            hmm(vPositions, vLog2Ratios, vCNStates, m_vPriorProb, m_vMu, m_vSigma, m_dTransMatDecay, m_strStateEstMethod, vStateMargProb, m_uiEMIterations, m_dEMConvergentThreshold, false, m_iNormalState);
            copyNumberStatePostProcessing(iChromosome, vLog2Ratios, vCNStates);
        }
    }
//    Verbose::progressEnd(1, "Done");
}

/**
 * @brief Anneal the CN calls and load them into the CNProbeSet vector.
 * @param int - The chromosome to process.
 * @param std::vector<float>& - The log2 ration vector
 * @param std::vector<int>& - The CN calls vector
 */
void CNAnalysisMethodCN::copyNumberStatePostProcessing(int iChromosome, std::vector<float>& vLog2Ratios, std::vector<int>& vCNStates)
{
    if (m_iPostCNFitMaxOutlierRemoveRunSize == 1)
    {
        std::vector<int> vDiff((int)vLog2Ratios.size() - 1);
        for (int iIndex = 1; (iIndex < (int)vLog2Ratios.size()); iIndex++)
        {
            vDiff[iIndex - 1] = (vCNStates[iIndex] - vCNStates[iIndex - 1]);
        }
        for (int iIndex = 1; (iIndex < (int)vDiff.size()); iIndex++)
        {
            if ((vDiff[iIndex - 1] == -vDiff[iIndex]) && (vDiff[iIndex - 1] != 0))
            {
                vCNStates[iIndex] = vCNStates[iIndex - 1];
            }
        }
    }
    int iIndex = 0;
    for (int iProbeSetIndex = 0; (iProbeSetIndex < (int)getProbeSets()->size()); iProbeSetIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iProbeSetIndex);
        if (pobjProbeSet->getChromosome() == iChromosome)
        {
            pobjProbeSet->setCNState(vCNStates[iIndex]);
            iIndex++;
        }
    }
}

/**
 * @brief Gives the CNAnalysisMethodCN object the diagonal probabilities it
 * should use.
 * @param dpvec  - Diagonal probabilities.
 */
void
CNAnalysisMethodCN::diagProb(const vector<double>& dpvec) {
  m_diagProb.clear();
  m_diagProb = dpvec;
}



/**
 * @brief Compute the observed log-likelihoods 'aka' the emission probabilities.
 * @param samplr         - Sample log2 ratios.
 * @param sumState       - Number of states.
 * @param updated_mu     - Mean values per state.
 * @param updated_sigma  - standard deviation per state.
 * @param log_lik_obs    - Observed log likelihood.
 */
void CNAnalysisMethodCN::obsloglik(
    //Input
    const VecCNLR&          samplr,
    const int             numState,
    const vector<double>&         updated_mu,
    const vector<double>&         updated_sigma,
    //Output
    vector< vector<double> >&     log_lik_obs
    )
{
  // Initialization
  int numobs =     (int)samplr.size();         // Number of observed symbols
  double ln_sqrt_2pi =    log(sqrt(2*M_PI));

  //Verbose::out(3,"CNAnalysisMethodCN::hmm:obsloglik - Estimating Log Likelihood... ");

  for (int n = 0; n < numobs; n++) {
      for (int stateI = 0; stateI < numState; stateI++) {
        double dTMP = (samplr[n] - updated_mu[stateI])/updated_sigma[stateI];
        log_lik_obs[stateI][n] = -(ln_sqrt_2pi + log(updated_sigma[stateI]) +0.5 * (dTMP * dTMP));
        //Verbose::out(10,"CNAnalysisMethodCN::hmm:obsloglik - Observation: " + ToStr(n) + " state: " + ToStr(stateI) + " log_lik_obs : " + ToStr(log_lik_obs[stateI][n]));
      } //stateI
  } //numobs
}

/**
 * @brief Compute the forward log likelihood.
 * @param dist_covariate
 * @param log_priorProb         -
 * @param updated_log_priorProb -
 * @param log_lik_obs           -
 * @param transMatDecay         -
 * @param trans_type_stat       -
 * @param log_transmat_stat     -
 * @param alpha                 -
 * @param beta                  -
 * @param gamma                 -
 * @param gamma_state_sum       -
 * @param xi_sum                -
 * @param previous_loglik       -
 * @return The forward log likelihood.
 */
double CNAnalysisMethodCN::fwdback(
    //Input
    const vector<double>&         dist_covariate,
    const vector<double>&         log_priorProb,
    const vector<double>&         updated_log_priorProb,
    const vector<vector<double> >&     log_lik_obs,
    const double                    transMatDecay,
    const bool                      trans_type_stat,
    const vector< vector<double> >& log_transmat_stat,
    //Output
    vector< vector<double> >&     alpha,
    vector< vector<double> >&     beta,
    vector< vector<double> >&     gamma,
    vector<double>&                 gamma_state_sum,
    vector< vector<double> >&     xi_sum,
    const double             previous_loglik
    )
{
  // Initialization
  double decay_factor, log_decay_factor, prior_mult_exp_fn, trans_elem;
  int numobs =   (int)log_lik_obs[0].size();
  int numState = (int)log_lik_obs.size();


  //Verbose::out(3,"CNAnalysisMethodCN::hmm:fwdback - Computing forward-back coefficients... ");
  // Verbose::out(10,"CNAnalysisMethodCN::hmm:fwdback - fwd_only: " + ToStr(m_fwdOnly) +  " scaled: " + ToStr(m_scaled));

  //Error checking: Make sure sum(priors) =1
  double prior_sum = -FLT_MAX;
  for(int stateI = 0; stateI < numState; stateI++) {
    prior_sum = ADD_LOG_PROB(prior_sum, updated_log_priorProb[stateI]);
  }
  CNERR_CHECK((fabs(exp(prior_sum)-1)<1e-5),"CNAnalysisMethodCN::hmm:fwdback - Sum of priors: " + ToStr(exp(prior_sum)) + " does not equal 1");

  //Error checking: Make sure sum(transitions) = 1: only necessary in case of stationary transition matrix
  if (trans_type_stat) {
    for(int stateI = 0; stateI < numState; stateI++) {
      double trans_sum = -FLT_MAX;
      for(int stateJ = 0; stateJ < numState; stateJ++) {
        trans_sum = ADD_LOG_PROB(trans_sum, log_transmat_stat[stateI][stateJ]);
      }
      if(fabs(exp(trans_sum)-1)>1e-5) {
        for(int stateJ = 0; stateJ < numState; stateJ++) {
          //Verbose::out(3,"CNAnalysisMethodCN::hmm:fwdback: log_transmat_stat[" + ToStr(stateI) + "][" + ToStr(stateJ) + "] = " + ToStr(log_transmat_stat[stateI][stateJ]));
        }
        CNERR_CHECK((fabs(exp(trans_sum)-1)<1e-5),"CNAnalysisMethodCN::hmm:fwdback - Sum of transitions from state " + ToStr(stateI) + " is " + ToStr(exp(trans_sum)) + " and does not equal 1");
      }
    }
  }

  // Step A: FORWARD of fwdback
  // NOTE: Alpha/Beta/Gamma are stored in log space
  for(int stateI = 0; stateI < numState; stateI++) {
    // Step 1: Initialization of alpha
    alpha[stateI][0] = ADD_LOG_PROB(updated_log_priorProb[stateI], log_lik_obs[stateI][0]);
  }

  //Step 2-3: Induction & Termination
  for(int n = 1; n < numobs; n++) {
    log_decay_factor = -(1/transMatDecay)*dist_covariate[n-1];
    decay_factor = exp(log_decay_factor);
    for(int stateI = 0; stateI < numState; stateI++) {
      alpha[stateI][n] = -FLT_MAX;
      prior_mult_exp_fn = (m_updateTransMat) ? (updated_log_priorProb[stateI]+log(1-decay_factor)) : (log_priorProb[stateI]+log(1-decay_factor));
      for(int stateJ = 0; stateJ < numState; stateJ++) {
        if (trans_type_stat) {
          trans_elem = log_transmat_stat[stateJ][stateI];
        } else {
          trans_elem = prior_mult_exp_fn;
          if(stateJ == stateI) {
            trans_elem = ADD_LOG_PROB(trans_elem, log_decay_factor);
          }
        } // if-else trans_type_stat
        alpha[stateI][n] = ADD_LOG_PROB(alpha[stateI][n], alpha[stateJ][n-1] + trans_elem);
      } // stateJ
      alpha[stateI][n] += log_lik_obs[stateI][n];  // This term is a multiplier, so can happen outside of the loop over stateJ
    } // stateI
  } //numobs

  // Compute full log likelihood
  double fwd_loglik = -FLT_MAX;
  for(int stateI = 0; stateI < numState; stateI++) {
    fwd_loglik = ADD_LOG_PROB(fwd_loglik, alpha[stateI][numobs-1]);
  }
  //Verbose::out(10,"CNAnalysisMethodCN::hmm:fwdback - Log(Probability of the observation sequence given a model): " + ToStr(fwd_loglik));


  // Step B: BACKWARD of fwdback
  // Step 1: Initialization of beta
  for(int stateI = 0; stateI < numState; stateI++) {
    beta[stateI][numobs-1] = 0;
  } // stateI

  for(int stateI = 0; stateI < numState; stateI++) {
      //Verbose::out(10,"CNAnalysisMethodCN::hmm:fwdback(back) - Observation: " + ToStr(numobs-1) + ", state " + ToStr(stateI)  +  ", beta(scaled) = " + ToStr(beta[stateI][numobs-1]));
      //Verbose::out(10,"CNAnalysisMethodCN::hmm:fwdback(back) - Observation: " + ToStr(numobs-1) + ", state " + ToStr(stateI)  +  ", gamma(scaled) = " + ToStr(gamma[stateI][numobs-1]));
  } //stateI

  // Step 2 -3: Induction + Termination
  for(int n = numobs-2; n >=0; n--) {
    log_decay_factor = -(1/transMatDecay)*dist_covariate[n];
    decay_factor = exp(log_decay_factor);
    for(int stateI = 0; stateI < numState; stateI++) {
      beta[stateI][n] = -FLT_MAX;
      prior_mult_exp_fn = (m_updateTransMat) ? (updated_log_priorProb[stateI]+log(1-decay_factor)) : (log_priorProb[stateI]+log(1-decay_factor));
      for(int stateJ = 0; stateJ < numState; stateJ++) {
        if (trans_type_stat) {
          trans_elem = log_transmat_stat[stateI][stateJ];
        } else {
          trans_elem = prior_mult_exp_fn;
          if(stateJ == stateI) {
            trans_elem = ADD_LOG_PROB(trans_elem, log_decay_factor);
          }
        }
        beta[stateI][n] = ADD_LOG_PROB(beta[stateI][n], beta[stateJ][n+1] + trans_elem + log_lik_obs[stateJ][n]);
      } //stateJ
    } //stateI
  } //nobs

  // Compute gamma, gamma_sum & xi_sum
  vector< vector<double> > xi_temp(numState);
  for(int stateI = 0; stateI < numState; stateI++) {
    gamma_state_sum[stateI] = -FLT_MAX;
    xi_temp[stateI].resize(numState);
    for(int stateJ = 0; stateJ < numState; stateJ++) {
      xi_sum[stateI][stateJ] = -FLT_MAX;
    }
  }
  double min_gamma_sum=-FLT_MAX;
  double max_gamma_sum=-FLT_MAX;
  for(int n = 0; n < numobs; n++) {
    double gamma_sum = -FLT_MAX;
    for(int stateI = 0; stateI < numState; stateI++) {
      gamma[stateI][n] = alpha[stateI][n] + beta[stateI][n];
      gamma_sum = ADD_LOG_PROB(gamma_sum, gamma[stateI][n]);
    }
    if(n==0) {
      min_gamma_sum = max_gamma_sum = gamma_sum;
    } else {
      // A sanity-check on forward & backward algorithms - sum of gamma should be the same for every n
      if(gamma_sum > max_gamma_sum) {
        max_gamma_sum = gamma_sum;
      }
      if(gamma_sum < min_gamma_sum) {
        min_gamma_sum = gamma_sum;
      }
    }
    // Scale and sum the gammas
    for(int stateI = 0; stateI < numState; stateI++) {
      gamma[stateI][n] -= gamma_sum;
      gamma_state_sum[stateI] = ADD_LOG_PROB(gamma_state_sum[stateI], gamma[stateI][n]);
    }

    // Computation of xi is for the re-estimation component
    if(n < (numobs-1)) {
      double xi_scale_factor = -FLT_MAX;
      log_decay_factor = -(1/transMatDecay)*dist_covariate[n];
      decay_factor = exp(log_decay_factor);
      for(int stateI = 0; stateI < numState; stateI++) {
        prior_mult_exp_fn = (m_updateTransMat) ? (updated_log_priorProb[stateI]+log(1-decay_factor)) : (log_priorProb[stateI]+log(1-decay_factor));
        for(int stateJ = 0; stateJ < numState; stateJ++) {
          if (trans_type_stat) {
            trans_elem = log_transmat_stat[stateI][stateJ];
          } else {
            trans_elem = prior_mult_exp_fn;
            if(stateJ == stateI) {
              trans_elem = ADD_LOG_PROB(trans_elem, log_decay_factor);
            }
          }
          xi_temp[stateI][stateJ] = alpha[stateI][n] + trans_elem + beta[stateJ][n+1] + log_lik_obs[stateJ][n+1];
          xi_scale_factor = ADD_LOG_PROB(xi_scale_factor, xi_temp[stateI][stateJ]);
        }
      }
      // scale and sum
      for(int stateI = 0; stateI < numState; stateI++) {
        for(int stateJ = 0; stateJ < numState; stateJ++) {
          xi_temp[stateI][stateJ] -= xi_scale_factor;
          xi_sum[stateI][stateJ] = ADD_LOG_PROB(xi_sum[stateI][stateJ], xi_temp[stateI][stateJ]);
        }
      }
    }
  } //n
  //Verbose::out(3,"CNAnalysisMethodCN::fwdback: max delta between gamma_sum is " + ToStr(max_gamma_sum - min_gamma_sum));

  return fwd_loglik;
}

/**
 * @brief
 * @param previous_loglik     - Previous log likelihood.
 * @param loglik              - Log likelihood.
 * @param emConvergenceThresh - The EM convergence threshold.
 * @param bEMconverged        - NOT NEEDED (see todo)
 * @return The status of EM convergence.
 * @todo Remove the "bEMconverged" parameter since it is not used.
 */
bool CNAnalysisMethodCN::EMconverge(
    const double             previous_loglik,
    const double            loglik,
    const double            emConvergenceThresh,
    bool                bEMconverged
    )
{

  bEMconverged = false;

  // eps: is the distance from 1 -> next larger dbl precision #
  // eps is a small number set to 2^-52
  double    eps = 2e-52;

  double delta_loglik = fabs(loglik - previous_loglik);
  double loglik_avg = 0.5*(fabs(loglik) + fabs(previous_loglik) + eps);

  //Verbose::out(3,"CNAnalysisMethodCN::hmm:EMconverge - Testing for EM convergence... ");
  if ( fabs(loglik_avg) !=HUGE_VAL ) {
    bEMconverged = (((delta_loglik/loglik_avg) < emConvergenceThresh));
  }
  //Verbose::out(10,"CNAnalysisMethodCN::hmm:EMconverge - previous_loglik: " + ToStr(previous_loglik) + ", loglik: " + ToStr(loglik)  +  ", emConvergenceThresh " + ToStr(emConvergenceThresh) + ", delta_loglik: " + ToStr(delta_loglik) + ", loglik_avg " + ToStr(loglik_avg) + ", bEMconverged: " + ToStr(bEMconverged));

  return bEMconverged;
}

/**
 * @brief
 * @param samplr                - Sample log2 ratios
 * @param gamma                 -
 * @param gamma_state_sum       -
 * @param xi_sum                -
 * @param log_transmat_stat     -
 * @param updated_log_priorProb - Updated prior probabilities.
 * @param updated_mu            - Updated mu.
 * @param updated_sigma         - Updated sigma.
 * @param compute_obsloglik     -
 * @param bEMconverged          -
 * @param normalState           - The normal copy number state (user specified).
 */
void CNAnalysisMethodCN::reEstimate (
    const VecCNLR&                     samplr,
    const vector< vector<double> >&     gamma,
    const vector<double>&                   gamma_state_sum,
    const vector< vector<double> >&     xi_sum,
    vector< vector<double> >&             log_transmat_stat,
    vector<double>&                           updated_log_priorProb,
    vector<double>&                           updated_mu,
    vector<double>&                           updated_sigma,
    bool&                                   compute_obsloglik,
    bool&                        bEMconverged,
    int                                     normalState
    )
{
  int numState = (int)updated_log_priorProb.size(),
    numobs = (int)samplr.size();
  bool ok_to_update_emat = true;
  vector<double> sigma_update(numState,0), mu_update(numState,0), gamma_update(numState,0);

  //Verbose::out(3,"CNAnalysisMethodCN::hmm:reEstimate - Performing re-estimation... ");
  //Re-estimate and update the prior
  double prior_sum = -FLT_MAX;
  if (m_updatePrior) {
    for(int stateI = 0; stateI < numState; stateI++) {
      updated_log_priorProb[stateI] = gamma[stateI][0];
      prior_sum = ADD_LOG_PROB(prior_sum, updated_log_priorProb[stateI]);
    }
    if(fabs(exp(prior_sum)-1.0)>1e-5) {
     // for(int stateI = 0; stateI < numState; stateI++)
        //Verbose::out(3,"  Prior[" + ToStr(stateI) + "]\t" + ToStr(updated_log_priorProb[stateI]) + "\n");
    }
    CNERR_CHECK((fabs(exp(prior_sum)-1)<1e-5),"CNAnalysisMethodCN::hmm:fwdback - Sum of priors: " + ToStr(exp(prior_sum)) + " does not equal 1");
  }


  // Re-estimate & update the emission parameters
  if (!m_updateMu || !m_updateSigma) {
    compute_obsloglik = false;
  }

  if (m_updateMu || m_updateSigma || m_updateTransMat) {
    for (int stateI = 0; stateI < numState; stateI++) {
      mu_update[stateI] = 0;
      sigma_update[stateI] = 0;
      for ( int n = 0; n < numobs; n++) {
        double weight = exp(gamma[stateI][n] - gamma_state_sum[stateI]);
    mu_update[stateI] += weight*samplr[n];
    sigma_update[stateI] += weight*((samplr[n] - updated_mu[stateI]) * (samplr[n] - updated_mu[stateI]));
      } //numobs
      if(sigma_update[stateI] < MIN_SIGMA) {
        sigma_update[stateI] = MIN_SIGMA;
      }
    } // stateI

    // CHECK:  Check for some possible illegal parameter values before making the update.
    for (int stateI = 0; stateI < numState; stateI++) {
      // The order of mu should always be maintained- sometimes the order can be disrupted on account of negligible datapoints in a given state
      if((stateI > 0) && (mu_update[stateI] < mu_update[stateI-1])) {
        //Verbose::out(3,"CNAnalysisMethodCN::hmm:reEstimate - mu[" + ToStr(stateI) + "] is less than mu[" + ToStr(stateI-1) + "]");
        ok_to_update_emat = false;
      }
      //Check to ensure updated mu for deletion states do not switch to a positive  value
      if((stateI < normalState) && (mu_update[stateI] > updated_mu[normalState])) {
        //Verbose::out(3,"CNAnalysisMethodCN::hmm:reEstimate - mu[" + ToStr(stateI) + "] should be less than " + ToStr(updated_mu[normalState]));
        ok_to_update_emat = false;
      }
      //Check to ensure updated mu for amplification states do not switch to negative value
      else if ((stateI > normalState) && (mu_update[stateI] < updated_mu[normalState])) {
        //Verbose::out(3,"CNAnalysisMethodCN::hmm:reEstimate - mu[" + ToStr(stateI) + "] should be greater than " + ToStr(updated_mu[normalState]));
        ok_to_update_emat = false;
      }
     } //stateI

    compute_obsloglik = ok_to_update_emat;

    //Verbose::out(10,"CNAnalysisMethodCN::hmm:reEstimate - Update Emission Matrix? " + ToStr(ok_to_update_emat));
    //Verbose::out(10,"CNAnalysisMethodCN::hmm:reEstimate - Comput_obslik? " + ToStr(compute_obsloglik));

    if (ok_to_update_emat){
      for (int stateI = 0; stateI < numState; stateI++) {
        if(m_updateMu) {
      updated_mu[stateI] = mu_update[stateI];
      //Verbose::out(10,"CNAnalysisMethodCN::hmm:reEstimate - Updated Mu: " + ToStr(updated_mu[stateI]));
        }
        if(m_updateSigma) {
      updated_sigma[stateI] = sqrt(sigma_update[stateI]);
      //Verbose::out(10,"CNAnalysisMethodCN::hmm:reEstimate - Updated Sigma: " + ToStr(updated_sigma[stateI]));
        }
      } //stateI
    }
  } // check update status


  //NOTE: The non-stationary form gets updated via the prior
  //Re-estimate the transition matrix - the stationary form
  if ( m_updateTransMat & ok_to_update_emat ) {
    for (int stateI = 0; stateI < numState; stateI++) {
      // Need to have a gamma sum over all but last element for transmat update
      double denominator = SUBTRACT_LOG_PROB(gamma_update[stateI], gamma[stateI][numobs-1]);
      for (int stateJ = 0; stateJ < numState; stateJ++) {
    log_transmat_stat[stateI][stateJ] = xi_sum[stateI][stateJ] - denominator;
           //Verbose::out(10,"CNAnalysisMethodCN::hmm:reEstimate - Updated Stat Transmat: " + ToStr(log_transmat_stat[stateJ][stateI]) + "state : " + ToStr(stateI) + ToStr(stateJ));
      } //stateJ
    } //stateI
  } //update_log_transmat_stat
}

/**
 * @brief
 * @param log_lik_obs           -
 * @param dist_covariate        -
 * @param log_priorProb         -
 * @param updated_log_priorProb -
 * @param transMatDecay         -
 * @param trans_type_stat       -
 * @param log_transmat_stat     -
 * @param state                 - in/out: CN states.
 * @return The logP
 */
double CNAnalysisMethodCN::viterbi (
        // Inputs
    const vector< vector<double> >& log_lik_obs,
    const vector<double>&         dist_covariate,
    const vector<double>&           log_priorProb,
    const vector<double>&            updated_log_priorProb,
    const double                    transMatDecay,
    const bool                      trans_type_stat,
    const vector< vector<double> >& log_transmat_stat,
    // Outputs
    vector<int>&                    state
    )
{
  double max_val,
  trans_elem,
  prior_mult_exp_fn;
  int max_id = 0,
  numobs  = (int)log_lik_obs[0].size(),
  numState = (int)log_lik_obs.size();
  vector< vector<double> > delta(numState, vector<double>(numobs));
  vector< vector<int> >    psi(numState, vector<int>(numobs));

  //Verbose::out(3,"CNAnalysisMethodCN::hmm:viterbi - Computing Viterbi path... ");
  CNERR_CHECK((int)state.size()==numobs, "CNAnalysisMethodCN::viterbi - The number of states does NOT equal the number of observations.");

  // Step 1:Initialization
  // Note: log_lik_obs ~ log(emission probability matrix)
  for(int stateI = 0; stateI < numState; stateI++) {
    delta[stateI][0] = (m_updateTransMat) ? (log_lik_obs[stateI][0] + updated_log_priorProb[stateI]) : (log_lik_obs[stateI][0] + log_priorProb[stateI]);
    psi[stateI][0] = max_id;
  } //stateI

  // Step 2: Recursion
  double log_decay_factor, decay_factor;
  for(int n = 1; n < numobs; n++) {
    log_decay_factor = -(1/transMatDecay)*dist_covariate[n-1];
    decay_factor = exp(log_decay_factor);
    for(int stateI = 0; stateI < numState; stateI++) {
      prior_mult_exp_fn = (m_updateTransMat) ? (updated_log_priorProb[stateI]+log(1-decay_factor)) : (log_priorProb[stateI]+log(1-decay_factor));
      max_id = 0;
      // Choose a very small # -Infinity to start w/
      max_val = -1 * HUGE_VAL;
      for(int stateJ = 0; stateJ < numState; stateJ++) {
        if (trans_type_stat) {
          trans_elem = log_transmat_stat[stateJ][stateI];
        } else {
          trans_elem = prior_mult_exp_fn;
          if(stateJ == stateI) {
            trans_elem = ADD_LOG_PROB(trans_elem, log_decay_factor);
          }
        }

        double dTMP = delta[stateJ][n-1] + trans_elem;
        if ( dTMP > max_val ) {
          max_val = dTMP;
          max_id = stateJ;
        }
      } //stateJ
      // Save the transition information for later backtracking
      psi[stateI][n] = max_id;
      // Update delta
      delta[stateI][n] = log_lik_obs[stateI][n] + max_val;
    } //stateI
  }// numobs

  // Step 3: Termination: Determine the final step(stateF) w/ the max probability(logP)
  double logP = delta[0][numobs-1];
  int stateF = 0;
  for(int stateI = 0; stateI < numState; stateI++) {
    if ( delta[stateI][numobs-1] > logP ) {
      logP = delta[stateI][numobs-1];
      stateF = stateI;
    } // delta[stateI][numobs-1] > logP
  } //stateI

  // Step 4: Path(state sequence) backtracking
  state[numobs-1] = stateF;
  //Verbose::out(10,"CNAnalysisMethodCN::hmm:viterbi - Observation: " + ToStr(numobs-1) + " state: " + ToStr(state[numobs-1]));
  for(int n = numobs-2; n >=0; n--) {
    state[n] = psi[state[n+1]][n+1];
    // Verbose::out(10,"CNAnalysisMethodCN::hmm:viterbi - Observation: " + ToStr(n) + " state: " + ToStr(state[n]));
    if (state[n]==0) {
      // Verbose::out(10,"CNAnalysisMethodCN::hmm:viterbi - Observation: " + ToStr(n) + " state: " + ToStr(state[n]));
    }
  } //numobs

   return logP;
}

/**
 * @brief Performs the CopyNumber algorithm on log2 ratios of a sample on
 * a per chromosome basis.
 *
 * Note:  Currently only states will be used after calling this but in the
 *        future this may change.
 * @param position            - SNP positions.
 * @param samplr              - Sample log2 ratios.
 * @param state               - in/out: CN states.
 * @param priorProb           - in/out: Prior probabilities.
 * @param mu                  - in/out: Mu values.
 * @param sigma               - in/out: Sigma values.
 * @param transMatDecay       - Transition decay matrix.
 * @param stateEstMethod      - State estimation method.
 * @param stateMargProb       - State marginal probability (may be modified in future).
 * @param emIterations        - Number of EM iterations.
 * @param emConvergenceThresh - EM convergence threshold.
 * @param b_allelespecific    -
 * @param normalState         - The normal copy number state (user specified).
 */
void
CNAnalysisMethodCN::hmm(
           const vector<int>&  position,
           const VecCNLR&               samplr,
           vector<int>&                 state,
           const vector<double>&        priorProb,
           const vector<double>&        mu,
           const vector<double>&        sigma,
           double                       transMatDecay,
           const string&                stateEstMethod,
           vector<double>&              stateMargProb,
           const unsigned int           emIterations,
           const double                 emConvergenceThresh,
           bool                         b_allelespecific,
           int                          normalState
           ) {


  // Error Handling
  CNERR_CHECK(position.size()==samplr.size(),
         "CNAnalysisMethodCN::hmm - the number of positions " + ToStr(position.size()) + " does not equal the number of SNPs being evaluated " + ToStr(samplr.size()));
  CNERR_CHECK(position.size()==state.size(),
         "CNAnalysisMethodCN::hmm - the number of positions " + ToStr(position.size()) + " does not equal the number of SNPs whose states are being evaluated " + ToStr(state.size()));
  CNERR_CHECK(priorProb.size()==mu.size(),
         "CNAnalysisMethodCN::hmm - the number of states determined from Priors " + ToStr(priorProb.size()) + " and mu " + ToStr(mu.size()) + " do not match");
  CNERR_CHECK(priorProb.size()==sigma.size(),
         "CNAnalysisMethodCN::hmm - the number of states determined from Priors " + ToStr(priorProb.size()) + " and variance " + ToStr(sigma.size()) + " do not match");
  CNERR_CHECK(mu.size()==sigma.size(),
         "CNAnalysisMethodCN::hmm - the number of states in the emission matrix components(mu and sigma) do not match; #States(mu): " + ToStr(mu.size()) + " #States(variance)" + ToStr(sigma.size()));
  CNERR_CHECK(emIterations>=1,
         "CNAnalysisMethodCN::hmm - Number of EM Iterations must be at least 1");

  // Need at least 2datapoints: corresponding to current and prior markers
  if (position.size() == 0 ) {
    return;
  }

  // Initialization
  int numState = (int)priorProb.size(),
  numobs = (int)samplr.size();
  double previous_loglik = -1*HUGE_VAL,
  loglik = 0.0,
  logP = 0.0;
  bool bEMconverged = false,
  trans_type_stat = false,
  compute_obsloglik=true;
  vector< vector<double> > log_lik_obs(numState, vector<double>(numobs)),
  alpha(numState, vector<double>(numobs)),
  beta(numState, vector<double>(numobs)),
  gamma(numState, vector<double>(numobs)),
  xi_sum(numState, vector<double>(numState,0)),
  log_transmat_stat(numState, vector<double>(numState,-FLT_MAX));
  vector<double> dist_covariate(numobs-1,0),
  log_priorProb(numState,-FLT_MAX),
  updated_log_priorProb(numState,-FLT_MAX),
  updated_mu(numState,0),
  updated_sigma(numState,0),
  gamma_state_sum(numState,0);

  if ( trans_type_stat ) {
    CNERR_CHECK((int)m_diagProb.size()==numState,
       "CNAnalysisMethodCN::hmm - number of elements in the Diagonal Probability Matrix " + ToStr(m_diagProb.size()) + " does not match number of states in the model " + ToStr(numState));
  }

  //Verbose::out(10,"CNAnalysisMethodCN::hmm:Number of EM iterations: " + ToStr(emIterations));

  // Do CN HMM algorithm.
  // 0: Initialize the HMM parameters
  for (int stateI = 0; stateI < numState; stateI++) {
    log_priorProb[stateI] = (priorProb[stateI] > 0) ? log(priorProb[stateI]) : -FLT_MAX;
    updated_log_priorProb[stateI] = log_priorProb[stateI];
    // For allele specific case
    if (b_allelespecific) {
      updated_mu[stateI] = mu[stateI] + 1 ;
    }
    else {
      updated_mu[stateI] = mu[stateI];
    }
    updated_sigma[stateI] = sigma[stateI];
    // if stationary transition matrix is preferred:
    for (int stateJ = 0; stateJ < numState; stateJ++) {
      if(stateJ==stateI) {
        log_transmat_stat[stateI][stateJ] = (m_diagProb[stateJ] > 0.0) ? log(m_diagProb[stateJ]) : -FLT_MAX;
      } else {
        log_transmat_stat[stateI][stateJ] = log(1.0-m_diagProb[stateJ]) - log(numState-1.0);
      }
    } //stateJ
  } //stateI

  // 1: Generate the inter-SNP distance - covariate for non-stationary transition matrix
  for(int n = 0; n < (numobs - 1); n++) {
    dist_covariate[n] = position[n+1] - position[n];
  }

  // EM step
  for (int EMiter = 0; EMiter < (int)emIterations; EMiter++) {
    //Verbose::out(3,"Performing EM iteration... " + ToStr(EMiter));

    // E Step: 2-3:
    // 2: Generate the observed log likelihood(emission probabilities)
    if (compute_obsloglik) {
      obsloglik(samplr,numState,updated_mu,updated_sigma,log_lik_obs);
    }

    // 3: The fwdback computation
    loglik = fwdback(dist_covariate,log_priorProb,updated_log_priorProb,log_lik_obs,transMatDecay,trans_type_stat,log_transmat_stat,alpha,beta,gamma,gamma_state_sum,xi_sum,previous_loglik);
    //Verbose::out(3,"CNAnalysisMethodCN::hmm:Updated Log Likelihood post fwdback: " + ToStr(loglik));

    if ((EMiter == 0)) {
      bEMconverged = (numobs == 1) ? true : false;
    } else {
      bEMconverged = (numobs == 1) ? true : EMconverge(previous_loglik,loglik,emConvergenceThresh,bEMconverged);
      if(loglik < previous_loglik) {
        // This shouldn't actually happen, since EM should guarantee a non-decreasing likelihood
        //Verbose::out(3,"CNAnalysisMethodCN::hmm: Terminating EM, log likelihood decreased from " + ToStr(previous_loglik) + " to " + ToStr(loglik));
        bEMconverged = true;
      }
    }
    previous_loglik = loglik;

    // 5: M step: reEstimation & updates
    if (!bEMconverged) {
      reEstimate(samplr,gamma,gamma_state_sum,xi_sum,log_transmat_stat,updated_log_priorProb,updated_mu,updated_sigma,compute_obsloglik,bEMconverged,normalState);
    } else {
      break;
    }
    //Verbose::out(3,"  CNAnalysisMethodCN::hmm: End of EM iteration " + ToStr(EMiter) + ":");
    for (int stateI = 0; stateI < numState; stateI++) {
      //Verbose::out(3,"    CNAnalysisMethodCN::hmm:  priorProb[" + ToStr(stateI) + "] is " + ToStr(exp(updated_log_priorProb[stateI])));
      //Verbose::out(3,"    CNAnalysisMethodCN::hmm:  mu[" + ToStr(stateI) + "] is " + ToStr(updated_mu[stateI]));
      //Verbose::out(3,"    CNAnalysisMethodCN::hmm:  sigma[" + ToStr(stateI) + "] is " + ToStr(updated_sigma[stateI]));
    }
  }//EMiter

  for (int stateI = 0; stateI < numState; stateI++) {
    //Verbose::out(3, "CNAnalysisMethodCN::state:  " + ToStr(stateI) + " - Updated prior probability: " + ToStr(exp(updated_log_priorProb[stateI])) + " - Updated mu: " + ToStr(updated_mu[stateI]) + " - Updated sigma: " + ToStr(updated_sigma[stateI]));
  }

  // 6: Selection of optimal path sequence via Viterbi
  logP = viterbi(log_lik_obs,dist_covariate,log_priorProb,updated_log_priorProb,transMatDecay,trans_type_stat,log_transmat_stat,state);

  //Verbose::out(3,"Done performing CopyNumber HMM.");
} // CNAnalysisMethodCN::hmm(...)



/**
 * @brief Smoothing outliers is done on a per chromosome basis for a given sample.
 * @param state      - in/out: CN states.
 * @param pos        - Position of SNP: used to compute distances.
 * @param distThresh - Distance threshold.
 */
void
CNAnalysisMethodCN::smoothOutliers(
                  vector<int>&          state,
                  const vector<double>& pos,
                  double                distThresh
                  ) {
  //Verbose::out(3,"CopyNumber HMM: Smoothing outliers...");
  CNERR_CHECK(state.size()>0,"CNAnalysisMethodCN::smoothOutliers - no states were provided");
  CNERR_CHECK(state.size()==pos.size(), "CNAnalysisMethodCN::smoothOutliers - The number of chromosome states (" + ToStr(state.size()) + ") and positions (" + ToStr(pos.size()) + ") don't match");

  unsigned int n = (unsigned int)state.size();
  if(n > 1) {
    vector<bool> change(n+1);  // change[i] is true if state[i] is different to state[i-1]
    vector<double> dist(n+1);  // dist[i] is pos[i]-pos[i-1]
    change[0] = change[n] = true;
    dist[0] = dist[n] = pos[n-1]-pos[0]+1+distThresh;

    // Pass through once to locate all state transitions
    for(unsigned int i=1; i<n; i++) {
      dist[i] = pos[i]-pos[i-1];
      if(state[i] != state[i-1]) {
        change[i] = true;
      } else {
        change[i] = false;
      }
    }

    // Pass through a second time to identify states of duration 1 within distThresh of adjacent states.
    // Any such states found are assigned to the nearest surrounding state.
    vector<int> newstate(n);
    for(unsigned int i=0; i<n; i++) {
      if(change[i] && change[i+1] && ((dist[i] <= distThresh) || (dist[i+1] <= distThresh))) {
        if(dist[i] < dist[i+1]) {
          newstate[i] = state[i-1];
        } else {
          newstate[i] = state[i+1];
        }
      } else {
        newstate[i] = state[i];
      }
    }

    // Copy final results
    for(unsigned int i=0; i<n; i++) {
      state[i] = newstate[i];
    }
  }

  //Verbose::out(3,"CopyNumber HMM: Done Smoothing outliers.");
}

/**
 * @brief This is done across all chromosomes for a given sample.
 * The indexes for the "vector< vector<TYPE> >" items are
 * [chromosome][state/log2ratio/medianCN/pValue]
 * Example:  For chromosome index 2 the 5th log2 ratio values would
 * be "samplr[2][4]"
 * @param state                      - in/out: CN states.
 * @param logRatio                   - Log2 ratios for all SNPs and all chromosomes.
 * @param normalState                - The normal copy number state (user specified).
 * @param normalStateMinObservations - Minimum number of observations for the normal state.
 * @param estCorrectionFactor        -
 * @param correctionFactor           -
 * @param medianCN                   - Median copy number values for all chromosomes.
 */
void
CNAnalysisMethodCN::postProcessHMM(
                  const vector< vector<int> >&    state,
                  const VecVecCNLR&               logRatio,
                  int                             normalState,
                      int                             normalStateMinObservations,
                  bool                            estCorrectionFactor,
                  double&                         correctionFactor,
                  VecVecCNLR&                     medianCN
                  ) {
  //Verbose::out(3,"CopyNumber HMM: Post processing HMM...");
  CNERR_CHECK(state.size()==logRatio.size(), "CNAnalysisMethodCN::postProcessHMM - The number of states and log2 ratios don't match");
  if(state.size() == 0) {
    return;
  }


  //
  // Loop over positions in chromosomes to determine segment bounds and
  // collect summary statistics per segment from which we can compute
  // t-stats etc.
  //
  unsigned int nChrom = (unsigned int)state.size();
  //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm:# of chromosome processed (nChrom): " + ToStr(nChrom));
  //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm:What is the normal state? " + ToStr(normalState));
  //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm:normalStateMinObservations: " + ToStr(normalStateMinObservations));
  //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm:CorrectionFactor: " + ToStr(correctionFactor));

  vector< vector<int> >    segment_state(nChrom);            // estimated state for each segment
  vector< vector<int> >    segment_end_index(nChrom);        // index of first position after end of segment
  VecVecCNLR segment_average_x(nChrom);        // the average log ratio for each segment
  VecVecCNLR segment_average_xsquared(nChrom); // the average squared log ratio for each segment
  VecVecCNLR segment_median(nChrom);           // the median log ratio for each segment

  // The follwoing normal_state variables hold stats relating to
  // observations of the normal state across all chromosomes.
  VecCNLR normal_state_average(N_NORMAL_SEGMENTS);    // the averages for each normal segment
  VecCNLR normal_state_sd(N_NORMAL_SEGMENTS);         // the sds for each normal segment
  vector<U_INT> normal_state_nObs(N_NORMAL_SEGMENTS); // the number of SNPs in each normal segment
  U_INT normal_state_n=0;                             // the number of normal segments

  for(unsigned int chromIndex=0; chromIndex<nChrom; chromIndex++) {
    CNERR_CHECK( state[chromIndex].size()==logRatio[chromIndex].size(),
      "CNAnalysisMethodCN::postProcessHMM - The number of snp states and log2 ratios don't match for seq " + ToStr(chromIndex)
    );
    //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm:Processing chromsome... " + ToStr(chromIndex));
    collect_segment_stats(
      state[chromIndex],
      logRatio[chromIndex],
      normalState,
      normalStateMinObservations,
      segment_state[chromIndex],
      segment_end_index[chromIndex],
      segment_average_x[chromIndex],
      segment_average_xsquared[chromIndex],
      segment_median[chromIndex],
      normal_state_average,
      normal_state_sd,
      normal_state_nObs,
      normal_state_n
    );
  }

  //
  // Process collection of all normal copy number state states to
  // get a weighted average of the normal copy number state mean and sd,
  // which will be used in two-sample t-test measuring significance
  //
  unsigned int nTotal=0;
  for(unsigned int i=0; i<normal_state_n; i++) {
    nTotal += normal_state_nObs[i];
  }
  //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm: # of normal_state segments "    + ToStr(normal_state_n) + " with nTotal SNPs " + ToStr(nTotal));
  double baseline_average = 0;
  double baseline_sd = 0;
  for(unsigned int i=0; i<normal_state_n; i++) {
    double weight = (normal_state_nObs[i]/(double)nTotal);
    baseline_average +=  weight * normal_state_average[i];
    baseline_sd +=  weight * normal_state_sd[i];
  }
  //The baseline average & sd reflects the weighted average; where the weighting is based on the #of SNPs per normal segment
  //Verbose::out(3,
  //  "copyNumber HMM: post processing HMM: normal state (average,sd) are (" + ToStr(baseline_average) +
  //  "," + ToStr(baseline_sd) + ") based on " + ToStr(nTotal) + " SNPs in " + ToStr(normal_state_n) + " segments"
  //);

  //
  // Loop over all chomosomes and all segments to write out t-stats and median log ratio
  //
  medianCN.resize(nChrom);
//  double term1 = (nTotal > 0) ? baseline_sd * baseline_sd / (double) nTotal : 0;
  for(unsigned int chromIndex=0; chromIndex<nChrom; chromIndex++) {
    unsigned int nState = (unsigned int)state[chromIndex].size();
    medianCN[chromIndex].resize(nState);
    unsigned int segIndex,stateIndex;
//  int prev_segment_end_index=0;
    for(segIndex=0,stateIndex=0; segIndex< segment_state[chromIndex].size(); segIndex++) {
 //     double numerator = segment_average_x[chromIndex][segIndex] - baseline_average;
//      double sd2 = segment_average_xsquared[chromIndex][segIndex] - segment_average_x[chromIndex][segIndex] * segment_average_x[chromIndex][segIndex];
//      int n2 = segment_end_index[chromIndex][segIndex] - prev_segment_end_index;
 //     double term2 = sd2 * sd2 / (double) (n2);
//      double denominator = sqrt(term1 + term2);
      //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm: # chromindex: " + ToStr(chromIndex) + " segIndex " + ToStr(segIndex) + " stateIndex " + ToStr(stateIndex));
      //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm: # term1: " + ToStr(term1));
      //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm: # term2: " + ToStr(term2) + " sd2 " + ToStr(sd2) + " n2: " + ToStr(n2));
      //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm: # Numerator: " + ToStr(numerator) + " Denominator: " + ToStr(denominator));
      double med = segment_median[chromIndex][segIndex];

      for( ; (int)stateIndex < segment_end_index[chromIndex][segIndex]; stateIndex++) {
        medianCN[chromIndex][stateIndex] = (float)med;
      }
    }
  }

  //Verbose::out(3,"copyNumber HMM: Done post processing HMM.");
}

/**
 * @brief Collects segment stats.
 * @param state                      - in/out: CN states.
 * @param logRatio                   - Log2 ratios for all SNPs and all chromosomes.
 * @param normalState                - The normal copy number state (user specified).
 * @param normalStateMinObservations - Minimum number of observations for the normal state.
 * @param segment_state              -
 * @param segment_end_index          -
 * @param segment_average_x          -
 * @param segment_average_xsquared   -
 * @param segment_median             -
 * @param normal_state_average       - Normal state log2 ratio average.
 * @param normal_state_sd            - Normal state log2 ratio standard deviation.
 * @param normal_state_nObs          - Normal state number of observations.
 * @param normal_state_n             -
 */
void CNAnalysisMethodCN::collect_segment_stats(
          const vector<int>&     state,
          const VecCNLR&         logRatio,
          int                    normalState,
          int                    normalStateMinObservations,
          vector<int>&           segment_state,
          vector<int>&           segment_end_index,
          VecCNLR&               segment_average_x,
          VecCNLR&               segment_average_xsquared,
          VecCNLR&               segment_median,
          VecCNLR&               normal_state_average,
          VecCNLR&               normal_state_sd,
          vector<U_INT>&         normal_state_nObs,
          U_INT&                 normal_state_n
) {
  unsigned int nState = (unsigned int)state.size();
  segment_state.resize(nState);
  segment_end_index.resize(nState);
  segment_average_x.resize(nState);
  segment_average_xsquared.resize(nState);
  segment_median.resize(nState);

  if(nState==0) {
    return;
  }

  VecCNLR segment_logRatio;       // Buffer for holding logRatios within segment to allow taking medians
  segment_logRatio.resize(nState);
  unsigned int j;
  unsigned int n=0;
  int k = 0, prevStateIndex = 0;
  for(j=0,k=-1; j<nState; j++) {
    if(j==0 || state[j] != segment_state[k]) {
      // We're either at the first state, or we just finished a segment
      ++k;
      segment_state[k] = state[j];
      if(k > 0) {  // just finished a segment
        prevStateIndex = k-1;
        segment_end_index[prevStateIndex] = j;
        segment_median[prevStateIndex] = median_in_place(segment_logRatio.begin(),segment_logRatio.begin()+n);
        // if the finished segment is a normal state with sufficient observations, record its stats
        if((segment_state[prevStateIndex] == normalState) && ((int)n >= normalStateMinObservations)) {
          // Check if we need to expand the normal state vectors
          if(normal_state_average.size() == normal_state_n) {
            unsigned int new_size = (unsigned int) ceil(normal_state_average.size() * N_NORMAL_SEGMENTS_GROWTH_FACTOR);
            normal_state_average.resize(new_size);
            normal_state_sd.resize(new_size);
            normal_state_nObs.resize(new_size);
          }
          // Write in new normal state stats
      //Verbose::out(10,"CNAnalysisMethodCN::postProcesshmm:collect_segment_stats:# of normal states in a segment: " + ToStr(n));
          normal_state_average[normal_state_n] = segment_average_x[prevStateIndex];
          CNLRT sd = (CNLRT)((segment_average_xsquared[prevStateIndex] - segment_average_x[prevStateIndex] * segment_average_x[prevStateIndex]) * n / (double) (n-1));
          normal_state_sd[normal_state_n] = sd;
          normal_state_nObs[normal_state_n] = n;
          normal_state_n++;
        }
      }
      segment_average_x[k] = logRatio[j];
      segment_average_xsquared[k] = logRatio[j] * logRatio[j];
      segment_logRatio[0] = logRatio[j];
      n=1;
    } else {
      segment_logRatio[n] = logRatio[j];
      ++n;
      segment_average_x[k] = (float)(((n-1)*segment_average_x[k]/(double)n) + (logRatio[j]/(double)n));
      segment_average_xsquared[k] = (float)(((n-1)*segment_average_xsquared[k]/(double)n) + ((logRatio[j] * logRatio[j])/(double)n));
    }
  }
  // Will have just finished the final segment
  prevStateIndex = k;
  if((segment_state[prevStateIndex] == normalState) && ((int)n >= normalStateMinObservations)) {
    // Check if we need to expand the normal state vectors
    if(normal_state_average.size() == normal_state_n) {
      unsigned int new_size = (unsigned int) ceil(normal_state_average.size() * N_NORMAL_SEGMENTS_GROWTH_FACTOR);
      normal_state_average.resize(new_size);
      normal_state_sd.resize(new_size);
      normal_state_nObs.resize(new_size);
    }
    // Write in new normal state stats
    normal_state_average[normal_state_n] = segment_average_x[prevStateIndex];
    CNLRT sd = (CNLRT)((segment_average_xsquared[prevStateIndex] - segment_average_x[prevStateIndex] * segment_average_x[prevStateIndex]) * n / (double) (n-1));
    normal_state_sd[normal_state_n] = sd;
    normal_state_nObs[normal_state_n] = n;
    normal_state_n++;
  }
  segment_end_index[prevStateIndex] = j;
  segment_median[prevStateIndex] = median_in_place(segment_logRatio.begin(),segment_logRatio.begin()+n);

  ++k;
  segment_state.resize(k);
  segment_end_index.resize(k);
  segment_average_x.resize(k);
  segment_average_xsquared.resize(k);
  segment_median.resize(k);

  return;
}

int CNAnalysisMethodCN::getLastAutosomeChromosome()
{
    int iLastAutosomeChromosome = 1;
    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        if (pobjProbeSet->getChromosome() < m_iXChromosome)
        {
            iLastAutosomeChromosome = Max(iLastAutosomeChromosome, (int)pobjProbeSet->getChromosome());
        }
    }
    return iLastAutosomeChromosome;
}
