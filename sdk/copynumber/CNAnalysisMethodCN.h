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

#ifndef _CNAnalysisMethodCN_H_
#define _CNAnalysisMethodCN_H_
/**
 * @file CNAnalysisMethodCN.h
 *
 * @brief This header contains the CNAnalysisMethodCN class definition.
 */

#include "copynumber/CNAnalysisMethod.h"

#define CNAnalysisMethodCN_Version "1.0"

#define CNERR_CHECK(x, y) if ( !(x) ) throw(Except(y))
typedef unsigned int U_INT; ///< Unsigned int.

///////////////////////////////////////////////////////////////////////////
// Basic types.
// NOTE: CNLRT should be used since CNLR is deprecated.
///////////////////////////////////////////////////////////////////////////
typedef float                               CNLRT;    ///< CN log2 ratio.
typedef CNLRT                               CNLR;     ///< CN log2 ratio.
typedef long                                FragLenT; ///< Fragment length.
typedef int                                    PosT;     ///< Position.
typedef double                              GCCT;     ///< GC content.
typedef int                                 HmmStateT;///< State.
typedef double                              PValT;    ///< P-Value.

/// Genotype Type
enum geno_type_t { AA, AB, BB, N};

///////////////////////////////////////////////////////////////////////////
// Collections: this is for convenience so one does not have to write
// out something like "std::map<std::string,std::string>"
///////////////////////////////////////////////////////////////////////////
typedef std::map<std::string,std::string>   StrStrMap; ///< Map of string to string.
typedef std::map<std::string,int>           StrIntMap; ///< Map of String to int.
typedef std::map<int,std::string>           IntStrMap; ///< Map of int to string.
typedef std::map<int,int>                   IntIntMap; ///< Map of int to int.
typedef std::set<int>                       IntSet;    ///< Set of ints.
typedef std::set<std::string>               StrSet;    ///< Set of Strings.
typedef std::vector<int>                    VecInt;    ///< Vector of ints.
typedef std::vector<double>                 VecD;      ///< Vector of doubles.
typedef std::vector<float>                  VecF;      ///< Vector of floats.
typedef std::vector< std::vector<int> >     VecVecInt; ///< Vector of vector of ints.
typedef std::vector< std::vector<int> >     VecVecI;   ///< Vector of vector of ints.
typedef std::vector< VecD >                 VecVecD;   ///< Vector of vector of doubles.
typedef std::vector< VecF >                 VecVecF;   ///< Vector of vector of floats.
typedef std::vector<std::string>            VecStr;    ///< Vector of strings.
typedef VecStr                              StrVec;    ///< Vector of strings.
typedef std::vector<VecStr>                 VecVecStr; ///< Vector of vector of strings.

typedef std::vector<CNLRT>                  VecCNLR;     ///< Vector of CNLRT.
typedef std::vector<HmmStateT>              VecHmmState; ///< Vector of HmmStateT.
typedef std::vector<PValT>                  VecProb;     ///< Vector of PValT.
typedef std::vector<PosT>                   VecPosition; ///< Vector of PosT.
typedef std::vector<FragLenT>               VecFraglen;  ///< Vector of FragLenT.
typedef std::vector<GCCT>                   VecGCC;      ///< Vector of GCCT.

typedef std::vector< VecCNLR >              VecVecCNLR;     ///< Vector of vector of CNLR.
typedef std::vector< VecHmmState >          VecVecHmmState; ///< Vector of vector of HMMStateT.
typedef std::vector< VecProb >              VecVecProb;     ///< Vector of vector of PValT.
typedef std::vector< VecPosition >          VecVecPosition; ///< Vector of vector of PosT.

typedef std::map<std::string,CNLRT>         MapStrCNLR;     ///< Map of string to CNLRT.

#define MIN_PVALUE 1e-16                    ///< The smallest p-value that we can report (since we'll be taking logs of p-values)
#define N_NORMAL_SEGMENTS 1000              ///< The initial buffer for storing normal state segment statistics
#define N_NORMAL_SEGMENTS_GROWTH_FACTOR 1.1 ///< Factor by which the normal state segment vectors will be grown when they overflow
#define MIN_SIGMA 1e-8                      ///< The smallest sigma that we allow for HMM state emissions
#define MAX_SIGMA_RATIO 3                   ///< The biggest allowable ratio between max & min sigma for HMM state emissions

/**
 * @brief Macro to return log(exp(x)+exp(y)) in a manner that avoids
 * underflow. Useful for HMM operations where we need to compute the log
 * of the sum of two probabilities which are stored in log scale and which
 * might be very small.
 */
#define ADD_LOG_PROB(X,Y) ( (Y) > (X) ) ? ((Y) + log(1 + exp((X)-(Y)))) : ((X) + log(1 + exp((Y)-(X))))
/**
 * @brief Macro to return log(exp(x)-exp(y)) in a manner that avoids
 * underflow. Useful for HMM operations where we need to compute the log
 * of the difference of two probabilities which are stored in log scale
 * and which might be very small.
 */
#define SUBTRACT_LOG_PROB(X,Y) ( (Y) > (X) ) ? ((Y) + log(exp((X)-(Y)) - 1)) : ((X) + log(1 - exp((Y)-(X))))

/**
 * @brief  The CN analysis method.
 *
 */
class CNAnalysisMethodCN : public CNAnalysisMethod
{
private:
    // Copy Number State Parameters
    std::vector<int> m_vCNState;
    std::vector<double> m_vPriorProb;
    std::vector<double> m_vMu;
    std::vector<double> m_vSigma;
    double m_dTransMatDecay;
    AffxString m_strStateEstMethod;
    unsigned int m_uiEMIterations;
    double m_dEMConvergentThreshold;
    int m_iNormalState;
    int m_iPostCNFitMaxOutlierRemoveRunSize;
    // Unused parameters.
  int   m_fwdOnly;
  int m_inormalStateMinObservations;
  int m_iSmoothOutliers;
  int   m_transTypeStat;


private:
  bool   m_updatePrior;
  bool   m_updateMu;
  bool   m_updateSigma;
  bool   m_updateTransMat;
  bool   m_scaled;

  std::vector<double> m_diagProb;

public:
    CNAnalysisMethodCN();
    virtual ~CNAnalysisMethodCN();

    static std::string getType();
    static std::string getDescription();
    static std::string getVersion();

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    virtual AffxString getName();

    virtual void setEngine(BaseEngine* p);

    virtual void run();
    virtual bool isSegmentTypeAnalysis();

protected:
    int getLastAutosomeChromosome();
    void copyNumberStatePostProcessing(int iChromosome, std::vector<float>& vLog2Ratios, std::vector<int>& vCNStates);

  /// @brief Compute the observed log-likelihoods 'aka' the emission probabilities.
  void obsloglik(
      const VecCNLR&                       samplr,
      const int                            numState,
      const std::vector<double>&           updated_mu,
      const std::vector<double>&           updated_sigma,
      std::vector< std::vector<double> >&  log_lik_obs
  );
  ///
  double fwdback(
      const std::vector<double>&            dist_covariate,
      const std::vector<double>&            log_priorProb,
      const std::vector<double>&            updated_log_priorProb,
      const std::vector< std::vector<double> >&   log_lik_obs,
      const double                          transMatDecay,
      const bool                            trans_type_stat,
      const std::vector< std::vector<double> >&  log_transmat_stat,
      std::vector< std::vector<double> >&   alpha,
      std::vector< std::vector<double> >&   beta,
      std::vector< std::vector<double> >&   gamma,
      std::vector<double>&                  gamma_state_sum,
      std::vector< std::vector<double> >&   xi_sum,
      const double   previous_loglik
  );
  ///
  bool EMconverge(
      const double   previous_loglik,
      const double   loglik,
      const double   emConvergenceThresh,
      bool           bEMconverged
  );
  ///
  void reEstimate (
      const VecCNLR&                                  samplr, // Sample log2 ratios
      const std::vector< std::vector<double> >&       gamma,
      const std::vector<double>&                      gamma_state_sum,
      const std::vector< std::vector<double> >&       xi_sum,
      std::vector< std::vector<double> >&             log_transmat_stat,
      std::vector<double>&                            updated_log_priorProb,
      std::vector<double>&                            updated_mu,
      std::vector<double>&                            updated_sigma,
      bool&                                           compute_obsloglik,
      bool&                                           bEMconverged,
      int                                             normalState
  );
  ///
  double viterbi (
      const std::vector< std::vector<double> >&       log_lik_obs,
      const std::vector<double>&                      dist_covariate,
      const std::vector<double>&                      log_priorProb,
      const std::vector<double>&                      updated_log_priorProb,
      const double                                    transMatDecay,
      const bool                                      trans_type_stat,
      const std::vector< std::vector<double> >&            log_transmat_stat,
      std::vector<int>&                               state
  );
  ///
  void collect_segment_stats(
      // Inputs
      const std::vector<int>&     state,                      // estimated HMM state sequence
      const VecCNLR&         logRatio,                          // observed log ratios
      int                    normalState,                // definition of which state corresponds to normal copy number
      int                    normalStateMinObservations, // min #observations in normal copy number segment to include in baseline for significance testing
      // Segment outputs - these vectors will be resized according to # segments found
      std::vector<int>&           segment_state,              // state for each segment
      std::vector<int>&           segment_end_index,          // index of state just after segment ends
      VecCNLR&               segment_average_x,          // average log ratio within each segment
      VecCNLR&               segment_average_xsquared,   // average squared log ratio within each segment
      VecCNLR&               segment_median,             // median log ratio within each segment
      // Normal segment outputs - these vectors will be written to starting at index normal_state_n and will be expanded in chunks as necessary
      VecCNLR&               normal_state_average,       // average log ratio per normal state segment
      VecCNLR&               normal_state_sd,            // sample sd of log ratio per normal state segment
      std::vector<U_INT>&         normal_state_nObs,          // number of observations of normal state per normal state segment
      U_INT&                 normal_state_n              // number of normal state segments
  );

  void setup(
      bool updatePrior,
      bool updateMu,
      bool updateSigma,
      bool updateTransMat,
      bool fwdOnly,
      bool scaled,
      bool transTypeStat
      );

  /**
   * @brief Gives the CopyNumberHMM object the diagonal probabilities
   * it should use.
   */
  void diagProb(const std::vector<double>& dpvec);
  /**
   * @brief Runs HMM algorithm on a set of log2 ratios for a given sample
   * on a given chromosome.
   * Note: currently, only "state" will be used outside this method.
   */
  void hmm(
      const std::vector<int>&      position,   // SNP positions
      const VecCNLR&              samplr,     // Sample log2 ratios
      std::vector<int>&           state,
      const std::vector<double>&  priorProb,
      const std::vector<double>&  mu,
      const std::vector<double>&  sigma,
      double                      transMatDecay,
      const std::string&          stateEstMethod,
      std::vector<double>&        stateMargProb,
      const unsigned int          emIterations,
      const double                emConvergenceThresh,
      bool                        b_allelespecific,
      int                         normalState
  );
  /**
   * @brief This is done on a per chromosome basis for a given sample.
   */
  void smoothOutliers(
                      std::vector<int>&           state,
                      const std::vector<double>&  pos,
                      double                      distThresh
  );
  /**
   * @brief This is done across all chromosomes for a given sample.
   * The indexes for the "vector< vector<TYPE> >" items are
   * [chromosome][state/log2ratio/medianCN/pValue]
   * Example:  For chromosome index 2 the 5th log2 ratio values would
   * be "samplr[2][4]"
   */
  void postProcessHMM(
                      const std::vector< std::vector<int> >&   state,
                      const VecVecCNLR&                        logRatio,
                      int                                      normalState,
                      int                                      normalStateMinObservations,
                      bool                                     estCorrectionFactor,
                      double&                                  correctionFactor,
                      VecVecCNLR&                              medianCN
  );


};

#endif


