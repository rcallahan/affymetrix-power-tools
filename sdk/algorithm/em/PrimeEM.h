////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
 * @file   PrimeEM.h
 * @author Xiaojun Di
 * @date   Tue Aug 14 2006
 * 
 * DRAFT: wraps the EM algorithmics for genotyping under SDK
 * 
 */

//
#ifndef __PRIMEEM_H__
#define __PRIMEEM_H__

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <vector>
//

/**
 * @brief Object to hold all prior information to constitute an EM seed 
 */
class CEMSeed
{
	friend class CPrimeEM;
public:
	/**
	* @brief, the default constructor
	* @param clusters - number of clusters to cluster data into 
	*/
	CEMSeed(const int clusters = 3); 
	/**
	* @brief set the number of cluters
	* @param clusters - number of clusters to cluster data into
	*/
	void set(const int clusters = 3); 
	/**
	* @brief assignment operator
	*/
	CEMSeed& operator = (const CEMSeed& rhs);
	/**
	* @brief equal operator
	*/
	bool operator == (const CEMSeed& rhs);
	/**
	* @brief Set the centers of clusters assuming upto three clusters in total
	* @param c1 - center for the first cluster
	* @param c2 - center for the second cluster
	* @param c3 - center for the third cluster
	* @param c[] - array of centers of all clusters
	* @param c - vector of centers of all clusters
	*/
	void setMu(const double c[], const int size = 3);
	void setMu(const double c1, const double c2,const double c3);
	void setMu(std::vector<double>& c);
	/**
	* @brief Set the minimum of centers of clusters assuming upto three clusters in total
	* @param c1 - minimum of center for the first cluster
	* @param c2 - minimum of center for the second cluster
	* @param c3 - minimum of center for the third cluster
	* @param c[] - array of minimum of centers for all clusters
	* @param c - vector of minimum of centers of all clusters
	*/	
	void setMinMu(const double c[], const int size = 3);
	void setMinMu(const double c1,const double c2,const double c3);
	void setMinMu(std::vector<double>& c);
	/**
	* @brief Set the maximum of centers of clusters assuming upto three clusters in total
	* @param c1 - maximum of center for the first cluster
	* @param c2 - maximum of center for the second cluster
	* @param c3 - maximum of center for the third cluster
	* @param c[] - array of maximum of centers for all clusters
	* @param c - vector of maximum of centers of all clusters
	*/	
	void setMaxMu(const double c[], const int size = 3);
	void setMaxMu(const double c1,const double c2,const double c3);
	void setMaxMu(std::vector<double>& c);
	/**
	* @brief Set the weight of clusters assuming upto three clusters in total
	* @param w1 - weight for the first cluster
	* @param w2 - weight of center for the second cluster
	* @param w3 - weight of center for the third cluster
	* @param w[] - array of weights for all clusters
	* @param w - vector of weights  of all clusters
	*/	
	void setWeight(const double w[], const int size = 3);
	void setWeight(const double w1,const double w2,const double w3);
	void setWeight(std::vector<double>& w);
	/**
	* @brief Set the sigmas of clusters assuming upto three clusters in total
	* @param s1 - sigma for the first cluster
	* @param s2 - sigma for the second cluster
	* @param s3 - sigma for the third cluster
	* @param s[] - array of sigmas for all clusters
	* @param s - vector of sigmas of all clusters
	*/		
	void setSigma(const double s[], const int size = 3);
	void setSigma(const double s1,const double s2,const double s3);
	void setSigma(std::vector<double>& s);
	/**
	* @brief Set the minimum of sigmas of clusters assuming upto three clusters in total
	* @param s1 - minimum of sigma for the first cluster
	* @param s2 - minimum of sigma for the second cluster
	* @param s3 - minimum of sigma for the third cluster
	* @param s[] - array of minimum sigmas for all clusters
	* @param s - vector of minimum sigmas of all clusters
	*/		
	void setMinSigma(const double s[], const int size = 3);
	void setMinSigma(const double s1,const double s2,const double s3);
	void setMinSigma(std::vector<double>& s);
	/**
	* @brief Set the maximum of sigmas of clusters assuming upto three clusters in total
	* @param s1 - maximum of sigma for the first cluster
	* @param s2 - maximum of sigma for the second cluster
	* @param s3 - maximum of sigma for the third cluster
	* @param s[] - array of maximum sigmas for all clusters
	* @param s - vector of maximum sigmas of all clusters
	*/	
	void setMaxSigma(const double s[], const int size = 3);
	void setMaxSigma(const double s1,const double s2,const double s3);
	void setMaxSigma(std::vector<double>& s);
	/** 
	* @brief get functions to all member variables
	* @return the contents of member variables
	*/
	std::vector<double> getMu() { return m_mu; }
	std::vector<double> getSigma() {return m_sigma; }
	std::vector<double> getWeight() { return m_weight; }
	std::vector<double> getMinMu() { return m_minMu; }
	std::vector<double> getMinSigma() { return m_minSigma; }
	std::vector<double> getMaxMu() { return m_maxMu; }
	std::vector<double> getMaxSigma() { return m_maxSigma; }
protected:
	// initial values for EM to start with
	/// centers of clusters
	std::vector<double> m_mu;
	/// sigmas of clusters
	std::vector<double> m_sigma;
	/// weights of clusters
	std::vector<double> m_weight;

	// constraints to restrict EM's move
	/// minimums of cluster centers
	std::vector<double> m_minMu;
	/// minimums of cluter sigmas
	std::vector<double> m_minSigma;
	/// maximums of cluster centers
	std::vector<double> m_maxMu;
	/// maximums of cluster sigmas
	std::vector<double> m_maxSigma;
};

/**
 * @brief Object to hold all EM posterior information including some classification quality stats
 */
class CEMEst 
{
	friend class CPrimeEM;
public:
	/**
	* @brief, the default constructor
	* @param clusters - number of clusters to cluster data into 
	*/
	CEMEst(const int clusters = 3); 
	/**
	* @brief set the number of cluters
	* @param clusters - number of clusters to cluster data into
	*/
	void set(const int clusters = 3); 
	/**
	* @brief assignment operator
	*/
	CEMEst& operator = (const CEMEst& rhs);
protected:
	/**
	* @brief reset the buffer size for posterior information
	* @param dim - dimension of the data to be classified
	*/
	void reset(const int dim);
public:
	/// cluster centers
	std::vector<double> m_mu;
	/// cluster sigmas
	std::vector<double> m_sigma;
	/// cluster weights
	std::vector<double> m_weight;
	/// minimum of cluster centers
	std::vector<double> m_minMu;
	/// maximums of cluster centers
	std::vector<double> m_maxMu;
	/// minimums of cluster sigmas
	std::vector<double> m_minSigma;
	/// maximums of cluster sigmas
	std::vector<double> m_maxSigma;
	/// number of clusters remaining
	std::vector<bool> m_cluster;
	/// EM cluster member assignment
	std::vector<int> m_class;
	/// probability to best fit cluster
	std::vector<double> m_prob;
	/// probability ratio, best vs. second best
	std::vector<double> m_probRatio;
	/// genotype confidence
	std::vector<double> m_confidence;
	/// total maximum likelihood
	double m_ML;
	// QC metrices, Fisher's linear discriminants
	/// Fisher's linear discriminant of the first two clusters
	double m_FLD1;
	/// Fisher's linear discriminant of the second two clusters
	double m_FLD2;
	/// separation of the clustering results
	double m_sep;
};

/////////////////////////////////////////////////////////////////////////////////////////////
// compare two pairs based only on first value
struct PairLess : public std::binary_function<std::pair<double,int>,std::pair<double,int>,bool>
{
	bool operator () (const std::pair<double,int>& x, const std::pair<double,int>& y) const
	{
		return (x.first < y.first);
	}
};

/**
* @brief different loss function options
*/
typedef enum {
	ML=1,
	AIC=2,
	BIC=4,
	BIC_SEP=8
} LOSSFUNC;

	/**
 * @brief Object to hold all EM parameters
 */
class CEMParam
{
public:
	/**
	* @brief defualt constructor
	* @param clusters - number of clusters
	*/
	CEMParam(const int clusters = 3); 
	/**
	* @brief assignment operator
	*/
	CEMParam& operator = (const CEMParam& rhs);
public:
	/// maximum number of clusters
	int m_nClusters;
	/// maximum number of iterations
	int m_nMaxIter;
	/// tolerance for insignificant change for convergence
	double m_fTolerance;
	/// to where we place the non-Gaussian tail
	double m_tailCut;
	/// flag to indicate weight is used or not
	bool m_useWeight;
	/// factor to balance the weights in confidence of probability and loglikelihood ratio
	double m_factor;
	/// parameter used for seed generation
	double m_delta;
	/// parameter used for seed generation
	double m_epsilon;
	/// the loss function to be used for EM evaluation
	LOSSFUNC m_lossFunction;
	/// the threshold to drop a cluster if total number of data points is less
	double m_clusterDropThres;
	/// the threshold on EM's confidence
	double m_threshold;
	/// factor for hom cluster variance shrinkage
	double m_HomReductionCoef;
	/// factor for het cluster variance shrinkage
	double m_HetReductionCoef;
	/// factor for both hom and het cluster variance shrinkage
	double m_SigmaPower;
    // factors to check if any cluster std have gone above x times of geometric mean of all stds
    std::vector<double> m_geo_ratio_multiplier;
    // factors to check if any cluster std have gone below x times of geometric mean of all stds
    std::vector<double> m_geo_ratio_divisor;
};

/**
 * @brief Object to hold all EM prior information
 */
class CEMPrior
{
	friend class CPrimeEM;
public:
	/**
	* @brief defualt constructor
	* @param clusters - number of clusters
	*/
	CEMPrior(const int clusters = 3);
	/**
	* @brief initialize prior information with defaults if clusters=3
	* @return void
	*/
	void initialize();
	/**
	* @brief set function to initialize buffers
	* @param clusters - number of clusters
	* @return void
	*/
	void set(const int clusters = 3); 
	/**
	* @brief assignment operator
	*/
	CEMPrior& operator = (const CEMPrior& rhs);
	/**
	* @brief function to set up the centers of clusters
	* @param c - a vector to hold all centers of clusters
	* @return void
	*/
	void setMu(const std::vector<double>& c);
	/**
	* @brief function to set up the stdev of cluster centers
	* @param c - a vector to hold all centers of clusters
	* @return void
	*/
	void setMuSD(const std::vector<double>& c);
	/**
	* @brief function to set up the sigmas of clusters
	* @param s - a vector to hold all of stdevs of clusters
	* @return void
	*/
	void setSigma(const std::vector<double>& s);
	/**
	* @brief function to set up the centers of clusters
	* @param s - a vector to hold all stdevs of sigmas
	* @return void
	*/
	void setSigmaSD(const std::vector<double>& s);
	/** 
	* @brief function to set the prior status, 
	*		must set true with valid information to proceed if initialzie()
	*		not called
	* @param good - flag to indicate the prior information status
	* @return void
	*/
	inline void setStatus(bool good) { m_good = good; }
protected:
	/// prior status
	bool m_good;
	/// cluster centers
	std::vector<double> m_mu;
	/// stdev of cluster centers
	std::vector<double> m_muSD;
	/// cluster sigmas
	std::vector<double> m_sigma;
	/// stdev of cluster sigmas
	std::vector<double> m_sigmaSD;
};

/**
 * @brief The primary class to hold all EM operations
 */
class CPrimeEM
{
public:
	/**
	* @brief default constructor
	* @param clusters - number of clusters
	*/
	CPrimeEM(); 
	/**
	* @brief generate variable number of EM seeds based on some prior information
	*		this function is designed specifically for 3-clusters scenario, hence only 
	*		works if three-clusters is specified, otherwise return false and do nothing
	* @return true if no error, false otherwise
	*/
	bool generateEmpiricalEMSeeds();
	/**
	* @brief The function to run EM for genotyping, function only works for three cluster cases
	* @return true if there is no error, false otherwise
	*/
	bool goPrimeEM();
	/**
	* @brief the function to get EM's estimate with a specific initial seed
	* @param seed - the specific initial seed
	* @return true if there is no error, false otherwise
	*/
	bool EMEstimate(CEMSeed& seed);
	/**
	* @brief get to posterior
	*/
	inline CEMEst* getEMEstimates() { return &m_est;}
	/**
	* @brief get seeds generated by generateEmpiricalEMSeeds()
	* @return void
	*/
	inline std::vector<CEMSeed> getSeeds() { return m_seeds; }
	/**
	* @brief reset the data dimension
	* @param dim - data dimenstion
	* @return true if no error, false otherwise
	*/
	bool reset(const int dim);
	/**
	* @brief set the data to be classified
	* @param x - the data to be classified
	* @return void
	*/
	void setData(const std::vector<double>& x); 
	/** 
	* @brief set the prior information
	* @param prior - the prior information
	* @return void
	*/
	inline void setPrior(const CEMPrior& prior){ m_prior = prior;}
	/** 
	* @brief set EM parameters
	* @param param - the parameters
	* @return void
	*/
	void setParam(const CEMParam& param);
	/**
	* @brief check integrety of settings
	* @return true if settings are valid and consistent, false otherwise
	*/
	bool checkPrior();
	/**
	* @brief confidence threshold set function, for backward compatibility
	*/
	void setThreshold(const double thresh) { m_param.m_threshold = thresh; }

public:
	/// the data to be classified
	std::vector<double> m_data;
	/**
	* @brief the Expectation step of EM algorithm
	* @param final - flag to indicate it's step to calculate posterior
	* @return true if no error, false otherwise
	*/
	bool doEstep(bool final=false);
	/**
	* @brief the Maximization step of EM algorithm
	* @return true if no error, false otherwise
	*/
	bool doMstep();
	/** 
	* @brief check the convergence status
	* @param it - current iteration
	* @return true if not converged, false otherwise
	*/
	bool checkStatus(const int it);
	/**
	* @brief compute all posterior estimates
	* @return true if no error, false otherwise
	*/
	bool computeOptimalEstimates();
	/**
	* @brief compare the current estimates to the best by far
	* @param bestEst - best estimates so far
	* @return true if current is better
	*/
	bool evaluateEstimates(CEMEst& bestEst);
	/**
	* @brief roll back to the initial values
	* @param seed - the seed which contains the initial info to roll back to
	* @return true if no error, false otherwise
	*/
	bool useInitialEstimate(CEMSeed& seed);

	/**
	* @brief update cluster centers 
	* @return true if no error, false otherwise
	*/
	bool updateMu();
	/**
	* @brief update cluster variances 
	* @return true if no error, false otherwise
	*/
	bool updateSigma();
	/**
	* @brief update cluster weights 
	* @return true if no error, false otherwise
	*/
	bool updateWeight();
	/**
	* @brief adjust cluster centers according to constraints
	* @return true if no error, false otherwise
	*/
	bool adjustMu();
	/**
	* @brief adjust cluster variances according to constraints
	* @return true if no error, false otherwise
	*/
	bool adjustSigma();
	/**
	* @brief validate and adjust seed according to constraints
	* @param seed - the seed to be adjusted and returned after adjustment
	* @return true if no error, false otherwise
	*/
	bool validate(CEMSeed& seed);
	// general functions
	/**
	* @brief initialize the type and likelihood pair vector, assign likelihood to infinity
	* @param typeLikelihoodPair - the type and likelihood pair vector to be initialized
	* @return void
	*/
	void initializeTypeLikelihoodPair(std::vector<std::pair<double,int> >& typeLikelihoodPair);
	/**
	* @brief assign likelihood to each cluster for data point i
	* @param typeLikelihoodPair - the type and likelihood pair container
	* @param i - data point index
	* @return void
	*/
	void assignLogLikelihood(std::vector<std::pair<double,int> >& typeLikelihoodPair, const int i);
	/**
	* @brief compute Fisher's Linear discrimants
	*		 only computed if clusters are 2 or three
	* @param FLD1 - FLD for the first two clusters
	* @param FLD2 - FLD for the second two clusters
	* @return void
	*/
	void computeFLD();
	/**
	* @brief assign class membership and associated confidence to data point i
	* @param typeLikelihoodPair - the type and likelihood pair container
	* @param i - data point index
	* @return void
	*/	
	void assignClassAndConfidence(const int i, std::vector<std::pair<double,int> >& typeLikelihoodPair);
	/**
	* @brief adjust class membership and associated confidence to data point i
	* @param i - data point index
	* @return void
	*/	
	void adjustClassAndConfidence(const int i);

public:
	/// container of probabilities
	std::vector<std::vector<double> > m_p;
	/// container of normalized probabilities
	std::vector<std::vector<double> > m_np;
	/// container of total cluster probabilities
	std::vector<double> m_clusterP;
	/// container of total cluster normalized probabilities
	std::vector<double> m_clusterNP;
	/// container of total probabilities
	double m_sumNP;
	/// current estimates
	CEMEst m_est;
	/// pointer to the current seed in use
	CEMSeed* m_pSeed;
protected:
	/// number of data points
	int m_nDim;
	/// parameter to the scale of FLD penalty of a classification, calculated from m_nDim on the fly
	double m_lamda;
	/// container to hold all EM seeds 
	std::vector<CEMSeed> m_seeds;
	/// the generic prior, kept in the level to be adjustable case specific intead of generic
	CEMPrior m_prior;
	/// EM parameters
	CEMParam m_param;
	/// best estimates so far
	CEMEst m_preEst;
};

/**
* Utility function to calculate value from Gaussian function 
* with an exponential tail
* @param z - the z value
* @param zCut - where the exponential part starts
* @param factor - a scale factor
* @return the value of z
*/
double calcProbability(const double x,const double mu,const double sigma,const double weight,const double zCut,bool useWeight=false);
/**
* Utility function to calculate probability from Gaussian distibution 
* with an exponential tail given the mean and standard deviation
*
* @param x - the observation
* @param mu - mean of the gaussian distribution
* @param sigma - standard deviation of the gaussian distribution
* @param weight - the weight
* @param zCut - where the exponential part starts
* @param useWeight - flag to indicate whether weight should be used or not
* @return the probability of x from the given Gaussian distribution
*/
double calcGaussianWithExpTail(const double z,const double zCut,const double factor=1) ;

#endif //__PRIMEEM_H__
