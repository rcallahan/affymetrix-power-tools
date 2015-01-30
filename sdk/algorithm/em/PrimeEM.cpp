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
 * @file   PrimeEM.cpp
 * @author Xiaojun Di
 * @date   Tue Aug 14 11:56 2006
 * 
 * DRAFT: wraps the EM algorithmics for genotyping under SDK
 * 
 */

#include "algorithm/em/PrimeEM.h"
//
#include <iostream>
//

using namespace std;

// some global constants
/// NA
#define NA	-1*numeric_limits<double>::max()
/// maximum of doubles
#define MAX_DOUBLE	numeric_limits<double>::max()
/// infinity of doubles
#define INF	numeric_limits<double>::infinity()
/// very small positive double number
#define EPSILON	numeric_limits<double>::epsilon()

/**
* Utility function to calculate value from Gaussian function 
* with an exponential tail
* @param z - the z value
* @param zCut - where the exponential part starts
* @param factor - a scale factor
* @return the value of z
*/
double calcGaussianWithExpTail(const double z,const double zCut,const double factor) 
{
	double absz = fabs(z);
	double exponent = 0;
	if(absz < zCut) 
	{
		exponent = -0.5f*absz*absz*factor; 
	}
	else 
	{
		exponent = -0.5f*zCut*absz*factor;
	} 

	if(exponent < -200) {
		exponent = -200;
	}
	return exp(exponent);
}
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
double calcProbability(const double x,
                       const double mu,
                       const double sigma,
                       const double weight,
                       const double zCut,
                       bool useWeight)
{
	double z = (mu - x)/sigma;
	double p = 0;
	double exponent = calcGaussianWithExpTail(z,zCut,1);
	if(useWeight) 
	{
		p = (weight * exponent)/sigma;
	}
	else 
	{
		p = exponent/sigma;
	}
	return p;
}
/**
* @brief, the default constructor
* @param clusters - number of clusters to cluster data into 
*/
CEMSeed::CEMSeed(const int clusters) 
{
	m_mu.resize(clusters);
	m_sigma.resize(clusters);
	m_weight.resize(clusters);
	m_minMu.resize(clusters);
	m_minSigma.resize(clusters);
	m_maxMu.resize(clusters);
	m_maxSigma.resize(clusters);
}
/**
* @brief set the number of cluters
* @param clusters - number of clusters to cluster data into
*/
void CEMSeed::set(const int clusters) 
{
	m_mu.resize(clusters);
	m_sigma.resize(clusters);
	m_weight.resize(clusters);
	m_minMu.resize(clusters);
	m_minSigma.resize(clusters);
	m_maxMu.resize(clusters);
	m_maxSigma.resize(clusters);
}
/**
* @brief assignment operator
*/
CEMSeed& CEMSeed::operator = (const CEMSeed& rhs)
{
	copy(rhs.m_mu.begin(),rhs.m_mu.end(),m_mu.begin());
	copy(rhs.m_sigma.begin(),rhs.m_sigma.end(),m_sigma.begin());
	copy(rhs.m_weight.begin(),rhs.m_weight.end(),m_weight.begin());
	copy(rhs.m_minMu.begin(),rhs.m_minMu.end(),m_minMu.begin());
	copy(rhs.m_minSigma.begin(),rhs.m_minSigma.end(),m_minSigma.begin());
	copy(rhs.m_maxMu.begin(),rhs.m_maxMu.end(),m_maxMu.begin());
	copy(rhs.m_maxSigma.begin(),rhs.m_maxSigma.end(),m_maxSigma.begin());

	return *this;
}
/**
* @brief equal operator
*/
bool CEMSeed::operator == (const CEMSeed& rhs)
{
	bool ret = true;
	ret &= equal(rhs.m_mu.begin(),rhs.m_mu.end(),m_mu.begin());
	ret &= equal(rhs.m_sigma.begin(),rhs.m_sigma.end(),m_sigma.begin());
	ret &= equal(rhs.m_weight.begin(),rhs.m_weight.end(),m_weight.begin());
	ret &= equal(rhs.m_maxMu.begin(),rhs.m_maxMu.end(),m_maxMu.begin());
	ret &= equal(rhs.m_minMu.begin(),rhs.m_minMu.end(),m_minMu.begin());
	ret &= equal(rhs.m_minSigma.begin(),rhs.m_minSigma.end(),m_minSigma.begin());
	ret &= equal(rhs.m_maxSigma.begin(),rhs.m_maxSigma.end(),m_maxSigma.begin());
	return ret;
}
/**
* @brief Set the centers of clusters assuming upto three clusters in total
* @param c1 - center for the first cluster
* @param c2 - center for the second cluster
* @param c3 - center for the third cluster
* @param c[] - array of centers of all clusters
* @param c - vector of centers of all clusters
*/
void CEMSeed::setMu(const double c[], const int size)
{
	m_mu.resize(size);
	copy(c,&c[size],m_mu.begin());
}
void CEMSeed::setMu(const double c1, const double c2,const double c3)
{
	m_mu.resize(3);
	m_mu[0] = c1;
	m_mu[1] = c2;
	m_mu[2] = c3;
}
void CEMSeed::setMu(std::vector<double>& c)
{
	m_mu.resize(c.size());
	copy(c.begin(),c.end(),m_mu.begin());
}
/**
* @brief Set the minimum of centers of clusters assuming upto three clusters in total
* @param c1 - minimum of center for the first cluster
* @param c2 - minimum of center for the second cluster
* @param c3 - minimum of center for the third cluster
* @param c[] - array of minimum of centers for all clusters
* @param c - vector of minimum of centers of all clusters
*/	
void CEMSeed::setMinMu(const double c[], const int size)
{
	m_minMu.resize(size);
	copy(c,&c[size],m_minMu.begin());
}
void CEMSeed::setMinMu(const double c1,const double c2,const double c3)
{
	m_minMu.resize(3);
	m_minMu[0] = c1;
	m_minMu[1] = c2;
	m_minMu[2] = c3;
}
void CEMSeed::setMinMu(std::vector<double>& c)
{
	m_minMu.resize(c.size());
	copy(c.begin(),c.end(),m_minMu.begin());
}	
/**
* @brief Set the maximum of centers of clusters assuming upto three clusters in total
* @param c1 - maximum of center for the first cluster
* @param c2 - maximum of center for the second cluster
* @param c3 - maximum of center for the third cluster
* @param c[] - array of maximum of centers for all clusters
* @param c - vector of maximum of centers of all clusters
*/	
void CEMSeed::setMaxMu(const double c[], const int size)
{
	m_maxMu.resize(size);
	copy(c,&c[size],m_maxMu.begin());
}
void CEMSeed::setMaxMu(const double c1,const double c2,const double c3)
{
	m_maxMu.resize(3);
	m_maxMu[0] = c1;
	m_maxMu[1] = c2;
	m_maxMu[2] = c3;
}
void CEMSeed::setMaxMu(std::vector<double>& c)
{
	m_maxMu.resize(c.size());
	copy(c.begin(),c.end(),m_maxMu.begin());
}	
/**
* @brief Set the weight of clusters assuming upto three clusters in total
* @param w1 - weight for the first cluster
* @param w2 - weight of center for the second cluster
* @param w3 - weight of center for the third cluster
* @param w[] - array of weights for all clusters
* @param w - vector of weights  of all clusters
*/	
void CEMSeed::setWeight(const double w[], const int size)
{
	m_weight.resize(size);
	copy(w,&w[size],m_weight.begin());
}
void CEMSeed::setWeight(const double w1,const double w2,const double w3)
{
	m_weight.resize(3);
	m_weight[0] = w1;
	m_weight[1] = w2;
	m_weight[2] = w3;
}
void CEMSeed::setWeight(std::vector<double>& w)
{
	m_weight.resize(w.size());
	copy(w.begin(),w.end(),m_weight.begin());
}		
/**
* @brief Set the sigmas of clusters assuming upto three clusters in total
* @param s1 - sigma for the first cluster
* @param s2 - sigma for the second cluster
* @param s3 - sigma for the third cluster
* @param s[] - array of sigmas for all clusters
* @param s - vector of sigmas of all clusters
*/		
void CEMSeed::setSigma(const double s[], const int size)
{
	m_sigma.resize(size);
	copy(s,&s[size],m_sigma.begin());
}
void CEMSeed::setSigma(const double s1,const double s2,const double s3)
{
	m_sigma.resize(3);
	m_sigma[0] = s1;
	m_sigma[1] = s2;
	m_sigma[2] = s3;
}
void CEMSeed::setSigma(std::vector<double>& s)
{
	m_sigma.resize(s.size());
	copy(s.begin(),s.end(),m_sigma.begin());
}		
/**
* @brief Set the minimum of sigmas of clusters assuming upto three clusters in total
* @param s1 - minimum of sigma for the first cluster
* @param s2 - minimum of sigma for the second cluster
* @param s3 - minimum of sigma for the third cluster
* @param s[] - array of minimum sigmas for all clusters
* @param s - vector of minimum sigmas of all clusters
*/		
void CEMSeed::setMinSigma(const double s[], const int size)
{
	m_minSigma.resize(size);
	copy(s,&s[size],m_minSigma.begin());
}
void CEMSeed::setMinSigma(const double s1,const double s2,const double s3)
{
	m_minSigma.resize(3);
	m_minSigma[0] = s1;
	m_minSigma[1] = s2;
	m_minSigma[2] = s3;
}
void CEMSeed::setMinSigma(std::vector<double>& s)
{
	m_minSigma.resize(s.size());
	copy(s.begin(),s.end(),m_minSigma.begin());
}		
/**
* @brief Set the maximum of sigmas of clusters assuming upto three clusters in total
* @param s1 - maximum of sigma for the first cluster
* @param s2 - maximum of sigma for the second cluster
* @param s3 - maximum of sigma for the third cluster
* @param s[] - array of maximum sigmas for all clusters
* @param s - vector of maximum sigmas of all clusters
*/	
void CEMSeed::setMaxSigma(const double s[], const int size)
{
	m_maxSigma.resize(size);
	copy(s,&s[size],m_maxSigma.begin());
}
void CEMSeed::setMaxSigma(const double s1,const double s2,const double s3)
{
	m_maxSigma.resize(3);
	m_maxSigma[0] = s1;
	m_maxSigma[1] = s2;
	m_maxSigma[2] = s3;
}
void CEMSeed::setMaxSigma(std::vector<double>& s)
{
	m_maxSigma.resize(s.size());
	copy(s.begin(),s.end(),m_maxSigma.begin());
}		

/**
* @brief, the default constructor
* @param clusters - number of clusters to cluster data into 
*/
CEMEst::CEMEst(const int clusters) 
{
	m_mu.resize(clusters);
	m_sigma.resize(clusters);
	m_weight.resize(clusters);
	m_minMu.resize(clusters);
	m_minSigma.resize(clusters);
	m_maxMu.resize(clusters);
	m_maxSigma.resize(clusters);
	m_cluster.resize(clusters);
	m_ML = -1;
	m_FLD1 = NA;
	m_FLD2 = NA;
	m_sep = NA;
}
/**
* @brief set the number of cluters
* @param clusters - number of clusters to cluster data into
*/
void CEMEst::set(const int clusters) 
{
	m_mu.resize(clusters);
	m_sigma.resize(clusters);
	m_weight.resize(clusters);
	m_minMu.resize(clusters);
	m_minSigma.resize(clusters);
	m_maxMu.resize(clusters);
	m_maxSigma.resize(clusters);
	m_cluster.resize(clusters);
	m_ML = -1;
	m_FLD1 = NA;
	m_FLD2 = NA;
	m_sep = NA;
}
/**
* @brief assignment operator
*/
CEMEst& CEMEst::operator = (const CEMEst& rhs)
{
	copy(rhs.m_mu.begin(),rhs.m_mu.end(),m_mu.begin());
	copy(rhs.m_sigma.begin(),rhs.m_sigma.end(),m_sigma.begin());
	copy(rhs.m_weight.begin(),rhs.m_weight.end(),m_weight.begin());
	copy(rhs.m_minMu.begin(),rhs.m_minMu.end(),m_minMu.begin());
	copy(rhs.m_minSigma.begin(),rhs.m_minSigma.end(),m_minSigma.begin());
	copy(rhs.m_maxMu.begin(),rhs.m_maxMu.end(),m_maxMu.begin());
	copy(rhs.m_maxSigma.begin(),rhs.m_maxSigma.end(),m_maxSigma.begin());
	copy(rhs.m_cluster.begin(),rhs.m_cluster.end(),m_cluster.begin());
	copy(rhs.m_class.begin(),rhs.m_class.end(),m_class.begin());
	copy(rhs.m_prob.begin(),rhs.m_prob.end(),m_prob.begin());
	copy(rhs.m_probRatio.begin(),rhs.m_probRatio.end(),m_probRatio.begin());
	copy(rhs.m_confidence.begin(),rhs.m_confidence.end(),m_confidence.begin());
	m_ML = rhs.m_ML;
	m_FLD1 = rhs.m_FLD1;
	m_FLD2 = rhs.m_FLD2;
	m_sep = rhs.m_sep;

	return *this;
}
/**
* @brief reset the buffer size for posterior information
* @param dim - dimension of the data to be classified
*/
void CEMEst::reset(const int dim)
{
	m_class.resize(dim);
	m_confidence.resize(dim);
	m_probRatio.resize(dim);
	m_prob.resize(dim);
}
/**
* @brief defualt constructor
* @param clusters - number of clusters
*/
CEMParam::CEMParam(const int clusters) 
{
	m_fTolerance = 0.001f;
	m_nMaxIter = 50;
	m_factor = 0.01f;
	m_nClusters = clusters;
	m_lossFunction = BIC_SEP;
	m_tailCut = 2.5f;
	m_useWeight = false;
	m_threshold = 0.05f;
	m_clusterDropThres = 0.8f;
	m_delta = 0.003f;
	m_epsilon = 0.005f;
	m_HomReductionCoef = 2.0f;
	m_HetReductionCoef = 1.2f;
	m_SigmaPower = 0.8f;
	m_geo_ratio_multiplier.push_back(2.5f);
	m_geo_ratio_multiplier.push_back(3.0f);
	m_geo_ratio_multiplier.push_back(2.5f);
	m_geo_ratio_divisor.push_back(2.5f);
	m_geo_ratio_divisor.push_back(2.0f);
	m_geo_ratio_divisor.push_back(2.5f);
}
/**
* @brief assignment operator
*/
CEMParam& CEMParam::operator = (const CEMParam& rhs)
{
	m_nMaxIter = rhs.m_nMaxIter;
	m_factor = rhs.m_factor;
	m_nClusters = rhs.m_nClusters;
	m_lossFunction = rhs.m_lossFunction;
	m_tailCut = rhs.m_tailCut;
	m_useWeight = rhs.m_useWeight;
	m_threshold = rhs.m_threshold;
	m_delta = rhs.m_delta;
	m_epsilon = rhs.m_epsilon;
	m_HomReductionCoef = rhs.m_HomReductionCoef;
	m_HetReductionCoef = rhs.m_HetReductionCoef;
	m_geo_ratio_multiplier.resize(rhs.m_geo_ratio_multiplier.size());
	m_geo_ratio_divisor.resize(rhs.m_geo_ratio_divisor.size());
	copy(rhs.m_geo_ratio_multiplier.begin(),rhs.m_geo_ratio_multiplier.end(),m_geo_ratio_multiplier.begin());
	copy(rhs.m_geo_ratio_divisor.begin(),rhs.m_geo_ratio_divisor.end(),m_geo_ratio_divisor.begin());
	m_SigmaPower = rhs.m_SigmaPower;
	m_clusterDropThres = rhs.m_clusterDropThres;
	return *this;
}
/**
* @brief defualt constructor
* @param clusters - number of clusters
*/
CEMPrior::CEMPrior(const int clusters) 
{
	m_mu.resize(clusters);
	m_sigma.resize(clusters);
	m_muSD.resize(clusters);
	m_sigmaSD.resize(clusters);
	m_good = false; // default to invalid
}
/**
* @brief initialize prior information with defaults if clusters=3
* @return void
*/
void CEMPrior::initialize()
{
	// number of clusters has been set, if equal to 3 initialize the generic
	if(m_mu.size() == 3)
	{
		// default prior
		std::vector<double> x(3);
		// the static default prior
		double c[] = {-0.66f,0,0.66f}; // centers
		copy(&c[0], &c[3],x.begin());
		setMu(x);
		double c_s[] = {0.125f,0.175f,0.125f}; // stdev of centers
		copy(&c_s[0], &c_s[3],x.begin());
		setMuSD(x);
		double s[] = {0.065f,0.085f,0.065f}; // stdev
		copy(&s[0], &s[3],x.begin());
		setSigma(x);
		double s_s[] = {0.035f,0.015f,0.035f}; // stdev of stdev
		copy(&s_s[0], &s_s[3],x.begin());
		setSigmaSD(x);
		m_good = true;
	}
	else 
	{
		m_good = false; // don't know how to set up defaults
	}
}

/**
* @brief set function to initialize buffers
* @param clusters - number of clusters
* @return void
*/
void CEMPrior::set(const int clusters) 
{
	m_mu.resize(clusters);
	m_sigma.resize(clusters);
	m_muSD.resize(clusters);
	m_sigmaSD.resize(clusters);
}
/**
* @brief assignment operator
*/
CEMPrior& CEMPrior::operator = (const CEMPrior& rhs)
{
	int clusters = rhs.m_mu.size();
	m_mu.resize(clusters);
	m_sigma.resize(clusters);
	m_muSD.resize(clusters);
	m_sigmaSD.resize(clusters);

	copy(rhs.m_mu.begin(),rhs.m_mu.end(),m_mu.begin());
	copy(rhs.m_muSD.begin(),rhs.m_muSD.end(),m_muSD.begin());
	copy(rhs.m_sigma.begin(),rhs.m_sigma.end(),m_sigma.begin());
	copy(rhs.m_sigmaSD.begin(),rhs.m_sigmaSD.end(),m_sigmaSD.begin());
	return *this;
}
/**
* @brief function to set up the centers of clusters
* @param c - a vector to hold all centers of clusters
* @return void
*/
void CEMPrior::setMu(const std::vector<double>& c)
{
	m_mu.resize(c.size());
	copy(c.begin(),c.end(),m_mu.begin());
}
/**
* @brief function to set up the stdev of cluster centers
* @param c - a vector to hold all centers of clusters
* @return void
*/
void CEMPrior::setMuSD(const std::vector<double>& c)
{
	m_muSD.resize(c.size());
	copy(c.begin(),c.end(),m_muSD.begin());
}
/**
* @brief function to set up the sigmas of clusters
* @param s - a vector to hold all of stdevs of clusters
* @return void
*/
void CEMPrior::setSigma(const std::vector<double>& s)
{
	m_sigma.resize(s.size());
	copy(s.begin(),s.end(),m_sigma.begin());
}
/**
* @brief function to set up the centers of clusters
* @param c - a vector to hold all stdevs of sigmas
* @return void
*/
void CEMPrior::setSigmaSD(const std::vector<double>& s)
{
	m_sigmaSD.resize(s.size());
	copy(s.begin(),s.end(),m_sigmaSD.begin());
}
/**
* @brief default constructor
* @param clusters - number of clusters
*/
CPrimeEM::CPrimeEM() 
{
	m_nDim = 0;
	m_sumNP = -1;
	m_pSeed = 0;
	m_est.set(m_param.m_nClusters);
	m_prior.set(m_param.m_nClusters);
	m_clusterP.resize(m_param.m_nClusters);
	m_clusterNP.resize(m_param.m_nClusters);
	m_prior.initialize();
}
/**
* @brief reset the data dimension
* @param dim - data dimenstion
* @return true if no error, false otherwise
*/
bool CPrimeEM::reset(const int dim)
{
	// only reset if data dim changes
	if(m_nDim != dim)
	{
		m_nDim = dim;
		m_est.reset(dim);
		m_preEst.reset(dim);
		m_p.resize(dim);
		m_np.resize(dim);
		m_data.resize(dim);
		int i = 0;
		for(i=0; i < dim; i ++)
		{
			m_p[i].resize(m_param.m_nClusters);
			m_np[i].resize(m_param.m_nClusters);
		}
		m_lamda = sqrt((double) dim);
	}
	return true;
}
/**
* @brief check integrety of settings
* @return true if settings are valid and consistent, false otherwise
*/
bool CPrimeEM::checkPrior()
{
	bool status = (m_prior.m_mu.size() == m_prior.m_muSD.size());
	status &= (m_prior.m_mu.size() == m_prior.m_sigma.size());
	status &= (m_prior.m_mu.size() == m_prior.m_sigmaSD.size());
	status &= (m_prior.m_mu.size() == m_param.m_nClusters);
	if(status) m_prior.setStatus(true);
	else m_prior.setStatus(false);
	return status;
}
/** 
* @brief set EM parameters
* @param param - the parameters
* @return void
*/
void CPrimeEM::setParam(const CEMParam& param)
{
	m_param = param;
	m_est.set(m_param.m_nClusters);
	m_prior.set(m_param.m_nClusters);
	m_clusterP.resize(m_param.m_nClusters);
	m_clusterNP.resize(m_param.m_nClusters);
	m_prior.initialize();
}
/**
* @brief set the data to be classified
* @param x - the data to be classified
* @return void
*/
void CPrimeEM::setData(const std::vector<double>& x) 
{ 
	reset(x.size());
	copy(x.begin(),x.end(),m_data.begin());
}
/**
* Run EM for genotyping, first generate variable number of EM seeds
* based some prior information and the data to be classified, then
* use EM to get estimates starting from each seeds
* this function only works for three cluster cases as it is designed for 
* genotyping only.
*
* @return true if there is no error, false otherwise
*/
bool CPrimeEM::goPrimeEM()
{
	// check consistency
	if(!checkPrior()) return false;

	// if we are here we need generate seeds first
	if(!generateEmpiricalEMSeeds()) return false;

	size_t nSeeds = m_seeds.size();
//	if(nSeeds == 1) return false;
	if(nSeeds == 1)
	{
		// get EM estimates for this seed
		bool bret = EMEstimate(m_seeds[0]);
		accumulate(m_est.m_cluster.begin(),m_est.m_cluster.end(),0);
		return bret;
	}

	size_t it = 0;
	CEMEst bestEst;
	bestEst.set(m_param.m_nClusters);
	bestEst.reset(m_nDim);
	for(it = 0; it < nSeeds; it ++)
	{
		// get EM estimates for this seed
		EMEstimate(m_seeds[it]);

		// evaluate estimates
		if(it == 0)
			bestEst = m_est;
		else
			evaluateEstimates(bestEst);
	}

	// let's keep the best
	m_est = bestEst;

	return true;
}

/** 
* @brief Roll back to initial values for all estimates
* @param seed - the container with the initial values
* @return true is no error, false otherwise
*/
bool CPrimeEM::useInitialEstimate(CEMSeed& seed)
{
	// validate seed, modify if necessary
	validate(seed);

	// final updates
	doEstep(true);

	// get optimal estimates
	computeOptimalEstimates();

	// Fisher's discriminants based on EM posterior, only do this when three cluster is requested
	m_est.m_FLD1 = NA, m_est.m_FLD2 = NA;
	computeFLD();
     
	return true;
}

/**
* Start EM from the initial seed to get estimates for a given data set,
* all the posterior information will be calculated
*
* @param seed - the initial seed for EM to start with
* @return true if no error, false otherwise
*/
bool CPrimeEM::EMEstimate(CEMSeed& seed)
{
	// validate seed, modify if necessary
	validate(seed);

	// keep a global pointer of the current seed
	m_pSeed = &seed;

	int it = 0;
	while(it < m_param.m_nMaxIter)
	{
		// E step
		doEstep();

		// M step
		doMstep();

		if(!checkStatus(it))
			break;

		m_preEst = m_est;
		it++;
	}

	// final updates
	doEstep(true);

	// get optimal estimates
	computeOptimalEstimates();

	// Fisher's discriminants based on EM posterior, only do this when three cluster is requested
	m_est.m_FLD1 = NA, m_est.m_FLD2 = NA;
	computeFLD();
     
	return true;
}

/**
* @brief the Expectation step of EM algorithm
* @param final - flag to indicate it's step to calculate posterior
* @return true if no error, false otherwise
*/
bool CPrimeEM::doEstep(bool final)
{
	int j =0, i = 0;
	// sum of normalized probabilities for each cluster over all samples
	vector<double> sum_p;
	sum_p.resize(m_nDim);
	fill(sum_p.begin(),sum_p.end(),0.0f);
	double z = 0;

	// calculating probabilities for each sample in each cluster// and sum up for each sample
	for(j= 0; j < m_param.m_nClusters; j ++) 
	{
		if(m_est.m_cluster[j]) 
		{
			for(i = 0; i < m_nDim; i ++) 
			{
				z = fabs(m_data[i] - m_est.m_mu[j])/m_est.m_sigma[j];
				if(final || z < 50)
				{
					m_p[i][j] = calcProbability(m_data[i],m_est.m_mu[j],m_est.m_sigma[j],m_est.m_weight[j],m_param.m_tailCut,m_param.m_useWeight);
				}
				else
				{
				   m_p[i][j] = 0;
				}
				sum_p[i] = sum_p[i] + m_p[i][j];
			}
		}
	}

	// normalize probabilities for each sample in each cluster// and sum up for each cluster
	m_sumNP = 0;
	fill(m_clusterNP.begin(),m_clusterNP.end(),0.0f);
	for(j=0; j < m_param.m_nClusters; j ++) 
	{
		if(m_est.m_cluster[j]) 
		{
			for(i = 0; i < m_nDim; i ++) 
			{
				if(sum_p[i] > 0) {
					m_np[i][j] = m_p[i][j]/sum_p[i];
				}
				else 
				{
					m_np[i][j] = 0;
				}
				m_clusterNP[j] = m_clusterNP[j] + m_np[i][j];
			}
			m_sumNP += m_clusterNP[j];
		}
	}
   
    // check if we need to drop any clusters
    if(!final) 
	{
		for(j = 0; j < m_param.m_nClusters; j ++) 
		{
			if(m_clusterNP[j] < m_param.m_clusterDropThres)
    			m_est.m_cluster[j] = false;
		}
    }

	return true;
}
/**
* @brief the Maximization step of EM algorithm
* @return true if no error, false otherwise
*/
bool CPrimeEM::doMstep()
{
	bool bret = true;
	// we have the Estep, now we are ready to update mu and sigma
	bret &= updateMu();

	// adjust mu
	bret &= adjustMu();

	// update weights
	bret &= updateWeight();

	// update sigmas
	bret &= updateSigma();

	// adjust sigmas
	bret &= adjustSigma();

	return bret;
}
/** 
* @brief check the convergence status
* @param it - current iteration
* @return true if not converged, false otherwise
*/
bool CPrimeEM::checkStatus(const int it)
{
	if(it == 0) return true;

	bool changed = false;

	// check changes
	int j = 0;
	for(j=0; j < m_param.m_nClusters; j ++) 
	{
		if(m_est.m_cluster[j]) 
		{
			double w_diff = fabs(m_preEst.m_weight[j] - m_est.m_weight[j]);
			double m_diff = fabs(m_preEst.m_mu[j] - m_est.m_mu[j]);
			double s_diff = fabs(m_preEst.m_sigma[j] - m_est.m_sigma[j]);
			if(w_diff > m_param.m_fTolerance || m_diff > m_param.m_fTolerance || s_diff > m_param.m_fTolerance) 
			{
				changed = true;
				break;
			}
		}
	}

	return changed;
}
/**
* @brief initialize the type and likelihood pair vector, assign likelihood to infinity
* @param typeLikelihoodPair - the type and likelihood pair vector to be initialized
* @return void
*/
void CPrimeEM::initializeTypeLikelihoodPair(std::vector<std::pair<double,int> >& typeLikelihoodPair)
{
	int j = 0;
	for(j = 0; j < m_param.m_nClusters; j ++)
	{
		typeLikelihoodPair[j].first = MAX_DOUBLE;
		typeLikelihoodPair[j].second = j;
	}
}
/**
* @brief assign likelihood to each cluster for data point i
* @param typeLikelihoodPair - the type and likelihood pair container
* @param i - data point index
* @return void
*/
void CPrimeEM::assignLogLikelihood(std::vector<std::pair<double,int> >& typeLikelihoodPair, const int i)
{
	int j = 0;
	for(j = 0; j < m_param.m_nClusters; j ++)
	{
		if(m_est.m_cluster[j])
		{
			if(m_np[i][j] > 0)
				typeLikelihoodPair[j].first = -log(m_np[i][j]);
			else 
				typeLikelihoodPair[j].first = MAX_DOUBLE;
			typeLikelihoodPair[j].second = j;
		}
		else 
		{
			typeLikelihoodPair[j].first = MAX_DOUBLE;
			typeLikelihoodPair[j].second = j;
		}
	}
}
/**
* @brief compute Fisher's Linear discrimants
*		 only computed if clusters are 2 or three
* @param FLD1 - FLD for the first two clusters
* @param FLD2 - FLD for the second two clusters
* @return void
*/
void CPrimeEM::computeFLD()
{
	if(m_param.m_nClusters > 3 || m_param.m_nClusters < 2)
	{
		m_est.m_FLD1 = NA;
		m_est.m_FLD2 = NA;
		m_est.m_sep = NA;
		return;
	}

	if(m_est.m_cluster[0] && m_est.m_cluster[1]) 
	{
		m_est.m_FLD1 = (m_est.m_mu[1]-m_est.m_mu[0])/sqrt(m_est.m_sigma[1]*m_est.m_sigma[1]+m_est.m_sigma[0]*m_est.m_sigma[0]);
	}
	
	if(m_param.m_nClusters == 3)
	{
		if(m_est.m_cluster[1] && m_est.m_cluster[2]) 
		{
			m_est.m_FLD2 = (m_est.m_mu[2]-m_est.m_mu[1])/sqrt(m_est.m_sigma[2]*m_est.m_sigma[2]+m_est.m_sigma[1]*m_est.m_sigma[1]);
		}

		if(m_est.m_cluster[0] && m_est.m_cluster[2] && !m_est.m_cluster[1]) 
		{
			m_est.m_FLD2 = (m_est.m_mu[2]-m_est.m_mu[0])/sqrt(m_est.m_sigma[2]*m_est.m_sigma[2]+m_est.m_sigma[0]*m_est.m_sigma[0]);
			m_est.m_FLD2 = m_est.m_FLD2/4;
		}
	}
	else
	{
		m_est.m_FLD2 = NA;
	}

	if(m_est.m_FLD1 != NA && m_est.m_FLD2 != NA)
		m_est.m_sep = sqrt(m_est.m_FLD1*m_est.m_FLD2);
	else if(m_est.m_FLD1 != NA)
		m_est.m_sep = m_est.m_FLD1;
	else if(m_est.m_FLD2 != NA)
		m_est.m_sep = m_est.m_FLD2;
	else m_est.m_sep = NA;
}
/**
* @brief adjust class membership and associated confidence to data point i
* @param i - data point index
* @return void
*/
void CPrimeEM::adjustClassAndConfidence(const int i)
{
	int nValidClusters = accumulate(m_est.m_cluster.begin(),m_est.m_cluster.end(),0);

	// reassign so that 0==AA,1==AB,2==BB,-1==NC
	m_est.m_class[i] = 2 - m_est.m_class[i];
	if(m_est.m_class[i] > 2)
		m_est.m_class[i] = -1;

	// calculate confidence scores
	m_est.m_confidence[i] = 1 - m_est.m_prob[i];
	if(nValidClusters > 1)
	{
		double pRatio = m_est.m_probRatio[i]/(m_est.m_probRatio[i] + 1);
		m_est.m_confidence[i] = 1- m_est.m_prob[i] * pRatio;
	}
	// thresholding
	if(m_est.m_confidence[i] > m_param.m_threshold)
		m_est.m_class[i] = -1;
}
/**
* @brief assign class membership and associated confidence to data point i
* @param typeLikelihoodPair - the type and likelihood pair container
* @param i - data point index
* @return void
*/	
void CPrimeEM::assignClassAndConfidence(const int i, std::vector<std::pair<double,int> >& typeLikelihoodPair)
{
	// results holders
	int best = 0, second = -1, third = -1, temp = -1;

	// number of valid clusters
	int nValidClusters = accumulate(m_est.m_cluster.begin(),m_est.m_cluster.end(),0);

	if(nValidClusters > 1) 
	{
		best = typeLikelihoodPair[0].second, second = typeLikelihoodPair[1].second;
		double dotProd = (m_data[i] - m_est.m_mu[second])*(m_est.m_mu[best]-m_est.m_mu[second]);
		if(dotProd < 0) 
		{ //inconsistency found
            temp = best;
            best = second;
            second = temp;
		}
		if(nValidClusters > 2) 
		{
			third = typeLikelihoodPair[2].second;
            dotProd = (m_data[i] - m_est.m_mu[third])*(m_est.m_mu[best]-m_est.m_mu[third]);
            if(dotProd < 0) 
			{ //inconsistency found
				temp = second;
				second = best;
				best = third;
				third = temp;
            }
		}
	}
	m_est.m_class[i] = best;

	if(nValidClusters > 1) 
	{
		// check that the second most likely cluster is much less probable
		if(m_np[i][second] == 0) 
            m_est.m_probRatio[i] = 0;
		else 
            m_est.m_probRatio[i] = m_np[i][best]/m_np[i][second];
		
		// need to modify if the ratio is too large, do this later
		if(m_est.m_probRatio[i] == INF) 
			m_est.m_probRatio[i] = MAX_DOUBLE;
		else if(m_est.m_probRatio[i] > MAX_DOUBLE) 
			m_est.m_probRatio[i] = MAX_DOUBLE;
	}			
	
	double z = (m_data[i] - m_est.m_mu[best])/m_est.m_sigma[best];
	m_est.m_prob[i] = calcGaussianWithExpTail(z,m_param.m_tailCut,m_param.m_factor);
	
	// the probability to be used for seed quality
	double prob = calcProbability(m_data[i],m_est.m_mu[best],m_est.m_sigma[best],m_est.m_weight[best],m_param.m_tailCut);	
	prob = prob*m_est.m_sigma[best];
	double logProb = -log(prob);

	// sum up the best negative log probabilities over all samples
	// there is tweek waiting for Chad's response
	if(logProb > 20) 
		m_est.m_ML = m_est.m_ML + 20;
	else 
		m_est.m_ML = m_est.m_ML + logProb;   
	
	// adjust call and compute confidence
    // reassign so that 0==AA,1==AB,2==BB,-1==NC
	adjustClassAndConfidence(i);
}

/**
* @brief compute all posterior estimates
* @return true if no error, false otherwise
*/
bool CPrimeEM::computeOptimalEstimates()
{
	// number of valid clusters
	int nValidClusters = accumulate(m_est.m_cluster.begin(),m_est.m_cluster.end(),0);
		
	int i = 0;
	vector<pair<double,int> > typeLikelihoodPair(m_param.m_nClusters);
	if(nValidClusters > 0) 
	{ 
		// reset likelihood for each data point
		initializeTypeLikelihoodPair(typeLikelihoodPair);

		// calculate negative log likelihoods
		m_est.m_ML = 0;
		for(i=0; i < m_nDim; i ++) 
		{
			// calculate and assign log likelihood for each data point
			assignLogLikelihood(typeLikelihoodPair,i);

			// order class index based on likelihood
			sort(typeLikelihoodPair.begin(),typeLikelihoodPair.end(),PairLess());
			
			// build call and confidence for data point i
			assignClassAndConfidence(i,typeLikelihoodPair);
		}
	}
		
	return true;
}
/**
* @brief compare the current estimates to the best by far
* @param bestEst - best estimates so far
* @return true if current is better
*/
bool CPrimeEM::evaluateEstimates(CEMEst& bestEst)
{
	bool test = false;
	int nValidClusters = accumulate(m_est.m_cluster.begin(),m_est.m_cluster.end(),0);
	switch(m_param.m_lossFunction)
	{
	case ML:
		test = (m_est.m_ML < bestEst.m_ML);
		break;
	case AIC:
		test = (2*m_est.m_ML + 2*nValidClusters 
			< 2*bestEst.m_ML + 2*nValidClusters);
		break;
	case BIC:
		test = (2*m_est.m_ML + 2*m_lamda*nValidClusters 
			< 2*bestEst.m_ML + 2*m_lamda*nValidClusters);
		break;
	case BIC_SEP:
		test = (2*m_est.m_ML + 2*m_lamda*(nValidClusters-m_est.m_sep) 
			< 2*bestEst.m_ML + 2*m_lamda*(nValidClusters-bestEst.m_sep));
		break;
	default:
		test = (2*m_est.m_ML + 2*m_lamda*(nValidClusters-m_est.m_sep) 
			< 2*bestEst.m_ML + 2*m_lamda*(nValidClusters-bestEst.m_sep));
		break;
	}
	if(test)
		bestEst = m_est;

	return true;
}

/**
* @brief update cluster centers 
* @return true if no error, false otherwise
*/
bool CPrimeEM::updateMu()
{
    // calculate mu
	int j = 0, i = 0;
	fill(m_est.m_mu.begin(),m_est.m_mu.end(),0.0f);
	for(j = 0; j < m_param.m_nClusters; j++) 
	{
		// only work on valid clusters
		if(m_est.m_cluster[j]) 
		{
			for(i = 0; i < m_nDim; i ++)
				m_est.m_mu[j] = m_est.m_mu[j] + m_np[i][j] * m_data[i];
	
			m_est.m_mu[j] = m_est.m_mu[j]/m_clusterNP[j];
		}
	}

	return true;
}
/**
* @brief update cluster variances 
* @return true if no error, false otherwise
*/
bool CPrimeEM::updateSigma()
{
    double s_2 = 0;
	int i = 0, j = 0;
    for(j = 0; j < m_param.m_nClusters; j ++) 
	{
		// only work on valid clusters
		if(m_est.m_cluster[j]) 
		{
			s_2 = 0;
			for(i=0; i<m_nDim; i ++) 
			{
				s_2 += m_np[i][j] * (m_data[i]-m_est.m_mu[j])*(m_data[i]-m_est.m_mu[j]);
			}
			s_2 = s_2/m_clusterNP[j];
			m_est.m_sigma[j] = sqrt(s_2);
		}
    }

	return true;
}
/**
* @brief update cluster weights 
* @return true if no error, false otherwise
*/
bool CPrimeEM::updateWeight()
{
	int j = 0;
	for(j = 0; j < m_param.m_nClusters; j ++)
	{
		if(m_est.m_cluster[j])
			m_est.m_weight[j] = m_clusterNP[j]/m_sumNP;
	}

	return true;
}

/**
* @brief validate and adjust seed according to constraints
* @param seed - the seed to be adjusted and returned after adjustment
* @return true if no error, false otherwise
*/
bool CPrimeEM::validate(CEMSeed& seed)
{
	// consistency
	bool consistent = (seed.m_maxMu.size() == seed.m_maxSigma.size());
	consistent &= (seed.m_maxMu.size() == seed.m_minMu.size());
	consistent &= (seed.m_maxMu.size() == seed.m_minSigma.size());
	consistent &= (seed.m_maxMu.size() == seed.m_mu.size());
	consistent &= (seed.m_maxMu.size() == seed.m_sigma.size());
	consistent &= (seed.m_maxMu.size() == seed.m_weight.size());
	consistent &= (seed.m_maxMu.size() == m_param.m_nClusters);
	if(!consistent) return false;

	int i = 0;
	for(i = 0; i < m_param.m_nClusters; i ++)
	{
		if(seed.m_weight[i] > 0)
		{
			m_est.m_cluster[i] = true;

			if(seed.m_mu[i] > seed.m_maxMu[i])
				seed.m_mu[i] = seed.m_maxMu[i];
			else if(seed.m_mu[i] < seed.m_minMu[i])
				seed.m_mu[i] = seed.m_minMu[i];

			if(seed.m_sigma[i] > seed.m_maxSigma[i])
				seed.m_sigma[i] = seed.m_maxSigma[i];
			else if(seed.m_sigma[i] < seed.m_minSigma[i])
				seed.m_sigma[i] = seed.m_minSigma[i];
		}
		else m_est.m_cluster[i] = false;
	}

	// set initial step
	copy(seed.m_mu.begin(),seed.m_mu.end(),m_est.m_mu.begin());
	copy(seed.m_sigma.begin(),seed.m_sigma.end(),m_est.m_sigma.begin());
	copy(seed.m_weight.begin(),seed.m_weight.end(),m_est.m_weight.begin());
	copy(seed.m_minMu.begin(),seed.m_minMu.end(),m_est.m_minMu.begin());
	copy(seed.m_minSigma.begin(),seed.m_minSigma.end(),m_est.m_minSigma.begin());
	copy(seed.m_maxMu.begin(),seed.m_maxMu.end(),m_est.m_maxMu.begin());
	copy(seed.m_maxSigma.begin(),seed.m_maxSigma.end(),m_est.m_maxSigma.begin());

	// the results holder
	m_est.set(m_param.m_nClusters);
	m_preEst.set(m_param.m_nClusters);

	return true;
}
/**
* @brief adjust cluster centers according to constraints
* @return true if no error, false otherwise
*/
bool CPrimeEM::adjustMu()
{
	int i = 0;
	for(i = 0; i < m_param.m_nClusters; i ++)
	{
		if(m_est.m_cluster[i])
		{
			if(m_est.m_mu[i] > m_est.m_maxMu[i])
				m_est.m_mu[i] = m_est.m_maxMu[i];
			else if(m_est.m_mu[i] < m_est.m_minMu[i])
				m_est.m_mu[i] = m_est.m_minMu[i];

			if(m_est.m_sigma[i] > m_est.m_maxSigma[i])
				m_est.m_sigma[i] = m_est.m_maxSigma[i];
			else if(m_est.m_sigma[i] < m_est.m_minSigma[i])
				m_est.m_sigma[i] = m_est.m_minSigma[i];
		}
	}

	return true;
}
/**
* @brief adjust cluster variances according to constraints
* @return true if no error, false otherwise
*/
bool CPrimeEM::adjustSigma()
{
	// if not three clusters just return
	if(m_param.m_nClusters != 3) return true;

    // adjust sigmas
    bool bHet_adjusted = false;
	double reducedSigma = 0;

	double weights = accumulate(m_est.m_weight.begin(),m_est.m_weight.end(),0.0f);
    if(weights > 0) 
	{
		double diff1 = 0, diff2 = 0;
		if(m_est.m_weight[0] != 0 && m_est.m_weight[1] != 0)
			diff1 = fabs(m_est.m_mu[1] - m_est.m_mu[0]);
		if(m_est.m_weight[1] != 0 && m_est.m_weight[2] != 0)
			diff2 = fabs(m_est.m_mu[2] - m_est.m_mu[1]);
		vector<double> sigmaMax = m_pSeed->m_maxSigma;

		// do sigma ajustment
		if(diff1 != 0 && diff1 < 1) 
		{
			reducedSigma = m_param.m_HomReductionCoef * sigmaMax[0] * (1-pow(diff1,m_param.m_SigmaPower));
			if(reducedSigma > sigmaMax[0])
				reducedSigma = sigmaMax[0];
			reducedSigma = reducedSigma*(m_est.m_weight[1]/(m_est.m_weight[0]+m_est.m_weight[1]));
			sigmaMax[0] = sigmaMax[0] - reducedSigma;	
		}

		if(diff2 != 0 && diff2 < 1) 
		{
			reducedSigma = m_param.m_HomReductionCoef * sigmaMax[2] * (1-pow(diff2,m_param.m_SigmaPower));
			if(reducedSigma > sigmaMax[2])
				reducedSigma = sigmaMax[2];
			reducedSigma = reducedSigma*(m_est.m_weight[1]/(m_est.m_weight[1]+m_est.m_weight[2]));
			sigmaMax[2] = sigmaMax[2] - reducedSigma;
		}

		// give hom more weights
		if((diff1!=0 && diff1<diff2) || diff2==0 || diff2 >= 1) 
		{
			if(diff1 != 0 && diff1 < 1) 
			{
				reducedSigma = m_param.m_HetReductionCoef * sigmaMax[1] * (1-pow(diff1,m_param.m_SigmaPower));
				if(reducedSigma > sigmaMax[1])
					reducedSigma = sigmaMax[1];
				reducedSigma = reducedSigma*(m_est.m_weight[0]/(m_est.m_weight[0]+m_est.m_weight[1]));
				sigmaMax[1] = sigmaMax[1] - reducedSigma;
				bHet_adjusted = true;
			}
		}

		if(!bHet_adjusted) 
		{
			if(diff2 != 0 && diff2 < 1) 
			{
				reducedSigma = m_param.m_HetReductionCoef * sigmaMax[1] * (1-pow(diff2,m_param.m_SigmaPower));
				if(reducedSigma > sigmaMax[1])
					reducedSigma = sigmaMax[1];
				reducedSigma = reducedSigma*(m_est.m_weight[2]/(m_est.m_weight[2]+m_est.m_weight[1]));
				sigmaMax[1] = sigmaMax[1] - reducedSigma;			
			}
		}

		// sigmaMax should not be larger than the simgaMax or smaller than sigmaMin seeded
		int j = 0;
		for(j =0; j < m_param.m_nClusters; j ++) 
		{
			if(sigmaMax[j] > m_pSeed->m_maxSigma[j])
				sigmaMax[j] = m_pSeed->m_maxSigma[j];
			if(sigmaMax[j] < m_pSeed->m_minSigma[j])
				sigmaMax[j] = m_pSeed->m_minSigma[j];	
		}

		m_est.m_maxSigma = sigmaMax;
    }

    // check any unusual sigmas
	int j = 0;
	for(j =0; j < m_param.m_nClusters; j ++)
	{
		if(m_est.m_cluster[j]) 
		{
			if(m_est.m_sigma[j] < m_est.m_minSigma[j])
				m_est.m_sigma[j] = m_est.m_minSigma[j];
			else if(m_est.m_sigma[j] > m_est.m_maxSigma[j])
				m_est.m_sigma[j] = m_est.m_maxSigma[j];
		}
    }
    
    // check if any sigmas have gone below or above 2.5/3.0 times of 
    // geometric mean of all sigmas
    double sigmaGeoMean = 1;

	for(j =0; j < m_param.m_nClusters; j ++)
	{
		if(m_est.m_cluster[j]) 
        	sigmaGeoMean = sigmaGeoMean*m_est.m_sigma[j];
    }
	int nValidClusters = accumulate(m_est.m_cluster.begin(),m_est.m_cluster.end(),0);
    if(nValidClusters > 0)
    	sigmaGeoMean = pow(sigmaGeoMean,1/nValidClusters);
    
	for(j =0; j < m_param.m_nClusters; j ++)
	{
		if(m_est.m_cluster[j])
		{
			if(m_est.m_sigma[j] > m_param.m_geo_ratio_multiplier[j]*sigmaGeoMean) 
			{
				m_est.m_sigma[j] = m_param.m_geo_ratio_multiplier[j]*sigmaGeoMean;
				if(m_est.m_sigma[j] > m_est.m_maxSigma[j])
					m_est.m_sigma[j] = m_est.m_maxSigma[j];
			}
	    }
	    else if(m_est.m_sigma[j] < (sigmaGeoMean/m_param.m_geo_ratio_divisor[j])) 
		{
			m_est.m_sigma[j] = (sigmaGeoMean/m_param.m_geo_ratio_divisor[j]);
			if(m_est.m_sigma[j] < m_est.m_minSigma[j])
		   		m_est.m_sigma[j] = m_est.m_minSigma[j];
	    }
    }
	
	return true;
}

/**
* @brief generate variable number of EM seeds based on some prior information
*		this function is designed specifically for 3-clusters scenario, hence only 
*		works if three-clusters is specified, otherwise return false and do nothing
* @return true if no error, false otherwise
*/
bool CPrimeEM::generateEmpiricalEMSeeds()
{
	// this only works for three clusters or fewer
	if(m_param.m_nClusters != 3) return false;

	// clear buffer
	m_seeds.clear();

	vector<double> sigmaMax(3),sigmaMin(3), mu(3), weight(3);
	transform(m_prior.m_sigma.begin(),m_prior.m_sigma.end(),m_prior.m_sigmaSD.begin(),sigmaMax.begin(),plus<double>());
	transform(m_prior.m_sigma.begin(),m_prior.m_sigma.end(),m_prior.m_sigmaSD.begin(),sigmaMin.begin(),minus<double>());

	int a = 0, b = 0, h = 0,k = 0;
	CEMSeed bestSeed;
	double bestFLD= 0, hardshell = 0.25;
	double ratio1 = m_prior.m_muSD[0]/(m_prior.m_muSD[0]+m_prior.m_muSD[1]);
	double ratio2 = m_prior.m_muSD[2]/(m_prior.m_muSD[2]+m_prior.m_muSD[1]);

	vector<double> data(m_nDim);
	vector<double>::iterator uit, lit,it;
	copy(m_data.begin(),m_data.end(),data.begin());
	sort(data.begin(),data.end(),less<double>());
	for(a = 1; a <= 3; a ++)
	{
		for(b = 1; b <= 3; b ++)
		{
			for(h = -2; h <= 2; h ++)
			{
				bool skip = false;
				CEMSeed seed;
				copy(sigmaMax.begin(),sigmaMax.end(),seed.m_maxSigma.begin());
				copy(sigmaMin.begin(),sigmaMin.end(),seed.m_minSigma.begin());
				copy(m_prior.m_sigma.begin(),m_prior.m_sigma.end(),seed.m_sigma.begin());

				double factor[] = {b-1,h,1-a};
				transform(m_prior.m_muSD.begin(),m_prior.m_muSD.end(),&factor[0],mu.begin(),multiplies<double>());
				transform(m_prior.m_mu.begin(),m_prior.m_mu.end(),mu.begin(),seed.m_mu.begin(),plus<double>());
				// mu minimum
				seed.m_minMu[0] = -1.25f;
				double ext = (1-ratio1)*(seed.m_mu[1]-seed.m_mu[0]);
				seed.m_minMu[1] = seed.m_mu[1] - ext + m_param.m_epsilon;
				ext = ratio2*(seed.m_mu[2]-seed.m_mu[1]);
				seed.m_minMu[2] = seed.m_mu[2] - ext + m_param.m_epsilon;
				// mu maximum
				ext = ratio1*(seed.m_mu[1]-seed.m_mu[0]);
				seed.m_maxMu[0] = seed.m_mu[0] + ext - m_param.m_epsilon;
				ext = (1-ratio2)*(seed.m_mu[2]-seed.m_mu[1]);
				seed.m_maxMu[1] = seed.m_mu[1] + ext - m_param.m_epsilon;
				seed.m_maxMu[2] = 1.25f;

				//heck consistencey
				if(seed.m_maxMu[0] > seed.m_minMu[1])
					skip = true;
				if(seed.m_maxMu[1] > seed.m_minMu[2])
					skip = true;

				if(skip) continue;

				int nTotal = 0;
				double percent = 0;
				double FLD = NA, FLD12 = NA, FLD23 = NA, FLD13 = NA;
				vector<double> temp(m_nDim);
				copy(m_data.begin(),m_data.end(),temp.begin());	
				sort(temp.begin(),temp.end(),less<double>());
				for(k = 0; k < m_param.m_nClusters; k ++)
				{
					uit = lower_bound(temp.begin(),temp.end(),seed.m_maxMu[k]);
					lit = upper_bound(temp.begin(),temp.end(),seed.m_minMu[k]);
					int nCount = uit - lit;
					if(nCount > 0)
					{
						seed.m_mu[k] = accumulate(lit,uit,0.0f)/nCount;
						seed.m_weight[k] = nCount;

						vector<double> var(nCount),mean(nCount);
						fill(mean.begin(),mean.end(),seed.m_mu[k]);
						transform(lit,uit,mean.begin(),var.begin(),minus<double>());
						copy(var.begin(),var.end(),mean.begin());
						seed.m_sigma[k] = sqrt(inner_product(var.begin(),var.end(),mean.begin(),0.0f)/double(nCount));
						if(seed.m_sigma[k] < EPSILON)
							seed.m_sigma[k] = sigmaMin[k];
					}
					else 
					{
						seed.m_weight[k] = 0;
						seed.m_mu[k] = m_prior.m_mu[k];;
						seed.m_sigma[k] = m_prior.m_sigma[k];
					}

					nTotal += nCount;
				}

				// did not cover any data points
				if(nTotal <= 0) continue;

				// calculate weights
				for(it=seed.m_weight.begin(); it != seed.m_weight.end(); it++) *it /= double(nTotal);

				// readjust the minmax of the centers
				seed.m_minMu[1] = seed.m_mu[1] - (1-ratio1)*(seed.m_mu[1]-seed.m_mu[0]) + m_param.m_epsilon;
				seed.m_minMu[2] = seed.m_mu[2] - ratio2*(seed.m_mu[2]-seed.m_mu[1]) + m_param.m_epsilon;
				seed.m_maxMu[0] = seed.m_mu[0] + ratio1*(seed.m_mu[1]-seed.m_mu[0]) - m_param.m_epsilon;
				seed.m_maxMu[1] = seed.m_mu[1] + (1-ratio2)*(seed.m_mu[2]-seed.m_mu[1]) - m_param.m_epsilon;

				percent = double(nTotal)/data.size();
				if(seed.m_weight[0] > 0 && seed.m_weight[1] > 0)
				{
					if(seed.m_maxMu[0] > seed.m_minMu[1] || (seed.m_mu[1] - seed.m_mu[0]) < hardshell)
						skip = true;
					FLD12 = (seed.m_mu[1] - seed.m_mu[0])/sqrt(seed.m_sigma[1]*seed.m_sigma[1]+seed.m_sigma[0]*seed.m_sigma[0]);
				}
				if(seed.m_weight[1] > 0 && seed.m_weight[2] > 0)
				{
					if(seed.m_maxMu[1] > seed.m_minMu[2] || (seed.m_mu[2] - seed.m_mu[1]) < hardshell)
						skip = true;
					FLD23 = (seed.m_mu[2] - seed.m_mu[1])/sqrt(seed.m_sigma[1]*seed.m_sigma[1]+seed.m_sigma[2]*seed.m_sigma[2]);
				}
				if(seed.m_weight[0] > 0 && seed.m_weight[2] > 0 && seed.m_weight[1] == 0)
				{
					if(seed.m_mu[2] - seed.m_mu[0] < 3*hardshell)
						skip = true;
					FLD13 = (seed.m_mu[2] - seed.m_mu[0])/sqrt(seed.m_sigma[0]*seed.m_sigma[0]+seed.m_sigma[2]*seed.m_sigma[2]);
				}				

				if(skip) continue;

				// lets look at the stats first to see whether we need more seeds or not
				if(FLD12 != NA && FLD23 != NA) 
					FLD = FLD12 < FLD23? FLD12:FLD23;
				else if(FLD12 != NA)
					FLD = FLD12;
				else if(FLD23 != NA)
					FLD = FLD23;
				else if(FLD13 != NA)
					FLD = FLD13;

				double HET = seed.m_weight[1];
				double HOM = (seed.m_weight[0] < seed.m_weight[2])? seed.m_weight[0] : seed.m_weight[2];
				if(FLD != NA && FLD > bestFLD && percent > 0.9f && HET >= 0.5*HOM && HOM > 0)
				{
					bestSeed = seed;
					bestFLD = FLD;
				}
			
				// collect this seed
				m_seeds.push_back(seed);
			}
		}
	}

	// single seed conditions satisfied
	if(bestFLD > 3)
	{
		m_seeds[0] = bestSeed;
		m_seeds.resize(1);
		return true;
	}

	// add the generic prior seed and two two-cluster generic seeds to cover rare cases
	CEMSeed seed;
	copy(sigmaMax.begin(),sigmaMax.end(),seed.m_maxSigma.begin());
	copy(sigmaMin.begin(),sigmaMin.end(),seed.m_minSigma.begin());
	copy(m_prior.m_mu.begin(),m_prior.m_mu.end(),seed.m_mu.begin());
	copy(m_prior.m_sigma.begin(),m_prior.m_sigma.end(),seed.m_sigma.begin());
	transform(m_prior.m_mu.begin(),m_prior.m_mu.end(),m_prior.m_muSD.begin(),seed.m_maxMu.begin(),plus<double>());
	transform(m_prior.m_mu.begin(),m_prior.m_mu.end(),m_prior.m_muSD.begin(),seed.m_maxMu.begin(),minus<double>());

	// three cluster seed with equal weights
	seed.setWeight(0.33f,0.34f,0.33f);
	m_seeds.push_back(seed);
	// BB and AB cluster seed
	seed.setWeight(0.75f,0.25f,0);
	m_seeds.push_back(seed);
	// AB and AA cluster seed
	seed.setWeight(0,0.25f,0.75f);
	m_seeds.push_back(seed);

	// remove all redundant seeds
	vector<CEMSeed>::iterator sit = m_seeds.begin(), tit;
	while(sit != m_seeds.end())
	{
		for(tit=sit+1; tit != m_seeds.end();)
		{
			if(*sit == *tit) tit = m_seeds.erase(tit);
			else tit ++;
		}
		sit ++;
	}
//	unique(m_seeds.begin(),m_seeds.end());

	return true;
}
