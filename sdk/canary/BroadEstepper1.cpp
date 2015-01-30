////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include "canary/BroadEstepper1.h"
//
#include "canary/canary_utils.h"
//
#include <cfloat>
#include <cmath>
//


BroadEstepper1::BroadEstepper1(CanaryOptions& opts,
                               valarray<double> vals_in,
                               CanaryPrior & prior_in,
                               vector<int> & clusters_in)
  : _opts(opts)
  {
  N = vals_in.size();
  G = clusters_in.size();
  _vals.resize(vals_in.size()); _vals = vals_in;
  _prior = prior_in;
  _cvec.resize(G); _cvec = clusters_in;
  _prop.resize(G,1.0/G);
  _mean.resize(G);
  _var.resize(G);

  for (unsigned int k=0; k<G; k++)
    {
    int cpos = _cvec[k];
//  The algorithm converges to a cycle dependent on initial values.  If it
//  converged to a point - as it should - it would make sense to uncomment.
//    _prop[k] = _prior.prop()[cpos];
    _mean[k] = _prior.mean()[cpos];
    _var[k] = _prior.var()[cpos];
    }

  _prob.ReSize(N,G);
  _prob = 0.0;

  update_prob();
  }

BroadEstepper1::BroadEstepper1(CanaryOptions& opts)
  : _opts(opts)
  {
  N=0;
  G=0;
  _prob.ReSize(0,0);
  }

BroadEstepper1::BroadEstepper1(const BroadEstepper1 & BE)
  : _opts(BE._opts)
  {
  N = BE.N;
  G = BE.G;
  _vals.resize(N); _vals = BE._vals;
  _prior = BE._prior;
  _cvec.resize(G); _cvec = BE._cvec;
  _prop.resize(G); _prop = BE._prop;
  _mean.resize(G); _mean = BE._mean;
  _var.resize(G); _var = BE._var;
  //
  _prob.ReSize(N,G);
  _prob = BE._prob;
  }

BroadEstepper1 & BroadEstepper1::operator=(const BroadEstepper1 & BE)
  {
  if (this == &BE) return *this;
  //
  this->_opts=BE._opts;
  //
  this->N = BE.N;
  this->G = BE.G;
  this->_vals.resize(this->N); this->_vals = BE._vals;
  this->_prior = BE._prior;
  this->_cvec.resize(this->G); this->_cvec = BE._cvec;
  this->_prop.resize(this->G); this->_prop = BE._prop;
  this->_mean.resize(this->G); this->_mean = BE._mean;
  this->_var.resize(this->G); this->_var = BE._var;
  this->_prob.ReSize(N,G);
  this->_prob = BE._prob;

  return *this;
  }


void BroadEstepper1::update_broad()
  {
  update_prop_broad();
  update_mean_broad();
  update_var_broad();
  update_prob();
  }


void BroadEstepper1::update_prop_broad()
  {
  double sum_prob = _prob.Sum();

  for (int k=0; k<G; k++)
    {
    double sum_k=0.0;
    for (int i=0; i<N; i++) sum_k += _prob.element(i,k);
    _prop[k] = sum_k/sum_prob;
    }
  }


void BroadEstepper1::update_mean_broad()
  {
  int PPF=_opts.tune_pseudopoint_factor; 
  valarray<double>pmean(G);
  for (int k=0; k<G; k++) pmean[k] = _prior.mean(_cvec[k]);
  int num_pseudopoints = int(N/PPF);
  if (!num_pseudopoints) num_pseudopoints = 1;
  valarray<double>xmean(0.0,G);

  for (int k=0; k<G; k++)
    {
    double kwsum=0.0, wsum=0.0;
    for (int i=0; i<N; i++)
      {
      kwsum += _vals[i]*_prob.element(i,k);
      wsum += _prob.element(i,k);
      if (wsum > 0.0) xmean[k] = kwsum/wsum;
      }
    }

  for (int k=0; k<G; k++)
    {
    _mean[k] = xmean[k]*N*_prop[k] + pmean[k]*num_pseudopoints;
    _mean[k] /= N*_prop[k] + num_pseudopoints;
    }
  }


void BroadEstepper1::update_var_broad()
  {
  valarray<double>xvar(0.0,G);

  for (int k=0; k<G; k++)
    {
    double kwsum=0.0, wsum=0.0;
    for (int i=0; i<N; i++)
      {
      double diff = _vals[i] - _mean[k];
      kwsum += _prob.element(i,k)*diff*diff;
      wsum += _prob.element(i,k);
      }
    if (wsum > 0.0) xvar[k] = kwsum/wsum;
    else xvar[k] = _opts.tune_min_cluster_variance;
    }

  double mean_var = xvar.sum()/xvar.size();
  double reg_weight = _opts.tune_regularize_variance_factor;
  _var = reg_weight*mean_var + (1.0 - reg_weight)*xvar;
  }


void BroadEstepper1::update_prob()
  {
  for (int i=0; i<N; i++)
    {
    valarray<double>xprobs(0.0,G);
    for (int k=0; k<G; k++)
      {
      double diff = _vals[i] - _mean[k];
      xprobs[k] = _prop[k]*exp(-0.5*diff*diff/_var[k])/sqrt(_var[k]);
      }

    double sum_probs = xprobs.sum();

    // just in case all clusters are very unlikely
    if (sum_probs < 1e-10)
      for (int k=0; k<G; k++) _prob.element(i,k) = 1.0/G;

    else for (int k=0; k<G; k++)
        _prob.element(i,k) = xprobs[k]/sum_probs;
    }
  }


double BroadEstepper1::log_likelihood()
  {
  double log_likelihood_sum = 0.0;

  for (int i=0; i<N; i++)
    {
    double this_likelihood = FLT_MIN;
    for (int k=0; k<G; k++)
      {
      double diff = _vals[i] - _mean[k];
	// here is an example of how sensitive the algorithm is
	// uncommenting this line for numerical stability has a
	// large effect on over 1 percent of the CNVs
//      if (diff*diff/_var[k] > 100.0) continue;

      double prob = exp(-0.5*(diff*diff/_var[k]));
      prob /= sqrt(2.0*M_PI*_var[k]);
      prob *= _prop[k];
      this_likelihood += prob;
      }
    log_likelihood_sum += log(this_likelihood);
    }

  return log_likelihood_sum;
  }

//vector<pair<int,double> > BroadEstepper1::confidences()
vector<double> BroadEstepper1::confidences()
  {
  double norm_const = 1.0/sqrt(2.0*M_PI);
//  vector<pair<int,double> >ret_vec;
  vector<double> ret_vec;
  for (unsigned int i=0; i<N; i++)
    {
//    int best_clust = -1;
    double best_prob = 0.0;
    double total_prob = 0.0;
    for (unsigned int k=0; k<_cvec.size(); k++)
      {
      double dev = (_vals[i] - _mean[k])/sqrt(_var[k]);

      // look into how to better handle small probabilities
      double prob = FLT_MIN;
      if (abs(dev) > FLT_MIN) prob = norm_const*exp(-0.5*dev*dev);
      total_prob += prob;
      if (prob > best_prob)
        {
        best_prob = prob;
//        best_clust = _cvec[k];
        }
      }
//      pair<int,float> val(best_clust,best_prob/total_prob);
//      ret_vec.push_back(val);
      ret_vec.push_back(best_prob/total_prob);
    }
  return ret_vec;
  }


// Just assigns the sample point to the cluster of highest probability
vector<int> BroadEstepper1::assignments()
  {
  vector<int> ret_vec;
  for (unsigned int i=0; i<N; i++)
    {
    int best_clust = -1;
    double best_prob = 0.0;

    for (unsigned int k=0; k<_cvec.size(); k++)
      {
      if (_prob.element(i,k) > best_prob)
        {
        best_clust = _cvec[k];
        best_prob = _prob.element(i,k);
        }
      }
      ret_vec.push_back(best_clust);
    }
  return ret_vec;
  }

