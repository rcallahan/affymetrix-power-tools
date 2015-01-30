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

#include "canary/CanaryWithPrior.h"
//
#include "canary/CanaryModel.h"
#include "canary/InterferenceModel.h"
#include "canary/canary_utils.h"
//
#include <cmath>
//

CanaryWithPrior::CanaryWithPrior(CanaryOptions& opts, 
                                 valarray<double> & vals_in,
                                 CanaryPrior & prior_in)
  : _opts(opts)
  {
  _vals.resize(vals_in.size()); _vals = vals_in;
  _prior = prior_in;
  valarray<double>vars = _prior.var();
  for (int k=0; k<vars.size(); k++)
    if (vars[k] < _opts.tune_min_cluster_variance)
      vars[k] = _opts.tune_min_cluster_variance;
  _prior.var(vars);

  if ((prior_in.prop(0) + prior_in.prop(1)) == 0.0)
    {
    _canary_model.erase("model0");
	_canary_model.erase("model01");
	_canary_model.erase("model12");
	_canary_model.erase("model012");
	_canary_model.erase("model123");
	_canary_model.erase("model0123");
	_canary_model.erase("model1234");
	_canary_model.erase("model01234");
	}

  if ((prior_in.prop(3) + prior_in.prop(4)) == 0.0)
    {
    _canary_model.erase("model23");
	_canary_model.erase("model34");
	_canary_model.erase("model234");
	_canary_model.erase("model123");
	_canary_model.erase("model0123");
	_canary_model.erase("model1234");
	_canary_model.erase("model01234");
	}

//  leave these print statements in for future debugging 
//  printf("\n");
//  printf("model      loglik  bic     close     over   af   hwe                  fit\n");
  for (unsigned long k=0; k<_canary_model.size(); k++)
    {
    string model_name = _canary_model(k);
    vector<int>cvec = _canary_model[model_name];
    BroadEstepper1 BE = gmm_withprior(_vals,cvec);
    BE_map.insert(std::pair<std::string,BroadEstepper1>(model_name,BE));

    CanaryPenalty CP = assay_fit(BE,cvec);
//	printf("%-8s %8.1f %5.1f %8.1f %8.1f %5.3f %20.17f %7.1f\n",
//		model_name.c_str(), CP.loglikelihood(), CP.bic(), CP.closeness(),
//		CP.overlap(), CP.af(), CP.hwe(), CP.fit());
	
    CP_map[model_name] = CP;
    } 
//  printf("\n");
//  printf("\n");
  fflush(stdout);
  }

string CanaryWithPrior::best_fit()
  {
  vector<pair<string,double> >fits;

  for (unsigned int k=0; k<_canary_model.size(); k++)
    {
    string model = _canary_model(k);
    double fit = CP_map[model].fit();
    fits.push_back(pair<string,double>(model,fit));
    }

  string best_model = fits[0].first;
  double best_fit = fits[0].second;

  for (unsigned int k=1; k<fits.size(); k++)
    {
    if (fits[k].second > best_fit)
      {
      best_model=fits[k].first;
      best_fit=fits[k].second;
      }
    }

  return best_model;
  }
  


BroadEstepper1 CanaryWithPrior::gmm_withprior(valarray<double> vals,
      vector<int> cluster_vec)
  {
  int MIN_ITER = 10;
  int MAX_ITER = 100;
  double LOGLIK_THRESHOLD = 0.01;

  BroadEstepper1 BE(_opts,vals,_prior,cluster_vec);

  double old_loglik=0.0;
  int iterations = 0;
  while(1)
    {
    iterations++;
    double loglik = BE.log_likelihood();

    BE.update_broad();

    double delta_loglik = fabs(loglik - old_loglik);
    old_loglik = loglik;
	

    if (iterations < MIN_ITER) continue;
    if (iterations == MAX_ITER) break;
    if (delta_loglik < LOGLIK_THRESHOLD) break;
    }

  return BE;
  }


CanaryPenalty CanaryWithPrior::assay_fit(BroadEstepper1 BE,
                                         vector<int> cluster_vec)
  {
  unsigned long N = BE.vals().size();
  unsigned long G = cluster_vec.size();

  double loglik = BE.log_likelihood();
  double bic = 2.0*loglik - G*log((double)N);
  double reward_closeness=0.0;

  for (unsigned int k=0; k<G; k++)
    {
    double closeness = abs(BE.mean(k) - _prior.mean(cluster_vec[k]));
    if (closeness < _opts.tune_TOL) closeness = _opts.tune_TOL;
    reward_closeness += BE.prop(k)/closeness;
    }


  // OPTIONS
  reward_closeness *= 0.1*G*abs(bic);

  valarray<double>frequencies = BE.prop();
  valarray<double>delta_frequencies = _prior.prop();

  for (unsigned int k=0; k<G; k++)
    delta_frequencies[cluster_vec[k]] -= frequencies[k];

  double af = _opts.tune_af_weight*(delta_frequencies*delta_frequencies).sum();

  double hwe = hardy_weinberg(BE.cvec(),(double)BE.vals().size()*BE.prop());

  double hwe_penalty=0.0;
  if (hwe < _opts.tune_hwe_tol) hwe_penalty = -log10(hwe);
  if (hwe < _opts.tune_hwe_tol2) hwe_penalty = -50*log10(hwe);

  double overlap = InterferenceModel(_opts,BE.mean(),BE.var()).interference();
  overlap *= 1000;

  double fit = 5.0*bic - abs(af*bic);
  fit +=  reward_closeness/(overlap + af + 1);
  fit -= overlap + 30*hwe_penalty;

  return CanaryPenalty(loglik,bic,reward_closeness,af,fit,hwe,overlap);
  }
