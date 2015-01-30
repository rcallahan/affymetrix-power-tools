////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
#include "copynumber/BinarySegmenter.h"


#include <cmath>
#include <cstdio>
#include <vector>

using std::pair;
using std::valarray;
using std::vector;

static void noprior_prob(std::vector<BinaryCNode>::iterator bcn_iter,
		 const double mu, const double S);

static void noprior_prob_normal(vector<BinaryCNode>::iterator bcn_iter,
		const pair<double,double> mu, const pair<double,double> var,
		const double bias);

static void carry_forward(std::vector<BinaryCNode>::iterator bcn_prev,
	std::vector<BinaryCNode>::iterator bcn, double prior_prob);

typedef std::vector<BinaryCNode>::iterator BCN_iter;


void binary_segmenter(std::vector<BinaryCNode>::iterator bcn_start,
		std::vector<BinaryCNode>::iterator bcn_end,
		std::vector<BinaryCNodeParam> bcn_params,
		int level)
	{
	if (bcn_params.empty()) {
		for (BCN_iter Bit = bcn_start; Bit != bcn_end; Bit++) {
			Bit->prob = 1.0;
			Bit->level = level;
			}
		return;
		}

	double mu = bcn_params[0].mu;
	double S = bcn_params[0].S;
	double tprob = bcn_params[0].tprob;

	for (BCN_iter Bit = bcn_start; Bit != bcn_end; Bit++) {
		noprior_prob(Bit,mu,S);
		}

	bcn_start->path_back = std::make_pair(0,1);

	for (BCN_iter Bit = bcn_start+1; Bit != bcn_end; Bit++) {
		carry_forward(Bit-1,Bit,tprob);
		}

	std::deque<int> binary_path;

	BCN_iter iter = bcn_end - 1;
	int back_to = iter->prob > 0.5 ? 0 : 1;

	binary_path.push_front(back_to);

	 while (iter != bcn_start) {
		if (!back_to) back_to = iter->path_back.first;
		else back_to = iter->path_back.second;
		binary_path.push_front(back_to);
		iter--;
		}


	iter = bcn_start;
	int idx = 0;

	while (iter != bcn_end) {
		if (!binary_path[idx]) {
			iter->level = level;
			idx++;
			iter++;
			continue;
			}

		BCN_iter seg_start = iter;

		while (iter != bcn_end) {
			if (!binary_path[idx]) break;
			idx++;
			iter++;
			}

		BCN_iter seg_end = iter;

		std::vector<BinaryCNodeParam> bcn_params_short;
		bcn_params_short.assign(bcn_params.begin()+1,bcn_params.end());

		binary_segmenter(seg_start,seg_end,bcn_params_short,level+1);
		}
	}



void binary_segmenter_normal(BCN_iter bcn_start, BCN_iter bcn_end,
		std::vector<BinaryCNodeParamNormal> bcn_params,
		int level)
	{
	if (bcn_params.empty()) {
		for (BCN_iter Bit = bcn_start; Bit != bcn_end; Bit++) {
			Bit->prob = 1.0;
			Bit->level = level;
			}
		return;
		}

	pair<double,double> mu = bcn_params[0].mu;
	pair<double,double> var = bcn_params[0].var;
	double tprob = bcn_params[0].tprob;
	double bias = bcn_params[0].bias;

	for (BCN_iter Bit = bcn_start; Bit != bcn_end; Bit++) {
		noprior_prob_normal(Bit,mu,var,bias);
		}


	bcn_start->path_back = std::make_pair(0,1);

	for (BCN_iter Bit = bcn_start+1; Bit != bcn_end; Bit++) {
		carry_forward(Bit-1,Bit,tprob);
		}

	std::deque<int> binary_path;

	BCN_iter iter = bcn_end - 1;
	int back_to = iter->prob > 0.5 ? 0 : 1;

	binary_path.push_front(back_to);

	 while (iter != bcn_start) {
		if (!back_to) back_to = iter->path_back.first;
		else back_to = iter->path_back.second;
		binary_path.push_front(back_to);
		iter--;
		}


	iter = bcn_start;
	int idx = 0;

	while (iter != bcn_end) {
		if (!binary_path[idx]) {
			iter->level = level;
			idx++;
			iter++;
			continue;
			}

		BCN_iter seg_start = iter;

		while (iter != bcn_end) {
			if (!binary_path[idx]) break;
			idx++;
			iter++;
			}

		BCN_iter seg_end = iter;

		std::vector<BinaryCNodeParamNormal> bcn_params_short;
		bcn_params_short.assign(bcn_params.begin()+1,bcn_params.end());

		binary_segmenter_normal(seg_start,seg_end,bcn_params_short,level+1);
		}
	}



//---------------------------------------------------------------------


void noprior_prob(std::vector<BinaryCNode>::iterator bcn_iter,
		 const double mu, const double S)
	{
	double x = (bcn_iter->value - mu)/S;
	bcn_iter->prob = 1.0 - 1.0/(1.0 + exp(-x));
	}


void noprior_prob_normal(std::vector<BinaryCNode>::iterator bcn_iter,
		const pair<double,double> mu, const pair<double,double> var,
		const double bias)
	{
	double x = (bcn_iter->value - mu.first);
	double p1 = bias*exp(-0.5*(x*x)/var.first)/sqrt(var.first);

	x = (bcn_iter->value - mu.second);
	double p2 = exp(-0.5*(x*x)/var.second)/sqrt(var.second);

	bcn_iter->prob = p1/(p1 + p2);
	}


void carry_forward(std::vector<BinaryCNode>::iterator bcn_prev,
	std::vector<BinaryCNode>::iterator bcn, double prior_prob)
	{
	double prev_prob = bcn_prev->prob;

	double p0 = bcn->prob;
	p0 *= prev_prob*prior_prob +  (1.0 - prev_prob)*(1.0 - prior_prob);

	double p1 = 1.0 - bcn->prob;
	p1 *= prev_prob*(1.0 - prior_prob) + (1.0 - prev_prob)*prior_prob;

	bcn->prob = p0/(p0 + p1);

	bcn->path_back = std::make_pair(0,0);

	if (prev_prob*prior_prob < (1.0 - prev_prob)*(1.0 - prior_prob)) {
		bcn->path_back.first = 1;
		}

	if (prev_prob*(1.0 - prior_prob) < (1.0 - prev_prob)*prior_prob) {
		bcn->path_back.second = 1;
		}
	}


std::vector<BinaryCNodeParam> make_bcn_params(std::valarray<double> mu,
        double prob0, double tprob)
    {
    std::vector<BinaryCNodeParam> bcn_params;

    for (int i=1; i<mu.size(); i++) {
        double m1 = exp(mu[i-1]);
        double m2 = exp(mu[i]);
        double m = log(0.5) + log(m1 + m2);
        double d = m - mu[i-1];
        double S = d/log(prob0/(1.0 - prob0));
        bcn_params.push_back(BinaryCNodeParam(m,S,tprob));

		printf("-- make_bcn_params m,S,tprob = %f,%f,%f\n", m,S,tprob);
        }
    return bcn_params;
	}



vector<BinaryCNodeParamNormal> make_bcn_params_normal( valarray<double> mu,
		valarray<double> var, valarray<double> tprob, valarray<double> bias)
    {
    vector<BinaryCNodeParamNormal> bcn_params;

	for (int i=1; i<mu.size(); i++) {
		std::pair<double,double>mu1 =
			std::make_pair<double,double>(mu[i-1],mu[i]);
		std::pair<double,double>var1 =
			std::make_pair<double,double>(var[i-1],var[i]);

		bcn_params.push_back(BinaryCNodeParamNormal(mu1,var1,
			tprob[i-1],bias[i-1]));

		}
    return bcn_params;
	}

