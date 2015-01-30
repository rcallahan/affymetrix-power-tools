////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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

#include "copynumber/WaveletShrink.h"

#include <cmath>
#include <cstdio>


using std::valarray;
using std::vector;

static void wavelet_sandbox_iter(valarray<double> & data_vec, Wavelet * W,
		vector<double>likelihood_prec, double prior_prec, MODE mode,
		double converge);

static void shrink_mrf1(valarray<double> & vec, double likelihood_prec,
		double prior_prec, double converge);


void wavelet_sandbox(valarray<double> & data_vec, Wavelet * W,
		vector<double>likelihood_prec, double prior_prec, MODE mode,
		double converge, bool outlier, double nu, int maxiter)
    {
	if (!outlier) {
		wavelet_sandbox_iter(data_vec, W, likelihood_prec, prior_prec,
			mode, converge);
		return;
		}

	int N = data_vec.size();
	valarray<double> data(data_vec);
	valarray<double> resids(N);

	double S = 1.0;
	valarray<double> prec_vec(1.0,N);

	int count=0;
	while (count < maxiter) {
		wavelet_sandbox_iter(data, W, likelihood_prec, prior_prec,
			mode, converge);

		resids = data_vec - data;

		S *= (resids*resids*prec_vec).sum()/N;

		prec_vec = (nu + 1.0)/(nu*S + resids*resids);

		data += resids*(prec_vec/prec_vec.sum());

		count++;
		}

	data_vec = data;
    }


void wavelet_sandbox_iter(valarray<double> & data_vec, Wavelet * W,
		vector<double>likelihood_prec, double prior_prec, MODE mode,
		double converge)
    {
	// nothing to do when levels of resolution are exhausted
	if (likelihood_prec.empty()) {
	 return;
	 }

	// length of the data vector
    index_t data_len = (index_t)data_vec.size();

	// length of vector needed to hold the deconstruction / wavelet
	// expansion of the data vector
    index_t dec_len = dwt_buffer_length(data_len,W->dec_len,mode);

	// for holding the expansions
    valarray<double> approximation(0.0, dec_len);
    valarray<double> detail(0.0, dec_len);

	// compute the expansions
    double_dec_a(&data_vec[0], data_len, W, &approximation[0], dec_len, mode);
    double_dec_d(&data_vec[0], data_len, W, &detail[0], dec_len, mode);

	// shrink the details of the expansion
	shrink_mrf1(detail, likelihood_prec[0], prior_prec, converge);

	// send the recursive call a front popped version of the 
	// likelihood weights.  Note that an empty vector will immediately
	// return.
	vector<double> lp_next(likelihood_prec.begin() + 1, likelihood_prec.end());

	// send the next level of approximations to have their details shrunk
	wavelet_sandbox_iter(approximation, W, lp_next, prior_prec, mode, converge);

	// the length of the vector and vector to reconstruct the approximation
    index_t rec_len = reconstruction_buffer_length(dec_len,W->rec_len);
    valarray<double> rec_vec(0.0, rec_len);

	// add the approximation to the reconstruction
    double_upsampling_convolution_valid_sf(&approximation[0],dec_len,
            W->rec_lo_double, W->rec_len,
            &rec_vec[0], rec_len, mode);

	// add the shrunk details to the reconstruction
    double_upsampling_convolution_valid_sf(&detail[0],dec_len,
            W->rec_hi_double, W->rec_len,
            &rec_vec[0], rec_len, mode);

	// write over the input approximation based on shrunk details
	// the copy must be made because rec_vec could be longer than
	// the original data vec due to padding
    for (int i=0; i<data_len; i++) data_vec[i] = rec_vec[i];
    }



// need to comment how shrinkage works
void shrink_mrf1(valarray<double> & vec,
		double likelihood_prec, double prior_prec,
		double converge) {

	int N = vec.size();
	double tau=likelihood_prec;
	double beta=prior_prec;

	double x;

	double prec1 = 1.0/(tau + fabs(beta));
	double prec2 = 1.0/(tau + 2.0*fabs(beta));

	valarray<double>vec1(vec.size());
	valarray<double>old_vec(vec.size());

	vec1 = vec;

	while (1) {
		old_vec = vec1;
		x = beta*vec1[1] + tau*vec[0];
		vec1[0] = x*prec1;


		for (int i=1; i<N-1; i++) {
			x = beta*(vec1[i-1] + vec1[i+1]) + tau*vec[i];
			vec1[i] = x*prec2;
			}

		x = beta*vec1[N-2] + tau*vec[N-1];
		vec1[N-1] = x*prec1;

		double diff_max = abs(vec1 - old_vec).max();
		
		if (diff_max < converge) break;
		}
	vec = vec1;
	}

