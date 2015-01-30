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
#ifndef BINARYCOPYNUMBERNODE_H_
#define BINARYCOPYNUMBERNODE_H_

#include <cmath>

#include <deque>
#include <utility>
#include <valarray>
#include <vector>

typedef struct BinaryCNode_t {
	BinaryCNode_t(const double val): value(val)
		{
		prob = 0.0;
		level = 0;
		path_back = std::make_pair(0,1);
		}

	double value;
	double prob;
	int level;
	std::pair<int,int> path_back;
	} BinaryCNode;



typedef struct BinaryCNodeParam_t {
	BinaryCNodeParam_t(double mu_in, double S_in, double tprob_in):
		mu(mu_in), S(S_in), tprob(tprob_in) { }
	double mu;
	double S;
	double tprob;
	} BinaryCNodeParam;


typedef struct BinaryCNodeParamNormal_t {
	BinaryCNodeParamNormal_t(std::pair<double,double>mu_in,
		std::pair<double,double>var_in, double tprob_in, double bias_in):
			mu(mu_in), var(var_in), tprob(tprob_in), bias(bias_in) {
		}
	std::pair<double,double> mu;
	std::pair<double,double> var;
	double tprob;
	double bias;
	} BinaryCNodeParamNormal;


#endif
