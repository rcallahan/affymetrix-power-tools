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

#ifndef _canary_utils_H
#define _canary_utils_H

#include "canary/BroadEstepper1.h"
#include "canary/CanaryPrior.h"
//
#include <cstring>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <valarray>
#include <vector>
//

pair<list<string>,map<string,valarray<double> > >
  load_broad_intensity_map(string ifname);

map<string,CanaryPrior> load_broad_prior_map(string ifname);

double hardy_weinberg(vector<int> cluster_vec,valarray<double> props);
double hardy_weinberg2(double Naa, double Nab, double Nbb);
double hardy_weinberg3(double N0,double N1,double N2,double N3,double N4);
double hardy_weinberg3_bad(double N0,double N1,double N2,double N3,double N4);

BroadEstepper1 directional_impute(CanaryOptions& opts,
                                  pair<string,BroadEstepper1> model_pair,
                                  CanaryPrior CP);

BroadEstepper1 directional_impute_01234(CanaryOptions& opts,
                                        pair<string,BroadEstepper1> model_pair,
                                        CanaryPrior CP);

valarray<double> fill_prop(CanaryOptions& opts,
                           valarray<double>prop_in,vector<int>cvec_in,
vector<int>cvec_forced);

valarray<double> median_polish(Matrix & M, int niter=1000);
double compute_median(valarray<double> valvec);

string swap_suffix(string name1, string new_suffix);

#endif
