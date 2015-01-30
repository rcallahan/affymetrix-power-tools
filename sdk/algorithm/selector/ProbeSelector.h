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
 * @file   ProbeSelector.h
 * @author Earl Hubbell
 * @date   Wed Dec 6 2006
 * 
 */

#ifndef _PROBESELECTOR_H 
#define _PROBESELECTOR_H

//
#include <iostream>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
//
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>
//

using namespace std;

// header: parameters, functions, etc for generating covariate normalization
// predictors and making predictions for correction factors

double LogisticRegression(
	vector<double> &Summaries,
	vector<double> &Predictor,
	const vector<int> &TargetGenotypes,
	const vector< vector<double> > &Covariates,
				  vector<double> PriorMeans,
			  vector<double> PriorWeights,
	int loopMax,
	double Converged);

#endif /* _PROBESELECTOR_H_ */
