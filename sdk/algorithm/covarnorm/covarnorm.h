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
 * @file   covarnorm.h
 * @author Earl Hubbell
 * @date   Fri 14 Sept 2006
 */

#ifndef _COVARNORM_H_
#define _COVARNORM_H_

//
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>
//
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
//

// Dont use this in a ".h" file, as it will alter the namespace of the including file.
// using namespace std;

// header: parameters, functions, etc for generating covariate normalization
// predictors and making predictions for correction factors

void td(const std::vector<double> &v);

/** this class holds information required to make predictions on new snps for normalization */
class NormalizationPredictor{
public:
	// information required to make predictions on new SNPs in a given experiment
	// vector for aa,ab,bb genotypes
	// sigma for each genotype
	std::vector<double> fitAA; ///< prediction function for AA
	std::vector<double> fitAB; ///< prediction function for AB
	std::vector<double> fitBB; ///< prediction function for BB
	// detect genotype by distance scaled by sigma
	std::vector<double> sigma; ///< sigma within each genotype
	// offset for covariates
	// covariates are centered for predictor
	std::vector<double> CenterCovariates; ///< mean values for covariates
	// where are the true centers supposed to be
	std::vector<double> TargetCenter; ///< where are the true centers.

	void MakeCov(std::vector<double> &cov,const double Strength);
	double MakePrediction(const double Contrast,const std::vector<double> &cov);
	void Copy(NormalizationPredictor &Source);
	void Dump(std::ostream &out);
};

/** Holds some useful values in relation to each other */
class NormalizationFrame{
public:
	std::vector<double> Contrast; ///< want to predict this value
	// covariates
	std::vector<double> Strength; ///< using this value

	// vector of contrast
	// matrix of covar
	NormalizationPredictor ThisFit; ///< functions for predicting
	// fill in contrasts
	// fill in covariates
	// build normalization predictor
	void FillInCovariates(const std::vector<double> &S);
	void FillInContrast(const std::vector<double> &v);	
	void FitEM();
};

/** this is where the magic works to fit a cubic regression */
void
FitWeightedCubic(const std::vector<double> &contrast,
                 const std::vector<double> &strength,
                 const std::vector<double> &weights,
                 std::vector<double> &Predictor,
                 std::vector<double> &Predicted); 

#endif /* _COVARNORM_H_ */
