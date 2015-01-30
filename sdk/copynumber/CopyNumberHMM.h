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
#ifndef COPYNUMBERHMM_H_
#define COPYNUMBERHMM_H_

#include "TransitionMatrix.h"
#include "CopyNumberNode.h"

#include <valarray>

class CopyNumberHMM {

public:
CopyNumberHMM(const TransitionMatrix TransMat_in);

void add_data(const std::valarray<double> & data_in,
    const std::valarray<double> & mu_in,
    const std::valarray<double> & var_in);

// Issues the instructions to solve.  Probably used for each chromosome.
std::valarray<int> solve_map(bool constantMeansFlag = true);

// How many states there are.  Maybe a better name could be used.
unsigned int size() const { return TransMat.size(); }

// Return the probabilities for any of the states given.
std::valarray<double> state_prob(const std::valarray<int> states);


private:
unsigned int nstate;
TransitionMatrix TransMat;
std::valarray<double> _prec;  // use precisions instead of variances
std::valarray<double> _mu;    // means
std::valarray<double> _data;  // observed values
};

#endif // COPYNUMBERHMM_H_
