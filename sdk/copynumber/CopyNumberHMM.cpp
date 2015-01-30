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
#include "CopyNumberHMM.h"

#include <vector>

using namespace std;


// Need to throw exceptions when real users come along.
CopyNumberHMM::CopyNumberHMM(
        const TransitionMatrix TransMat_in)
        : TransMat(TransMat_in) { }


// Write over the top of any old data with new data.
void CopyNumberHMM::add_data(
        const valarray<double> & data_in,
        const valarray<double> & mu_in,
        const valarray<double> & var_in)
    {
    _data.resize(data_in.size());
    _data = data_in;
    _mu.resize(mu_in.size());
    _mu = mu_in;
    _prec.resize(this->size());
    _prec = 1.0/var_in;
    }


valarray<int> CopyNumberHMM::solve_map(bool constantMeansFlag)
    {
    const unsigned int N = _data.size();  // number of nodes in the HMM
    const unsigned int P = this->size();  // number of states per node
    vector<CopyNumberNode>cnn;
    cnn.reserve(N);

    // This spares the need to allocate copies of many little slices
    valarray<double>mu_marker(P);
    mu_marker = _mu[slice(0,P,1)];

    cnn.push_back(CopyNumberNode(TransMat,_data[0],mu_marker,_prec));

    for (unsigned int t=1; t<_data.size(); t++) {
        if (!constantMeansFlag) mu_marker = _mu[slice(t*P,P,1)];
            cnn.push_back(CopyNumberNode(cnn[t-1],TransMat,_data[t],
            mu_marker, _prec));
        }

    valarray<int>copy_number(N);
    copy_number[N-1] = cnn[N-1].backward();

    for (unsigned int t=1; t<N; t++)
        copy_number[N-t-1] = cnn[N-t].backward(copy_number[N-t]);

    return copy_number;
    }

