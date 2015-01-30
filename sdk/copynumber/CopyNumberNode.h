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
#ifndef COPYNUMBERNODE_H_
#define COPYNUMBERNODE_H_


#include <float.h>
#include "TransitionMatrix.h"
#include <valarray>


class CopyNumberNode {
public:

CopyNumberNode() { }

CopyNumberNode(const TransitionMatrix & TM, const double obs_val,
        const std::valarray<double>mu, const std::valarray<double>prec);

CopyNumberNode(const CopyNumberNode & CNN, const TransitionMatrix & TM,
        const double obs_val, const std::valarray<double>mu,
        const std::valarray<double>prec)
    {
    fill_in(CNN,TM,obs_val,mu,prec);
    }

// In case someone needs to know how many of X can be got from *this.
unsigned long size() const { return prob_forward.size(); }

double forward(const unsigned int k) const { return prob_forward[k]; }
int backward(const int k) const { return path_backward[k]; }

int backward() const
    {
    double mval=0.0;
    int mpos=0;
    for (int k=0; k<this->size(); k++)
        if (prob_forward[k] > mval) { mval=prob_forward[k]; mpos=k; }
    return path_backward[mpos];
    }


// Assignment operator
CopyNumberNode & operator=(const CopyNumberNode & CNN)
    {
    if (this == &CNN) return *this;

    prob_forward.resize(CNN.prob_forward.size());
    prob_forward = CNN.prob_forward;

    path_backward.resize(CNN.path_backward.size());
    path_backward = CNN.path_backward;

    return *this;
    }


// Copy constructor
CopyNumberNode(const CopyNumberNode & CNN)
    {
    prob_forward.resize(CNN.prob_forward.size());
    prob_forward = CNN.prob_forward;

    path_backward.resize(CNN.path_backward.size());
    path_backward = CNN.path_backward;
    }


private:
// The relative probabilities of MAP estimates leading to each copy
std::valarray<double>prob_forward;

// Where to look in the node before
std::valarray<int>path_backward;

void fill_in(const CopyNumberNode & CNN, const TransitionMatrix & TM,
    const double obs_val, const std::valarray<double>mu,
    const std::valarray<double>prec);

};

inline CopyNumberNode::CopyNumberNode(const TransitionMatrix & TM,
        const double obs_val, const std::valarray<double>mu,
        const std::valarray<double>prec)
    {
    prob_forward.resize(TM.size());
    for (unsigned int k=0; k<TM.size(); k++)
        {
        double diff = obs_val - mu[k];
        // Transition probability times likelihood.
        prob_forward[k] = TM(k)*sqrt(prec[k])*exp(-0.5*prec[k]*diff*diff);
        }

    // For numerical stability rescale the probabilities
    prob_forward /= prob_forward.sum()/prob_forward.size();

    // Not that it matters because this is the first node but the path
    // back is consistent with the use of the transition matrix.
    path_backward.resize(this->size());
    for (int k=0; k<this->size(); k++) path_backward = k;
    }



// Use this to fill in all constructors that conform.  Any constructor
// used after the first node is in place will probably use this one.
inline void CopyNumberNode::fill_in(const CopyNumberNode & CNN,
        const TransitionMatrix & TM, const double obs_val,
        const std::valarray<double>mu, const std::valarray<double>prec)
    {
    prob_forward.resize(TM.size());
    path_backward.resize(TM.size());

    // For each level find the most probable path to it.  The most probable
    // path is the product of previous probabilities and their respective
    // transition probabilities.
    // mval keeps track of the highest probability path at level k
    // mpos keeps track of the previous level for this level k
    for (unsigned int k=0; k<this->size(); k++)
        {
        int mpos=0;
        double mval=0.0;
        for (int j=0; j<this->size(); j++)
            {
            // previous probability times transition
            double this_val = TM(j,k)*CNN.forward(j);
            if (this_val > mval) { mval=this_val; mpos=j; }
            }

        double diff = obs_val - mu[k];
        // multiply by the likelihood or emmision in HMM terminology
        prob_forward[k] = mval*sqrt(prec[k])*exp(-0.5*prec[k]*diff*diff);

        // the highest probability path to backtrack
        path_backward[k] = mpos;
        }

    // if the likelihood says nothing due to an outlier, carry the previous
    // probe forward so the old probes acts as a proxy for this bad data probe.
    if (prob_forward.sum() < FLT_MIN)
        {
        for (int k=0; k<this->size(); k++)
            {
            prob_forward[k] = CNN.forward(k);
            path_backward[k] = k;
            }
        }

    // rescale the probabilities for numerical stability
    prob_forward /= prob_forward.sum()/prob_forward.size();
    }

#endif // COPYNUMBERNODE_H_
