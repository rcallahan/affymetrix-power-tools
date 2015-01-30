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

#include "copynumber/ArrayLevel.h"
//
#include <new>
//
using namespace std;

// row,col givs the size of the ArrayLevel
ArrayLevel::ArrayLevel(int nrow, int ncol)
    {
    m_size = make_pair<int,int>(nrow,ncol);
    m_array.resize(m_size.first*m_size.second,ArrayElement());

    // Assume the next level will not be made and if incorrect allocate it.
    // sz1 and sz2 are computed with the effect of rounding up this array
    // to an even number in size.
    m_next_level = 0;
    if (nrow > 1 && ncol > 1) {
        int sz1 = (nrow - 1)/2 + 1;
        int sz2 = (ncol - 1)/2 + 1;
        m_next_level = new ArrayLevel(sz1,sz2);
        }
    }


// Recursive construction so recursive destruction
ArrayLevel::~ArrayLevel()
    {
    delete m_next_level;
    m_next_level = 0;
    }


pair<unsigned int,unsigned int> ArrayLevel::size() { return m_size; }


// Local phi affects short range smoothing.  Global phi affects long
// range smoothing. Increasing the phis increases their effects.  Increasing
// the phis decreases the local effect of the data.  Convergence terminates
// the method.  An alternative convergence criteria could/should be implemented.
int ArrayLevel::mrf_smooth(double local_phi,double global_phi,double converge)
    {
    // Always smooth the level above first
    if (m_next_level) {
        m_next_level->mrf_smooth(local_phi,global_phi,converge);
        }

    // beat on the solution until one interation at a time.
    int counter = 0;
    double this_converge;
    do {
        counter += 1;
        this_converge = mrf_smooth_iter(local_phi,global_phi);
        } while (this_converge > converge);

    return counter;
    }


// This is the heart of the smoothing.  All the other code is plumbing
// to get the data here.
//
// This algorithm is based on a Markov random field.  You need to
// understand a second order neighborhood in a MRF and how MRFs are
// applied to image analysis to understand this algorithm.  Which is
// beyond the scope of comment lines.

double ArrayLevel::mrf_smooth_iter(double local_phi,double global_phi)
    {
    double converge = 0.0;
    double wt_edge = local_phi*1.0/sqrt(2.0);  // weight in the corners
    double wt_corner = local_phi;              // weight along the probes

    for (int i=0; i<m_size.first; i++) {
        for (int j=0; j<m_size.second; j++) {
            // This is the part where observed data weigh in.
            // Grab a reference to save white space.
            ArrayElement & AE = m_array[i*m_size.second + j];

            // Proportion of the block loaded with observed values.
            double prop = double(AE.count)/double(AE.size);

            // Accumulating a weighted average here.
            double total = AE.mean*prop;
            double sum_wt = prop;

            // If the next level in the linked lists is there use its
            // summary to effect long range smoothing.
            if (m_next_level) {
                total += global_phi*(*m_next_level)(i/2,j/2).mean;
                sum_wt += global_phi;
                }

            // average in the neighbors to fill out the clique.  If you
            // understand an MRF in the context of image analysis this
            // is fairly simple.
            for (int row=i-1; row<i+2; row++) {
                if (row<0 || row==m_size.first) continue;
                for (int col=j-1; col<j+2; col++) {
                    if (col<0 || col==m_size.second) continue;
                    if (row==i && col==j) continue;
                    double wt = (row==i || col==j) ? wt_edge : wt_corner;
                    total += wt*m_array[row*m_size.second + col].mean;
                    sum_wt += wt;
                    }
                }

            // A not so elegant use of a relaxation factor.  Not pretty
            // but it works.
            double old_back = AE.back;
            AE.back = total/sum_wt;
            double diff = AE.back - old_back;
            AE.back += MRF_ITER_RELAX*diff;
            diff += MRF_ITER_RELAX*diff;
            converge += diff*diff;
            }
        }

    return converge;
    }


ArrayElement & ArrayLevel::operator()(int i, int j)
    {
    return m_array[i*m_size.second + j];
    }


// Once the data is loaded in at the head of the linked list, the
// 2x2 blocks of data must be aggregated and passed on to the next
// level.  This is how long range smoothing is achieved.  Think
// wavelets if it helps.
void ArrayLevel::load_next_level()
    {
    // maybe at the end of the list
    if (!m_next_level) return;

    // This is all about adding up information in 2x2 arrays, summarising
    // it by weight and passing it on.  Also making sure that boundaries
    // are not exceeded.
    for (int i=0; i<m_next_level->m_size.first; i++) {
        for (int j=0; j<m_next_level->m_size.second; j++) {
            double mean=0.0;
            int count = 0;
            int size = 0;
            for (int row=2*i; row<2*i + 2; row++) {
                if (row==m_size.first) continue;
                for (int col=2*j; col<2*j + 2; col++) {
                    if (col==m_size.second) continue;
                    ArrayElement & AE = m_array[row*m_size.second + col];
                    size += AE.size;
                    if (!AE.count) continue;
                    mean += AE.count*AE.mean;
                    count += AE.count;
                    }
                }

            if (count) mean /= count;
            (*m_next_level)(i,j) = ArrayElement(mean,count,size);
            }
        }

    m_next_level->load_next_level();
    }
