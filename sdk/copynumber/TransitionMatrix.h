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
#ifndef TRANSITIONMATRIX_H_
#define TRANSITIONMATRIX_H_

#include <valarray>

class TransitionMatrix {
public:

TransitionMatrix(const unsigned int N, const double val) : _size(N)
    {
    double off_diag = (1.0 - val)/(_size - 1);
    data.resize(_size*_size,off_diag);
    for (unsigned int j=0; j<_size; j++) {
        data[j*_size + j] = val;
        }
    }


TransitionMatrix(const TransitionMatrix & T)
    {
    _size = T._size;
    data.resize(T.data.size());
    data = T.data;
    }


TransitionMatrix & operator=(const TransitionMatrix & T)
    {
    if (this == &T) return *this;

    _size = T._size;
    data.resize(T.data.size());
    data = T.data;
    return *this;
    }


double operator()(const unsigned int k1, const unsigned int k2) const;
double operator()(const unsigned int k) const;
unsigned long size() const { return _size; }

private:
std::valarray<double> data;
unsigned long _size;
};


inline double TransitionMatrix::
    operator()(const unsigned int k1, const unsigned int k2) const
    {
    return data[k1*_size + k2];
    }

inline double TransitionMatrix::operator()(const unsigned int k) const
    {
    return data[k*_size + k];
    }

#endif // TRANSITIONMATRIX_H_
