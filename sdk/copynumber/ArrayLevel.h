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
#ifndef ARRAYLEVEL_H
#define ARRAYLEVEL_H

#include "copynumber/ArrayElement.h"
//
#include <cmath>
#include <utility>
#include <vector>
//

#define MRF_ITER_RELAX 0.45


class ArrayLevel {
public:
    ArrayLevel();
    ArrayLevel(int nrow, int ncol);
    ~ArrayLevel();
    std::pair<unsigned int, unsigned int> size();
    ArrayElement & operator()(int i, int j);
    int mrf_smooth(double local_phi, double global_phi, double converge);
    double mrf_smooth_iter(double local_phi, double global_phi);
    void load_next_level();

private:
    ArrayLevel * m_next_level;
    std::pair<unsigned int,unsigned int> m_size;
    std::vector<ArrayElement> m_array;
};


#endif
