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
#ifndef DATABLOCK_H
#define DATABLOCK_H

#include "copynumber/ArrayElement.h"
#include "copynumber/ArrayLevel.h"
//
#include <utility>
#include <valarray>
#include <vector>
//



class DataBlock {
public:
    DataBlock(unsigned int nrow, unsigned int ncol, int brows, int bcols);
    ~DataBlock();
    bool missing(int i, int j) const;
    void set_value(const int i, const int j, const double value);
    double get_value(const int i, const int j) const;
    double residual(const int i, const int j) const;
    double background(const int i, const int j) const;
    std::pair<unsigned int, unsigned int> size() { return m_size; }
    std::pair<unsigned int, unsigned int> block_size() { return m_block_size; }
    ArrayElement get_block(const int i, const int j);
    void load_array_level();
    void mrf_smooth(const double local_phi, const double global_phi,
            const double converge);


private:
    std::pair<unsigned int, unsigned int> m_size;
    std::pair<unsigned int, unsigned int> m_block_size;
    std::pair<unsigned int, unsigned int> m_nblocks;
    std::valarray<double>m_data;
    std::valarray<bool>m_missing;
    ArrayLevel * m_array;
};


#endif
