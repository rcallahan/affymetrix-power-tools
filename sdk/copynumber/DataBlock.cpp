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

#include "copynumber/DataBlock.h"
//

using namespace std;

// The constructor creates a block of data a multiple of the smaller
// block sizes that are summarised for ArrayLevel
DataBlock::DataBlock(unsigned int nrow, unsigned int ncol,
        int brows, int bcols)
    {
    nrow = brows*((nrow - 1)/brows + 1);
    ncol = bcols*((ncol - 1)/bcols + 1);

    // Here is where the numeric data are stored.
    m_data.resize(nrow*ncol,0.0);
    m_size = make_pair<unsigned int,unsigned int>(nrow,ncol);

    // Assume the data are missing until seen
    m_missing.resize(nrow*ncol,true);

    // The background is estimated in blocks to save computation time.
    m_block_size = make_pair<unsigned int,unsigned int>(brows,bcols);

    int sz1 = (nrow - 1)/brows + 1;
    int sz2 = (ncol - 1)/bcols + 1;
    m_nblocks = make_pair<unsigned int,unsigned int>(sz1,sz2);

    // This is where the background will be estimated.
    m_array = new ArrayLevel(sz1,sz2);
    }


DataBlock::~DataBlock()
    {
    delete m_array;
    m_array = 0;
    }


// After construction, the data are added
void DataBlock::set_value(const int i, const int j, const double value)
    {
    m_data[i*m_size.second + j] = value;
    m_missing[i*m_size.second + j] = false;
    }


// For the user to handle by missing values.
bool DataBlock::missing(const int i, const int j) const
    {
    return m_missing[i*m_size.second + j];
    }


double DataBlock::get_value(const int i, const int j) const
    {
    return m_data[i*m_size.second + j];
    }


// The background is pulled out of the head of the ArrayLevel
double DataBlock::background(const int i, const int j) const
    {
    return (*m_array)(i/m_block_size.first,j/m_block_size.second).back;
    }


// Observed minus background.
double DataBlock::residual(const int i, const int j) const
    {
    if (m_missing[i*m_size.second + j]) return 0;
    double value = m_data[i*m_size.second + j];
    return value - (*m_array)(i/m_block_size.first,j/m_block_size.second).back;
    }


// After all the data are loaded by set_value() summaries are propagated
// through the ArrayLevel list.
void DataBlock::load_array_level()
    {
    int block_size = m_block_size.first*m_block_size.second;

    for (unsigned int i=0; i<m_nblocks.first; i++) {
        for (unsigned int j=0; j<m_nblocks.second; j++) {
            int row_start = i*m_block_size.first;
            int row_end = (i+1)*m_block_size.first;
            double mean = 0.0;
            int count = 0;
            for (int row=row_start; row<row_end; row++)
                {
                int col_start = j*m_block_size.second;
                int col_end = (j+1)*m_block_size.second;
                for (int col=col_start; col<col_end; col++)
                    {
                    if (m_missing[row*m_size.second + col]) continue;
                    mean += m_data[row*m_size.second + col];
                    count++;
                    }
                }
            if (count) mean /= count;
            (*m_array)(i,j) = ArrayElement(mean,count,block_size);
            }
        }

    m_array->load_next_level();
    }


// Tell the ArrayLevel list to estimate the background
void DataBlock::mrf_smooth(const double local_phi, const double global_phi,
        const double converge)
    {
    m_array->mrf_smooth(local_phi, global_phi, converge);
    }


// This may not be needed.
ArrayElement DataBlock::get_block(const int i, const int j)
    {
    return (*m_array)(i,j);
    }
