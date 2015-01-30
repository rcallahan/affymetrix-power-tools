////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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

/*
 * FILE MultiDimensionalArray.h
 */

#include "broadutil/MultiDimensionalArray.h"
//
#include <cassert>
#include <vector>
//

using namespace std;

MultiDimensionalArraySpec::MultiDimensionalArraySpec(const vector<size_t> &dimensions):
    m_dimensions(dimensions)
{
  assert(m_dimensions.size() > 1);
  m_numElements = 1;
  for (vector<size_t>::const_iterator it = m_dimensions.begin();
       it != m_dimensions.end(); ++it) {
    m_numElements *= *it;
  }
}

size_t MultiDimensionalArraySpec::getIndex(const vector<size_t> &indices) const
{
  assert(m_dimensions.size() == indices.size());
  size_t index = indices[0];
  for (size_t i = 1; i < indices.size(); ++i) {
    index = (index * m_dimensions[i]) + indices[i];
  }
  return index;
}

vector<size_t> MultiDimensionalArraySpec::getIndices(size_t index) const
{
  vector<size_t> indices(m_dimensions.size());
  for (size_t i = m_dimensions.size() - 1; i > 0; --i) {
    indices[i] = index % m_dimensions[i];
    index /= m_dimensions[i];
  }
  indices[0] = index;
  return indices;
}


/******************************************************************/
/**************************[END OF MultiDimensionalArray.cpp]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
