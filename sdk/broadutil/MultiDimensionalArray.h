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

#ifndef _MULTIDIMENSIONALARRAY_H
#define _MULTIDIMENSIONALARRAY_H

#include <cassert>
#include <vector>
//

using namespace std;

/*
 * Holds the dimensions of a dynamically-sized multi-dimensional array,
 * and converts between linear index and set of indices.
 * The convention is that the indices will be [Y][X].
 */
class MultiDimensionalArraySpec
{
private:
  vector<size_t> m_dimensions;
  size_t m_numElements;

public:
  MultiDimensionalArraySpec(const vector<size_t> &dimensions);

  size_t getNumElements() const {
    return m_numElements;
  }

  size_t getIndex(const vector<size_t> &indices) const;

  vector<size_t> getIndices(size_t index) const;

  const vector<size_t> &getDimensions() const {
    return m_dimensions;
  }
};

template <typename T>
class MultiDimensionalArray
{
private:
  MultiDimensionalArraySpec m_spec;
  vector<T> m_values;

public:
  MultiDimensionalArray(const vector<size_t> &dimensions):
      m_spec(dimensions) {
    m_values.resize(m_spec.getNumElements());
  }

  MultiDimensionalArray(const MultiDimensionalArraySpec &spec):
      m_spec(spec) {
    m_values.resize(m_spec.getNumElements());
  }

  T &operator[](const vector<size_t> &indices) {
    return m_values[m_spec.getIndex(indices)];
  }

  const T &operator[](const vector<size_t> &indices) const {
    return m_values[m_spec.getIndex(indices)];
  }


  T &operator[](size_t index) {
    return m_values[index];
  }

  const T &operator[](size_t index) const {
    return m_values[index];
  }

  // These are handy for arrays of bools, where you can't return a reference,
  // because vector is smart enough to represent it as a bit vector.
  void setVal(const vector<size_t> &indices, T val) {
    setVal(m_spec.getIndex(indices), val);
  }

  void setVal(size_t index, T val) {
    m_values[index] = val;
  }

  T getVal(const vector<size_t> &indices)  const {
    return getVal(m_spec.getIndex(indices));
  }

  T getVal(size_t index) const {
    return m_values[index];
  }

  const MultiDimensionalArraySpec &getSpec() const {
    return m_spec;
  }
};


#endif /* _MULTIDIMENSIONALARRAY_H */

/******************************************************************/
/**************************[END OF MultiDimensionalArray.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
