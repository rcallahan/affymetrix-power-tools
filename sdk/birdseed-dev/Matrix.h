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
 * FILE Matrix.h
 * This is a very rough matrix library that has just enough functionality
 * to implement FitSNPGaussiansPriors3().
 */

#ifndef _BIRDSEED_MATRIX_H
#define _BIRDSEED_MATRIX_H

//
#include <cassert>
#include <sstream>
#include <vector>

// uncomment below to get the assert statements back.
// #define _MATRIX_DEBUG_ASSERT_ 1

namespace birdseed
{
namespace dev
{
    template<class T, int LEN>
    class FixedVector
    {
      public:
        typedef T value_type;
        enum {kLen = LEN};

        T m_data[LEN];

        FixedVector() {}
        FixedVector(size_t size)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(size == LEN);
#endif
        }
        
        inline T &operator[](size_t i) 
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < size());
#endif
            return m_data[i];
        }

        inline const T &operator[](size_t i) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < size());
#endif
            return m_data[i];
        }

        inline size_t size() const 
        {
            return LEN;
        }

        // This is not really right, because the dimension of the output vector
        // is not the length of the input vector, but is the # of columns in the input matrix.
        template<class MATRIX>
        FixedVector operator*(const MATRIX &m) const
        {
            FixedVector ret;
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == m.numRows());
            assert(ret.size() == m.numCols());
#endif
            for (size_t i = 0; i < m.numCols(); ++i) {
                typename MATRIX::value_type acc = 0;
                for (size_t j = 0; j < this->size(); ++j) {
                    acc += (*this)[j] * m[j][i];
                }
                ret[i] = acc;
            }
            return ret;
        }

        inline value_type operator*(const FixedVector &other) const
        {
            value_type ret = 0;
            for (size_t i = 0; i < this->size(); ++i) {
                ret += (*this)[i] * other[i];
            }
            return ret;
        }
        
        template<class VECTOR>
        FixedVector operator-(const VECTOR &other) const
        {
            FixedVector ret;
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] - other[i];
            }
            return ret;
        }

        inline FixedVector operator+(const FixedVector &other) const
        {
            FixedVector ret;
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] + other[i];
            }
            return ret;
        }

        inline FixedVector operator*(value_type scalar) const
        {
            FixedVector ret;
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] * scalar;
            }
            return ret;
        }

        inline FixedVector operator-() const
        {
            FixedVector ret;
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = -((*this)[i]);
            }
            return ret;
        }

        inline FixedVector operator/(value_type scalar) const
        {
            FixedVector ret;
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i]/scalar;
            }
            return ret;
        }

        inline void setAllElements(value_type scalar) 
        {
            for (size_t i = 0; i < this->size(); ++i) {
                (*this)[i] = scalar;
            }
        }
    };
    

    template<class T, int ROWS, int COLS>
    class FixedMatrix
    {
      public:
        typedef T value_type;
        typedef FixedVector<T, ROWS> ColumnType;
        typedef FixedVector<T, COLS> RowType;
        typedef RowType ROW;
        enum {kRows = ROWS, kCols = COLS};
        
        
        RowType m_data[ROWS];

        FixedMatrix() {}
        FixedMatrix(size_t len)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(len == ROWS);
#endif
        }

        template<class MATRIX>
        FixedMatrix(const MATRIX &other)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(other.numRows() == ROWS);
          assert(other.numCols() == COLS);
#endif
            for (size_t row = 0; row < other.numRows(); ++row) {
                for (size_t col = 0; col < other.numCols(); ++col) {
                    (*this)[row][col] = other[row][col];
                }
            }
        }
        
        inline RowType &operator[](size_t i)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < numRows());
#endif
            return m_data[i];
        }

        inline const RowType &operator[](size_t i) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < numRows());
#endif
            return m_data[i];
        }
        
        inline size_t numRows() const 
        {
            return ROWS;
        }

        inline size_t numCols() const 
        {
            return COLS;
        }

        template<class VECTOR>
        ColumnType operator*(const VECTOR &v) const
        {
            ColumnType ret;
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(v.size() == this->numCols());
            assert(ret.size() == this->numRows());
#endif
            for (size_t i = 0; i < this->numRows(); ++i) {
                value_type acc = 0;
                for (size_t j = 0; j < v.size(); ++j) {
                    acc += v[j] * (*this)[i][j];
                }
                ret[i] = acc;
            }
            return ret;
        }

        FixedMatrix operator*(const T mult) const
        {
            FixedMatrix ret;
            for (size_t row = 0; row < numRows(); ++row) {
                for (size_t col = 0; col < numCols(); ++col) {
                    ret[row][col] = (*this)[row][col] * mult;
                }
            }
            return ret;
        }

        template<class MATRIX>
        FixedMatrix operator-(const MATRIX &minus) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(numRows() == minus.numRows());
          assert(numCols() == minus.numCols());
#endif
            FixedMatrix ret;
            for (size_t row = 0; row < numRows(); ++row) {
                for (size_t col = 0; col < numCols(); ++col) {
                    ret[row][col] = (*this)[row][col] - minus[row][col];
                }
            }
            return ret;
        }
            
        template<class MATRIX>
        FixedMatrix operator+(const MATRIX &matrix2) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(numRows() == matrix2.numRows());
          assert(numCols() == matrix2.numCols());
#endif
            FixedMatrix ret;
            for (size_t row = 0; row < this->numRows(); ++row) {
                for (size_t col = 0; col < this->numCols(); ++col) {
                    ret[row][col] = (*this)[row][col] + matrix2[row][col];
                }
            }
            return ret;
        }

        
        ColumnType columnSlice(size_t col) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(col < numCols());
#endif
            ColumnType ret;
            for (size_t i = 0; i < numRows(); ++i) {
                ret[i] = (*this)[i][col];
            }
            return ret;
        }
    };

    template<class T>
    class VarVector: public std::vector<T>
    {
      public:
        typedef std::vector<T> base_type;
        typedef typename base_type::value_type value_type;
        VarVector(size_t len):
            base_type(len)
        {}
        
        VarVector()
        {}

        inline T &operator[](size_t i)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < this->size());
#endif
            return base_type::operator[](i);
        }
        
        inline const T &operator[](size_t i) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < this->size());
#endif
            return base_type::operator[](i);
        }
        
        inline VarVector operator*(value_type scalar) const
        {
            VarVector ret(this->size());
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] * scalar;
            }
            return ret;
        }

        inline VarVector operator+(value_type scalar) const
        {
            VarVector ret(this->size());
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] + scalar;
            }
            return ret;
        }

        template<class VECTOR>
        VarVector operator-(const VECTOR &other) const
        {
            VarVector ret(this->size());
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] - other[i];
            }
            return ret;
        }

        template<class VECTOR>
        VarVector operator+(const VECTOR &other) const
        {
            VarVector ret(this->size());
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] + other[i];
            }
            return ret;
        }

        template<class VECTOR>
        VarVector multiplyElementwise(const VECTOR &other) const
        {
            VarVector ret(this->size());
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] * other[i];
            }
            return ret;
        }
        template<class VECTOR>
        VarVector divideElementwise(const VECTOR &other) const
        {
            VarVector ret(this->size());
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] / other[i];
            }
            return ret;
        }

    };
    
    template<class T, int MAXLENGTH>
    class VarVectorWithReservedLength
    {
      public:
        enum{kMaxLength = MAXLENGTH};
        typedef T value_type;

        size_t m_len;
        T m_data[kMaxLength];
        
        VarVectorWithReservedLength(size_t len):
            m_len(len)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(len <= kMaxLength);
#endif
        }
        
        inline T &operator[](size_t i)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < this->size());
#endif
            return m_data[i];
        }
        
        inline const T &operator[](size_t i) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < this->size());
#endif
            return m_data[i];
        }
        
        inline VarVectorWithReservedLength operator*(value_type scalar) const
        {
            VarVectorWithReservedLength ret(this->size());
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] * scalar;
            }
            return ret;
        }

        inline VarVectorWithReservedLength operator+(value_type scalar) const
        {
            VarVectorWithReservedLength ret(this->size());
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] + scalar;
            }
            return ret;
        }

        template<class VECTOR>
        VarVectorWithReservedLength operator+(const VECTOR &other) const
        {
            VarVectorWithReservedLength ret(this->size());
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] + other[i];
            }
            return ret;
        }

        template<class VECTOR>
        VarVectorWithReservedLength operator-(const VECTOR &other) const
        {
            VarVectorWithReservedLength ret(this->size());
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] - other[i];
            }
            return ret;
        }

        template<class VECTOR>
        VarVectorWithReservedLength multiplyElementwise(const VECTOR &other) const
        {
            VarVectorWithReservedLength ret(this->size());
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] * other[i];
            }
            return ret;
        }
        template<class VECTOR>
        VarVectorWithReservedLength divideElementwise(const VECTOR &other) const
        {
            VarVectorWithReservedLength ret(this->size());
#ifdef _MATRIX_DEBUG_ASSERT_
            assert(this->size() == other.size());
#endif
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = (*this)[i] / other[i];
            }
            return ret;
        }

        // Return a vector in which each element is val/this[i]
        inline VarVectorWithReservedLength divisor(T val) const
        {
            VarVectorWithReservedLength ret(this->size());
            for (size_t i = 0; i < this->size(); ++i) {
                ret[i] = val/(*this)[i];
            }
            return ret;
        }
        
        inline void setAllElements(T val)
        {
            for (size_t i = 0; i < size(); ++i){
                m_data[i] = val;
            }
        }

        inline size_t size() const
        {
            return m_len;
        }

        inline void resize(size_t newlen)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(newlen <= kMaxLength);
#endif
            m_len = newlen;
        }
    };
    
    
    template<class T, int COLS>
    class VarMatrix
    {
      public:
        typedef T value_type;
        enum {kCols = COLS};

        typedef VarMatrix self;
        typedef FixedVector<T, COLS> ROW;

      private:
        std::vector<ROW> m_data;

      public:

        VarMatrix(size_t len):
            m_data(len)
        {
        }

        VarMatrix(size_t len, size_t reserve):
            m_data()
        {
            m_data.reserve(reserve);
            m_data.resize(len);
        }

        // To avoid picking up the template below if passing an int constant
        VarMatrix(int len):
            m_data(len)
        {
        }

        template<class MATRIX>
        VarMatrix(const MATRIX &other):
            m_data(other.numRows())
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(other.numCols() == COLS);
#endif
            for (size_t row = 0; row < other.numRows(); ++row) {
                m_data[row] = other[row];
            }
        }

        template<class MATRIX>
        VarMatrix &operator=(const MATRIX &other)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(other.numCols() == COLS);
#endif
            resize(other.numRows());
            for (size_t row = 0; row < other.numRows(); ++row) {
                m_data[row] = other[row];
            }
        }

        inline ROW &operator[](size_t i)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < numRows());
#endif
            return m_data[i];
        }
        
        inline const ROW &operator[](size_t i) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < numRows());
#endif
            return m_data[i];
        }

        inline size_t numRows() const 
        {
            return m_data.size();
        }

        inline size_t numCols() const 
        {
            return COLS;
        }

        inline VarMatrix operator*(const T mult) const
        {
            VarMatrix ret(numRows());
            for (size_t row = 0; row < numRows(); ++row) {
                for (size_t col = 0; col < numCols(); ++col) {
                    ret[row][col] = (*this)[row][col] * mult;
                }
            }
            return ret;
        }

        inline VarVector<T> columnSlice(size_t col) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(col < numCols());
#endif
            VarVector<T> ret(numRows());
            for (size_t i = 0; i < numRows(); ++i) {
                ret[i] = (*this)[i][col];
            }
            return ret;
        }

        inline void resize(size_t newNumRows) 
        {
            m_data.resize(newNumRows);
        }

        inline void push_back(const ROW &row)
        {
            m_data.push_back(row);
        }
    };

    template<class T, int MAXROWS, int COLS>
    class VarMatrixWithReservedRows
    {
      public:
        typedef T value_type;
        enum {kMaxRows = MAXROWS, kCols = COLS};

        typedef VarMatrixWithReservedRows  self;
        typedef FixedVector<T, COLS> ROW;

      private:
        size_t m_len;
        ROW m_data[MAXROWS];

      public:

        VarMatrixWithReservedRows(size_t len):
            m_len(len),
            m_data()
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(m_len <= kMaxRows);
#endif
        }

        // To avoid picking up the template below if passing an int constant
        VarMatrixWithReservedRows(int len):
            m_len(len),
            m_data()
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(m_len < kMaxRows);
#endif
        }

        template<class MATRIX>
        VarMatrixWithReservedRows &operator=(const MATRIX &other)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(other.numCols() == COLS);
          assert(other.numRows() <= kMaxRows);
#endif
            resize(other.numRows());
            for (size_t row = 0; row < other.numRows(); ++row) {
                m_data[row] = other[row];
            }
            return *this;
        }

        inline ROW &operator[](size_t i)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < numRows());
#endif
            return m_data[i];
        }
        
        inline const ROW &operator[](size_t i) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < numRows());
#endif
            return m_data[i];
        }

        inline size_t numRows() const 
        {
            return m_len;
        }

        inline size_t numCols() const 
        {
            return COLS;
        }

        inline VarMatrixWithReservedRows operator*(const T mult) const
        {
            VarMatrixWithReservedRows ret(numRows());
            for (size_t row = 0; row < numRows(); ++row) {
                for (size_t col = 0; col < numCols(); ++col) {
                    ret[row][col] = (*this)[row][col] * mult;
                }
            }
            return ret;
        }

        VarVectorWithReservedLength<T, MAXROWS> columnSlice(size_t col) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(col < numCols());
#endif
            VarVectorWithReservedLength<T, MAXROWS> ret(numRows());
            for (size_t i = 0; i < numRows(); ++i) {
                ret[i] = (*this)[i][col];
            }
            return ret;
        }

        inline void resize(size_t newNumRows) 
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(newNumRows <= kMaxRows);
#endif
            m_len = newNumRows;
        }

        template<class MATRIX>
        VarMatrixWithReservedRows operator-(const MATRIX &minus) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(numRows() == minus.numRows());
          assert(numCols() == minus.numCols());
#endif
            VarMatrixWithReservedRows ret(numRows());
            for (size_t row = 0; row < numRows(); ++row) {
                for (size_t col = 0; col < numCols(); ++col) {
                    ret[row][col] = (*this)[row][col] - minus[row][col];
                }
            }
            return ret;
        }
    };

    template<class T>
    class VarVarMatrix
    {
      public:
        typedef T value_type;
        typedef VarVector<T> ROW;

      private:
        std::vector<ROW> m_data;

      public:
        VarVarMatrix(size_t rows, size_t cols)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(rows > 0);
          assert(cols > 0);
#endif
            m_data.resize(rows);
            for (typename std::vector<ROW>::iterator it = m_data.begin();
                 it != m_data.end(); ++it) {
                (*it).resize(cols);
            }
        }
        
        VarVarMatrix(size_t rows, size_t cols, T initializer)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(rows > 0);
          assert(cols > 0);
#endif
            m_data.resize(rows);
            for (typename std::vector<ROW>::iterator it = m_data.begin();
                 it != m_data.end(); ++it) {
                (*it).resize(cols, initializer);
            }
        }

        inline ROW &operator[](size_t i)
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < numRows());
#endif
            return m_data[i];
        }
        
        inline const ROW &operator[](size_t i) const
        {
#ifdef _MATRIX_DEBUG_ASSERT_
          assert(i < numRows());
#endif
            return m_data[i];
        }

        inline size_t numRows() const 
        {
            return m_data.size();
        }

        inline size_t numCols() const 
        {
            return m_data[0].size();
        }

        inline void resizeRows(size_t newrows)
        {
            m_data.resize(newrows);
        }
    };
    
    template<class MATRIX>
    typename MATRIX::value_type sumColumn(const MATRIX &matrix, size_t column)
    {
        typename MATRIX::value_type total = 0.0;
        for (size_t i = 0; i < matrix.numRows(); ++i) {
            total += matrix[i][column];
        }
        return total;
    }

    template<class MATRIX>
    typename MATRIX::value_type mean(const MATRIX &matrix, size_t column)
    {
        return sumColumn(matrix, column)/matrix.numRows();
    }

    template<class MATRIX>
    typename MATRIX::value_type variance(const MATRIX &matrix, typename MATRIX::value_type theMean, size_t column)
    {
        if (matrix.numRows() < 2) {
            return 1.0;
        }
#ifdef _MATRIX_DEBUG_ASSERT_
        assert(matrix.numRows() > 1);
#endif
        typename MATRIX::value_type acc = 0.0;
        for (size_t i = 0; i < matrix.numRows(); ++i) {
            acc += (matrix[i][column] - theMean) * (matrix[i][column] - theMean);
        }
        if (acc == 0.0) {
            return 1.0;
        }
        return acc / (matrix.numRows() - 1);
    }

    template<class MATRIX>
    typename MATRIX::value_type minInColumn(const MATRIX &matrix, size_t column) 
    {
#ifdef _MATRIX_DEBUG_ASSERT_
      assert(matrix.numRows() > 0);
#endif
        typename MATRIX::value_type theMin = matrix[0][column];
        for (size_t i = 1; i < matrix.numRows(); ++i) {
            theMin = min(theMin, matrix[i][column]);
        }
        return theMin;
    }

    template<class MATRIX>
    typename MATRIX::value_type maxInColumn(const MATRIX &matrix, size_t column) 
    {
#ifdef _MATRIX_DEBUG_ASSERT_
      assert(matrix.numRows() > 0);
#endif
        typename MATRIX::value_type theMax = matrix[0][column];
        for (size_t i = 1; i < matrix.numRows(); ++i) {
            theMax = max(theMax, matrix[i][column]);
        }
        return theMax;
    }

    template<class VECTOR>
    typename VECTOR::value_type sumVector(const VECTOR &vec)
    {
        typename VECTOR::value_type acc = 0;
        for (size_t i = 0; i < vec.size(); ++i) {
            acc += vec[i];
        }
        return acc;
    }

    template<class MATRIX>
    typename MATRIX::value_type sumRow(const MATRIX &matrix, size_t row)
    {
        return sumVector(matrix[row]);
        /*
        typename MATRIX::value_type total = 0.0;
        for (size_t i = 0; i < matrix.numCols(); ++i) {
            total += matrix[row][i];
        }
        return total;
        */
    }
    
    
    template<class VECTOR>
    typename VECTOR::value_type meanVector(const VECTOR &vec)
    {
        return sumVector(vec)/vec.size();
    }

    template<class VECTOR>
    typename VECTOR::value_type maxInVector(const VECTOR &vec) 
    {
#ifdef _MATRIX_DEBUG_ASSERT_
      assert(vec.size() > 0);
#endif
        typename VECTOR::value_type theMax = vec[0];
        for (size_t i = 1; i < vec.size(); ++i) {
            theMax = max(theMax, vec[i]);
        }
        return theMax;
    }

    template<class VECTOR>
    typename VECTOR::value_type minInVector(const VECTOR &vec) 
    {
#ifdef _MATRIX_DEBUG_ASSERT_
      assert(vec.size() > 0);
#endif
        typename VECTOR::value_type theMin = vec[0];
        for (size_t i = 1; i < vec.size(); ++i) {
            theMin = min(theMin, vec[i]);
        }
        return theMin;
    }

    template<class MATRIX>
    MATRIX maxMatrixElementwise(const MATRIX &matrix, typename MATRIX::value_type val)
    {
        MATRIX ret(matrix.numRows());
        for (size_t row = 0; row < matrix.numRows(); ++row) {
            for (size_t col = 0; col < matrix.numCols(); ++col) {
                ret[row][col] = max(val, matrix[row][col]);
            }
        }
        return ret;
    }

    template<class VECTOR>
    std::string vectorToString(const VECTOR &vec)
    {
        std::stringstream strm;
        strm << "[";
        for (size_t i = 0; i < vec.size(); ++i) {
            if (i != 0) strm << ", ";
            strm << vec[i];
        }
        strm << "]";
        return strm.str();
    }
    
    
    template<class MATRIX>
    std::string matrixToString(const MATRIX &matrix)
    {
        std::stringstream strm;
        strm << "[";
        for (size_t row = 0; row < matrix.numRows(); ++row) {
            if (row > 0) {
                strm << ", ";
            }
            strm << "[";
            for (size_t col = 0; col < matrix.numCols(); ++col) {
                if (col > 0) {
                    strm << ", ";
                }
                strm << matrix[row][col];
            }
            strm << "]";
        }
        strm << "]";
        return strm.str();
    }
};
};



#endif /* _BIRDSEED_MATRIX_H */

/******************************************************************/
/**************************[END OF Matrix.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
