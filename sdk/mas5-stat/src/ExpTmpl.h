////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#ifndef _EXPTMPL_H_
#define _EXPTMPL_H_

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>
//

using namespace std;

template <class T> class matrix {
private: vector<vector<T> > array;
public:
	matrix(): array(0) {}
	matrix(int rows, int cols) : array(rows) {
		for (int i = 0; i < rows; ++i)
			array[i].resize(cols);
	}
	void resize(int rows, int cols) {
		array.resize(rows);
		for (int i=0; i < rows; ++i)
			array[i].resize(cols);
	}
	matrix (int rows, int cols, T value) : array(rows) {
		for (int i = 0; i < rows; ++i) {
			array[i].resize(cols);
			for (int j = 0; j < cols; ++j)
				array[i][j] = value;
		}
	}
	const vector<T> & operator[](int row) const {
		return array[row];}
	vector<T> & operator[](int row) {
		return array[row];}
	int nRows() const {return array.size();}
	int nCols() const {return nRows() > 0 ? array[0].size() : 0;}
	double average() const {double sum = 0;
		for (int i=0; i < nRows(); ++i) 
			sum += mean(array[i]);
		return sum / (double) nRows();
	}
	double average(const int minr, const int maxr,
		const int minc, const int maxc) const {
		if (minr <0 || minr > maxr || maxr >= nRows()) {
			exit(1);
		}
		if (minc <0 || minc > maxc || maxc >= nCols()) {
			exit(1);
		}
		double sum=0;
		for (int i = minr; i <= maxr; ++i) {
			double rowSum = 0;
			for (int j = minc; j <= maxc; ++j) {
				rowSum += array[i][j];
			}
			sum += rowSum / (double) (maxc - minc + 1);
		}
		return sum / (double)(maxr - minr + 1);
	}
	const matrix<T> & operator=(const matrix<T> & rhs) {
		if (this != &rhs) {
			this->resize(rhs.nRows(), rhs.nCols());
			for (int i = 0; i < rhs.nRows(); ++i) {
				for (int j = 0; j < rhs.nCols(); ++j) {
					(*this)[i][j] = rhs[i][j];
				}
			}
		}
		return *this;
	}
};

struct ExpResults{
	double p_value;
	int call;
	ExpResults(): p_value(0), call(0) {}
	ExpResults(double p, int c): p_value(p), call(c) {}
	bool operator<(const struct ExpResults & rhs) const
	{return (p_value < rhs.p_value);}
	bool operator>(const struct ExpResults & rhs) const
	{return (p_value > rhs.p_value);}
};

struct FloatPair{
	float value1;
	float value2;
	FloatPair(): value1(0.0f), value2(0.0f) {}
	FloatPair(float v1, float v2): value1(v1), value2(v2) {}
	bool operator<(const struct FloatPair & rhs) const
	{return (value1 < rhs.value1);}
};

template<class T> ExpResults oneSidedSignRank2(const vector<T> & x, const double alpha);
template<class T> float mean(const vector<T> & v);
template<class T> float median(const vector<T> & v);
template<class T> float stddev(const vector<T> & v);
template<class T> ExpResults newSignRank(const vector<T> & dif, 
	const double alpha1, const double alpha2);
template<class T> float medianAbsoluteDeviation(const vector<T> & x);
template <class T> FloatPair trimMeanAndStd(vector<T> & v, const double p1, const double p2);
template <class T> double trimMean(const vector<T> & vec, const double p1, const double p2);

#endif /*_EXPTMPL_H_*/
