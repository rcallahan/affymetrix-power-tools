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

#ifndef  _MATHLIB_H
#define  _MATHLIB_H


#define REAL double
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
//

using namespace std;




class DI {
public:
	double d;
	int i;
	DI(): d(0), i(0) {}
	DI(double dv, int iv): d(dv), i(iv) {}
	bool operator<(const DI & rhs) const
		{return (d < rhs.d);}
};

class FMATRIX {
private: std::vector<std::vector<float> > array;
public:
	FMATRIX(): array() {}
	FMATRIX(int rows, int cols) : array(rows) {
		for (int i = 0; i < rows; ++i)
			array[i].resize(cols);
	}
	void resize(int rows) {
		array.resize(rows);
	}
	void resize(int rows, int cols) {
		array.resize(rows);
		for (int i=0; i < rows; ++i)
			array[i].resize(cols);
	}
	FMATRIX (int rows, int cols, float value) : array(rows) {
		for (int i = 0; i < rows; ++i) {
			array[i].resize(cols);
			for (int j = 0; j < cols; ++j)
				array[i][j] = value;
		}
	}
	FMATRIX (int rows, int cols, float value, float value2) : array(rows) {
		for (int i = 0; i < rows; ++i) {
			array[i].resize(cols);
			for (int j = 0; j < cols; ++j)
				array[i][j] = (i == j) ? value : value2;
		}
	}
	FMATRIX(const FMATRIX & mat) : array(mat.nRows()) {
		int n = mat.nCols();
		for (int i = 0; i < mat.nRows(); ++i) {
			array[i].resize(n);
			for (int j = 0; j < n; ++j) {
				array[i][j] = mat[i][j];
			}
		}
	}
	~FMATRIX() {
		array.clear();
	}

	const std::vector<float> & operator[](int row) const {
		return array[row];}
	std::vector<float> & operator[](int row) {
		return array[row];}
	int nRows() const {return (int) array.size();}
	int nCols() const {return nRows() > 0 ? (int) array[0].size() : 0;}
	const FMATRIX & operator=(const FMATRIX & rhs) {
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
	std::vector<float> getCol(const int jj) const {
		std::vector<float> col(nRows());
		for (int i = 0; i < nRows(); ++i) {
			col[i] = array[i][jj];
		}
		return col;
	}
};

class RMATRIX {
private: std::vector<std::vector<REAL> > array;
public:
	RMATRIX(): array() {}
	RMATRIX(int rows, int cols) : array(rows) {
		for (int i = 0; i < rows; ++i)
			array[i].resize(cols);
	}
	void resize(int rows) {
		array.resize(rows);
	}
	void resize(int rows, int cols) {
		array.resize(rows);
		for (int i=0; i < rows; ++i)
			array[i].resize(cols);
	}
	RMATRIX (int rows, int cols, REAL value) : array(rows) {
		for (int i = 0; i < rows; ++i) {
			array[i].resize(cols);
			for (int j = 0; j < cols; ++j)
				array[i][j] = value;
		}
	}
	RMATRIX (int rows, int cols, REAL value, REAL value2) : array(rows) {
		for (int i = 0; i < rows; ++i) {
			array[i].resize(cols);
			for (int j = 0; j < cols; ++j)
				array[i][j] = (i == j) ? value : value2;
		}
	}
	RMATRIX(const RMATRIX & mat) : array(mat.nRows()) {
		int n = mat.nCols();
		for (int i = 0; i < mat.nRows(); ++i) {
			array[i].resize(n);
			for (int j = 0; j < n; ++j) {
				array[i][j] = mat[i][j];
			}
		}
	}
	~RMATRIX() {
		array.clear();
	}

	const std::vector<REAL> & operator[](int row) const {
		return array[row];}
	std::vector<REAL> & operator[](int row) {
		return array[row];}
	int nRows() const {return (int) array.size();}
	int nCols() const {return nRows() > 0 ? (int) array[0].size() : 0;}
	const RMATRIX & operator=(const RMATRIX & rhs) {
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
	std::vector<REAL> getCol(const int jj) const {
		std::vector<REAL> col(nRows());
		for (int i = 0; i < nRows(); ++i) {
			col[i] = array[i][jj];
		}
		return col;
	}
};

/// @todo These should be in a cpp file.
const double PI = 3.1415926535897932384626433832795;
const double TINY = 1e-20;
const double Eps = numeric_limits<double>::epsilon();
const double RelativeBound = 5 * Eps;
const int MaxIterations = 5000;
const int MeaningLess = -9999;

float median(const std::vector<float> & v);
DI oneSidedSignRank3(const std::vector<float> & x, const double alpha);
double cFraction(double x, double a, double b);
double incompleteBeta(const double x, const double a, const double b);
double tCDF(double t, double df);
double tCDFinversed(double p, double df);
double logGamma(const double x);
double erf(const double x);
double normalCDF(const double x);
double t_distribution_lookup(const int df);

#endif
