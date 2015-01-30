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

/**
 * @file   SpectClust.h
 * @author Chuck Sugnet
 * @date   Fri Dec  8 13:46:28 2006
 * 
 * @brief Some utility functions for doing spectral clustering, very similar in
 * spirit to kernel PCA.
 * 
 */


#ifndef SPECTCLUST_H
#define SPECTCLUST_H

#include "portability/affy-base-types.h"
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Util.h"
//
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
//
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <vector>
//
typedef double Numeric;

/** Use dot product between two vectors as distance metric. */
class DistanceMetric {

public: 

  /** 
   * @brief Virtual destructor for a virtual class.
   */
  virtual ~DistanceMetric();

  /** 
   * Compute a distance metric between two columns of a
   * matrix. <b>Note that the indexes are *1* based (not 0) as that is
   * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
   * 
   * @param M - Matrix whose columns represent individual items to be clustered.
   * @param col1Ix - Column index to be compared (1 based).
   * @param col2Ix - Column index to be compared (1 based).
   * 
   * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
   */
  virtual double dist(const Matrix &M, int col1Ix, int col2Ix) const;

};

/** Use correlation between two vectors as distance metric. */
class CorrelationMetric : public DistanceMetric {

public: 
  
  /** 
   * Compute a distance metric between two columns of a
   * matrix. <b>Note that the indexes are *1* based (not 0) as that is
   * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
   * 
   * @param M - Matrix whose columns represent individual items to be clustered.
   * @param col1Ix - Column index to be compared (1 based).
   * @param col2Ix - Column index to be compared (1 based).
   * 
   * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
   */
  virtual double dist(const Matrix &M, int col1Ix, int col2Ix) const;
};

/** Use correlation between two vectors as distance metric. */
class GuassianRadial : public DistanceMetric {

public: 

  GuassianRadial(double sigma);
  
  /** 
   * Compute a distance metric between two columns of a
   * matrix. <b>Note that the indexes are *1* based (not 0) as that is
   * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
   * 
   * @param M - Matrix whose columns represent individual items to be clustered.
   * @param col1Ix - Column index to be compared (1 based).
   * @param col2Ix - Column index to be compared (1 based).
   * 
   * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
   */
  virtual double dist(const Matrix &M, int col1Ix, int col2Ix) const;

protected:
  /// Scaling factor.
  double m_Sigma;
};


/** Use angle between two vectors as a distance metric. */
class AngleMetric : public DistanceMetric {

public: 

  AngleMetric(const Matrix &M);

  /** 
   * Set the matrix. This precalculates the norms for each column so
   * they are only calculated once.
   * @param M - Matrix to be used.
   */
  void setMatrix(const Matrix &M);
  
  /** 
   * Compute a distance metric between two columns of a
   * matrix. <b>Note that the indexes are *1* based (not 0) as that is
   * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
   * 
   * In this case the distnance metric is the cosine of the angle between the
   * two vectors, this is also sometimes known as the uncentered correlation coefficient
   * and is similar to correlation except that it requires the length of the vectors to 
   * be the same for cos(x,y) to be 1.
   * 
   * If x and y are vectors then:
   * \f$ angleDist(\vec{x},\vec{y}) = cos(\theta) = \frac{ \vec{x}^t \cdot \vec{y}}{\|\vec{x}\|\|\vec{y}\|} \f$
   * Where the numerator is the dot product of <b>x</b> and <b>y</b>.
   * @param M - Matrix whose columns represent individual items to be clustered.
   * @param col1Ix - Column index to be compared (1 based).
   * @param col2Ix - Column index to be compared (1 based).
   * 
   * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
   */
  virtual double dist(const Matrix &M, int col1Ix, int col2Ix) const;
protected:
  /// Precalculated norms for the matrix we're working with.
  std::vector<double> m_Norms;
};

/** Use angle between two vectors as a distance metric. */
class ExpAngleMetric : public AngleMetric {

public: 

  ExpAngleMetric(const Matrix &M, double sigma);
  
  /** 
   * Compute a distance metric between two columns of a
   * matrix. <b>Note that the indexes are *1* based (not 0) as that is
   * Newmat's convention</b>. Note that dist(M,i,j) must equal dist(M,j,i);
   * 
   * In this case the distnance metric is the exponated 1 - cosine of the angle between the
   * two vectors, this is also sometimes known as the uncentered correlation coefficient
   * and is similar to correlation except that it requires the length of the vectors to 
   * be the same for cos(x,y) to be 1.
   * 
   * If x and y are vectors then:
   * \f$ ExpAngleDist(\vec{x},\vec{y}) = e^{-1 * (1 - cos(\theta))/2\sigma^2} \f$
   * Where:
   * \f$ cos(\theta) = \frac{ \vec{x}^t \cdot \vec{y}}{\|\vec{x}\|\|\vec{y}\|} \f$
   * Where the numerator is the dot product of <b>x</b> and <b>y</b>.
   * @param M - Matrix whose columns represent individual items to be clustered.
   * @param col1Ix - Column index to be compared (1 based).
   * @param col2Ix - Column index to be compared (1 based).
   * 
   * @return - "Distance" or "dissimilarity" metric between two columns of matrix.
   */
  virtual double dist(const Matrix &M, int col1Ix, int col2Ix) const;

  double m_Sigma;
};


class SpectClust {

public:

  static void multByMatrix(ColumnVector &N, ColumnVector &O, const Matrix &M);

  static double fast_corr(std::vector<double> &x, std::vector<double> &y);      

  static void fillInDistance(SymmetricMatrix &A, const Matrix &M, const DistanceMetric &dMetric, bool expon=true);

  static double rowMedian(const Matrix &M, int rowIx);

  static void subColAvg(Matrix &M);
  
  static void MatrixScatter(Matrix &M, Matrix &C);

  static void MatrixCov(Matrix &M, Matrix &C);

  static void MatrixCor(Matrix &M, Matrix &C);

  static void rowMedianDivide(Matrix &M);

  static void normalizeSum(SymmetricMatrix &A);

  static bool findNLargestSymEvals(const SymmetricMatrix &W, int numLamda, std::vector<Numeric> &eVals, Matrix &EVec);

  static bool MaxEigen(const Matrix &M, double &maxValue, ColumnVector &MaxVec, int maxIterations = 75);

  static bool findNLargestEvals(const Matrix &M, int numLamda, std::vector<Numeric> &eVals, Matrix &EVec, int maxIterations=75);

  static double calcNormalizedCut(const Matrix &D, 
                                  const std::vector<std::pair<double,int> > &indices, 
                                  int cut);


  static void calcNormalizedCutValues(const Matrix &D, const std::vector<std::pair<double,int> > &indices, 
                                      std::vector<double> &cutVals);

  static void colToVector(const Matrix &M, int colIx , std::vector<double> &v);

  static void partitionClusters(const Matrix &D, const Matrix &E, std::vector<Numeric> &eVals, int numClusters,
                                std::vector<int> &clusters, int hardMin, double cutVal = 0, double margin=1);

  static void orderIntraClusterCorr(const Matrix &M, vector<int> &clusters, 
                                    int numClusters, std::vector<double> &clusterCorr);

  static void AICcRmaMedianCluster(const Matrix &PM, int numClusters, const std::vector<int> &clusters, 
                                   std::vector<double> &AICc, std::vector<double> &BIC, bool doLog);

  struct PairLess : public std::binary_function<std::pair<double,int>, std::pair<double,int>, bool> {
    bool operator () (const std::pair<double,int>& x, const std::pair<double,int>& y) const
    {
      return (x.first < y.first);
    }
  };

};

#endif /* SPECTCLUST_H */
