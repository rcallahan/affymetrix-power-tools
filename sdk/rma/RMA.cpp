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
 * @file   RMA.cpp
 * @author Chuck Sugnet
 * @date   Mon May 16 13:30:37 2005
 * 
 * @brief Functions to implement Robust Multiarray Average. 
 *
 * Implements the method described in "Exploration, normalization, and summaries
 * of high density oligonucleotide array probe level data", Irizarry RA, Hobbs
 * B, Collin F, Beazer-Barclay YD, Antonellis KJ, Scherf U, Speed TP.,
 * Biostatistics. 2003 Apr;4(2):249-64.
 * 
 */

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

//
#include "rma/RMA.h"
//
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#define M_1_SQRTPI      (0.5*M_2_SQRTPI)

using namespace std;

/** 
 * Subtract the median off of each row.
 * 
 * @param psMatrix - matrix.
 * @param numRow - number of rows.
 * @param numCol - number of columns.
 * @param rowEffect - values to subtract from each row.
 */
void RMA::subRowEffect(vector< vector<float> > &psMatrix, int numRow, int numCol, vector<float> &rowEffect) {
  int colIx = 0, rowIx = 0;
  /* some sanity checks. */
  assert(numRow > 0);
  assert(numCol > 0);
  
  for(rowIx = 0; rowIx < numRow; rowIx++) {
    for(colIx = 0; colIx < numCol; colIx++) {
      psMatrix[rowIx][colIx] = psMatrix[rowIx][colIx] - rowEffect[rowIx];
    }
  }
}

/** 
 * Find the median for the rows.
 * 
 * @param psMatrix - matrix.
 * @param numRow - number of rows in matrix.
 * @param numCol - number of rows in matrix.
 * @param colEffect - fill in row effect.
 */
void RMA::getRowEffect(vector< vector<float> > &psMatrix, int numRow, int numCol, vector<float> &rowEffect) {
  vector<float> row(numCol, 0.0);
  int colIx = 0, rowIx = 0;
  
  /* some sanity checks. */
  assert(numRow > 0);
  assert(numCol > 0);
  rowEffect.clear();

  /* loop through rows calculating effect for each one. */
  for(rowIx = 0; rowIx < numRow; rowIx++) {
    /* fill in row for effect calculation. */
    row.clear();
    for(colIx = 0; colIx < numCol; colIx++) {
      row.push_back(psMatrix[rowIx][colIx]);
    }
    /* median of row is effect. */
    rowEffect.push_back((float)median(row.begin(), row.end()));
  }
}

/** 
 * Subtract the median off of each row.
 * 
 * @param psMatrix - matrix.
 * @param numRow - number of rows.
 * @param numCol - number of columns.
 * @param rowEffect - values to subtract from each row.
 */
void RMA::subColEffect(vector< vector<float> > &psMatrix, int numRow, int numCol, vector<float> &colEffect) {
  int rowIx = 0, colIx = 0;

  /* some sanity checks. */
  assert(numCol > 0);
  assert(numRow > 0);
  
  for(colIx = 0; colIx < numCol; colIx++) {
    for(rowIx = 0; rowIx < numRow; rowIx++) {
      psMatrix[rowIx][colIx] = psMatrix[rowIx][colIx] - colEffect[colIx];
    }
  }
}

/** 
 * Find the median for the columns.
 * 
 * @param psMatrix - matrix.
 * @param numRow - number of rows in matrix.
 * @param numCol - number of columns in matrix.
 * @param colEffect - fill in column effect.
 */
void RMA::getColEffect(vector< vector<float> > &psMatrix, int numRow, int numCol,  vector<float> &colEffect) {
  vector<float> col(numRow, 0.0);
  int rowIx = 0, colIx = 0;
  
  /* some sanity checks. */
  assert(numCol > 0);
  assert(numRow > 0);
  colEffect.clear();

  /* loop through cols calculating effect for each one. */
  for(colIx = 0; colIx < numCol; colIx++) {
    /* fill in col for effect calculation. */
    col.clear();
    for(rowIx = 0; rowIx < numRow; rowIx++) {
      col.push_back(psMatrix[rowIx][colIx]);
    }
    /* median of col is effect. */
    colEffect.push_back((float)median(col.begin(), col.end()));
  }
}

/** 
 * Return the sum of a matrix.
 * @param psMatrix - matrix.
 * @param numRow - number of rows in matrix.
 * @param numCol - number of columns in matrix.
 * 
 * @return sum of matrix.
 */
double RMA::matrixSum(vector< vector<float> > &psMatrix, int numRow, int numCol) {
  double sum = 0.0;
  int rowIx = 0, colIx = 0;
  for(rowIx = 0; rowIx < numRow; rowIx++) {
    for(colIx = 0; colIx < numCol; colIx++) {
      sum += fabs(psMatrix[rowIx][colIx]);
    }
  }
  return sum;
}

/** 
 * Static version of medianPolishSet for when it isn't convenient to 
 * have an RMA instance.
 * 
 * @param PM - Matrix of perfect match probe intensities. Values
 *             will be modified in place, residuals will be filled in
 *             here.
 * @param numRow - Number of rows in the matrix.
 * @param numCol - Number of columns in the matrix.
 * @param indexes - Indexes of the probes for probe set of interest
 *        into PM matrix.
 * @param results - Vector in which probe sets summaries, one per GeneChip,
 *        will be returned.
 *
**/
void RMA::medianPolishPsetFromMatrix(vector< vector<float> > &pmMatrix, int numRow, int numCol, 
                                     vector<float> &colEstimates, vector<float> &rowEstimates, bool doLog) {
  vector<float> colEffect(numCol, 0.0);
  vector<float> rowEffect(numRow, 0.0);
  vector<float> workBuff(max(numCol, numRow), 0.0);
  float totalEffect = 0;
  int maxIter = 10; // How many iterations to do. 
  int iterIx = 0;   // Current iteration number.
  int i = 0, j = 0; // Indexes into arrays.
  double epsilon = 0.01;
  double oldSum = 0.0, newSum = 0.0;
  float delta = 0.0;
//  static bool first = true;
//   if(first) {
//     cerr.setf(ios::fixed, ios::floatfield);
//     cerr.precision(10);
//     for(i = 0; i < numRow; i++) {
//       for(j = 0; j < numCol; j++) {
//         cerr << pmMatrix[i][j] << "\t";
//       }
//       cerr << endl;
//     }
//  }
  /* sanity checks. */
  assert(numRow > 0 && numCol > 0);
  assert(numRow <= (int)pmMatrix.size());
  assert((int)pmMatrix.size() > 0 && numCol <= (int)pmMatrix[0].size());

  /* Clear out the estimates. */
  colEstimates.clear();
  rowEstimates.clear();
  colEstimates.insert(colEstimates.begin(), numCol, 0);
  rowEstimates.insert(rowEstimates.begin(), numRow, 0);

  if(doLog) {
    /* Log transform. */
    float log_2=logf(2.0);
    for(i = 0; i < numRow; i++) {
      for(j = 0; j < numCol; j++) {
          if(pmMatrix[i][j] <= 0.0){
            Err::errAbort("RMA::medianPolishPSetFromMatrix. Values must be strictly positive. "); 
          }
        pmMatrix[i][j] = logf(pmMatrix[i][j])/log_2;
      }
    }
  }

  // Here goes the polishing...
  for(iterIx = 0; iterIx < maxIter; iterIx++) {

  int k;
  if(maxIter == 9)
    k = 3;


    // Do the columns
    getColEffect(pmMatrix, numRow, numCol, workBuff);
    subColEffect(pmMatrix, numRow, numCol, workBuff);
    for(i = 0; i < numCol; i++) 
      colEffect[i] += workBuff[i];

    delta = (float)median(rowEffect.begin(), rowEffect.end());

    for(i = 0; i < numRow; i++) 
      rowEffect[i] -= delta;

    totalEffect += delta;

    // Do the rows    
    getRowEffect(pmMatrix, numRow, numCol, workBuff);
    subRowEffect(pmMatrix, numRow, numCol, workBuff);
    
    for(i = 0; i < numRow; i++) 
      rowEffect[i] += workBuff[i];

    delta = (float)median(colEffect.begin(), colEffect.end());

    for(i = 0; i < numCol; i++) 
      colEffect[i] -= delta;
    totalEffect += delta;

    // Check for convergance.
    newSum = matrixSum(pmMatrix, numRow, numCol);
    if(newSum == 0 || fabs(newSum - oldSum) < epsilon  * newSum) 
      break;
    oldSum = newSum;
  }

  /* data to return. */
  for(i = 0; i < numCol; i++) {
    colEstimates[i] = totalEffect + colEffect[i];
  }
  for(i = 0; i < numRow; i++) {
    rowEstimates[i] = totalEffect + rowEffect[i];
  }
//   if(first) {
//     cerr.setf(ios::fixed, ios::floatfield);
//     cerr.precision(10);
//     first = false;
//     for(i = 0; i < numRow; i++) {
//       cerr << rowEstimates[i] << "\t";
//     }
//     cerr << endl;
//   }
}


/** 
 * Static version of medianPolish.  In this version precomputed values are used as seed values to the usual 
 * median-polish algorithm.  It is allowed to run until the usual convergence criterion are met. 
 * 
 * 
 * @param PM - Matrix of perfect match probe intensities. Values
 *             will be modified in place, residuals will be filled in
 *             here.
 * @param numRow - Number of rows in the matrix.
 * @param numCol - Number of columns in the matrix.
 * @param indexes - Indexes of the probes for probe set of interest
 *        into PM matrix.
 * @param results - Vector in which probe sets summaries, one per GeneChip,
 *        will be returned.
 */

void RMA::medianPolishWithPrecomputedEffectsUsedAsSeedValues(vector< vector<float> > &pmMatrix, int numRow, int numCol, 
                                     vector<float> &colEstimates, vector<float> &rowEstimates, bool doLog) {
  vector<float> colEffect(numCol, 0.0);
  vector<float> rowEffect(numRow, 0.0);
  vector<float> workBuff(max(numCol, numRow), 0.0);
  float totalEffect = 0;
  int maxIter = 10; // How many iterations to do. 
  int iterIx = 0;   // Current iteration number.
  int i = 0, j = 0; // Indexes into arrays.
  double epsilon = 0.01;
  double oldSum = 0.0, newSum = 0.0;
  float delta = 0.0;
  /* sanity checks. */
  assert(numRow > 0 && numCol > 0);
  assert(numRow <= (int)pmMatrix.size());
  assert((int)pmMatrix.size() > 0 && numCol <= (int)pmMatrix[0].size());

  /* Clear out the estimates. */
  rowEstimates.clear();
  rowEstimates.insert(rowEstimates.begin(), numRow, 0);

  if(doLog) {
    /* Log transform. */
    float log_2=logf(2.0);
    for(i = 0; i < numRow; i++) {
      for(j = 0; j < numCol; j++) {
          if(pmMatrix[i][j] <= 0.0){
            Err::errAbort("RMA::medianPolishWithPrecomputedEffectsUsedAsSeedValues. Values must be strictly positive. "); 
          }
        pmMatrix[i][j] = logf(pmMatrix[i][j])/log_2;
      }
    }
  }
  bool firstTime = true;
  // Here goes the polishing...
  for(iterIx = 0; iterIx < maxIter; iterIx++) {

    // Do the columns
    if(firstTime == true){
      firstTime = false;
      for(i=0; i<numCol; i++)
        workBuff[i] = colEstimates[i];
    } else{ 
      getColEffect(pmMatrix, numRow, numCol, workBuff);
    }
    subColEffect(pmMatrix, numRow, numCol, workBuff);
    for(i = 0; i < numCol; i++) 
      colEffect[i] += workBuff[i];

    delta = (float)median(rowEffect.begin(), rowEffect.end());

    for(i = 0; i < numRow; i++) 
      rowEffect[i] -= delta;

    totalEffect += delta;

    // Do the rows    
    getRowEffect(pmMatrix, numRow, numCol, workBuff);
    subRowEffect(pmMatrix, numRow, numCol, workBuff);
    
    for(i = 0; i < numRow; i++) 
      rowEffect[i] += workBuff[i];

    delta = (float)median(colEffect.begin(), colEffect.end());

    for(i = 0; i < numCol; i++) 
      colEffect[i] -= delta;
    totalEffect += delta;

    // Check for convergance.
    newSum = matrixSum(pmMatrix, numRow, numCol);
    if(newSum == 0 || fabs(newSum - oldSum) < epsilon  * newSum) 
      break;
    oldSum = newSum;
  }

  /* data to return. */
  for(i = 0; i < numCol; i++) {
    colEstimates[i] = totalEffect + colEffect[i];
  }
  for(i = 0; i < numRow; i++) {
    rowEstimates[i] = totalEffect + rowEffect[i];
  }
}




/**
 * Static version of medianPolishWithPrecomputedEffects.  This version does the last steps of the usual
 * median-polish algorithm but does not iterate until convergence. If the feature effects that are used 
 * have been generated from the same data on which the analysis is done, the chip effect values will not
 * change.  
 *
 * @param pmMatrix -            Matrix of perfect match probe intensities. Values
 *                              will be modified in place, residuals will be filled in
 *                              here.
 * @param numRow -              Number of rows in the matrix.
 * @param numCol -              Number of columns in the matrix.
 * @param colEstimates -        This is the precomputed vector of probe (feature) effects.
 * @param rowEstimates -        This is the output vector of chip effects.
 * @param doLog -		This is a boolean switch to turn on or off the computation of logs of the input matrix.
 */

  void RMA::medianPolishWithPrecomputedEffectsOnePass(vector< vector<float> > &pmMatrix, int numChips, int numProbes,
                                vector<float> &colEstimates, vector<float> &rowEstimates, bool doLog){

    // In this function we subtract off the precalculated feature effects (colEstimates) from the pmMatrix,
    // calculate the medians of the remaining rows, and to these medians, add the median of the precalculated
    // feature effects.  These are the last steps in the original medianPolish algorithm.  

    if(doLog) {
      float log_2=logf(2.0);
      for(int i = 0; i < numChips; i++) {
        for(int j = 0; j < numProbes; j++) {
          if(pmMatrix[i][j] <= 0.0){
            Err::errAbort("RMA::medianPolishWithPrecomputedEffectsOnePass. Values must be strictly positive. "); 
          }
          pmMatrix[i][j] = logf(pmMatrix[i][j])/log_2;
        }
      }
    }

    rowEstimates.clear();

    subColEffect(pmMatrix, numChips, numProbes, colEstimates);
    getRowEffect(pmMatrix, numChips, numProbes, rowEstimates);
    float delta = (float)median(colEstimates.begin(), colEstimates.end());

    vector<float>::iterator begin = rowEstimates.begin();
    vector<float>::iterator end = rowEstimates.end();
    for(; begin!=end; begin++){
      *begin += delta;
    } 
  }


/** 
 * Compute the mean of a given vector.
 * sum(column)/n
 * @param dat - vector of doubles to take average of.
 * @return - mean value of vector
 */
double RMA::vectorMean(vector<float> &dat) {
  double sum = 0;
  unsigned int i = 0;
  for(i = 0; i < dat.size(); i++) {
    sum += dat[i];
  }
  return sum / dat.size();
}


/** 
 * Compute the variance of a given column. 
 * (sum(average - actual))^2 / n-1
 * @param dat - vector of doubles to take variance of.
 * @return - variance of vectorn
 */
double RMA::vectorVariance(vector<float> &dat) {
  double mean = 0;
  double sumOfSquares = 0;
  unsigned int i = 0;
  
  assert(dat.size() > 0);

  mean = vectorMean(dat);
  for(i = 0; i < dat.size(); i++) {
    sumOfSquares += (mean - dat[i]) * (mean - dat[i]);
  }
  return sumOfSquares / (dat.size() - 1);
}

/** 
 * Return the value at a particular quantile, Example quantile
 * would be .5 which is equivalent to the median.
 * 
 * @param dat - Data of interest.
 * @param quantile - Quantile of interest to report between [0,1];
 * @return value at quantile.
 */
double RMA::vectorQuantile(vector<float> &dat, double quantile) {
  vector<float> copy(dat);
  int i = 0;
  double r, f, result;

  /* Check to make sure there is something in the vector. */
  assert(dat.size() > 0);

  /* Check for size == 1 where can't do i+1. */
  if(dat.size() == 1)
    return dat[0];

  std::sort(copy.begin(), copy.end());

  /* Find the indexes of interest. */
  r = (dat.size() - 1) * quantile;
  i = (int)floor(r);
  f = r - i;

  result = ((1 - f) * copy[i]) + (f * copy[i+1]);
  return result;
}
  
/** 
 * Return the range of values from .75 quantile to .25 quantile.
 * 
 * @param dat - Data of interest.
 * @return - difference between .75 and .25 quantiles.
 */
double RMA::vectorIqr(vector<float> &dat) {
  double lowQuant = vectorQuantile(dat, .25);
  double highQuant = vectorQuantile(dat, .75);
  return highQuant - lowQuant;
}

/** 
 * Find the minimum value in a vector.
 * 
 * @param dat - vector of data.
 * 
 * @return minimum value.
 */
double RMA::vectorMin(vector <float> &dat) {
  unsigned int i = 0;
  float minimum = -1;
  assert(dat.size() > 0);
  minimum = dat[0];
  for(i = 0; i < dat.size(); i++) {
    minimum = min(minimum, dat[i]);
  }
  return minimum;
}

/** 
 * Find the maximum value in a vector.
 * 
 * @param dat - vector of data.
 * 
 * @return maximum value.
 */
double RMA::vectorMax(vector <float> &dat) {
  unsigned int i = 0;
  float maximum = -1;
  assert(dat.size() > 0);
  maximum = dat[0];
  for(i = 0; i < dat.size(); i++) {
    maximum = max(maximum, dat[i]);
  }
  return maximum;
}

/** 
 * Implementation of the bw.nrd0() function in R to decide kernel
 * bandwidth.
 * 
 * @param dat - vector of data of interest.
 * @return - Bandwidth of kernel to use.
 */
double RMA::findBandWidth(vector<float> &dat) {
  double bandWidth = 0;
  double high = 0, low = 0;
  assert(dat.size() > 1);
  high = sqrt(vectorVariance(dat));
  low = min(high, vectorIqr(dat)/1.34);
  if(!(fabs(low) > 0)) {
      Verbose::out(1, "In RMA::findBandWidth, boundary condition hit. Possible problem with input data.");
      low = high;
  }
  if(!(fabs(low) > 0)) {
      Verbose::out(1, "In RMA::findBandWidth, boundary condition hit. Possible problem with input data.");
      low = fabs(dat[0]);
  }
  if(!(fabs(low) > 0)) {
      Verbose::out(1, "In RMA::findBandWidth, boundary condition hit. Possible problem with input data.");
      low = 1.0;
  }
  assert(low);
  bandWidth = .9 * low * pow(1.0*dat.size(), -0.2);
  return bandWidth;
}

/** 
 * Spread the mass around our histogram bins.
 * 
 * @param dat - data of interest.
 * @param bandWidth - Bandwidth of mass distribution.
 * @param yVec - vector of bins to spread mass into. 
 * @param numBins - number of bins in yVec. 
 */
void RMA::colMassDist(vector<float> &dat, double bandWidth, double from, 
                      double to, vector<float> &yVec, int numBins) {
  double xDelta = 0;
  double ixMin = 0;
  double ixMax = numBins - 2;
  double ixMass = 0;
  double minimum, maximum;
  int i = 0;
  assert(dat.size() > 0);

  ixMass = 1.0 / dat.size();
  
  maximum = vectorMax(dat);
  minimum = vectorMin(dat);
  from = minimum - 7 * bandWidth;
  to = maximum + 7 * bandWidth;

  xDelta = (to - from) / (numBins -1);
  for(i = 0; i < 2*numBins; i++)
    yVec[i] = 0;

  for(i = 0; i < (int)dat.size(); i++) {
    double xpos = (dat[i] - from) / xDelta;
    int ix = (int)floor(xpos);
    double fx = xpos - ix;
    if(ixMin <= ix && ix <= ixMax) {
      yVec[ix] += (float)(1 - fx);
      yVec[(ix+1)] += (float)fx;
    }
    else if(ix == -1) {
      yVec[0] += (float)fx;
    }
    else if(ix == ixMax + 1) {
      yVec[ix] += (float)(1 - fx);
    }
  }

  for(i = 0; i < numBins; i++)
    yVec[i] *= (float)ixMass;
}

/** 
 * Convolves the real numbers in two vectors. Fill in kerelwords with
 * the convolved data. Note that yVec and kernel words must be a size
 * which is a power of two (i.e. 256, 512, 1024, etc. )
 * 
 * @param yVec - binned up data.
 * @param kernelWords - kernel function.
 * @param count - Number of entries in yVec and kernelWords. 
 */
void RMA::convolveDataKernel(vector<float> &yVec, vector<float> &kernelWords, int count) {
  vector<float> y(2*count,0.0);
  vector<float> k(2*count,0.0);
  vector<float> r(2*count,0.0);
  int cCount = 2 * count; /* Count for storing as complex type for fft. */
  int i = 0;

  /* Check to make sure we have a power of two. */
  assert(yVec.size() > 0);
  /// @todo This check looks ok, but on some hardware very small deviations make the floor and ceil fail
  ///  assert(floor(log((double)yVec.size()) / log(2.0)) == ceil(log((double)yVec.size()) / log(2.0)));

  /* Get data ready for fft call. */
  for(i = 0; i < count; i++) {
    y[i*2] = yVec[i];
    y[i*2+1] = 0;
    k[i*2] = kernelWords[i];
    k[i*2+1] = 0;
  }

  fft(y, false);
  fft(k, false);

  /* multiply results to convolve. */
  for(i = 0; i < cCount; i+=2) {
    r[i] = y[i] * k[i] + y[i+1] * k[i+1];
    r[i+1] = y[i+1]*k[i] - y[i]*k[i+1];
  }
  
  fft(r, true);

  for(i = 0; i < count; i++) {
    kernelWords[i] = r[2*i]; /* Pull out the real numbers as result. */
  }
}

/** 
 * Find the interpolation for value d between x and y. Assumes
 * that x and y are sorted.
 * 
 * @param d - value to interpolate.
 * @param x - vector 1
 * @param y - vector 2
 * 
 * @return interpolation of d from x into y. 
 */
double RMA::linearInterpolate(double d, vector<float> &x, vector<float> &y) {
  int start = 0, end = 0, mid = -1;
  double value = 0, yRange = 0, xRange = 0;
  
  assert(x.size() > 0 && y.size() > 0);
  end = int(x.size()) - 1;
  /* Check edge cases. */
  if(d < x[start]) 
    return y[start];
  if(d > x[end])
    return y[end];
  
  /* Do a binary search for values. */
  while( start + 1 < end ) {
    mid = (start + end) / 2;
    if(d > x[mid]) {  /* pick larger half of vector for further search */
      start = mid;
    }
    else { /* pick smaller half of vector for further search. */
      end = mid;
    }
  }
  assert(start+1 == end);

  /* Exact match? */
  if(d == x[start]) 
    return y[start];
  if(d == x[end])
    return y[end];
  
  yRange = y[end] - y[start];
  xRange = x[end] - x[start];
  value = y[start] + yRange * ((d - x[start])/xRange);
  return value;
}

/** 
 * Linearly interpolate x and y using values in x out. 
 * 
 * @param x - vector 1 for interpolation from.
 * @param y - vector 2 for interpolation to.
 * @param xOut - vector with data to interpolate.
 * @param yOut - output vector for interpolated data.
 */
void RMA::linearApprox(vector<float> &x, vector<float> &y, 
                       vector<float> &xOut, vector<float> &yOut) {
  const size_t size = x.size();
  size_t i = 0;

  assert(x.size() > 0 && y.size() > 0 && xOut.size() > 0);
  assert(x.size() == y.size());

  for(i = 0; i < size; i++) {
    yOut[i] = (float)linearInterpolate(xOut[i], x, y);
  }
}

void RMA::printVector(vector<float> &dat) {
  unsigned int i = 0;
  for(i = 0; i < dat.size(); i++) {
    cerr << dat[i] << ", ";
  }
  cerr << endl;
}

void RMA::printArray(float *dat, int size, float thresh){
  int i = 0;
  for(i = 0; i < size; i++) {
    if(fabs(dat[i]) >= thresh) 
      cerr << dat[i] << ",";
  }
  cerr << endl;
}

void RMA::printArray(double *dat, int size, double thresh){
  int i = 0;
  for(i = 0; i < size; i++) {
    if(fabs(dat[i]) >= thresh) 
      cerr << dat[i] << ",";
  }
  cerr << endl;
}


double RMA::arrayMax(double *d, int count) {
  double m = d[0];
  int i =0;
  for(i = 0; i < count; i++) {
    m = max(m, d[i]);
  }
  return  m;
}

/** 
 * Fill in a vector from a given column of data. 
 * 
 * @param PM - Matrix to fill in from.
 * @param numRows - Number of rows in matrix, will be vector length.
 * @param colIx - Column number to create vector from.
 * @param toFill - Vector to fill in with column data. 
 */
void RMA::columnToVector(double **PM, int numRows, int colIx, vector<float> &toFill) {
  int i = 0;
  assert(PM != NULL && numRows > 0);
  assert(toFill.size() == 0); // This vector should be empty when we get it. 
  toFill.reserve(numRows);    // avoid reallocation for efficiency.
  for(i = 0; i < numRows; i++) {
    toFill.push_back((float)PM[i][colIx]);
  }
}
    
/** 
 * Here we are recreating the density() function in R with an
 * Epanechnikov kernel. There are three
 * steps 1) Compute a Epanechikov kernel. 2) Convolve the mass
 * distribution with a fast fourier transform. 3) Linearly interpolate
 * the results.
 * @param dat - data of interest.
 * @param numBins - number of bins to spread data over (16384 for original
 *                  RMA method).
 * 
 * @return - Mode of the data.
 */
double RMA::modeOfData(vector<float> &dat, 
                       int numBins) {
  double bandWidth = 0;
  vector<float> yVec(2 * numBins, 0.0);
  vector<float> kernelWords(2 * numBins, 0.0);
  double from, to, maximum, minimum, upper, downer;
  vector<float> xords(numBins, 0.0), xOut(numBins, 0.0), 
    output(numBins, 0.0), kernelWordsVec(numBins);
  int i = 0;
  int maxIx = 0;
  double width = 0;
  double mode = 0;

  assert(dat.size() > 0);
  bandWidth = findBandWidth(dat);
  maximum = RMA::vectorMax(dat);
  minimum = RMA::vectorMin(dat);

  /* Set the limits. */
  from = minimum - 7 * bandWidth;
  to = maximum + 7 * bandWidth;
  downer = from + 4 * bandWidth;
  upper = to - 4 * bandWidth;

  /* Spread the mass of the distribution around. */
  colMassDist(dat, bandWidth, from, to, yVec, numBins);

  for(i = 0; i < numBins * 2; i++) {
    kernelWords[i] = float(i * (2*(to - from))/(numBins*2-1));
  }
  for(i = numBins+1; i < numBins * 2; i++) {
    kernelWords[i] = float(-1.0 * kernelWords[numBins - (i - numBins)]);
  }

  /* Apply the epanechnikov kernel. */
  width = bandWidth * sqrt(5.0);
  for(i = 0; i < numBins * 2; i++) {
    double abKw = fabs(kernelWords[i]);
    double abKwSq = (1.0 - (abKw/width * abKw/width));
    if(width > abKw) 
      kernelWords[i] = float(.75 * abKwSq / width);
    else
      kernelWords[i] = 0;
  }

  /* Use a fast fourier tranform to combine kernel and data. */
  convolveDataKernel(yVec, kernelWords, 2*numBins);

  /* Set the range for the interpolation. */
  for(i = 0; i < numBins; i++) {
    xords[i] =  float(from + i * ((to - from) / (numBins -1)));
    xOut[i] = float(downer + i * ((upper - downer) / (numBins-1)));
    kernelWordsVec[i] = kernelWords[i] / numBins;
  }

  /* interpolate. */
  linearApprox(xords, kernelWordsVec, xOut, output);

  /* Find bin with the most weight, i.e. the mode. */
  maximum = output[0];
  maxIx = 0;
  for(i = 0; i < numBins; i++) {
    if(maximum <= output[i]) {
      maximum = output[i];
      maxIx = i;
    }
  }
  
  mode = xOut[maxIx];
  return mode;
}

/** 
 * Calculate the parameters necesary to do the RMA background 
 * subtraction.
 * @param data - all of the data for a given cel file.
 * @param mu - pointer to parameter mu to be filled in by function.
 * @param sigma - pointer to parameter sigma to be filled in by function.
 * @param alpha - pointer to paramter alpha to be filled in.
 * @param numBins - number of bins to use in denstiy. RMA traditionally uses
 *                  16,384. It is essential that this number be a power of two.
 */
void RMA::estimateBgParam(vector<float> &data, double *mu, 
                          double *sigma, double *alpha, int numBins) {

  vector<float> upper, lower; // For upper and lower quantiles.
  double mode = 0, expMean = 0;
  double upperSum = 0, lowerSum = 0;
  unsigned int i =0, upperCount = 0, lowerCount = 0;
  assert(data.size() > 0 && mu != NULL && sigma != NULL && alpha != NULL);
  assert( (int)Util::round(pow(2.0, Util::round(log((double)numBins) / log(2.0)))) == numBins);
  
  /* Get the mode of all the data. */
  mode = modeOfData(data, numBins);

  /* Get the data in the lower .5 quantile. */
  for(i = 0; i < data.size(); i++) {
    if(data[i] < mode) {
      lower.push_back(data[i]);
    }
  }
  
  mode = modeOfData(lower, numBins);
  lower.clear();
  
  /* Divide out into above background and below background. */
  for(i = 0; i < data.size(); i++) {
    if(data[i] < mode) {
      double tmp = data[i] - mode;
      lower.push_back((float)tmp);
      lowerSum += (tmp * tmp);
      lowerCount++;
    }
    else {
      double tmp = data[i] - mode;
      upper.push_back((float)tmp);
      upperSum += (tmp * tmp);
      upperCount++;
    }
  }
  expMean = modeOfData(upper, numBins);
  (*alpha) = 1 / expMean;
  (*mu) = mode;
  (*sigma) = sqrt((double)lowerSum/(lowerCount - 1)) * sqrt(2.0);
}

/** 
 * The cdf function for standard normal distribution.
 * @param x - value.
 * @return - probability that standard normal value has
 * probability <= x.
 */
double RMA::PHI(double x) {
  return pnorm(x);
}

/** 
 * Standard normal probability density function.
 * @param x - value of interest.
 * @return - density at value supplied
 */
double RMA::phi(double x){
  double pi = 3.14159265358979323846; 
  return 1 / sqrt(2 * pi)* exp(-0.5 * x * x);
}

/** 
 * Cumulative normal density function.
 * @param u - value for which to calculate the distribution function.
 * @return - probability that value is < u.
 */
double RMA::pnorm(double u)
{
  const double a[5] = {
    1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001,
    1.887426188426510e+002,3.209377589138469e+003
  };
  const double b[5] = {
    1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002,
    1.813893686502485e+003,8.044716608901563e+003
  };
  const double c[9] = {
    2.15311535474403846e-8,5.64188496988670089e-1,8.88314979438837594e00,
    6.61191906371416295e01,2.98635138197400131e02,8.81952221241769090e02,
    1.71204761263407058e03,2.05107837782607147e03,1.23033935479799725E03
  };
  const double d[9] = {
    1.00000000000000000e00,1.57449261107098347e01,1.17693950891312499e02,
    5.37181101862009858e02,1.62138957456669019e03,3.29079923573345963e03,
    4.36261909014324716e03,3.43936767414372164e03,1.23033935480374942e03
  };
  const double p[6] = {
    1.63153871373020978e-2,3.05326634961232344e-1,3.60344899949804439e-1,
    1.25781726111229246e-1,1.60837851487422766e-2,6.58749161529837803e-4
  };
  const double q[6] = {
    1.00000000000000000e00,2.56852019228982242e00,1.87295284992346047e00,
    5.27905102951428412e-1,6.05183413124413191e-2,2.33520497626869185e-3
  };
  register double y, z;

  y = fabs(u);
  if (y <= 0.46875*M_SQRT2) {
    /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
    z = y*y;
    y = u*((((a[0]*z+a[1])*z+a[2])*z+a[3])*z+a[4])
      /((((b[0]*z+b[1])*z+b[2])*z+b[3])*z+b[4]);
    return 0.5+y;
  }
  z = exp(-y*y/2)/2;
  if (y <= 4.0) {
    /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
    y = y/M_SQRT2;
    y = ((((((((c[0]*y+c[1])*y+c[2])*y+c[3])*y+c[4])*y+c[5])*y+c[6])*y+c[7])*y+c[8])
        /((((((((d[0]*y+d[1])*y+d[2])*y+d[3])*y+d[4])*y+d[5])*y+d[6])*y+d[7])*y+d[8]);
    y = z*y;
  } else {
    /* evaluate erfc() for |u| > sqrt(2)*4.0 */
    z = z*M_SQRT2/y;
    y = 2/(y*y);
    y = y*(((((p[0]*y+p[1])*y+p[2])*y+p[3])*y+p[4])*y+p[5])
      /(((((q[0]*y+q[1])*y+q[2])*y+q[3])*y+q[4])*y+q[5]);
    y = z*(M_1_SQRTPI-y);
  }
  return (u < 0.0 ? y : 1-y);
};

/** 
 * RMA background correct subtraction. The RMA bacground subtraction
 * is a little complicated. The details are available in Ben Bolstad's
 * thesis dissertation:
 * http://stat-www.berkeley.edu/users/bolstad/Dissertation/Bolstad_2004_Dissertation.pdf
 * The basic idea is that the signal is exponential and the background
 * is additive gaussian. There are a lot of fun integrals to show that
 * the joint distribution can be found to estimate the true signal
 * given the observed signal on the chip and the parameters mu and
 * sigma for the background and alpha for the exponential. Estimating
 * mu, sigma and alpha explicitly is difficult though so RMA uses an
 * ad-hoc method where mu is estimated as the mode of the data, sigma
 * is estimated by the variation in the lower 50% quantile tail and
 * alpha is estimated as an exponential fit to the upper 50% quantile
 * tail. The "true" signal is then estimated as:
 * 
 * X = a + b(phi(a/b)/PHI(a/b))
 * where a = signal - mu - (sigma^2 * alpha)
 *       b = sigma
 
 *       phi and PHI are standard normal density function and
 *       distribution function respectively given as:
 *
 *       phi(z) = 1 / sqrt(2*pi) * exp(-.5*z^2)
 *       PHI(z) = \int_-inf^z 1 / sqrt(2*pi) exp(-.5w^2) dw
 * 
 * Note that this all happens on a per chip basis.
 * @param data - vector of data for a chip. 
 * @param numBins - number of bins to use for estimating mode
 *    of data for background parameters. Must be a power of 2.
 */
void RMA::backgroundCorrect(vector<float> &data, int numBins) {
  double mu = 0, sigma = 0, alpha = 0;
  unsigned int i;
  assert(data.size() > 0);
  estimateBgParam(data, &mu, &sigma, &alpha, numBins);
  bool error = false;
  for(i = 0; i < data.size(); i++) {
    data[i] = bgSubData(data[i], mu, sigma, alpha, error);
    if(error) {
       Err::errAbort("RMA::backgroundCorrect() - error in background calculation");
    }
  }
}

/** 
 * Do a fast fourier transform on data using the Danielson-Lanczos
 * algorithm. Values in data vector will be overwritten with the
 * results. The format of data must be that of real and imaginary
 * numbers interleaved together. So if the real portion of the number
 * is at index 0, the imaginary portion will be at index 1. The next
 * real number is at index 2 and imaginary portion is at index 3,
 * etc. The input data vector must be a power of two (i.e. 32, 64,
 * 128, etc) if your data is not a power of two in length you must pad
 * it with zeros. Be warned that this bit of code is not terribly user
 * friendly and can be tough to work with.
 * 
 * @param data - data of interest with interleaving real and imaginary
 * portions.
 * @param doInverse - should we perform the inverse fourier transform?
 */
void RMA::fft(vector<float> &data, bool doInverse)
{
  int length = 0, next = 0, step = 0;
  int i = 0, j = 0, k = 0;
  double sinHalfTemp = 0, sinTemp = 0, theta = 0;
  double wr = 0, wi = 0, wpr = 0;
  double realTmp = 0, imagTmp = 0;
  int sign = 1;

  /* Sanity checks */
  assert(data.size() > 0);
/// @todo This check looks ok, but on some hardware very small deviations make the floor and ceil fail
//   if(!(floor(log((double)data.size()) / log(2.0)) == ceil(log((double)data.size()) / log(2.0))))
//     Err::errAbort("RMA::fft() - data must be a power of two!");

  /* Are we calculating the inverse transform? */
  if(doInverse)
    sign = -1;
  
  /* Bit reversal. */
  j = 1;
  length = data.size();
  for (i = 1; i < length; i+= 2) {
    if (j > i) {
      realTmp = float(data[j-1]); data[j-1] = data[i-1]; data[i-1] = float(realTmp);
      realTmp = float(data[j]); data[j] = data[i]; data[i] = float(realTmp);
    }
    k = length / 2;
    while (k >= 2 && j > k) {
      j -= k;
      k /= 2;
    }
    j += k;
  }

  /* Danielson-Lanczos algorithm */
  next = 2;
  while (length > next) {
    step = next * 2;
    theta= sign * (6.28318530717959 / next); // Note that 6.28318530717959 is 2 * pi
    sinTemp = sin(theta);
    sinHalfTemp = sin(0.5 * theta);
    wpr = -2.0 * sinHalfTemp * sinHalfTemp;
    wr = 1.0;
    wi = 0.0;
    for (i = 1; i < next; i += 2) { 
      for (j = i;j <= length; j += step) {
        k = j + next;
        realTmp = wr * data[k-1] - wi * data[k];
        imagTmp= wr * data[k] + wi * data[k-1];
        data[k-1] = float(data[j-1] - realTmp);
        data[k] = float(data[j] - imagTmp);
        data[j-1] += (float)realTmp;
        data[j] += (float)imagTmp;
      }
      sinHalfTemp = wi;
      wi = wi * wpr + wr * sinTemp + wi;
      wr = wr * wpr - sinHalfTemp * sinTemp + wr;
    }
    next = step;
  }
}

struct _rmaPair 
/** Utility structure for sorting by data, but
    remembering rank. */
{
  double data; ///< intensity value.
  int rank;    ///< rank compared to other intensity values in chip.
};

/** 
 * Sort based on data in _rmaPair structure.
 * @param v1 pointer to _rmaPair structure.
 * @param v2 pointer to _rmaPair structure.
 * 
 * @return -1 if v1 < v2, 1 if v1 > v2, 0 if equal.
 */
static int pairCompare(const void *v1, const void *v2){
  struct _rmaPair *p1, *p2;
  p1 = (struct _rmaPair *)v1;
  p2 = (struct _rmaPair *)v2;

  if (p1->data < p2->data)
    return -1;
  else if (p1->data > p2->data)
    return 1;
  else
    return 0;
}


/** 
 * Get the ranks of _rmaPair structures using the same algorthim as R
 * to break ties.
 * @param pairs - array of _rmaPair strucutres.
 * @param count - size of array.
 * @param ranks - vector of ranks.
 */
void RMA::getRanks(struct _rmaPair *pairs, int count, vector<int> &ranks) {
  int current = 0, mark = 0;
  ranks.clear();
  while(current < count) {
    double total = 0;
    int i = 0;
    mark = current;
    while(mark < count -1 && pairs[mark].data == pairs[mark+1].data) {
      mark++;
    }
    total = floor((double)(mark + current + 2) / 2);
    for(i = 0; i < mark - current + 1; i++) {
      ranks.push_back((int)total);
    }
    current = mark + 1;
  }
}
        

/** 
 * Quantile Norm data using same algorithm to break ties as RMA
 * This function is deprecated as normalization sdk module now supports
 * bioconductor style quantile normalization.
 * @param data - matrix of data.
 */
void RMA::quantileNorm(vector< vector<float> > &data) {
  struct _rmaPair **rmaMatrix = NULL;
  vector<float> avgs;
  vector<int> ranks;
  int i = 0, rowIx = 0, colIx = 0;
  int numRow = (int)data.size();
  assert(numRow > 0);
  assert(data[0].size() > 0);

  /* There are issues with how affy sdk and bioconductor 
     break ties in the sorting. Here we are going to step into
     c-style qsort() and R version rank() to break ties to ensure
     compatability with bioconductor version of code. */

  /* Allocate memory. */
  rmaMatrix = new struct _rmaPair*[numRow];
  if(rmaMatrix == NULL) 
    Err::errAbort("RMA::quantileNorm() - Out of memory.");
  
  for(i = 0; i < numRow; i++) {
    rmaMatrix[i] = new struct _rmaPair[data[i].size()];
    if(rmaMatrix[i] == NULL) 
      Err::errAbort("RMA::quantileNorm() - Out of memory.");
  }

  /* Fill in data. */
  for(rowIx = 0; rowIx < numRow; rowIx++) {
    for(colIx = 0; colIx < (int)data[rowIx].size(); colIx++) {
      rmaMatrix[rowIx][colIx].data = data[rowIx][colIx];
      rmaMatrix[rowIx][colIx].rank = colIx;
    }
  }

  /* Sort everything. */
  for(rowIx = 0; rowIx < numRow; rowIx++) {
    qsort(rmaMatrix[rowIx], data[rowIx].size(), sizeof(struct _rmaPair), pairCompare);
  }

  /* Roll through and take the averages. */
  for(colIx = 0; colIx < (int)data[0].size(); colIx++) {
    double avg = 0.0;
    for(rowIx = 0; rowIx < numRow; rowIx++) {
      avg += rmaMatrix[rowIx][colIx].data / numRow;
    }
    avgs.push_back((float)avg);
  }
  for(rowIx = 0; rowIx < numRow; rowIx++) {
    getRanks(rmaMatrix[rowIx], data[rowIx].size(), ranks);
    for(colIx = 0; colIx < (int)data[0].size(); colIx++) {
      data[rowIx][rmaMatrix[rowIx][colIx].rank] = avgs[ranks[colIx]-1];
    }
  }

  for(i = 0; i < numRow; i++) {
    delete [] rmaMatrix[i];
  }
  delete [] rmaMatrix;
}
