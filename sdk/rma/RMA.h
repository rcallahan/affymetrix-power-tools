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
 * @file   RMA.h
 * @author Chuck Sugnet
 * @date   Mon May 16 13:34:44 2005
 * 
 * @brief  Functions to implement Robust Multiarray Average. 
 */

#ifndef RMA_H
#define RMA_H

//
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <vector>

/**
 * @class RMA 
 * @brief Holds methods for computing the popular Robust Muliarray Average of
 * Microarrays (RMA). 
 *
 * Implements the method described in "Exploration,
 * normalization, and summaries of high density oligonucleotide array probe
 * level data", Irizarry RA, Hobbs B, Collin F, Beazer-Barclay YD, Antonellis
 * KJ, Scherf U, Speed TP., Biostatistics. 2003 Apr;4(2):249-64.  Core routines
 * are two static functions for manipulating matrices. This class is mainly a
 * collection of static methods to manipulate vectors and matrices to perform an
 * RMA analysis. RMA really consists of three separate steps:
 *
 * 1) Background estimation which occurs at the chip level as implemented by
 *    backgroundCorrect(). The background subtraction is more complicated than
 *    a simple subtraction and accounts for a substantial amount of the code.

 * 2) Quantile normalization. As the original normalization module of the affy
 *    sdk used a different manner of resolving ties than the bioconductor version
 *    there are methods to normalize a matrix of chip data implemented in
 *    quantileNorm(). The normalization module has since been extended to reproduce
 *    the bioconductor results and the code here will be deprecated soon.
 *
 * 3) Probeset summary. Use the median polish algorithm to estimate the target
 *    response and feature response using medianPolishPsetFromMatrix(). Target
 *    response is used as the signal level estimates. While this is where the real
 *    "work" of doing probeset level estimates occurs, it takes up very little code
 *    as the median polish algorithm is simple to implement.
 */
class RMA {
  
public: 
  
/** 
 * Static version of medianPolishSet for when it isn't convenient to have an RMA
 * instance. Note that both the the colEstimates and rowEstimates by convention
 * have the global median added back to them. The identifiability function used
 * by RMA is median(log(f)) = 0 for both feature and target effects. By
 * convention the global median is added back onto target estimates to avoid
 * negative expression summaries.  The same is done for feature effects although
 * subtracting off the median for a probeset will give the median(log(f)) = 0
 * adjusted feature effects.
 * 
 * @param pmMatrix - Matrix of perfect match probe intensities. Values
 *             will be modified in place, residuals will be filled in
 *             here.
 * @param numRow - Number of rows in the matrix.
 * @param numCol - Number of columns in the matrix.
 * @param colEstimates - Feature effects.
 * @param rowEstimates - Target effects.
 */
static void medianPolishPsetFromMatrix(std::vector< std::vector<float> > &pmMatrix, 
                                       int numRow, int numCol, 
                                       std::vector<float> &colEstimates,
                                       std::vector<float> &rowEstimates,
                                       bool doLog = true);

/** 
 * Static version of medianPolishSet for when it isn't convenient to 
 * have an RMA instance.
 * 
 * @param pmMatrix - Matrix of perfect match probe intensities. Values
 *             will be modified in place, residuals will be filled in
 *             here.
 * @param numRow - Number of rows in the matrix.
 * @param numCol - Number of columns in the matrix.
 * @param indexes - Indexes of the probes for probe set of interest
 *        into PM matrix.
 * @param results - Vector in which probe sets summaries, one per GeneChip,
 *        will be returned as log2 values.
 */
static inline void medianPolish(std::vector< std::vector<float> > &pmMatrix,
                                int numRow, int numCol, 
                                std::vector<float> &colEstimates,
                                std::vector<float> &rowEstimates) {
  medianPolishPsetFromMatrix(pmMatrix, numRow, numCol, 
                             colEstimates, rowEstimates);
}


/**
 * Static version of medianPolishWithPrecomputedEffectsOnePass
 * This version will use precomputed feature effects as input values and perform the usual median polish 
 * algorithm until convergence 
 *
 * @param pmMatrix - 		Matrix of perfect match probe intensities. Values
 *             			will be modified in place, residuals will be filled in
 *             			here.
 * @param numRow - 		Number of rows in the matrix.
 * @param numCol - 		Number of columns in the matrix.
 * @param colEstimates - 	This is the precomputed vector of probe (feature) effects.
 * @param rowEstimates -	This is the output vector of chip effects. 
 * @param doLog -               This is a boolean switch to turn on or off the computation of logs of the input matrix.

 */
  static void medianPolishWithPrecomputedEffectsOnePass(std::vector< std::vector<float> > &pmMatrix,
                                                        int numChips, int numProbes, 
                                                        std::vector<float> &colEstimates,
                                                        std::vector<float> &rowEstimates, 
                                                        bool doLog = true);

/**
 * Static version of medianPolishWithPrecomputedEffectsUsedAsSeedValues
 * This version will use precomputed feature effects as input values and perform the usual median polish 
 * algorithm until convergence 
 *
 * @param pmMatrix - 		Matrix of perfect match probe intensities. Values
 *             			will be modified in place, residuals will be filled in
 *             			here.
 * @param numRow - 		Number of rows in the matrix.
 * @param numCol - 		Number of columns in the matrix.
 * @param colEstimates - 	This is the precomputed vector of probe (feature) effects.
 * @param rowEstimates -	This is the output vector of chip effects. 
 * @param doLog -               This is a boolean switch to turn on or off the computation of logs of the input matrix.

 */
  static void medianPolishWithPrecomputedEffectsUsedAsSeedValues(std::vector< std::vector<float> > &pmMatrix,
                                                                 int numChips, int numProbes, 
                                                                 std::vector<float> &colEstimates,
                                                                 std::vector<float> &rowEstimates, 
                                                                 bool doLog = true);

  /** 
   * Find the minimum value in a vector.
   * 
   * @param dat - vector of data.
   * 
   * @return minimum value.
   */
  static double vectorMin(std::vector <float> &dat);
  
  /** 
   * Find the maximum value in a vector.
   * 
   * @param dat - vector of data.
   * 
   * @return maximum value.
   */
  static double vectorMax(std::vector <float> &dat);

   /** 
    * Compute the mean of a given vector.
    * sum(column)/n
    * @param dat - vector of doubles to take average of.
    * @return - mean value of vector
    */
  static double vectorMean(std::vector<float> &dat);

   /** 
    * Compute the variance of a given column. 
    * (sum(average - actual))^2 / n-1
    * @param dat - vector of doubles to take variance of.
    * @return - variance of vectorn
    */
  static double vectorVariance(std::vector<float> &dat);

  /** 
   * Return the value at a particular quantile, Example quantile
   * would be .5 which is equivalent to the median.
   * 
   * @param dat - Data of interest.
   * @param quantile - Quantile of interest to report between [0,1];
   * @return value at quantile.
   */
  static double vectorQuantile(std::vector<float> &dat, double quantile);

  /** 
   * Return the range of values from .75 quantile to .25 quantile.
   * 
   * @param dat - Data of interest.
   * @return - difference between .75 and .25 quantiles.
   */
  static double vectorIqr(std::vector<float> &dat);

  /** 
   * Implementation of the bw.nrd0() function in R to decide kernel
   * bandwidth.
   * 
   * @param dat - vector of data of interest.
   * @return - Bandwidth of kernel to use.
   */
  static double findBandWidth(std::vector<float> &dat);

  /** 
   * Spread the mass around our histogram bins.
   * 
   * @param dat - data of interest.
   * @param bandWidth - Bandwidth of mass distribution.
   * @param from - begininning of dist.
   * @param to - end of values.
   * @param yVec - vector of bins to spread mass into. 
   * @param numBins - number of bins in yVec. 
   */
  static void colMassDist(std::vector<float> &dat, double bandWidth, double from, 
                          double to, std::vector<float> &yVec, int numBins);
  /** 
   * Fill in a vector from a given column of data. 
   * 
   * @param PM - Matrix to fill in from.
   * @param numRows - Number of rows in matrix, will be vector length.
   * @param colIx - Column number to create vector from.
   * @param toFill - Vector to fill in with column data. 
   */
  static void columnToVector(double **PM, int numRows, int colIx, std::vector<float> &toFill);

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
  static double linearInterpolate(double d, std::vector<float> &x, std::vector<float> &y);

  /** 
   * Linearly interpolate x and y using values in x out. 
   * 
   * @param x - vector 1 for interpolation from.
   * @param y - vector 2 for interpolation to.
   * @param xOut - vector with data to interpolate.
   * @param yOut - output vector for interpolated data.
   */
  static void linearApprox(std::vector<float> &x, std::vector<float> &y, 
                           std::vector<float> &xOut, std::vector<float> &yOut);
   /** 
    * Here we are recreating the density() function in R with an Epanechnikov
    * kernel. There are three steps 1) Compute a Epanechikov kernel. 2) Convolve
    * the mass distribution with a fast fourier transform. 3) Linearly
    * interpolate the results.
    * @param dat - data of interest.
    * @param numBins - number of bins to spread data over (16384 for original
    *                  RMA method).
    * 
    * @return - Mode of the data.
    */
  static double modeOfData(std::vector<float> &dat, int numBins);
  
  static void printArray(float *dat, int size, float thresh=0);
  static void printArray(double *dat, int size, double thresh=0);
  static double arrayMax(double *d, int count);
  static void printVector(std::vector<float> &dat);

  /** 
   * Do a fast fourier transform on data. Values in data vector will be
   * overwritten with the results. The format of data must be that of
   * real and imaginary numbers interleaved together. So if the real
   * portion of the number is at index 0, the imaginary portion will be
   * at index 1. The next real number is at index 2 and imaginary
   * portion is at index 3, etc. The input data vector must be a power
   * of two (i.e. 32, 64, 128, etc) if your data is not a power of two
   * in length you must pad it with zeros. Be warned that this bit of
   * code is not terribly user friendly and can be tough to work with.
   * 
   * @param dat - data of interest.
   * @param inverse - should we perform the inverse fourier transform?
   */
  static void fft(std::vector<float> &dat, bool inverse);

  /** 
   * Convolves the real numbers in two vectors. Fill in kerelwords with
   * the convolved data.
   * 
   * @param yVec - binned up data.
   * @param kernelWords - kernel function.
   * @param count - Number of entries in yVec and kernelWords. 
   */
  static void convolveDataKernel(std::vector<float> &yVec, std::vector<float> &kernelWords, int count);

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
  static void estimateBgParam(std::vector<float> &data, double *mu, 
                              double *sigma, double *alpha, int numBins= 16384);

  /**
   * The "true" signal is then estimated as:
   * 
   * X = a + b(phi(a/b)/PHI(a/b))
   * where a = signal - mu - (sigma^2 * alpha)
   *       b = sigma
   
   *       phi and PHI are standard normal density function and
   *       distribution function respectively given as:
   *
   *       phi(z) = 1 / sqrt(2*pi) * exp(-.5*z^2)
   *       PHI(z) = \int_-inf^z 1 / sqrt(2*pi) exp(-.5w^2) dw
   * @param value original value.
   * @param mu - mode of "background" signal.
   * @param sigma - stdev of "background" signal
   * @param alpha - mean of exponential "true" signal
   * @return bg adjusted value.
   */
  static inline float bgSubData(float value, double mu, double sigma, double alpha, bool &error) {
    //Verbose::out(5,"bgSubData: value " + ToStr(value) + " [in]");
    //Verbose::out(5,"bgSubData: mu    " + ToStr(mu) + " [in]");
    //Verbose::out(5,"bgSubData: sigma " + ToStr(sigma) + " [in]");
    //Verbose::out(5,"bgSubData: alpha " + ToStr(alpha) + " [in]");
    value = (float)(value - mu - alpha*sigma*sigma);
    //Verbose::out(5,"bgSubData: value " + ToStr(value) + " [post transformation]");
    float denom = float(pnorm(value/sigma));
    //Verbose::out(5,"bgSubData: denom " + ToStr(denom));
    if(denom == 0) {
        //Verbose::out(5,"bgSubData: error is set");
        error = true;
    } else 
        error = false;
    value = (float)(value + sigma * phi(value/sigma)/denom);
  return value;
  }

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
  static void backgroundCorrect(std::vector<float> &data, int numBins = 16384);

  /** 
   * Find the median for the rows.
   * 
   * @param psMatrix - matrix.
   * @param numRow - number of rows in matrix.
   * @param numCol - number of rows in matrix.
   * @param rowEffect - fill in row effect.
   */
  static void getRowEffect(std::vector< std::vector<float> > &psMatrix,
                           int numRow, int numCol,
                           std::vector<float> &rowEffect);

  /** 
   * Subtract the median off of each row.
   * 
   * @param psMatrix - matrix.
   * @param numRow - number of rows.
   * @param numCol - number of columns.
   * @param rowEffect - values to subtract from each row.
   */
  static void subRowEffect(std::vector< std::vector<float> > &psMatrix, 
                           int numRow, int numCol,
                           std::vector<float> &rowEffect);

  /** 
   * Find the median for the columns.
   * 
   * @param psMatrix - matrix.
   * @param numRow - number of rows in matrix.
   * @param numCol - number of columns in matrix.
   * @param colEffect - fill in column effect.
   */
  static void getColEffect(std::vector< std::vector<float> > &psMatrix,
                           int numRow, int numCol,
                           std::vector<float> &colEffect);

  /** 
   * Subtract the median off of each column.
   * 
   * @param psMatrix - matrix.
   * @param numRow - number of rows.
   * @param numCol - number of columns.
   * @param colEffect - values to subtract from each column.
   */
  static void subColEffect(std::vector< std::vector<float> > &psMatrix,
                           int numRow, int numCol,
                           std::vector<float> &colEffect);

  /** 
   * Return the sum of a matrix.
   * @param psMatrix - matrix.
   * @param numRow - number of rows in matrix.
   * @param numCol - number of columns in matrix.
   * 
   * @return sum of matrix.
   */
  static double matrixSum(std::vector< std::vector<float> > &psMatrix, int numRow, int numCol);

  /** 
   * Use median polish to get a summary of probe set concentrations
   * using probes that are at the row positions in PM specified by 
   * indexes. results will be a vector with one summary value for each
   * GeneChip for the probe set requested.
   * 
   * @param indexes - Rows that are probes in the probe set of interest.
   * @param results - Vector in which probe sets summaries, one per GeneChip,
   *        will be returned as log2 values.
   */
  void medianPolishPSet(const std::vector<int> &indexes, std::vector<float> &results);

  /** 
   * Standard normal probability density function.
   * @param x - value of interest.
   * @return - density at value supplied. Specifically: 1/sqrt(2*pi)*exp(-0.5*x^2)
   */
  static double phi(double x);

  /** 
   * The cdf function for standard normal distribution consistant with Bolstad's
   * RMA notatation. Also see pnorm.
   * @param x - value.
   * @return - probability that standard normal value has
   * probability <= x.
   */
  static double PHI(double x);

  /** 
   * Cumulative standard normal density function.
   * @param z - value for which to calculate the distribution function.
   * @return - probability that value is < z.
   */
  static double pnorm(double z);

  /** 
   * Get the ranks of _rmaPair structures using the same algorthim as R
   * to break ties.
   * @param pairs - array of _rmaPair strucutres.
   * @param count - size of array.
   * @param ranks - vector of ranks.
   */
  static void getRanks(struct _rmaPair *pairs, int count, std::vector<int> &ranks);
  
  /** 
   * Quantile Norm data using same algorithm to break ties as RMA
   * This function is deprecated as normalization sdk module now supports
   * bioconductor style quantile normalization.
   * @param data - matrix of data.
   */
  static void quantileNorm(std::vector< std::vector<float> > &data);

};
#endif /* RMA_H */
