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

/*! \file statfun.h This file statistical analysis functions.
 */

#ifndef __STATFUN_H_
#define __STATFUN_H_

#include <cmath>
#include <cstdio>
#include <vector>

//

// This should be defined, however it isnt on solaris 5.8
#ifndef NAN
// if NAN is not defined then dont use the "nan()" function
// as it is not likely to be defined either.
// passing NULL to nan on redhat-8 causes a SEGV
// #define NAN nan(NULL)
// Instead we generate our own NAN value.
#ifdef _WIN32 
  #include <limits>
  #define NAN numeric_limits<double>::quiet_NaN()
#else
 #define NAN (0.0/0.0)
#endif
#endif



namespace affxstat {

/*! The square root of 2. */
#define SQRT_TWO 1.4142135623730951455

/*! The log of 2 */
#define LOG_TWO  0.69314718055994528623

/*! Used for comparing two floating point numbers */
#define AFFX_STAT_CONTINUED_FRACTION_EPSILON 1e-15

/*! The maximum number of iterations */
#define AFFX_STAT_CONTINUED_FRACTION_MAX_ITERATION 1000

/*! Used for comparing a floating point value to zero */
#define ZERO_COMPARE_EPSILON 1e-8

/*! Max n for which we can compute psignrank() */
#define PSIGNRANK_MAX_N 1000

/*! Max n for which we can compute pwilcox() */
#define PWILCOX_MAX_N 1000

/*! For n exceeding this the normal approximation to the signed rank test will be used. */
#define APPROX_SIGNRANK_CUTOFF 50   

/*! For sample size exceeding this the normal approximation will be used for signed rank tests */
#define APPROX_RANKSUM_CUTOFF 50   

/*! Tail types for statistical tests. */
typedef enum _TAIL_TYPE {
	/*! One sided lower tail */
	ONE_SIDED_LOWER,

	/*! Once sided upper tail */
	ONE_SIDED_UPPER,

	/*! Two sided tail */
	TWO_SIDED
} TAIL_TYPE;

/*! This function ...
 *
 * @param xx The ...
 * @return The ...
 */
double gammln(double xx);

/*! This function ...
 *
 * @param aa The ...
 * @return The ...
 */
double log_gamma(double aa);

/*! This function ...
 *
 * @param a The ...
 * @param b The ...
 * @return The ...
 */
double log_beta(double a, double b);

/*! This function ...
 *
 * @param a The ...
 * @param b The ...
 * @param x The ...
 * @return The ...
 */
double Incomplete_Beta(double a, double b, double x);

/*! This function ...
 *
 * @param a The ...
 * @param z The ...
 * @param doLog The ...
 * @return The ...
 */
double Incomplete_Gamma(double a, double z, bool doLog);

/*! Calculate the hardy weinburg equilibrium
 * @param nAA The number of AA calls
 * @param nAB The number of AB calls
 * @param nBB The number of BB calls
 * @return The  HW equilibrium
 */
double CalcHWEqPValue(int nAA, int nAB, int nBB);

/*! This function ...
 *
 * @param F The ...
 * @param vone The ...
 * @param vtwo The ...
 * @return The ...
 */
double Ftest(double F, double vone, double vtwo);

/*! This function ...
 *
 * @param x The ...
 * @return The ...
 */
double erf(double x);

/*! This function ...
 *
 * @param x The ...
 * @param mu The ...
 * @param sigma The ...
 * @param lower_tail Flag indicating if lower tail should be computed.
 * @param log_p Flag indicating if the p-value should be logged.
 * @return The ...
 */
double pnorm(double x, double mu, double sigma, bool lower_tail, bool log_p);

/*! This function ...
 *
 * @param w The ...
 * @param nx The ...
 * @param xy The ...
 * @param lower_tail Flag indicating if lower tail should be computed.
 * @param log_p Flag indicating if the p-value should be logged.
 * @return The ...
 */
double pwilcox(unsigned int w, unsigned int nx, unsigned int ny, bool lower_tail, bool log_p);

/*! The implementation of wilcoxon signed rank statistic.
 *
 * @param w The ...
 * @param n The ...
 * @param lower_tail Flag indicating if lower tail should be computed.
 * @param log_p Flag indicating if the p-value should be logged.
 * @return The signed rank statistic.
 */
double psignrank(unsigned int t, unsigned int n, bool lower_tail, bool log_p);

/*! This function compute the N choose K value.
 *
 * @param n The N value
 * @param k The K value
 * @return N choose K
 */
double nChooseK(unsigned int n, unsigned int k);

/*! Performs a Wilcoxon rank sum test.  This is the simplest interface.
 *
 * @param x1 The array of values for the first sample.
 * @param n1 The number of points in the x1 array.
 * @param x2 The array of values for the second sample.
 * @param n2 The number of points in the x2 array.
 * @param tail_type The tail type for the test (one of ONE_SIDED_UPPER, ONE_SIDED_LOWER, TWO_SIDED)
 * @param log_p Boolean flag indicating if the p-value should be logged.
 * @return The computed p-value.
 */
double ranksumTest(double *x1, int n1, double *x2, int n2, int tail_type, bool log_p);


/*! Performs a Wilcoxon rank sum test.  This is the full interface.
 *
 * @param x1 The array of values for the first sample.
 * @param n1 The number of points in the x1 array.
 * @param x2 The array of values for the second sample.
 * @param n2 The number of points in the x2 array.
 * @param tied Boolean flag, value will reflect whether or not there were ties in the input data.
 * @param pval The computed p-value.
 * @param tail_type The tail type for the test (one of ONE_SIDED_UPPER, ONE_SIDED_LOWER, TWO_SIDED)
 * @param log_p Boolean flag indicating if the p-value should be logged.
 */
void ranksumTest(double *x1, int n1, double *x2, int n2, bool *tied, double *pval, int tail_type, bool log_p);


/*! Performs a Wilcoxon rank sum test.  This is the legacy interface.
 *
 * @param x1 The array of values for the first sample.
 * @param n1 The number of points in the x1 array.
 * @param x2 The array of values for the second sample.
 * @param n2 The number of points in the x2 array.
 * @param nTied Value will be incremented by one if there were ties in the input data.
 * @param pval The computed p-value.
 * @param signal The median of all pairwise differences in x1 and x2 (Hodges-Lehmann estimator)
 * @param tail_type The tail type for the test (one of ONE_SIDED_UPPER, ONE_SIDED_LOWER, TWO_SIDED)
 */
void ranksumTest(double *x1, int n1, double *x2, int n2, int *nTied, double *pval, double *signal, int tail_type);


/*! Performs a Wilcoxon signed rank test.  This is the simplest interface.
 *
 * @param x The array of values to analyze.
 * @param n The number of points in the x array.
 * @param tail_type The tail type for the test (one of ONE_SIDED_UPPER, ONE_SIDED_LOWER, TWO_SIDED)
 * @param log_p Flag indicating if the p-value should be logged.
 * @return The computed p-value.
 */
double signedRankTest(double *x, int n, int tail_type, bool log_p);

/*! Performs a signed rank test.  This is the full interface.
 *
 * @param x The array of values to analyze.
 * @param n The number of points in the x array.
 * @param tied The number of observations that were tied.
 * @param nZero The number of observations that were zero.
 * @param tiedAndZero The number of observations that were tied and zero.
 * @param pval The p-value.
 * @param signal The pseudo median.
 * @param tail_type The tail type for the test (one of ONE_SIDED_UPPER, ONE_SIDED_LOWER, TWO_SIDED)
 * @param log_p Flag indicating if the p-value should be logged.
 */
void signedRankTest(double *x, int n, bool *tied, bool *nZero, bool *tiedAndZero, double *pval, int tail_type, bool log_p);

/*! Performs a signed rank test, returns logged lower-tail p-value and computes pseudoMedian.  This is the legacy interface.
 *
 * @param x The array of values to analyze.
 * @param n The number of points in the x array.
 * @param tied The number of observations that were tied.
 * @param nZero The number of observations that were zero.
 * @param tiedAndZero The number of observations that were tied and zero.
 * @param pval The p-value.
 * @param signal The pseudo median.
 */
void signedRankTest(double *x, int n, int *tied, int *nZero, int *tiedAndZero, double *pval, double *signal);

/*! Computes the pseudo median.
 *
 * @param x The array of values for which to compute the pseudomedian.
 * @param n The number of points in the x array.
 * @return The pseudo median.
 */
double pseudoMedian(double *x, int n);

/*! Computes the median over all pairwise differences in two arrays of values.
 *
 * @param x  The first array of values.
 * @param nX The number of points in the x1 array.
 * @param y  The second array of values.
 * @param nY The number of points in the x2 array.
 * @return The median of all n1*n2 pairwise differences.
 */
double medianDifference(double *x, int nX, double *y, int nY);

/*! This function calculates Pearson's correlation coefficient between 2 arrays of floats
 *
 * @param x1 The first array of floats
 * @param x2 The second array of floats
 * @param nPoints The size of each array. assumed the same
 * @param transform A function to transform the data
 * @return Pearson's correlation
 */
double PearsonCorrelation(float *x1, float *x2, int nPoints, float(*transform)(float x) = NULL);

}

#endif // __STATFUN_H_
