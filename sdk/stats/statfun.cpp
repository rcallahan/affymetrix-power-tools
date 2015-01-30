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

// C99 defines NAN.  However if NAN is not defined, then generate one.
#include "stats/statfun.h"
//
#include "stats/stats-distributions.h"
//
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
//
#ifdef WIN32
#include <math.h>
#endif
#ifndef NAN
#define NAN (0.0/0.0)
#endif

using namespace std;
using namespace affxstat;

double continued_fraction_recurrence(double *a, double *b, double n, double d);
double Incomplete_Beta_Continued_Fraction(double a, double b, double x);
double Incomplete_Gamma_Continued_Fraction(double a, double z, bool doLog);
double Incomplete_Gamma_Power_Series(double a, double z, bool doLog);


typedef struct {
  double val;
  int index;
  int vector;
  double rank;
} Rank;
int rankCmp(Rank *a, Rank *b);
int doubleCmp(double *a, double *b);
double ranksum_nways(unsigned int n, unsigned int t);
double ranksum_nways(int n, int m, int t);
double approxSignrank(double wilcoxStat, int n, int sumTies, int tail_type, bool log_p);

double approxRankSum(double rankSumStat, int n, int m, int sumTies, int tail_type, bool log_p);
void rankSum(double *x, int nX, double *y, int nY, double *w, int *sumTies);
int multi_rank(int m, double **f, int *n, double **r);
int computeRanks(Rank *rs, int n);
double choose(unsigned int n, unsigned int k);


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

		if ((minr <0 || minr > maxr || maxr >= nRows()) || 
				(minc <0 || minc > maxc || maxc >= nCols())) {
			throw "average";
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
	void show() {
		for (int i = 0; i < nRows(); ++i) {
			printVector(array[i]);
		}
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

template<class T> double mean(const vector<T> & v) {
	if (v.empty()) {
		throw "mean";
	}
	double sum = 0;
	for (int i = 0; i < v.size(); ++i) {
		sum += v[i];
	}
	return sum / (double) v.size();
}

template<class T> double sumSquareDiff(const vector<T> & v, const double mv) {
	if (v.empty()) {
		throw "sumSquareDiff";
	}
	double sum = 0;
	for (int i = 0; i < v.size(); ++i) {
		double difference = v[i] - mv;
		sum += difference * difference;
	}
	return sum;
}

double affxstat::gammln(double xx)
{
        double x,y,tmp,ser;

        y = x = xx;
        tmp = x+5.5;
        tmp -= (x+.5)*log(tmp);
        ser = 1.000000000190015;
        ser += 76.18009172947146/++y;
        ser += -86.5053/++y;
        ser += 24.014/++y;
        ser += -1.23/++y;
        ser += -0.0012/++y;
        return((double)(-1*tmp+log(2.506628*ser/x)));
}

double affxstat::log_gamma(double aa)
{
        double retval;
        double series;
        double a = aa;
        double g = 5;
        double z;
        double t;

        // compute log_gamma(aa)) by Lanczos approximation
        // order 6
        // log_gamma(a) = log((sqrt(2*pi)/a)*(x_0+x_1/(a+1)+...)(a+g+.5)^(a+.5)exp(-(a+g+.5))
        // is this even correct?
        series = 1.000000000190015;
        series +=  76.18009172947146/(a+1);
        series += -86.50532032941677/(a+2);
        series +=  24.01409824083091/(a+3);
        series += -1.231739516/(a+4);
        series +=  0.0012058003/(a+5);
        series += -0.00000536382/(a+6);
        series *= 2.50662827465/a; // sqrt(2*pi)/a
        z = a+g+.5;
        t = z-(a+.5)*log(z);
        retval = log(series)-t;
        return(retval);
}

double continued_fraction_recurrence(double *a, double *b, double n, double d)
{
        long i;
        double tmp;
        // do steps of continued fraction
        // shift right
        a[2] = a[1];
        a[1] = a[0];
        b[2] = b[1];
        b[1] = b[0];
        // n/d
        a[0] = a[1]*d + a[2]*n;  // numerator times second in past, denominator times past
        b[0] = b[1]*d + b[2]*n;  // numerator times second in past, denominator times past
        // what if b[0] = 0?  how do I fix?
        tmp = b[0];
        for (i=0; i<3; i++)
        {
                a[i] /= tmp;
                b[i] /= tmp;
        }
        return(a[0]); // current value of function
}

double affxstat::log_beta(double a, double b)
{
        double t, s;
        t = affxstat::gammln(a)+ affxstat::gammln(b)- affxstat::gammln(a+b);
        s = affxstat::log_gamma(a)+ affxstat::log_gamma(b)- affxstat::log_gamma(a+b);
        return(s);
}

double Incomplete_Gamma_Continued_Fraction(double a, double z, bool doLog) {
  // compute Q_(a,z)
  // as [z^a exp(-z) / Gamma(a)] * 1/(z+(1-a)/1+1/z+(2-a)/1+2/z+(3-a)/...)
  double Pn[3], Qn[3];
  int j;
  double old,numerator,denominator;

  // start one step into the recurrence
  Pn[1] = 1;
  Pn[0] = 0;
  Qn[1] = 0;
  Qn[0] = 1;

  numerator   = 1;
  denominator = z;
  continued_fraction_recurrence(Pn, Qn, numerator, denominator);
  for(j=2; j < AFFX_STAT_CONTINUED_FRACTION_MAX_ITERATION; j++) {
    old = 1/Pn[0]; // old value
    if(j % 2 == 0) {
      numerator   = (j/2)-a;
      denominator = 1;
    } else {
      numerator   = (j-1)/2;
      denominator = z;
    }
    continued_fraction_recurrence(Pn, Qn, numerator, denominator);
    if (fabs(1-Pn[0]*old) < AFFX_STAT_CONTINUED_FRACTION_EPSILON)
      break;
  }

  // now we've done the continued fraction, multiply in the initial term
  double ans = a*log(z)-z- affxstat::log_gamma(a)+log(Pn[0]);
  if(!doLog)
    ans = exp(ans);

  return(ans);
}

double Incomplete_Gamma_Power_Series(double a, double z, bool doLog) {
  // compute Q_(a,z)
  // as 1 - [ z^a / Gamma(a) * (1/a - z/(a+1) + z^2/(2(a+2)) + ...) ]

  double oldsum,sum,numerator,denominator,multiplier,log_z;
  numerator = 0;
  multiplier = 1;
  log_z = log(z);
  sum=log(1/a);
  int i;
  for(i=1; i < AFFX_STAT_CONTINUED_FRACTION_MAX_ITERATION; i++) {
    oldsum = sum;
    numerator += log_z;
    denominator = log((double)i)+log(a+i);
    multiplier *= -1;
    sum += log(1 + multiplier * exp(numerator-denominator-sum)); // trying to do it in an overflow-avoiding manner
    if (fabs(sum-oldsum) < AFFX_STAT_CONTINUED_FRACTION_EPSILON)
      break;
  }
  
  double ans = 1-exp(a*log(z) -  affxstat::log_gamma(a) + sum);
  if(doLog)
    ans = log(ans);

  return(ans);
}

double affxstat::Incomplete_Gamma(double a, double z, bool doLog)
{
  double ans;

  if (z<=0) {
    ans = 1;
  } else if(z > 0.03) {
    // use continued fraction
    return(Incomplete_Gamma_Continued_Fraction(a,z,doLog));
  } else {
    // use power series
    return(Incomplete_Gamma_Power_Series(a,z,doLog));
  }

  return(ans);
}

double Incomplete_Beta_Continued_Fraction(double a, double b, double x)
{
        // compute I_x(a,b)
        // as [(x^a)*(1-x)^b /(a*Beta(a,b))]*
        // 1/(1+c1/1+c2/...) - continued fraction expansion
        // c_2n = n(b-n)x/(a+2n-1)*(a+2n)
        // c_2n+1 = -(a+n)(a+b+n)x/(a+2n-1)(a+2n)
        // take log to fix up...

        double an[3], bn[3];
        long j;
        long MaxJ = 1000;
        double epsilon = AFFX_STAT_CONTINUED_FRACTION_EPSILON;
        double old;
        double coef;
        // start one step into the recurrence
        an[1] = 1;
        bn[1] = 1;
        an[0] = 1-((a+b)*x)/(a+1);
        bn[0] = 1;
        // use the continued fraction approximation to incomplete beta function
        // note - computing 1/(1+(c1)/1+c2/1+...)
        // where d (2m+1) = -(a+m)(a+b+m)x/(a+2m)(a+2m+1)
        // where d (2m)   = m(b-m)x/(a+2m-1)(a+2m)

        j=1;
        while (j<MaxJ)
        {
                old = 1/an[0]; // old value
                // 2j j = 1 continues
                // even

                coef = j*(b-j)*x/((a+2*j-1)*(a+2*j));
                continued_fraction_recurrence(an, bn, coef, 1);
                // 2j+1 j= 0 starts
                // odd

                coef = -1*(a+j)*(a+b+j)*x/((a+2*j)*(a+2*j+1));
                continued_fraction_recurrence(an, bn, coef , 1);
                j++;  // update j here
                if (fabs(1-an[0]*old)<epsilon)
                        break;
        }      
        // now we've done the continued fraction
        // do the coefficient
        // exp(a*log(x)+b*log(1-x) +log_gamma(a+b)-log_gamma(a)-log_gamma(b))/a

        return(exp(a*log(x)+b*log(1-x)- affxstat::log_beta(a,b))/(a*an[0]));
}

double affxstat::Incomplete_Beta(double a, double b, double x)
{
        if (x<=0)
                return(0);
        if (x>=1)
                return(1); // limiting cases!
        // check for convergence criteria
        if (x*(a+b+2)<(a+1))
                return(Incomplete_Beta_Continued_Fraction(a,b,x));
        else
                return(1-Incomplete_Beta_Continued_Fraction(b,a,1-x));
}

double affxstat::Ftest(double F, double vone, double vtwo)
{
        double a,b,x,s;
        x = vtwo/(vtwo+vone*F);
        a = vtwo/2;
        b = vone/2;
        s= affxstat::Incomplete_Beta(a,b,x);
        return(s);
}


// Implementation of erf taken directly from MAS5 codebase.
double affxstat::erf(double x) {
	static const int erfCoeffSize[] = {4, 8, 5};
	static const double erfA1[] = {1.85777706184603153e-1,
		3.16112374387056560,    1.13864154151050156e2,
		3.77485237685302021e2,  3.20937758913846947e3};
    static const double erfB1[] = {2.36012909523441209e1,
		2.44024637934444173e2,  1.28261652607737228e3,
		2.84423683343917062e3};
	static const double erfA2[] = {2.15311535474403846e-8,
		5.64188496988670089e-1, 8.88314979438837594,
		6.61191906371416295e1,	2.98635138197400131e2,
		8.81952221241769090e2,	1.71204761263407058e3,
		2.05107837782607147e3,	1.23033935479799725e3};
    static const double erfB2[] = {1.57449261107098347e1,
		1.17693950891312499e2,  5.37181101862009858e2,
		1.62138957456669019e3,  3.29079923573345963e3,
		4.36261909014324716e3,  3.43936767414372164e3,
		1.23033935480374942e3};
	static const double erfA3[] ={1.63153871373020978e-2,
		3.05326634961232344e-1, 3.60344899949804439e-1,
		1.25781726111229246e-1, 1.60837851487422766e-2,
		6.58749161529837803e-4};
    static const double erfB3[] = {2.56852019228982242,
		1.87295284992346047,    5.27905102951428412e-1,
		6.05183413124413191e-2, 2.33520497626869185e-3};
	static const double erfXBounds[] = {0.46875, 4.0};
	static const double invSqrtPi = 0.56418958354775627928; // 1 / sqrt(pi)
	
	double absX = fabs(x);
	double xSquared = absX * absX;
	double temp = 0;
	

	if (absX <= erfXBounds[0]) {    
		double num = erfA1[0] * xSquared;
		double den = xSquared;
		int last = erfCoeffSize[0];
		for (int i = 1; i < last; i++) {
			num = (num + erfA1[i]) * xSquared;
			den = (den + erfB1[i - 1]) * xSquared;
		}
		return x * (num + erfA1[last]) / (den + erfB1[last - 1]);
	}
	else {
		if (absX <= erfXBounds[1]) {
			double num = erfA2[0] * absX;
			double den = absX;
			int last = erfCoeffSize[1];
			for (int i = 1; i < last; i++) {
				num = (num + erfA2[i]) * absX;
				den = (den + erfB2[i - 1]) * absX;
			}
			temp = (num + erfA2[last]) / (den + erfB2[last - 1]);
		}
		else {
			double xInvSquared = 1.0 / xSquared;
			double num = erfA3[0] * xInvSquared;
			double den = xInvSquared;
			int last = erfCoeffSize[2];
			for (int i = 1; i < last; i++) {
				num = (num + erfA3[i]) * xInvSquared;
				den = (den + erfB3[i - 1]) * xInvSquared;
			}
			temp = xInvSquared * (num + erfA3[last]) / (den + erfB3[last - 1]);
			temp = (invSqrtPi - temp) / absX;
		}
		temp = 1 - exp(-xSquared) * temp;
		if (x > 0)      // in fact, we may use if (x > erfXBounds[0])
			return temp;
		else
			return -temp;
	}               
}


double affxstat::pnorm(double x, double mu, double sigma, bool lower_tail, bool doLog) {
  double prob,e;

  x = (x-mu)/(SQRT_TWO*sigma);
  e = ((x > 0) ?  affxstat::erf(x) :  affxstat::erf(-x))/2.0;

  if((lower_tail && x>0) || (!lower_tail && x<0))
    prob = 0.5 + e;
  else
    prob = 0.5 - e;
  if(doLog)
    prob = log(prob);
  return(prob);
}


//
// Returns nways[n][t] which is defined as the number of ways you can get a rank sum of t from integers (1,...,n).
// Computed using the recursion nways[n][t] = nways[n-1][t] + nways[n-t][t]
//
double ranksum_nways(unsigned int n, unsigned int t) {
	static bool first=true;
	static double **nways=NULL;

	if(first) {
		first=false;
		nways = new double*[PSIGNRANK_MAX_N];
		for(int i=0; i<PSIGNRANK_MAX_N; i++)
			nways[i] = NULL;
	}

	unsigned int ranksum_max = n*(n+1)/2;
        unsigned int ranksum_halfmax = (ranksum_max/2);

	if( (n < 0) || (n >= PSIGNRANK_MAX_N) ) {
		throw("Invalid n specified for ranksum_nways\n");
		return(0);
	}
	if( (t < 0) || (t > ranksum_max) )
		return(0);

	if(nways[n] == NULL) {
		nways[n] = new double[1+ranksum_halfmax];
		for(unsigned int i=0; i <= ranksum_halfmax; i++)
			nways[n][i] = -1;
	}
	if(t > ranksum_halfmax)
		t = ranksum_max-t;

	if(nways[n][t] < 0) {
		// need to initialize it
		if(n==0)
			nways[n][t] = 1;
		else
			nways[n][t] = ranksum_nways(n-1,t) + ranksum_nways(n-1,t-n);
	}

	return(nways[n][t]);
}

//
// Returns nways[n1][n2][t] which is defined as the number of ways you can get a rank sum
// of t from a sample of size n1 sorted with a sample of size n2.
//
// Computed using the recursion nways[n1][n2][t] = nways[n1-1][n2][t-n1] + nways[n-t][t]
//
double ranksum_nways(int n1, int n2, int t) {
	static bool first=true;
	static double ***nways=NULL;

	if(first) {
		first=false;
		nways = new double**[PWILCOX_MAX_N];
		for(int i=0; i<PWILCOX_MAX_N; i++) {
			nways[i] = new double*[PWILCOX_MAX_N];
			for(int j=0; j<PWILCOX_MAX_N; j++)
				nways[i][j] = NULL;
		}
	}

	int ranksum_max = n1*n2 + n1*(n1+1)/2;

	if( (n1 < 0) || (n1 >= PSIGNRANK_MAX_N) )
		throw("Invalid n1 specified for ranksum_nways\n");
	if( (n2 < 0) || (n2 >= PSIGNRANK_MAX_N) )
		throw("Invalid n2 specified for ranksum_nways\n");
	if( (t < 0) || (t > ranksum_max) )
		return(0);

	if(nways[n1][n2] == NULL) {
		nways[n1][n2] = new double[1+ranksum_max];
		for(int i=0; i <= ranksum_max; i++)
			nways[n1][n2][i] = -1;
	}

	if(nways[n1][n2][t] < 0) {
		// need to initialize it
		if(n1==0) {
			nways[n1][n2][t] = (t==0);
		} else if(n2==0) {
			int x = n1*(n1+1)/2;
			nways[n1][n2][t] = (t==x);
		} else {
			nways[n1][n2][t] = ranksum_nways(n1,n2-1,t) + ranksum_nways(n1-1,n2,t-n1-n2);
		}
	}

	return(nways[n1][n2][t]);
}

//
// Efficient implementation of wilcoxon rank sum statistic.
//
double affxstat::pwilcox(unsigned int t, unsigned int n1, unsigned int n2, bool lower_tail, bool log_p)
{
	if( n1 < 1) {
		throw("n1 must be positive for pwilcox\n");
                //		return(NAN);
	}
	if( n1 >= PWILCOX_MAX_N) {
		throw("n1 too large for psignrank\n");
                //		return(NAN);
	}
	if( n2 < 1) {
		throw("n2 must be positive for pwilcox\n");
                //		return(NAN);
	}
	if( n2 >= PWILCOX_MAX_N) {
		throw("n2 too large for psignrank\n");
                //		return(NAN);
	}
		
	unsigned int ranksum_max = n1*n2 + (n1*(n1+1))/2;
	if(t <= 0) {
		double pval = lower_tail ? 0 : 1;
		if(log_p)
			pval = log(pval);
		return(pval);
	}
	else if (t >= ranksum_max) {
		double pval = lower_tail ? 1 : 0;
		if(log_p)
			pval = log(pval);
		return(pval);
	}

	double divisor = choose(n1+n2,n1);
	double pval = 0;
	if(t <= (unsigned int) ranksum_max/2) {
		// shorter to compute lower tail, flip if necessary
		for(int i=1; i <= (int)t; i++)
			pval += ranksum_nways((int)n1,(int)n2,i) / divisor;
		if(!lower_tail)
			pval = 1-pval;
	} else {
		// shorter to compute upper tail, flip if necessary
		for(unsigned int i=t+1; i <= ranksum_max; i++)
			pval += ranksum_nways((int)n1,(int)n2,i) / divisor;
		if(lower_tail)
			pval = 1-pval;
	}

	if(log_p)
		pval = log(pval);
		
	return(pval);
}


//
// Efficient implementation of wilcoxon signed rank statistic.
//
double affxstat::psignrank(unsigned int t, unsigned int n, bool lower_tail, bool log_p)
{
	if( n < 1) {
		throw("n must be positive for psignrank\n");
                //		return(NAN);
	}
	if( n >= PSIGNRANK_MAX_N) {
		throw("n too large for psignrank\n");
		// return(NAN);
	}
		
	unsigned int ranksum_max = n*(n+1)/2;
	if(t < 0)
		return(0);
	else if (t >= ranksum_max)
		return(1);

	double divisor = exp(n*LOG_TWO);    // 2^n
	double pval = 0;
	if(t <= (unsigned int) ranksum_max/2) {
		// shorter to compute lower tail, flip if necessary
		for(unsigned int i=0; i <= t; i++)
			pval += ranksum_nways(n,i) / divisor;
		if(!lower_tail)
			pval = 1-pval;
	} else {
		// shorter to compute upper tail, flip if necessary
		for(unsigned int i=t+1; i <= ranksum_max; i++)
			pval += ranksum_nways(n,i) / divisor;
		if(lower_tail)
			pval = 1-pval;
	}

	if(log_p)
		pval = log(pval);
		
	return(pval);
}


double affxstat::nChooseK(unsigned int n, unsigned int k)
{
  if ((n<0) || (k < 0 || k > n))
    return 0;
  
  unsigned int k1 = (k <= n/2)? k : n - k;
  double ans = 0;
  for (unsigned int i = 0; i < k1; ++i)
    ans += log((double)(n-i))-log((double)(i+1));
  
  return(exp(ans));
}



int rank(double *f, int n, double *r);
void signrank(double *f, int n, double *w, int *m, int *sumTies);

//
//  This interface is being preserved for legacy reasons.  You are probably better off using one
//  of the other interfaces to signedRankTest.
//
//  Performs a signed rank test, returns logged lower-tail p-value and computes pseudoMedian
//
void affxstat::signedRankTest(double *x, int n, int *tied, int *zero, int *tiedAndZero, double *pval, double *signal) {
  bool tied_b,zero_b,tiedAndZero_b;
  signedRankTest(x,n,&tied_b,&zero_b,&tiedAndZero_b,pval,ONE_SIDED_LOWER,true);
  if(tied_b)
    (*tied)++;
  if(zero_b)
    (*zero)++;
  if(tiedAndZero_b)
    (*tiedAndZero)++;
  *signal = pseudoMedian(x,n);
  return;
}

//
// This interface is for those who don't care to find out how many observations were tied, zero or both.
//
double affxstat::signedRankTest(double *x, int n, int tail_type, bool log_p) {
  double pval;
  bool tied,zero,tiedAndZero;
  signedRankTest(x,n,&tied,&zero,&tiedAndZero,&pval,tail_type,log_p);
  return(pval);
}

//
// Wilcoxon signed rank test (with continuity correction)
//
// Inputs:
//   x: pointer to an array of n doubles
//   n: size of the array of doubles
//   tail_type: one of (ONE_SIDED_UPPER, ONE_SIDED_LOWER, TWO_SIDED)
//   log_p: should the p-value be logged?
// Outputs:
//   tied: will be true if any of the observations are tied
//   zero: will be true if any of the observations are zero
//   tiedAndZero: will be true if at least one observation is zero and at least two observations are tied
//   pval: the p-value, with lower/upper status determined by lower_tail and log/linear status by log_p
//
void affxstat::signedRankTest(double *x, int n, bool *tied, bool *zero, bool *tiedAndZero, double *pval, int tail_type, bool log_p) {
  double wilcoxStat;
  int nonZero,sumTies;

  signrank(x,n,&wilcoxStat,&nonZero,&sumTies);
  *tiedAndZero = *tied = *zero = false;
  if(sumTies && (nonZero < n))
    *tiedAndZero = true;
  else if(sumTies)
    *tied = true;
  else if(nonZero < n)
    *zero = true;

  if(*zero || *tied || *tiedAndZero || n >= APPROX_SIGNRANK_CUTOFF) {
    // approximate test
    *pval = approxSignrank(wilcoxStat,nonZero,sumTies,tail_type,log_p);
  } else {
    // exact test
    switch(tail_type) {
      case ONE_SIDED_LOWER:
        *pval = psignrank((unsigned int) wilcoxStat, (unsigned int) nonZero,true,log_p);
        break;
      case ONE_SIDED_UPPER:
        *pval = psignrank((unsigned int) (wilcoxStat-1.0), (unsigned int) nonZero,false,log_p);
        break;
      case TWO_SIDED:
	double median;
	median = (nonZero*(nonZero+1)/4.0);
        if(wilcoxStat > median)
          *pval = psignrank((unsigned int) (wilcoxStat-1.0), (unsigned int) nonZero,false,log_p);
        else if(wilcoxStat < median)
          *pval = psignrank((unsigned int) wilcoxStat, (unsigned int) nonZero,true,log_p);
        else
          *pval = (log_p) ? -LOG_TWO : 0.5;
        if(log_p)
          *pval += LOG_TWO;
        else
          *pval *= 2.0;
        break;
      default:
        throw("Invalid tail type specified!\n");
        break;
    }
  }

  return;
}

/*
** Computes wilcoxon signed rank sum.  The answer is stored in *w and the number of
** data which went into its calculation is stored in *m.  sumTies will contain the sum
** of nTies^3 - nTies over all tied groups.
*/
void signrank(double *f,int n, double *w, int *m, int *sumTies) {
  int i,j,*isPos;
  double *absData,*ranks;

  absData = (double *) malloc(n * sizeof(double));
  isPos = (int *) malloc(n * sizeof(int));

  /* Get absolute values of non-zero entries. */
  for(i=0,j=0; i<n; i++) {
    if(fabs(f[i]) > ZERO_COMPARE_EPSILON) {
      absData[j] = fabs(f[i]);
      isPos[j] = (f[i] > 0);
      j++;
    }
  }
  *m = j;

  /* Compute ranks. */
  ranks = (double *) malloc(*m * sizeof(double));
  *sumTies = ::rank(absData,*m,ranks);

  /* Get rank sum for positive entries. */
  for(j=0,*w=0; j<*m; j++) {
    if(isPos[j])
      *w += ranks[j];
  }

  /* Tidy up. */
  free(isPos);
  free(absData);
  free(ranks);

  return;
}


//
// Computes the ranks of an array of doubles.  Returns sumTies, which is the sum of nTies^3-nTies for each tied group. */
//
int rank(double *f,int n,double *r) {
  double rankSum,thisRank,oldVal;
  int i,j,nRanks,last;
  int sumTies=0;
  Rank *rs;

  /* Initialise the rank structure. */
  rs = (Rank *) malloc(n * sizeof(Rank));
  for(i=0; i<n; i++) {
    rs[i].val = f[i];
    rs[i].index = i;
  }

  /* Sort the doubles. */
  qsort((void *)rs, n, sizeof(Rank), (int (*)(const void*, const void*))(rankCmp));

  /* Assign ranks. */
  for(i=0, rankSum=0, last=0, nRanks=0, oldVal=0; i<n; i++) {
    if((i == 0) || (fabs(rs[i].val - oldVal) < ZERO_COMPARE_EPSILON)) {
      rankSum += i+1;
      nRanks++;
    } else {
      sumTies += (nRanks * nRanks * nRanks) - nRanks;
      // printf("nRanks is %d, sumTies is now %d\n",nRanks,sumTies);
      thisRank = rankSum/nRanks;
      for(j=last; j<i; j++)
        rs[j].rank = thisRank;
      last = i;
      rankSum = i+1;
      nRanks = 1;
    }
    oldVal = rs[i].val;
  }
  sumTies += (nRanks * nRanks * nRanks) - nRanks;
  thisRank = rankSum/nRanks;
  while(last < n)
    rs[last++].rank = thisRank;

  /* Copy ranks to output array. */
  for(i=0; i<n; i++)
    r[rs[i].index] = rs[i].rank;

  free(rs);
  return(sumTies);
}


double approxSignrank(double wilcoxStat, int n, int sumTies, int tail_type, bool log_p) {
  double z,sigma;
  double logTerm1,logTerm2;

  z = wilcoxStat - n*(n+1)/4.0;
  if(n <= 1000) {
    sigma = sqrt(n*(n+1)*(2*n+1)/24.0 - sumTies/48.0);
  } else {
    /* For larger n, we compute sigma the following more complex way to avoid overflow */
    if(sumTies > 0) {
      logTerm1 = log((double)n)+log((double)(n+1))+log((double)(2*n+1)) - 3.1780538303479457518;
      logTerm2 = log((double)sumTies) - 3.8712010109078911491;
      sigma = exp(0.5*(logTerm1 + log(1-exp(logTerm2-logTerm1))));
    } else {
      sigma = exp(0.5*(log((double)n)+log((double)(n+1))+log((double)(2*n+1))-3.1780538303479457518));
    }
  }
  
  // apply continuity correction
  double correction = 0;
  switch(tail_type) {
    case ONE_SIDED_LOWER:
      correction = 0.5;
      break;
    case ONE_SIDED_UPPER:
      correction = -0.5;
      break;
    case TWO_SIDED:
      correction = (z > 0) ? -0.5 : 0.5;
      break;
    default:
      throw("Invalid tail_type specified!");
      break;
  }

  double pval;
  z = z + correction;
  if(tail_type == ONE_SIDED_UPPER) {
    pval = pnorm(z,0,sigma,false,log_p);
  } else if (tail_type == ONE_SIDED_LOWER) {
    pval = pnorm(z,0,sigma,true,log_p);
  } else if (tail_type == TWO_SIDED) {
    pval = (z > 0) ? pnorm(z,0,sigma,false,log_p) : pnorm(z,0,sigma,true,log_p);
    if(log_p)
      pval += LOG_TWO;
    else
      pval *= 2.0;
  } else {
    throw("Invalid tail_type specified!");
  }

  return(pval);
}


double affxstat::pseudoMedian(double *x, int n) {
  double *walshAvg=NULL;
  double pMed;
  int nW,i,j,k;

  if(n <= 0) {
    fprintf(stderr,"ERROR: pseudoMedian: empty dataset supplied.\n");
    return(0);
  } else if(n == 1) {
    return(x[0]);
  } else {
    nW = n * (n+1) / 2;
    walshAvg = (double *) malloc(nW * sizeof(double));
    for(i=0,k=0; i<n; i++)
      for(j=i; j<n; j++)
        walshAvg[k++] = (x[i] + x[j]) / 2.0;
    qsort((void *)walshAvg, nW, sizeof(double), (int (*)(const void*, const void*))(doubleCmp));
    if((nW % 2) == 0) {
      pMed = 0.5 * (walshAvg[nW/2] + walshAvg[nW/2 - 1]);
    } else {
      pMed = walshAvg[nW/2];
    }
    free(walshAvg);
    return(pMed);
  }
}

int rankCmp(Rank *a, Rank *b) {
  if( a->val > b->val)
    return(1);
  else if( a->val < b->val)
    return(-1);
  else
    return(0);
}

int doubleCmp(double *a, double *b) {
  if( *a > *b)
    return(1);
  else if( *a < *b)
    return(-1);
  else
    return(0);
}


//
//  This interface is being preserved for legacy reasons.  You are probably better off using one
//  of the other interfaces to ranksumTest.
//
//  Performs a rank sum test, returns logged lower-tail p-value and computes median of all pairwise differences.
//
void affxstat::ranksumTest(double *x1, int n1, double *x2, int n2, int *nTied, double *pval, double *signal, int tail_type) {
  bool tied;

  // Do the test
  ranksumTest(x1,n1,x2,n2,&tied,pval,tail_type,false);

  // Check if there were ties.
  if(tied)
    (*nTied)++;

  // compute HL estimator
  *signal = medianDifference(x1,n1,x2,n2);

  return;
}



//
// Performs a Wilcoxon signed rank test, returns the computed p-value.
//
double affxstat::ranksumTest(double *x1, int n1, double *x2, int n2, int tail_type, bool log_p) {
  double pval=0;
  bool tied;
  ranksumTest(x1,n1,x2,n2,&tied,&pval,tail_type,log_p);
  return(pval);
}



//
// Full interface to Wilcoxon rank sum test.  Returns computed p-value, reports back if there
// were any ties in the data.
//
void affxstat::ranksumTest(double *x, int nX, double *y, int nY, bool *tied, double *pval, int tail_type, bool log_p) {
  double rankSumStat;
  int sumTies=0;

  rankSum(x,nX,y,nY,&rankSumStat,&sumTies);
  *tied = (sumTies != 0) ? true : false;
  if(*tied || nX >= APPROX_RANKSUM_CUTOFF || nY >= APPROX_RANKSUM_CUTOFF) {
    /* approximate test */
    *pval = approxRankSum(rankSumStat,nX,nY,sumTies,tail_type,log_p);
  } else {
    /* exact test */
    switch(tail_type) {
      case ONE_SIDED_LOWER:
        *pval = pwilcox((unsigned int) rankSumStat, (unsigned int) nX, (unsigned int) nY, true, log_p);
        break;
      case ONE_SIDED_UPPER:
        *pval = pwilcox((unsigned int) (rankSumStat-1), (unsigned int) nX, (unsigned int) nY, false, log_p);
        break;
      case TWO_SIDED:
	double median;
	median = ((nX*nY + nX*(nX+1))/2.0);
        if(rankSumStat > median)
          *pval = pwilcox((unsigned int) (rankSumStat-1), (unsigned int) nX, (unsigned int) nY, false, log_p);
        else if(rankSumStat < median)
          *pval = pwilcox( (unsigned int) rankSumStat, (unsigned int) nX, (unsigned int) nY, true, log_p);
	else
	  *pval = (log_p) ? -LOG_TWO : 0.5;
        if(log_p)
          *pval += LOG_TWO;
        else
          *pval *= 2.0;
        break;
      default:
        throw("Invalid tail type specified!\n");
        break;
    }
  }

  return;
}

double approxRankSum(double rankSumStat, int n, int m, int sumTies, int tail_type, bool log_p) {
  double nm,npm,z,pval,sigma,logTerm1,logTerm2;

  nm = n*m;
  npm = n+m;

  z = rankSumStat - (nm + n*(n+1))/2.0;

  /* compute sigma */
  if(n <= 1000 && m <= 1000) {
    sigma = sqrt(nm * (npm + 1 - sumTies/(npm * (npm-1))) / 12.0);
  } else {
    if(sumTies > 0) {
      logTerm1 = log(nm)+log(npm+1) - 2.4849066497880003546;
      logTerm2 = log(nm)+log((double)sumTies) - log(npm) - log(npm+1) - 2.4849066497880003546;
      sigma = exp(0.5*(logTerm1 + log(1-exp(logTerm2-logTerm1))));
    } else {
      sigma = exp(0.5*(log(nm)+log(npm+1) - 2.4849066497880003546));
    }
  }

//printf("z is %f, sigma is %f, sumTies is %d\n",z,sigma,sumTies);
  /* apply continuity correction */
  if(tail_type == ONE_SIDED_UPPER) {
    z -= 0.5;
  } else if(tail_type == ONE_SIDED_LOWER) {
    z += 0.5;
    pval = pnorm(z,0,sigma,ONE_SIDED_LOWER,log_p);
  } else {
    if(z < 0) {
      z += 0.5;
    } else if (z > 0) {
      z -= 0.5;
    }
  }
//printf("corrected z is %f\n",z);

  if(tail_type == ONE_SIDED_UPPER) {
    pval = pnorm(z,0,sigma,false,log_p);
  } else if (tail_type == ONE_SIDED_LOWER) {
    pval = pnorm(z,0,sigma,true,log_p);
  } else if (tail_type == TWO_SIDED) {
    pval = (z > 0) ? pnorm(z,0,sigma,false,log_p) : pnorm(z,0,sigma,true,log_p);
    if(log_p)
      pval += LOG_TWO;
    else
      pval *= 2.0;
  } else {
    throw("Invalid tail_type specified!");
  }
//printf("pval is %f\n",pval);

  return(pval);
}

//
// Computes wilcoxon rank sum statistic.  The answer is stored in *w and sumTies
// will contain nTies^3 - nTies summed over all tied groups.
//
void rankSum(double *x, int nX, double *y, int nY, double *w, int *sumTies) {
  double *f[2];
  double *r[2];
  int n[2];
  int i;

  f[0] = x;
  f[1] = y;
  r[0] = (double *) malloc(nX * sizeof(double));
  r[1] = (double *) malloc(nY * sizeof(double));
  n[0] = nX;
  n[1] = nY;
  
  *sumTies = multi_rank(2,f,n,r);
  for(*w=0, i=0; i<nX; i++)
    *w += r[0][i];

  free(r[0]);
  free(r[1]);

  return;
}


//
// Computes the ranks of an array of n lists of doubles, when all are sorted together
//
int multi_rank(int m, double **f, int *n, double **r) {
  int N,i,j,l,ties;
  Rank *rs;

  /* How many data in total? */
  for(N=0,i=0; i<m; i++)
    N += n[i];

  /* Initialise the rank structure. */
  rs = (Rank *) malloc(N * sizeof(Rank));
  for(i=0,l=0; i<m; i++) {
    for(j=0; j<n[i]; j++,l++) {
      rs[l].val    = f[i][j];
      rs[l].index  = j;
      rs[l].vector = i;
    }
  }

  ties = computeRanks(rs,N);

  /* Copy ranks to output array. */
  for(i=0; i<N; i++)
    r[rs[i].vector][rs[i].index] = rs[i].rank;

  free(rs);
  return(ties);
}


int computeRanks(Rank *rs, int n) {
  double rankSum,thisRank,oldVal;
  int i,j,nRanks,last,sumTies;

  /* Sort the doubles. */
  qsort((void *)rs, n, sizeof(Rank), (int (*)(const void*, const void*))(rankCmp));

  /* Assign ranks. */
  for(i=0, rankSum=0, last=0, nRanks=0, sumTies=0, oldVal=0; i<n; i++) {
    if((i == 0) || (fabs(rs[i].val - oldVal) < ZERO_COMPARE_EPSILON)) {
      rankSum += i+1;
      nRanks++;
    } else {
      sumTies += (nRanks * nRanks * nRanks) - nRanks;
      thisRank = rankSum/nRanks;
      /* printf("thisRank is %f, nRanks is %d, sumTies is now %d\n",thisRank,nRanks,sumTies); */
      for(j=last; j<i; j++)
        rs[j].rank = thisRank;
      last = i;
      rankSum = i+1;
      nRanks = 1;
    }
    oldVal = rs[i].val;
  }
  sumTies += (nRanks * nRanks * nRanks) - nRanks;
  thisRank = rankSum/nRanks;
  while(last < n)
    rs[last++].rank = thisRank;

  return(sumTies);
}

double choose(unsigned int n, unsigned int k) {
  unsigned int half_n;

  half_n = n/2;
  if(k > half_n)
    k = n-k;

  double ans = (double) n / (double) k;
  while(--k > 0) {
    --n;
    ans *= (double) n / (double) k;
  }

  return(ans);

}


//
// Returns the median of all pairwise differences between a two arrays of numbers.
// Used in conjunction with the rank sum test to return the Hodges Lehmann estimator.
//
double affxstat::medianDifference(double *x, int nX, double *y, int nY) {
  double *diffs,*diffsP,med;
  int nDiffs,i,j;

  if((nX <= 0) || (nY <= 0)) {
    throw("medianDifference: neither dataset should be empty.\n");
  } else {
    nDiffs = nX * nY;
    diffs = (double *) malloc(nDiffs * sizeof(double));
    diffsP = diffs;
    for(diffsP=diffs, i=0; i<nX; i++)
      for(j=0; j<nY; j++)
        *diffsP++ = x[i] - y[j];
    qsort((void *)diffs, nDiffs, sizeof(double), (int (*)(const void*, const void*))(doubleCmp));
    if((nDiffs % 2) == 0) {
      med = 0.5 * (diffs[nDiffs/2] + diffs[nDiffs/2 - 1]);
    } else {
      med = diffs[nDiffs/2];
    }
    free(diffs);
    return(med);
  }
}

double affxstat::PearsonCorrelation(float *x1, float *x2, int nPoints, float(*transform)(float x))
{
 if (nPoints < 2){
  return(numeric_limits<double>::quiet_NaN());
 }
 // Compute the means of the data
 double mean1=0;
 double mean2=0;
 for(int i=0; i<nPoints; i++)
 {
	 mean1 = mean1 + (transform == NULL ? x1[i] : transform(x1[i]));
	 mean2 = mean2 + (transform == NULL ? x2[i] : transform(x2[i]));
 }
 mean1 = mean1/nPoints;
 mean2 = mean2/nPoints;
 
 // compute the cross product sums and the sum of squares
 double crossproduct=0;
 double ssq1=0;
 double ssq2=0;
 for(int i=0; i<nPoints; i++)
 {
  crossproduct = crossproduct+((transform == NULL ? x1[i] : transform(x1[i]))-mean1)*((transform == NULL ? x2[i] : transform(x2[i]))-mean2);
  ssq1 = ssq1 + ((transform == NULL ? x1[i] : transform(x1[i]))-mean1)*((transform == NULL ? x1[i] : transform(x1[i]))-mean1);
  ssq2 = ssq2 + ((transform == NULL ? x2[i] : transform(x2[i]))-mean2)*((transform == NULL ? x2[i] : transform(x2[i]))-mean2);
 }
 if ((ssq1 ==0) || (ssq2==0)){
  return(numeric_limits<double>::quiet_NaN());
 }

 double correlation = crossproduct/(sqrt(ssq1)*sqrt(ssq2));
 return( correlation );
}

 

/* The following functions are taken from the GTYPE code base to calculate the Hardy-Weinburg equilibrium value. */

const double PI = 3.1415926535897932384626433832795;
const double TINY = 1e-20;
const double Eps = numeric_limits<double>::epsilon();
const double RelativeBound = 5 * Eps;
const int MaxIterations = 5000;
const int MeaningLess = -9999;

//----------------------------------//
static double logGamma( double x)
//----------------------------------//
{
	static double a[8] = {5.7083835261e-03, -1.910444077728e-03,
		8.4171387781295e-04, -5.952379913043012e-04,
		7.93650793500350248e-04, -2.777777777777681622553e-03,
		8.333333333333333331554247e-02, 0.9189385332046727417803297};

	static double c1[9] = {4.945235359296727046734888e0,
		2.018112620856775083915565e2, 2.290838373831346393026739e3,
		1.131967205903380828685045e4, 2.855724635671635335736389e4,
		3.848496228443793359990269e4, 2.637748787624195437963534e4,
		7.225813979700288197698961e3, -5.772156649015328605195174e-1};
	
	static double d1[8] = {6.748212550303777196073036e1,
		1.113332393857199323513008e3, 7.738757056935398733233834e3,
		2.763987074403340708898585e4, 5.499310206226157329794414e4,
		6.161122180066002127833352e4, 3.635127591501940507276287e4,
		8.785536302431013170870835e3};
     
	static double c2[9] = {4.974607845568932035012064e0,
		5.424138599891070494101986e2, 1.550693864978364947665077e4,
		1.847932904445632425417223e5, 1.088204769468828767498470e6,
		3.338152967987029735917223e6, 5.106661678927352456275255e6,
		3.074109054850539556250927e6, 4.227843350984671393993777e-1};

	static double d2[8] = {1.830328399370592604055942e2,
		7.765049321445005871323047e3, 1.331903827966074194402448e5,
		1.136705821321969608938755e6, 5.267964117437946917577538e6,
		1.346701454311101692290052e7, 1.782736530353274213975932e7,
		9.533095591844353613395747e6};
		
	static double c4[9] = {-1.474502166059939948905062e4,
		-2.426813369486704502836312e6, -1.214755574045093227939592e8,
		-2.663432449630976949898078e9, -2.940378956634553899906876e10,
		-1.702665737765398868392998e11, -4.926125793377430887588120e11,
		-5.606251856223951465078242e11, 1.791759469228055000094023e0};

	static double d4[8] = {-2.690530175870899333379843e3,
		-6.393885654300092398984238e5, -4.135599930241388052042842e7,
		-1.120872109616147941376570e9, -1.488613728678813811542398e10,
		-1.016803586272438228077304e11, -3.417476345507377132798597e11,
		-4.463158187419713286462081e11};

   	double numerator = 0;
	double denominator = 1;
	int n = 8;

	if (x <= Eps) {
		return - log(x);
	}
	else if (x <= 0.5) {
		for (int i = 0; i < n; ++i) {
			numerator = x * numerator + c1[i];
			denominator = x * denominator + d1[i];
		}
		return - log(x) + x * (x * numerator / denominator + c1[n]);
	}
	else if (x <= 0.6796875) {
		double z = x - 1;
		for (int i = 0; i < n; ++i) {
			numerator = z * numerator + c2[i];
			denominator = z * denominator + d2[i];
		}
		return -log(x) + z * (z * numerator / denominator + c2[n]);
	}
	else if (x <= 1.5) {
		double z = x - 1;
		for (int i = 0; i < n; ++i) {
			numerator = z * numerator + c1[i];
			denominator = z * denominator + d1[i];
		}
		return z * (z * numerator / denominator + c1[n]);
	}
	else if (x <= 4) {
		double z = x - 2;
		for (int i = 0; i < n; ++i) {
			numerator = z * numerator + c2[i];
			denominator = z * denominator + d2[i];
		}
		return z * (z * numerator / denominator + c2[n]);
	}
	else if (x <= 12) {
		double z = x - 4;
		for (int i = 0; i < n; ++i) {
			numerator = z * numerator + c4[i];
			denominator = z * denominator + d4[i];
		}
		return z * numerator / denominator + c4[n];
	}
	else {
		double z = x * x;
		numerator = a[0];
		for (int i = 1; i <= 6; ++i) {
			numerator = numerator / z + a[i];
		}
		numerator /= x;
		return numerator + log(x) * (x - 0.5) - x + a[7];
	}	
}

//----------------------------------------------------//
static double incompleteGamma( double x,  double a) 
//----------------------------------------------------//
{
	if (x == 0) {
		return 0;
	}
	
	bool convergence = false;
	double logGammaA = (a < TINY) ? logGamma(a+Eps) : logGamma(a);

	if (x < a + 1) {
		double temp = 1.0 / a;
		double update = temp;
		double a0 = a;	
		for (int i = 0; i < MaxIterations; ++i) {
			if (fabs(update) < RelativeBound * fabs(temp)) {
				convergence = true;
				break;
			}
			a0 = a0 + 1;
			update = x * update / a0;
			temp += update;
		}
		if (!convergence) 
		{
			return -1;
		}
		return temp * exp(a * log(x)  - x - logGammaA);
	}
	else {
		double p0 = 1,	p1 = x,	q0 = 0, q1 = 1;
		double c = 1, fold = q0, fnew = q1;
		for (int i = 1; i < MaxIterations; ++ i) {
			if (fabs(fnew - fold) < RelativeBound * fabs(fnew)) {
				convergence = true;
				break;
			}
			fold = fnew;
			double iMinusA = i - a;
			p0 = (p0 * iMinusA + p1) * c;
			q0 = (q0 * iMinusA + q1) * c;
			double iTimesC = i * c;
			p1 = p1 * iTimesC + p0 * x;
			q1 = q1 * iTimesC + q0 * x;
			c = 1.0 / p1;
			fnew = q1 * c;
		}
		if (!convergence) {
			return -1;
		}
		return 1 - fnew * exp(a * log(x) - x - logGammaA);
	}
}

//----------------------------------------------//
static double chiSquareCDF( double x,  int df) 
//----------------------------------------------//
{
	if (x <= 0.0) 
	{
		return 0.0;
	}
	else 
	{
		double y = incompleteGamma(x / 2.0, df / 2.0);
		if (y > 1.0)
			y = 1.0;
		if (y < 0.0)
			y = 0.0;
		return y;
	}
}

double affxstat::CalcHWEqPValue(int nAA, int nAB, int nBB)
{
    if (nAA + nAB + nBB <= 0)
        return 1.0;

	double fPHW = 0.0f;
	double fa = 0.0f, fb = 0.0f, fab = 0.0f, fchi = 0.0f, faa = 0.0f, fbb = 0.0f;
	double fea = 0.0f, feb = 0.0f, feab = 0.0f;
	double fTotal = nAA + nAB + nBB;

	fa = nAA/fTotal;
	fb = nBB/fTotal;
	fab = nAB/fTotal;

	faa = fa + 0.5*fab;
	fbb = fb + 0.5*fab;

	if(faa <= 0.0 || fbb <= 0.0)
		return 1.0;

	fea = fTotal*faa*faa;
	feab = 2*fTotal*faa*fbb;
	feb = fTotal*fbb*fbb;

	fchi = (fea - nAA)*(fea - nAA)/fea + 
			(feab - nAB)*(feab - nAB)/feab + 
			(feb - nBB)*(feb - nBB)/feb;

	fPHW = 1- chiSquareCDF(fchi, 1);
	fPHW = max(fPHW, 0.0);
	fPHW = min(fPHW, 1.0);

	return fPHW;
}
