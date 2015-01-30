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
 * @file   snp.label.cpp
 * @author Earl Hubbell
 * @date   Jan 22 2007
 * 
 * Does the labeling algorithmics.
 * 
 */

//
#include "label/snp.label.h"
//
#include "chipstream/BioTypes.h"
#include "chipstream/MatrixUtil.h"
//
#include "newmat.h"
//
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <memory>
#include <new>
#include <utility>
#include <vector>
//

// draft of plain-vanilla code for doing labeling numerics
// just match the code blindly

using namespace std;

/**
 * hash indexer
 * 
 * @param i - first index
 * @param j - second index
 * @param length - size of side of array
 */
inline int hax(int i, int j, int length)
{
    return(i*length+j);
}

/**
 * copies one float vector to another
 * 
 * @param final - destination vector
 * @param source - source vector
 * @param length - length of data
 */
void transfer(double *final, double *source, int length)
{
    for (int i=0; i<length; i++)
        final[i] = source[i];
}

/**
 * cumulative sum of source vector
 *
 * @param final - destination for cumulative sum
 * @param source - source of individual values
 * @param length - length of data
 */
void cumsum(double *final, double *source, int length)
{
    final[0] = source[0];
    for (int i=1; i<length; i++)
        final[i] = source[i]+final[i-1];
}

/**
 * cumulative sum of source vector using STL
 * 
 * @param final - destination for cumulative sum
 * @param source - source of individual values
 * @param length - length of data to sum over [carried over]
 */
void vcumsum(vector<double> &final,const vector<double> &source, int length)
{
    final.resize(length);
    final[0] = source[0];
    for (int i=1; i<length; i++)
        final[i] = source[i] + final[i-1];
}

/**
 * difference in cumulative sums: i.e. sum from b to a
 *
 * @param source - vector of cumulative sums
 * @param length - length of data
 * @param a - larger cumulative point
 * @param b - smaller cumulative point
 */
double delta(double *source,int length, int a, int b)
{
    double tmp;
    // difference in cumuldtive sums
    tmp = source[a]-source[b];
    if (a<b)
        tmp = 0;
    return(tmp);
}

/**
 * sum up the source array along one axis
 *
 * @param final - output, marginal sums of array
 * @param source - input, all the individual values
 * @param length - dimension of array
 * @param flag - which margin to sum to
 */
void apply_sum(vector<double> &final,const vector<double> &source, int length, int flag)
{
    int i,j;
    double sum;

    for (i=0; i<length; i++)
    {
        sum=0;
        for (j=0; j<length; j++)
        {
            if (flag)
                sum += source[hax(i,j,length)];
            else
                sum += source[hax(j,i,length)];
        }
        final[i] = sum;
    }
}

/**
 * takes a cumulative sum and makes it cumulate from the other end
 *
 * @param final - output, reverse cumulative sum
 * @param source - input, forward direction cumulative sum
 * @param length - length of data
 */
void reversecumsum(double *final,const double *source, int length)
{
    for (int i=0; i<length; i++)
        final[i] = source[length-1]-source[i];
}

/**
 * takes a cumulative sum and makes it cumulate from the other end
 *
 * @param final - output, reverse cumulative sum
 * @param source - input, forward direction cumulative sum
 * @param length - length of data (not needed for STL)
 */
void vreversecumsum(vector<double> &final, const vector<double> &source, int length)
{
    double tmp;
    final.resize(length);
    tmp = source[length-1];
    for (int i=0; i<length; i++)
        final[i] = tmp-source[i];
}

/**
 * takes outer difference of two vectors
 * 
 * @param final - output matrix
 * @param a - input vector
 * @param b - input vector
 * @param length - length of data
 */
void outer_minus(double *final, const double *a,const double *b, int length)
{
    int i,j;

    for (i=0; i<length; i++)
        for (j=i; j<length; j++)
            final[hax(i,j,length)] = a[j]-b[i];
}

/**
 * mask out a portion of an array
 *
 * @param final - array to be masked
 * @param length - size of margin of array
 * @param val - value to be placed in masked area
 */
void mask(double *final, int length, double val)
{
    int i,j;
    for (i=0; i<length; i++)
        for (j=0; j<i; j++)
            final[hax(i,j,length)] = val;
}

/**
 * mask out a portion of an array
 *
 * @param final - array to be masked
 * @param length - size of margin of array
 * @param val - value to be placed in masked area
 */
void vmask(vector<double> &final, int length, double val)
{
    int i,j;
    for (i=0; i<length; i++)
        for (j=0; j<i; j++)
            final[hax(i,j,length)] = val;
    
}

/**
 * sum two vectors
 *
 * @param final - output
 * @param a - input vector
 * @param b - input vector
 * @param length - length of input vectors
 */
void vsum(double *final,const double *a,const double *b, int length)
{
    for (int i=0; i<length; i++)
        final[i] = a[i]+b[i];
}

/**
 * add a constant to each element of a vector
 * 
 * @param final - output vector
 * @param a - input vector
 * @param k - value to be added
 * @param length - length of data
 */
void vplusscaldr(double *final, double *a, double k, int length)
{
    for (int i=0; i<length; i++)
        final[i] = a[i] + k;
}

/**
 * multiply each element of a vector by a value
 * 
 * @param X - vector
 * @param factor - factor to multiply each element by
 */
void vmult(vector<double> &X, double factor)
{
    for (int i=0; i<X.size(); i++)
        X[i] *= factor;
}

/**
 * inverts 3x3 matrix
 * specialty code, avoid calling external routines
 *
 * @param w - output inverse of matrix
 * @param v - input matrix
 */
void inversethree(vector<double> &w, vector<double> &v)
{
    w.resize(9);
    double det;
    // v = 9 element matrix 11,12,...33
    w[0] = v[4]*v[8]-v[5]*v[7];
    w[1] = v[2]*v[7]-v[1]*v[8];
    w[2] = v[1]*v[5]-v[2]*v[4];
    w[3] = v[5]*v[6]-v[3]*v[8];
    w[4] = v[0]*v[8]-v[2]*v[6];
    w[5] = v[2]*v[3]-v[0]*v[5];
    w[6] = v[3]*v[7]-v[4]*v[6];
    w[7] = v[1]*v[6]-v[0]*v[7];
    w[8] = v[0]*v[4]-v[1]*v[3];
    
    det = v[0]*w[0]+v[3]*w[1]+v[6]*w[2];
    for (int i=0; i<9; i++)
        w[i]/=det;
}

/**
 * sum two vectors
 *
 * @param z - output
 * @param a - input
 * @param b - input
 */
void sumthree(vector<double> &z, vector<double> &a, vector<double> &b)
{
    z.resize(a.size());
    for (int i=0; i<a.size(); i++)
        z[i] = a[i]+b[i];
}

/**
 * multiply two 3x3 matrices
 *
 * @param z - output
 * @param a - input 3x3 matrix (as 9 length vector)
 * @param b - input 3x3 matrix (as 9 length vector)
 */
void timesthree(vector<double> &z, vector<double> &a, vector<double> &b)
{
    z.resize(9);
    z[0] = a[0]*b[0]+a[1]*b[3]+a[2]*b[6];
    z[1] = a[0]*b[1]+a[1]*b[4]+a[2]*b[7];
    z[2] = a[0]*b[2]+a[1]*b[5]+a[2]*b[8];
    z[3] = a[3]*b[0]+a[4]*b[3]+a[5]*b[6];
    z[4] = a[3]*b[1]+a[4]*b[4]+a[5]*b[7];
    z[5] = a[3]*b[2]+a[4]*b[5]+a[5]*b[8];
    z[6] = a[6]*b[0]+a[7]*b[3]+a[8]*b[6];
    z[7] = a[6]*b[1]+a[7]*b[4]+a[8]*b[7];
    z[8] = a[6]*b[2]+a[7]*b[5]+a[8]*b[8];
}

/**
 * multiply a 1x3 vector by a 3x3 matrix
 *
 * @param z - output
 * @param a - 3x3 matrix (as 9 length vector)
 * @param b - 1x3 vector (as 3 length vector)
 */
void timesone(vector<double> &z, vector<double> &a, vector<double> &b)
{
    z.resize(3);
    z[0] = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    z[1] = a[3]*b[0]+a[4]*b[1]+a[5]*b[2];
    z[2] = a[6]*b[0]+a[7]*b[1]+a[8]*b[2];
}


/**
 * bayesian mean - prior + observed
 * 
 * @param k - prior pseudo-observations
 * @param m - prior mean
 * @param mm - observed sum of data
 * @param n - number of observations
 */
double b_mean(double k, double m,double mm, double n)
{
    // mean is prior + observed
    return (k*m+mm)/(k+n);
}

/**
 * bayesian variance - prior + observed + shift in means
 *
 * @param v - prior pseudo-observations
 * @param pv - prior variance
 * @param ddx - sum of squares of observations
 * @param cc - sum of observations
 * @param nc - number of observations
 * @param k - prior pseudo-observations for mean
 * @param mcc - observed mean
 * @param pm - prior location for mean
 */
double b_var(double v, double pv, double ddx, double cc, double nc, double k, double mcc, double pm)
{
    double tmp;
    // variance is prior
    tmp = v*pv;
    // plus observed ssq sum(x^2)-sum(x)^2/n
    tmp += ddx-cc*cc/(nc+.0001);
    // plus change in means squared
    tmp += (k/(k+nc))*nc*(mcc-pm)*(mcc-pm);
    // divided by effective number of observations
    tmp = tmp/(v+nc);
    return (tmp);
}

/**
 * gaussian likelihood, in convenient running sums form
 *
 * @param xcx - sum of squares of observations
 * @param lc  - posterior mean
 * @param xc - sum of observations
 * @param yc - number of observations
 * @param dvc - posterior variance
 */
double l_posterior(double xcx, double lc, double xc, double yc, double dvc)
{
    double tmp;
    tmp = (xcx-2*lc*xc+lc*lc*yc)/dvc + yc*log(dvc);
    return(tmp);
}

/**
 * gaussian likelihood, for just one piece of data
 * 
 * @param m - mean
 * @param x - observed value
 * @param v - variance
 */
double l_normal(double m, double x, double v)
{
    return ((m-x)*(m-x)/v + log(v));
}

double l_two_normal(double x, double y, cluster_data &cl, double inflatePRA)
{
    double dx,dy,cs,r,z;
    double TWOPI = 6.28318531;
    double up;

    // increase variance for calls
    up = 1 + inflatePRA/cl.k;
    
    dx = x-cl.m;
    dy = y-cl.ym;

    // compute correlation
    cs = sqrt(cl.ss*cl.yss);
    r = cl.xyss/cs;
    
    // log-likelihood
    z = dx*dx/(up*cl.ss) - 2*r*dx*dy/(up*cs) + dy*dy/(up*cl.yss);
    z /= 2*(1-r*r);
    z += log(TWOPI*up*cs*sqrt(1-r*r));    
    return(z);
}

/**
 * cheap inverse-gamma likelihood, for one piece of data
 * 
 * @param av - prior expected variance
 * @param dv - observed variance
 * @param v - pseudo-observations for prior
 */
double l_inverse(double av, double dv, double v)
{
    return(av/dv + (v+1)*log(dv));
}

/**
 * check if centers too close, return bignum if true
 *
 * @param a - cluster aa
 * @param b - cluster bb
 * @param h - cluster ab
 * @param shell - minimum distance squared for penalty
 */
double fix_shell(double a, double b, double h, double shell)
{
    // hard shell between entries
    if ((a-h)*(a-h)<shell)
        return(100);
    if ((h-b)*(h-b)<shell)
        return(100);
    if ((a-b)*(a-b)<6*shell)
        return(100);
    return(0);
}

/**
 * prevents hom clusters from approaching midline
 * not useful on B chip due to asymmetric probe choice
 *
 * @param a - cluster aa mean
 * @param b - cluster bb mean
 * @param h - cluater ab mean
 * @param shell - homs can't be closer to center than this or penalty
 */
double fix_hom_shell(double a, double b, double h, double shell)
{
    // hard shell between homs and midline
    double tmp;
    tmp = sqrt(shell);
    if (a< tmp)
        return(100);
    if (b> -1*tmp)
        return(100);
    return(0);
}

/**
 * Given 3 genotype probabilities, make a call/confidence output
 * 
 * @call - call of final genotype 2,1,0
 * @conf - confidence of final genotypes 1-likelihood
 * @txd  - relative probability of d = b allele dominant
 * @txh  - relative probability of a+b allele
 * @txc  - relative probability of c = a allele dominant
 */
void make_a_call(affx::GType &call, double &conf, double txd, double txh, double txc)
{
    double tmp;

  call = affx::BB;  // 2;
  tmp = txc;
  if (txh>tmp) {
    call = affx::AB; // 1
    tmp = txh;
  }
  if (txd>tmp) {
    call = affx::AA; // 0
    tmp =txd;    
  }
  conf = 1-tmp;
}
    
/**
 *  Given cluster parameters for each genotype, call data as genotypes
 *
 * @param call - output called genotypes
 * @param conf - output confidences for each call
 * @param x    - input data in contrast space
 * @param length - number of data points
 * @param rbma - aa cluster mean
 * @param rbva - aa cluster variance
 * @param rbmb - bb cluster mean
 * @param rbvb - bb cluster variance
 * @param rbmh - ab cluster mean
 * @param rbvh - ab cluster variance
 * @param cn - copy number - are we 1 copy (all homs) or two allowing hets
 */
void make_calls(vector<affx::GType> &call,
                vector<double> &conf,
                vector<double> &x,
                int length,
                cluster_data &pa,
                cluster_data &ph,
                cluster_data &pb,
                int cn,
                int freqflag,
                vector<double> &vTxa,
                vector<double> &vTxh,
                vector<double> &vTxb,
                bool store_probabilities
    )
{
    int i;
    double tmp,sm;
    double txa,txb,txh;
    double lfreqa, lfreqb,lfreqh;

    // precompute modification for frequency of call
    // i.e. rare alleles less likely than high frequency alleles
    // change this to use frequency instead of precision of mean
    lfreqa=0;
    lfreqb=0;
    lfreqh=0;
    if (freqflag)
    {
        // log frequency
        // will be normalized out when computing relative probability
        // so no need to adjust here
        // also note posterior is always >0 
        // note negation later
        lfreqa = -log(pa.k);
        lfreqb = -log(pb.k);
        lfreqh = -log(ph.k);
    }

    for (i=0; i<length; i++)
    {
        txa=l_normal(x[i],pa.m,pa.ss)/2 + lfreqa;
        txb=l_normal(x[i],pb.m,pb.ss)/2 + lfreqb;
        txh=l_normal(x[i],ph.m,ph.ss)/2 + lfreqh;
        tmp = txa;
        if (txh<tmp && cn>1)
            tmp=txh;
        if (txb<tmp)
            tmp =txb;
        txa = txa-tmp;
        txb = txb-tmp;
        txh = txh-tmp;
        txa = exp(-txa);
        txb = exp(-txb);
        txh = exp(-txh);
        if (cn<2)
            txh = 0; // no chance of hets if 1 copy
        sm = txa+txb+txh;
        txa /= sm;
        txb /= sm;
        txh /= sm;
        make_a_call(call[i],conf[i],txa,txh,txb);
        if (store_probabilities) {
            vTxa[i]=txa;
            vTxh[i]=txh;
            vTxb[i]=txb;            
        }
    }

  // JHG -debug
  // printVec("conf",conf);
  // printVec("x",x);
}

//make_two_calls(call,conf,x,y,x.size(),sp.posterior, sp.copynumber,sp.mix);

/**
 *  Given cluster parameters for each genotype, call data as genotypes
 *
 * @param call - output called genotypes
 * @param conf - output confidences for each call
 * @param x    - input data in contrast space
 * @param length - number of data points
 * @param rbma - aa cluster mean
 * @param rbva - aa cluster variance
 * @param rbmb - bb cluster mean
 * @param rbvb - bb cluster variance
 * @param rbmh - ab cluster mean
 * @param rbvh - ab cluster variance
 * @param cn - copy number - are we 1 copy (all homs) or two allowing hets
 */
void make_two_calls(vector<affx::GType> &call,
                    vector<double> &conf,
                    vector<double> &x,
                    vector<double> &y,
                    int length,
                    snp_distribution &p,
                    int cn,
                    int freqflag,
                    double ocean,
                    double inflatePRA,
                    vector<double> &vTxa,
                    vector<double> &vTxh,
                    vector<double> &vTxb,
                    bool store_probabilities
    )
{
    int i;
    double tmp,sm;
    double txa,txb,txh,txt;
    double lfreqa, lfreqb,lfreqh, ltot;


    // precompute modification for frequency of call
    // i.e. rare alleles less likely than outliers to high frequency alleles
    // Use frequency here, no longer reuse posterior precision of mean
    lfreqa=0;
    lfreqb=0;
    lfreqh=0;
    if (freqflag)
    {
        // log frequency
        // also note posterior is always >0 
        // f not working?
        // note negation later
        lfreqa = -log(p.aa.k);
        lfreqb = -log(p.bb.k);
        lfreqh = -log(p.ab.k);
        ltot = log(p.aa.k+p.bb.k+p.ab.k);
        lfreqa += ltot;
        lfreqb += ltot;
        lfreqh += ltot;
    }

    // want to maintain absolute scale in case comparing to 'absolute' issues
    for (i=0; i<length; i++)
    {
        txa=l_two_normal(x[i],y[i],p.aa,inflatePRA) + lfreqa;
        txb=l_two_normal(x[i],y[i],p.bb,inflatePRA) + lfreqb;
        txh=l_two_normal(x[i],y[i],p.ab,inflatePRA) + lfreqh;

        //cout << txa << "\t" << txh << "\t" << txb << endl;
        // find the best of three
        tmp = txa;
        if (txh<tmp && cn>1)
            tmp=txh;
        if (txb<tmp)
            tmp =txb;
        txa = txa-tmp;
        txb = txb-tmp;
        txh = txh-tmp;
        txa = exp(-txa);
        txb = exp(-txb);
        txh = exp(-txh);
        txt = ocean*exp(tmp); // how deep is the uniform ocean relative to these guys
        //cout << txa << "\t" << txh << "\t" << txb << endl;
        if (cn<2)
            txh = 0; // no chance of hets if 1 copy
        // rescale by ocean to reduce confidence
        sm = txa+txb+txh+txt;
        txa /= sm;
        txb /= sm;
        txh /= sm;
        
        //cout << txa << "\t" << txh << "\t" << txb << endl;
        make_a_call(call[i],conf[i],txa,txh,txb);
        if (store_probabilities) {
            vTxa[i]=txa;
            vTxh[i]=txh;
            vTxb[i]=txb;
        }
    }

  // JHG -debug
  // printVec("conf",conf);
  // printVec("x",x);
}

/** make utility for comparing pairs */
struct first_less : public binary_function<pair<int,double>,pair<int,double>,bool>
{
    bool operator() (const pair<int,double>& x, const pair<int,double>&y) const
    {
        return(x.first <y.first);
    }
};

/** utility for comparing pairs */
struct second_less : public binary_function<pair<int,double>,pair<int,double>,bool>
{
    bool operator()(const pair<int,double>& x,const pair<int,double>&y) const
    {
        return(x.second<y.second);
    }
};    

/**
 * make calls using the final labelings - primarily a matter of unsorting the contrast values
 *
 * @param call - output called genotypes
 * @param conf - confidence of call
 * @param x - raw contrast values (unsorted)
 * @param length - number of data points
 * @param ed - estimated relative probability for d = bb genotype (sorted)
 * @param ec - estimated relative probability for c = aa genotype (sorted)
 * @param eh - estimated relative probability for h = ab genotype (sorted)
 * @param cn - copynumber 1 = homs only, copynumber 2 = hets allowed
 */
void make_label_calls(vector<affx::GType> &call, vector<double> &conf, vector<double> &x, 
                      int length, vector<double> &ed, vector<double> &ec, vector<double> &eh, int cn)
{
    // z = sorted values
    // ed, ec, eh = cluster memberships
    // x = values to be matched with z
    // note: for future work, this should be "bin" compatible

    int i,j;
    vector< pair<int,double> > data;

    for (i=0; i<length; i++)
    {
        pair<int,double> p;
        p.first = i;
        p.second = x[i];
        data.push_back(p);
    }

    sort(data.begin(),data.end(),second_less());

    for (i=0; i<length; i++)
    {
        // for each value in data (sorted x)
        // find the corresponding z value
        // place the corresponding call/conf for x
        j = data[i].first;
        if (cn>1)
            make_a_call(call[j],conf[j],ec[i],eh[i],ed[i]);    
        else
            make_a_call(call[j],conf[j],ec[i],0,ed[i]);
    }
}

/**
 * assigns the sorted constrast data into "bins" for fast processing
 *
 * @param zall - the bins of the data, number per bin, sums, sum-squares, hints
 * @param z - the raw data, sorted by constrast
 * @param zerohints - what are the reference genotypes supplied
 * @param zlength - how many values in the data
 * @param bins - how many bins at most for the data
 * @param hok - can homs be flipped in hints?
 */
int setup_bins(binned_data &zall, vector<double> &z, vector<double> &w, vector<int> &zerohints, vector<double> &zeroInbred, int zlength, int bins, int hok, double Cpenalty)
{

    double delta,boundary;
    int ti,i;

    delta = (z[zlength-1]-z[0])/(bins+1);

    zall.zz[0] = 0;
    zall.zzx[0]= 0;
    zall.zy[0] = 0;
    zall.zzy[0] = 0;
    zall.zxy[0] = 0;
    zall.nz[0] = 0;
    zall.nod[0] = 0;
    zall.noh[0] = 0;    
    zall.noc[0] = 0;
    ti = 0;
    boundary = z[0]-1;
    for (i=0; i<zlength; i++)
    {
        // when do we go to the next thing?
        // when a) bins is turned off
        // b) when fewer points than bins
        // c) when point outside domain
        if (z[i]>boundary || bins==0 || zlength<bins)
        {
            ti++;
            boundary = z[i]+delta;
            // initialize the new bin
            zall.nz[ti] = 0;
            zall.zz[ti] = 0;
            zall.zzx[ti] = 0;
            zall.zy[ti] = 0;
            zall.zzy[ti] = 0;
            zall.zxy[ti] = 0;
            // hints
            zall.nod[ti] = 0;    
            zall.noh[ti] = 0;
            zall.noc[ti] = 0;
        }
        // add current point to bins
        zall.nz[ti]  += 1;
        zall.zz[ti]  += z[i];
        zall.zzx[ti] += z[i]*z[i];
        zall.zy[ti]  += w[i];
        zall.zzy[ti] += w[i]*w[i];
        zall.zxy[ti] += z[i]*w[i];
        // count contradictions
        // if no information, no contradictions
        if (zerohints[i]==0)
        {
            if (hok<1)
                zall.nod[ti] +=Cpenalty;
            zall.noh[ti] +=Cpenalty;
        }
        if (zerohints[i]==1)
        {
            zall.nod[ti] +=Cpenalty;
            zall.noc[ti] +=Cpenalty;
        }
        if (zerohints[i]==2)
        {
            zall.noh[ti] +=Cpenalty;
            if (hok<1)
                zall.noc[ti] +=Cpenalty;
        }
        // bin ti contains item [i] which has a bias for or against hets by being inbred.
        zall.noh[ti] += zeroInbred[i];
    }
    
     // ti indicates last point here
    // how much did we actually use
    return(ti+1);
}

/**
 * generate sums of data based on weighted cluster membership
 * 
 * @param n - effective number of data points 
 * @param m - effective sum of data points
 * @param s - effective sum of squares of data points
 * @param zall - binned contrast data for points
 * @param wi - weight of bins
 * @param ZS - number of bins
 */
void compute_binned_cluster(double &n, double &m, double &s, binned_data &zall,
         vector<double> &wi, int ZS)
{
    n=0;
    m=0;
    s=0;

    for (int i=0; i<ZS-1; i++)
    {
        n +=  zall.nz[i+1]*wi[i];
        m +=  zall.zz[i+1]*wi[i];
        s +=  zall.zzx[i+1]*wi[i];
    }
}

/**
 * generate sums of data based on weighted cluster membership
 * 
 * @param n - effective number of data points 
 * @param m - effective sum of data points
 * @param s - effective sum of squares of data points
 * @param zall - binned contrast data for points
 * @param wi - weight of bins
 * @param ZS - number of bins
 */
void compute_two_binned_cluster(vector<double> &six, binned_data &zall,
         vector<double> &wi, int ZS)
{
    for (int j=0; j<6; j++)
        six[j] = 0;
    /*n=0;
    m=0;
    s=0;
    my=0;
    sy=0;
    xy=0;*/
    for (int i=0; i<ZS-1; i++)
    {
        six[0] +=  zall.nz[i+1]*wi[i];
        six[1] +=  zall.zz[i+1]*wi[i];
        six[2] +=  zall.zzx[i+1]*wi[i];
        six[3] +=  zall.zy[i+1]*wi[i];
        six[4] +=  zall.zzy[i+1]*wi[i];
        six[5] +=  zall.zxy[i+1]*wi[i];
    }
}

/**
 * take probabilities of bins and turn them into probabilities of individual points
 *
 * @param tex - output per point
 * @param length - number of real data points
 * @param ex - binned data
 * @param nz - number of points per bin
 * @param ZS - number of bins
 */
void expand_data(vector<double> &tex, int length,vector<double> &ex,vector<double> &nz,int ZS){
    int ti,tj,i;
    i=0;
    // for each bin
    for (ti=0; ti<ZS-1; ti++)
    {
        // for each data point in bin
        for (tj=0; tj<floor(nz[ti+1]+.0001); tj++)
        {
            // expanded value is value for bin
            tex[i] = ex[ti];
            i++;
        }
    }
}

/**
 *  convert log-likelihood of all labelings into relative probability
 * 
 * @param q - matrix to be converted
 * @param ZS - edge size of matrix
 */
void qtoeq(vector<double> &q, int ZS)
{
    double qmin;
    int tin;
    int i,j;

// compute minimum of useful values
    // so as to rescale log-likelihood to relative levels
    qmin = q[0];
    for (i=0; i<ZS; i++)
        for (j=i; j<ZS; j++)
        {
            tin = hax(i,j,ZS);
            if (qmin>q[tin])
                qmin=q[tin];
        }
    // q, eq
    for (i=0; i<q.size(); i++)
    {
        q[i] = qmin-q[i];
        q[i] = exp(q[i]);
    }
    // mask out useless values
    vmask(q,ZS,0);
}

/**
 * turn the overall relative probability of labelings into relative probability of individual genotypes
 *
 * @param ec - relative probability of c = a allele dominant
 * @param eh - relative probability of h = ab allele about equal
 * @param ed - relative probability of d = b allele dominant
 * @param ZS - marginal size of matrix
 * @param q  - matrix of relative probability of labelings
 */
void eqtoproblabel(vector<double> &ec, vector<double> &eh, vector<double> &ed, int ZS, vector<double> q)
{
    int i;
    double total;

    total =0; 
    for (i=0; i<q.size(); i++)
        total += q[i];
    
    // get cluster memberships
    //  note 0 here is entry "0" in the original list of length l
    // entry j is in ec[l] iff j<=l
    // so straight cumulative sum
    vector<double> eqt;
    eqt.resize(ZS);
    apply_sum(eqt,q,ZS,0);
    vcumsum(ec,eqt,ZS);
    // entry i is in ed[l] iff i>l
    // do cumsum to get i<=l, then subtract
    apply_sum(eqt,q,ZS,1); 
    vcumsum(ed,eqt,ZS);
    
    for (i=0; i<ZS-1; i++)
    {
        ed[i] = total-ed[i];
        // entry i is in eh[l] iff not in a or b
        eh[i] = total-ed[i]-ec[i];
        // and therefore the sum is always total
        //sm     = total;
        ed[i] /= total;
        ec[i] /= total;
        eh[i] /= total;
    }
}

/**
 *  setup the matrixes for doing individual BRLMM steps
 *
 * @param m - prior mean vector
 * @param Minv - inverse matrix of covariances (pseudo-observations) for means
 * @param Sinv - inverse matrix of variances of observations within genotypes
 * @param snp_param - parameters of algorithm
 */
void setupBRLMM(vector<double> &m, vector<double> &Minv, vector<double> &Sinv, snp_param &sp)
{
            m.assign(3,0);
            m[0] = sp.prior.aa.m;
            m[1] = sp.prior.ab.m;
            m[2] = sp.prior.bb.m;
            // only work with Minv
            //inversethree(Minv,M);
            Minv.assign(9,0);
            // set up allowed variation in cluster center
            Minv[0] = sp.prior.aa.k/sp.prior.aa.ss;
            Minv[4] = sp.prior.ab.k/sp.prior.ab.ss;
            Minv[8] = sp.prior.bb.k/sp.prior.bb.ss;
            // handle cross-terms here
            Minv[1] = sp.prior.xah/sqrt(sp.prior.aa.ss*sp.prior.ab.ss);
            Minv[2] = sp.prior.xab/sqrt(sp.prior.aa.ss*sp.prior.bb.ss);
            Minv[3] = Minv[1];
            Minv[5] = sp.prior.xhb/sqrt(sp.prior.ab.ss*sp.prior.bb.ss);
            Minv[6] = Minv[2];
            Minv[7] = Minv[5];
            // only work with Sinv;
            //inversethree(Sinv,S);
            Sinv.assign(9,0);
            Sinv[0] = 1/sp.prior.aa.ss;
            Sinv[4] = 1/sp.prior.ab.ss;
            Sinv[8] = 1/sp.prior.bb.ss;
}


void force_isotonic(double &mb, double &mh, double &ma, double wb, double wh, double wa,double delta)
{
    double gamma, tmp;
    // forces means to be isotonic, separated by delta
    // moves low weight centers most
    gamma = delta*(wb-wa)/(wb+wh+wa);
    // shift means towards center of mass by appropriate values
    mb += delta-gamma;
    mh += 0-gamma;
    ma += -1*delta-gamma;
    // Pool Adjacent Violators
    // possible that mb<=mh<=ma not true
    if (mb>mh)
    {
        tmp = (wb*mb+wh*mh)/(wb+wh);
        mb = tmp;
        mh = tmp;
    }
    if (mh>ma)
    {
        tmp = (wh*mh+wa*ma)/(wh+wa);
        mh = tmp;
        ma = tmp;
        if (mb>mh)
        {
            tmp = (wb*mb+wh*mh+wa*ma)/(wb+wh+wa);
            mb = tmp;
            mh = tmp;
            ma = tmp;
        }
    }
    // now mb<=mh<=ma is true and average data hasn't changed
    // now we make sure that they are separated by delta
    mb -= delta-gamma;
    mh -= 0-gamma;
    ma -= -1*delta-gamma;
    // now mb+delta<=mh, mh+delta<ma
}

void force_IsoHetY(snp_distribution &post, double IsoHetWeight)
{
    // ensure that the Y means of clusters are in the right order
    // Het mean > Hom mean average
    double ystarab, ystar;
    double wab,wh,wa,wb;
    double delta,deltab,deltaa;
    double dxah,dxhb,dxab;
    

    // interpolated Y value from hom clusters
    dxhb = post.ab.m-post.bb.m;
    dxah = post.aa.m-post.ab.m;
    dxab = post.aa.m-post.bb.m;
    
    ystarab = (dxhb*post.aa.ym + dxah*post.bb.ym)/dxab;
    if (ystarab>post.ab.ym)
    {
        wh = post.ab.k*IsoHetWeight;
        wa = post.aa.k;
        wb = post.bb.k;
        wab = wa+wb;
        // violation of constraint
        ystar = (ystarab*wab+post.ab.ym*wh)/(wab+wh);
        post.ab.ym = ystar; // here's where I've moved the het

        // now move homs to make their center work
        // on average must move this much
        delta = (ystar-ystarab);
        // move each by weight of corresponding entity
        deltab = delta/((wb*dxhb+wa*dxah)/(wa*dxab));
        deltaa = delta/((wb*dxhb+wa*dxah)/(wb*dxab));
        post.bb.ym += deltab;
        post.aa.ym += deltaa;
    }
}

void bayes_two_variance(cluster_data &post, vector<double> &sixc, cluster_data &prior)
{
    post.v = prior.v + sixc[0];

    post.ss = prior.v*prior.ss;  // prior
        post.ss += sixc[2]-sixc[1]*sixc[1]/(sixc[0]+0.001); // observed ssq, avoid dividebyzero
        post.ss += ((prior.k*sixc[0])/(prior.k+sixc[0])) * (post.m-prior.m)*(post.m-prior.m);  
        post.ss /= post.v;

    post.yss = prior.v*prior.yss; // prior
    post.yss += sixc[4]-sixc[3]*sixc[3]/(sixc[0]+0.001);
    post.yss += ((prior.k*sixc[0])/(prior.k+sixc[0])) * (post.ym-prior.ym)*(post.ym-prior.ym);  
    post.yss /= post.v;

    //cout << prior.v*prior.yss << "\t" << sixc[4]-sixc[3]*sixc[3]/(sixc[0]+0.001) << "\t" << ((prior.k*sixc[0])/(prior.k+sixc[0])) * (post.ym-prior.ym)*(post.ym-prior.ym) << endl;

    post.xyss = prior.v*prior.xyss; // prior
    post.xyss += sixc[5]-sixc[1]*sixc[3]/(sixc[0]+0.001);
    post.xyss += ((prior.k*sixc[0])/(prior.k+sixc[0])) * (post.ym-prior.ym)*(post.m-prior.m);
    post.xyss /= post.v;
}



void shrink_cluster_variance(cluster_data &out, cluster_data &a, cluster_data &b, cluster_data &c, double lambda)
{
    double main,cross;
    
    main = 3-2*lambda;
    cross = lambda;

    out.ss = (main*a.ss*a.v+cross*b.ss*b.v+cross*c.ss*c.v)/(main*a.v+cross*b.v+cross*c.v);
    out.yss = (main*a.yss*a.v+cross*b.yss*b.v+cross*c.yss*c.v)/(main*a.v+cross*b.v+cross*c.v);
    out.xyss = sqrt(out.ss*out.yss)*a.xyss/(sqrt(a.ss*a.yss));  // covariance shrinkage, same correlation
    // what should degree of confidence in variance be???
}

void bayes_two_shrink_variance(snp_param &sp)
{
    // keep the correlation the same within each cluster
    // move the variances towards each other?
    // analogue of 1-d common variance
    
    snp_distribution tmp;
    tmp.Copy(sp.posterior);

    shrink_cluster_variance(tmp.aa,sp.posterior.aa,sp.posterior.ab,sp.posterior.bb,sp.lambda);
    shrink_cluster_variance(tmp.ab,sp.posterior.ab,sp.posterior.bb,sp.posterior.aa,sp.lambda);
    shrink_cluster_variance(tmp.bb,sp.posterior.bb,sp.posterior.aa,sp.posterior.ab,sp.lambda);
    
    sp.posterior.Copy(tmp); // copy shrunk values
}

void Dump_Matrix(Matrix &M)
{
    cout << endl;
    
    for (int i=0; i<M.Nrows(); i++)
    {
        for (int j=0; j<M.Ncols(); j++)
        {
            cout << M.element(i,j) << "\t";
        }
        cout << endl;
    }
}

void Dump_Column(ColumnVector &V)
{
    cout << endl;
    for (int i=0; i<V.Nrows(); i++)
            cout << V.element(i) << "\t";
    cout << endl;
}

void labels_two_posterior(vector<double> &ec, vector<double> &eh, vector<double> &ed, int ZS, binned_data &zall, snp_param &sp)
{
    // Means = (K+N)^-1 (K*u+N*X)
    // Fix isotonic X dimension means
    // K<- K+N
    // V1<-V + SSQ + [kn/(k+n)] * (u-Means)^2
    // Vshrink = common var    
    
    vector<double> sixd,sixc,sixh;
    sixd.resize(6);
    sixc.resize(6);
    sixh.resize(6);
    // compute basic stats, raw mean/sum-square/number of entries
    // for each genotype assignment mixed fraction
    compute_two_binned_cluster(sixd,zall,ed,ZS);
    compute_two_binned_cluster(sixc,zall,ec,ZS);
    compute_two_binned_cluster(sixh,zall,eh,ZS);

    // now we set up the N matrix
    // why is this banded in BRLMM?
    Matrix N(6,6);
      N = 0;
      N.element(0,0) = sixc[0];
      N.element(1,1) = sixc[0];
     N.element(2,2) = sixh[0];
      N.element(3,3) = sixh[0];
      N.element(4,4) = sixd[0];
      N.element(5,5) = sixd[0];

    Matrix K(6,6);
    K=0;
    K.element(0,0) = sp.prior.aa.k;
    K.element(1,1) = sp.prior.aa.k;
    K.element(2,2) = sp.prior.ab.k;
    K.element(3,3) = sp.prior.ab.k;
    K.element(4,4) = sp.prior.bb.k;
    K.element(5,5) = sp.prior.bb.k;
    // fill in off-diagonal elements for x to x
    K.element(0,2) = sp.prior.xah;
    K.element(2,0) = sp.prior.xah;
    K.element(2,4) = sp.prior.xhb;
    K.element(4,2) = sp.prior.xhb;
    K.element(0,4) = sp.prior.xab;
    K.element(4,0) = sp.prior.xab;
    // fill in off diagonal elements for y to y
    K.element(1,3) = sp.prior.yah;
    K.element(3,1) = sp.prior.yah;
    K.element(3,5) = sp.prior.yhb;
    K.element(5,3) = sp.prior.yhb;
    K.element(1,5) = sp.prior.yab;
    K.element(5,1) = sp.prior.yab;

    // temporary!!!
    // Need 4x3 covariance pairs

    // left hand side
    Matrix Left = K + N;
    //Dump_Matrix(Left);
    //Dump_Matrix(K);
    //Dump_Matrix(N);
    Left = Left.i();
    //Dump_Matrix(Left);    
    
    
    ColumnVector Mu;
    Mu.ReSize(6);
    Mu.element(0) = sp.prior.aa.m;
    Mu.element(1) = sp.prior.aa.ym;
    Mu.element(2) = sp.prior.ab.m;
    Mu.element(3) = sp.prior.ab.ym;
    Mu.element(4) = sp.prior.bb.m;
    Mu.element(5) = sp.prior.bb.ym;
    //Dump_Column(Mu);


    Matrix KMu = K*Mu;
    //Dump_Matrix(KMu);
    
    // premultipy, prevent div/0 in computing mean
    
    ColumnVector Nv;
    Nv.ReSize(6);
    Nv.element(0) = sixc[1];
    Nv.element(1) = sixc[3];
    Nv.element(2) = sixh[1];
    Nv.element(3) = sixh[3];
    Nv.element(4) = sixd[1];
    Nv.element(5) = sixd[3];

    //Dump_Column(Nv);
    
    Matrix Right = KMu + Nv;

    //Dump_Matrix(Right);

    Matrix Final = Left*Right;

    //Dump_Matrix(Final);
    // update means
    sp.posterior.aa.m  = Final.element(0,0);
    sp.posterior.aa.ym = Final.element(1,0);
    sp.posterior.ab.m  = Final.element(2,0);
    sp.posterior.ab.ym = Final.element(3,0);
    sp.posterior.bb.m  = Final.element(4,0);
    sp.posterior.bb.ym = Final.element(5,0); 

    //cout << "Executed through" << endl;

    // update posterior strength
    sp.posterior.aa.k = sixc[0] + sp.prior.aa.k;
    sp.posterior.ab.k = sixh[0] + sp.prior.ab.k;
    sp.posterior.bb.k = sixd[0] + sp.prior.bb.k;

    // update posterior frequency
    // Omit currently
    
    // fix up positions of posterior in contrast dimension
    if (sp.hardshell==3)
        force_isotonic(sp.posterior.bb.m,sp.posterior.ab.m,sp.posterior.aa.m,sp.posterior.bb.k,sp.posterior.ab.k,sp.posterior.aa.k, sp.shellbarrier);

    // put het above hom average
    // this is to avoid unrealistic situations in the face of copy number
    // the numerical control is the relative weight of Hets: near zero means move het, near infinity means move hom
    // 1.0 means weight hets and homs equally for the final line.
    if (sp.IsoHetY>0)
        force_IsoHetY(sp.posterior,sp.IsoHetY);

    // now I have the means computed!!!
    // on to compute variances
    // Technically, these are variances computed from the non-isotonic means
    bayes_two_variance(sp.posterior.aa,sixc,sp.prior.aa);
    bayes_two_variance(sp.posterior.ab,sixh,sp.prior.ab);
    bayes_two_variance(sp.posterior.bb,sixd,sp.prior.bb);
    // how do I shrink variances in 2-D correctly?
    bayes_two_shrink_variance(sp);
    
    // cluster covariances in pseudo-obs
    sp.posterior.xah = sp.prior.xah;
    sp.posterior.xhb = sp.prior.xhb;
    sp.posterior.xab = sp.prior.xab;
    sp.posterior.yah = sp.prior.yah;
    sp.posterior.yhb = sp.prior.yhb;
    sp.posterior.yab = sp.prior.yab;
}


/**
 * convert binned label relative probabilities into cluster parameters
 *
 * @param ec - relative probability of c = a allele dominates
 * @param eh - relative probability of h = ab about equal
 * @param ed - relative probability of d = b allele dominates
 * @param ZS - length of bins
 * @param zall - all binned data, sums, sums of squares
 * @param sp - snp parameters
 */
void labelstoposterior(vector<double> &ec, vector<double> &eh, vector<double> &ed, int ZS, binned_data &zall, snp_param &sp)
{
    // compute posterior mean/variance based on model average labels
    double rnd,rnc,rnh,rmd,rmc,rmh,rsd,rsc,rsh;
    double rbma,rbmb,rbmh,rbva,rbvb,rbvh;
    //double tmp;
    
    // compute basic stats, raw mean/sum-square/number of entries
    // for each genotype assignment mixed fraction
    compute_binned_cluster(rnd,rmd,rsd,zall,ed,ZS);
    compute_binned_cluster(rnc,rmc,rsc,zall,ec,ZS);
    compute_binned_cluster(rnh,rmh,rsh,zall,eh,ZS);

    // Compute the posterior mean/variance
    // i.e. do BRLMM step here
    //rbmb = b_mean(sp.bb.k,sp.bb.m,rmd,rnd);
    //rbma = b_mean(sp.aa.k,sp.aa.m,rmc,rnc);
    //rbmh = b_mean(sp.ab.k,sp.ab.m,rmh,rnh);
    
// BRLMM
// setup BRLMM
            vector<double> m;
            vector<double> Minv;
            vector<double> Sinv;
            setupBRLMM(m,Minv,Sinv,sp);
            // only work with Sinv;
            // scratch vectors
            vector<double> Tz;
            vector<double> Zt;
            vector<double> Xinv;
            vector<double> Xt;
            vector<double> Mz;
            // vectors that change
            // sum of x
            vector<double> Nv;
            // number of data points
            vector<double> N;
    // brlmm here
            Nv.resize(3);
            Nv[0] = rmc;
            Nv[1] = rmh;
            Nv[2] = rmd;
            N.assign(9,0);
            N[0] = rnc;
            N[4] = rnh;
            N[8] = rnd;
            timesthree(Tz,N,Sinv);
            sumthree(Zt, Minv,Tz);
            inversethree(Xinv,Zt);
            timesone(Tz,Sinv, Nv);
            timesone(Zt,Minv, m);
            sumthree(Xt,Tz,Zt);
            timesone(Mz,Xinv,Xt);
            rbma = Mz[0];
            rbmh = Mz[1];
            rbmb = Mz[2];
    // fix up positions of posterior
    if (sp.hardshell==3)
        force_isotonic(rbmb,rbmh,rbma,rnd+sp.prior.bb.k,rnh+sp.prior.ab.k,rnd+sp.prior.aa.k, sp.shellbarrier);

    // do within-cluster variance as before
    rbvb = b_var(sp.prior.bb.v,sp.prior.bb.ss,rsd,rmd,rnd,sp.prior.bb.k,rbmb,sp.prior.bb.m);
    rbva = b_var(sp.prior.aa.v,sp.prior.aa.ss,rsc,rmc,rnc,sp.prior.aa.k,rbma,sp.prior.aa.m);
    rbvh = b_var(sp.prior.ab.v,sp.prior.ab.ss,rsh,rmh,rnh,sp.prior.ab.k,rbmh,sp.prior.ab.m);
    
    // fix up variances if needed
    if (sp.comvar==1 )
    {
            double tc,td,th;
            double tnc,tnd,tnh;
            double tl,tlm;
            tnd = rnd+sp.prior.bb.v;
            tnc = rnc+sp.prior.aa.v;
            tnh = rnh+sp.prior.ab.v;
            td = rbvb*tnd;
            tc = rbva*tnc;
            th = rbvh*tnh;
            // cross-talk between variances 0 = independent, 1 = common variance
            tl = sp.lambda;
            tlm = 3-2*sp.lambda;
            rbvb = (tlm*td+tl*tc+tl*th)/(tlm*tnd+tl*tnc+tl*tnh);
            rbva = (tl*td+tlm*tc+tl*th)/(tl*tnd+tlm*tnc+tl*tnh);
            rbvh = (tl*td+tl*tc+tlm*th)/(tl*tnd+tl*tnc+tlm*tnh);
    }

    // output posterior + estimated fraction of data used for each cluster
    sp.posterior.aa.Set(rbma,rbva,rnc+sp.prior.aa.k,rnc+sp.prior.aa.v,sp.prior.aa.ym,sp.prior.aa.yss,sp.prior.aa.xyss);
    sp.posterior.ab.Set(rbmh,rbvh,rnh+sp.prior.ab.k,rnh+sp.prior.ab.v,sp.prior.ab.ym,sp.prior.ab.yss,sp.prior.ab.xyss);
    sp.posterior.bb.Set(rbmb,rbvb,rnd+sp.prior.bb.k,rnd+sp.prior.bb.v,sp.prior.bb.ym,sp.prior.bb.yss,sp.prior.bb.xyss);
    // set up posterior covariances to be prior covariances
    // since working with inverse, not actually changed by increased observations
    sp.posterior.xah = sp.prior.xah;
    sp.posterior.xab = sp.prior.xab;
    sp.posterior.xhb = sp.prior.xhb;
}


/**
 * This penalizes unusual frequency distributions
 * 
 * @param a - first cluster freq
 * @param b - second cluster freq
 * @param c - third cluster freq
 * @param oa - number of observations
 * @param ob - number of observations
 * @param oc - number of observations
 * @param lambda - parameterized safety value to avoid zero frequency
 */
double mixture_penalty(double a,double b, double c,double oa, double ob, double oc,double lambda)
{
    double total, entropy;

    // penalty for mixing fraction
    // want to avoid splitting clusters
    // so 20/0/0 is good, 7/7/7 less good
    // i.e. need more evidence to provoke a new cluster
    total = a+b+c+3*lambda;
    entropy = oa * log((a+lambda)/total)+ob*log((b+lambda)/total)+oc*log((c+lambda)/total);
    entropy = -1*entropy;
    return(entropy);
}

/**
 * This penalizes unusual cluster distributions (violating HW extremely)
 * Expected to be useful for screening/finding initial priors
 *
 * @param a - hom cluster freq
 * @param b - hom cluster freq
 * @param c - het cluster freq
 */
double HW_penalty(double a, double b, double c, double lambda)
{
    double total, penalty, p, q;

    total = 2*a+2*b+2*c+2*lambda; // total number of chromosomes +1 for always potential errors
    p = (2*a+c+lambda)/total; // hom chromosomes + het chromosomes belonging to the allele
    q = 1-p;              // other hom
    // a* log(p*p)+c*log(2*p*q) + b*log(q*q) = (2a+c)*log(p)+(2b+c)*log(q)+c*log(2)
    penalty = (2*a+c)*log(p) + (2*b+c)*log(q)+c*log((double)2);
    penalty = -1*penalty;
    return(penalty);
}


/**
 * This is the "labeling" routine that tries out all plausible genotypes for best likelihood
 * Basic goal: efficiently evaluate likelihood of three clusters based on seed
 * Loops over BB,...BB (i times), AB,...AB (j-i), AA,...AA (n-j+1)
 *
 * @param ec - relative probability of c = a allele dominates
 * @param eh - relative probability of h = ab allele similar
 * @param ec - relative probability of d = b allele dominates
 * @param ZS - number of bins + 1 = size of data
 * @param zall - data in bins, sum, sum of squares
 * @param sp - snp parameters
 */
void integratebrlmmoverlabelings(vector<double> &ec, vector<double> &eh, vector<double> &ed, int ZS, binned_data &zall,snp_param &sp)
{
    int ZSS;
    int NS;
    int i,j;
    //double tmp;
    // allocate memory for these puppies
    // over allocate for most of these, but temporary
    ZSS = ZS*ZS;

    // likelihood of data, given 
    vector<double> q;
    q.resize(ZSS);
    q.assign(ZSS,0);
    vector<double> ddx;
    vector<double> ccx;
    vector<double> nd;
    vector<double> nc;
    vector<double> dd;
    vector<double> cc;
    ddx.resize(ZS);
    ccx.resize(ZS);
    nd.resize(ZS);
    nc.resize(ZS);
    dd.resize(ZS);
    cc.resize(ZS);

    vector<double> dnod;
    vector<double> dnoh;
    vector<double> dnoc,cnoc;
    dnod.resize(ZS);
    //hnod.resize(ZS);
    //cnod.resize(ZS);
    dnoh.resize(ZS);
    //hnoh.resize(ZS);
    //cnoh.resize(ZS);
    dnoc.resize(ZS);
    //hnoc.resize(ZS);
    cnoc.resize(ZS);

    // count up the number of contradictions
    // cumulative contradictions
    vcumsum(dnod,zall.nod,ZS);
    vcumsum(dnoh,zall.noh,ZS);
    vcumsum(dnoc,zall.noc,ZS);
    //vreversecumsum(cnod,dnod,ZS);
    //vreversecumsum(cnoh,dnoh,ZS);
    vreversecumsum(cnoc,dnoc,ZS);
    // multiply by the penalty to reduce computational load
    // obsolete: contradictions are directly penalized in the bins
    //vmult(cnoc,sp.contradictionpenalty);
    //vmult(dnod,sp.contradictionpenalty);
    //vmult(dnoh,sp.contradictionpenalty);

    

    // cumulative number - weight for bins
    vcumsum(nd,zall.nz,ZS);
    NS = (int)floor(nd[ZS-1]+.00001); // number of points total
    vreversecumsum(nc,nd,ZS);
    //cumulative sum of values
    vcumsum(dd,zall.zz,ZS);
    vreversecumsum(cc,dd,ZS);
    // cumulative sum squares
    vcumsum(ddx,zall.zzx,ZS);
    vreversecumsum(ccx,ddx,ZS);

    int tin;

    // temporaries
    double ycij,xcij,xcxij;
    double ydij,xdij,xdxij;
    double yhij,xhij,xhxij;
    double ldij,lcij,lhij;
    double dvdij,dvcij,dvhij;

    // setup brlmm bits that don't change
    // speed
            vector<double> m;
            vector<double> Minv;
            vector<double> Sinv;
            setupBRLMM(m,Minv,Sinv,sp);
            // only work with Sinv;
            // scratch vectors
            vector<double> Tz;
            vector<double> Zt;
            vector<double> Xinv;
            vector<double> Xt;
            vector<double> Mz;
            // vectors that change
            // sum of x
            vector<double> Nv;
            // number of data points
            vector<double> N;
            Nv.resize(3);
            N.assign(9,0);

    // do some more setup
    // so that things that aren't executed don't get repeated conditionals
    if (sp.copynumber<2)
    {
        // no hets
        for (i=0; i<ZS; i++)
            for (j=i+1; j<ZS; j++)
            {
                tin = hax(i,j,ZS);
                q[tin] += NS*1000;
            }
    }

    if (sp.hints==1)
    {
        // penalty for contradiction of reference
        for (i=0; i<ZS; i++)
            for (j=i; j<ZS; j++)
            {
                tin = hax(i,j,ZS);
                q[tin] += (dnod[i]+cnoc[j]+(dnoh[j]-dnoh[i]));
            }

    }
    
    if (sp.mix==1)
    {
        // mixture penalty under observed data
        for (i=0; i<ZS; i++)
            for (j=i; j<ZS; j++)
            {
                tin = hax(i,j,ZS);
                q[tin] += mixture_penalty(nc[j],nd[i],nd[j]-nd[i],nc[j],nd[i],nd[j]-nd[i],sp.SafetyFrequency);
            }
    }
    if (sp.mix==2)
    {
        // mixture penalty under observed+prior number of observations
        // change to use frequency not "prior precision"
        // make room for better code
        // i.e. posterior should be grouped correctly
        for (i=0; i<ZS; i++)
            for (j=i; j<ZS; j++)
            {
                tin = hax(i,j,ZS);
                q[tin] += mixture_penalty(nc[j]+sp.prior.aa.k,nd[i]+sp.prior.bb.k,nd[j]-nd[i]+sp.prior.ab.k,nc[j],nd[i],nd[j]-nd[i],sp.SafetyFrequency);
            }
    }
    if (sp.mix==3)
    {
        for (i=0; i<ZS; i++)
            for (j=i; j<ZS; j++)
            {
                tin = hax(i,j,ZS);
                q[tin] += HW_penalty(nc[j],nd[i],nd[j]-nd[i], 0.5*sp.SafetyFrequency);
            }
    }
    if (sp.bic>0)
    {
        // if useful to do some bic correction
        // bic * log(n) - empirical correction
        for (i=0; i<ZS; i++)
        {
            // subtract off once for 2 cluster
            // [0,0],[0,ZS-1],[ZS-1,ZS-1] occur twice
            // so subtract off 2 times for 1 cluster
            tin = hax(i,i,ZS);
                    q[tin] -= sp.bic*log((double)NS);
            tin = hax(0,i,ZS);
            q[tin] -= sp.bic*log((double)NS);
            tin = hax(i,ZS-1,ZS);
            q[tin] -= sp.bic*log((double)NS);
            for (j=i; j<ZS; j++)
            {
                tin = hax(i,j,ZS);
                // add 3 times for everyone
                q[tin] += 3*sp.bic*log((double)NS);
            }
        }
    }
    // likelihood of data under posterior
    // (x-m)^2 = x^2-2*x*m+m^2
    // running sums sum(x^2)-2*(sum(x)*m+sum(1)*m^2
    // big loop: try all reasonable labelings
    // bb^i->ab^(j-i)->aa^(n-j+1)
    for (i=0; i<ZS; i++)
        for (j=i; j<ZS; j++)
        {
            tin = hax(i,j,ZS);
            // initialize global matrix
            // c
            ycij = nc[j];
            xcij = cc[j];
            xcxij = ccx[j];
            // d
            ydij = nd[i];
            xdij = dd[i];
            xdxij = ddx[i];
            // h
            yhij = nd[j]-nd[i];
            xhij = dd[j]-dd[i];
            xhxij = ddx[j]-ddx[i];


            // brlmm here
            // setup centers
            Nv[0] = xcij;
            Nv[1] = xhij;
            Nv[2] = xdij;
            // number of data points
            N[0] = ycij;
            N[4] = yhij;
            N[8] = ydij;
            // do some inversions
            timesthree(Tz,N,Sinv);
            sumthree(Zt, Minv,Tz);
            inversethree(Xinv,Zt);
            timesone(Tz,Sinv, Nv);
            timesone(Zt,Minv, m);
            sumthree(Xt,Tz,Zt);
            timesone(Mz,Xinv,Xt);
            // Tz = three means
            //printf("%d %f %f %f\n", tin, Mz[0],Mz[1],Mz[2]);
            
            
            // figure out the triple mean/variance combo
            ldij = Mz[2];
            lhij = Mz[1];
            lcij = Mz[0];
            // force means apart by shellbarrier
            if (sp.hardshell==3)
                force_isotonic(ldij,lhij,lcij,ydij+sp.prior.bb.k,yhij+sp.prior.ab.k,ydij+sp.prior.aa.k, sp.shellbarrier);
            // compute variance given means
            // within-cluster variance + prior + distance shifted
            dvcij = b_var(sp.prior.aa.v,sp.prior.aa.ss,xcxij,xcij,ycij,sp.prior.aa.k,lcij,sp.prior.aa.m);
            dvdij = b_var(sp.prior.bb.v,sp.prior.bb.ss,xdxij,xdij,ydij,sp.prior.bb.k,ldij,sp.prior.bb.m);
            dvhij = b_var(sp.prior.ab.v,sp.prior.ab.ss,xhxij,xhij,yhij,sp.prior.ab.k,lhij,sp.prior.ab.m);

            // fix up common variance
            if (sp.comvar==1)
            {
                double tc,td,th;
                double tnc,tnd,tnh;
                double tl,tlm;
                tnd = ydij+sp.prior.bb.v;
                tnc = ycij+sp.prior.aa.v;
                tnh = yhij+sp.prior.ab.v;
                td = dvdij*tnd;
                tc = dvcij*tnc;
                th = dvhij*tnh;
                // cross-talk between variances 0 = independent, 1 = common variance
                tl = sp.lambda;
                tlm = 3-2*sp.lambda;
                dvdij = (tlm*td+tl*tc+tl*th)/(tlm*tnd+tl*tnc+tl*tnh);
                dvcij = (tl*td+tlm*tc+tl*th)/(tl*tnd+tlm*tnc+tl*tnh);
                dvhij = (tl*td+tl*tc+tlm*th)/(tl*tnd+tl*tnc+tlm*tnh);
            }
        ////
        ////
        //// now do the global likelihood with calculated values
            // likelihood of the data gvien mean/variance    
            //q[tin] = 0;
            q[tin] += l_posterior(xdxij,ldij,xdij,ydij,dvdij);
            q[tin] += l_posterior(xhxij,lhij,xhij,yhij,dvhij);
            q[tin] += l_posterior(xcxij,lcij,xcij,ycij,dvcij);    
            //printf("%d %lf\n", i, q[i]);    
    
            // likelihood of mean under prior
            //q[tin] += l_normal(sp.bb.m,ldij,sp.bb.ss/sp.bb.k);
            //q[tin] += l_normal(sp.ab.m,lhij,sp.ab.ss/sp.ab.k);
            //q[tin] += l_normal(sp.aa.m,lcij,sp.aa.ss/sp.aa.k);
            q[tin] += l_normal(sp.prior.bb.m,ldij,1/Minv[8]);
            q[tin] += l_normal(sp.prior.ab.m,lhij,1/Minv[4]);
            q[tin] += l_normal(sp.prior.aa.m,lcij,1/Minv[0]);
            // likelihood of variance under prior
            q[tin] += l_inverse(sp.prior.bb.ss,dvdij,sp.prior.bb.v);
            q[tin] += l_inverse(sp.prior.ab.ss,dvhij,sp.prior.ab.v);
            q[tin] += l_inverse(sp.prior.aa.ss,dvcij,sp.prior.aa.v);
            // half
            q[tin] /= 2;
            
            // penalize assorted peculiar features of the data
            // i.e. too close, too strange
            if (sp.hardshell==2)
                q[tin] += NS*fix_shell(ldij,lcij,lhij,sp.shellbarrier);
            if (sp.hardshell==1)
                q[tin] += NS*fix_hom_shell(ldij,lcij,lhij,sp.shellbarrier);
            if (sp.CSepPen>0)
            {
                // avoid cluster splitting by ad-hoc penalty for FLD too small
                // rather ad-hoc favor for having FLD large!
                double flddh,fldhc,flddc;
                // compute FLD for cluster separation
                flddh = (ldij-lhij)*(ldij-lhij)/(dvdij+dvhij);
                fldhc = (lhij-lcij)*(lhij-lcij)/(dvhij+dvcij);
                flddc = (ldij-lcij)*(ldij-lcij)/(dvdij+dvcij);
                // threshold FLD to avoid too much repulsion
                flddh = flddh/(1+flddh/sp.CSepThr);
                fldhc = fldhc/(1+fldhc/sp.CSepThr);
                // might have only homs, so keep this
                flddc = flddc/(1+flddc/(2*sp.CSepThr));
                // apply data points at risk
                flddh *= (ydij+yhij);
                fldhc *= (yhij+ycij);
                flddc *= (ydij+ycij);
                // favor! [larger FLD better, up to a point]
                q[tin] -= sp.CSepPen * (flddh+fldhc+flddc);
            }
        }

    
    qtoeq(q,ZS);
    // convert q to label
    eqtoproblabel(ec,eh,ed,ZS,q);
}

/**
 * tunes up the prior to allow for cluster shifts
 *
 */
void apply_wobble(snp_param &sp)
{
    // apply "wobble" to account for possible cluster shifts
    // i.e. even if infinite amount of prior data
    // this experimental set could be different
    //
    // // change to "f" here to preserve prior relative frequency after changing k
    //
    // note that this is currently "undone" in the actual 1-D seeding routine to preserve continuity of behavior
    // flag will be installed to allow better behavior
    
    double tmp;
    tmp = 1/sp.wobble;
    if (tmp<sp.prior.aa.k)
        sp.prior.aa.k = tmp;
    if (tmp<sp.prior.ab.k)
        sp.prior.ab.k = tmp;
    if (tmp<sp.prior.bb.k)
        sp.prior.bb.k = tmp;
    // note that the posterior will reflect the "lower" preknowledge
    // of the prior data
    
}

/*
 * Do QC check to see if values are weird, indicative of error
 *
 * @param call - genotype per data point
 * @param conf - confidence per data point
 * @param y - strength per data point
 */
void qc_modify_confidences(vector<affx::GType> &call, vector<double> &conf, vector<double> &y, double perr)
{
    // use scatter in y direction to find outlier data points
    // outlier = pr(under model for y)<pr(error)
    // assume normal distribution of residuals
    // assume common variance (improve power for minor cluster)
    // don't assume common mean
    int cIx, i;
    double tmp,num,ss;
    vector<double> residual;
    double wgt;


    residual.resize(y.size());
    for (cIx=0; cIx<3; cIx++) // valid calls
    {
        tmp = 0;
        num = 0.00001; // tiny
        for (i=0; i<call.size(); i++)
            if (call[i]==cIx)
            {
                tmp += y[i];
                num +=1;
            }
        tmp/=num; // now have mean
        for (i=0; i<call.size(); i++)
            if (call[i]==cIx)
            {
                residual[i] = y[i]-tmp;
            }
    }
    ss = 0;
    for (i=0; i<residual.size(); i++)
        ss+= residual[i]*residual[i];
    ss /= residual.size(); // degrees of freedom close enough
    ss += 0.001; // avoid zero
    for (i=0; i<residual.size(); i++)
    {
        // "standardized" probability (don't care about variance)
        tmp = exp(-1*residual[i]*residual[i]/(2*ss))/(2.51);
        // generate a weight
        if (tmp>0)
            wgt = perr/tmp;
        else
            wgt = perr/(0.000001);
        conf[i] = (conf[i] + wgt)/(1+wgt); // perr = "relative likelihood of error"
    }
}

void qc_compute_copy_posterior(vector<affx::GType> call, vector<double> &y, snp_distribution &posterior)
{
    vector<double> cmean;
    vector<double> cnum;
    vector<double> residual;
    double shrinkage = 0.00001;
    double vshrinkage = 0.1;
    double minvar = 0.1;
    int cIx, i;
    double cy,ss;
    
    // old_style computation
    residual.resize(call.size());
    cmean.resize(3);
    cnum.resize(3);
    
    for (cIx=0; cIx<3; cIx++)
    {
        cmean[cIx] = 0;
        cnum[cIx] = 0;
        for (i=0; i<call.size(); i++)
        {
            if (call[i]==cIx)
            {
                cmean[cIx] += y[i];
                cnum[cIx] +=1;
            }
        }
    }
    cy = (cmean[0] + cmean[1] + cmean[2])/(cnum[0]+cnum[1]+cnum[2]);
    
    // shrink towards common mean to fill in nonexistent calls
    cmean[0] = (cy*shrinkage+cmean[0])/(cnum[0]+shrinkage);
    cmean[1] = (cy*shrinkage+cmean[1])/(cnum[1]+shrinkage);
    cmean[2] = (cy*shrinkage+cmean[2])/(cnum[2]+shrinkage);
    // all three exist
    for (i=0; i<call.size(); i++)
    {
        // if (call[i]>-1 && call[i]<3) {
    if (affx::GType_called(call[i])) { // aa,ab,bb
            residual[i] = y[i] - cmean[call[i]];
        }
        else {
            residual[i] = 0.0;
    }
    }
    ss = 0;
    for (i=0; i<residual.size(); i++)
        ss+= residual[i]*residual[i];
    ss = (ss + minvar*vshrinkage)/(residual.size()+vshrinkage); // degrees of freedom close enough for large data sets

    // output means (shrunk towards center)
    posterior.aa.ym = cmean[0];
    posterior.ab.ym = cmean[1];
    posterior.bb.ym = cmean[2];
    // output common variance
    posterior.aa.yss = ss;
    posterior.ab.yss = ss;
    posterior.bb.yss = ss;
}


/*
 * Do QC check to see if values are weird, indicative of error
 *
 * @param call - genotype per data point
 * @param conf - confidence per data point
 * @param y - strength per data point
 */
void qc_posterior_modify_confidences(vector<affx::GType> &call, vector<double> &conf, vector<double> &y, double perr, snp_distribution dist)
{
    // use scatter in y direction to find outlier data points
    // outlier = pr(under model for y)<pr(error)
    // assume normal distribution of residuals
    // assume common variance (improve power for minor cluster)
    // don't assume common mean
    int cIx, i;
    double tmp;
  // double num,ss;
    double residual;
    vector<double> cmean,cvar;
    double wgt;

    cmean.resize(3);
    cvar.resize(3);
    cmean[0] = dist.aa.ym;
    cmean[1] = dist.ab.ym;
    cmean[2] = dist.bb.ym;

    cvar[0] = dist.aa.yss;
    cvar[1] = dist.ab.yss;
    cvar[2] = dist.bb.yss;
    
    for (cIx=0; cIx<3; cIx++) // valid calls
    {
        for (i=0; i<call.size(); i++)
            if (call[i]==cIx)
            {
                residual = y[i]-cmean[cIx];
                tmp = exp(-1*residual*residual/(2*cvar[cIx]))/(2.51);
                // generate a weight
                if (tmp>0)
                {
                    wgt = perr/tmp;
                    conf[i] = (conf[i] + wgt)/(1+wgt); // perr = "relative likelihood of error"
                }
                else
                    conf[i] = 1;            
            }
    }
}

/**
 * Set up the binned data from the regular data
 * compute running sums, allocate sorted bins
 *
 * @param zall - overall structure to keep track of data + soft genotype labels by bin
 * @param x    - initial data points (unsorted)
 * @param genohints - hints associated with data points
 * @param numbins - how many bins to make, if bins>points, go with points
 *
 */
void initialize_bins(binned_data &zall, vector<double> &x, vector<double> &y, vector<int> &genohints, vector<double> &SubsetInbred, int numbins, int hok, double CP)
{
    int i;

    // takes a sorted vector of length z and labels it
    // this will be different when using bins



    // call sides "d" and "c" to avoid genotype confusion
    // negative genotypes are "d", positive are "c", middle are "h"

    // divide sorted points into bins by z
    // division by range(points)/(bins+1)
    // if bins==0, no change in behavior, exactly 1 point/bin
    // min 1 point per bin (nz)
    // sum(z) in zz
    // sum(z^2) in zzx
    // count the bins
    //     
    zall.length = x.size();
    zall.NS = zall.length+1;
    
    vector<double> z,w;
    vector<int> zerohints;
    vector<double> zeroInbred;
    // local hints always null unless assigned otherwise
    zerohints.assign(zall.length,-1);
    zeroInbred.assign(zall.length,0);
    /// store sorted data
    z.resize(zall.length);
    w.resize(zall.length); // y values stored in sorted x order
    ///
    zall.nz.resize(zall.NS);
    zall.zz.resize(zall.NS);
    zall.zzx.resize(zall.NS);
    // y terms
    zall.zy.resize(zall.NS);
    zall.zzy.resize(zall.NS);
    zall.zxy.resize(zall.NS);
    // hints
    zall.nod.resize(zall.NS);
    zall.noh.resize(zall.NS);
    zall.noc.resize(zall.NS);

    vector< pair<int,double> > data;

    for (i=0; i<zall.length; i++)
    {
        pair<int,double> p;
        /*
        if (genohints.size()>0)
            p.first = genohints[i];
        else
            p.first = -1; // no hint here */
        p.first = i; // keep index
        p.second = x[i];
        data.push_back(p);
    }

    // sort the data keeping the index to refer to others
    sort(data.begin(),data.end(),second_less());

    for (i=0; i<zall.length; i++)
    {
        z[i] = x[data[i].first];
        w[i] = y[data[i].first];
        if (genohints.size()>0)
            zerohints[i] = genohints.at(data[i].first);
        else
            zerohints[i] = -1; //no hint here
        if (SubsetInbred.size()>0)
            zeroInbred[i] = SubsetInbred.at(data[i].first);
    }
    zall.ZS = setup_bins(zall,z,w,zerohints, zeroInbred, zall.length,numbins,hok,CP);
    zall.ec.resize(zall.ZS);
    zall.ed.resize(zall.ZS);
    zall.eh.resize(zall.ZS);    
}


/**
 * The master routine: takes data and turns it into genotypes with no seeds needed
 *
 * @param call - output of genotypes per data point
 * @param conf - output of confidences (1-relative probability)
 * @param x    - contrast values per data point
 * @param y     - size values per data point - used as QC check for copy number
 * @param genohints - any reference genotypes you have lying around as hints
 * @param sp - snp parameters: prior information, algorithm parameters, etc.
 */
void bayes_label(vector<affx::GType> &call, 
                 vector<double> &conf, 
                 vector<double> &x, 
                 vector<double> &y, 
                 vector<int> &genohints, 
                 vector<double> &SubsetInbred, 
                 snp_param &sp)
{
  vector<double> ted,teh,tec;
  bayes_label(call, conf, x, y, genohints, SubsetInbred, sp, ted, teh, tec, false);
}
void bayes_label(vector<affx::GType> &call, 
                 vector<double> &conf, 
                 vector<double> &x, 
                 vector<double> &y, 
                 vector<int> &genohints, 
                 vector<double> &SubsetInbred, 
                 snp_param &sp, 
                 vector<double> &ted,
                 vector<double> &teh,
                 vector<double> &tec, 
                 bool store_probabilities)
{
    binned_data zall;

    // set up posterior in case no processing of data ("single sample mode")
    sp.priortoposterior();
    // weaken prior if cluster shifting expected
    apply_wobble(sp);

    // if we're not running single-sample mode
    // we want to fit the binned data
    if (sp.callmethod<2)
    {
        // format the data & hints into bins
        initialize_bins(zall, x, y, genohints, SubsetInbred, sp.bins,sp.hok,sp.contradictionpenalty);

        // take the setup-bins, integrate over all possible labelings
        // computing the relative likelihood of the data    
        integratebrlmmoverlabelings(zall.ec,zall.eh,zall.ed,zall.ZS,zall,sp);

        // take the estimated mixture fraction of labels
        // turn into approximate "posterior" distribution 
        if (sp.clustertype==1)
            labelstoposterior(zall.ec,zall.eh,zall.ed,zall.ZS,zall,sp);
        else
            labels_two_posterior(zall.ec,zall.eh,zall.ed,zall.ZS,zall,sp);
    }
    // make calls based on mean/variance/confidence
    // x >unsorted< matrix of calls
    if (sp.callmethod>0)
    {
        if (store_probabilities) {
            ted.resize(x.size());
            teh.resize(x.size());
            tec.resize(x.size());
        }
        // make calls only by posterior distribution
        if (sp.clustertype==1)
            make_calls(call,conf,x,x.size(),sp.posterior.aa, sp.posterior.ab, sp.posterior.bb, sp.copynumber,sp.mix,ted,teh,tec,store_probabilities);
        else
            make_two_calls(call,conf,x,y,x.size(),sp.posterior, sp.copynumber,sp.mix,sp.ocean,sp.inflatePRA,ted,teh,tec, store_probabilities);
    }
    else
    {
        // output the raw labels matched to data points
        ted.resize(zall.length);
        tec.resize(zall.length);
        teh.resize(zall.length);
        // expand data from bins into 1-1 data (still sorted)
        expand_data(ted,zall.length,zall.ed,zall.nz,zall.ZS);
        expand_data(tec,zall.length,zall.ec,zall.nz,zall.ZS);
        expand_data(teh,zall.length,zall.eh,zall.nz,zall.ZS);
        // make calls by expanded data (unsorting)
        make_label_calls(call,conf,x,x.size(),ted,tec,teh,sp.copynumber);
        if (!store_probabilities) {
            ted.clear();
            teh.clear();
            tec.clear();
        }
    }
    // do copy number QC test
    // only if not running in single-sample mode
    if (sp.copytype==0)
    {
        if (sp.copyqc>0 && sp.callmethod<2)
        {
            qc_modify_confidences(call,conf,y,sp.copyqc);
        }
    }
    // new style copyqc
    // uses prior information if single-sample mode to modify confidence
    // Note possible low confidence for unobserved clusters!!!
    // must have 2-dimensional clusters stored/allowed
    if (sp.copytype==1)
    {
        if (sp.copyqc>0)
        {
            // do not compute posterior if single-sample mode
            if (sp.callmethod<2)
                qc_compute_copy_posterior(call,y,sp.posterior);
            // modify confidences if either mode
            qc_posterior_modify_confidences(call,conf,y,sp.copyqc,sp.posterior);
        }
    }
}    
