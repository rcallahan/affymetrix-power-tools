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
 * @file   snp.label.h
 * @author Earl Hubbell 
 * @date   Jan 22 2007
 * 
 * 
 */

#ifndef _SNPLABEL_H_
#define _SNPLABEL_H_

//
#include "chipstream/BioTypes.h"
#include "chipstream/SelfDoc.h"
//
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>

/**
 *  this class summarizes a single genotype
 */
class cluster_data {

public:
	// mean/variance
	// strength of both
	double m; ///< mean of cluster
	double ss; ///< variance of cluster
	double k; ///< strength of mean (pseudo-observations)
	double v; ///< strength of variance (pseudo-observations)
	double ym; ///< mean of cluster in other dimension
	double yss; ///< variance of cluster in other dimension
	double xyss; ///< covariance of cluster in both directions

	/** set the values */
	void Set(double im, double iss, double ik, double iv, double iym, double iyss, double ixyss){
		m=im;
		ss = iss;
		k=ik;
		v = iv;
		ym = iym;
		yss = iyss;
		xyss = ixyss;
	};
  	/** copy from another cluster */
	void Copy(cluster_data &source){
		m = source.m;
		ss = source.ss;
		k = source.k;
		v = source.v;
		ym = source.ym;
		yss = source.yss;
		xyss = source.xyss;
	};
 	/** do a diagnostic dump */
	void Dump(){
		printf("%lf %lf %lf %lf %lf %lf %lf\n",m,ss,k,v, ym, yss, xyss);
	}
  void Clear(){
		m    = 0.0;
		ss   = 0.0;
		k    = 0.0;
		v    = 0.0;
		ym   = 0.0;
		yss  = 0.0;
		xyss = 0.0;
  };
};

/** this class summarizes binned data */
class binned_data{
public:
	// length of data
	int length;
	int NS;
	int ZS;
	//bin summaries
	std::vector<double> nz; ///< number of things in each bin
	std::vector<double> zz; ///< sum of things in each bin
	std::vector<double> zzx; ///< sum of squares in each bin
	// bin summaries y
	std::vector<double> zy; ///< sum of y values in each bin
	std::vector<double> zzy; ///< sum of squares of y values in each bin
	// cross terms
	std::vector<double> zxy; ///< sum of x*y values in each bin
	// hints: counts of genotypes per bin
	std::vector<double> nod; ///< penalty*number of hint genotypes not d
	std::vector<double> noh; ///< penalty*number of hint genotypes not h
	std::vector<double> noc; ///< penalty*number of hint genotypes not c
	// soft genotypes by bin
	std::vector<double> ec;
	std::vector<double> ed;
	std::vector<double> eh;
};

/** this class holds all three genotype clusters */
class snp_distribution{
public:
	// prior
	// note: a-b for transform
	// therefore bb allele is negative
	// despite allegations to the contrary
	cluster_data aa; ///< aa allele cluster
	cluster_data ab; ///< ab genotype cluster
	cluster_data bb; ///< bb genotype cluster

	// cross-correlation terms
	// Like everything else, in pseudo-obs
	// Can be asymmetrical because of rotation between variances within clusters
	double xah,xab,xhb; ///< cross-correlation terms between clusters x
	double yah,yab,yhb; ///< cross-correlation terms between clusters y
	double xyah,xyab,xyhb; ///< cross-correlation terms between clusters x and y
	double yxah,yxab,yxhb; ///< cross-correlation terms between clusters y and x

	void Copy(snp_distribution &s)
	{
		aa.Copy(s.aa);
		ab.Copy(s.ab);
		bb.Copy(s.bb);
		xah = s.xah;
		xab = s.xab;
		xhb = s.xhb;
		yah = s.yah;
		yab = s.yab;
		yhb = s.yhb;
		xyah = s.xyah;
		xyab = s.xyab;
		xyhb = s.xyhb;
		yxah = s.yxah;
		yxab = s.yxab;
		yxhb = s.yxhb;
	}
  void Clear_Cross() {
		xah  = 0.0;
		xab  = 0.0;
		xhb  = 0.0;
		yah  = 0.0;
		yab  = 0.0;
		yhb  = 0.0;
		xyah = 0.0;
		xyab = 0.0;
		xyhb = 0.0;
		yxah = 0.0;
		yxab = 0.0;
		yxhb = 0.0;
  };
  void Clear() {
    aa.Clear();
    ab.Clear();
    bb.Clear();
    Clear_Cross();
  };
};

/** this holds all the parameters for a snp, including alg choices */
class snp_param{
public:

	snp_distribution prior;
	snp_distribution posterior;
	
	// hints
	int comvar; ///< do we use common variances for all clusters?
	double lambda; ///< control mixing of variances if common
	int callmethod; ///< do we use raw labels or posterior distributions?
	int hardshell;  ///< do we stop cluster centers from being too close?
	double shellbarrier; ///< how close can cluster centers be?
	int bins;  ///< use quick method
	int copynumber; ///< only one copy, for example, XY snps
	int hints; ///< use hinted genotypes as reference
	double contradictionpenalty; ///< penalty for contradicting a hint
	int hok; ///< allow hints to be flipped in genotype 0,1,2<->2,1,0
	int mix; ///< turn on mixture distribution penalty
	double bic; ///< bic penalty level bic*nparam*log(n)
	double CSepPen; ///< penalty for FLD separation too low
	double CSepThr; ///< stop penalizing using Geman-McClure
	double wobble; ///< how much do we allow clusters to shift from prior
	double copyqc; ///< check for possible outlier "size" values
	int copytype; ///< how to handle copy number outlier detection
	int clustertype; ///< do I handle multiple cluster dimensions
	double ocean; ///< uniform density to compare weird data points against
	double inflatePRA; ///< increase variance by uncertainty in mean location when making calls
	double IsoHetY; ///< ensure Het > average Hom Y by at least epsilon

	// safetyparameter
	double SafetyFrequency; ///< ensure avoid taking log of zero frequencies

	/** initialize to some useful defaults */
	void Initialize(){
		prior.aa.Set(.66,.005,4,10, 9, 0.1, 0);
		prior.ab.Set(0,.01,.2,10,9,0.1,0);
		prior.bb.Set(-.66,.005,4,10,9,0.1,0);
		prior.xah=prior.xab=prior.xhb=0;
		prior.yah=prior.yab=prior.yhb=0;
		prior.xyah=prior.xyab=prior.xyhb=0;
		prior.yxah=prior.yxab=prior.yxhb=0;
		posterior.Copy(prior);

		comvar=1;
		lambda = 1;
		callmethod = 0; // labeling
		hardshell=2;
		shellbarrier = .05;
		bins = 0; // no bins used
		copynumber = 2; // default two copies
		hints = 0;
		contradictionpenalty = 0; // no penalty
		hok = 0;
		mix = 0; // no mixture penalty
		bic = 0; // no bic turned on
		wobble = 0.00001; // max prior 100,000 data points
		copyqc = 0; // no test for copy qc
		copytype = 0; // standard copy qc method
		clustertype = 1; // standard 1-d clustering only
		CSepPen = 0;
		CSepThr = 16;
		ocean = 0; // similar to copyqc, turn this off normally
		inflatePRA = 0; // no increase in uncertainty
		IsoHetY = 0; // nothing here
		SafetyFrequency= 1; // at least one count per genotype safety
	};
 	/** copy parameters */
	void copy(snp_param s){
		// copy prior
		prior.Copy(s.prior);

		// copy posterior
		posterior.Copy(s.posterior);


		// copy algorithm tweaks/parameters
		comvar=s.comvar;
		lambda=s.lambda;
		hardshell=s.hardshell;
		shellbarrier=s.shellbarrier;
		callmethod = s.callmethod;
		copynumber = s.copynumber;
		bins = s.bins;
		hints = s.hints;
		contradictionpenalty = s.contradictionpenalty;
		hok = s.hok;
		mix = s.mix;
		bic = s.bic;
		wobble = s.wobble;
		copyqc = s.copyqc;
		copytype = s.copytype;
		clustertype = s.clustertype;
		CSepPen = s.CSepPen;
		CSepThr = s.CSepThr;
		ocean = s.ocean;
		inflatePRA = s.inflatePRA;
		IsoHetY = s.IsoHetY;
		SafetyFrequency = s.SafetyFrequency;
	};
	/** get only the prior value from some source*/
	void getprior(snp_param s){
		// just copy the prior data over
		// keep the current methods
		prior.Copy(s.prior);
	};
	/** transfer prior directly to posterior */
	void priortoposterior(){
		posterior.Copy(prior);
	}
    void setDocValues(SelfDoc &doc) {
        // prior strength
        doc.setOptValue("KX",ToStr(prior.aa.k));
        doc.setOptValue("KH",ToStr(prior.ab.k));
        doc.setOptValue("V",ToStr(prior.ab.v));

        // prior centers
        doc.setOptValue("BBM",ToStr(prior.bb.m));
        doc.setOptValue("ABM",ToStr(prior.ab.m));
        doc.setOptValue("AAM",ToStr(prior.aa.m));

	doc.setOptValue("BBY",ToStr(prior.bb.ym));
	doc.setOptValue("ABY",ToStr(prior.ab.ym));
	doc.setOptValue("AAY",ToStr(prior.aa.ym));

        // prior variances
        doc.setOptValue("AAV",ToStr(prior.aa.ss));
        doc.setOptValue("BBV",ToStr(prior.bb.ss));
        doc.setOptValue("ABV",ToStr(prior.ab.ss));

	doc.setOptValue("AAYV",ToStr(prior.aa.yss));
	doc.setOptValue("BBYV",ToStr(prior.bb.yss));
	doc.setOptValue("ABYV",ToStr(prior.ab.yss));

	doc.setOptValue("AAXY",ToStr(prior.aa.xyss));
	doc.setOptValue("ABXY",ToStr(prior.ab.xyss));
	doc.setOptValue("BBXY",ToStr(prior.bb.xyss));

	//Cryptic between-cluster terms
        doc.setOptValue("KXX",ToStr(prior.xab));
        doc.setOptValue("KAH",ToStr(prior.xah));
        doc.setOptValue("KHB",ToStr(prior.xhb));
 	// Y versions
	doc.setOptValue("KYAB",ToStr(prior.yab));
	doc.setOptValue("KYAH",ToStr(prior.yah));
	doc.setOptValue("KYHB",ToStr(prior.yhb));

	// other parameters
        doc.setOptValue("COMVAR",ToStr(comvar));
        doc.setOptValue("HARD",ToStr(hardshell));
        doc.setOptValue("SB",ToStr(shellbarrier));
        doc.setOptValue("CM",ToStr(callmethod));
        doc.setOptValue("bins",ToStr(bins));
        doc.setOptValue("hints",ToStr(hints));
        doc.setOptValue("CP",ToStr(contradictionpenalty));
	doc.setOptValue("Hok",ToStr(hok));
        doc.setOptValue("mix",ToStr(mix));
        doc.setOptValue("bic",ToStr(bic));
	doc.setOptValue("CSepPen",ToStr(CSepPen));
	doc.setOptValue("CSepThr",ToStr(CSepThr));
        doc.setOptValue("lambda",ToStr(lambda));
        doc.setOptValue("wobble",ToStr(wobble));
        doc.setOptValue("copyqc",ToStr(copyqc));
	doc.setOptValue("copytype",ToStr(copytype));
	doc.setOptValue("clustertype",ToStr(clustertype));
	doc.setOptValue("ocean",ToStr(ocean));
	doc.setOptValue("inflatePRA",ToStr(inflatePRA));
	doc.setOptValue("IsoHetY",ToStr(IsoHetY));
    }
};

/** This is the magic function that does the clustering for the outside world */
void bayes_label(std::vector<affx::GType> &call,
                 std::vector<double> &conf,
                 std::vector<double> &x,
                 std::vector<double> &y,
                 std::vector<int> &genohints,
                 std::vector<double> &SubsetInbred,
                 snp_param &sp);

void bayes_label(std::vector<affx::GType> &call, 
                 std::vector<double> &conf, 
                 std::vector<double> &x, 
                 std::vector<double> &y, 
                 std::vector<int> &genohints, 
                 std::vector<double> &SubsetInbred, 
                 snp_param &sp, 
                 std::vector<double> &ted,
                 std::vector<double> &teh,
                 std::vector<double> &tec, 
                 bool store_probabilites);

#endif /* _SNPLABEL_H_ */
