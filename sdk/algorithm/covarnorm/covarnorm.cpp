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

//
#include "algorithm/covarnorm/covarnorm.h"
//

/*
 * Just dump the vector out
 *
 * @param v - some vector we want to see dumped out
 */
void td(const std::vector<double> &v)
{
	cout.precision(12);
  cout.setf(ios_base::showpoint|ios_base::fixed|ios_base::showpos);
  //
	for(unsigned int i=0; i<v.size(); i++) {
		cout << v[i] << "\t";
  }
	cout << endl;
}

/*
 * comma delimited dump
 *
 * @param v - vector to be dumped
 * @param out - stream to dump it to
 */
void ctd(const std::vector<double> &v,std::ostream &out)
{
	for(unsigned int i=0; i<v.size(); i++)
		out << v[i] << ",";
	out<< "\t";
}


// scratch utility functions

/*
 * mean of a vector
 *
 * @param v - vector to be averaged
 */
double compute_mean(const std::vector<double> &v) {
	double sum;
	sum =0;
	for (unsigned int i=0; i<v.size(); i++)
		sum += v[i];
	sum /= v.size();
	return(sum);
}

/*
 * normal likelihood density without constants
 *
 * @param x - observed value
 * @param m - mean
 * @param s - standard deviation
 */
double dnorm(const double x,const double m,const double s)
{
	// don't care about constants
	return( exp(-(x-m)*(x-m)/(2*s*s))/s);
}

/*
 * multiply element-wise two vectors and sum
 *
 * @param a - first vector
 * @param b - second vector
 */
double vprod(const std::vector<double> &a,const std::vector<double> &b)
{
	double tmp=0.0;
	for (unsigned int i=0; i<a.size(); i++)
	{
		tmp += a[i]*b[i];
	}
	return(tmp);
}

/*
 * Fill vector with value x
 *
 * @param v - my vector
 * @param x - the value
 */
void FillVector(std::vector<double> &v,const double x)
{
	for (unsigned int i=0; i<v.size(); i++)
		v[i] = x;
}

/*
 * compute weights of data points for a genotype
 *
 * @param w - output weights
 * @param m - vector of estimated means for each snp for this genotype
 * @param sigma - within-snp variance
 * @param Contrast - actual observed values for each snp
 */
void Compute_Weights(std::vector<double> &w,
                     const std::vector<double> &m,
                     const double sigma,
                     const std::vector<double> &Contrast)
{
	// compute dnorm for everything
	for (unsigned int i=0; i<Contrast.size(); i++)
	{
		w[i] = dnorm(Contrast[i],m[i],sigma);
	}
}

/*
 *  Likelihood of the snp
 *
 * @param lik_per_snp - vector of likelihoods, one per snp
 * @param wAA - weight of AA genotype  model
 * @param wAB - weight of AB genotype model
 * @param wBB - weight of BB genotype model
 * @param p   - probability of a genotype overall
 */
void LikelihoodPerSNP(std::vector<double> &lik_per_snp,
                      const std::vector<double> &wAA,
                      const std::vector<double> &wAB,
                      const std::vector<double> &wBB,
                      const std::vector<double> &p)
{
	for (unsigned int i=0; i<wAA.size(); i++)
	{
		lik_per_snp[i] = p[2]*wAA[i]+p[1]*wAB[i]+p[0]*wBB[i];
	}
}

/*
 * Odds of this genotype per snp
 *
 * @param geno - relative prob of genotype
 * @param w - weight of snp
 * @param d - overall sum
 * @param p - probability of this genotype
 */
void GenotypePerSNP(std::vector<double> &geno,
                    const std::vector<double> &w,
                    const std::vector<double> &d,
                    const double p)
{
	for (unsigned int i=0; i<w.size(); i++)
		geno[i] = p*w[i]/d[i];
}

/*
 * computes the log-likelihood of the whole model
 *
 * @param lik_per_snp - likeihood per snp to be added up on log scale
 */
double ComputeLL(const std::vector<double> &lik_per_snp)
{
	double sum=0;
	for (unsigned int i=0; i<lik_per_snp.size(); i++)
		sum += log(lik_per_snp[i]);
	return(sum);
}

/*
 * Computes the within-genotype variation
 *
 * @param c - contrast value per snp
 * @param m - expected mean per snp
 * @param z - weight of this snp for computing sigma
 */
double WeightedSigma(const std::vector<double> &c,
                     const std::vector<double> &m,
                     const std::vector<double> &z)
{
	double sum,denom;
	sum=denom=0;
	for (unsigned int i=0; i<m.size(); i++)
	{
		sum += (c[i]-m[i])*(c[i]-m[i])*z[i];
		denom += z[i];
	}
	return(sqrt(sum/denom));
}

/*
 * In case variation too small, set a floor
 *
 * @param s - sigma values
 * @param m - minimal value
 */
void MinSigma(std::vector<double> &s,const double m)
{
	for (unsigned int i=0; i<s.size(); i++)
		if (s[i]<m)
			s[i] = m;
}

/*
 * Initialize the variables driving EM fit
 *
 * @param mu - mean for each genotype
 * @param sigma - variance within each genotype
 * @param probs - probabilities of each genotype
 * @pram Contrast - contrast for each snp
 */
void InitializeEMVars(std::vector<double> &mu,
                      std::vector<double> &sigma,
                      std::vector<double> &probs,
                      const std::vector<double> Contrast)
{
	vector<double> quantiles;
  // @todo fix this
	quantiles.resize(Contrast.size());
	for (unsigned int i=0; i<Contrast.size(); i++)
		quantiles[i] = Contrast[i];
	sort(quantiles.begin(),quantiles.end(),less<double>());
  //
	int qlevel;
	qlevel = quantiles.size()/6;
	mu.push_back(quantiles[qlevel]);
	mu.push_back(0);
	mu.push_back(quantiles[quantiles.size()-qlevel]);
	double tmpsigma;
	tmpsigma = mu[0]-quantiles[qlevel/2];
	tmpsigma += quantiles[quantiles.size()-qlevel/2]-mu[2];
	tmpsigma /=2;
	sigma.push_back(tmpsigma);
	sigma.push_back(tmpsigma/2);
	sigma.push_back(tmpsigma);
	probs.push_back(1/3.0);
	probs.push_back(1/3.0);
	probs.push_back(1/3.0);
}

/*
 * Does a "null" fit to the data, all estimates same for all snps
 *
 * @param Contrast - contrast values for each snp
 * @param Strength - strength values for each snp
 * @param zG - weight for each snp
 * @param Predictor - predictor function for this snp
 * @param Fitted - fitted values
 */
void Dummy_Fit(const std::vector<double> &Contrast,
               const std::vector<double> &Strength,
               const std::vector<double> &zG,
               std::vector<double> &Predictor,
               std::vector<double> &Fitted)
{
	// fit Contrast to covariates with weights zG
	// returning predicted values mG
	// currently uses no covariates at all and just takes means
	double sum, denom,tmp;
	sum=denom=0;

	for (unsigned int i=0; i<Contrast.size(); i++)
	{
		sum += Contrast[i]*zG[i];
		denom += zG[i];
	}
	tmp = sum/denom;
	Predictor.resize(4);
	Predictor[0]=tmp;
	Predictor[1]=Predictor[2]=Predictor[3]=0;
	FillVector(Fitted,tmp);
}

/*
 * Fits a weighted cubic regression on predictor(s)
 *
 * @param contrast - want to predict this value per snp
 * @param strength - covariate of choice
 * @param weights - weight of data points for this genotype
 * @param Predictor - output, prediction function coefficients
 * @param Predicted - output, predicted contrast per snp
 */
void
FitWeightedCubic(const std::vector<double> &contrast,
                 const std::vector<double> &strength,
                 const std::vector<double> &weights,
                 std::vector<double> &Predictor,
                 std::vector<double> &Predicted) {

  	// Singular value decomposition method
	unsigned int i;
	unsigned int nobs;
	unsigned int npred;
	npred = 3+1;
	nobs= contrast.size();

	// convert double into doubles to match newmat
  vector<Real> tmp_vec(nobs);
  Real* tmp_ptr = &tmp_vec[0];
	vector<Real> obs_vec(nobs);
	Real *obs_ptr = &obs_vec[0];
	vector<Real> weight_vec(nobs);

  Matrix covarMat(nobs,npred);
	ColumnVector observedVec(nobs);

	// fill in the data
	// modified by weights
	for (i=0; i<nobs; i++)
		weight_vec[i] = sqrt(weights[i]);

  	// load data - 1s into col 1 of matrix
	for (i=0; i<nobs; i++)
		tmp_vec[i] = weight_vec[i];
  	covarMat.Column(1) << tmp_ptr;
	for (i=0; i<nobs; i++)
		tmp_vec[i] *= strength[i];
  	covarMat.Column(2) << tmp_ptr;
	for (i=0; i<nobs; i++)
		tmp_vec[i] *= strength[i];
  	covarMat.Column(3) << tmp_ptr;
	for (i=0; i<nobs; i++)
		tmp_vec[i] *= strength[i];
  	covarMat.Column(4) << tmp_ptr;

  	for (i=0; i<nobs; i++)
		obs_vec[i] = contrast[i]*weight_vec[i];
  	observedVec << obs_ptr;

  	// do SVD
  	Matrix U, V;
  	DiagonalMatrix D;
    ColumnVector Fitted(nobs);
    ColumnVector A(npred);

  	SVD(covarMat,D,U,V);

  	Fitted = U.t() * observedVec;
  	A = V * ( D.i() * Fitted );

	// this predicts "0" for low weights
	// because of weighted regression
  	Fitted = U * Fitted;

	// this is the predictor
	Predictor.resize(npred);
	for (i=0; i<npred; i++)
		Predictor[i] = A.element(i);


  // export data back to doubles
	// and therefore this predicts "0" for low-weighted points
	// which is >not< the desired outcome!!!!
	// instead we need to predict all points at once
	// >unweighted< as output
	vector<double> Goofy;
	Predicted.resize(nobs);
  	for (i = 0; i < nobs; ++i) {
      Goofy.resize(npred);
      Goofy[0] = 1;
      Goofy[1] = strength[i];
      Goofy[2] = strength[i]*Goofy[1];
      Goofy[3] = strength[i]*Goofy[2];
      Predicted[i] = vprod(Goofy,Predictor);
  	}
}

/*
 * Sets a hard shell so that genotypes can't cross the center
 *
 * @param Contrast - observed contrast value
 * @param wAA - weight for the case of this point being AA
 * @param wBB - weight for the case of this point being BB
 */
void HardShell(const std::vector<double> &Contrast,
               std::vector<double> &wAA,
               std::vector<double> &wBB)
{
	for (unsigned int i=0; i<Contrast.size(); i++) {
    if (Contrast[i]<0)
      wAA[i] = 0;
    else
      wBB[i] = 0;
  }
}

/*
 *  Do EM repeatedly to fit normalization functions for each genotype
 *
 * @param Contrast - contrast values for each snp
 * @param Strength - strength values for each snp
 * @param max_iterations - don't go to infinity
 * @param log_lik_convergence - stopped making progress
 * @param NP - output prediction functions for each genotype
 */
void EMContrast(const std::vector<double> &Contrast,
                const std::vector<double> &Strength,
                const int max_iterations,
                const double log_lik_convergence,
                NormalizationPredictor &NP)
{
	// this is the only routine that needs newmat
	// so I can isolate it!
	vector<double> mu;
	vector<double> sigma;
	vector<double> probs;
	vector<double> wAA,wAB,wBB;
	vector<double> mAA,mAB,mBB;
	vector<double> zAA,zAB,zBB;
	vector<double> lik_per_snp;

	double last_log_lik,log_lik_diff,curr_log_lik;
	unsigned int n_iterations,nSNP;

	InitializeEMVars(mu,sigma,probs,Contrast);
	NP.TargetCenter.resize(3);
	NP.TargetCenter[0]=-.66;
	NP.TargetCenter[1]=0;
	NP.TargetCenter[2]= .66;
	NP.CenterCovariates.resize(3);
	NP.CenterCovariates[0] = 0;
	NP.CenterCovariates[1] = 0;
	NP.CenterCovariates[2] = 0;

	nSNP = Contrast.size();
	mAA.resize(nSNP);
	mAB.resize(nSNP);
	mBB.resize(nSNP);
	wAA.resize(nSNP);
	wAB.resize(nSNP);
	wBB.resize(nSNP);
	zAA.resize(nSNP);
	zAB.resize(nSNP);
	zBB.resize(nSNP);
	lik_per_snp.resize(nSNP);

	FillVector(mBB,mu[0]);
	FillVector(mAB,mu[1]);
	FillVector(mAA,mu[2]);
	Compute_Weights(wAA,mAA,sigma[2],Contrast);
	Compute_Weights(wAB,mAB,sigma[1],Contrast);
	Compute_Weights(wBB,mBB,sigma[0],Contrast);

	last_log_lik = -1000000;
	log_lik_diff = log_lik_convergence+1;
	n_iterations = 0;

	while (log_lik_diff > log_lik_convergence && n_iterations < max_iterations){
		n_iterations += 1;
		// E step
		// update probs for each group
		// relative genotype probability per SNP
		// relative likelihood steps
		LikelihoodPerSNP(lik_per_snp,wAA,wAB,wBB,probs);
		//cout<< "likpersnp\t";
		//td(lik_per_snp);
		GenotypePerSNP(zAA,wAA,lik_per_snp,probs[2]);
		GenotypePerSNP(zAB,wAB,lik_per_snp,probs[1]);
		GenotypePerSNP(zBB,wBB,lik_per_snp,probs[0]);
		/*td(zAA);
		td(zAB);
		td(zBB);*/
		probs[0] = compute_mean(zBB);
		probs[1] = compute_mean(zAB);
		probs[2] = compute_mean(zAA);
		//cout << "probs:\t";
		//td(probs);
		curr_log_lik = ComputeLL(lik_per_snp);
		//cout << "curr_lik:\t" << curr_log_lik << endl;
		log_lik_diff = curr_log_lik-last_log_lik;
		last_log_lik = curr_log_lik;
		// done with E step

		// now do the big M step
		// do magic: generate predictors for each group
		// how to do magic:  SVD of of design matrix outside of loop
		// do magic within loop
		// fitAA = predictor weighted by zAA
		// mAA = current predictions for this
		// fitAB = predictor weighted by zAB
		// mAB = current predictions for this
		// fitBB = predictor weighted by zBB
		// mBB = current predictions for this genotype by snp
		//Dummy_Fit(Contrast,mAA,zAA);
		Dummy_Fit(Contrast,Strength,zAB,NP.fitAB,mAB);
		//Dummy_Fit(Contrast,mBB,zBB);
		// try the real one
		FitWeightedCubic(Contrast,Strength,zAA,NP.fitAA,mAA);
		//FitWeightedCubic(Contrast,Strength,zAB,NP.fitAB,mAB);
		FitWeightedCubic(Contrast,Strength,zBB,NP.fitBB,mBB);

		/*td(zAA);
		td(mAA);
		td(zAB);
		td(mAB);
		td(zBB);
		td(mBB);*/

		sigma[0] = WeightedSigma(Contrast,mBB,zBB);
		sigma[1] = WeightedSigma(Contrast,mAB,zAB);
		sigma[2] = WeightedSigma(Contrast,mAA,zAA);
		MinSigma(sigma, .001);
		Compute_Weights(wAA,mAA,sigma[2],Contrast);
		Compute_Weights(wAB,mAB,sigma[1],Contrast);
		Compute_Weights(wBB,mBB,sigma[0],Contrast);
		HardShell(Contrast,wAA,wBB);
		/*td(wAA);
		td(wAB);
		td(wBB);
		cout << "sigma\t";
		td(sigma);*/

	}
	/*cout << "Fitted values" << endl;
	unsigned int dummy;
	for (dummy=0; dummy<Contrast.size(); dummy++)
	{
		cout << dummy << "\t";
		cout << Contrast[dummy] << "\t";
		cout << Strength[dummy] << "\t";
		cout << mAA[dummy] << "\t";
		cout << mAB[dummy] << "\t";
		cout << mBB[dummy] << "\t";
		cout << zAA[dummy] << "\t";
		cout << zAB[dummy] << "\t";
		cout << zBB[dummy] << "\t";
		cout << endl;
	}*/
	//td(Contrast);
	//td(mAA);
	//td(mAB);
	//td(mBB);
	// transfer useful data out of program
	// need: predictor vectors fAA, fAB, fBB
	// need: sigma
	NP.sigma.resize(3);
	NP.sigma[0]=sigma[0];
	NP.sigma[1]=sigma[1];
	NP.sigma[2]=sigma[2];
}

/*
 * Make covariates
 *
 * @param cov - covariate listing
 * @param Strength - the one initial covariate
 */
void NormalizationPredictor::MakeCov(std::vector<double> &cov,const double Strength)
{
	// refactor later to include vector<double> othervariables
	// also strength-target offset
	cov.resize(4);
	cov[0] = 1;
	cov[1] = Strength*cov[0];
	cov[2] = Strength*cov[1];
	cov[3] = Strength*cov[2];
}

/*
 * Make a prediction, given the covariate for a snp
 *
 * @param Contrast - input contrast
 * @param cov - covariates
 * @return new normalized contrast
 */
double NormalizationPredictor::MakePrediction(const double Contrast,
                                              const std::vector<double> &cov)
{
	// make a prediction given the covariate
	vector<double> predicted;
	double offset, AAd,ABd,BBd,sum;

	predicted.push_back(0);
	predicted.push_back(0);
	predicted.push_back(0);
	// get predicted center
	predicted[2] = vprod(cov,fitAA);
	predicted[1] = vprod(cov,fitAB);
	predicted[0] = vprod(cov,fitBB);

	// estimate odds of each genotype
	AAd = dnorm(Contrast,predicted[2],sigma[2]);
	ABd = dnorm(Contrast,predicted[1],sigma[1]);
	BBd = dnorm(Contrast,predicted[0],sigma[0]);
	if (Contrast>predicted[1])
		BBd = 0;
	if (Contrast<predicted[1])
		AAd = 0;
	sum = AAd+ABd+BBd;
	AAd /= sum;
	ABd /= sum;
	BBd /= sum;

	// weight needed offset by probability of each genotype
	offset = (TargetCenter[2]-predicted[2])*AAd;
	offset += (TargetCenter[1]-predicted[1])*ABd;
	offset += (TargetCenter[0]-predicted[0])*BBd;
	return(Contrast+offset); // estimated correction factor for contrast
}

/** copy off a previous set of functions */
void NormalizationPredictor::Copy(NormalizationPredictor &Source)
{
	fitAA.resize(Source.fitAA.size());
	fitAB.resize(Source.fitAB.size());
	fitBB.resize(Source.fitBB.size());
	sigma.resize(Source.sigma.size());
	CenterCovariates.resize(Source.CenterCovariates.size());
	TargetCenter.resize(Source.TargetCenter.size());
	copy(Source.fitAA.begin(),Source.fitAA.end(),fitAA.begin());
	copy(Source.fitAB.begin(),Source.fitAB.end(),fitAB.begin());
	copy(Source.fitBB.begin(),Source.fitBB.end(),fitBB.begin());
	copy(Source.sigma.begin(),Source.sigma.end(),sigma.begin());
	copy(Source.CenterCovariates.begin(),Source.CenterCovariates.end(),CenterCovariates.begin());
	copy(Source.TargetCenter.begin(),Source.TargetCenter.end(),TargetCenter.begin());
}

/** dump out some diagnostics */
void NormalizationPredictor::Dump(std::ostream &out)
{
	ctd(fitBB,out);
	ctd(fitAB,out);
	ctd(fitAA,out);
	ctd(sigma,out);
	ctd(CenterCovariates,out);
	ctd(TargetCenter,out);
	out << endl;
}

/** just the one covariate */
void NormalizationFrame::FillInCovariates(const std::vector<double> &S)
{
	Strength.resize(S.size());
	for (unsigned int i=0; i<S.size(); i++)
	{
		Strength[i] = S[i];
	}
}

/** just the dependent variable */
void NormalizationFrame::FillInContrast(const std::vector<double> &v) {
  	Contrast.resize(v.size());
  	for(unsigned int i = 0; i < v.size(); i++)
	{
    		Contrast[i] = v[i];
	}
}

/** wrapper for numeric function */
void NormalizationFrame::FitEM() {
	EMContrast(Contrast,Strength,50, .05,ThisFit);
}
