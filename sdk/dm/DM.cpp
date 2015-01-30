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
 * @file   DM.cpp
 * @author Xiaojun Di
 * @date   March 8 2006
 * 
 * DRAFT: wraps the DM algorithmics for genotyping under SDK
 * 
 */

//
#include "dm/DM.h"
//
#include "dm/DmPTable.h" // renamed to avoid clash with mas5-stat/src/pTable.h
//
#include <string.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////
// min max macros
#define MIN(x, y)	((x)<(y))? (x) : (y)
#define MAX(x, y)	((x)>(y))? (x) : (y)


  /**
  * @brief get the probe quartet given the allele type and match type
  * @param allele - allele type: A or B
  * @param matchtype - the match type: perfect or mis
  * @return the probe index
  */
int CQuartet::GetIndex(const int& allele, const int& matchtype)
  {
    if(allele == A_Allele && matchtype == PM) return 0;
    else if(allele == A_Allele && matchtype == MM) return 1;
    else if(allele == B_Allele && matchtype == PM) return 2;
    else if(allele == B_Allele && matchtype == MM) return 3;
    else return -1; // should never be here
  }

  /**
  * @brief default constructor
  */
CQuartet::CQuartet()
  {
    memset(pixels, 1, sizeof(int)*NUMCELLS);
    memset(intensity, 0, sizeof(float)*NUMCELLS);
    memset(variance, 0, sizeof(float)*NUMCELLS);
  }

/////////////////////////////////////////////////////////////////////////////////////////////
// function for convenience
/**
* @brief converting function between model and state
* @param model - model index
* @return the state the model is mapped to
*/
int GetState(const int& model) 
{
	if(model == 0) return 1;
	else if(model == 1) return 5;
	else if(model == 2) return 4;
	else return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// function for convenience
/**
* @brief A simple utility function
* 
* @param mu - observed mean
* @param var - variance
* @param estmu - estimated mean
* @return the omega value
*/
double Omega(const double& mu, const double& var, const double& estmu)
{	
	return (var + (mu - estmu)*(mu - estmu));
}

/////////////////////////////////////////////////////////////////////////////////////////////
// calculate the estimated means
/**
* @brief calculate the mean for a given model
* @param quartet - the probe quartet
* @param model - the current model
* @param mu - the mean to be returned
* @return void
*/
void CalcMu(const CQuartet &quartet, const int& model, double mu[])
{
	double bgNum = 0, bgDen = 0, fgNum = 0, fgDen = 0, bgmu = 0, fgmu = 0;
	int state = GetState(model), idx = 0;
	for(idx = 0; idx < NUMCELLS; idx++) 
	{	
		if(!(state & (1 << idx)))
		{ 
			bgNum += quartet.pixels[idx]*quartet.intensity[idx];
			bgDen += quartet.pixels[idx];
		}
		else 
		{ 
			fgNum += quartet.pixels[idx] * quartet.intensity[idx];
			fgDen += quartet.pixels[idx];
		}
	}

	// Calculate estimated intensities
	if(bgDen > 0) bgmu = bgNum/bgDen;
	if(fgDen>0) fgmu = fgNum/fgDen;

	for(idx = 0; idx < NUMCELLS; idx++) 
	{	
		if(!(state & (1 << idx)))
		{
			mu[idx] = bgmu;
		}
		else 
		{
			mu[idx] = fgmu;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
// calculate the estimated variances
/**
* @brief calculate the sigma for a given model with given means
* @param quartet - the probe quartet
* @param model - the current model
* @param mu - the mean to this model
* @param sigma - the sigmas to be returned
* @return void
*/
void CalcSigma(const CQuartet &quartet, double mu[], double sigma[], const int& model)
{
	double omega = 0, bgNum = 0, bgDen = 0, fgNum = 0, fgDen = 0, bgParm,fgParm = 0;
	int state = GetState(model), ib = 0;
	for(ib = 0; ib < NUMCELLS; ib ++) 
	{	
		if(!(state & (1 << ib)))
		{ 
			omega = Omega(quartet.intensity[ib],quartet.variance[ib],mu[ib]);
			bgNum += quartet.pixels[ib]*omega;
			bgDen += quartet.pixels[ib];
		}
		else 
		{ 
			omega = Omega(quartet.intensity[ib],quartet.variance[ib],mu[ib]);
			fgNum += quartet.pixels[ib] * omega;
			fgDen += quartet.pixels[ib];
		}
	}

	// Calculate estimated variances	
	bgParm = bgNum/bgDen;
	if(fgDen > 0) fgParm = fgNum/fgDen;
	for(ib = 0; ib < NUMCELLS; ib++) 
	{	
		if(!(state & ( 1 << ib)))
		{
			sigma[ib] = bgParm;
		}
		else 
		{
			sigma[ib] = fgParm;
		}
	}
}

/** 
 * Round the double and return it.
 * @param d - The double to round.
 * @param iDecimalPlaces - The number of decimal places to round to.
 * @return - rounded double value.
 */ 
double round(double val,int places)
{
  double v_int,v_frac;
  double place_10,v_int10,v_frac10;

  // split into parts.
  v_frac=modf(val,&v_int);

  // create our shift
  place_10=pow(10.0,places);
  // shift the frac part up and truncate
  v_frac10=modf(v_frac*place_10,&v_int10);
  // std-c++ round says that it rounds away from zero.
  // apply the rounding of the last digit.
  if (v_frac10>=0.5) {
    v_int10=v_int10+1.0;
  }
  else if (v_frac10<=-0.5) {
    v_int10=v_int10-1.0;
  }

  // shift fract back down.
  v_frac=v_int10/place_10;

  // put back together.
  return (v_int+v_frac);
}


/////////////////////////////////////////////////////////////////////////////////////////////
// calculate the loglikelihood for a quartet
/**
* @brief calculate loglikelihood of a probe quartet given model mean and sigma
* @param quartet - the probe quartet
* @param mu - the mean to the model
* @param sigma - the mean to the model
* @return the loglikelihood
*/
double CalcLogLikelihood(const CQuartet& quartet, double mu [], double sigma[])
{	
	double omega;
	double sum = 0;

	for( int ib=0; ib < NUMCELLS; ib++) 
	{	
		if(sigma[ib]>0)
		{ 
			omega = Omega(quartet.intensity[ib],quartet.variance[ib],mu[ib]);
			omega = omega/sigma[ib];
			double f = quartet.pixels[ib] * (omega + log10(sigma[ib]) + LOG2PI);
			sum += round(f, 10);
		}
	}
	sum *= -.5;

	return sum;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// calculate the loglikelihood for a quartet
/**
* @brief calculate loglikelihood of a probe quartet for a given model
* @param quartet - the probe quartet
* @param model - the model in consideration
* @return the loglikelihood
*/
float CalcLogLikelihood(const CQuartet& quartet, const int& model)
{
	double mu[NUMCELLS], sigma[NUMCELLS], fVal;
	
	// calculate estimated mean
	CalcMu(quartet, model, mu);

	// calculate estimated variances
	CalcSigma(quartet, mu, sigma, model);

	// calculate loglikelihood
	fVal = CalcLogLikelihood(quartet, mu, sigma);

	return (float)fVal;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// check large background, true if background is greater than foreground
/**
* @brief check whether backgroud is higher than foreground or not for a given model
*
* @param quartet - the probe quartet to be checked
* @param model - the model to be tested
* @return true if yes, false otherwise
*/
bool CheckLargeBackground(const CQuartet &quartet, const int& model)
{
	float fg = 0, bg = 0;

	// this is the original one
	if(model == 0) // AA
	{
		fg = quartet.intensity[0];
		bg = MAX(MAX(quartet.intensity[1],quartet.intensity[2]),quartet.intensity[3]);
	}
	else if(model == 2) // BB
	{
		fg = quartet.intensity[2];
		bg = MAX(MAX(quartet.intensity[0],quartet.intensity[1]),quartet.intensity[3]);
	}
	else if(model == 1) // AB
	{
		fg = MIN(quartet.intensity[0],quartet.intensity[2]);
		bg = MAX(quartet.intensity[1],quartet.intensity[3]);
	}
	else return false;

/*
	// this would be a better one in term of ballancing het and hom
	// may add a flag to switch the two later
	if(model == 0) // AA
	{
		fg = quartet.intensity[0];
		bg = (quartet.intensity[1] + quartet.intensity[2] + quartet.intensity[3])/3;
	}
	else if(model == 2) // BB
	{
		fg = quartet.intensity[2];
		bg = (quartet.intensity[0] + quartet.intensity[1] + quartet.intensity[3])/3;
	}
	else if(model == 1) // AB
	{
		fg = 0.5*(quartet.intensity[0] + quartet.intensity[2]);
		bg = 0.5*(quartet.intensity[1] + quartet.intensity[3]);
	}
	else return false;
*/

	return (bg >= fg);
}

/////////////////////////////////////////////////////////////////////////////////////////////
// get p-value through wilcoxon signed rank test
/**
* Apply Wilcoxon signed rank test on the loglikelihood difference vector to compute a p-value
*
* @param x - the loglikelihood difference vector
* @return the p-value
*/
float GetPValue(const vector<float>& x)
{
  //	int len = x.size();

	// 1. Ignore all zero differences.
	vector <float> newdiff = x;
	int n = 0, i = 0;
	for (i = 0; i < x.size(); ++i) 
	{
		if (x[i] != 0.0) 
		{
			newdiff[n] = x[i];
			n++;
		}
	}
	newdiff.resize(n);

	// 2.  Assign integer ranks to the differences.
	vector <float> ranks(n);
	for (i=0; i<n; ++i) 
	{
		ranks[i] = i + 1.0f;
	}

	float p_value = 1.0f;
	if (n == 0) // No non-zero differences.  Output 0.5 as the one-sided p-value and detection is absent.
		return 0.5f;
	else 
	{
		// 3. Convert differences to absolute values and sort in ascending order.
		vector <pair<float,int> > absdiff(n);
		for (i = 0; i < n; ++i) 
		{
			absdiff[i].first = fabs(newdiff[i]);
			absdiff[i].second = i;
		}
		sort(absdiff.begin(), absdiff.end());

		// 4. If there are ties among absolute differences, all differences in a tie
		//    group are assigned to a rank equal to the average of the integer ranks.
		int nTies = 0;	
		for (i = 0; i < n - 1; ++i) 
		{
			if (absdiff[i].first == absdiff[i+1].first) 
			{
				nTies ++;
				break;
			}
		}

		vector <float> invr(n);
		for (i = 0; i < n; ++i) 
		{
			invr[absdiff[i].second] = ranks[i];
		}
		
		float w = 0;
		for (i = 0; i < n; ++i) 
		{
			if (newdiff[i] > 0)
				w += invr[i];
		}

		if (nTies == 0)
		{
			int iCode = 0;
			for (i=0; i < n; ++i)
			{
				if (newdiff[i] > 0)
				{
					iCode += 1 << ((int) invr[i] - 1);
				}
			}
			p_value = dmFGetPValue(n-1, iCode);
		}
		else
		{
			// p-value = sum(u(Sj > S) + 0.5 * u(Sj = S)) / 2 to power n
			// where j = 1 to 2 to power n.
			// and u(Sj > S) = 1 if Sj > S and 0 otherwise.
			// and u(Sj = S) = 1 if Sj = S and 0 otherwise.
			int twoToN = 1 << n;
			vector<int> mask(n);
			for (i = 0; i < n; ++i) 
			{
				mask[i] = 1 << i;
			}
			vector <int> posRanks(twoToN);
			for (i = 0; i < twoToN; ++i) 
			{
				float sum = 0;
				for (int j = 0; j < n; ++j) 
				{
					if (i & mask[j])
						sum += ranks[j];
				}
				posRanks[i] = (int)sum; /// @todo Why are we summing
			}
			float tail = 0;
			for (i = 0; i < twoToN; ++i) 
			{
				if (posRanks[i] > w) 
				{
					tail ++;
				}
				else if (posRanks[i] == w) 
				{
					tail += 0.5;
				}
			}
			p_value = tail / (float) twoToN;
		}
	}

	return p_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// compare two pairs based only on p-values
struct PairLess : public binary_function<pair<float,int>,pair<float,int>,bool>
{
	bool operator () (const pair<float,int>& x, const pair<float,int>& y) const
	{
		return (x.first < y.first);
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////
// Get DM genotype call with a p-value as confidence
// Convention for models: NC-> -1, AA->0, AB->1, BB->2
/**
* the main entry function of DM algorithm, which starts with a vector of probe quartets and some parameters
* then run through the whole DM flow, return the genotype along with a p-value associated with
* 
* @param vquartets - a vector of probe quartets
* @param results - a result container to hold the genotype and p-values and return
* @param hetMult - the heterozygote multiplier to be multiplied on AB model
* @return true if no error, false otherwise
*/
bool CallDM(const vector<CQuartet> &vquartets, pair<float,int>& results, const float& hetMult)
{
	int nQ = vquartets.size();

	// no less than 1, no more than 14
	if(nQ <= 0 || nQ > 14) return false;

	// the loglikelihood matrix for all quartets
	vector<vector<float> > LLMatrix(nQ);
	int iq = 0, model = -1;
	for(iq = 0; iq < nQ; iq ++)
	{
		LLMatrix[iq].resize(NUMCELLS);
	}

	// there are multiple quartets, we make call base on all of them
	// calculate probe quartet level loglikelihoods for all available quartets
	for(iq = 0; iq < nQ; iq ++)
	{
		LLMatrix[iq][0] = CalcLogLikelihood(vquartets[iq], -1);
		for(model = 1; model < 4; model ++)
		{
			// check large background for models other than NC
			if(CheckLargeBackground(vquartets[iq], model-1))
			{
				LLMatrix[iq][model] = LLMatrix[iq][0];
			}
			else
			{
				LLMatrix[iq][model] = CalcLogLikelihood(vquartets[iq], model - 1);
			}
		}

		// apply the het multiplier to AB model by adding the desired constant to AB loglikelihood
		LLMatrix[iq][2] += hetMult;
	}

	// the loglikelihood ratio vector, prepare for SNP level aggregation
	vector<float> LLRatioVector(nQ);
	vector<pair<float,int> > PvalModelVector(4);

	for(model = -1; model < 3; model ++)
	{
		for(iq = 0; iq < nQ; iq ++)
		{
			float fRatio = FLOATINFINITY;
			for(int ib = 0; ib < NUMCELLS; ib ++)
			{
				if(ib == model + 1) continue;
				if(fRatio > (LLMatrix[iq][model+1] - LLMatrix[iq][ib]))
				{
					fRatio = LLMatrix[iq][model+1] - LLMatrix[iq][ib];
				}
			}
			LLRatioVector[iq] = fRatio;
		}
		PvalModelVector[model+1].first = GetPValue(LLRatioVector);
		PvalModelVector[model+1].second = model;
	}

	// to make it backward compatible, swap AB and BB to keep the original order for breaking ties 
	pair<float,int> temp = PvalModelVector[2];
	PvalModelVector[2] = PvalModelVector[3];
	PvalModelVector[3] = temp;

	// SNP level aggregation, the best is the one with the smallest p-value
	sort(PvalModelVector.begin(), PvalModelVector.end(), PairLess());

	// the results, p-value first, genotype call second
	results.first = round(PvalModelVector[0].first, 6);
	results.second = PvalModelVector[0].second;

	return true;
}
