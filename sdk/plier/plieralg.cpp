////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

/*
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The PLIER (Probe Logarithmic Error Intensity Estimate) method produces
an improved signal by accounting for experimentally observed patterns 
in feature behavior and handling error at the appropriately at low and 
high signal values.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/

#ifndef __APPLE__
	#include <malloc.h>
#endif
#include "plier/plieralg.h"
//
#include "plier/affyheapsort.h"
#include "plier/error.h"
//
#include "stats/stats.h"
//
#include <cmath>
#include <cstring>
#include <memory>
#include <new>
#include <string.h>
//
//////////////////////////////////////////////////////////////////////

enum NormErrorCodes {
	NoError,
	MemoryError,
	NumErrors
};

//////////////////////////////////////////////////////////////////////

// static int g_ErrorCode=NoError; unused pk

//////////////////////////////////////////////////////////////////////
#define FPMIN	    0.00000001f
#define EPS		 1e-2
#define IQR_DIVIDER   1.34f

//////////////////////////////////////////////////////////////////////


inline
void TransferVector(double *A, double *B, long Length)
{
	for (long i=0; i<Length; i++)
		A[i] = B[i];
}

//////////////////////////////////////////////////////////////////////
inline
void ScrambleTransferVector(double *A, double *B, long *Twist, long Length)
{
	long i;
	for (i=0; i<Length; i++)
		A[i] = B[Twist[i]];
}

//////////////////////////////////////////////////////////////////////
inline
void InitializeVector(double *v, long n, double fValue)
{
	for (long i=0; i<n; i++)
		v[i] = fValue;
}
//////////////////////////////////////////////////////////////////////

double glog(double Value, double Offset)
{
	return( log((Value + sqrt (Value*Value+Offset))/2));
}

/////////////////////////////////////////////////////////////////////

double slog(double Value, double Start)
{
	return( log (Value + Start));
}


//////////////////////////////////////////////////////////////////////

void
LogVector(double *Vector, long Length)
{
	long i;
	for (i=0; i<Length; i++)
		Vector[i] = log(Vector[i]);
}
//////////////////////////////////////////////////////////////////////

void
sLogVector(double *Vector, long Length,double start)
{
	long i;
	for (i=0; i<Length; i++)
		Vector[i] = slog(Vector[i],start);
}

//////////////////////////////////////////////////////////////////////

void
ExpVector(double *Vector, long Length)
{
	long i;
	for (i=0; i<Length; i++)
		Vector[i] = exp(Vector[i]);
}

//////////////////////////////////////////////////////////////////////
inline
void ShrinkVector(double *Vector, double *ShiftVector, long Length, double dropmax)
{
	long i;
	for (i=0; i<Length; i++)
		if (Vector[i]+ShiftVector[i] < Vector[i]/dropmax)
			ShiftVector[i] = Vector[i]/dropmax - Vector[i];
}

//////////////////////////////////////////////////////////////////////

void StepVector(double *Vector, double *ShiftVector, long Length, double lambda)
{
	long i;
	for (i=0; i<Length; i++)
		Vector[i] += lambda*ShiftVector[i];
}

//////////////////////////////////////////////////////////////////////

inline
void AugmentData(double **Matrix, long NumExperiments, long Numfeatures, double Safety)
{
	long i,j;
	// ensure data is positive by adding small constant everywhere
	for (i=0; i<NumExperiments; i++)
		for (j=0; j<Numfeatures; j++)
			Matrix[i][j] += Safety; // small positive number
}

//////////////////////////////////////////////////////////////////////

void AugmentFloatData(float* data, long nSize, float safety)
{
	for (long i=0; i<nSize; i++)
		data[i] += safety;
}

//////////////////////////////////////////////////////////////////////

void AugmentDoubleData(double* data, long nSize, double safety)
{
	for (long i=0; i<nSize; i++)
		data[i] += safety;
}

//////////////////////////////////////////////////////////////////////

long FitAdditiveModel(plier_data *pData, double *TargetResponse, double *FeatureResponse,
			   double **TestDataTable, double **Weight, double ssqlimit, bool bNoFeatureResponse)
{
	long iter;
	long i,j, index, ti;
	//long medxval, medxvallo, medpval, medpvallo;
	double ssq;

	// do median polish by heapindex
	double median;
  // double *Val;
	// long *Rank;
	long ValSize;

	long nSize = pData->m_nAnalyses*pData->m_nFeatures;
  // this wont be freed at the end of the program.
  static double *Val = NULL;
  static int memSize = 0;

  if(nSize > memSize) {
    // free before reallocating
    delete[] Val;
    // now realloc scratch space
    Val = new double [nSize];
    if (Val == 0) {
      return NO_DATAMEM;
    }
    // remember the size.
    memSize = nSize;
  }

// 	Rank = new long [nSize]; // scratch space for rank
// 	if (Rank == 0)
// 	{
// 		delete[] Val;
// 		return NO_DATAMEM;
// 	}

	ssq = 10;
	for (iter=0; iter<pData->m_algParams->seaiteration && ssq>ssqlimit; iter++)
	{
		ssq = 0;
		for (i=0; i<pData->m_nAnalyses;)
		{
			index = 0;
			for (ti=i; ti<pData->m_nReplicates[i]; ti++)
			{
				// if replicates, would repeat this operation to load all values into Val
				for (j=0; j<pData->m_nFeatures; j++)
				{
					if (Weight[ti][j]>0)
					{
						Val[index] = TestDataTable[ti][j]-FeatureResponse[j]; // copy the values
						index++;
						//Val[index] = TestDataTable[ti*pData->m_nFeatures + j]-FeatureResponse[j]; // copy the values
					}
				}
			}
			// Valsize is number of things loaded
			ValSize = index;
			if (index>0)
			{
// 				HeapIndex(Val, 0, Rank, ValSize);

// 				// get median, subtract from column
// 				medpval = (ValSize-ValSize%2)/2;
// 				medpvallo = (ValSize-1-(ValSize-1)%2)/2;
// 				median = (Val[Rank[medpval]]+Val[Rank[medpvallo]])/2;
					 median = percentile_in_place(&Val[0], &Val[ValSize], 50);
			}
			else
				median = 0; // no data for this column, null model in output

			// if replicates, would set TargetResponse[i] for each replicate
			for (ti=i; ti<pData->m_nReplicates[i]; ti++)
			{
				ssq += (TargetResponse[ti]-median)*(TargetResponse[ti]-median);
				TargetResponse[ti] = median;
			}
			i = pData->m_nReplicates[i]; //jump forwards by an appropriate amount
		}
		// WBS 04/02/2008 - Calculate the ssq first so we can break if needed, then calculate the feature reposnse.
		for (j=0; j<pData->m_nFeatures && !bNoFeatureResponse; j++) // only do this if updating FeatureResponses
		{
			// no replicate features
			// but weights can zero out some features
			index = 0;
			for (i=0; i<pData->m_nAnalyses; i++)
			{
				if (Weight[i][j]>0)
				{
					Val[index] = TestDataTable[i][j]- TargetResponse[i];
					index++;
				}
			}
			ValSize = index;
			if (index>0)
			{
				median = percentile_in_place(&Val[0], &Val[ValSize], 50);
			}
			else {median = 0;} // no data, null model this column
			ssq += (FeatureResponse[j]-median)*(FeatureResponse[j]-median);
		}
		if ((pData->m_algParams->FixFeatureEffect) && (ssq <= ssqlimit)) {break;}
		for (j=0; j<pData->m_nFeatures && !bNoFeatureResponse; j++) // only do this if updating FeatureResponses
		{
			// no replicate features
			// but weights can zero out some features
			index = 0;
			for (i=0; i<pData->m_nAnalyses; i++)
			{
				if (Weight[i][j]>0)
				{
					Val[index] = TestDataTable[i][j]- TargetResponse[i];
					index++;
					//Val[i] = TestDataTable[i*pData->m_nFeatures+j]- TargetResponse[i];
				}
			}
			ValSize = index;
			if (index>0)
			{
// 				HeapIndex(Val, 0, Rank, ValSize);
// 				// get median, subtract from column
// 				medxval = (ValSize-ValSize%2)/2;
// 				medxvallo = (ValSize-1-(ValSize-1)%2)/2;
// 				median = (Val[Rank[medxval]]+Val[Rank[medxvallo]])/2;
						  median = percentile_in_place(&Val[0], &Val[ValSize], 50);
			}
			else {median = 0;} // no data, null model this column
			FeatureResponse[j] = median;
		}
	}

	   //	delete[] Val;
	   //	delete[] Rank;

	if (iter==pData->m_algParams->seaiteration && ssq>ssqlimit)
		return MAXIT_SEA_REACHED;
	else
		return NO_PLIER_ERROR;
}

//////////////////////////////////////////////////////////////////////

void BalanceFeatureResponse(double *TargetResponse, double *FeatureResponse,
					 long NumExperiments, long Numfeatures)
{
	double total;
	long i,j;

	// balance FeatureResponses for identifiability - crude - should build in at top?
	total = 0;
	for (j=0; j<Numfeatures; j++)
		total += FeatureResponse[j];
	total /= Numfeatures;

	for (i=0; i<NumExperiments; i++)
		TargetResponse[i] *= total;
	for (j=0; j<Numfeatures; j++)
		FeatureResponse[j] /=total;
}

//////////////////////////////////////////////////////////////////////

void BalanceMedianFeatureResponse(double *TargetResponse, double *FeatureResponse, long NumExperiments, long Numfeatures)
{
	// does this once at end if selected flag
  //	long *TmpRank;
	double med; 
	   //	long v, vlo,i,j;
	   int i,j;
// 	TmpRank = new long [Numfeatures];
// 	HeapIndex(FeatureResponse, 0,TmpRank, Numfeatures);
// 	v = (Numfeatures-Numfeatures%2)/2;
// 	vlo = (Numfeatures-1-(Numfeatures-1)%2)/2;
//	    median = (FeatureResponse[TmpRank[v]]+FeatureResponse[TmpRank[vlo]])/2;
	med = median(&FeatureResponse[0], &FeatureResponse[Numfeatures]);

	for (i=0; i<NumExperiments; i++)
		TargetResponse[i] *= med;
	for (j=0; j<Numfeatures; j++)
		FeatureResponse[j] /=med;
	   //	delete[] TmpRank;
}

//////////////////////

double ComputeSignalSize(double *TargetResponse, long NumExperiments)
{
	// compute reasonable values
	double sum;
	long i;
	sum = 0;
	for (i=0; i<NumExperiments; i++)
		sum += TargetResponse[i];
	sum /= NumExperiments;
	return(sum);
}

/////////////////////
void BalanceSignalSize(double *TargetResponse, double *FeatureResponse, long NumExperiments, long Numfeatures, double SignalSize)
{
	// does this once at end if selected flag
	double sum,med; 
	int i,j;

	sum = ComputeSignalSize(TargetResponse, NumExperiments);

	med = SignalSize/sum;
	// never increase the average signal based on robust measure
	// only decrease to maintain approximate balance
	if (med>1)
		med =1;

	for (i=0; i<NumExperiments; i++)
		TargetResponse[i] *= med;
	for (j=0; j<Numfeatures; j++)
		FeatureResponse[j] /=med;
}


//////////////////////////////////////////////////////////////////////
// Y = I-B
// U = glog(Y,H)
// SEA: :"additive model"
// SEA trick: H = 4*I*B*L (attenuated background)
// V = log(f) + log(t)
// R = U-V
// minimize abs(R) by Median Polish

long doSEA(plier_data *pData, double* TargetResponse,
		   double* FeatureResponse, double **U, double **Weight, bool bNoFeatureResponse)
{
	long i,j;
	double fFourTimesAttenuation = 4*pData->m_algParams->attenuation;
	// initialize data using safety factor added to the logarithm
	// when using input feature responses, could be zero due to rounding
	// better than nothing!
	//double safetyzeta = .000001;

	// set up the (glog) transformed intensity minus background matrix 
	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nFeatures; j++)
		{
			U[i][j] = glog(pData->m_fPM[i][j] - pData->m_fMM[i][j],fFourTimesAttenuation*pData->m_fPM[i][j]*pData->m_fMM[i][j]);
		}

	// additive model on log scale
	sLogVector(TargetResponse, pData->m_nAnalyses, pData->m_algParams->safetyZero);
	sLogVector(FeatureResponse, pData->m_nFeatures, pData->m_algParams->safetyZero);

	long lRes = FitAdditiveModel(pData, TargetResponse, FeatureResponse, U, Weight, pData->m_algParams->seaconvergence, bNoFeatureResponse);
	if (lRes != NO_PLIER_ERROR)
		return lRes;

	// return to regular data scale
	// note: explicitly ignoring the fact that safetyzero has been added
	ExpVector(TargetResponse, pData->m_nAnalyses);
	ExpVector(FeatureResponse, pData->m_nFeatures);

	// balance our data
	BalanceFeatureResponse(TargetResponse, FeatureResponse, pData->m_nAnalyses, pData->m_nFeatures);

	return NO_PLIER_ERROR;
}

//////////////////////////////////////////////////////////////////////
inline
double GenerateReferenceForX(double *X, long Length, double safetyZero)
{
	long i;
	// use slog in case reference contains zero values
	// due to rounding in input values
	// note: may add twice, but still fine
	double total=0;
	for (i=0; i<Length; i++)
		total += slog(X[i],safetyZero); // logarithm of X
	total /= Length;
	total = exp(total); // geometric mean
	return(total); // geometric mean
}

//////////////////////////////////////////////////////////////////////
double JustError(double TargetResponse,
		   double FeatureResponse,
		   double Hash, double PM, double MM, bool bUseGlog)
{
	double f, t, y, e;

	if (bUseGlog)
	{
	   f = FeatureResponse;
	 t = TargetResponse;
	 y = f * t;
	 //q = sqrt (y * y + Hash);
	 //r = (y + q) / (2 * PM);
	 //e = log (r);
	  // note glog(PM-MM, Hash) = log(PM)
	  e = glog(y,Hash)-log(PM);
	}
	else  // slog
	{
	 f = FeatureResponse;
	 t = TargetResponse;
	 y = f * t;
	 //r = (y + MM) / PM;
	 //e = log (r);
	  // note slog(PM-MM,MM) = log(PM)
	  e = slog(y,MM) - log(PM);
	}
	return(e);
}

//////////////////////////////////////////////////////////////////////

double JustFit(double TargetResponse, double FeatureResponse, double Hash, double PM, double MM, double z, bool bUseGlog)
{
	double e;

	e = JustError(TargetResponse, FeatureResponse, Hash, PM, MM, bUseGlog);
	return(e*e/(1+e*e/z));
}

//////////////////////////////////////////////////////////////////////
// drop derivatives for compute speed
// do not reformulate glog in case of rounding issues
inline
void JustPLIERMMLikelihood(double *Likelihood,
				  double TargetResponse,
				  double FeatureResponse,
				  double Hash,
				  double PM,
				  double MM,
				  double z, bool bUseGlog)
{
	double f,t,y,q,r,e,x,esq;

	if (bUseGlog)
	{
		f = FeatureResponse;
		t = TargetResponse;
		y = f*t;
		// unrolling the computation for efficiency, since I need q later 
		q = sqrt(y*y+Hash);
		r = (y+q)/(2*PM);
		// e = glog(y,Hash)-glog(PM-MM,Hash) = glog(y,Hash)-log(PM)
		// but nicer to take logs of things near 1 (when we fit)
		e = log(r);

		// Geman-McClure
		esq = e*e;
		x = 1+(esq)/z;
		// h = g^2/(1+g^2/z)
		*Likelihood = esq/x;
	}
	else  // use slog
	{
		f = FeatureResponse;
		t = TargetResponse;
		y = f*t;
		r = (y+MM)/PM;

		// e = slog(y,MM)-slog(PM-MM,MM) = slog(y,MM)-log(PM)
		e = log(r);

		// Geman-McClure
		esq = e*e;
		x = 1+esq/z;

		*Likelihood = esq/x;
	}
}


//////////////////////////////////////////////////////////////////////
inline
void PLIERMMLikelihood(double *Likelihood,
				  double *TDeriv,
				  double *FDeriv,
				  double *TGrad,
				  double *FGrad,
				  double TargetResponse,
				  double FeatureResponse,
				  double Hash,
				  double PM,
				  double MM,
				  double z, bool bUseGlog)
{
	double f,t,y,q,r,e,x,dh,dgT,dgF,ddh, xsq, esq;

	if (bUseGlog)
	{
		f = FeatureResponse;
		t = TargetResponse;
		y = f*t;
		// unrolling the computation for efficiency, since I need q later 
		q = sqrt(y*y+Hash);
		r = (y+q)/(2*PM);
		// e = glog(y,Hash)-glog(PM-MM,Hash) = glog(y,Hash)-log(PM)
		e = log(r);

		// Geman-McClure
		esq = e*e;
		x = 1+(esq)/z;
		xsq = x*x;
		// h = g^2/(1+g^2/z)
		*Likelihood = esq/x;
	// give me the likelihood I need
		dh = (2*e)/xsq;

		// dh = 2g dg
		//df = 2*e;

		// dg = 1/r dr
		// dr = (f/q)r
		dgT = f/q;
		dgF = t/q;
		*TDeriv = dh*dgT;
		*FDeriv = dh*dgF;

		// d^2h = 2g/(1+g^2/z)^2 d^2g + 2[3*(1+g^2/z)-4]/(1+g^2/z)^3 (dg)^2
		// drop negatives to get a useful gradient
		ddh = 2/xsq;

		*TGrad = ddh*(dgT*dgT);
		*FGrad = ddh*(dgF*dgF);
	}
	else  // use slog
	{
		f = FeatureResponse;
		t = TargetResponse;
		y = f*t;
		r = (y+MM)/PM;

		// e = slog(y,MM)-slog(PM-MM,MM) = slog(y,MM)-log(PM)
		e = log(r);

		// Geman-McClure
		esq = e*e;
		x = 1+esq/z;
		xsq = x*x;

		*Likelihood = esq/x;
		dh = (2*e)/xsq;

		// 1/r dr = PM/(y+MM)*f/PM = f/(y+MM)
		dgT = f/(y+MM);
		dgF = t/(y+MM);

		*TDeriv = dh*dgT;
		*FDeriv = dh*dgF;

		ddh = 2/xsq;
		*TGrad = ddh*(dgT*dgT);
		*FGrad = ddh*(dgF*dgF);
	}
}

//////////////////////////////////////////////////////////////////////
inline
void RoughnessPenaltyForX(double *LocalLikelihood, double *Deriv, double *Grad, double X, double Reference, double alpha, double alpha_times_2, double safetyZero)
{
	// some function of a TargetResponse against a reference level
	double green;
	double Y;
	Y = X + safetyZero;
	//double alpha_times_2=alpha*2;
	// reference as computed should already have safety factor
	green = Y/Reference;
	double loggreen = log(green);
	*LocalLikelihood = alpha*loggreen*loggreen;
	*Deriv = alpha_times_2*loggreen/Y;
	*Grad = alpha_times_2/(Y*Y);
}

//////////////////////////////////////////////////////////////////////
inline
double UpdateLikelihoodForRoughness(double *X, double *XDeriv, double *XGrad, long Length, double penalty, double Reference, double safetyZero)
{
	double  LocalLikelihood, Deriv, Grad, LogLikelihood;
	long i;
	LogLikelihood = 0;

	// do not compute reference here, to uncouple likelihood
	double penalty_times_2 = penalty*2;
	for (i=0; i<Length; i++)
	{
		RoughnessPenaltyForX(&LocalLikelihood,&Deriv, &Grad,X[i],Reference,penalty,penalty_times_2, safetyZero);
		LogLikelihood += LocalLikelihood;
		XDeriv[i] += Deriv;
		XGrad[i] += Grad;
	}
	return(LogLikelihood);
}

//////////////////////////////////////////////////////////////////////
inline
double ComputeGlobalLikelihood(plier_data* pData, double *TargetResponse, double *FeatureResponse,
							   double *NewtonTDeriv, double *NewtonFDeriv, double *NewtonTGrad, double *NewtonFGrad,
							   double **IHash, double **Weight, int NoDerivative)
{
	double LogLikelihood, LocalLikelihood;
	double TDeriv, FDeriv, TGrad, FGrad;
	long i,j;
	double ReferenceFeatureResponse, ReferenceTargetResponse;

	LogLikelihood = 0;
	memset(NewtonTDeriv, 0, pData->m_nAnalyses*sizeof(double));
	memset(NewtonTGrad, 0,  pData->m_nAnalyses*sizeof(double));
	memset(NewtonFDeriv, 0, pData->m_nFeatures*sizeof(double));
	memset(NewtonFGrad, 0, pData->m_nFeatures*sizeof(double));

	if (!NoDerivative)
	{
	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nFeatures; j++)
		{
			PLIERMMLikelihood(&LocalLikelihood, &TDeriv, &FDeriv, &TGrad, &FGrad,
				TargetResponse[i], FeatureResponse[j], IHash[i][j],
				pData->m_fPM[i][j], pData->m_fMM[i][j],
				pData->m_algParams->gmcutoff, pData->m_algParams->usemm);
			LogLikelihood += Weight[i][j]*LocalLikelihood;
			NewtonTDeriv[i] += Weight[i][j]*TDeriv;
			NewtonFDeriv[j] += Weight[i][j]*FDeriv;
			NewtonTGrad[i] += Weight[i][j]*TGrad;
			NewtonFGrad[j] += Weight[i][j]*FGrad;
		}
	}
	else
	{
	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nFeatures; j++)
		{
			JustPLIERMMLikelihood(&LocalLikelihood,
				TargetResponse[i], FeatureResponse[j], IHash[i][j],
				pData->m_fPM[i][j], pData->m_fMM[i][j],
				pData->m_algParams->gmcutoff, pData->m_algParams->usemm);
			LogLikelihood += Weight[i][j]*LocalLikelihood;
		}
	}

	ReferenceFeatureResponse = GenerateReferenceForX(FeatureResponse, pData->m_nFeatures, pData->m_algParams->safetyZero);
	ReferenceTargetResponse = GenerateReferenceForX(TargetResponse, pData->m_nAnalyses, pData->m_algParams->safetyZero);

	LogLikelihood += UpdateLikelihoodForRoughness(FeatureResponse, NewtonFDeriv, NewtonFGrad,
						pData->m_nFeatures, pData->m_algParams->differentialfeaturepenalty, ReferenceFeatureResponse, pData->m_algParams->safetyZero);
	LogLikelihood += UpdateLikelihoodForRoughness(TargetResponse, NewtonTDeriv, NewtonTGrad,
						pData->m_nAnalyses, pData->m_algParams->differentialtargetpenalty, ReferenceTargetResponse, pData->m_algParams->safetyZero);

	return(LogLikelihood);
}

//////////////////////////////////////////////////////////////////////

void Join_Replicates(long *Replicate, double *NewtonTDeriv, double *NewtonTGrad, long NumExperiments)
{
	long i, ti, NumRep;
	double TDeriv;
	double TGrad;

	for (i=0; i<NumExperiments;)
	{
		TDeriv = TGrad = 0;
		for (ti=i; ti<Replicate[i]; ti++)
		{
			TDeriv += NewtonTDeriv[ti];
			TGrad += NewtonTGrad[ti];
		}
		NumRep = Replicate[i]-i;
		TDeriv /= NumRep;
		TGrad /= NumRep;
		for (ti=i; ti<Replicate[i]; ti++)
		{
			NewtonTDeriv[ti] = TDeriv;
			NewtonTGrad[ti] = TGrad;
		}
		i = Replicate[i];
	}
}

//////////////////////////////////////////////////////////////////////

double ComputeExperimentLogLikelihood(plier_data *pData,
									  double *TargetResponse, double *FeatureResponse,
									  double *NewtonTDeriv, double *NewtonTGrad,
									  double **IHash, double **Weight, long WhichExp)
{
	double ReferenceTargetResponse;
	double LocalLikelihood;
	double LogLikelihood;
	long j, ti;


	ReferenceTargetResponse = GenerateReferenceForX(TargetResponse, pData->m_nAnalyses, pData->m_algParams->safetyZero);

	// for each experimental replicate
	LogLikelihood = 0;
	for (ti=WhichExp; ti<pData->m_nReplicates[WhichExp]; ti++)
	{
		for (j=0; j<pData->m_nFeatures; j++)
		{
			JustPLIERMMLikelihood(&LocalLikelihood, 
				TargetResponse[ti], FeatureResponse[j], IHash[ti][j],
				pData->m_fPM[ti][j], pData->m_fMM[ti][j],
				pData->m_algParams->gmcutoff, pData->m_algParams->usemm);
			LogLikelihood += Weight[ti][j]*LocalLikelihood;
		}
	}
	// roughness penalty over all experiments
	LogLikelihood += UpdateLikelihoodForRoughness(TargetResponse, NewtonTDeriv, NewtonTGrad,
						pData->m_nAnalyses, pData->m_algParams->differentialtargetpenalty, ReferenceTargetResponse, pData->m_algParams->safetyZero);
	return(LogLikelihood);
}

//////////////////////////////////////////////////////////////////////

double ComputefeatureLogLikelihood(plier_data* pData,
								 double *TargetResponse, double *FeatureResponse,
								 double *NewtonFDeriv, double *NewtonFGrad,
								 double **IHash, double **Weight, long Whichfeature)
{
	double ReferenceFeatureResponse;
	double  LocalLikelihood;
	double LogLikelihood;
	long i;

	ReferenceFeatureResponse = GenerateReferenceForX(FeatureResponse, pData->m_nFeatures, pData->m_algParams->safetyZero);
	// for each feature generate default values
	LogLikelihood = 0;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		JustPLIERMMLikelihood(&LocalLikelihood, 
				TargetResponse[i], FeatureResponse[Whichfeature], IHash[i][Whichfeature],
				pData->m_fPM[i][Whichfeature], pData->m_fMM[i][Whichfeature],
				pData->m_algParams->gmcutoff, pData->m_algParams->usemm);
		LogLikelihood += Weight[i][Whichfeature]*LocalLikelihood;
	}
	LogLikelihood += UpdateLikelihoodForRoughness(FeatureResponse, NewtonFDeriv, NewtonFGrad,
						pData->m_nFeatures, pData->m_algParams->differentialfeaturepenalty, ReferenceFeatureResponse, pData->m_algParams->safetyZero);
	return(LogLikelihood);
}

//////////////////////////////////////////////////////////////////////
long SearchExperimentGrid(plier_data* pData, double *TargetResponse, double *FeatureResponse,
					   double *NewtonTDeriv, double *NewtonTGrad, double *NewtonFDeriv, double *NewtonFGrad,
					   double **IHash, double **Weight,
					   double epsilon, bool bNoFeatureResponse)
{
	long i,j, ti, tti;
	double trialval, oldval;
	double oldLogLikelihood;
	double LogLikelihood;
	long converged;

	converged = 1;
	for (i=0; i<pData->m_nAnalyses;)
	{
		// for each replicate group compute best old likelihood
		oldLogLikelihood = ComputeExperimentLogLikelihood(pData, TargetResponse, FeatureResponse,
								NewtonTDeriv, NewtonTGrad, IHash, Weight, i);
		for (ti=i; ti<pData->m_nReplicates[i]; ti++)
		{
			for (j=0; j<pData->m_nFeatures; j++) // generate trial values for each replicate - denser attempts
			{
				if (FeatureResponse[j]>0)
					trialval = (pData->m_fPM[ti][j]-pData->m_fMM[ti][j])/(FeatureResponse[j]);
					//trialval = (pData->m_fPM[nIndex]-pData->m_fMM[nIndex])/(FeatureResponse[j]);
				else
					trialval = -1;
				if (trialval>0)
				{
					// set each replicate within the group containing experiment i to the same trial value
					oldval = TargetResponse[i];
					for (tti=i; tti<pData->m_nReplicates[i]; tti++)
						TargetResponse[tti] = trialval;
					LogLikelihood = ComputeExperimentLogLikelihood(pData, TargetResponse, FeatureResponse,
										NewtonTDeriv, NewtonTGrad, IHash, Weight, i);
					if (LogLikelihood<oldLogLikelihood)
					{
						converged = 0; // obviously found something better
						oldLogLikelihood = LogLikelihood; // the new champion
					}
					else // reset each replicate to the same old value
						for (tti=i; tti<pData->m_nReplicates[i]; tti++)
							TargetResponse[tti] = oldval;
				}
			}
		}
		i = pData->m_nReplicates[i]; // update to next value
	}
	return(converged);
}

///////////////////////////////////////////////////////////////////////////

long SearchFeatureGrid(plier_data* pData, double *TargetResponse, double *FeatureResponse,
					   double *NewtonTDeriv, double *NewtonTGrad, double *NewtonFDeriv, double *NewtonFGrad,
					   double **IHash, double **Weight,
					   double epsilon, bool bNoFeatureResponse)
{
	long i,j;
  //long ti, tti;
	double trialval, oldval;
	double oldLogLikelihood;
	double LogLikelihood;
	long converged;
	// now do the same for FeatureResponses
	converged = 1;
	for (j=0; j<pData->m_nFeatures && !bNoFeatureResponse; j++)
	{
		oldLogLikelihood = ComputefeatureLogLikelihood(pData, TargetResponse, FeatureResponse,
								NewtonFDeriv, NewtonFGrad, IHash, Weight, j);
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			// run thru experiments generating trial values
			if (TargetResponse[i]>0)
				trialval = (pData->m_fPM[i][j]-pData->m_fMM[i][j])/(TargetResponse[i]);
				//trialval = (pData->m_fPM[nIndex]-pData->m_fMM[nIndex])/(TargetResponse[i]);
			else
				trialval = -1;
			if (trialval>0)
			{
				oldval = FeatureResponse[j];
				FeatureResponse[j] = trialval;
				LogLikelihood = ComputefeatureLogLikelihood(pData, TargetResponse, FeatureResponse,
									NewtonFDeriv, NewtonFGrad, IHash, Weight, j);
				if (LogLikelihood<oldLogLikelihood)
				{
					converged = 0; // obviously found something better
					oldLogLikelihood = LogLikelihood; // the new champion
				}
				else
					FeatureResponse[j] = oldval;
			}
		}
	}
	return(converged);

}

//////////////////////////////////////////////////////////////////////

long SearchGridOptimum(plier_data* pData, double *TargetResponse, double *FeatureResponse,
					   double *NewtonTDeriv, double *NewtonTGrad, double *NewtonFDeriv, double *NewtonFGrad,
					   double **IHash, double **Weight,
					   double epsilon, bool bNoFeatureResponse)
{
	//long i,j, ti, tti;
	//double trialval, oldval;
	//double oldLogLikelihood;
	//double LogLikelihood;
	long converged,xconverged;

	converged = 1; // assume at optimum
	// search plausible grid for better possible optima
	// obvious grid placement is determined by TargetResponse values
	converged = SearchExperimentGrid(pData,TargetResponse,FeatureResponse,NewtonTDeriv,NewtonTGrad,NewtonFDeriv,NewtonFGrad,IHash,Weight,epsilon, bNoFeatureResponse);
	xconverged =converged;
	converged = SearchFeatureGrid(pData,TargetResponse,FeatureResponse,NewtonTDeriv,NewtonTGrad,NewtonFDeriv,NewtonFGrad,IHash,Weight,epsilon, bNoFeatureResponse);
	if (!xconverged || !converged)
		converged = 0;
	return(converged);
}

//////////////////////////////////////////////////////////////////////

long CorrectReplicatesSlow(long *NewOrder, long *Replicate, long NumExperiments)
{
	// but, what about replicates?
	long *TOrder = 0;
	long *TRep = 0;
	long i,j, ti, tj;

	TOrder = new long [NumExperiments];
	if (TOrder == 0)
		return NO_DATAMEM;

	TRep = new long [NumExperiments];
	if (TRep == 0)
	{
		delete[] TOrder;
		return NO_DATAMEM;
	}

	// do n^2 correction for replicates
	// which basically is slow, but compared with PLIER iterations is still quick
	// quickest writing of code to get something functional
	// use OldOrder as scratch space here
	j = 0; // start off with nothing
	for (i=0; i<NumExperiments; i++)
	{
		// start with the current least element
		if (NewOrder[i]>-1)
		{
			// add to new order
			tj =j;
			TOrder[j] = NewOrder[i];
			j++;
			// search rest of list
			for (ti=i+1; ti<NumExperiments; ti++)
			{
				if (NewOrder[ti]>-1)
				{
					// if part of the same replicate group, add here
					if (Replicate[NewOrder[ti]]==Replicate[NewOrder[i]])
					{
						TOrder[j] = NewOrder[ti]; // next replicate in order
						j++;
						NewOrder[ti] = -1;
					}
				}
			}
			NewOrder[i] = -1; // done this, turn off
			for (; tj<j; tj++)
				TRep[tj] = j; // set replicate number
		}
	}
	// replace NewOrder
	for (i=0; i<NumExperiments; i++)
	{
		NewOrder[i] = TOrder[i];
		Replicate[i] = TRep[i];
	}

	delete[] TOrder;
	delete[] TRep;
	return NO_PLIER_ERROR;
}

///////////////////////////////////////////////////////////////////////

long UnScrambleMatrix(plier_data *pData, long *Twist)
{
	long lRes = NO_PLIER_ERROR;
	double *TmpValue;
	long i,j;

	TmpValue = new double [pData->m_nAnalyses];
	if (TmpValue==0)
		return(NO_DATAMEM);

	// unscramble PM/MM values to original order
	for (j=0; j<pData->m_nFeatures; j++)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			TmpValue[i] = pData->m_fPM[Twist[i]][j];
		}
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			pData->m_fPM[i][j] = TmpValue[i];
		}
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			TmpValue[i] = pData->m_fMM[Twist[i]][j];
		}
		for (i=0; i<pData->m_nAnalyses; i++)
		{
			pData->m_fMM[i][j] = TmpValue[i];
		}
		if (pData->m_fWT!=0)
		{
			for (i=0; i<pData->m_nAnalyses; i++)
			{
				TmpValue[i] = pData->m_fWT[Twist[i]][j];
			}
			for (i=0; i<pData->m_nAnalyses; i++)
			{
				pData->m_fWT[i][j] = TmpValue[i];
			}
		}
	}

	delete[] TmpValue;
	return(lRes);
}

/////////////////////////////////////////////////////////////////////////

long UnScrambleReplicates(plier_data *pData, long *Twist)
{
	long lRes = NO_PLIER_ERROR;
	long *TmpValue;
	long i;

	TmpValue = new long [pData->m_nAnalyses];
	if (TmpValue==0)
		return(NO_DATAMEM);

	for (i=0; i<pData->m_nAnalyses; i++)
		TmpValue[i] = pData->m_nReplicates[Twist[i]];
	for (i=0; i<pData->m_nAnalyses; i++)
		pData->m_nReplicates[i] = TmpValue[i];

	delete[] TmpValue;
	return(lRes);
}



//////////////////////////////////////////////////////////////////////

long SortInputs(plier_data* pData, long *OldOrder)
{
	// sort the input data into a more useful format, keeping the old order in the appropriate position
	// Treats PMdata and MMdata columns as a single, large number
	// must handle ties, due to quantile normalization and possible duplicate values
	// i.e. A vs B is handled by A[0]<>B[0], ties broken by A[1]<>B[1], etc
	double **TmpMatrix = 0; // hold both PM and MM together
	long *NewOrder = 0;
	long i,j;
	long TmpLength;
	long lRes = NO_PLIER_ERROR;

	TmpLength = 2*pData->m_nFeatures;
	TmpMatrix = new double * [pData->m_nAnalyses];
	if (TmpMatrix == 0)
		return NO_DATAMEM;

	NewOrder = new long [pData->m_nAnalyses];
	if (NewOrder == 0)
	{
		delete[] TmpMatrix;
		return NO_DATAMEM;
	}

	long nFailedAt = -1;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		TmpMatrix[i] = new double [TmpLength];
		if (TmpMatrix[i] == 0)
		{
			nFailedAt = i;
			break;
		}
	}
	if (nFailedAt != -1)
	{
		for (i=0; i<nFailedAt; i++)
			delete[] TmpMatrix[i];
		delete[] TmpMatrix;
		delete[] NewOrder;
	}

	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nFeatures; j++)
		{
			TmpMatrix[i][j] = pData->m_fPM[i][j];
			TmpMatrix[i][j+pData->m_nFeatures] = pData->m_fMM[i][j];
		}

	for (i=0; i<pData->m_nAnalyses; i++)
		NewOrder[i] = i;

	// okay, have set up TmpMatrix
	HeapIndexMatrix(TmpMatrix, NewOrder, pData->m_nAnalyses, TmpLength);
	// NewOrder now contains a sort of the experiments

	// but the replicates need to be put in correct order!
	// i.e. adjacent so that the iterators work
	lRes = CorrectReplicatesSlow(NewOrder, pData->m_nReplicates, pData->m_nAnalyses);
	if (lRes == NO_PLIER_ERROR)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
			OldOrder[NewOrder[i]] = i; // inverse map

/*		for (i=0; i<pData->m_nAnalyses; i++)
			for (j=0; j<pData->m_nFeatures; j++)
			{
				pData->m_fPM[i][j] = TmpMatrix[NewOrder[i]][j];
				pData->m_fMM[i][j] = TmpMatrix[NewOrder[i]][j+pData->m_nFeatures];
			}*/
		UnScrambleMatrix(pData,NewOrder); // do the sorting
	}
	for (i=0; i<pData->m_nAnalyses; i++)
		delete[] TmpMatrix[i];
	delete[] TmpMatrix;
	delete[] NewOrder;
	return lRes;
}


////////////////////////////////////////////////////////////

void Zero_DataSheet(plier_datasheet *DataSheet)
{
	DataSheet->TargetResponse = 0;
	DataSheet->OldTargetResponse = 0;
	DataSheet->TDeriv = 0;
	DataSheet->TGrad = 0;
	DataSheet->TStep =0;
	DataSheet->FeatureResponse = 0;
	DataSheet->OldFeatureResponse = 0;
	DataSheet->FDeriv = 0;
	DataSheet->FGrad = 0;
	DataSheet->FStep = 0;
	DataSheet->IHash = 0;
	DataSheet->OldOrder = 0;

	DataSheet->IMinusB = 0;
	DataSheet->U = 0;
	DataSheet->Weight = 0;
}

//////////////////////////////////////////////////////////////////////////////////

long Allocate_DataSheet(plier_data *pData, plier_datasheet *DataSheet)
{
	long i, nFailedAt;

	DataSheet->TargetResponse = new double [pData->m_nAnalyses]; // experiments
	if (DataSheet->TargetResponse == 0)
		return NO_DATAMEM;

	DataSheet->OldTargetResponse = new double [pData->m_nAnalyses];
	if (DataSheet->OldTargetResponse == 0)
		return NO_DATAMEM;

	DataSheet->TDeriv = new double [pData->m_nAnalyses]; // direction
	if (DataSheet->TDeriv == 0)
		return NO_DATAMEM;

	DataSheet->TGrad = new double [pData->m_nAnalyses]; // how far to go
	if (DataSheet->TGrad == 0)
		return NO_DATAMEM;

	DataSheet->TStep = new double [pData->m_nAnalyses];
	if (DataSheet->TStep == 0)
		return NO_DATAMEM;

	DataSheet->FeatureResponse = new double [pData->m_nFeatures];
	if (DataSheet->FeatureResponse == 0)
		return NO_DATAMEM;

	DataSheet->OldFeatureResponse = new double [pData->m_nFeatures];
	if (DataSheet->OldFeatureResponse == 0)
		return NO_DATAMEM;

	DataSheet->FDeriv = new double [pData->m_nFeatures];
	if (DataSheet->FDeriv == 0)
		return NO_DATAMEM;

	DataSheet->FGrad = new double [pData->m_nFeatures];
	if (DataSheet->FGrad == 0)
		return NO_DATAMEM;

	DataSheet->FStep = new double [pData->m_nFeatures];
	if (DataSheet->FStep == 0)
		return NO_DATAMEM;


	DataSheet->IHash = new double * [pData->m_nAnalyses];
	if (DataSheet->IHash == 0)
		return NO_DATAMEM;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		DataSheet->IHash[i] = 0;
	}

	nFailedAt = -1;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		DataSheet->IHash[i] = new double [pData->m_nFeatures];
		if (DataSheet->IHash[i] == 0)
		{
			nFailedAt = i;
			break;
		}
	}
	if (nFailedAt != -1)
		return NO_DATAMEM;

	DataSheet->IMinusB = new double * [pData->m_nAnalyses];
	if (DataSheet->IMinusB == 0)
		return NO_DATAMEM;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		DataSheet->IMinusB[i] = 0;
	}

	nFailedAt = -1;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		DataSheet->IMinusB[i] = new double [pData->m_nFeatures];
		if (DataSheet->IMinusB[i] == 0)
		{
			nFailedAt = i;
			break;
		}
	}
	if (nFailedAt != -1)
		return NO_DATAMEM;

	DataSheet->U = new double * [pData->m_nAnalyses];
	if (DataSheet->U == 0)
		return NO_DATAMEM;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		DataSheet->U[i] = 0;
	}

	nFailedAt = -1;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		DataSheet->U[i] = new double [pData->m_nFeatures];
		if (DataSheet->U[i] == 0)
		{
			nFailedAt = i;
			break;
		}
	}
	if (nFailedAt != -1)
		return NO_DATAMEM;

	DataSheet->Weight = new double * [pData->m_nAnalyses];
	if (DataSheet->Weight == 0)
		return NO_DATAMEM;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		DataSheet->Weight[i] = 0;
	}

	nFailedAt = -1;
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		DataSheet->Weight[i] = new double [pData->m_nFeatures];
		if (DataSheet->Weight[i] == 0)
		{
			nFailedAt = i;
			break;
		}
	}
	if (nFailedAt != -1)
		return NO_DATAMEM;

	DataSheet->OldOrder = new long [pData->m_nAnalyses];
	if (DataSheet->OldOrder == 0)
		return NO_DATAMEM;

	return (0);
}

////////////////////////////////////////////////////////////////////////

int Delete_DataSheet(plier_data *pData, plier_datasheet *DataSheet)
{
	long i;

	if (DataSheet->TargetResponse!=0)
	{
		delete[] DataSheet->TargetResponse;
		DataSheet->TargetResponse = 0;
	}
	if (DataSheet->OldTargetResponse!=0)
		delete[] DataSheet->OldTargetResponse;
	if (DataSheet->TDeriv!=0)
		delete[] DataSheet->TDeriv;
	if (DataSheet->TGrad!=0)
		delete[] DataSheet->TGrad;
	if (DataSheet->TStep!=0)
		delete[] DataSheet->TStep;
	if (DataSheet->FeatureResponse!=0)
		delete[] DataSheet->FeatureResponse;
	if (DataSheet->OldFeatureResponse!=0)
		delete[] DataSheet->OldFeatureResponse;
	if (DataSheet->FDeriv!=0)
		delete[] DataSheet->FDeriv;
	if (DataSheet->FGrad!=0)
		delete[] DataSheet->FGrad;
	if (DataSheet->FStep!=0)
		delete[] DataSheet->FStep;
	if (DataSheet->IHash!=0)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
			if (DataSheet->IHash[i]!=0)
				delete[] DataSheet->IHash[i];
		delete[] DataSheet->IHash;
	}
	if (DataSheet->IMinusB!=0)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
			if (DataSheet->IMinusB[i]!=0)
				delete[] DataSheet->IMinusB[i];
		delete[] DataSheet->IMinusB;
	}
	if (DataSheet->U!=0)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
			if (DataSheet->U[i]!=0)
				delete[] DataSheet->U[i];
		delete[] DataSheet->U;
	}
	if (DataSheet->Weight!=0)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
			if (DataSheet->Weight[i]!=0)
				delete[] DataSheet->Weight[i];
		delete[] DataSheet->Weight;
	}
	if (DataSheet->OldOrder!=0)
		delete[] DataSheet->OldOrder;
	return(0);
}

//////////////////////////////////////////////////////////////////////////////////


void SetWeights(plier_data *pData, plier_datasheet &pSheet)
{
	long i,j;

	// initialize weights
	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nFeatures; j++)
			pSheet.Weight[i][j] = 1.0;
	// if weights have been passed in, use them
	if (pData->m_fWT!=0)
	{
		for (i=0; i<pData->m_nAnalyses; i++)
			for (j=0; j<pData->m_nFeatures; j++)
				pSheet.Weight[i][j] = pData->m_fWT[i][j];
	}
	// conceivably, there's other weighting that could happen
}

//////////////////////////////////////////////////////////////////////////////////

void SetPLIERHash(plier_data *pData, plier_datasheet &pSheet)
{
	long i,j;
	// set up hash for glog

	for (i=0; i<pData->m_nAnalyses; i++)
		for (j=0; j<pData->m_nFeatures; j++)
			pSheet.IHash[i][j] = 4*pData->m_fPM[i][j]*pData->m_fMM[i][j]; // prepped data for distance
}

///////////////////////////////////////////////////////////////////////////////////

void FindDescentDirection(plier_data *pData, plier_datasheet &pSheet)
{
	long i,j;

			// handle replicates:  Average Deriv and Grad for all replicates, substitute back
			// then everything steps the same amount
			// i.e. treat replicates as identical copies
// initialize step to zero
memset(pSheet.TStep, 0, pData->m_nAnalyses*sizeof(double));
memset(pSheet.FStep, 0, pData->m_nFeatures*sizeof(double)); 

			Join_Replicates(pData->m_nReplicates, pSheet.TDeriv, pSheet.TGrad, pData->m_nAnalyses);

			for (i=0; i<pData->m_nAnalyses; i++)
			{
				if (pSheet.TGrad[i]>0) // cut off bad values
					pSheet.TStep[i] = -pSheet.TDeriv[i]/(pSheet.TGrad[i]); // approximate inverse matrix by inverse diagonal
			}
			for (j=0; j<pData->m_nFeatures; j++)
			{
				if (pSheet.FGrad[j]>0) // cut off bad values
					pSheet.FStep[j] = -pSheet.FDeriv[j]/(pSheet.FGrad[j]);
			}
	// use logarithmic barrier to avoid actually hitting zero quickly
			ShrinkVector(pSheet.TargetResponse, pSheet.TStep, pData->m_nAnalyses, pData->m_algParams->dropmax);
			ShrinkVector(pSheet.FeatureResponse, pSheet.FStep, pData->m_nFeatures, pData->m_algParams->dropmax);
			if (!pData->m_algParams->fitFeatureResponse)
				memset(pSheet.FStep, 0, pData->m_nFeatures*sizeof(double)); // no change to feature responses based on data
}

////////////////////////////////////////////////////////////////////////////////////

void TryDescent(plier_data *pData, plier_datasheet &pSheet)
{
	double lambda;

	// back up current location
			pSheet.oldLogLikelihood = pSheet.LogLikelihood;
			TransferVector(pSheet.OldTargetResponse, pSheet.TargetResponse, pData->m_nAnalyses);
			TransferVector(pSheet.OldFeatureResponse, pSheet.FeatureResponse, pData->m_nFeatures);

			// try decreasing moves in current direction until succeed
			pSheet.LogLikelihood +=1; //not yet evaluated, but probably worse
			lambda = 2;
			while (pSheet.oldLogLikelihood<pSheet.LogLikelihood && lambda>pData->m_algParams->lambdalimit)
			{
				lambda = lambda/2;

				// restore from backup
				TransferVector(pSheet.TargetResponse, pSheet.OldTargetResponse, pData->m_nAnalyses);
				TransferVector(pSheet.FeatureResponse, pSheet.OldFeatureResponse, pData->m_nFeatures);

				// take a trial step
				StepVector(pSheet.TargetResponse, pSheet.TStep, pData->m_nAnalyses, lambda);
				StepVector(pSheet.FeatureResponse, pSheet.FStep, pData->m_nFeatures, lambda);
				
				// see if it did good!
				pSheet.LogLikelihood = ComputeGlobalLikelihood(pData, pSheet.TargetResponse, pSheet.FeatureResponse,
					pSheet.TDeriv, pSheet.FDeriv, pSheet.TGrad, pSheet.FGrad, pSheet.IHash, pSheet.Weight, 1);
			}

	// if didn't do any good, put back the way it was before
			if (pSheet.oldLogLikelihood<pSheet.LogLikelihood)
			{
				// do no harm
				TransferVector(pSheet.TargetResponse, pSheet.OldTargetResponse, pData->m_nAnalyses);
				TransferVector(pSheet.FeatureResponse, pSheet.OldFeatureResponse, pData->m_nFeatures);
				pSheet.LogLikelihood = pSheet.oldLogLikelihood;
			}
}
/////////////////////////////////////////////////////////////////////////////////////////////

double SumConvergence(double *X, double *Y, long length)
{
	double ssq = 0;
	long i;

	for (i=0; i<length; i++)
		ssq += (X[i]-Y[i])*(X[i]-Y[i]);
	return(ssq);
}

/////////////////////////////////////////////////////////////////////////////////////////////

long TestNumericalConvergence(plier_data *pData,plier_datasheet &pSheet)
{
	long converged;
	//double NumericalLimit = 0.1;
	double ssq;
	converged = 0;

	// test for numerical convergence between old/new values
	// only used when feature responses are fixed
	// because of stability issues
	ssq = SumConvergence(pSheet.TargetResponse,pSheet.OldTargetResponse,pData->m_nAnalyses);
	if (ssq< pData->m_algParams->NumericalConvergenceLimit)
	{
		converged = 1; 	
	}
    return converged;
}


///////////////////////////////////////////////////////////////////////////////////////////////

long TryGrid(plier_data *pData, plier_datasheet &pSheet)
{
	long converged;
	TransferVector(pSheet.OldTargetResponse, pSheet.TargetResponse, pData->m_nAnalyses);
	TransferVector(pSheet.OldFeatureResponse, pSheet.FeatureResponse, pData->m_nFeatures);
	pSheet.oldLogLikelihood = pSheet.LogLikelihood;
	// found at least a local optimum - try other possible minima
	converged = SearchGridOptimum(pData, pSheet.TargetResponse, pSheet.FeatureResponse,
	pSheet.TDeriv, pSheet.TGrad, pSheet.FDeriv, pSheet.FGrad, pSheet.IHash, pSheet.Weight,
	pData->m_algParams->plierconvergence, !pData->m_algParams->fitFeatureResponse);

	pSheet.LogLikelihood = ComputeGlobalLikelihood(pData, pSheet.TargetResponse, pSheet.FeatureResponse,
			pSheet.TDeriv, pSheet.FDeriv, pSheet.TGrad, pSheet.FGrad, pSheet.IHash, pSheet.Weight, 1);

	if (pSheet.oldLogLikelihood<pSheet.LogLikelihood+pData->m_algParams->plierconvergence)
	{
		converged = 1; // didn't do any good to find another point
		// test again, because roundoff error can build up if epsilon too small
		if (pSheet.oldLogLikelihood<pSheet.LogLikelihood)
		{
			TransferVector(pSheet.TargetResponse, pSheet.OldTargetResponse, pData->m_nAnalyses);
			TransferVector(pSheet.FeatureResponse, pSheet.OldFeatureResponse, pData->m_nFeatures);
		}
	}
	else
		converged =0; // it did work to find another point
	return(converged);
}

///////////////////////////////////////////////////////////////////////////////////////////////
long TestConvergence(plier_data *pData, plier_datasheet &pSheet)
{
	long converged;
	converged = 0;

	// are we within the zone where we don't care about likelihood?
	if (pSheet.oldLogLikelihood<pSheet.LogLikelihood+pData->m_algParams->plierconvergence)
	{
		converged++;
		// can we stop iterating here because numerical convergence reached?
		if (!pData->m_algParams->fitFeatureResponse && pData->m_algParams->fixPrecomputed)
			converged = TestNumericalConvergence(pData,pSheet);
		// if we've reached numerical convergence, or are ignoring it because feature responses still shifting
		// see if we're trapped in a local minimum anywhere
		//computationally expensive, so don't do it if don't need it.
		if (converged)
			converged=TryGrid(pData,pSheet);
	}
	else
		converged = 0;
	return(converged);
}

///////////////////////////////////////////////////////////////////////////////
// Y = I-B
// U = glog(Y,H)
// PLIER: "same transformation both sides"
// V = glog (f*t,H)
// PLIER trick: H = 4*I*B
// R = U-V
// minimize Geman-McClure(R) for f>=0, t>=0

long doPlier(plier_data *pData, plier_datasheet &pSheet, double &output)
{
	long count, icount, converged;

	long lRes = NO_PLIER_ERROR;

	// actually do the PLIER iterations

		count=0;
		icount = 0;
		converged = 0;
		pSheet.LogLikelihood=0;
		pSheet.oldLogLikelihood = 1;

		SetPLIERHash(pData, pSheet);

		while(count<pData->m_algParams->plieriteration && converged<1)
		{
			count++;

			pSheet.LogLikelihood = ComputeGlobalLikelihood(pData, pSheet.TargetResponse, pSheet.FeatureResponse,
				pSheet.TDeriv, pSheet.FDeriv, pSheet.TGrad, pSheet.FGrad, pSheet.IHash, pSheet.Weight, 0);

			FindDescentDirection(pData,pSheet);

			TryDescent(pData,pSheet);

			converged = TestConvergence(pData,pSheet);
			BalanceFeatureResponse(pSheet.TargetResponse, pSheet.FeatureResponse, pData->m_nAnalyses, pData->m_nFeatures);
		}

		if (count == pData->m_algParams->plieriteration && converged<1)
		{
			lRes = MAXIT_PLIER_REACHED;
		}
	return(lRes);
}

//////////////////////////////////////////////////////////////////////
// Y = I-B
// U = glog(Y,H)
// SEA: :"additive model"
// SEA trick: H = 4*I*B*L (attenuated background)
// V = log(f) + log(t)
// R = U-V
// minimize abs(R) by Median Polish
// Y = I-B
// U = glog(Y,H)
// PLIER: "same transformation both sides"
// V = glog (f*t,H)
// PLIER trick: H = 4*I*B
// R = U-V
// minimize Geman-McClure(R) for f>=0, t>=0

long NewtonPlier(plier_data* pData, double& output)
{
	long lRes;
	plier_datasheet pSheet;
	double SignalSize;
	
	output = -1;
	Zero_DataSheet(&pSheet);
	lRes = Allocate_DataSheet(pData, &pSheet);
	if (lRes==NO_DATAMEM)
	{
		Delete_DataSheet(pData, &pSheet);
		return(lRes);
	}

	// preprocess for stability
	SortInputs(pData, pSheet.OldOrder);

	// initialize inputs or retrieve preset feature responses
	InitializeVector(pSheet.TargetResponse, pData->m_nAnalyses, pData->m_algParams->defaultTargetResponse);
	if (pData->m_bUseModel == false)
		InitializeVector(pSheet.FeatureResponse, pData->m_nFeatures, pData->m_algParams->defaultFeatureResponse);
	else
		TransferVector(pSheet.FeatureResponse, pData->m_fFeatureResponse, pData->m_nFeatures);

	// set weights on the individual observations
	// perhaps based on feature quality, etc.
	SetWeights(pData,pSheet);

	AugmentData(pData->m_fPM, pData->m_nAnalyses, pData->m_nFeatures, pData->m_algParams->augmentation);
	AugmentData(pData->m_fMM, pData->m_nAnalyses, pData->m_nFeatures, pData->m_algParams->augmentation);

	// start with SEA for stability!
	lRes = doSEA(pData, pSheet.TargetResponse, pSheet.FeatureResponse, pSheet.U, pSheet.Weight, !pData->m_algParams->fitFeatureResponse);
	if (lRes != NO_PLIER_ERROR)
	{
		Delete_DataSheet(pData, &pSheet);
		return lRes;
	}
	SignalSize = ComputeSignalSize(pSheet.TargetResponse, pData->m_nAnalyses);

	if (pData->m_algParams->optimization == FULL_MLE)
	{
		lRes = doPlier(pData,pSheet, output);
		if (pData->m_algParams->fitFeatureResponse && pData->m_algParams->fixPrecomputed)
		{
		  // pretend we're reading in the data just
		  // as if we were not fitting feature responses
		  // for max-compatibility with one-cel-at-a time 
		  pData->m_algParams->fitFeatureResponse=0;
		  // do SEA to provide same initial value as if reading in
		  InitializeVector(pSheet.TargetResponse, pData->m_nAnalyses, pData->m_algParams->defaultTargetResponse);
		  // do not re-initialize feature response, of course!
		  lRes = doSEA(pData, pSheet.TargetResponse, pSheet.FeatureResponse, pSheet.U, pSheet.Weight, !pData->m_algParams->fitFeatureResponse);
		  	lRes = doPlier(pData,pSheet, output);
		  pData->m_algParams->fitFeatureResponse=1;

		}
	}

	if (pData->m_algParams->balancetype == IDENTIFY_MEDIAN)
		BalanceMedianFeatureResponse(pSheet.TargetResponse, pSheet.FeatureResponse,pData->m_nAnalyses,pData->m_nFeatures);
	if (pData->m_algParams->balancetype == IDENTIFY_SIZE)
	{
		// SignalSize comes from SEA
		// amplify low end
		BalanceSignalSize(pSheet.TargetResponse, pSheet.FeatureResponse,pData->m_nAnalyses,pData->m_nFeatures,SignalSize);
	}

	//TransferVector(pData->m_fTargetResponse, TargetResponse, pData->m_nAnalyses);
	ScrambleTransferVector(pData->m_fTargetResponse, pSheet.TargetResponse, pSheet.OldOrder, pData->m_nAnalyses);
	TransferVector(pData->m_fFeatureResponse, pSheet.FeatureResponse, pData->m_nFeatures);

	// clean up PM/MM matrix to be polite to further analysis
	// note that "data augmentation" is still in effect so everything will have augmentation added
	UnScrambleMatrix(pData, pSheet.OldOrder);
	UnScrambleReplicates(pData, pSheet.OldOrder);

	// track the overall fit
	output = pSheet.LogLikelihood/(pData->m_nAnalyses*pData->m_nFeatures);

	Delete_DataSheet(pData,&pSheet);
	return lRes;
}

//////////////////////////////////////////////////////////////////////////////////////////

long Compute_Signed_Residuals(plier_data* pData, long ResidualType)
{
	long i,j;

	// just compute signed error for each PM/MM pair currently
	// do not count TargetResponse penalty (goes with "TargetResponse" related values)
	// do not count FeatureResponse penalty (goes with "FeatureResponse" related values)
	// do not do Geman-McClure transformation 
	// user can compute weights from raw error values if they like
	// signed residuals more useful than unsigned
	for (i=0; i<pData->m_nAnalyses; i++)
	{
		for (j=0; j<pData->m_nFeatures; j++)
		{
			pData->m_fRS[i][j] = JustError(pData->m_fTargetResponse[i],
				pData->m_fFeatureResponse[j], 
				4*pData->m_fPM[i][j]*pData->m_fMM[i][j], 
				pData->m_fPM[i][j], 
				pData->m_fMM[i][j], 
				pData->m_algParams->usemm);

		}
	}
	return(NO_PLIER_ERROR);
}

