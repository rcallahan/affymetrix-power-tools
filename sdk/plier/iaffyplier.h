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
in probe behavior and handling error at the appropriately at low and 
high signal values.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/

/*
 * iaffyplier.h: iaffyplier interface class.
 *
 * NOTE: Users should use "caffyplier" as defined in "affyplier.h"
 *       This is an internal class.
 *
 */

#if !defined(PLIER_INTERFACE_HEADER)
#define PLIER_INTERFACE_HEADER

#ifdef _NO_PLIER_DLL_
#define plierdll
#else
#ifdef _MSVC
#define plierdll __declspec( dllexport )
#else
#ifdef _GCC
#define plierdll
#else
#ifdef __GNUC__
#define plierdll
#else
#define plierdll __declspec( dllimport )
#endif /* __GNUC__ */
#endif /* _GCC */
#endif /* _MSVC */
#endif /* _NO_PLIER_DLL_ */

#include "plier/affy_ptr.h"
//

#define FULL_MLE 0
#define SEA 1
#define MAX_IMPLNAME_LENGTH 255
#define MAX_IMPLVER_LENGTH 32
#define MAX_ERROR_LENGTH 1024
#define IDENTIFY_SIZE 2
#define IDENTIFY_MEDIAN 1
#define IDENTIFY_SUM 0

/*
 * PLIER Interfaces
 * An abstract C++ class
 */

class iaffyplier
{
protected:
	iaffyplier() { refcnt = 0; }
	virtual ~iaffyplier() {}

public:

	/*
	 * Increment reference count
	 */
	void addref() { ++refcnt; }

	/*
	 * Decrement reference count, if zero, delete the current object
	 */
	void release() 
	{ 
		if (refcnt > 0 && --refcnt == 0)
			delete this;
	}

	/* PLIER parameters - Initialization */
	virtual void setInitAugmentation(double) = 0;          /* avoid zero values in input data */
	virtual void getInitAugmentation(double*) = 0;         /* return byref the value of augmentation */
	virtual void setInitDefaultFeatureResponse(double) = 0;       /* default FeatureResponses if not supplied */
	virtual void getInitDefaultFeatureResponse(double*) = 0;      /* return byref the value of default FeatureResponse */
	virtual void setInitDefaultTargetResponse(double) = 0;  /* default TargetResponse if not supplied */
	virtual void getInitDefaultTargetResponse(double*) = 0; /* return byref the value of default TargetResponse*/

	/* PLIER parameters - SEA */
	virtual void setSeaAttenuation(float) = 0;             /* set up attenuation of background */
	virtual void getSeaAttenuation(float*) = 0;            /* return byref the value of SEA attenuation */

	/* PLIER parameters - SEA Optimization */
	virtual void setSeaOptConvergence(double) = 0;        /* change in log-value at which to stop */
	virtual void getSeaOptConvergence(double*) = 0;       /* return byref the value of SEA convergence */
	virtual void setSeaOptIteration(long) = 0;            /* max number of SEA iteration to avoid infinite loops in SEA */
	virtual void getSeaOptIteration(long*) = 0;           /* return byref the value of SEA iteration */

	/* PLIER parameters - PLIER */
	virtual void setPlierGmCutoff(float) = 0;              /* set up robustness */
	virtual void getPlierGmCutoff(float*) = 0;             /* return byref the value of GMCutoff */
	virtual void setPlierDifferentialFeaturePenalty(float) = 0;          /* Bayes penalty for peculiar features */
	virtual void getPlierDifferentialFeaturePenalty(float*) = 0;         /* return byref the value of feature penalty */
	virtual void setPlierDifferentialTargetPenalty(float) = 0;           /* Bayes penalty for really peculiar TargetResponses */
	virtual void getPlierDifferentialTargetPenalty(float*) = 0;          /* return byref the value of TargetResponse penalty */
	virtual void setPlierUseMMLikelihood(bool) = 0;                 /* use mm or background based likelihood */
	virtual void getPlierUseMMLikelihood(bool*) = 0;                /* return byref the value of usemm */
	virtual void setPlierUseInputModel(bool) = 0;              /* Use provided values as the initial model of Feature Responses */
	virtual void getPlierUseInputModel(bool*) = 0;             /* return byref the value of usemodel */
	virtual void setPlierFitFeatureResponse(bool) = 0;           /* Fit Feature Response dynamically or don't update from initial values */
	virtual void getPlierFitFeatureResponse(bool*) = 0;          /* return byref the value of fitFeatureResponse */

	/* PLIER parameters - PLIER Optimization */
	virtual void setPlierOptConvergence(double) = 0;      /* change in likelihood */
	virtual void getPlierOptConvergence(double*) = 0;     /* return byref the value of PLIER convergence */
	virtual void setPlierOptIteration(long) = 0;          /* max number of PLIER iteration to avoid infinite loops in PLIER */
	virtual void getPlierOptIteration(long*) = 0;         /* return byref the value of PLIER iteration */
	virtual void setPlierOptDropMax(double) = 0;          /* "logarithmic" barrier near zero */
	virtual void getPlierOptDropMax(double*) = 0;         /* return byref the value of dropmax */
	virtual void setPlierOptLambdaLimit(double) = 0;      /* minimum step multiplier in method */
	virtual void getPlierOptLambdaLimit(double*) = 0;     /* return byref the value of lambdalimit */
	virtual void setPlierOptOptimizationMethod(long) = 0;       /* optimization method. either FULL_OPTIMIZATION or SEA */
	virtual void getPlierOptOptimizationMethod(long*) = 0;      /* return byref the value of optimization */
	virtual void setPlierOptBalanceMethod(long) = 0;       /* identifiability method: MEDIAN or SUM */
	virtual void getPlierOptBalanceMethod(long*) = 0;      /* return byref the value of identifiability method */

	/* PLIER Inputs */
	virtual void setNumExp(long) = 0;                      /* number of input experiments */
	virtual void getNumExp(long*) = 0;                     /* return byref the number of input experiments */
	virtual void setNumFeature(long) = 0;                    /* number of input features in each experiment */
	virtual void getNumFeature(long*) = 0;                   /* return byref the number of input features in each experiment */
	virtual void setReplicate(long*) = 0;                   /* replicate groups */
	virtual void getReplicate(long**) = 0;                  /* return byref the pointer that set by using "set_replicate" method */
	virtual void setPM(double**) = 0;                       /* matrix of num_exp by num_probe stored by experiments (i.e. by row) */
	virtual void getPM(double***) = 0;                      /* return byref the pointer that set by using "set_pm" method */
	virtual void setMM(double**) = 0;                       /* matrix of mismatch or background values, corresponds to structure of pm */
	virtual void getMM(double***) = 0;                      /* return byref the pointer that set by using "set_mm" method */
	
	/*PLIER OPTIONAL INPUT*/
	virtual void setWT(double ** val) =0;					/*matrix of weights for input intensities*/
	virtual void getWT(double *** val)=0 ;					/*return by ref the pointer that set the matrix of weights*/

	/* PLIER Inputs/Outputs */
	virtual void setTargetResponse(double*) = 0;             /* array of length num_exp for estimated target responses */
	virtual void getTargetResponse(double**) = 0;            /* return byref the pointer that set by using "set_targetresponse" method */
	virtual void setFeatureResponse(double*) = 0;                  /* array of length num_probe for estimated feature responses*/
	virtual void getFeatureResponse(double**) = 0;                 /* return byref the pointer that set by using "set_featureresponse" method */

	/* PLIER Parameters (RESERVED FOR FUTURE) */
	virtual void setUseTargetResponse(bool) = 0;            /* 1 if initial estimates of TargetResponse are provided */
	virtual void getUseTargetResponse(bool*) = 0;           /* return byref the value of the use_TargetResponse flag */
	virtual void setUseBg(bool) = 0;                       /* 0 if the bg data is not initialized */
	virtual void getUseBg(bool*) = 0;                      /* return byref the value of usebg */
	virtual void setUseTargetResponseSd(bool) = 0;                  /* 1 if initial estimates of TargetResponses' standard deviation are provided */
	virtual void getUseTargetResponseSd(bool*) = 0;                 /* return byref the value of the use_TargetReponse_sd flag */
	virtual void setUseFeatureResponseSd(bool) = 0;              /* 1 if estimates of FeatureResponse standard deviation are provided */
	virtual void getUseFeatureResponseSd(bool*) = 0;             /* return byref the value of the use_FeatureResponse_sd flag */
	
	/* PLIER Parameters {added} */
	virtual void setFixPrecomputed(bool) = 0; /* 1 to recompute values as though supplying feature responses*/
	virtual void getFixPrecomputed(bool *) = 0; /* return byref precomputed flag */
	virtual void setNumericalTolerance(double ) = 0; /*converge to what numerical level*/
	virtual void getNumericalTolerance(double *) = 0; /*return byref numerical convergence level*/
	virtual void setSafetyZero(double ) = 0; /* ward off log(0) errors when reading in possible zero values*/
	virtual void getSafetyZero(double *) =0; /* get the tiny value by ref*/
	virtual void setFixFeatureEffect(bool b) = 0; /* Break when limit reached in plierapg::FitAdditiveModel() */
	virtual void getFixFeatureEffect(bool* b) = 0; /* Break when limit reached in plierapg::FitAdditiveModel() */


	

	/* PLIER Inputs (RESERVED FOR FUTURE) */
	virtual void setBg(double**) = 0;                       /* currently NULL, intended to allow for provision of background, corresponds to structure of pm */
	virtual void getBg(double***) = 0;                      /* return byref the pointer that set by using "set_bg" method */

	/* PLIER Inputs/Outputs (RESERVED FOR FUTURE) */
	virtual void setTargetResponseSd(double*) = 0;                   /* array of length num_exp for standard deviation of estimated intensities */
	virtual void getTargetResponseSd(double**) = 0;                  /* return byref the pointer that set by using "set_signal_sd" method */
	virtual void setFeatureResponseSd(double*) = 0;               /* array of length num_probe for standard deviation of estimated probe effects */
	virtual void getFeatureResponseSd(double**) = 0;              /* return byref the pointer that set by using "set_FeatureResponse_sd" method */

	/* PLIER Outputs (RESERVED FOR FUTURE) */
	virtual void setResiduals(double**) = 0;                /* Used to return residuals, structure corresponds to pm */
	virtual void getResiduals(double***) = 0;               /* return byref the pointer that set by using "set_residuals" method */
	virtual void setCallPvalue(double*) = 0;               /* array of length num_exp to return detection pvalues */
	virtual void getCallPvalue(double**) = 0;              /* return byref the pointer that set by using "set_call_pvalue" method */
	virtual void setLikelihood(double*) = 0;                /* The likelihood, or some similarly-interpretable quantity */
	virtual void getLikelihood(double**) = 0;               /* return byref the pointer that set by using "set_likelihood" method */

	/*
	 * Set to default PLIER parameters (usable for HG-U133A, for example)
	 * Input parameters :
	 *       buffer - [in] Reserved; must be 0
	 */
	virtual void setDefault(void* buffer=0) = 0;

	/*
	 * Interface to run the PLIER algorithm to obtain outputs
	 * Input parameter :
	 *       error_code - [out] Pointer to a long integer which receives
	 *                    the error code value, set to 0 if no error
	 */
	virtual void run(long* error_code) = 0;

	/*
	 * Get the plier object implementation name and version number
	 * Input parameters :
	 *       plier_impl_name - [out] Pointer to a buffer that receives
	 *                         a null-terminated string containing the implementation name.
	 *                         The buffer size should be large enough to contain
	 *                         MAX_IMPLNAME_LENGTH + 1 characters
	 *
	 *       plier_impl_ver  - [out] Pointer to a buffer that receives
	 *                         a null-terminated string containing the 
	 *                         The buffer size should be large enough to contain
	 *                         MAX_IMPLVER_LENGTH + 1 characters
	 */
	virtual void getImplInfo(char* plier_impl_name, char* plier_impl_ver) = 0;


private:
	int refcnt;  // reference count
};

/*
 * PLIER object creation function
 * Input parameters :
 *       plier_impl_name - [in] Pointer to a buffer that stores the plier implementation name.
 *                         It can be null and then default to PLIER v1 implementation
 *                         The buffer size should be at most MAX_IMPLNAME_LENGTH + 1 characters
 *
 *       plier           - [out] Pointer to a pointer of iaffyplier object to be returned
 */
plierdll void createPlierObject(const char* plier_impl_name, iaffyplier** plier);

/*
 * Get the error string
 * Input parameters :
 *       error_code - [in] Specifies the error code returned in run method
 *
 *       error      - [out] Pointer to a buffer that receives
 *                    a null-terminated string containing the error string.
 *                    The buffer size should be large enough to contain
 *                    MAX_ERROR_LENGTH + 1 characters
 */
plierdll void getPlierError(long error_code, char* error);

#endif /* defined(PLIER_INTERFACE_HEADER) */
