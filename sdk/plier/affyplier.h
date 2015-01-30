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

/*
 * affyplier.h: interface for the caffyplier class.
 */

#if !defined(PLIERCLASS_HEADER)
#define PLIERCLASS_HEADER

#include "plier/iaffyplier.h"
//

#define PLIER_IMPL_NAME "PLIER"
#define PLIER_IMPL_VER "1.1"

/*
 * structures for storing PLIER parameters
 */

typedef struct
{
	double augmentation;
	double defaultFeatureResponse;
	double defaultTargetResponse;
	double seaconvergence;
	double plierconvergence;
	double dropmax;
	double lambdalimit;
	float attenuation;
	float gmcutoff;
	float differentialfeaturepenalty;
	float differentialtargetpenalty;
	bool usemm;
	bool usemodel;
	bool fitFeatureResponse;
	long seaiteration;
	long plieriteration;
	long optimization;
	long balancetype;
	// new parameters
	// flag to use new tricks
	long fixPrecomputed;
	// when zero values are present take 'safe' started logs
	double safetyZero;
	// check convergence of TargetResponse to this limit
	double NumericalConvergenceLimit;
	bool FixFeatureEffect;
} plier_params;

/*
 * PLIER implementation class which inherits from iplier interface class
 */

class caffyplier : public iaffyplier
{
  // plier doesnt own any of this memory - dont try to free it.
	long num_exp;
	long num_feature;
	double** pm;
	double** mm;
	double** wt; // weights, if there
	double** residuals;
	long* replicate;
	double* TargetResponse;
	double* FeatureResponse;
	plier_params params;

public:
	caffyplier();
	virtual ~caffyplier();

public:

	/* PLIER parameters - Initialization */
	virtual void setInitAugmentation(double val) { params.augmentation=val; }
	virtual void getInitAugmentation(double* outval) { *outval=params.augmentation; }
	virtual void setInitDefaultFeatureResponse(double val) { params.defaultFeatureResponse=val; }
	virtual void getInitDefaultFeatureResponse(double* outval) { *outval=params.defaultFeatureResponse; }
	virtual void setInitDefaultTargetResponse(double val) {params.defaultTargetResponse=val; }
	virtual void getInitDefaultTargetResponse(double* outval) { *outval=params.defaultTargetResponse; }

	/* PLIER parameters - SEA */
	virtual void setSeaAttenuation(float val) { params.attenuation=val; }
	virtual void getSeaAttenuation(float* outval) { *outval=params.attenuation; }

	/* PLIER parameters - SEA Optimization */
	virtual void setSeaOptConvergence(double val) { params.seaconvergence=val; }
	virtual void getSeaOptConvergence(double* outval) { *outval=params.seaconvergence; }
	virtual void setSeaOptIteration(long val) { params.seaiteration=val; }
	virtual void getSeaOptIteration(long* outval) { *outval=params.seaiteration; }

	/* PLIER parameters - PLIER */
	virtual void setPlierGmCutoff(float val) { params.gmcutoff=val; }
	virtual void getPlierGmCutoff(float* outval) { *outval=params.gmcutoff; }
	virtual void setPlierDifferentialFeaturePenalty(float val) { params.differentialfeaturepenalty=val; }
	virtual void getPlierDifferentialFeaturePenalty(float* outval) { *outval=params.differentialfeaturepenalty; }
	virtual void setPlierDifferentialTargetPenalty(float val) { params.differentialtargetpenalty=val; }
	virtual void getPlierDifferentialTargetPenalty(float* outval) { *outval=params.differentialtargetpenalty; }
	virtual void setPlierUseMMLikelihood(bool val) { params.usemm=val; }
	virtual void getPlierUseMMLikelihood(bool* outval) { *outval=params.usemm; }
	virtual void setPlierUseInputModel(bool val) { params.usemodel=val; }
	virtual void getPlierUseInputModel(bool* outval) { *outval=params.usemodel; }
	virtual void setPlierFitFeatureResponse(bool val) { params.fitFeatureResponse=val; }
	virtual void getPlierFitFeatureResponse(bool* outval) { *outval=params.fitFeatureResponse; }
	
	/* PLIER parameters - PLIER Optimization */
	virtual void setPlierOptConvergence(double val) { params.plierconvergence=val; }
	virtual void getPlierOptConvergence(double* outval) { *outval=params.plierconvergence; }
	virtual void setPlierOptIteration(long val) { params.plieriteration=val; }
	virtual void getPlierOptIteration(long* outval) { *outval=params.plieriteration; }
	virtual void setPlierOptDropMax(double val) { params.dropmax=val; }
	virtual void getPlierOptDropMax(double* outval) { *outval=params.dropmax; }
	virtual void setPlierOptLambdaLimit(double val) { params.lambdalimit=val; }
	virtual void getPlierOptLambdaLimit(double* outval) { *outval=params.lambdalimit; }
	virtual void setPlierOptOptimizationMethod(long val) { params.optimization=val; }
	virtual void getPlierOptOptimizationMethod(long *outval) { *outval=params.optimization; }
	virtual void setPlierOptBalanceMethod(long val) { params.balancetype=val; }
	virtual void getPlierOptBalanceMethod(long *outval) { *outval=params.balancetype; }

	/* PLIER Inputs */
	virtual void setNumExp(long val) { num_exp=val; }
	virtual void getNumExp(long* val) { *val=num_exp; }
	virtual void setNumFeature(long val) { num_feature=val; }
	virtual void getNumFeature(long* val) { *val=num_feature; }
	virtual void setReplicate(long* val) { replicate=val; }
	virtual void getReplicate(long** val) { *val=replicate; }
	virtual void setPM(double** val) { pm=val; }
	virtual void getPM(double*** val) { *val=pm; }
	virtual void setMM(double** val) { mm=val; }
	virtual void getMM(double*** val) { *val=mm; }

	/*PLIER OPTIONAL INPUT*/
	virtual void setWT(double ** val){wt = val;}
	virtual void getWT(double *** val){*val = wt;}

	/* PLIER Inputs/Outputs */
	virtual void setTargetResponse(double* val) { TargetResponse=val; }
	virtual void getTargetResponse(double** val) { *val=TargetResponse; }
	virtual void setFeatureResponse(double* val) { FeatureResponse=val; }
	virtual void getFeatureResponse(double** val) { *val=FeatureResponse; }

	/* PLIER Parameters (RESERVED FOR FUTURE) */
	virtual void setUseTargetResponse(bool) {};
	virtual void getUseTargetResponse(bool*outval) {};
	virtual void setUseBg(bool) {};
	virtual void getUseBg(bool*outval) {};
	virtual void setUseTargetResponseSd(bool) {};
	virtual void getUseTargetResponseSd(bool*outval) {};
	virtual void setUseFeatureResponseSd(bool) {};
	virtual void getUseFeatureResponseSd(bool*outval) {};

	/* PLIER Parameters {added} */
	virtual void setFixPrecomputed(bool val){params.fixPrecomputed = val;};
	virtual void getFixPrecomputed(bool *val){*val = (params.fixPrecomputed!=0);};
	virtual void setNumericalTolerance(double inval){params.NumericalConvergenceLimit = inval;};
	virtual void getNumericalTolerance(double *outval){*outval = params.NumericalConvergenceLimit;};
	virtual void setSafetyZero(double val){params.safetyZero = val;};
	virtual void getSafetyZero(double *val){*val = params.safetyZero;};
	virtual void setFixFeatureEffect(bool b) {params.FixFeatureEffect = b;}
	virtual void getFixFeatureEffect(bool* b) {*b = params.FixFeatureEffect;}

	/* PLIER Inputs (RESERVED FOR FUTURE) */
	virtual void setBg(double**) {};
	virtual void getBg(double***) {};

	/* PLIER Inputs/Outputs (RESERVED FOR FUTURE) */
	virtual void setTargetResponseSd(double*) {};
	virtual void getTargetResponseSd(double**) {};
	virtual void setFeatureResponseSd(double*) {};
	virtual void getFeatureResponseSd(double**) {};

	// IMPLEMENTED FOR UTILITY
	virtual void setResiduals(double** val) {residuals=val;};
	virtual void getResiduals(double*** val) { *val=residuals;};
	/* PLIER Outputs (RESERVED FOR FUTURE) */
	virtual void setCallPvalue(double*) {};
	virtual void getCallPvalue(double**) {};
	virtual void setLikelihood(double*) {};
	virtual void getLikelihood(double**) {};

	/*
	 * Set to default PLIER parameters. This is called by the constructor.
     * in general you probably do not want to call any of these unless you
     * are explicitly trying to model a particular PLIER behavior.
     * 
	 * Input parameters :
	 *       buffer - [in] Reserved; must be 0
	 */
	virtual void setDefault(void* buffer=0);

    /*
     * Set the parameters to the AA set. This is the older PLIER behavior
     */
	virtual void setDefaultAA();

    /*
     * Set the parameters to the AB set. This is the default set. This
     * parameter set improves the stability of plier results when using
     * pre-computed feature responses. 
     */
	virtual void setDefaultAB();

	/*
	 * Interface to run the PLIER algorithm to obtain outputs
	 * Input parameter :
	 *       error_code - [out] Pointer to a long integer which receives
	 *                    the error code value, set to 0 if no error
	 */
	virtual void run(long* error_code);

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
	virtual void getImplInfo(char* plier_impl_name, char* plier_impl_ver);

protected:

	/* 
	 * if input replicate is NULL, set it to default replicate
	 * i.e. each experiment is in its own replicate group
	 */
	virtual void setDefaultReplicate(long size, long* replicate);

	/*
	 * validate the input parameters before running the algorithm
	 * return 0 if validation ok
	 */
	virtual long validateParams();

	/*
	 * validate the inputs before running the algorithm
	 * return 0 if validation ok
	 */
	virtual long validateInputs();

};

#endif /* !defined(PLIERCLASS_HEADER) */
