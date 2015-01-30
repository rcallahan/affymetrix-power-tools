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
an improved signal (estimate of target response) by accounting for experimentally observed patterns 
in probe behavior and handling error at the appropriately at low and 
high signal values.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/

/*
 * plier.cpp: implementation of the caffyplier class.
 */

#include "plier/affyplier.h"
//
#include "plier/error.h"
#include "plier/plieralg.h"
//
#include <cstring>
#include <memory.h>
#include <string.h>
#include <string>
//

#ifdef _MSC_VER
#pragma warning(disable: 4996) // don't show deprecated warnings.
#endif

/* PLIER 1.0 like parameters -- should give the same results as earlier plier versions */
#define DefaultParameterSettingAAaugmentation	0.1
#define DefaultParameterSettingAAdefaultFeatureResponse	1.0
#define DefaultParameterSettingAAdefaultTargetResponse	1.0
#define DefaultParameterSettingAAattenuation	0.005f
#define DefaultParameterSettingAAseaconvergence	0.000001
#define DefaultParameterSettingAAseaiteration	2000
#define DefaultParameterSettingAAgmcutoff	0.15f
#define DefaultParameterSettingAAdifferentialfeaturepenalty	0.001f
#define DefaultParameterSettingAAdifferentialtargetpenalty	0.000001f
#define DefaultParameterSettingAAusemm 1
#define DefaultParameterSettingAAusemodel 0
#define DefaultParameterSettingAAfitFeatureResponse 1
#define DefaultParameterSettingAAplierconvergence	0.000001
#define DefaultParameterSettingAAplieriteration	3000
#define DefaultParameterSettingAAdropmax	3.0
#define DefaultParameterSettingAAlambdalimit	0.01
#define DefaultParameterSettingAAoptimization 0
#define DefaultParameterSettingAAbalancetype 0
// These give the original PLIER 1.0 behavior
#define DefaultParameterSettingAAfixPrecomputed 0
#define DefaultParameterSettingAAsafetyZero 0
#define DefaultParameterSettingAANumericalConvergenceLimit 0
#define DefaultParameterSettingAAFixFeatureEffect false


/* PLIER 1.1: Changed 2006-06-26 -- stability changes related to using precomputed feature responses */
#define DefaultParameterSettingABaugmentation	0.1
#define DefaultParameterSettingABdefaultFeatureResponse	1.0
#define DefaultParameterSettingABdefaultTargetResponse	1.0
#define DefaultParameterSettingABattenuation	0.005f
#define DefaultParameterSettingABseaconvergence	0.000001
#define DefaultParameterSettingABseaiteration	2000
#define DefaultParameterSettingABgmcutoff	0.15f
#define DefaultParameterSettingABdifferentialfeaturepenalty	0.001f
#define DefaultParameterSettingABdifferentialtargetpenalty	0.000001f
#define DefaultParameterSettingABusemm 1
#define DefaultParameterSettingABusemodel 0
#define DefaultParameterSettingABfitFeatureResponse 1
#define DefaultParameterSettingABplierconvergence	0.000001
#define DefaultParameterSettingABplieriteration	3000
#define DefaultParameterSettingABdropmax	3.0
#define DefaultParameterSettingABlambdalimit	0.01
#define DefaultParameterSettingABoptimization 0
#define DefaultParameterSettingABbalancetype 0
// These give the PLIER 1.1 behavior with improved stability when using precomputed feature responses
#define DefaultParameterSettingABfixPrecomputed 1
#define DefaultParameterSettingABsafetyZero 0.000001
#define DefaultParameterSettingABNumericalConvergenceLimit .1
#define DefaultParameterSettingABFixFeatureEffect false


/***************************************************************/
/** caffyplier class implementation                           **/
/***************************************************************/

/*
 * Constructor
 */
caffyplier::caffyplier()
{
	memset(&params, 0, sizeof(plier_params));
	num_exp = 0;
	num_feature = 0;
	pm = NULL;
	mm = NULL;
	wt = NULL;
	residuals = NULL;
	replicate = NULL;
	TargetResponse = NULL;
	FeatureResponse = NULL;
	setDefault();
}

/*
 * Destructor
 */
caffyplier::~caffyplier()
{
}

/*
 * set to default parameters
 */
void caffyplier::setDefault(void* buffer)
{
    setDefaultAB();
}

// Set the AA parameter set -- PLIER 1.0 like behavior
void caffyplier::setDefaultAA() {
	params.augmentation =               DefaultParameterSettingAAaugmentation;
	params.defaultFeatureResponse =     DefaultParameterSettingAAdefaultFeatureResponse;
	params.defaultTargetResponse =      DefaultParameterSettingAAdefaultTargetResponse;
	params.attenuation =                DefaultParameterSettingAAattenuation;
	params.seaconvergence =             DefaultParameterSettingAAseaconvergence;
	params.seaiteration =               DefaultParameterSettingAAseaiteration;
	params.gmcutoff =                   DefaultParameterSettingAAgmcutoff;
	params.differentialfeaturepenalty = DefaultParameterSettingAAdifferentialfeaturepenalty;
	params.differentialtargetpenalty =  DefaultParameterSettingAAdifferentialtargetpenalty;
	params.usemm =                      DefaultParameterSettingAAusemm;
	params.usemodel =                   DefaultParameterSettingAAusemodel;
	params.fitFeatureResponse =         DefaultParameterSettingAAfitFeatureResponse;
	params.plierconvergence =           DefaultParameterSettingAAplierconvergence;
	params.plieriteration =             DefaultParameterSettingAAplieriteration;
	params.dropmax =                    DefaultParameterSettingAAdropmax;
	params.lambdalimit =                DefaultParameterSettingAAlambdalimit;
	params.optimization =               DefaultParameterSettingAAoptimization;
	params.balancetype =                DefaultParameterSettingAAbalancetype;
	params.fixPrecomputed =             DefaultParameterSettingAAfixPrecomputed;
	params.safetyZero =                 DefaultParameterSettingAAsafetyZero;
	params.NumericalConvergenceLimit =  DefaultParameterSettingAANumericalConvergenceLimit;
	params.FixFeatureEffect =			DefaultParameterSettingAAFixFeatureEffect;
}

// Set the default parameter set -- PLIER 1.1
void caffyplier::setDefaultAB() {
	params.augmentation =               DefaultParameterSettingABaugmentation;
	params.defaultFeatureResponse =     DefaultParameterSettingABdefaultFeatureResponse;
	params.defaultTargetResponse =      DefaultParameterSettingABdefaultTargetResponse;
	params.attenuation =                DefaultParameterSettingABattenuation;
	params.seaconvergence =             DefaultParameterSettingABseaconvergence;
	params.seaiteration =               DefaultParameterSettingABseaiteration;
	params.gmcutoff =                   DefaultParameterSettingABgmcutoff;
	params.differentialfeaturepenalty = DefaultParameterSettingABdifferentialfeaturepenalty;
	params.differentialtargetpenalty =  DefaultParameterSettingABdifferentialtargetpenalty;
	params.usemm =                      DefaultParameterSettingABusemm;
	params.usemodel =                   DefaultParameterSettingABusemodel;
	params.fitFeatureResponse =         DefaultParameterSettingABfitFeatureResponse;
	params.plierconvergence =           DefaultParameterSettingABplierconvergence;
	params.plieriteration =             DefaultParameterSettingABplieriteration;
	params.dropmax =                    DefaultParameterSettingABdropmax;
	params.lambdalimit =                DefaultParameterSettingABlambdalimit;
	params.optimization =               DefaultParameterSettingABoptimization;
	params.balancetype =                DefaultParameterSettingABbalancetype;
	params.fixPrecomputed =             DefaultParameterSettingABfixPrecomputed;
	params.safetyZero =                 DefaultParameterSettingABsafetyZero;
	params.NumericalConvergenceLimit =  DefaultParameterSettingABNumericalConvergenceLimit;
	params.FixFeatureEffect =			DefaultParameterSettingABFixFeatureEffect;
}

/*
 * validate the input parameters before calling PLIER algorithm
 */
long caffyplier::validateParams()
{
	if (params.augmentation < 0)
		return INVALID_AUGMENTATION;

	/* gmcutoff should not be 0, otherwise, division by zero will occur */
	if (params.gmcutoff == 0)
		return INVALID_GMCUTOFF;

	/* dropmax should be > 0, otherwise, division by zero will occur */
	if (params.dropmax <= 0)
		return INVALID_DROPMAX;

	/* concpenalty should be non-zero */
	if (params.differentialtargetpenalty == 0)
		return INVALID_DIFFERENTIALTARGETPENALTY;

	/* probepenalty should be non-zero */
	if (params.differentialfeaturepenalty == 0)
		return INVALID_DIFFERENTIALFEATUREPENALTY;

	if (params.optimization != FULL_MLE && params.optimization != SEA)
		return INVALID_OPTIMIZATION;

	/* seaiteration should be > 0 */
	if (params.seaiteration <= 0)
		return INVALID_SEAITERATION;

	/* plieriteration should be > 0	*/
	if (params.optimization == FULL_MLE && params.plieriteration <= 0)
		return INVALID_PLIERITERATION;

	if (params.safetyZero<0)
		return INVALID_SAFETYZERO;
	
	if (params.NumericalConvergenceLimit<0)
		return INVALID_NUMERICALCONVERGENCE;
	
	return NO_PLIER_ERROR;
}

/*
 * validate the inputs before calling PLIER algorithm
 */
long caffyplier::validateInputs()
{
	if (num_exp <= 0)
		return INVALID_NUM_EXP;

	if (num_feature <= 0)
		return INVALID_NUM_FEATURE;

	if (pm == 0)
		return INVALID_PM;

	if (mm == 0)
		return INVALID_MM;

	// if weights null, assume all weights 1.0

	// if residuals null, assume no return required.

	if (TargetResponse == 0)
		return INVALID_TARGETRESPONSE;

	if (FeatureResponse == 0)
		return INVALID_FEATURERESPONSE;

	return NO_PLIER_ERROR;
}

/*
 * run the plier algorithm
 */

void caffyplier::run(long* error_code)
{
	*error_code = validateParams();
	if (*error_code != NO_PLIER_ERROR)
		return;

	*error_code = validateInputs();
	if (*error_code != NO_PLIER_ERROR)
		return;

	bool create_def_replicate = false;
	if (!replicate)
	{
		replicate = new long[num_exp];
		if (replicate == 0)
		{
			*error_code = NO_DATAMEM;
			return;
		}
		setDefaultReplicate(num_exp, replicate);
		create_def_replicate = true;
	}

	plier_data inputs;
	inputs.m_nAnalyses = num_exp;
	inputs.m_nFeatures = num_feature;
	inputs.m_nReplicates = replicate;
	inputs.m_fPM = pm;
	inputs.m_fMM = mm;
	inputs.m_fRS = residuals;
	inputs.m_fWT = wt; // set weights if present
	inputs.m_fFeatureResponse = FeatureResponse;
	inputs.m_fTargetResponse = TargetResponse;
	inputs.m_bUseModel = params.usemodel;
	inputs.m_algParams = &params;

	double output;
	*error_code = NewtonPlier(&inputs, output);

	// now we have affinities/TargetResponses fit, if we have residuals, fit them
	if (
		(
			*error_code == NO_PLIER_ERROR || 
			*error_code == MAXIT_SEA_REACHED ||
			*error_code == MAXIT_PLIER_REACHED
		) && residuals!=0
	)


	{
		// construct signed residuals & return 
		// assumptions - ignore TargetResponse/FeatureResponse penalties
		// want just the fit of the data
		// TargetResponses, affinities, >augmented< data, as returned from NewtonPlier
		// get raw signed residuals [no Geman-McClure] for each PM/MM value & return
		*error_code = Compute_Signed_Residuals(&inputs, 0);
	}

	if (create_def_replicate && replicate)
	{
		delete[] replicate;
		replicate = 0;
	}
}

/* 
 * set to default replicate groups
 */
void caffyplier::setDefaultReplicate(long size, long* replicate)
{
	if (replicate)
		for (long i=0; i<size; i++)
			replicate[i] = i;
}

/*
 * get the implementation name and version
 */
void caffyplier::getImplInfo(char* plier_impl_name, char* plier_impl_ver)
{
	try
	{
		strcpy(plier_impl_name, PLIER_IMPL_NAME);
		strcpy(plier_impl_ver, PLIER_IMPL_VER);
	}
	catch(...)
	{
	}
}
