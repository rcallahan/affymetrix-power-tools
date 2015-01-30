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

#if !defined(_PLIERALGORITHM_HEADER_)
#define _PLIERALGORITHM_HEADER_

#include "plier/affyplier.h"
//

typedef struct
{
	long					m_nAnalyses;
	long					m_nFeatures;
	long*					m_nReplicates;
	double*				m_fTargetResponse;
	double*				m_fFeatureResponse;
	double**			m_fPM;
	double**			m_fMM;
	double**			m_fRS; // residuals
	double**			m_fWT; // weights for data points
	bool					m_bUseModel;
	plier_params*	m_algParams;
} plier_data;

////////////////////////////////////////////////////////////////////////

typedef struct {
	double *TargetResponse;
	double *OldTargetResponse;
	double *TDeriv;
	double *TGrad; // diagonal only
	double *TStep;
	double *FeatureResponse;
	double *OldFeatureResponse;
	double *FDeriv;
	double *FGrad; // diagonal only
	double *FStep; // step size

	double **IMinusB; // I-Background
	double **U; // I-Background transformed
	double **IHash; // log/glog magic value for each probe
	double **Weight; // weight for each probe pair
	long *OldOrder;

	double LogLikelihood;
	double oldLogLikelihood;

} plier_datasheet;

////////////////////////////////////////////////////////////////

long NewtonPlier(plier_data* pData, double& output);
long Compute_Signed_Residuals(plier_data* pData, long ResidualType);

////////////////////////////////////////////////////////////////////////

#endif // !defined(_PLIERALGORITHM_HEADER_)
