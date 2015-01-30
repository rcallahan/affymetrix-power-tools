////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////


#include "stdafx.h"
#include "FusionGenotypeProbeSetResultsCOM.h"
#include "COMStringUtils.h"

// CFusionGenotypeProbeSetResultsCOM

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionGenotypeProbeSetResults
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::get_AlleleCall(CHAR* pVal)
{
	*pVal = results.GetAlleleCall();
	return S_OK;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::get_Confidence(float* pVal)
{
	*pVal = results.GetConfidence();
	return S_OK;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::get_RAS1(float* pVal)
{
	*pVal = results.GetRAS1();
	return S_OK;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::get_RAS2(float* pVal)
{
	*pVal = results.GetRAS2();
	return S_OK;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::get_PValueAA(float* pVal)
{
	*pVal = results.GetPValueAA();
	return S_OK;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::get_PValueAB(float* pVal)
{
	*pVal = results.GetPValueAB();
	return S_OK;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::get_PValueBB(float* pVal)
{
	*pVal = results.GetPValueBB();
	return S_OK;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::get_PValueNoCall(float* pVal)
{
	*pVal = results.GetPValueNoCall();
	return S_OK;
}

STDMETHODIMP CFusionGenotypeProbeSetResultsCOM::GetAlleleCallString(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(results.GetAlleleCallString());
	return S_OK;
}
