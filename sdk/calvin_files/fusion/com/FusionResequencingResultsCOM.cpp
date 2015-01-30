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
#include "FusionResequencingResultsCOM.h"
#include "FusionForceCallTypeCOM.h"
#include "FusionBaseCallTypeCOM.h"

// CFusionResequencingResultsCOM

STDMETHODIMP CFusionResequencingResultsCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionResequencingResults
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionResequencingResultsCOM::Clear(void)
{
	results.Clear();
	return S_OK;
}

STDMETHODIMP CFusionResequencingResultsCOM::GetCalledBase(int index, CHAR* pVal)
{
	*pVal = results.GetCalledBase(index);
	return S_OK;
}

STDMETHODIMP CFusionResequencingResultsCOM::get_CalledBasesSize(int* pVal)
{
	*pVal = results.GetCalledBasesSize();
	return S_OK;
}

STDMETHODIMP CFusionResequencingResultsCOM::GetScore(int index, float *pVal)
{
	*pVal = results.GetScore(index);
	return S_OK;
}

STDMETHODIMP CFusionResequencingResultsCOM::get_ScoresSize(int* pVal)
{
	*pVal = results.GetScoresSize();
	return S_OK;
}

STDMETHODIMP CFusionResequencingResultsCOM::GetForceCall(int index, IFusionForceCallType ** pVal)
{
	CComPtr<IFusionForceCallType> call;
	call.CoCreateInstance(CLSID_FusionForceCallType);
	CFusionForceCallTypeCOM *pcall = static_cast<CFusionForceCallTypeCOM *>(call.p);
	pcall->Force() = results.GetForceCall(index);
	call.CopyTo(pVal);
	return S_OK;
}

STDMETHODIMP CFusionResequencingResultsCOM::get_ForceCallsSize(int* pVal)
{
	*pVal = results.GetForceCallsSize();
	return S_OK;
}

STDMETHODIMP CFusionResequencingResultsCOM::GetOrigCall(int index, IFusionBaseCallType ** pVal)
{
	CComPtr<IFusionBaseCallType> call;
	call.CoCreateInstance(CLSID_FusionBaseCallType);
	CFusionBaseCallTypeCOM *pcall = static_cast<CFusionBaseCallTypeCOM *>(call.p);
	pcall->Call() = results.GetOrigCall(index);
	call.CopyTo(pVal);
	return S_OK;
}

STDMETHODIMP CFusionResequencingResultsCOM::get_OrigCallsSize(int* pVal)
{
	*pVal = results.GetOrigCallsSize();
	return S_OK;
}
