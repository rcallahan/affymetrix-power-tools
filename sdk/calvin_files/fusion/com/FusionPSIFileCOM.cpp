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
#include "FusionPSIFileCOM.h"
#include "COMStringUtils.h"

using namespace affxpsi;

// CFusionPSIFileCOM

STDMETHODIMP CFusionPSIFileCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionPSIFile
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionPSIFileCOM::get_FileName(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(psi.GetFileName());
	return S_OK;
}

STDMETHODIMP CFusionPSIFileCOM::put_FileName(BSTR newVal)
{
	psi.SetFileName(COMStringUtils::ConvertString(newVal).c_str());
	return S_OK;
}

STDMETHODIMP CFusionPSIFileCOM::get_ProbeSetCount(int* pVal)
{
	*pVal = psi.GetProbeSetCount();
	return S_OK;
}

STDMETHODIMP CFusionPSIFileCOM::GetProbeSetName(int index, BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(psi.GetProbeSetName(index));
	return S_OK;
}

STDMETHODIMP CFusionPSIFileCOM::GetProbePairs(int index, int * pVal)
{
	*pVal = psi.GetProbePairs(index);
	return S_OK;
}

STDMETHODIMP CFusionPSIFileCOM::Read(VARIANT_BOOL* pVal)
{
	*pVal = (psi.Read() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionPSIFileCOM::Exists(VARIANT_BOOL* pVal)
{
	*pVal = (psi.Exists() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionPSIFileCOM::Clear()
{
	psi.Clear();
	return S_OK;
}
