////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
#include "ProbeSetQuantificationDetectionDataCOM.h"
#include "COMStringUtils.h"

// CProbeSetQuantificationDetectionDataCOM

STDMETHODIMP CProbeSetQuantificationDetectionDataCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IProbeSetQuantificationDetectionData
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CProbeSetQuantificationDetectionDataCOM::get_name(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(data.name);
	return S_OK;
}

STDMETHODIMP CProbeSetQuantificationDetectionDataCOM::get_quantification(float* pVal)
{
	*pVal = data.quantification;
	return S_OK;
}

STDMETHODIMP CProbeSetQuantificationDetectionDataCOM::get_pvalue(float* pVal)
{
	*pVal = data.pvalue;
	return S_OK;
}

STDMETHODIMP CProbeSetQuantificationDetectionDataCOM::get_identifier(int* pVal)
{
	*pVal = data.id;
    return S_OK;
}
