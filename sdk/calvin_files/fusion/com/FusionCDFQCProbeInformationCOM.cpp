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
#include "FusionCDFQCProbeInformationCOM.h"


// CFusionCDFQCProbeInformationCOM

STDMETHODIMP CFusionCDFQCProbeInformationCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCDFQCProbeInformation
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCDFQCProbeInformationCOM::get_X(int* pVal)
{
	*pVal = probe.GetX();
	return S_OK;
}

STDMETHODIMP CFusionCDFQCProbeInformationCOM::get_Y(int* pVal)
{
	*pVal = probe.GetY();
	return S_OK;
}

STDMETHODIMP CFusionCDFQCProbeInformationCOM::get_PLen(int* pVal)
{
	*pVal = probe.GetPLen();
	return S_OK;
}

STDMETHODIMP CFusionCDFQCProbeInformationCOM::IsPerfectMatchProbe(VARIANT_BOOL* pVal)
{
	*pVal = (probe.IsPerfectMatchProbe() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCDFQCProbeInformationCOM::IsBackgroundProbe(VARIANT_BOOL* pVal)
{
	*pVal = (probe.IsBackgroundProbe() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}
