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
#include "FusionCDFProbeInformationCOM.h"

// CFusionCDFProbeInformationCOM

STDMETHODIMP CFusionCDFProbeInformationCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCDFProbeInformation
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCDFProbeInformationCOM::get_ListIndex(int* pVal)
{
	*pVal = probe.GetListIndex();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeInformationCOM::get_Expos(int* pVal)
{
	*pVal = probe.GetExpos();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeInformationCOM::get_X(int* pVal)
{
	*pVal = probe.GetX();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeInformationCOM::get_Y(int* pVal)
{
	*pVal = probe.GetY();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeInformationCOM::get_PBase(CHAR* pVal)
{
	*pVal = probe.GetPBase();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeInformationCOM::get_TBase(CHAR* pVal)
{
	*pVal = probe.GetTBase();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeInformationCOM::get_ProbeLength(USHORT* pVal)
{
	*pVal = probe.GetProbeLength();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeInformationCOM::get_ProbeGrouping(USHORT* pVal)
{
	*pVal = probe.GetProbeGrouping();
	return S_OK;
}
