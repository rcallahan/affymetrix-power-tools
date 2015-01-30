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
#include "BackgroundZoneInfoCOM.h"
#include "BackgroundZoneTypeCOM.h"


// CBackgroundZoneInfoCOM

STDMETHODIMP CBackgroundZoneInfoCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IBackgroundZoneInfo
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CBackgroundZoneInfoCOM::get_number_zones(int* pVal)
{
	*pVal = zone.number_zones;
	return S_OK;
}

STDMETHODIMP CBackgroundZoneInfoCOM::get_smooth_factor(float* pVal)
{
	*pVal = zone.smooth_factor;
	return S_OK;
}

STDMETHODIMP CBackgroundZoneInfoCOM::get_zones(VARIANT* pVal)
{
	VariantInit(pVal);
	pVal->vt = VT_ARRAY | VT_DISPATCH;
	pVal->parray = NULL;
	int npVal = (int)zone.zones.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = npVal;
	pVal->parray = SafeArrayCreate(VT_DISPATCH, 1, rgsaBound);
	long index=0;
	affxchp::BackgroundZoneTypeList::iterator it;
	for (it=zone.zones.begin(); it!=zone.zones.end(); ++it)
	{
		CComPtr<IBackgroundZoneType> z;
		z.CoCreateInstance(CLSID_BackgroundZoneType);
		CBackgroundZoneTypeCOM *pz = static_cast<CBackgroundZoneTypeCOM *>(z.p);
		pz->SetBg(*it);
		HRESULT hr = SafeArrayPutElement(pVal->parray, &index, z);
		++index;
	}
	return S_OK;
}
