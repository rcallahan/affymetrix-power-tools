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
#include "FusionCELFileEntryTypeCOM.h"


// CFusionCELFileEntryTypeCOM

STDMETHODIMP CFusionCELFileEntryTypeCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCELFileEntryType
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCELFileEntryTypeCOM::get_Intensity(float* pVal)
{
	*pVal = entry.Intensity;
	return S_OK;
}

STDMETHODIMP CFusionCELFileEntryTypeCOM::get_Stdv(float* pVal)
{
	*pVal = entry.Stdv;
	return S_OK;
}

STDMETHODIMP CFusionCELFileEntryTypeCOM::get_Pixels(short* pVal)
{
	*pVal = entry.Pixels;
	return S_OK;
}
