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
#include "FusionTagValuePairTypeCOM.h"
#include "COMStringUtils.h"

// CFusionTagValuePairTypeCOM

STDMETHODIMP CFusionTagValuePairTypeCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionTagValuePairType
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionTagValuePairTypeCOM::get_Tag(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(param.Tag);
	return S_OK;
}

STDMETHODIMP CFusionTagValuePairTypeCOM::get_Value(BSTR* pVal)
{
	std::wstring s = param.DetailedType().ToString();
	if (s.length() > 0)
		*pVal = COMStringUtils::ConvertWString(s);
	else
		*pVal = COMStringUtils::ConvertWString(param.Value);
	return S_OK;
}
