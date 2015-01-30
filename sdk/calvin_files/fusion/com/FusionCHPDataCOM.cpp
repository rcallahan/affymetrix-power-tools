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
#include "FusionCHPDataCOM.h"
#include "COMStringUtils.h"


// CFusionCHPDataCOM

STDMETHODIMP CFusionCHPDataCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCHPData
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCHPDataCOM::get_FileTypeIdentifier(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertString(chp->FileTypeIdentifier());
	return S_OK;
}

STDMETHODIMP CFusionCHPDataCOM::get_FileId(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(chp->FileId());
	return S_OK;
}

STDMETHODIMP CFusionCHPDataCOM::get_FileTypeIdentifiers(VARIANT* pVal)
{
	if (chp == NULL)
		return E_FAIL;

	VariantInit(pVal);
	pVal->vt = VT_ARRAY | VT_BSTR;
	pVal->parray = NULL;
	affymetrix_calvin_utilities::AffymetrixGuidTypeList &ids = chp->FileTypeIdentifiers();
	int nids = (int)ids.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nids;
	pVal->parray = SafeArrayCreate(VT_BSTR, 1, rgsaBound);
	long index=0;
	affymetrix_calvin_utilities::AffymetrixGuidTypeList::iterator it;
	for (it=ids.begin(); it!=ids.end(); ++it)
	{
		BSTR id = COMStringUtils::ConvertString(*it);
		HRESULT hr = SafeArrayPutElement(pVal->parray, &index, id);
		++index;
	}
	return S_OK;
}
