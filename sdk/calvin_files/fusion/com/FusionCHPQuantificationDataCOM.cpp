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
#include "FusionCHPDataCOM.h"
#include "FusionCHPQuantificationDataCOM.h"
#include "ProbeSetQuantificationDataCOM.h"
#include "calvin_files/fusion/src/FusionCHPQuantificationData.h"
#include "FusionTagValuePairTypeCOM.h"
#include "COMStringUtils.h"

// CFusionCHPQuantificationDataCOM

STDMETHODIMP CFusionCHPQuantificationDataCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCHPQuantificationData
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::get_FileId(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(chp->FileId());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::get_FileTypeIdentifiers(VARIANT* pVal)
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

STDMETHODIMP CFusionCHPQuantificationDataCOM::FromBase(IFusionCHPData * baseChp, VARIANT_BOOL* pVal)
{
	CFusionCHPDataCOM *pBaseChp = static_cast<CFusionCHPDataCOM*>(baseChp);
	chp = affymetrix_fusion_io::FusionCHPQuantificationData::FromBase(pBaseChp->GetChip());
	if (chp == NULL)
	{
		*pVal = VARIANT_FALSE;
		return S_OK;
	}
	*pVal = VARIANT_TRUE;
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::get_AlgName(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetAlgName());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::get_AlgVersion(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetAlgVersion());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::get_ArrayType(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetArrayType());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::get_AlgorithmParameters(VARIANT* pVal)
{
	if (chp == NULL)
		return E_FAIL;

	VariantInit(pVal);
	pVal->vt = VT_ARRAY | VT_DISPATCH;
	pVal->parray = NULL;
	affymetrix_calvin_parameter::ParameterNameValueTypeList params = chp->GetAlgParams();
	int nparams = (int)params.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nparams;
	pVal->parray = SafeArrayCreate(VT_DISPATCH, 1, rgsaBound);
	long index=0;
	affymetrix_calvin_parameter::ParameterNameValueTypeList::iterator it;
	for (it=params.begin(); it!=params.end(); ++it)
	{
		CComPtr<IFusionTagValuePairType> param;
		affymetrix_fusion_io::FusionTagValuePairType fparam;
		param.CoCreateInstance(CLSID_FusionTagValuePairType);
		CFusionTagValuePairTypeCOM *pParam = static_cast<CFusionTagValuePairTypeCOM *>(param.p);
		fparam.Tag = it->GetName();
		fparam.Value = it->ToString();
		pParam->SetParam(fparam);
		HRESULT hr = SafeArrayPutElement(pVal->parray, &index, param);
		++index;
	}
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::get_SummaryParameters(VARIANT* pVal)
{
	if (chp == NULL)
		return E_FAIL;

	VariantInit(pVal);
	pVal->vt = VT_ARRAY | VT_DISPATCH;
	pVal->parray = NULL;
	affymetrix_calvin_parameter::ParameterNameValueTypeList params = chp->GetSummaryParams();
	int nparams = (int)params.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nparams;
	pVal->parray = SafeArrayCreate(VT_DISPATCH, 1, rgsaBound);
	long index=0;
	affymetrix_calvin_parameter::ParameterNameValueTypeList::iterator it;
	for (it=params.begin(); it!=params.end(); ++it)
	{
		CComPtr<IFusionTagValuePairType> param;
		affymetrix_fusion_io::FusionTagValuePairType fparam;
		param.CoCreateInstance(CLSID_FusionTagValuePairType);
		CFusionTagValuePairTypeCOM *pParam = static_cast<CFusionTagValuePairTypeCOM *>(param.p);
		fparam.Tag = it->GetName();
		fparam.Value = it->ToString();
		pParam->SetParam(fparam);
		HRESULT hr = SafeArrayPutElement(pVal->parray, &index, param);
		++index;
	}
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::get_EntryCount(int* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = chp->GetEntryCount();
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::GetQuantificationEntry(int index, IProbeSetQuantificationData * pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CProbeSetQuantificationDataCOM *e = static_cast<CProbeSetQuantificationDataCOM *>(pVal);
	chp->GetQuantificationEntry(index, e->Data());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDataCOM::Close()
{
    delete chp;
    chp = NULL;
    return S_OK;
}
