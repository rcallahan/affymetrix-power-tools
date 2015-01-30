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
#include "FusionCHPQuantificationDetectionDataCOM.h"
#include "FusionCHPDataCOM.h"
#include "ProbeSetQuantificationDetectionDataCOM.h"
#include "calvin_files/fusion/src/FusionCHPQuantificationDetectionData.h"
#include "FusionTagValuePairTypeCOM.h"
#include "COMStringUtils.h"


// CFusionCHPQuantificationDetectionData

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCHPQuantificationDetectionData
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}


STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::get_FileId(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(chp->FileId());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::get_FileTypeIdentifiers(VARIANT* pVal)
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

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::FromBase(IFusionCHPData * baseChp, VARIANT_BOOL* pVal)
{
	CFusionCHPDataCOM *pBaseChp = static_cast<CFusionCHPDataCOM*>(baseChp);
	chp = affymetrix_fusion_io::FusionCHPQuantificationDetectionData::FromBase(pBaseChp->GetChip());
	if (chp == NULL)
	{
		*pVal = VARIANT_FALSE;
		return S_OK;
	}
	*pVal = VARIANT_TRUE;
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::get_AlgName(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetAlgName());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::get_AlgVersion(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetAlgVersion());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::get_ArrayType(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetArrayType());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::get_AlgorithmParameters(VARIANT* pVal)
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

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::get_SummaryParameters(VARIANT* pVal)
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

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::get_EntryCount(int* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = chp->GetEntryCount();
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::GetQuantificationDetectionEntry(int index, IProbeSetQuantificationDetectionData * pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CProbeSetQuantificationDetectionDataCOM *e = static_cast<CProbeSetQuantificationDetectionDataCOM *>(pVal);
	chp->GetQuantificationDetectionEntry(index, e->Data());
	return S_OK;
}

STDMETHODIMP CFusionCHPQuantificationDetectionDataCOM::Close()
{
    delete chp;
    chp = NULL;
    return S_OK;
}
