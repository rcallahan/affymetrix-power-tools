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
#include "FusionCHPHeaderCOM.h"
#include "FusionTagValuePairTypeCOM.h"
#include "COMStringUtils.h"
#include "BackgroundZoneInfoCOM.h"
#include "BackgroundZoneTypeCOM.h"
#include "calvin_files/fusion/src/FusionTagValuePairType.h"

// CFusionCHPHeaderCOM

STDMETHODIMP CFusionCHPHeaderCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCHPHeader
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_Cols(int* pVal)
{
	*pVal = header->GetCols();
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_Rows(int* pVal)
{
	*pVal = header->GetRows();
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_NumProbeSets(int* pVal)
{
	*pVal = header->GetNumProbeSets();
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_AssayType(AssayType* pVal)
{
	*pVal = (AssayType) header->GetAssayType();
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_ChipType(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(header->GetChipType());
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_AlgName(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(header->GetAlgName());
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_AlgVersion(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(header->GetAlgVersion());
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_AlgorithmParameters(VARIANT* pVal)
{
	VariantInit(pVal);
	pVal->vt = VT_ARRAY | VT_DISPATCH;
	pVal->parray = NULL;
	affymetrix_fusion_io::FusionTagValuePairTypeList values;
	header->AlgorithmParameters(values);
	int nparams = (int)values.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nparams;
	pVal->parray = SafeArrayCreate(VT_DISPATCH, 1, rgsaBound);
	long index=0;
	affymetrix_fusion_io::FusionTagValuePairTypeList::iterator it;
	for (it=values.begin(); it!=values.end(); ++it)
	{
		CComPtr<IFusionTagValuePairType> param;
		param.CoCreateInstance(CLSID_FusionTagValuePairType);
		CFusionTagValuePairTypeCOM *pParam = static_cast<CFusionTagValuePairTypeCOM *>(param.p);
		pParam->SetParam((*it));
		HRESULT hr = SafeArrayPutElement(pVal->parray, &index, param);
		++index;
	}
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_AlgorithmParameterCount(int* pVal)
{
	*pVal = header->AlgorithmParameterCount();
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_SummaryParameterCount(int* pVal)
{
	*pVal = header->SummaryParameterCount();
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_SummaryParameters(VARIANT* pVal)
{
	VariantInit(pVal);
	pVal->vt = VT_ARRAY | VT_DISPATCH;
	pVal->parray = NULL;
	affymetrix_fusion_io::FusionTagValuePairTypeList values;
	header->SummaryParameters(values);
	int nparams = (int)values.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nparams;
	pVal->parray = SafeArrayCreate(VT_DISPATCH, 1, rgsaBound);
	long index=0;
	affymetrix_fusion_io::FusionTagValuePairTypeList::iterator it;
	for (it=values.begin(); it!=values.end(); ++it)
	{
		CComPtr<IFusionTagValuePairType> param;
		param.CoCreateInstance(CLSID_FusionTagValuePairType);
		CFusionTagValuePairTypeCOM *pParam = static_cast<CFusionTagValuePairTypeCOM *>(param.p);
		pParam->SetParam((*it));
		HRESULT hr = SafeArrayPutElement(pVal->parray, &index, param);
		++index;
	}
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_ParentCellFile(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(header->GetParentCellFile());
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_ProgID(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(header->GetProgID());
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::GetAlgorithmParameter(BSTR tag, BSTR * pVal)
{
	*pVal = COMStringUtils::ConvertWString(header->GetAlgorithmParameter(COMStringUtils::ConvertWString(tag).c_str()));
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::GetSummaryParameter(BSTR tag, BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(header->GetSummaryParameter(COMStringUtils::ConvertWString(tag).c_str()));
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::GetAlgorithmParameterName(int index, BSTR * pVal)
{
	if (algParamNames.size() == 0)
	{
		affymetrix_fusion_io::FusionTagValuePairTypeList params;
		header->AlgorithmParameters(params);
		algParamNames.resize(params.size());
		affymetrix_fusion_io::FusionTagValuePairTypeList::iterator it;
		int i=0;
		for (it=params.begin(); it!=params.end(); ++it)
		{
			algParamNames[i++] = it->Tag;
		}
	}
	*pVal = COMStringUtils::ConvertWString(algParamNames[index]);
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::GetSummaryParameterName(int index, BSTR* pVal)
{
	if (sumParamNames.size() == 0)
	{
		affymetrix_fusion_io::FusionTagValuePairTypeList params;
		header->SummaryParameters(params);
		sumParamNames.resize(params.size());
		affymetrix_fusion_io::FusionTagValuePairTypeList::iterator it;
		int i=0;
		for (it=params.begin(); it!=params.end(); ++it)
		{
			sumParamNames[i++] = it->Tag;
		}
	}
	*pVal = COMStringUtils::ConvertWString(sumParamNames[index]);
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::Clear(void)
{
	header->Clear();
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::get_BackgroundZoneInfo(IBackgroundZoneInfo ** pVal)
{
	CComPtr<IBackgroundZoneInfo> info;
	info.CoCreateInstance(CLSID_BackgroundZoneInfo);
	CBackgroundZoneInfoCOM *pinfo = static_cast<CBackgroundZoneInfoCOM *>(info.p);
	header->GetBackgroundZoneInfo(pinfo->Zone());
	info.CopyTo(pVal);
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::GetBackgroundZones(VARIANT *zones)
{
	VariantInit(zones);
	zones->vt = VT_ARRAY | VT_DISPATCH;
	zones->parray = NULL;
	affxchp::BackgroundZoneTypeList czones;
	header->GetBackgroundZones(czones);
	int nzones = (int)czones.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nzones;
	zones->parray = SafeArrayCreate(VT_DISPATCH, 1, rgsaBound);
	long index=0;
	affxchp::BackgroundZoneTypeList::iterator it;
	for (it=czones.begin(); it!=czones.end(); ++it)
	{
		CComPtr<IBackgroundZoneType> z;
		z.CoCreateInstance(CLSID_BackgroundZoneType);
		CBackgroundZoneTypeCOM *pz = static_cast<CBackgroundZoneTypeCOM *>(z.p);
		pz->SetBg(*it);
		HRESULT hr = SafeArrayPutElement(zones->parray, &index, z);
		++index;
	}
	return S_OK;
}

STDMETHODIMP CFusionCHPHeaderCOM::GetBackgroundZone(int x, int y, IBackgroundZoneType ** pVal)
{
	affxchp::BackgroundZoneType b;
	header->GetBackgroundZone(b, x, y);
	CComPtr<IBackgroundZoneType> zone;
	zone.CoCreateInstance(CLSID_BackgroundZoneType);
	CBackgroundZoneTypeCOM *pzone = static_cast<CBackgroundZoneTypeCOM *>(zone.p);
	pzone->SetBg(b);
	zone.CopyTo(pVal);
	return S_OK;
}
