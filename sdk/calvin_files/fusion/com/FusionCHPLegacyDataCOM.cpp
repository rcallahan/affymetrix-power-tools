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
#include "FusionCHPLegacyDataCOM.h"
#include "FusionCHPDataCOM.h"
#include "FusionExpressionProbeSetResultsCOM.h"
#include "FusionGenotypeProbeSetResultsCOM.h"
#include "FusionUniversalProbeSetResultsCOM.h"
#include "FusionResequencingResultsCOM.h"
#include "FusionCHPHeaderCOM.h"
#include "COMStringUtils.h"

// CFusionCHPLegacyDataCOM

STDMETHODIMP CFusionCHPLegacyDataCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCHPLegacyData
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::get_FileTypeIdentifier(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertString(chp->FileTypeIdentifier());
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::get_FileId(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(chp->FileId());
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::get_FileTypeIdentifiers(VARIANT* pVal)
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

STDMETHODIMP CFusionCHPLegacyDataCOM::FromBase(IFusionCHPData * baseChp, VARIANT_BOOL* pVal)
{
	CFusionCHPDataCOM *pBaseChp = static_cast<CFusionCHPDataCOM*>(baseChp);
	chp = affymetrix_fusion_io::FusionCHPLegacyData::FromBase(pBaseChp->GetChip());
	if (chp == NULL)
	{
		*pVal = VARIANT_FALSE;
		return S_OK;
	}
	*pVal = VARIANT_TRUE;
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::GetExpressionResults(int index, IFusionExpressionProbeSetResults * pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CFusionExpressionProbeSetResultsCOM *e = static_cast<CFusionExpressionProbeSetResultsCOM *>(pVal);
	chp->GetExpressionResults(index, e->Results());
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::GetGenotypingResults(int index, IFusionGenotypeProbeSetResults * pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CFusionGenotypeProbeSetResultsCOM *e = static_cast<CFusionGenotypeProbeSetResultsCOM *>(pVal);
	chp->GetGenotypingResults(index, e->Results());
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::GetUniversalResults(int index, IFusionUniversalProbeSetResults * pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CFusionUniversalProbeSetResultsCOM *e = static_cast<CFusionUniversalProbeSetResultsCOM *>(pVal);
	chp->GetUniversalResults(index, e->Results());
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::GetReseqResults(IFusionResequencingResults * pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CFusionResequencingResultsCOM *e = static_cast<CFusionResequencingResultsCOM *>(pVal);
	chp->GetReseqResults(e->Results());
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::Clear(void)
{
	if (chp)
	{
		chp->Clear();
		delete chp;
		chp = NULL;
	}
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::GetHeader(IFusionCHPHeader ** pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CComPtr<IFusionCHPHeader> header;
	header.CoCreateInstance(CLSID_FusionCHPHeader);
	CFusionCHPHeaderCOM *pHeader = static_cast<CFusionCHPHeaderCOM *>(header.p);
	pHeader->SetHeader(&chp->GetHeader());
	header.CopyTo(pVal);
	return S_OK;
}

STDMETHODIMP CFusionCHPLegacyDataCOM::GetProbeSetName(int index, BSTR *pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertString(chp->GetProbeSetName(index));
    return S_OK;
}
