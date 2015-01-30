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
#include "FusionCHPMultiDataDataCOM.h"
#include "ProbeSetMultiDataGenotypeDataCOM.h"
#include "ProbeSetMultiDataExpressionDataCOM.h"
#include "ProbeSetMultiDataCopyNumberDataCOM.h"
#include "ProbeSetMultiDataCytoRegionDataCOM.h"
#include "ProbeSetMultiDataCopyNumberVariationRegionDataCOM.h"
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"
#include "FusionTagValuePairTypeCOM.h"
#include "COMStringUtils.h"
#include <vector>
#include <string>
#include <comdef.h>

using namespace std;

// CFusionCHPMultiDataDataCOM

STDMETHODIMP CFusionCHPMultiDataDataCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCHPMultiDataData
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::get_FileId(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(chp->FileId());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::get_FileTypeIdentifiers(VARIANT* pVal)
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

STDMETHODIMP CFusionCHPMultiDataDataCOM::FromBase(IFusionCHPData * baseChp, VARIANT_BOOL* pVal)
{
	CFusionCHPDataCOM *pBaseChp = static_cast<CFusionCHPDataCOM*>(baseChp);
	chp = affymetrix_fusion_io::FusionCHPMultiDataData::FromBase(pBaseChp->GetChip());
	if (chp == NULL)
	{
		*pVal = VARIANT_FALSE;
		return S_OK;
	}
	*pVal = VARIANT_TRUE;
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::get_AlgName(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetAlgName());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::get_AlgVersion(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetAlgVersion());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::get_ArrayType(BSTR* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	*pVal = COMStringUtils::ConvertWString(chp->GetArrayType());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::get_AlgorithmParameters(VARIANT* pVal)
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

STDMETHODIMP CFusionCHPMultiDataDataCOM::get_SummaryParameters(VARIANT* pVal)
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

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetEntryCount(MultiDataType dataType, int* pVal)
{
	if (chp == NULL)
		return E_FAIL;
    *pVal = chp->GetEntryCount((affymetrix_calvin_io::MultiDataType)dataType);
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetGenotypeEntry(MultiDataType dataType, int index, IProbeSetMultiDataGenotypeData* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CProbeSetMultiDataGenotypeDataCOM *e = static_cast<CProbeSetMultiDataGenotypeDataCOM *>(pVal);
	chp->GetGenotypeEntry((affymetrix_calvin_io::MultiDataType)dataType, index, e->Data());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetGenoCall(MultiDataType dataType, int index, int * pVal)
{
	if (chp == NULL)
		return E_FAIL;
    *pVal = chp->GetGenoCall((affymetrix_calvin_io::MultiDataType)dataType, index);
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetGenoConfidence(MultiDataType dataType, int index, float * pVal)
{
	if (chp == NULL)
		return E_FAIL;
    *pVal = chp->GetGenoConfidence((affymetrix_calvin_io::MultiDataType)dataType, index);
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetProbeSetName(MultiDataType dataType, int index, BSTR * pVal)
{
	if (chp == NULL)
		return E_FAIL;
    *pVal = COMStringUtils::ConvertString(chp->GetProbeSetName((affymetrix_calvin_io::MultiDataType)dataType, index));
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetExpressionEntry(MultiDataType dataType, int index, IProbeSetMultiDataExpressionData * pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CProbeSetMultiDataExpressionDataCOM *e = static_cast<CProbeSetMultiDataExpressionDataCOM *>(pVal);
	chp->GetExpressionEntry((affymetrix_calvin_io::MultiDataType)dataType, index, e->Data());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetExpressionQuantification(MultiDataType dataType, int index, float * pVal)
{
	if (chp == NULL)
		return E_FAIL;
    *pVal = chp->GetExpressionQuantification((affymetrix_calvin_io::MultiDataType)dataType, index);
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetNumMetricColumns(MultiDataType dataType, int * pVal)
{
	if (chp == NULL)
		return E_FAIL;
    *pVal = chp->GetNumMetricColumns((affymetrix_calvin_io::MultiDataType)dataType);
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetMetricColumnName(MultiDataType dataType, int index, BSTR * pVal)
{
	if (chp == NULL)
		return E_FAIL;
    *pVal = COMStringUtils::ConvertWString(chp->GetMetricColumnName((affymetrix_calvin_io::MultiDataType)dataType, index));
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::InitializeGetGenotypeEntries(MultiDataType dataType, BSTR* forcedName, BSTR *signalAName, BSTR *signalBName)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetNumMetricColumns(dt);
    forceIndex = -1;
    sigAIndex = -1;
    sigBIndex = -1;
    for (int i=0; i<n; i++)
    {
        wstring name = chp->GetMetricColumnName(dt, i);
        if (name.find(L"Forced") == 0)
        {
            forceIndex = i;
            *forcedName = COMStringUtils::ConvertWString(name);
        }
        else if (sigAIndex == -1)
        {
            sigAIndex = i;
            *signalAName = COMStringUtils::ConvertWString(name);
        }
        else if (sigBIndex == -1)
        {
            sigBIndex = i;
            *signalBName = COMStringUtils::ConvertWString(name);
        }
    }
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetGenotypeEntries(MultiDataType dataType, VARIANT* names, VARIANT* calls, VARIANT* confidences, VARIANT* forcedcalls, VARIANT* signalsA, VARIANT* signalsB, int startIndex, int *pVal)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetEntryCount(dt);
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(names->parray, 1, &ubound);
    SafeArrayGetLBound(names->parray, 1, &lbound);
    int len = ubound - lbound + 1;
    int count = (len + startIndex <= n ? len : n - startIndex);
    *pVal = count;
    float fvalue;
    u_int8_t ubvalue;
    affymetrix_calvin_data::ProbeSetMultiDataGenotypeData entry;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
        chp->GetGenotypeEntry(dt, i+startIndex, entry);

        _bstr_t bstr = entry.name.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &index, bsstr);
        SysFreeString(bsstr);

        SafeArrayPutElement(calls->parray, &index, &entry.call);

        SafeArrayPutElement(confidences->parray, &index, &entry.confidence);

        ubvalue = entry.metrics[forceIndex].GetValueUInt8();
        SafeArrayPutElement(forcedcalls->parray, &index, &ubvalue);

        fvalue = entry.metrics[sigAIndex].GetValueFloat();
        SafeArrayPutElement(signalsA->parray, &index, &fvalue);

        fvalue = entry.metrics[sigBIndex].GetValueFloat();
        SafeArrayPutElement(signalsB->parray, &index, &fvalue);
    }

    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetGenotypeEntriesGivenIndicies(MultiDataType dataType, VARIANT* indicies, VARIANT* names, VARIANT* calls, VARIANT* confidences, VARIANT* forcedcalls, VARIANT* signalsA, VARIANT* signalsB)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(indicies->parray, 1, &ubound);
    SafeArrayGetLBound(indicies->parray, 1, &lbound);
    long len = ubound - lbound + 1;
    float fvalue;
    u_int8_t ubvalue;
    affymetrix_calvin_data::ProbeSetMultiDataGenotypeData entry;
    for (long i=0; i<len; i++)
    {
        ::SafeArrayGetElement(indicies->parray, &i, &index);
        chp->GetGenotypeEntry(dt, index, entry);

        _bstr_t bstr = entry.name.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &i, bsstr);
        SysFreeString(bsstr);

        SafeArrayPutElement(calls->parray, &i, &entry.call);

        SafeArrayPutElement(confidences->parray, &i, &entry.confidence);

        ubvalue = entry.metrics[forceIndex].GetValueUInt8();
        SafeArrayPutElement(forcedcalls->parray, &i, &ubvalue);

        fvalue = entry.metrics[sigAIndex].GetValueFloat();
        SafeArrayPutElement(signalsA->parray, &i, &fvalue);

        fvalue = entry.metrics[sigBIndex].GetValueFloat();
        SafeArrayPutElement(signalsB->parray, &i, &fvalue);
    }

    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetProbesetNamesAndGenotypeCalls(MultiDataType dataType, VARIANT* names, VARIANT* calls, int startIndex, int *pVal)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetEntryCount(dt);
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(names->parray, 1, &ubound);
    SafeArrayGetLBound(names->parray, 1, &lbound);
    int len = ubound - lbound + 1;
    int count = (len + startIndex <= n ? len : n - startIndex);
    *pVal = count;

    u_int8_t genotypeCall;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
        genotypeCall = chp->GetGenoCall(dt, i+startIndex);
        SafeArrayPutElement(calls->parray, &index, &genotypeCall);
    }

    string probeSetName;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
        probeSetName = chp->GetProbeSetName(dt, i+startIndex);
        _bstr_t bstr = probeSetName.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &index, bsstr);
        SysFreeString(bsstr);
    }

    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetProbesetNamesAndExpressionSignals(MultiDataType dataType, VARIANT* names, VARIANT* signals, int startIndex, int *pVal)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetEntryCount(dt);
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(names->parray, 1, &ubound);
    SafeArrayGetLBound(names->parray, 1, &lbound);
    int len = ubound - lbound + 1;
    int count = (len + startIndex <= n ? len : n - startIndex);
    *pVal = count;

    float expSignal;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
        expSignal = chp->GetExpressionQuantification(dt, i+startIndex);
        SafeArrayPutElement(signals->parray, &index, &expSignal);
    }

    string probeSetName;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
        probeSetName = chp->GetProbeSetName(dt, i+startIndex);
        _bstr_t bstr = probeSetName.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &index, bsstr);
        SysFreeString(bsstr);
    }

    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::Close()
{
    delete chp;
    chp = NULL;
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetCopyNumberEntry(MultiDataType dataType, int index, IProbeSetMultiDataCopyNumberData* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CProbeSetMultiDataCopyNumberDataCOM *e = static_cast<CProbeSetMultiDataCopyNumberDataCOM *>(pVal);
	chp->GetCopyNumberEntry((affymetrix_calvin_io::MultiDataType)dataType, index, e->Data());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetCopyNumberEntryLog2Ratio(MultiDataType dataType, int index, float* val)
{
	if (chp == NULL)
		return E_FAIL;
	chp->GetCopyNumberEntryLog2Ratio((affymetrix_calvin_io::MultiDataType)dataType, index, log2RatioIndex, val);
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetCytoRegionEntry(MultiDataType dataType, int index, IProbeSetMultiDataCytoRegionData* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CProbeSetMultiDataCytoRegionDataCOM *e = static_cast<CProbeSetMultiDataCytoRegionDataCOM *>(pVal);
	chp->GetCytoRegionEntry((affymetrix_calvin_io::MultiDataType)dataType, index, e->Data());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetCytoRegionEntries(MultiDataType dataType, VARIANT* names, VARIANT* calls, VARIANT* confidences, VARIANT* callsLOH, VARIANT* confidencesLOH, int startIndex, int *pVal)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetEntryCount(dt);
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(names->parray, 1, &ubound);
    SafeArrayGetLBound(names->parray, 1, &lbound);
    int len = ubound - lbound + 1;
    int count = (len + startIndex <= n ? len : n - startIndex);
    *pVal = count;
    affymetrix_calvin_data::ProbeSetMultiDataCytoRegionData entry;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
        chp->GetCytoRegionEntry(dt, i+startIndex, entry);
        _bstr_t bstr = entry.name.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &index, bsstr);
        SysFreeString(bsstr);
        SafeArrayPutElement(calls->parray, &index, &entry.call);
        SafeArrayPutElement(confidences->parray, &index, &entry.confidenceScore);

		for(ParameterNameValueTypeIt it = entry.metrics.begin(); it != entry.metrics.end(); it++)
		{
			wstring name = it->GetName();
			if(name == L"CallLOH")
			{
				int32_t callLOH = it->GetValueInt32();
				SafeArrayPutElement(callsLOH->parray, &index, &callLOH);
			}
			else if((name == L"ConfidenceLOH"))
			{
				float confidenceLOH = it->GetValueFloat();
				SafeArrayPutElement(confidencesLOH->parray, &index, &confidenceLOH);
			}
		}
    }
    return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetDataSetHeaderParameters(MultiDataType dataType, VARIANT* pVal)
{
	if (chp == NULL)
		return E_FAIL;

	VariantInit(pVal);
	pVal->vt = VT_ARRAY | VT_DISPATCH;
	pVal->parray = NULL;
	affymetrix_calvin_io::DataSetHeader *dsh = chp->GetDataSetHeader((affymetrix_calvin_io::MultiDataType)dataType);
	affymetrix_calvin_parameter::ParameterNameValueTypeConstIt begin;
	affymetrix_calvin_parameter::ParameterNameValueTypeConstIt end;
	dsh->GetNameValIterators(begin, end);
	int nparams = dsh->GetNameValParamCnt();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nparams;
	pVal->parray = SafeArrayCreate(VT_DISPATCH, 1, rgsaBound);
	long index=0;
	for (ParameterNameValueTypeConstIt it=begin; it!=end; ++it)
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

STDMETHODIMP CFusionCHPMultiDataDataCOM::InitializeGetCNStateEntries(MultiDataType dataType)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetNumMetricColumns(dt);
    cnStateIndex = -1;
    for (int i=0; i<n; i++)
    {
        wstring name = chp->GetMetricColumnName(dt, i);
        if (name == L"CNState")
        {
            cnStateIndex = i;
        }
    }
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetCopyNumberState(MultiDataType dataType, VARIANT* names, VARIANT* chromosomes, VARIANT* positions, VARIANT* cnStates, int startIndex, int *pVal)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetEntryCount(dt);
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(names->parray, 1, &ubound);
    SafeArrayGetLBound(names->parray, 1, &lbound);
    int len = ubound - lbound + 1;
    int count = (len + startIndex <= n ? len : n - startIndex);
    *pVal = count;
    u_int8_t ubvalue;
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData entry;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
		chp->GetCopyNumberEntry(dt, i+startIndex, entry);

        _bstr_t bstr = entry.name.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &index, bsstr);
        SysFreeString(bsstr);

		ubvalue = entry.chr;
        SafeArrayPutElement(chromosomes->parray, &index, &ubvalue);

		SafeArrayPutElement(positions->parray, &index, &entry.position);

		if (entry.metrics[cnStateIndex].GetParameterType() == ParameterNameValueType::FloatType)
			ubvalue = (u_int8_t)(entry.metrics[cnStateIndex].GetValueFloat() + 0.1f);
		else
			ubvalue = entry.metrics[cnStateIndex].GetValueUInt8();
        SafeArrayPutElement(cnStates->parray, &index, &ubvalue);
    }

	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::InitializeGetLOHStateEntries(MultiDataType dataType)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetNumMetricColumns(dt);
    cnStateIndex = -1;
    for (int i=0; i<n; i++)
    {
        wstring name = chp->GetMetricColumnName(dt, i);
        if (name == L"LOHState")
        {
            lohStateIndex = i;
        }
    }
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetLOHState(MultiDataType dataType, VARIANT* names, VARIANT* chromosomes, VARIANT* positions, VARIANT* lohStates, int startIndex, int *pVal)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetEntryCount(dt);
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(names->parray, 1, &ubound);
    SafeArrayGetLBound(names->parray, 1, &lbound);
    int len = ubound - lbound + 1;
    int count = (len + startIndex <= n ? len : n - startIndex);
    *pVal = count;
    u_int8_t ubvalue;
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData entry;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
		chp->GetCopyNumberEntry(dt, i+startIndex, entry);

        _bstr_t bstr = entry.name.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &index, bsstr);
        SysFreeString(bsstr);

		ubvalue = entry.chr;
        SafeArrayPutElement(chromosomes->parray, &index, &ubvalue);

		SafeArrayPutElement(positions->parray, &index, &entry.position);

		if (entry.metrics[lohStateIndex].GetParameterType() == ParameterNameValueType::FloatType)
			ubvalue = (u_int8_t)(entry.metrics[lohStateIndex].GetValueFloat() + 0.1f);
		else
			ubvalue = entry.metrics[lohStateIndex].GetValueUInt8();
        SafeArrayPutElement(lohStates->parray, &index, &ubvalue);
    }

	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::InitializeGetLog2RatioEntries(MultiDataType dataType)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetNumMetricColumns(dt);
    cnStateIndex = -1;
    for (int i=0; i<n; i++)
    {
        wstring name = chp->GetMetricColumnName(dt, i);
        if (name == L"Log2Ratio")
        {
            log2RatioIndex = i;
			break;
        }
    }
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::GetLog2Ratio(MultiDataType dataType, int startIndex, int count, VARIANT* names, VARIANT* positions, VARIANT* values)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData entry;
    for (int i=0; i<count; i++)
    {
        long index = (long) i;
		chp->GetCopyNumberEntry(dt, i+startIndex, entry);

		_bstr_t bstr = entry.name.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &index, bsstr);
        SysFreeString(bsstr);

		SafeArrayPutElement(positions->parray, &index, &entry.position);

		float val = entry.metrics[log2RatioIndex].GetValueFloat();
        SafeArrayPutElement(values->parray, &index, &val);
    }
	return S_OK;
}


STDMETHODIMP CFusionCHPMultiDataDataCOM::GetLog2Ratio2(MultiDataType dataType, int startIndex, int count, VARIANT* values)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberData entry;
    for (int i=0; i<count; i++)
    {
        long index = (long) i;
		float val = 0;
		chp->GetCopyNumberEntryLog2Ratio(dt, i+startIndex, log2RatioIndex, &val);
        SafeArrayPutElement(values->parray, &index, &val);
    }
	return S_OK;
}


STDMETHODIMP CFusionCHPMultiDataDataCOM::GetCopyNumberVariationRegionEntry(MultiDataType dataType, int index, IProbeSetMultiDataCopyNumberVariationRegionData* pVal)
{
	if (chp == NULL)
		return E_FAIL;
	CProbeSetMultiDataCopyNumberVariationRegionDataCOM *e = static_cast<CProbeSetMultiDataCopyNumberVariationRegionDataCOM *>(pVal);
	chp->GetCopyNumberVariationRegionEntry((affymetrix_calvin_io::MultiDataType)dataType, index, e->Data());
	return S_OK;
}

STDMETHODIMP CFusionCHPMultiDataDataCOM::InitializeGetCopyNumberVariationEntries(MultiDataType dataType)
{
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
    int n = chp->GetNumMetricColumns(dt);   
    signalIndex = -1;
    callIndex = -1;
    confidenceIndex = -1;
    for (int i=0; i<n; i++)
    {
        wstring name = chp->GetMetricColumnName(dt, i);
        if (name == L"Signal")
        {
            signalIndex = i;
        }
        else if (name == L"Call")
        {
            callIndex = i;
        }
        else if (name == L"Confidence")
        {
            confidenceIndex = i;
        }
    }
	return S_OK;
}

/*!
 *  returns the region names, calls, signals and confidenceScores in a cnvchp file
 */ 
STDMETHODIMP CFusionCHPMultiDataDataCOM::GetCopyNumberVariationRegionNamesAndResults(MultiDataType dataType, VARIANT* names, 
                                                                       VARIANT* calls, VARIANT* signals, VARIANT* confidenceScores, int *pVal) 
{
    if (chp == NULL)
		return E_FAIL;
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
   
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(names->parray, 1, &ubound);
    SafeArrayGetLBound(names->parray, 1, &lbound);
    int len = ubound - lbound + 1; 
    int count = chp->GetEntryCount(dt);
    if (len < count)
        return S_OK;
    
    *pVal = count;   
    string name;    
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberVariationRegionData entry;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
        chp->GetCopyNumberVariationRegionEntry(dt, i, entry);
        bstr_t bstr = entry.name.c_str();
        BSTR bsstr = bstr.copy();
        SafeArrayPutElement(names->parray, &index, bsstr);
        SysFreeString(bsstr);
        SafeArrayPutElement(calls->parray, &index, &entry.call);
        SafeArrayPutElement(signals->parray, &index, &entry.signal);
        SafeArrayPutElement(confidenceScores->parray, &index, &entry.confidenceScore);
    }

    return S_OK;
}

/*!
 *  returns the calls, signals and confidenceScores in a cnvchp file
 */ 
STDMETHODIMP CFusionCHPMultiDataDataCOM::GetCopyNumberVariationResults(MultiDataType dataType, VARIANT* calls, VARIANT* signals, 
                                                                       VARIANT* confidenceScores, int *pVal) 
{
    if (chp == NULL)
		return E_FAIL;
    affymetrix_calvin_io::MultiDataType dt = (affymetrix_calvin_io::MultiDataType) dataType;
   
    long ubound = 0;
    long lbound = 0;
    long index;
    SafeArrayGetUBound(calls->parray, 1, &ubound);
    SafeArrayGetLBound(calls->parray, 1, &lbound);
    int len = ubound - lbound + 1; 
    int count = chp->GetEntryCount(dt);
    if (len < count)
        return S_OK;
    
    *pVal = count;   
    string name;    
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberVariationRegionData entry;
    for (int i=0; i<count; i++)
    {
        index = (long) i;
        chp->GetCopyNumberVariationRegionEntry(dt, i, entry);        
        SafeArrayPutElement(calls->parray, &index, &entry.call);
        SafeArrayPutElement(signals->parray, &index, &entry.signal);
        SafeArrayPutElement(confidenceScores->parray, &index, &entry.confidenceScore);
    }

    return S_OK;
}


