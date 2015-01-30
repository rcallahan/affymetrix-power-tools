////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
// ProbeSetMultiDataCopyNumberVariationRegionDataCOM.h : Declaration of the CProbeSetMultiDataCopyNumberVariationRegionDataCOM

#pragma once
#include "resource.h"       // main symbols

#include "affx_fusion_com.h"
#include "calvin_files/data/src/ProbeSetMultiDataData.h"


#if defined(_WIN32_WCE) && !defined(_CE_DCOM) && !defined(_CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA)
#error "Single-threaded COM objects are not properly supported on Windows CE platform, such as the Windows Mobile platforms that do not include full DCOM support. Define _CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA to force ATL to support creating single-thread COM object's and allow use of it's single-threaded COM object implementations. The threading model in your rgs file was set to 'Free' as that is the only threading model supported in non DCOM Windows CE platforms."
#endif



// CProbeSetMultiDataCopyNumberVariationRegionDataCOM

class ATL_NO_VTABLE CProbeSetMultiDataCopyNumberVariationRegionDataCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CProbeSetMultiDataCopyNumberVariationRegionDataCOM, &CLSID_ProbeSetMultiDataCopyNumberVariationRegionData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IProbeSetMultiDataCopyNumberVariationRegionData, &IID_IProbeSetMultiDataCopyNumberVariationRegionData, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CProbeSetMultiDataCopyNumberVariationRegionDataCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_PROBESETMULTIDATACOPYNUMBERVARIATIONREGIONDATA)
//DECLARE_REGISTRY_RESOURCEID(IDR_PROBESETMULTIDATACOPYNUMBERVARIATIONREGIONDATACOM)


BEGIN_COM_MAP(CProbeSetMultiDataCopyNumberVariationRegionDataCOM)
	COM_INTERFACE_ENTRY(IProbeSetMultiDataCopyNumberVariationRegionData)
	COM_INTERFACE_ENTRY(IDispatch)
	COM_INTERFACE_ENTRY(ISupportErrorInfo)
END_COM_MAP()

// ISupportsErrorInfo
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);


	DECLARE_PROTECT_FINAL_CONSTRUCT()

	HRESULT FinalConstruct()
	{
		return S_OK;
	}

	void FinalRelease()
	{
	}

private:
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberVariationRegionData data;

public:
    affymetrix_calvin_data::ProbeSetMultiDataCopyNumberVariationRegionData &Data() { return data; }

public:
    STDMETHOD(get_name)(BSTR* pVal);
public:
    STDMETHOD(get_signal)(float* pVal);
public:
    STDMETHOD(get_call)(int* pVal);
public:
    STDMETHOD(get_confidenceScore)(float* pVal);
public:
    STDMETHOD(get_metricCount)(int* pVal);
public:
    STDMETHOD(GetMetric)(int index, VARIANT* pVal);
public:
    STDMETHOD(GetMetricName)(int index, BSTR* pVal);

};

OBJECT_ENTRY_AUTO(__uuidof(ProbeSetMultiDataCopyNumberVariationRegionData), CProbeSetMultiDataCopyNumberVariationRegionDataCOM)
