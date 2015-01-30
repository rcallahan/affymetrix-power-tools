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


#pragma once
#include "resource.h"       // main symbols

#include "affx_fusion_com.h"
#include "calvin_files/data/src/ProbeSetQuantificationDetectionData.h"

// CProbeSetQuantificationDetectionDataCOM

class ATL_NO_VTABLE CProbeSetQuantificationDetectionDataCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CProbeSetQuantificationDetectionDataCOM, &CLSID_ProbeSetQuantificationDetectionData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IProbeSetQuantificationDetectionData, &IID_IProbeSetQuantificationDetectionData, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CProbeSetQuantificationDetectionDataCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_PROBESETQUANTIFICATIONDETECTIONDATA)


BEGIN_COM_MAP(CProbeSetQuantificationDetectionDataCOM)
	COM_INTERFACE_ENTRY(IProbeSetQuantificationDetectionData)
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
	affymetrix_calvin_data::ProbeSetQuantificationDetectionData data;

public:
	affymetrix_calvin_data::ProbeSetQuantificationDetectionData &Data() { return data; }

public:
	STDMETHOD(get_name)(BSTR* pVal);
public:
	STDMETHOD(get_quantification)(float* pVal);
public:
    STDMETHOD(get_pvalue)(float* pVal);
public:
    STDMETHOD(get_identifier)(int* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(ProbeSetQuantificationDetectionData), CProbeSetQuantificationDetectionDataCOM)
