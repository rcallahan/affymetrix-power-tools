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


#pragma once
#include "resource.h"       // main symbols

#include "affx_fusion_com.h"
#include "file/CHPFileData.h"

// CBackgroundZoneInfoCOM

class ATL_NO_VTABLE CBackgroundZoneInfoCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CBackgroundZoneInfoCOM, &CLSID_BackgroundZoneInfo>,
	public ISupportErrorInfo,
	public IDispatchImpl<IBackgroundZoneInfo, &IID_IBackgroundZoneInfo, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CBackgroundZoneInfoCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_BACKGROUNDZONEINFO)


BEGIN_COM_MAP(CBackgroundZoneInfoCOM)
	COM_INTERFACE_ENTRY(IBackgroundZoneInfo)
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
	affxchp::BackgroundZoneInfo zone;

public:
	affxchp::BackgroundZoneInfo &Zone() { return zone; }

public:
	STDMETHOD(get_number_zones)(int* pVal);
public:
	STDMETHOD(get_smooth_factor)(float* pVal);
public:
	STDMETHOD(get_zones)(VARIANT* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(BackgroundZoneInfo), CBackgroundZoneInfoCOM)
