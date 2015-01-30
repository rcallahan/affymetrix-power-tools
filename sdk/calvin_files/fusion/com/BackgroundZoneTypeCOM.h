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

// CBackgroundZoneTypeCOM

class ATL_NO_VTABLE CBackgroundZoneTypeCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CBackgroundZoneTypeCOM, &CLSID_BackgroundZoneType>,
	public ISupportErrorInfo,
	public IDispatchImpl<IBackgroundZoneType, &IID_IBackgroundZoneType, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CBackgroundZoneTypeCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_BACKGROUNDZONETYPE)


BEGIN_COM_MAP(CBackgroundZoneTypeCOM)
	COM_INTERFACE_ENTRY(IBackgroundZoneType)
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
	affxchp::BackgroundZoneType bg;

public:
	void SetBg(affxchp::BackgroundZoneType b) { bg = b; }

public:
	STDMETHOD(get_centerx)(float* pVal);
public:
	STDMETHOD(get_centery)(float* pVal);
public:
	STDMETHOD(get_background)(float* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(BackgroundZoneType), CBackgroundZoneTypeCOM)
