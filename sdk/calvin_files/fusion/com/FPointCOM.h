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
#include "calvin_files/fusion/src/FusionCoords.h"


// CFPointCOM

class ATL_NO_VTABLE CFPointCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFPointCOM, &CLSID_FPoint>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFPoint, &IID_IFPoint, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFPointCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FPOINT)


BEGIN_COM_MAP(CFPointCOM)
	COM_INTERFACE_ENTRY(IFPoint)
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
	affymetrix_fusion_io::FPoint pt;

public:
	void SetPoint(affymetrix_fusion_io::FPoint &p) { pt = p; }

public:
	STDMETHOD(get_x)(float* pVal);
public:
	STDMETHOD(get_y)(float* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FPoint), CFPointCOM)
