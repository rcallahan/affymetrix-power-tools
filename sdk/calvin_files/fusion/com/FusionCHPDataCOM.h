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
#include "calvin_files/fusion/src/FusionCHPData.h"

// CFusionCHPDataCOM

class ATL_NO_VTABLE CFusionCHPDataCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCHPDataCOM, &CLSID_FusionCHPData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCHPData, &IID_IFusionCHPData, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCHPDataCOM()
	{
		chp = NULL;
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCHPDATA1)


BEGIN_COM_MAP(CFusionCHPDataCOM)
	COM_INTERFACE_ENTRY(IFusionCHPData)
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
	affymetrix_fusion_io::FusionCHPData *chp;

public:
	void SetChip(affymetrix_fusion_io::FusionCHPData *c) { chp = c; }

public:
	affymetrix_fusion_io::FusionCHPData *GetChip() { return chp; }

public:

public:
	STDMETHOD(get_FileTypeIdentifier)(BSTR* pVal);
public:
	STDMETHOD(get_FileTypeIdentifiers)(VARIANT* pVal);
public:
	STDMETHOD(get_FileId)(BSTR* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCHPData), CFusionCHPDataCOM)
