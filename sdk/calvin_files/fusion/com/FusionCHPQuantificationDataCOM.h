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
#include "calvin_files/fusion/src/FusionCHPQuantificationData.h"

// CFusionCHPQuantificationDataCOM

class ATL_NO_VTABLE CFusionCHPQuantificationDataCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCHPQuantificationDataCOM, &CLSID_FusionCHPQuantificationData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCHPQuantificationData, &IID_IFusionCHPQuantificationData, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCHPQuantificationDataCOM()
	{
		chp = NULL;
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCHPQUANTIFICATIONDATA)


BEGIN_COM_MAP(CFusionCHPQuantificationDataCOM)
	COM_INTERFACE_ENTRY(IFusionCHPQuantificationData)
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
		delete chp;
        chp = NULL;
	}

private:
	affymetrix_fusion_io::FusionCHPQuantificationData *chp;

public:
	STDMETHOD(get_FileTypeIdentifiers)(VARIANT* pVal);
public:
	STDMETHOD(FromBase)(IFusionCHPData * baseChp, VARIANT_BOOL* pVal);
public:
	STDMETHOD(get_AlgName)(BSTR* pVal);
public:
	STDMETHOD(get_AlgVersion)(BSTR* pVal);
public:
	STDMETHOD(get_ArrayType)(BSTR* pVal);
public:
	STDMETHOD(get_AlgorithmParameters)(VARIANT* pVal);
public:
	STDMETHOD(get_SummaryParameters)(VARIANT* pVal);
public:
	STDMETHOD(get_EntryCount)(int* pVal);
public:
	STDMETHOD(GetQuantificationEntry)(int index, IProbeSetQuantificationData * pVal);
public:
	STDMETHOD(get_FileId)(BSTR* pVal);
public:
    STDMETHOD(Close)();
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCHPQuantificationData), CFusionCHPQuantificationDataCOM)
