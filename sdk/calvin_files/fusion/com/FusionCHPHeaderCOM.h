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
#include "calvin_files/fusion/src/FusionCHPLegacyData.h"
#include <vector>
#include <string>


// CFusionCHPHeaderCOM

class ATL_NO_VTABLE CFusionCHPHeaderCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCHPHeaderCOM, &CLSID_FusionCHPHeader>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCHPHeader, &IID_IFusionCHPHeader, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCHPHeaderCOM()
	{
		header = NULL;
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCHPHEADER)


BEGIN_COM_MAP(CFusionCHPHeaderCOM)
	COM_INTERFACE_ENTRY(IFusionCHPHeader)
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
	affymetrix_fusion_io::FusionCHPHeader *header;
	std::vector<std::wstring> sumParamNames;
	std::vector<std::wstring> algParamNames;

public:
	void SetHeader(affymetrix_fusion_io::FusionCHPHeader *h) { header = h; sumParamNames.clear(); algParamNames.clear(); }

public:

public:
	STDMETHOD(get_Cols)(int* pVal);
public:
	STDMETHOD(get_Rows)(int* pVal);
public:
	STDMETHOD(get_NumProbeSets)(int* pVal);
public:
	STDMETHOD(get_AssayType)(AssayType* pVal);
public:
	STDMETHOD(get_ChipType)(BSTR* pVal);
public:
	STDMETHOD(get_AlgName)(BSTR* pVal);
public:
	STDMETHOD(get_AlgVersion)(BSTR* pVal);
public:
	STDMETHOD(get_AlgorithmParameters)(VARIANT* pVal);
public:
	STDMETHOD(get_AlgorithmParameterCount)(int* pVal);
public:
	STDMETHOD(get_SummaryParameters)(VARIANT* pVal);
public:
	STDMETHOD(get_SummaryParameterCount)(int* pVal);
public:
	STDMETHOD(get_ParentCellFile)(BSTR* pVal);
public:
	STDMETHOD(get_ProgID)(BSTR* pVal);
public:
	STDMETHOD(GetAlgorithmParameter)(BSTR tag, BSTR * pVal);
public:
	STDMETHOD(GetSummaryParameter)(BSTR tag, BSTR* pVal);
public:
	STDMETHOD(GetAlgorithmParameterName)(int index, BSTR * pVal);
public:
	STDMETHOD(GetSummaryParameterName)(int index, BSTR* pVal);
public:
	STDMETHOD(Clear)(void);
public:
	STDMETHOD(get_BackgroundZoneInfo)(IBackgroundZoneInfo ** pVal);
public:
	STDMETHOD(GetBackgroundZones)(VARIANT * zones);
public:
	STDMETHOD(GetBackgroundZone)(int x, int y, IBackgroundZoneType ** pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCHPHeader), CFusionCHPHeaderCOM)
