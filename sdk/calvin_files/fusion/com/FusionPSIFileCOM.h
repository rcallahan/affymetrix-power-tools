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
#include "file/PSIFileData.h"

// CFusionPSIFileCOM

class ATL_NO_VTABLE CFusionPSIFileCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionPSIFileCOM, &CLSID_FusionPSIFile>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionPSIFile, &IID_IFusionPSIFile, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionPSIFileCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONPSIFILE)


BEGIN_COM_MAP(CFusionPSIFileCOM)
	COM_INTERFACE_ENTRY(IFusionPSIFile)
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
	affxpsi::CPSIFileData psi;

public:

public:
	STDMETHOD(get_FileName)(BSTR* pVal);
public:
	STDMETHOD(put_FileName)(BSTR newVal);
public:
	STDMETHOD(get_ProbeSetCount)(int* pVal);
public:
	STDMETHOD(GetProbeSetName)(int index, BSTR* pVal);
public:
	STDMETHOD(GetProbePairs)(int index, int * pVal);
public:
	STDMETHOD(Read)(VARIANT_BOOL* pVal);
public:
	STDMETHOD(Exists)(VARIANT_BOOL* pVal);
public:
	STDMETHOD(Clear)();
};

OBJECT_ENTRY_AUTO(__uuidof(FusionPSIFile), CFusionPSIFileCOM)
