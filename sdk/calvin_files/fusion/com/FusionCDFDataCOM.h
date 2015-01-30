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
#include "calvin_files/fusion/src/FusionCDFData.h"


// CFusionCDFDataCOM

class ATL_NO_VTABLE CFusionCDFDataCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCDFDataCOM, &CLSID_FusionCDFData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCDFData, &IID_IFusionCDFData, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCDFDataCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCDFDATA)


BEGIN_COM_MAP(CFusionCDFDataCOM)
	COM_INTERFACE_ENTRY(IFusionCDFData)
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
	affymetrix_fusion_io::FusionCDFData cdf;

public:

public:
	STDMETHOD(get_FileName)(BSTR* pVal);
public:
	STDMETHOD(put_FileName)(BSTR newVal);
public:
	STDMETHOD(GetHeader)(IFusionCDFFileHeader ** pVal);
public:
	STDMETHOD(get_Error)(BSTR* pVal);
public:
	STDMETHOD(GetProbeSetName)(int index, BSTR * pVal);
public:
	STDMETHOD(get_ChipType)(BSTR* pVal);
public:
	STDMETHOD(Read)(VARIANT_BOOL* pVal);
public:
	STDMETHOD(ReadHeader)(VARIANT_BOOL* pVal);
public:
	STDMETHOD(Exists)(VARIANT_BOOL* pVal);
public:
	STDMETHOD(Close)(void);
public:
	STDMETHOD(IsXDACompatibleFile)(BSTR fileName, VARIANT_BOOL* pVal);
public:
	STDMETHOD(IsCalvinCompatibleFile)(BSTR fileName, VARIANT_BOOL* pVal);
public:
	STDMETHOD(GetProbeSetType)(int index, GeneChipProbeSetType * pVal);
public:
	STDMETHOD(GetProbeSetInformation)(int index, IFusionCDFProbeSetInformation * pVal);
public:
	STDMETHOD(GetQCProbeSetInformation)(int index, IFusionCDFQCProbeSetInformation * pVal);
public:
	STDMETHOD(GetQCProbeSetInformationByType)(GeneChipQCProbeSetType qcType, IFusionCDFQCProbeSetInformation * pVal);
public:
	STDMETHOD(get_ChipTypes)(VARIANT* pVal);
public:
	STDMETHOD(get_GUID)(BSTR* pVal);
public:
	STDMETHOD(get_IntegrityMd5)(BSTR* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCDFData), CFusionCDFDataCOM)
