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

// CFusionCDFProbeInformationCOM

class ATL_NO_VTABLE CFusionCDFProbeInformationCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCDFProbeInformationCOM, &CLSID_FusionCDFProbeInformation>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCDFProbeInformation, &IID_IFusionCDFProbeInformation, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCDFProbeInformationCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCDFPROBEINFORMATION)


BEGIN_COM_MAP(CFusionCDFProbeInformationCOM)
	COM_INTERFACE_ENTRY(IFusionCDFProbeInformation)
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
	affymetrix_fusion_io::FusionCDFProbeInformation probe;

public:
	affymetrix_fusion_io::FusionCDFProbeInformation &Probe() { return probe; }

public:

public:
	STDMETHOD(get_ListIndex)(int* pVal);
public:
	STDMETHOD(get_Expos)(int* pVal);
public:
	STDMETHOD(get_X)(int* pVal);
public:
	STDMETHOD(get_Y)(int* pVal);
public:
	STDMETHOD(get_PBase)(CHAR* pVal);
public:
	STDMETHOD(get_TBase)(CHAR* pVal);
public:
	STDMETHOD(get_ProbeLength)(USHORT* pVal);
public:
	STDMETHOD(get_ProbeGrouping)(USHORT* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCDFProbeInformation), CFusionCDFProbeInformationCOM)
