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


// CFusionCDFProbeSetInformationCOM

class ATL_NO_VTABLE CFusionCDFProbeSetInformationCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCDFProbeSetInformationCOM, &CLSID_FusionCDFProbeSetInformation>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCDFProbeSetInformation, &IID_IFusionCDFProbeSetInformation, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCDFProbeSetInformationCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCDFPROBESETINFORMATION)


BEGIN_COM_MAP(CFusionCDFProbeSetInformationCOM)
	COM_INTERFACE_ENTRY(IFusionCDFProbeSetInformation)
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
	affymetrix_fusion_io::FusionCDFProbeSetInformation set;

public:
	affymetrix_fusion_io::FusionCDFProbeSetInformation &Set() { return set; }

public:

public:
	STDMETHOD(get_ProbeSetType)(GeneChipProbeSetType* pVal);
public:
	STDMETHOD(get_Direction)(DirectionType* pVal);
public:
	STDMETHOD(get_NumLists)(int* pVal);
public:
	STDMETHOD(get_NumGroups)(int* pVal);
public:
	STDMETHOD(get_NumCells)(int* pVal);
public:
	STDMETHOD(get_NumCellsPerList)(int* pVal);
public:
	STDMETHOD(get_ProbeSetNumber)(int* pVal);
public:
	STDMETHOD(GetGroup)(int index, IFusionCDFProbeGroupInformation * group);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCDFProbeSetInformation), CFusionCDFProbeSetInformationCOM)
