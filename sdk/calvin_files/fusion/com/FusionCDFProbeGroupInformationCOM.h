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


// CFusionCDFProbeGroupInformationCOM

class ATL_NO_VTABLE CFusionCDFProbeGroupInformationCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCDFProbeGroupInformationCOM, &CLSID_FusionCDFProbeGroupInformation>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCDFProbeGroupInformation, &IID_IFusionCDFProbeGroupInformation, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCDFProbeGroupInformationCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCDFPROBEGROUPINFORMATION)


BEGIN_COM_MAP(CFusionCDFProbeGroupInformationCOM)
	COM_INTERFACE_ENTRY(IFusionCDFProbeGroupInformation)
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
	affymetrix_fusion_io::FusionCDFProbeGroupInformation group;

public:
	affymetrix_fusion_io::FusionCDFProbeGroupInformation &Group() { return group; }

public:

public:
	STDMETHOD(get_Direction)(DirectionType* pVal);
public:
	STDMETHOD(get_NumLists)(int* pVal);
public:
	STDMETHOD(get_NumCells)(int* pVal);
public:
	STDMETHOD(get_NumCellsPerList)(int* pVal);
public:
	STDMETHOD(get_Start)(int* pVal);
public:
	STDMETHOD(get_Stop)(int* pVal);
public:
	STDMETHOD(get_Name)(BSTR* pVal);
public:
	STDMETHOD(GetCell)(int index, IFusionCDFProbeInformation * info);
public:
	STDMETHOD(get_WobbleSituation)(USHORT *pVal);
public:
	STDMETHOD(get_AlleleCode)(USHORT *pVal);
public:
	STDMETHOD(get_Channel)(int *pVal);
public:
	STDMETHOD(get_RepType)(ReplicationType *pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCDFProbeGroupInformation), CFusionCDFProbeGroupInformationCOM)
