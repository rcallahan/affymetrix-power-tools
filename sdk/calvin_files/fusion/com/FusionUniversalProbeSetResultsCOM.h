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
#include "calvin_files/fusion/src/FusionProbeSetResults.h"


// CFusionUniversalProbeSetResultsCOM

class ATL_NO_VTABLE CFusionUniversalProbeSetResultsCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionUniversalProbeSetResultsCOM, &CLSID_FusionUniversalProbeSetResults>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionUniversalProbeSetResults, &IID_IFusionUniversalProbeSetResults, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionUniversalProbeSetResultsCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONUNIVERSALPROBESETRESULTS)


BEGIN_COM_MAP(CFusionUniversalProbeSetResultsCOM)
	COM_INTERFACE_ENTRY(IFusionUniversalProbeSetResults)
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
	affymetrix_fusion_io::FusionUniversalProbeSetResults results;

public:
	affymetrix_fusion_io::FusionUniversalProbeSetResults &Results() { return results; }
public:

public:
	STDMETHOD(get_Background)(float* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionUniversalProbeSetResults), CFusionUniversalProbeSetResultsCOM)
