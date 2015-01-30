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


// CFusionGenotypeProbeSetResultsCOM

class ATL_NO_VTABLE CFusionGenotypeProbeSetResultsCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionGenotypeProbeSetResultsCOM, &CLSID_FusionGenotypeProbeSetResults>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionGenotypeProbeSetResults, &IID_IFusionGenotypeProbeSetResults, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionGenotypeProbeSetResultsCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONGENOTYPEPROBESETRESULTS)


BEGIN_COM_MAP(CFusionGenotypeProbeSetResultsCOM)
	COM_INTERFACE_ENTRY(IFusionGenotypeProbeSetResults)
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
	affymetrix_fusion_io::FusionGenotypeProbeSetResults results;

public:
	affymetrix_fusion_io::FusionGenotypeProbeSetResults &Results() { return results; }
public:

public:
	STDMETHOD(get_AlleleCall)(CHAR* pVal);
public:
	STDMETHOD(get_Confidence)(float* pVal);
public:
	STDMETHOD(get_RAS1)(float* pVal);
public:
	STDMETHOD(get_RAS2)(float* pVal);
public:
	STDMETHOD(get_PValueAA)(float* pVal);
public:
	STDMETHOD(get_PValueAB)(float* pVal);
public:
	STDMETHOD(get_PValueBB)(float* pVal);
public:
	STDMETHOD(get_PValueNoCall)(float* pVal);
public:
	STDMETHOD(GetAlleleCallString)(BSTR* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionGenotypeProbeSetResults), CFusionGenotypeProbeSetResultsCOM)
