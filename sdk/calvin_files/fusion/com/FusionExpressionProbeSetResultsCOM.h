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

// CFusionExpressionProbeSetResultsCOM

class ATL_NO_VTABLE CFusionExpressionProbeSetResultsCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionExpressionProbeSetResultsCOM, &CLSID_FusionExpressionProbeSetResults>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionExpressionProbeSetResults, &IID_IFusionExpressionProbeSetResults, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionExpressionProbeSetResultsCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONEXPRESSIONPROBESETRESULTS)


BEGIN_COM_MAP(CFusionExpressionProbeSetResultsCOM)
	COM_INTERFACE_ENTRY(IFusionExpressionProbeSetResults)
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
	affymetrix_fusion_io::FusionExpressionProbeSetResults results;

public:
	affymetrix_fusion_io::FusionExpressionProbeSetResults &Results() { return results; }


public:

public:
	STDMETHOD(get_DetectionPValue)(float* pVal);
public:
	STDMETHOD(get_Signal)(float* pVal);
public:
	STDMETHOD(get_NumPairs)(int* pVal);
public:
	STDMETHOD(get_NumUsedPairs)(int* pVal);
public:
	STDMETHOD(get_Detection)(CHAR* pVal);
public:
	STDMETHOD(get_HasCompResults)(VARIANT_BOOL* pVal);
public:
	STDMETHOD(get_ChangePValue)(float* pVal);
public:
	STDMETHOD(get_SignalLogRatio)(float* pVal);
public:
	STDMETHOD(get_SignalLogRatioLow)(float* pVal);
public:
	STDMETHOD(get_SignalLogRatioHigh)(float* pVal);
public:
	STDMETHOD(get_NumCommonPairs)(int* pVal);
public:
	STDMETHOD(get_Change)(CHAR* pVal);
public:
	STDMETHOD(GetDetectionString)(BSTR * pVal);
public:
	STDMETHOD(GetChangeString)(BSTR* pVal);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionExpressionProbeSetResults), CFusionExpressionProbeSetResultsCOM)
