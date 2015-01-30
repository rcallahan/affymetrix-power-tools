////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

// MAS5.h : Declaration of the CMAS5

#pragma once
#include "resource.h"       // main symbols

#include "affx_mas5_com.h"
#include "MAS5Workflow.h"

// CMAS5

class ATL_NO_VTABLE CMAS5 :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CMAS5, &CLSID_MAS5>,
	public ISupportErrorInfo,
	public IDispatchImpl<IMAS5, &IID_IMAS5, &LIBID_affx_mas5_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CMAS5()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_MAS5)


BEGIN_COM_MAP(CMAS5)
	COM_INTERFACE_ENTRY(IMAS5)
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
	MAS5Workflow mas5;

public:
	STDMETHOD(RunAbsoluteAnalysis)(BSTR celFile, BSTR chpFile, BSTR cdfFile, VARIANT_BOOL* pVal);
public:
	STDMETHOD(RunComparisonAnalysis)(BSTR celFile, BSTR baselineFile, BSTR chpFile, BSTR cdfFile, VARIANT_BOOL* pVal);
public:
	STDMETHOD(get_Error)(BSTR* pVal);
public:
	STDMETHOD(put_SaveToLegacyFile)(VARIANT_BOOL newVal);
public:
    STDMETHOD(put_AlgorithmParameterFile)(BSTR newVal);
public:
    STDMETHOD(put_ReportControlsParameterFile)(BSTR newVal);
public:
    STDMETHOD(put_ProgramName)(BSTR newVal);
public:
    STDMETHOD(put_ProgramCompany)(BSTR newVal);
public:
    STDMETHOD(put_ProgramId)(BSTR newVal);
};

OBJECT_ENTRY_AUTO(__uuidof(MAS5), CMAS5)
