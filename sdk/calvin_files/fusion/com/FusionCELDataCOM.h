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
#include "calvin_files/fusion/src/FusionCELData.h"
#include "FusionCELFileEntryTypeCOM.h"


// CFusionCELDataCOM

class ATL_NO_VTABLE CFusionCELDataCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCELDataCOM, &CLSID_FusionCELData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCELData, &IID_IFusionCELData, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCELDataCOM()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCELDATA)


BEGIN_COM_MAP(CFusionCELDataCOM)
	COM_INTERFACE_ENTRY(IFusionCELData)
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
	affymetrix_fusion_io::FusionCELData cellData;

public:

public:
	STDMETHOD(get_FileName)(BSTR* pVal);
public:
	STDMETHOD(put_FileName)(BSTR newVal);
public:
	STDMETHOD(get_Error)(BSTR* pVal);
public:
	STDMETHOD(GetHeaderKey)(BSTR key, BSTR *pVal);
public:
	STDMETHOD(get_Version)(int* pVal);
public:
	STDMETHOD(get_Cols)(int* pVal);
public:
	STDMETHOD(get_Rows)(int* pVal);
public:
	STDMETHOD(get_NumCells)(int* pVal);
public:
	STDMETHOD(get_Header)(BSTR* pVal);
public:
	STDMETHOD(get_Alg)(BSTR* pVal);
public:
	STDMETHOD(get_Params)(BSTR* pVal);
public:
	STDMETHOD(GetAlgorithmParameter)(BSTR tag, BSTR * pVal);
public:
	STDMETHOD(GetAlgorithmParameterTag)(int index, BSTR *pVal);
public:
	STDMETHOD(get_NumberAlgorithmParameters)(int* pVal);
public:
	STDMETHOD(get_AlgorithmParameters)(BSTR* pVal);
public:
	STDMETHOD(get_DatHeader)(BSTR* pVal);
public:
	STDMETHOD(get_ChipType)(BSTR* pVal);
public:
	STDMETHOD(get_CellMargin)(int* pVal);
public:
	STDMETHOD(get_NumOutliers)(int* pVal);
public:
	STDMETHOD(get_NumMasked)(int* pVal);
public:
	STDMETHOD(IndexToX)(int index, int * x);
public:
	STDMETHOD(IndexToY)(int index, int * y);
public:
	STDMETHOD(XYToIndex)(int x, int y, int * index);
public:
	STDMETHOD(XYToIndexStatic)(int x, int y, int r, int c, int * index);
public:
	STDMETHOD(GetEntryByIndex)(int index, IFusionCELFileEntryType * entry);
public:
	STDMETHOD(GetEntryByCoordinate)(int x, int y, IFusionCELFileEntryType * entry);
public:
	STDMETHOD(GetIntensityByIndex)(int index, float * intensity);
public:
	STDMETHOD(GetIntensityByCoordinate)(int x, int y, float * intensity);
public:
	STDMETHOD(GetStdvByIndex)(int index, float * stdv);
public:
	STDMETHOD(GetStdvByCoordinate)(int x, int y, float * stdv);
public:
	STDMETHOD(GetPixelsByIndex)(int index, short * pixels);
public:
	STDMETHOD(GetPixelsByCoordinate)(int x, int y, short * pixels);
public:
	STDMETHOD(IsOutlierByIndex)(int index, VARIANT_BOOL * value);
public:
	STDMETHOD(IsOutlierByCoordinate)(int x, int y, VARIANT_BOOL * value);
public:
	STDMETHOD(IsMaskedByIndex)(int index, VARIANT_BOOL * value);
public:
	STDMETHOD(IsMaskedByCoordinate)(int x, int y, VARIANT_BOOL * value);
public:
	STDMETHOD(Close)(void);
public:
	STDMETHOD(GetFileSize)(LONG* size);
public:
	STDMETHOD(Exists)(VARIANT_BOOL* pVal);
public:
	STDMETHOD(ReadHeader)(VARIANT_BOOL* pVal);
public:
	STDMETHOD(Read)(VARIANT_BOOL includeMaskAndOutliers, VARIANT_BOOL* pVal);
public:
	STDMETHOD(Clear)(void);
public:
	STDMETHOD(GetParameters)(VARIANT* params);
public:
	STDMETHOD(GetGridCorners)(IFGridCoords ** grid);
public:
	STDMETHOD(get_FileId)(BSTR* pVal);
public:
    STDMETHOD(get_MasterFileName)(BSTR* pVal);
public:
    STDMETHOD(get_LibraryPackageName)(BSTR* pVal);
public:
	STDMETHOD(SetActiveDataGroup)(BSTR channel);
public:
	STDMETHOD(IsMultiColor)(VARIANT_BOOL *pVal);
public:
	STDMETHOD(GetChannels)(VARIANT *channels);
};

OBJECT_ENTRY_AUTO(__uuidof(FusionCELData), CFusionCELDataCOM)
