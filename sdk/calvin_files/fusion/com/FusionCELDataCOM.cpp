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


#include "stdafx.h"
#include "FusionCELDataCOM.h"
#include "COMStringUtils.h"
#include "FusionTagValuePairTypeCOM.h"
#include "FGridCoordsCOM.h"

// CFusionCELData

STDMETHODIMP CFusionCELDataCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCELData
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCELDataCOM::get_FileName(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(cellData.GetFileName());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::put_FileName(BSTR newVal)
{
	cellData.SetFileName(COMStringUtils::ConvertString(newVal).c_str());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_Error(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetError());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetHeaderKey(BSTR key, BSTR *pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetHeaderKey(COMStringUtils::ConvertWString(key).c_str()));
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_Version(int* pVal)
{
	*pVal = cellData.GetVersion();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_Cols(int* pVal)
{
	*pVal = cellData.GetCols();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_Rows(int* pVal)
{
	*pVal = cellData.GetRows();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_NumCells(int* pVal)
{
	*pVal = cellData.GetNumCells();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_Header(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetHeader());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_Alg(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetAlg());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_Params(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetParams());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetAlgorithmParameter(BSTR tag, BSTR * pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetAlgorithmParameter(COMStringUtils::ConvertWString(tag).c_str()));
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetAlgorithmParameterTag(int index, BSTR *pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetAlgorithmParameterTag(index));
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_NumberAlgorithmParameters(int* pVal)
{
	*pVal = cellData.GetNumberAlgorithmParameters();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_AlgorithmParameters(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetAlgorithmParameters());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_DatHeader(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetDatHeader());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_ChipType(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetChipType());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_CellMargin(int* pVal)
{
	*pVal = cellData.GetCellMargin();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_NumOutliers(int* pVal)
{
	*pVal = cellData.GetNumOutliers();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_NumMasked(int* pVal)
{
	*pVal = cellData.GetNumMasked();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::IndexToX(int index, int * x)
{
	*x = cellData.IndexToX(index);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::IndexToY(int index, int * y)
{
	*y = cellData.IndexToY(index);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::XYToIndex(int x, int y, int * index)
{
	*index = cellData.XYToIndex(x, y);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::XYToIndexStatic(int x, int y, int r, int c, int * index)
{
	*index = affymetrix_fusion_io::FusionCELData::XYToIndex(x, y, r, c);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetEntryByIndex(int index, IFusionCELFileEntryType * entry)
{
	CFusionCELFileEntryTypeCOM *e = static_cast<CFusionCELFileEntryTypeCOM*>(entry);
	cellData.GetEntry(index, e->Entry());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetEntryByCoordinate(int x, int y, IFusionCELFileEntryType * entry)
{
	CFusionCELFileEntryTypeCOM *e = static_cast<CFusionCELFileEntryTypeCOM*>(entry);
	cellData.GetEntry(x, y, e->Entry());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetIntensityByIndex(int index, float * intensity)
{
	*intensity = cellData.GetIntensity(index);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetIntensityByCoordinate(int x, int y, float * intensity)
{
	*intensity = cellData.GetIntensity(x,y);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetStdvByIndex(int index, float * stdv)
{
	*stdv = cellData.GetStdv(index);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetStdvByCoordinate(int x, int y, float * stdv)
{
	*stdv = cellData.GetStdv(x,y);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetPixelsByIndex(int index, short * pixels)
{
	*pixels = cellData.GetPixels(index);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetPixelsByCoordinate(int x, int y, short * pixels)
{
	*pixels = cellData.GetPixels(x,y);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::IsOutlierByIndex(int index, VARIANT_BOOL * value)
{
	*value = (cellData.IsOutlier(index) == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::IsOutlierByCoordinate(int x, int y, VARIANT_BOOL * value)
{
	*value = (cellData.IsOutlier(x,y) == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::IsMaskedByIndex(int index, VARIANT_BOOL * value)
{
	*value = (cellData.IsMasked(index) == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::IsMaskedByCoordinate(int x, int y, VARIANT_BOOL * value)
{
	*value = (cellData.IsMasked(x,y) == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::Close(void)
{
	cellData.Close();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetFileSize(LONG* size)
{
	*size = (LONG) cellData.GetFileSize();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::Exists(VARIANT_BOOL *pVal)
{
	*pVal = (cellData.Exists() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::ReadHeader(VARIANT_BOOL* pVal)
{
	*pVal = (cellData.ReadHeader() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::Read(VARIANT_BOOL includeMaskAndOutliers, VARIANT_BOOL* pVal)
{
	*pVal = (cellData.Read((includeMaskAndOutliers == VARIANT_TRUE ? true : false)) == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::Clear(void)
{
	cellData.Clear();
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetParameters(VARIANT* params)
{
	VariantInit(params);
	params->vt = VT_ARRAY | VT_DISPATCH;
	params->parray = NULL;
	affymetrix_fusion_io::FusionTagValuePairTypeList &cparams = cellData.GetParameters();
	int nparams = (int)cparams.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nparams;
	params->parray = SafeArrayCreate(VT_DISPATCH, 1, rgsaBound);
	long index=0;
	affymetrix_fusion_io::FusionTagValuePairTypeList::iterator it;
	for (it=cparams.begin(); it!=cparams.end(); ++it)
	{
		CComPtr<IFusionTagValuePairType> param;
		param.CoCreateInstance(CLSID_FusionTagValuePairType);
		CFusionTagValuePairTypeCOM *pParam = static_cast<CFusionTagValuePairTypeCOM *>(param.p);
		pParam->SetParam((*it));
		HRESULT hr = SafeArrayPutElement(params->parray, &index, param);
		++index;
	}
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetGridCorners(IFGridCoords ** pVal)
{
	CComPtr<IFGridCoords> grid;
	grid.CoCreateInstance(CLSID_FGridCoords);
	CFGridCoordsCOM *pGrid = static_cast<CFGridCoordsCOM*>(grid.p);
	pGrid->SetGrid(cellData.GetGridCorners());
	grid.CopyTo(pVal);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_FileId(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(cellData.GetFileId());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_MasterFileName(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetMasterFileName());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::get_LibraryPackageName(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertWString(cellData.GetLibraryPackageName());
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::SetActiveDataGroup(BSTR channel)
{
	cellData.SetActiveDataGroup(COMStringUtils::ConvertWString(channel));
	return S_OK;
}
STDMETHODIMP CFusionCELDataCOM::IsMultiColor(VARIANT_BOOL *pVal)
{
	*pVal = (cellData.IsMultiColor() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCELDataCOM::GetChannels(VARIANT *channels)
{
	VariantInit(channels);
	channels->vt = VT_ARRAY | VT_BSTR;
	channels->parray = NULL;
	WStringVector ch = cellData.GetChannels();
	int n = (int)ch.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = n;
	channels->parray = SafeArrayCreate(VT_BSTR, 1, rgsaBound);
	long index=0;
	for (WStringVector::iterator it=ch.begin(); it!=ch.end(); ++it)
	{
		BSTR c = COMStringUtils::ConvertWString(*it);
		SafeArrayPutElement(channels->parray, &index, c);
		++index;
	}
	return S_OK;
}
