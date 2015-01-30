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
#include "FusionCDFDataCOM.h"
#include "FusionCDFFileHeaderCOM.h"
#include "FusionCDFProbeSetInformationCOM.h"
#include "FusionCDFQCProbeSetInformationCOM.h"
#include "COMStringUtils.h"

// CFusionCDFDataCOM

STDMETHODIMP CFusionCDFDataCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCDFData
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCDFDataCOM::get_FileName(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(cdf.GetFileName());
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::put_FileName(BSTR newVal)
{
	cdf.SetFileName(COMStringUtils::ConvertString(newVal).c_str());
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::GetHeader(IFusionCDFFileHeader ** pVal)
{
	CComPtr<IFusionCDFFileHeader> header;
	header.CoCreateInstance(CLSID_FusionCDFFileHeader);
	CFusionCDFFileHeaderCOM *pHeader = static_cast<CFusionCDFFileHeaderCOM*>(header.p);
	pHeader->SetHeader(&cdf.GetHeader());
	header.CopyTo(pVal);
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::get_Error(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(cdf.GetError());
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::GetProbeSetName(int index, BSTR * pVal)
{
	*pVal = COMStringUtils::ConvertString(cdf.GetProbeSetName(index));
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::get_ChipType(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(cdf.GetChipType());
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::Read(VARIANT_BOOL* pVal)
{
	*pVal = (cdf.Read() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::ReadHeader(VARIANT_BOOL* pVal)
{
	*pVal = (cdf.ReadHeader() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::Exists(VARIANT_BOOL* pVal)
{
	*pVal = (cdf.Exists() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::Close(void)
{
	cdf.Close();
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::IsXDACompatibleFile(BSTR fileName, VARIANT_BOOL* pVal)
{
	*pVal = (affymetrix_fusion_io::FusionCDFData::IsXDACompatibleFile(COMStringUtils::ConvertString(fileName).c_str()) == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::IsCalvinCompatibleFile(BSTR fileName, VARIANT_BOOL* pVal)
{
	*pVal = (affymetrix_fusion_io::FusionCDFData::IsCalvinCompatibleFile(COMStringUtils::ConvertString(fileName).c_str()) == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::GetProbeSetType(int index, GeneChipProbeSetType *pVal)
{
	*pVal = (GeneChipProbeSetType) cdf.GetProbeSetType(index);
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::GetProbeSetInformation(int index, IFusionCDFProbeSetInformation * pVal)
{
	CFusionCDFProbeSetInformationCOM *set = (CFusionCDFProbeSetInformationCOM *) pVal;
	cdf.GetProbeSetInformation(index, set->Set());
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::GetQCProbeSetInformation(int index, IFusionCDFQCProbeSetInformation * pVal)
{
	CFusionCDFQCProbeSetInformationCOM *set = (CFusionCDFQCProbeSetInformationCOM *) pVal;
	cdf.GetQCProbeSetInformation(index, set->Set());
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::GetQCProbeSetInformationByType(GeneChipQCProbeSetType qcType, IFusionCDFQCProbeSetInformation * pVal)
{
	CFusionCDFQCProbeSetInformationCOM *set = (CFusionCDFQCProbeSetInformationCOM *) pVal;
	cdf.GetQCProbeSetInformation((affxcdf::GeneChipQCProbeSetType) qcType, set->Set());
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::get_ChipTypes(VARIANT* pVal)
{
	VariantInit(pVal);
	pVal->vt = VT_ARRAY | VT_BSTR;
	pVal->parray = NULL;
	std::vector<std::string> ids = cdf.GetChipTypes();
	int nids = (int)ids.size();
	SAFEARRAYBOUND  rgsaBound[1];
	rgsaBound[0].lLbound = 0;
	rgsaBound[0].cElements = nids;
	pVal->parray = SafeArrayCreate(VT_BSTR, 1, rgsaBound);
	for (long index=0; index<(long)ids.size(); index++)
	{
		BSTR id = COMStringUtils::ConvertString(ids[index]);
		HRESULT hr = SafeArrayPutElement(pVal->parray, &index, id);
	}
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::get_GUID(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(cdf.GetGUID());
	return S_OK;
}

STDMETHODIMP CFusionCDFDataCOM::get_IntegrityMd5(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(cdf.GetIntegrityMd5());
	return S_OK;
}

