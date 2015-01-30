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
#include "FusionCDFFileHeaderCOM.h"
#include "COMStringUtils.h"

// CFusionCDFFileHeaderCOM

STDMETHODIMP CFusionCDFFileHeaderCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCDFFileHeader
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCDFFileHeaderCOM::get_Cols(int* pVal)
{
	*pVal = header->GetCols();
	return S_OK;
}

STDMETHODIMP CFusionCDFFileHeaderCOM::get_Rows(int* pVal)
{
	*pVal = header->GetRows();
	return S_OK;
}

STDMETHODIMP CFusionCDFFileHeaderCOM::get_NumProbeSets(int* pVal)
{
	*pVal = header->GetNumProbeSets();
	return S_OK;
}

STDMETHODIMP CFusionCDFFileHeaderCOM::get_NumQCProbeSets(int* pVal)
{
	*pVal = header->GetNumQCProbeSets();
	return S_OK;
}

STDMETHODIMP CFusionCDFFileHeaderCOM::get_Reference(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(header->GetReference());
	return S_OK;
}

STDMETHODIMP CFusionCDFFileHeaderCOM::get_FormatVersion(int* pVal)
{
	*pVal = header->GetFormatVersion();
	return S_OK;
}

STDMETHODIMP CFusionCDFFileHeaderCOM::get_GUID(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(header->GetGUID());
	return S_OK;
}

STDMETHODIMP CFusionCDFFileHeaderCOM::get_IntegrityMd5(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(header->GetIntegrityMd5());
	return S_OK;
}
