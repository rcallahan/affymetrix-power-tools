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
#include "FusionCHPDataRegCOM.h"
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "FusionCHPDataCOM.h"
#include "COMStringUtils.h"
#include <string>

using namespace std;

// CFusionCHPDataRegCOM

STDMETHODIMP CFusionCHPDataRegCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCHPDataReg
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCHPDataRegCOM::Read(BSTR fileName, IFusionCHPData ** pVal)
{
	affymetrix_fusion_io::FusionCHPData *chpData =
		affymetrix_fusion_io::FusionCHPDataReg::Read(COMStringUtils::ConvertString(fileName));
	if (chpData == NULL)
		return E_FAIL;

	CComPtr<IFusionCHPData> chp;
	chp.CoCreateInstance(CLSID_FusionCHPData);
	CFusionCHPDataCOM *pchp = static_cast<CFusionCHPDataCOM *>(chp.p);
	pchp->SetChip(chpData);
	chp.CopyTo(pVal);

	return S_OK;
}

STDMETHODIMP CFusionCHPDataRegCOM::ReadHeader(BSTR fileName, IFusionCHPData ** pVal)
{
	affymetrix_fusion_io::FusionCHPData *chpData =
		affymetrix_fusion_io::FusionCHPDataReg::ReadHeader(COMStringUtils::ConvertString(fileName));
	if (chpData == NULL)
		return E_FAIL;

	CComPtr<IFusionCHPData> chp;
	chp.CoCreateInstance(CLSID_FusionCHPData);
	CFusionCHPDataCOM *pchp = static_cast<CFusionCHPDataCOM *>(chp.p);
	pchp->SetChip(chpData);
	chp.CopyTo(pVal);

	return S_OK;
}
