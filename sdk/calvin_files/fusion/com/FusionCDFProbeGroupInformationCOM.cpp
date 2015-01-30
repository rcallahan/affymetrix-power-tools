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
#include "FusionCDFProbeGroupInformationCOM.h"
#include "FusionCDFProbeInformationCOM.h"
#include "COMStringUtils.h"

// CFusionCDFProbeGroupInformationCOM

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCDFProbeGroupInformation
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_Direction(DirectionType* pVal)
{
	*pVal = (DirectionType)group.GetDirection();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_NumLists(int* pVal)
{
	*pVal = group.GetNumLists();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_NumCells(int* pVal)
{
	*pVal = group.GetNumCells();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_NumCellsPerList(int* pVal)
{
	*pVal = group.GetNumCellsPerList();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_Start(int* pVal)
{
	*pVal = group.GetStart();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_Stop(int* pVal)
{
	*pVal = group.GetStop();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_Name(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(group.GetName());
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::GetCell(int index, IFusionCDFProbeInformation * info)
{
	CFusionCDFProbeInformationCOM *c = static_cast<CFusionCDFProbeInformationCOM*>(info);
	group.GetCell(index, c->Probe());
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_WobbleSituation(USHORT* pVal)
{
	*pVal = group.GetWobbleSituation();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_AlleleCode(USHORT* pVal)
{
	*pVal = group.GetAlleleCode();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_Channel(int* pVal)
{
	*pVal = (int)group.GetChannel();
	return S_OK;
}
STDMETHODIMP CFusionCDFProbeGroupInformationCOM::get_RepType(ReplicationType *pVal)
{
	*pVal = (ReplicationType) group.GetRepType();
	return S_OK;
}
