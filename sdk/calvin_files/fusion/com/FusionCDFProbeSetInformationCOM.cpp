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
#include "FusionCDFProbeSetInformationCOM.h"
#include "FusionCDFProbeGroupInformationCOM.h"


// CFusionCDFProbeSetInformationCOM

STDMETHODIMP CFusionCDFProbeSetInformationCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionCDFProbeSetInformation
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionCDFProbeSetInformationCOM::get_ProbeSetType(GeneChipProbeSetType* pVal)
{
	*pVal = (GeneChipProbeSetType) set.GetProbeSetType();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeSetInformationCOM::get_Direction(DirectionType* pVal)
{
	*pVal = (DirectionType) set.GetDirection();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeSetInformationCOM::get_NumLists(int* pVal)
{
	*pVal = set.GetNumLists();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeSetInformationCOM::get_NumGroups(int* pVal)
{
	*pVal = set.GetNumGroups();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeSetInformationCOM::get_NumCells(int* pVal)
{
	*pVal = set.GetNumCells();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeSetInformationCOM::get_NumCellsPerList(int* pVal)
{
	*pVal = set.GetNumCellsPerList();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeSetInformationCOM::get_ProbeSetNumber(int* pVal)
{
	*pVal = set.GetProbeSetNumber();
	return S_OK;
}

STDMETHODIMP CFusionCDFProbeSetInformationCOM::GetGroup(int index, IFusionCDFProbeGroupInformation * group)
{
	CFusionCDFProbeGroupInformationCOM *g = static_cast<CFusionCDFProbeGroupInformationCOM*>(group);
	set.GetGroupInformation(index, g->Group());
	return S_OK;
}
