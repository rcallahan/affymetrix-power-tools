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
#include "FusionExpressionProbeSetResultsCOM.h"
#include "COMStringUtils.h"

// CFusionExpressionProbeSetResultsCOM

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFusionExpressionProbeSetResults
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_DetectionPValue(float* pVal)
{
	*pVal = results.GetDetectionPValue();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_Signal(float* pVal)
{
	*pVal = results.GetSignal();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_NumPairs(int* pVal)
{
	*pVal = results.GetNumPairs();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_NumUsedPairs(int* pVal)
{
	*pVal = results.GetNumUsedPairs();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_Detection(CHAR* pVal)
{
	*pVal = results.GetDetection();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_HasCompResults(VARIANT_BOOL* pVal)
{
	*pVal = (results.HasCompResults() == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_ChangePValue(float* pVal)
{
	*pVal = results.GetChangePValue();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_SignalLogRatio(float* pVal)
{
	*pVal = results.GetSignalLogRatio();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_SignalLogRatioLow(float* pVal)
{
	*pVal = results.GetSignalLogRatioLow();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_SignalLogRatioHigh(float* pVal)
{
	*pVal = results.GetSignalLogRatioHigh();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_NumCommonPairs(int* pVal)
{
	*pVal = results.GetNumCommonPairs();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::get_Change(CHAR* pVal)
{
	*pVal = results.GetChange();
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::GetDetectionString(BSTR * pVal)
{
	*pVal = COMStringUtils::ConvertString(results.GetDetectionString());
	return S_OK;
}

STDMETHODIMP CFusionExpressionProbeSetResultsCOM::GetChangeString(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(results.GetChangeString());
	return S_OK;
}
