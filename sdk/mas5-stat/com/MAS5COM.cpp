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

// MAS5.cpp : Implementation of CMAS5

#include "mas5-stat/com/MAS5COM.h"
//
#include "mas5-stat/com/COMStringUtils.h"
#include "mas5-stat/com/stdafx.h"
#include "mas5-stat/src/ExpressionAlgorithmImplementation.h"
#include "mas5-stat/workflow/MAS5ParameterExtraction.h"
//
#include "exp_report/src/ExpressionControlsParameterExtraction.h"
//


// CMAS5

STDMETHODIMP CMAS5::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IMAS5
	};

	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

/*
 * Run the algorithm and report and save the results to a CHP file.
 */
STDMETHODIMP CMAS5::RunAbsoluteAnalysis(BSTR celFile, BSTR chpFile, BSTR cdfFile, VARIANT_BOOL* pVal)
{
	std::string scelFile = COMStringUtils::ConvertString(celFile);
	std::string schpFile = COMStringUtils::ConvertString(chpFile);
	std::string scdfFile = COMStringUtils::ConvertString(cdfFile);
	std::string sbaselineFile;

	// Run the algorithm
	bool success = mas5.Execute(scelFile, sbaselineFile, scdfFile, schpFile);
	*pVal = (success == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

/*
 * Run the algorithm and report and save the results to a CHP file.
 */
STDMETHODIMP CMAS5::RunComparisonAnalysis(BSTR celFile, BSTR baselineFile, BSTR chpFile, BSTR cdfFile, VARIANT_BOOL* pVal)
{
	std::string scdfFile = COMStringUtils::ConvertString(cdfFile);
	std::string scelFile = COMStringUtils::ConvertString(celFile);
	std::string schpFile = COMStringUtils::ConvertString(chpFile);
	std::string sbaselineFile = COMStringUtils::ConvertString(baselineFile);

	// Run the algorithm
	bool success = mas5.Execute(scelFile, sbaselineFile, scdfFile, schpFile);
	*pVal = (success == true ? VARIANT_TRUE : VARIANT_FALSE);
	return S_OK;
}

/*
 * Get the error string.
 */
STDMETHODIMP CMAS5::get_Error(BSTR* pVal)
{
	*pVal = COMStringUtils::ConvertString(mas5.GetError());
	return S_OK;
}

STDMETHODIMP CMAS5::put_SaveToLegacyFile(VARIANT_BOOL newVal)
{
	mas5.SaveToLegacyFile() = ( newVal == VARIANT_TRUE ? true : false);
	return S_OK;
}

STDMETHODIMP CMAS5::put_AlgorithmParameterFile(BSTR newVal)
{
    std::string fileName = COMStringUtils::ConvertString(newVal);
    MAS5ParameterExtraction::ExtractParameters(fileName.c_str(), mas5.AlgParameters());
    return S_OK;
}

STDMETHODIMP CMAS5::put_ReportControlsParameterFile(BSTR newVal)
{
    std::string fileName = COMStringUtils::ConvertString(newVal);
    ExpressionControlsParameterExtraction::ExtractParameters(fileName.c_str(), mas5.ProbePairThreshold(), mas5.ReportControls());
    return S_OK;
}

STDMETHODIMP CMAS5::put_ProgramName(BSTR newVal)
{
    mas5.ProgramName() = newVal;
    return S_OK;
}

STDMETHODIMP CMAS5::put_ProgramCompany(BSTR newVal)
{
    mas5.ProgramCompany() = newVal;
    return S_OK;
}

STDMETHODIMP CMAS5::put_ProgramId(BSTR newVal)
{
    mas5.ProgramId() = newVal;
    return S_OK;
}
