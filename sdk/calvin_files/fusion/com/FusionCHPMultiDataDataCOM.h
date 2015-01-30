////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
#include "calvin_files/fusion/src/FusionCHPMultiDataData.h"

// CFusionCHPMultiDataDataCOM

class ATL_NO_VTABLE CFusionCHPMultiDataDataCOM :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFusionCHPMultiDataDataCOM, &CLSID_FusionCHPMultiDataData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFusionCHPMultiDataData, &IID_IFusionCHPMultiDataData, &LIBID_affx_fusion_comLib, /*wMajor =*/ 1, /*wMinor =*/ 0>
{
public:
	CFusionCHPMultiDataDataCOM()
	{
		chp = NULL;
        forceIndex = -1;
        sigAIndex = -1;
        sigBIndex = -1;
		cnStateIndex = -1;
		lohStateIndex = -1;
		log2RatioIndex = -1;
	}

DECLARE_REGISTRY_RESOURCEID(IDR_FUSIONCHPMULTIDATADATA)


BEGIN_COM_MAP(CFusionCHPMultiDataDataCOM)
	COM_INTERFACE_ENTRY(IFusionCHPMultiDataData)
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
		delete chp;
        chp = NULL;
	}

private:
	affymetrix_fusion_io::FusionCHPMultiDataData *chp;
    int forceIndex;
    int sigAIndex;
    int sigBIndex;
	int cnStateIndex;
	int lohStateIndex;
	int log2RatioIndex;
    int signalIndex;
    int callIndex;
    int confidenceIndex;

public:
	STDMETHOD(get_FileTypeIdentifiers)(VARIANT* pVal);
public:
	STDMETHOD(get_FileId)(BSTR* pVal);
public:
	STDMETHOD(FromBase)(IFusionCHPData * baseChp, VARIANT_BOOL* pVal);
public:
	STDMETHOD(get_AlgName)(BSTR* pVal);
public:
	STDMETHOD(get_AlgVersion)(BSTR* pVal);
public:
	STDMETHOD(get_ArrayType)(BSTR* pVal);
public:
	STDMETHOD(get_AlgorithmParameters)(VARIANT* pVal);
public:
	STDMETHOD(get_SummaryParameters)(VARIANT* pVal);
public:
	STDMETHOD(GetEntryCount)(MultiDataType dataType, int* pVal);
public:
	STDMETHOD(GetGenotypeEntry)(MultiDataType dataType, int index, IProbeSetMultiDataGenotypeData * pVal);
public:
    STDMETHOD(GetGenoCall)(MultiDataType dataType, int index, int * pVal);
public:
	STDMETHOD(GetGenoConfidence)(MultiDataType dataType, int index, float * pVal);
public:
    STDMETHOD(GetProbeSetName)(MultiDataType dataType, int index, BSTR * pVal);
public:
    STDMETHOD(GetExpressionEntry)(MultiDataType dataType, int index, IProbeSetMultiDataExpressionData * pVal);
public:
    STDMETHOD(GetExpressionQuantification)(MultiDataType dataType, int index, float * pVal);
public:
    STDMETHOD(GetNumMetricColumns)(MultiDataType dataType, int * pVal);
public:
    STDMETHOD(GetMetricColumnName)(MultiDataType dataType, int index, BSTR * pVal);
public:
    STDMETHOD(GetProbesetNamesAndGenotypeCalls)(MultiDataType dataType, VARIANT* names, VARIANT* calls, int startIndex, int *pVal);
public:
    STDMETHOD(GetProbesetNamesAndExpressionSignals)(MultiDataType dataType, VARIANT* names, VARIANT* signals, int startIndex, int *pVal);
public:
    STDMETHOD(InitializeGetGenotypeEntries)(MultiDataType dataType, BSTR* forcedName, BSTR *signalAName, BSTR *signalBName);
public:
    STDMETHOD(GetGenotypeEntries)(MultiDataType dataType, VARIANT* names, VARIANT* calls, VARIANT* confidences, VARIANT* forcedcalls, VARIANT* signalsA, VARIANT* signalsB, int startIndex, int *pVal);
public:
    STDMETHOD(GetGenotypeEntriesGivenIndicies)(MultiDataType dataType, VARIANT* indicies, VARIANT* names, VARIANT* calls, VARIANT* confidences, VARIANT* forcedcalls, VARIANT* signalsA, VARIANT* signalsB);
public:
    STDMETHOD(Close)();
public:
    STDMETHOD(GetCopyNumberEntry)(MultiDataType dataType, int index, IProbeSetMultiDataCopyNumberData* pVal);
public:
    STDMETHOD(GetCytoRegionEntry)(MultiDataType dataType, int index, IProbeSetMultiDataCytoRegionData* pVal);
public:
    STDMETHOD(GetCytoRegionEntries)(MultiDataType dataType, VARIANT* names, VARIANT* calls, VARIANT* confidences, VARIANT* callsLOH, VARIANT* confidencesLOH, int startIndex, int *pVal);
public:
	STDMETHOD(GetDataSetHeaderParameters)(MultiDataType dataType, VARIANT* pVal);
public:
	STDMETHOD(InitializeGetCNStateEntries)(MultiDataType dataType);
public:
	STDMETHOD(GetCopyNumberState)(MultiDataType dataType, VARIANT* names, VARIANT* chromosomes, VARIANT* positions, VARIANT* cnStates, int startIndex, int *pVal);
public:
    STDMETHOD(InitializeGetLOHStateEntries)(MultiDataType dataType);
public:
	STDMETHOD(GetLOHState)(MultiDataType dataType, VARIANT* names, VARIANT* chromosomes, VARIANT* positions, VARIANT* lohStates, int startIndex, int *pVal);
public:
    STDMETHOD(InitializeGetLog2RatioEntries)(MultiDataType dataType);
public:
	STDMETHOD(GetLog2Ratio)(MultiDataType dataType, int startIndex, int count, VARIANT* names, VARIANT* positions, VARIANT* values);
public:
	STDMETHOD(GetLog2Ratio2)(MultiDataType dataType, int startIndex, int count, VARIANT* values);
public:
    STDMETHOD(GetCopyNumberVariationRegionEntry)(MultiDataType dataType, int index, IProbeSetMultiDataCopyNumberVariationRegionData* pVal);
public:
    STDMETHOD(GetCopyNumberVariationRegionNamesAndResults)(MultiDataType dataType, VARIANT* regionNames, VARIANT* calls,  VARIANT* signals, VARIANT* confidenceScores, int *pVal);
public:
	STDMETHOD(InitializeGetCopyNumberVariationEntries)(MultiDataType dataType);
public:
    STDMETHOD(GetCopyNumberVariationResults)(MultiDataType dataType, VARIANT* calls,  VARIANT* signals, VARIANT* confidenceScores, int *pVal);
public:
    STDMETHOD(GetCopyNumberEntryLog2Ratio)(MultiDataType dataType, int index, float* pVal);

};

OBJECT_ENTRY_AUTO(__uuidof(FusionCHPMultiDataData), CFusionCHPMultiDataDataCOM)
