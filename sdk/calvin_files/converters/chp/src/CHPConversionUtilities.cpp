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


#include "calvin_files/converters/chp/src/CHPConversionUtilities.h"
//
#include "calvin_files/data/src/GenericDataTypes.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstring>
#include <string>
//

using namespace std;
using namespace affxchp;
using namespace affxcdf;
using namespace affxchpwriter;
using namespace affymetrix_calvin_io;
using namespace affymetrix_chp_converter;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_parameter;

#ifdef WIN32
#pragma warning(disable:4996) // Don't show deprecated messages.
#endif

/*
 * Copy the file contents.
 */
void CHPConversionUtilities::Copy(CHPData &inFile, CCHPFileWriter &outFile)
{
	GeneChipProbeSetType probeSetType=UnknownProbeSetType;
	string assayType = inFile.GetAssayType();
	if(assayType == CHP_EXPRESSION_ASSAY_TYPE)
		probeSetType = ExpressionProbeSetType;

	else if(assayType == CHP_GENOTYPING_ASSAY_TYPE)
		probeSetType = GenotypingProbeSetType;

	else if(assayType == CHP_RESEQUENCING_ASSAY_TYPE)
		probeSetType = ResequencingProbeSetType;

	else if(assayType == CHP_UNIVERSAL_ASSAY_TYPE)
		probeSetType = TagProbeSetType;

	outFile.InitializeForWriting(
		inFile.GetRows(), inFile.GetCols(), inFile.GetEntryCount(),
		StringUtils::ConvertWCSToMBS(inFile.GetArrayType()).c_str(), probeSetType);
	SetAlgName(outFile, StringUtils::ConvertWCSToMBS(inFile.GetAlgName()));
	outFile.SetAlgorithmVersion(StringUtils::ConvertWCSToMBS(inFile.GetAlgVersion()).c_str());
	outFile.SetParentCelFileName(StringUtils::ConvertWCSToMBS(inFile.GetParentCell()).c_str());
	outFile.SetProgID(StringUtils::ConvertWCSToMBS(inFile.GetProgId()).c_str());

	ParameterNameValueTypeVector params = inFile.GetAlgParams();
	int n = (int) params.size();
	for (int i=0; i<n; i++)
	{
		ParameterNameValueType &param = params[i];
		outFile.AddAlgorithmParameter(
			StringUtils::ConvertWCSToMBS(param.GetName()).c_str(),
			StringUtils::ConvertWCSToMBS(param.ToString()).c_str());
	}

	params = inFile.GetChipSums();
	n = (int) params.size();
	for (int i=0; i<n; i++)
	{
		ParameterNameValueType &param = params[i];
		outFile.AddChipSummaryParameter(
			StringUtils::ConvertWCSToMBS(param.GetName()).c_str(),
			StringUtils::ConvertWCSToMBS(param.ToString()).c_str());
	}

	n = inFile.GetBackgroundZoneCnt();
	CHPBackgroundZone zone;
	for (int i=0; i<n; i++)
	{
		inFile.GetBackgroundZone(i, zone);
		if (i == 0)
		{
			outFile.AddBackgroundInfo(n, zone.GetSmoothFactor());
		}
		outFile.AddBackgroundZone((int)zone.GetCenterX(), (int)zone.GetCenterY(), zone.GetBackground());
	}

	if (probeSetType == ExpressionProbeSetType)
	{
		CHPExpressionEntry entry;
		CExpressionProbeSetResults results;
		int n = inFile.GetEntryCount();
		for (int i=0; i<n; i++)
		{
			inFile.GetEntry(i, entry);
			results.Signal = entry.GetSignal();
			results.Detection = entry.GetDetection();
			results.DetectionPValue = entry.GetDetectionPValue();
			results.NumPairs = entry.GetNumPairs();
			results.NumUsedPairs = entry.GetNumPairsUsed();

			results.m_HasCompResults = entry.GetHasComparisonData();
			results.Change = entry.GetChange();
			results.ChangePValue = entry.GetChangePValue();
			results.SignalLogRatio = entry.GetSigLogRatio();
			results.SignalLogRatioLow = entry.GetSigLogRatioLo();
			results.SignalLogRatioHigh = entry.GetSigLogRatioHi();
			results.NumCommonPairs = entry.GetCommonPairs();

			outFile.SetExpressionEntry(i, &results);
		}
	}
	else if (probeSetType == GenotypingProbeSetType)
	{
		CHPGenotypeEntry entry;
		CGenotypeProbeSetResults results;
		int n = inFile.GetEntryCount();
		for (int i=0; i<n; i++)
		{
			inFile.GetEntry(i, entry);

			results.AlleleCall = entry.GetCall();
			results.Confidence = entry.GetConfidence();
			results.pvalue_AA = entry.GetAACall();
			results.pvalue_AB = entry.GetABCall();
			results.pvalue_BB = entry.GetBBCall();
			results.pvalue_NoCall = entry.GetNoCall();
			results.RAS1 = entry.GetRAS1();
			results.RAS2 = entry.GetRAS2();

			outFile.SetMappingEntry(i, &results);
		}
	}
	else if (probeSetType == ResequencingProbeSetType)
	{
		ForceCallType outF;
		CHPReseqForceCall inF;
		BaseCallType outO;
		CHPReseqOrigCall inO;
		CHPReseqEntry entry;
		CResequencingResults results;
		int n = inFile.GetEntryCount();
		results.ResizeCalledBases(n);
		results.ResizeScores(n);
		for (int i=0; i<n; i++)
		{
			inFile.GetEntry(i, entry);
			results.SetCalledBase(i, entry.call);
			results.SetScore(i, entry.score);
		}

		n = inFile.GetForceCnt();
		results.ResizeForceCalls(n);
		for (int i=0; i<n; i++)
		{
			inFile.GetForceCall(i, inF);
			outF.position = inF.position;
			outF.call = inF.call;
			outF.reason = inF.reason;
			results.SetForceCall(i, outF);
		}

		n = inFile.GetOrigCnt();
		results.ResizeOrigCalls(n);
		for (int i=0; i<n; i++)
		{
			inFile.GetOrigCall(i, inO);
			outO.position = inO.position;
			outO.call = inO.call;
			results.SetOrigCall(i, outO);
		}

		outFile.SetResequencingResults(&results);

	}
	else if (probeSetType == TagProbeSetType)
	{
		CHPUniversalEntry entry;
		CUniversalProbeSetResults results;
		int n = inFile.GetEntryCount();
		for (int i=0; i<n; i++)
		{
			inFile.GetEntry(i, entry);
			results.SetBackground(entry.GetBackground());
			outFile.SetUniversalEntry(i, &results);
		}
	}
}

void CHPConversionUtilities::SetAlgName(affxchpwriter::CCHPFileWriter& outFile, std::string algName)
{
	if (outFile.GetHeader().GetAssayType() == CCHPFileHeader::Expression ||
			outFile.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping)
	{
		// Change the AssayType before setting the algorithm name to side-step
		// CCHPFileWriter::SetAlgorithmName desire to prefix the algorithm name with
		// either "Expression" or "Genotyping" depending on the AssayType.
		CCHPFileHeader::GeneChipAssayType originalAssayType = outFile.GetHeader().GetAssayType();
		outFile.GetHeader().SetAssayType(CCHPFileHeader::Universal);
		outFile.SetAlgorithmName(algName.c_str());
		outFile.GetHeader().SetAssayType(originalAssayType);
	}

	else
	{
		outFile.SetAlgorithmName(algName.c_str());
	}
}

