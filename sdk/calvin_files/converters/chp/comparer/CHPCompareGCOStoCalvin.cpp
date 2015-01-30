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


#include "calvin_files/converters/chp/comparer/CHPCompareGCOStoCalvin.h"
//
#include "calvin_files/converters/chp/comparer/CHPCompareUtils.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cmath>
//

using namespace affxchp;
using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_comparer;
using namespace std;

/*
 * Clear the members.
 */
void CHPCompareGCOStoCalvin::Clear()
{
	differences = "";
	gcos.Clear();
	calvin.Clear();
}

/*
 * Open the files.
 */
bool CHPCompareGCOStoCalvin::OpenFiles(const char *gcosFile, const char *calvinFile)
{
	gcos.SetFileName(gcosFile);
	if (gcos.Read() == false)
	{
		differences = "Failed to read the GCOS file\r\n";
		return false;
	}

	CHPFileReader reader;
	try
	{
		reader.SetFilename(calvinFile);
		reader.Read(calvin);
	}
	catch (...)
	{
		differences = "Failed to read the Calvin file\r\n";
		return false;
	}
	return true;
}

/*
 * Compare the versions.
 */
bool CHPCompareGCOStoCalvin::CompareFiles(const char *gcosFile, const char *calvinFile)
{
	Clear();

	if (OpenFiles(gcosFile, calvinFile) == false)
		return false;

	CompareFiles();
	return (differences.length() == 0);
}

/*
 * Compare the versions.
 */
void CHPCompareGCOStoCalvin::CompareFiles()
{
	CompareHeader();
	CompareBackground();
	CompareData();
}

/*
 * Compare the header.
 */
void CHPCompareGCOStoCalvin::CompareHeader()
{
	CompareNonParameterHeader();
	CompareAlgParams();
	CompareChipSummary();
}

/*
 * Compare the header (not including the alg parameters or chip summary).
 */
void CHPCompareGCOStoCalvin::CompareNonParameterHeader()
{
	CCHPFileHeader &ghead = gcos.GetHeader();

	CHPCompareUtils::CompareInts(ghead.GetCols(), calvin.GetCols(), differences, "columns");
	CHPCompareUtils::CompareInts(ghead.GetRows(), calvin.GetRows(), differences, "rows");
	if (gcos.GetHeader().GetAssayType() != CCHPFileHeader::Resequencing)
		CHPCompareUtils::CompareInts(ghead.GetNumProbeSets(), calvin.GetEntryCount(), differences, "entry count");
	CHPCompareUtils::CompareStrings(ghead.GetChipType(), calvin.GetArrayType(), differences, "array type");
	CHPCompareUtils::CompareStrings(ghead.GetAlgName(), calvin.GetAlgName(), differences, "alg name");
	CHPCompareUtils::CompareStrings(ghead.GetAlgVersion(), calvin.GetAlgVersion(), differences, "alg version");
	CHPCompareUtils::CompareStrings(ghead.GetParentCellFile(), calvin.GetParentCell(), differences, "parent cel");
	CHPCompareUtils::CompareStrings(ghead.GetProgID(), calvin.GetProgId(), differences, "prog ID");
}

/*
 * Compare the algorithm parameters.
 */
void CHPCompareGCOStoCalvin::CompareAlgParams()
{
	TagValuePairTypeList &gparams = gcos.GetHeader().AlgorithmParameters();
	ParameterNameValueTypeVector cparams = calvin.GetAlgParams();
	CHPCompareUtils::CompareInts((int)gparams.size(), (int)cparams.size(), differences, "algorithm parameters size");
	if (gparams.size() != cparams.size())
	{
		return;
	}
	int index=0;
	for (TagValuePairTypeList::iterator it = gparams.begin(); it!=gparams.end(); ++it)
	{
		TagValuePairType &gparam = *it;
		ParameterNameValueType &cparam = cparams[index++];
		CHPCompareUtils::CompareStrings(gparam.Tag, cparam.GetName(), differences, "algorithm parameter name");
		CHPCompareUtils::CompareStrings(gparam.Value, cparam.ToString(), differences, "algorithm parameter value");
	}
}

/*
 * Compare the chip summary parameters.
 */
void CHPCompareGCOStoCalvin::CompareChipSummary()
{
	TagValuePairTypeList &gparams = gcos.GetHeader().SummaryParameters();
	ParameterNameValueTypeVector cparams = calvin.GetChipSums();
	CHPCompareUtils::CompareInts((int)gparams.size(), (int)cparams.size(), differences, "chip summary size");
	if (gparams.size() != cparams.size())
	{
		return;
	}
	int index=0;
	for (TagValuePairTypeList::iterator it = gparams.begin(); it!=gparams.end(); ++it)
	{
		TagValuePairType &gparam = *it;
		ParameterNameValueType &cparam = cparams[index++];
		CHPCompareUtils::CompareStrings(gparam.Tag, cparam.GetName(), differences, "chip summary name");
		CHPCompareUtils::CompareStrings(gparam.Value, cparam.ToString(), differences, "chip summary value");
	}
}

/*
 * Compare the background.
 */
void CHPCompareGCOStoCalvin::CompareBackground()
{
	BackgroundZoneTypeList &zones = gcos.GetHeader().GetBackgroundZones();
	CHPCompareUtils::CompareInts((int)zones.size(), calvin.GetBackgroundZoneCnt(), differences, "bg zone count");
	if ((int)zones.size() != calvin.GetBackgroundZoneCnt())
	{
		return;
	}

	int index=0;
	CHPBackgroundZone czone;
	for (BackgroundZoneTypeList::iterator it=zones.begin(); it!=zones.end(); ++it)
	{
		BackgroundZoneType &gzone = *it;
		calvin.GetBackgroundZone(index++, czone);
		CHPCompareUtils::CompareFloats(gzone.background, czone.GetBackground(), differences, "background");
		CHPCompareUtils::CompareFloats(czone.GetCenterX(), gzone.centerx, differences, "zone centerx");
		CHPCompareUtils::CompareFloats(czone.GetCenterY(), gzone.centery, differences, "zone centery");
	}
}

/*
 * Compare the data
 */
void CHPCompareGCOStoCalvin::CompareData()
{
	if (gcos.GetHeader().GetAssayType() == CCHPFileHeader::Expression)
		CompareExpression();

	else if (gcos.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping)
		CompareGenotyping();

	else if (gcos.GetHeader().GetAssayType() == CCHPFileHeader::Universal)
		CompareUniversal();

	else if (gcos.GetHeader().GetAssayType() == CCHPFileHeader::Resequencing)
		CompareResequencing();
}

/*
 * Compare the expression results
 */
void CHPCompareGCOStoCalvin::CompareExpression()
{
	int n=gcos.GetHeader().GetNumProbeSets();
	for (int i=0; i<n; i++)
	{
		CHPExpressionEntry centry;
		CExpressionProbeSetResults *gentry = gcos.GetExpressionResults(i);
		calvin.GetEntry(i, centry);
		CHPCompareUtils::CompareFloats(gentry->DetectionPValue, centry.GetDetectionPValue(), differences, "detection p-value");
		CHPCompareUtils::CompareFloats(gentry->Signal, centry.GetSignal(), differences, "signal"); 
		CHPCompareUtils::CompareInts(gentry->NumPairs, centry.GetNumPairs(), differences, "pairs");
		CHPCompareUtils::CompareInts(gentry->NumUsedPairs, centry.GetNumPairsUsed(), differences, "pairs used");
		CHPCompareUtils::CompareInts(gentry->Detection, centry.GetDetection(), differences, "detection");
		CHPCompareUtils::CompareInts(gentry->m_HasCompResults, centry.GetHasComparisonData(), differences, "comp flag");
		if (gentry->m_HasCompResults)
		{
			CHPCompareUtils::CompareFloats(gentry->ChangePValue, centry.GetChangePValue(), differences, "change p-value"); 
			CHPCompareUtils::CompareFloats(gentry->SignalLogRatio, centry.GetSigLogRatio(), differences, "slr"); 
			CHPCompareUtils::CompareFloats(gentry->SignalLogRatioLow, centry.GetSigLogRatioLo(), differences, "slr-low"); 
			CHPCompareUtils::CompareFloats(gentry->SignalLogRatioHigh, centry.GetSigLogRatioHi(), differences, "slr-high"); 
			CHPCompareUtils::CompareInts(gentry->NumCommonPairs, centry.GetCommonPairs(), differences, "common pairs");
			CHPCompareUtils::CompareInts(gentry->Change, centry.GetChange(), differences, "change");
		}
	}
}

/*
 * Compare the genotyping results.
 */
void CHPCompareGCOStoCalvin::CompareGenotyping()
{
	int n=gcos.GetHeader().GetNumProbeSets();
	for (int i=0; i<n; i++)
	{
		CHPGenotypeEntry centry;
		CGenotypeProbeSetResults *gentry = gcos.GetGenotypingResults(i);
		calvin.GetEntry(i, centry);
		CHPCompareUtils::CompareInts(gentry->AlleleCall, centry.GetCall(), differences, "call");
		CHPCompareUtils::CompareFloats(gentry->Confidence, centry.GetConfidence(), differences, "confidence"); 
		CHPCompareUtils::CompareFloats(gentry->pvalue_AA, centry.GetAACall(), differences, "aa p-value");
		CHPCompareUtils::CompareFloats(gentry->pvalue_AB, centry.GetABCall(), differences, "ab p-value");
		CHPCompareUtils::CompareFloats(gentry->pvalue_BB, centry.GetBBCall(), differences, "bb p-value");
		CHPCompareUtils::CompareFloats(gentry->pvalue_NoCall, centry.GetNoCall(), differences, "no call p-value");
		CHPCompareUtils::CompareFloats(gentry->RAS1, centry.GetRAS1(), differences, "ras1");
		CHPCompareUtils::CompareFloats(gentry->RAS2, centry.GetRAS2(), differences, "ras2");
	}
}

/*
 * Compare the universal results.
 */
void CHPCompareGCOStoCalvin::CompareUniversal()
{
	int n=gcos.GetHeader().GetNumProbeSets();
	for (int i=0; i<n; i++)
	{
		CHPUniversalEntry centry;
		CUniversalProbeSetResults *gentry = gcos.GetUniversalResults(i);
		calvin.GetEntry(i, centry);
		CHPCompareUtils::CompareFloats(gentry->GetBackground(), centry.GetBackground(), differences, "bg"); 
	}
}

/*
 * Compare resequencing results.
 */
void CHPCompareGCOStoCalvin::CompareResequencing()
{
	CResequencingResults *gresults = gcos.GetResequencingResults();
	CHPCompareUtils::CompareInts(gresults->GetCalledBasesSize(), calvin.GetEntryCount(), differences, "base count");
	if (gresults->GetCalledBasesSize() != calvin.GetEntryCount())
		return;

	int n=gresults->GetCalledBasesSize();
	CHPReseqEntry ce;
	for (int i=0; i<n; i++)
	{
		calvin.GetEntry(i, ce);
		CHPCompareUtils::CompareInts(gresults->GetCalledBase(i), ce.call, differences, "base call");
		CHPCompareUtils::CompareFloats(gresults->GetScore(i), ce.score, differences, "base score");
	}

	CHPCompareUtils::CompareInts(gresults->GetForceCallsSize(), calvin.GetForceCnt(), differences, "force size");
	if (gresults->GetForceCallsSize() != calvin.GetForceCnt())
		return;
	n=gresults->GetForceCallsSize();
	ForceCallType gf;
	CHPReseqForceCall cf;
	for (int i=0; i<n; i++)
	{
		gf = gresults->GetForceCall(i);
		calvin.GetForceCall(i, cf);
		CHPCompareUtils::CompareInts(gf.position, cf.position, differences, "force position");
		CHPCompareUtils::CompareInts(gf.call, cf.call, differences, "force call");
		CHPCompareUtils::CompareInts(gf.reason, cf.reason, differences, "force reason");
	}

	CHPCompareUtils::CompareInts(gresults->GetOrigCallsSize(), calvin.GetOrigCnt(), differences, "orig size");
	if (gresults->GetOrigCallsSize() != calvin.GetOrigCnt())
		return;
	n=gresults->GetOrigCallsSize();
	BaseCallType go;
	CHPReseqOrigCall co;
	for (int i=0; i<n; i++)
	{
		go = gresults->GetOrigCall(i);
		calvin.GetOrigCall(i, co);
		CHPCompareUtils::CompareInts(go.position, co.position, differences, "orig position");
		CHPCompareUtils::CompareInts(go.call, co.call, differences, "orig call");
	}
}
