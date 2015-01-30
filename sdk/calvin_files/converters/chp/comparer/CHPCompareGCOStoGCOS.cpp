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


#include "calvin_files/converters/chp/comparer/CHPCompareGCOStoGCOS.h"
//
#include "calvin_files/converters/chp/comparer/CHPCompareUtils.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cmath>
//

using namespace affxchp;
using namespace affymetrix_comparer;
using namespace std;

/*
 * Clear the members.
 */
void CHPCompareGCOStoGCOS::Clear()
{
	differences = "";
	gcos1.Clear();
	gcos2.Clear();
}

/*
 * Open the files.
 */
bool CHPCompareGCOStoGCOS::OpenFiles(const char *file1, const char *file2)
{
	gcos1.SetFileName(file1);
	if (gcos1.Read() == false)
	{
		differences = "Failed to read the GCOS file\r\n";
		return false;
	}
	gcos2.SetFileName(file2);
	if (gcos2.Read() == false)
	{
		differences = "Failed to read the GCOS file\r\n";
		return false;
	}
	return true;
}

/*
 * Compare the versions.
 */
bool CHPCompareGCOStoGCOS::CompareFiles(const char *file1, const char *file2)
{
	Clear();
	if (OpenFiles(file1, file2) == false)
		return false;
	CompareFiles();
	return (differences.length() == 0);
}

/*
 * Compare the versions.
 */
void CHPCompareGCOStoGCOS::CompareFiles()
{
	CompareHeader();
	CompareBackground();
	CompareData();
}

/*
 * Compare the header.
 */
void CHPCompareGCOStoGCOS::CompareHeader()
{
	CompareNonParameterHeader();
	CompareAlgParams();
	CompareChipSummary();
}

/*
 * Compare the header (not including the alg parameters or chip summary).
 */
void CHPCompareGCOStoGCOS::CompareNonParameterHeader()
{
	CCHPFileHeader &ghead1 = gcos1.GetHeader();
	CCHPFileHeader &ghead2 = gcos2.GetHeader();

	CHPCompareUtils::CompareInts(ghead1.GetCols(), ghead2.GetCols(), differences, "columns");
	CHPCompareUtils::CompareInts(ghead1.GetRows(), ghead2.GetRows(), differences, "rows");
	if (ghead1.GetAssayType() != CCHPFileHeader::Resequencing)
		CHPCompareUtils::CompareInts(ghead1.GetNumProbeSets(), ghead2.GetNumProbeSets(), differences, "entry count");
	CHPCompareUtils::CompareStrings(ghead1.GetChipType(), ghead2.GetChipType(), differences, "array type");
	CHPCompareUtils::CompareStrings(ghead1.GetAlgName(), ghead2.GetAlgName(), differences, "alg name");
	CHPCompareUtils::CompareStrings(ghead1.GetAlgVersion(), ghead2.GetAlgVersion(), differences, "alg version");
	CHPCompareUtils::CompareStrings(ghead1.GetParentCellFile(), ghead2.GetParentCellFile(), differences, "parent cel");
	CHPCompareUtils::CompareStrings(ghead1.GetProgID(), ghead2.GetProgID(), differences, "prog ID");
}

/*
 * Compare the algorithm parameters.
 */
void CHPCompareGCOStoGCOS::CompareAlgParams()
{
	TagValuePairTypeList &gparams1 = gcos1.GetHeader().AlgorithmParameters();
	TagValuePairTypeList &gparams2 = gcos2.GetHeader().AlgorithmParameters();
	CHPCompareUtils::CompareInts((int)gparams1.size(), (int)gparams2.size(), differences, "algorithm parameters size");
	if (gparams1.size() != gparams2.size())
	{
		return;
	}
	int index=0;
	TagValuePairTypeList::iterator it1 = gparams1.begin();
	TagValuePairTypeList::iterator it2 = gparams2.begin();
	while (it1!=gparams1.end())
	{
		TagValuePairType &gparam1 = *it1;
		TagValuePairType &gparam2 = *it2;
		CHPCompareUtils::CompareStrings(gparam1.Tag, gparam2.Tag, differences, "algorithm parameter name");
		CHPCompareUtils::CompareStrings(gparam1.Value, gparam2.Value, differences, "algorithm parameter value");

		++it1;
		++it2;
	}
}

/*
 * Compare the chip summary parameters.
 */
void CHPCompareGCOStoGCOS::CompareChipSummary()
{
	TagValuePairTypeList &gparams1 = gcos1.GetHeader().SummaryParameters();
	TagValuePairTypeList &gparams2 = gcos2.GetHeader().SummaryParameters();
	CHPCompareUtils::CompareInts((int)gparams1.size(), (int)gparams2.size(), differences, "chip summary size");
	if (gparams1.size() != gparams2.size())
	{
		return;
	}
	int index=0;
	TagValuePairTypeList::iterator it1 = gparams1.begin();
	TagValuePairTypeList::iterator it2 = gparams2.begin();
	while (it1!=gparams1.end())
	{
		TagValuePairType &gparam1 = *it1;
		TagValuePairType &gparam2 = *it2;
		CHPCompareUtils::CompareStrings(gparam1.Tag, gparam2.Tag, differences, "chip summary name");
		CHPCompareUtils::CompareStrings(gparam1.Value, gparam2.Value, differences, "chip summary value");

		++it1;
		++it2;
	}
}

/*
 * Compare the background.
 */
void CHPCompareGCOStoGCOS::CompareBackground()
{
	BackgroundZoneTypeList &zones1 = gcos1.GetHeader().GetBackgroundZones();
	BackgroundZoneTypeList &zones2 = gcos2.GetHeader().GetBackgroundZones();
	CHPCompareUtils::CompareInts((int)zones1.size(), (int)zones2.size(), differences, "bg zone count");
	if (zones1.size() != zones2.size())
	{
		return;
	}

	int index=0;
	BackgroundZoneTypeList::iterator it1=zones1.begin();
	BackgroundZoneTypeList::iterator it2=zones2.begin();
	while (it1 != zones1.end())
	{
		BackgroundZoneType &gzone1 = *it1;
		BackgroundZoneType &gzone2 = *it2;
		CHPCompareUtils::CompareFloats(gzone1.background, gzone2.background, differences, "background");
		CHPCompareUtils::CompareFloats(gzone1.centerx, gzone2.centerx, differences, "zone centerx");
		CHPCompareUtils::CompareFloats(gzone1.centery, gzone2.centery, differences, "zone centery");
		++it1;
		++it2;
	}
}

/*
 * Compare the data
 */
void CHPCompareGCOStoGCOS::CompareData()
{
	if (gcos1.GetHeader().GetAssayType() == CCHPFileHeader::Expression)
		CompareExpression();

	else if (gcos1.GetHeader().GetAssayType() == CCHPFileHeader::Genotyping)
		CompareGenotyping();

	else if (gcos1.GetHeader().GetAssayType() == CCHPFileHeader::Universal)
		CompareUniversal();

	else if (gcos1.GetHeader().GetAssayType() == CCHPFileHeader::Resequencing)
		CompareResequencing();
}

/*
 * Compare the expression results
 */
void CHPCompareGCOStoGCOS::CompareExpression()
{
	int n=gcos1.GetHeader().GetNumProbeSets();
	for (int i=0; i<n; i++)
	{
		CExpressionProbeSetResults *gentry1 = gcos1.GetExpressionResults(i);
		CExpressionProbeSetResults *gentry2 = gcos2.GetExpressionResults(i);

		CHPCompareUtils::CompareFloats(gentry1->DetectionPValue, gentry2->DetectionPValue, differences, "detection p-value");
		CHPCompareUtils::CompareFloats(gentry1->Signal, gentry2->Signal, differences, "signal"); 
		CHPCompareUtils::CompareInts(gentry1->NumPairs, gentry2->NumPairs, differences, "pairs");
		CHPCompareUtils::CompareInts(gentry1->NumUsedPairs, gentry2->NumUsedPairs, differences, "pairs used");
		CHPCompareUtils::CompareInts(gentry1->Detection, gentry2->Detection, differences, "detection");
		CHPCompareUtils::CompareInts(gentry1->m_HasCompResults, gentry2->m_HasCompResults, differences, "comp flag");
		if (gentry1->m_HasCompResults)
		{
			CHPCompareUtils::CompareFloats(gentry1->ChangePValue, gentry2->ChangePValue, differences, "change p-value"); 
			CHPCompareUtils::CompareFloats(gentry1->SignalLogRatio, gentry2->SignalLogRatio, differences, "slr"); 
			CHPCompareUtils::CompareFloats(gentry1->SignalLogRatioLow, gentry2->SignalLogRatioLow, differences, "slr-low"); 
			CHPCompareUtils::CompareFloats(gentry1->SignalLogRatioHigh, gentry2->SignalLogRatioHigh, differences, "slr-high"); 
			CHPCompareUtils::CompareInts(gentry1->NumCommonPairs, gentry2->NumCommonPairs, differences, "common pairs");
			CHPCompareUtils::CompareInts(gentry1->Change, gentry2->Change, differences, "change");
		}
	}
}

/*
 * Compare the genotyping results.
 */
void CHPCompareGCOStoGCOS::CompareGenotyping()
{
	int n=gcos1.GetHeader().GetNumProbeSets();
	for (int i=0; i<n; i++)
	{
		CGenotypeProbeSetResults *gentry1 = gcos1.GetGenotypingResults(i);
		CGenotypeProbeSetResults *gentry2 = gcos2.GetGenotypingResults(i);

		CHPCompareUtils::CompareInts(gentry1->AlleleCall, gentry2->AlleleCall, differences, "call");
		CHPCompareUtils::CompareFloats(gentry1->Confidence, gentry2->Confidence, differences, "confidence"); 
		CHPCompareUtils::CompareFloats(gentry1->pvalue_AA, gentry2->pvalue_AA, differences, "aa p-value");
		CHPCompareUtils::CompareFloats(gentry1->pvalue_AB, gentry2->pvalue_AB, differences, "ab p-value");
		CHPCompareUtils::CompareFloats(gentry1->pvalue_BB, gentry2->pvalue_BB, differences, "bb p-value");
		CHPCompareUtils::CompareFloats(gentry1->pvalue_NoCall, gentry2->pvalue_NoCall, differences, "no call p-value");
		CHPCompareUtils::CompareFloats(gentry1->RAS1, gentry2->RAS1, differences, "ras1");
		CHPCompareUtils::CompareFloats(gentry1->RAS2, gentry2->RAS2, differences, "ras2");
	}
}

/*
 * Compare the universal results.
 */
void CHPCompareGCOStoGCOS::CompareUniversal()
{
	int n=gcos1.GetHeader().GetNumProbeSets();
	for (int i=0; i<n; i++)
	{
		CUniversalProbeSetResults *gentry1 = gcos1.GetUniversalResults(i);
		CUniversalProbeSetResults *gentry2 = gcos2.GetUniversalResults(i);
		CHPCompareUtils::CompareFloats(gentry1->GetBackground(), gentry2->GetBackground(), differences, "bg"); 
	}
}

/*
 * Compare resequencing results.
 */
void CHPCompareGCOStoGCOS::CompareResequencing()
{
	CResequencingResults *gresults1 = gcos1.GetResequencingResults();
	CResequencingResults *gresults2 = gcos2.GetResequencingResults();
	CHPCompareUtils::CompareInts(gresults1->GetCalledBasesSize(), gresults2->GetCalledBasesSize(), differences, "base count");
	if (gresults1->GetCalledBasesSize() != gresults2->GetCalledBasesSize())
		return;

	int n=gresults1->GetCalledBasesSize();
	for (int i=0; i<n; i++)
	{
		CHPCompareUtils::CompareInts(gresults1->GetCalledBase(i), gresults2->GetCalledBase(i), differences, "base call");
		CHPCompareUtils::CompareFloats(gresults1->GetScore(i), gresults2->GetScore(i), differences, "base score");
	}

	CHPCompareUtils::CompareInts(gresults1->GetForceCallsSize(), gresults2->GetForceCallsSize(), differences, "force size");
	if (gresults1->GetForceCallsSize() != gresults2->GetForceCallsSize())
		return;
	n=gresults1->GetForceCallsSize();
	ForceCallType gf1;
	ForceCallType gf2;
	for (int i=0; i<n; i++)
	{
		gf1 = gresults1->GetForceCall(i);
		gf2 = gresults2->GetForceCall(i);
		CHPCompareUtils::CompareInts(gf1.position, gf2.position, differences, "force position");
		CHPCompareUtils::CompareInts(gf1.call, gf2.call, differences, "force call");
		CHPCompareUtils::CompareInts(gf1.reason, gf2.reason, differences, "force reason");
	}

	CHPCompareUtils::CompareInts(gresults1->GetOrigCallsSize(), gresults2->GetOrigCallsSize(), differences, "orig size");
	if (gresults1->GetOrigCallsSize() != gresults2->GetOrigCallsSize())
		return;
	n=gresults1->GetOrigCallsSize();
	BaseCallType go1;
	BaseCallType go2;
	for (int i=0; i<n; i++)
	{
		go1 = gresults1->GetOrigCall(i);
		go2 = gresults2->GetOrigCall(i);
		CHPCompareUtils::CompareInts(go1.position, go2.position, differences, "orig position");
		CHPCompareUtils::CompareInts(go1.call, go2.call, differences, "orig call");
	}
}
