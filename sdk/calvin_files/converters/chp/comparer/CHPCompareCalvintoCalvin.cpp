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


#include "calvin_files/converters/chp/comparer/CHPCompareCalvintoCalvin.h"
//
#include "calvin_files/converters/chp/comparer/CHPCompareUtils.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cmath>
//

using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_comparer;
using namespace std;

/*
 * Clear the members.
 */
void CHPCompareCalvintoCalvin::Clear()
{
	differences = "";
	calvin1.Clear();
	calvin2.Clear();
}

/*
 * Open the files.
 */
bool CHPCompareCalvintoCalvin::OpenFiles(const char *file1, const char *file2)
{
	CHPFileReader reader1;
	CHPFileReader reader2;
	try
	{
		reader1.SetFilename(file1);
		reader1.Read(calvin1);
		reader2.SetFilename(file2);
		reader2.Read(calvin2);
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
bool CHPCompareCalvintoCalvin::CompareFiles(const char *file1, const char *file2)
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
void CHPCompareCalvintoCalvin::CompareFiles()
{
	CompareHeader();
	CompareBackground();
	CompareData();
}

/*
 * Compare the header.
 */
void CHPCompareCalvintoCalvin::CompareHeader()
{
	CompareNonParameterHeader();
	CompareAlgParams();
	CompareChipSummary();
}

/*
 * Compare the header (not including the alg parameters or chip summary).
 */
void CHPCompareCalvintoCalvin::CompareNonParameterHeader()
{
	CHPCompareUtils::CompareInts(calvin2.GetCols(), calvin1.GetCols(), differences, "columns");
	CHPCompareUtils::CompareInts(calvin2.GetRows(), calvin1.GetRows(), differences, "rows");
	if (calvin1.GetAssayType() != CHP_RESEQUENCING_ASSAY_TYPE)
		CHPCompareUtils::CompareInts(calvin2.GetEntryCount(), calvin1.GetEntryCount(), differences, "entry count");
	CHPCompareUtils::CompareStrings(calvin2.GetArrayType(), calvin1.GetArrayType(), differences, "array type");
	CHPCompareUtils::CompareStrings(calvin2.GetAlgName(), calvin1.GetAlgName(), differences, "alg name");
	CHPCompareUtils::CompareStrings(calvin2.GetAlgVersion(), calvin1.GetAlgVersion(), differences, "alg version");
	CHPCompareUtils::CompareStrings(calvin2.GetParentCell(), calvin1.GetParentCell(), differences, "parent cel");
	CHPCompareUtils::CompareStrings(calvin2.GetProgId(), calvin1.GetProgId(), differences, "prog ID");
}

/*
 * Compare the algorithm parameters.
 */
void CHPCompareCalvintoCalvin::CompareAlgParams()
{
	ParameterNameValueTypeVector cparams1 = calvin1.GetAlgParams();
	ParameterNameValueTypeVector cparams2 = calvin2.GetAlgParams();
	CHPCompareUtils::CompareInts((int)cparams1.size(), (int)cparams2.size(), differences, "algorithm parameters size");
	if (cparams1.size() != cparams2.size())
	{
		return;
	}
	int n = (int)cparams1.size();
	for (int i=0; i<n; i++)
	{
		ParameterNameValueType &cparam1 = cparams1[i];
		ParameterNameValueType &cparam2 = cparams2[i];
		CHPCompareUtils::CompareStrings(cparam1.GetName(), cparam2.GetName(), differences, "algorithm parameter name");
		CHPCompareUtils::CompareStrings(cparam1.ToString(), cparam2.ToString(), differences, "algorithm parameter value");
	}
}

/*
 * Compare the chip summary parameters.
 */
void CHPCompareCalvintoCalvin::CompareChipSummary()
{
	ParameterNameValueTypeVector cparams1 = calvin1.GetChipSums();
	ParameterNameValueTypeVector cparams2 = calvin2.GetChipSums();
	CHPCompareUtils::CompareInts((int)cparams1.size(), (int)cparams2.size(), differences, "chip summary size");
	if (cparams1.size() != cparams2.size())
	{
		return;
	}
	int n = (int)cparams1.size();
	for (int i=0; i<n; i++)
	{
		ParameterNameValueType &cparam1 = cparams1[i];
		ParameterNameValueType &cparam2 = cparams2[i];
		CHPCompareUtils::CompareStrings(cparam1.GetName(), cparam2.GetName(), differences, "chip summary name");
		CHPCompareUtils::CompareStrings(cparam1.ToString(), cparam2.ToString(), differences, "chip summary value");
	}
}

/*
 * Compare the background.
 */
void CHPCompareCalvintoCalvin::CompareBackground()
{
	CHPCompareUtils::CompareInts(calvin1.GetBackgroundZoneCnt(), calvin2.GetBackgroundZoneCnt(), differences, "bg zone count");
	if (calvin1.GetBackgroundZoneCnt() != calvin2.GetBackgroundZoneCnt())
	{
		return;
	}

	int n = calvin1.GetBackgroundZoneCnt();
	CHPBackgroundZone czone1;
	CHPBackgroundZone czone2;
	for (int i=0; i<n; i++)
	{
		calvin1.GetBackgroundZone(i, czone1);
		calvin2.GetBackgroundZone(i, czone2);
		CHPCompareUtils::CompareFloats(czone1.GetBackground(), czone2.GetBackground(), differences, "background");
		CHPCompareUtils::CompareFloats(czone1.GetCenterX(), czone2.GetCenterX(), differences, "zone centerx");
		CHPCompareUtils::CompareFloats(czone1.GetCenterY(), czone2.GetCenterY(), differences, "zone centery");
	}
}

/*
 * Compare the data
 */
void CHPCompareCalvintoCalvin::CompareData()
{
	if (calvin1.GetAssayType() == CHP_EXPRESSION_ASSAY_TYPE)
		CompareExpression();

	else if (calvin1.GetAssayType() == CHP_GENOTYPING_ASSAY_TYPE)
		CompareGenotyping();

	else if (calvin1.GetAssayType() == CHP_UNIVERSAL_ASSAY_TYPE)
		CompareUniversal();

	else if (calvin1.GetAssayType() == CHP_RESEQUENCING_ASSAY_TYPE)
		CompareResequencing();
}

/*
 * Compare the expression results
 */
void CHPCompareCalvintoCalvin::CompareExpression()
{
	int n=calvin1.GetEntryCount();
	for (int i=0; i<n; i++)
	{
		CHPExpressionEntry centry1;
		CHPExpressionEntry centry2;
		calvin1.GetEntry(i, centry1);
		calvin2.GetEntry(i, centry2);
		CHPCompareUtils::CompareFloats(centry1.GetDetectionPValue(), centry2.GetDetectionPValue(), differences, "detection p-value");
		CHPCompareUtils::CompareFloats(centry1.GetSignal(), centry2.GetSignal(), differences, "signal"); 
		CHPCompareUtils::CompareInts(centry1.GetNumPairs(), centry2.GetNumPairs(), differences, "pairs");
		CHPCompareUtils::CompareInts(centry1.GetNumPairsUsed(), centry2.GetNumPairsUsed(), differences, "pairs used");
		CHPCompareUtils::CompareInts(centry1.GetDetection(), centry2.GetDetection(), differences, "detection");
		CHPCompareUtils::CompareInts(centry1.GetHasComparisonData(), centry2.GetHasComparisonData(), differences, "comp flag");
		if (centry1.GetHasComparisonData())
		{
			CHPCompareUtils::CompareFloats(centry1.GetChangePValue(), centry2.GetChangePValue(), differences, "change p-value"); 
			CHPCompareUtils::CompareFloats(centry1.GetSigLogRatio(), centry2.GetSigLogRatio(), differences, "slr"); 
			CHPCompareUtils::CompareFloats(centry1.GetSigLogRatioLo(), centry2.GetSigLogRatioLo(), differences, "slr-low"); 
			CHPCompareUtils::CompareFloats(centry1.GetSigLogRatioHi(), centry2.GetSigLogRatioHi(), differences, "slr-high"); 
			CHPCompareUtils::CompareInts(centry1.GetCommonPairs(), centry2.GetCommonPairs(), differences, "common pairs");
			CHPCompareUtils::CompareInts(centry1.GetChange(), centry2.GetChange(), differences, "change");
		}
	}
}

/*
 * Compare the genotyping results.
 */
void CHPCompareCalvintoCalvin::CompareGenotyping()
{
	int n=calvin1.GetEntryCount();
	for (int i=0; i<n; i++)
	{
		CHPGenotypeEntry centry1;
		CHPGenotypeEntry centry2;
		calvin1.GetEntry(i, centry1);
		calvin2.GetEntry(i, centry2);
		CHPCompareUtils::CompareInts(centry1.GetCall(), centry2.GetCall(), differences, "call");
		CHPCompareUtils::CompareFloats(centry1.GetConfidence(), centry2.GetConfidence(), differences, "confidence"); 
		CHPCompareUtils::CompareFloats(centry1.GetAACall(), centry2.GetAACall(), differences, "aa p-value");
		CHPCompareUtils::CompareFloats(centry1.GetABCall(), centry2.GetABCall(), differences, "ab p-value");
		CHPCompareUtils::CompareFloats(centry1.GetBBCall(), centry2.GetBBCall(), differences, "bb p-value");
		CHPCompareUtils::CompareFloats(centry1.GetNoCall(), centry2.GetNoCall(), differences, "no call p-value");
		CHPCompareUtils::CompareFloats(centry1.GetRAS1(), centry2.GetRAS1(), differences, "ras1");
		CHPCompareUtils::CompareFloats(centry1.GetRAS2(), centry2.GetRAS2(), differences, "ras2");
	}
}

/*
 * Compare the universal results.
 */
void CHPCompareCalvintoCalvin::CompareUniversal()
{
	int n=calvin1.GetEntryCount();
	for (int i=0; i<n; i++)
	{
		CHPUniversalEntry centry1;
		CHPUniversalEntry centry2;
		calvin1.GetEntry(i, centry1);
		calvin2.GetEntry(i, centry2);
		CHPCompareUtils::CompareFloats(centry1.GetBackground(), centry2.GetBackground(), differences, "bg"); 
	}
}

/*
 * Compare resequencing results.
 */
void CHPCompareCalvintoCalvin::CompareResequencing()
{
	CHPCompareUtils::CompareInts(calvin1.GetEntryCount(), calvin2.GetEntryCount(), differences, "base count");
	if (calvin1.GetEntryCount() != calvin2.GetEntryCount())
		return;

	int n=calvin1.GetEntryCount();
	CHPReseqEntry ce1;
	CHPReseqEntry ce2;
	for (int i=0; i<n; i++)
	{
		calvin1.GetEntry(i, ce1);
		calvin1.GetEntry(i, ce2);
		CHPCompareUtils::CompareInts(ce1.call, ce2.call, differences, "base call");
		CHPCompareUtils::CompareFloats(ce1.score, ce2.score, differences, "base score");
	}

	CHPCompareUtils::CompareInts(calvin1.GetForceCnt(), calvin2.GetForceCnt(), differences, "force size");
	if (calvin1.GetForceCnt() != calvin2.GetForceCnt())
		return;
	n=calvin1.GetForceCnt();
	CHPReseqForceCall cf1;
	CHPReseqForceCall cf2;
	for (int i=0; i<n; i++)
	{
		calvin1.GetForceCall(i, cf1);
		calvin2.GetForceCall(i, cf2);
		CHPCompareUtils::CompareInts(cf1.position, cf2.position, differences, "force position");
		CHPCompareUtils::CompareInts(cf1.call, cf2.call, differences, "force call");
		CHPCompareUtils::CompareInts(cf1.reason, cf2.reason, differences, "force reason");
	}

	CHPCompareUtils::CompareInts(calvin1.GetOrigCnt(), calvin2.GetOrigCnt(), differences, "orig size");
	if (calvin1.GetOrigCnt() != calvin2.GetOrigCnt())
		return;
	n=calvin1.GetOrigCnt();
	CHPReseqOrigCall co1;
	CHPReseqOrigCall co2;
	for (int i=0; i<n; i++)
	{
		calvin1.GetOrigCall(i, co1);
		calvin2.GetOrigCall(i, co2);
		CHPCompareUtils::CompareInts(co1.position, co2.position, differences, "orig position");
		CHPCompareUtils::CompareInts(co1.call, co2.call, differences, "orig call");
	}
}
