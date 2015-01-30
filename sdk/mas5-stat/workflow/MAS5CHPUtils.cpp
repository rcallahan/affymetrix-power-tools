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

#include "mas5-stat/workflow/MAS5CHPUtils.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPFileWriter.h"
//
#include <sstream>
#include <stdio.h>
//

using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_utilities;
using namespace ExpressionReport;
using namespace std;


/*! The algorith name. */
#define MAS5_ALG_NAME "ExpressionStat"

/*! The algorithm version. */
#define MAS5_ALG_VERSION "5.0"

/*! The prog ID of the GCOS component that implemented the MAS5 algorithm. */
#define MAS5_PROG_ID "GeneChip.CallGEBaseCall.1"

/*! The program name. */
#define PROGRAM_NAME_VALUE L"MAS5"

/*! The company name. */
#define PROGRAM_COMPANY_VALUE L"Affymetrix"

/*! The ID of the program. */
#define PROGRAM_ID_VALUE L"5.0"

/*
 * Convert to a string.
 */
string ToString(int value)
{
	ostringstream str;
	str << value;
	return str.str();
}

/*
 * Convert to a string.
 */
string ToString(float value)
{
    char buf[64];
    sprintf(buf, "%0.6f", value);
	return buf;
}

/*
 * Store the program information.
 */
void MAS5CHPUtils::SetProgramInformation(const std::wstring &name, const std::wstring &co, const std::wstring &id)
{
    programName = name;
    programCompany = co;
    programId = id;
}

/*
 * Add the parameter to the alg parameters section.
 */
void MAS5CHPUtils::AddAlgParam(const std::wstring &name, const std::wstring &value)
{
	if (cc_data != NULL)
		cc_data->AddAlgParam(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddAlgorithmParameter(StringUtils::ConvertWCSToMBS(name).c_str(), StringUtils::ConvertWCSToMBS(value).c_str());
}

/*
 * Add the parameter to the alg parameters section.
 */
void MAS5CHPUtils::AddAlgParam(const std::wstring &name, int value)
{
	if (cc_data != NULL)
		cc_data->AddAlgParam(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddAlgorithmParameter(StringUtils::ConvertWCSToMBS(name).c_str(), ToString(value).c_str());
}

/*
 * Add the parameter to the alg parameters section.
 */
void MAS5CHPUtils::AddAlgParam(const std::wstring &name, float value)
{
	if (cc_data != NULL)
		cc_data->AddAlgParam(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddAlgorithmParameter(StringUtils::ConvertWCSToMBS(name).c_str(), ToString(value).c_str());
}

/*
 * Clear the members.
 */
void MAS5CHPUtils::Clear()
{
	cc_data = NULL;
	gcos_data = NULL;
	mas5_data = NULL;
	report_data = NULL;
    programName.clear();
    programCompany.clear();
    programId.clear();
}

/*
 * Store the algorithm parameters to the chp object.
 */
void MAS5CHPUtils::StoreAlgParams()
{
	CExpStatAlgSettings &params = mas5_data->GetParameters();

	AddAlgParam(L"HZ", params.NumberHorZones);
	AddAlgParam(L"VZ", params.NumberVertZones);
	AddAlgParam(L"BG", (int) params.IntensityLowPercent);
	AddAlgParam(L"Alpha1", params.Alpha1);
	AddAlgParam(L"Alpha2", params.Alpha2);
	AddAlgParam(L"Tau", params.Tau);

	if (params.SFMethod != CExpStatAlgSettings::DEFINED_SCALING_FACTOR)
	{
		AddAlgParam(L"TGT", (int) params.TGT);
	}

	AddAlgParam(L"SF", params.ScaleFactor);

	AddAlgParam(L"NF", params.NormFactor);

	if (params.ProbeMaskFile.length() > 0)
		AddAlgParam(L"ProbeMask", StringUtils::ConvertMBSToWCS(params.ProbeMaskFile));

	if (params.SFMethod == CExpStatAlgSettings::SCALE_TO_ALL_PROBE_SETS)
		AddAlgParam(L"ScaleMask", L"All");
	else if (params.SFMethod == CExpStatAlgSettings::SCALE_TO_SELECTED_PROBE_SETS)
		AddAlgParam(L"ScaleMask", StringUtils::ConvertMBSToWCS(params.ScaleMaskFile));

	if (mas5_data->GetCompStatResult(0) != NULL)
	{
		AddAlgParam(L"Gamma1L", params.Gamma1L);
		AddAlgParam(L"Gamma1H", params.Gamma1H);
		AddAlgParam(L"Gamma2L", params.Gamma2L);
		AddAlgParam(L"Gamma2H", params.Gamma2H);
		AddAlgParam(L"Perturbation", params.Perturbation);

		if (params.NormMethod == CExpStatAlgSettings::NORM_TO_ALL_PROBE_SETS)
			AddAlgParam(L"NormMask", L"All");
		else if (params.NormMethod == CExpStatAlgSettings::NORM_TO_SELECTED_PROBE_SETS)
			AddAlgParam(L"NormMask", StringUtils::ConvertMBSToWCS(params.NormMaskFile));

		AddAlgParam(L"BaselineSF", params.BaseScaleFactor);
		AddAlgParam(L"BF", StringUtils::ConvertMBSToWCS(mas5_data->GetBaselineCellData().GetFileName()));
	}
}

/*
 * Store the summary parameters to the chp object.
 */
void MAS5CHPUtils::AddChipSum(const std::wstring &name, const std::wstring &value)
{
	if (cc_data != NULL)
		cc_data->AddChipSum(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddChipSummaryParameter(StringUtils::ConvertWCSToMBS(name).c_str(), StringUtils::ConvertWCSToMBS(value).c_str());
}

/*
 * Store the summary parameters to the chp object.
 */
void MAS5CHPUtils::AddChipSum(const std::wstring &name, int value)
{
	if (cc_data != NULL)
		cc_data->AddChipSum(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddChipSummaryParameter(StringUtils::ConvertWCSToMBS(name).c_str(), ToString(value).c_str());
}

/*
 * Store the summary parameters to the chp object.
 */
void MAS5CHPUtils::AddChipSum(const std::wstring &name, float value)
{
	if (cc_data != NULL)
		cc_data->AddChipSum(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddChipSummaryParameter(StringUtils::ConvertWCSToMBS(name).c_str(), ToString(value).c_str());
}

/*
 * Store the summary parameters to the chp object.
 */
void MAS5CHPUtils::StoreSummaryParams(bool separateEntries)
{
	wchar_t buf[256];
	AddChipSum(L"RawQ", mas5_data->GetRawQ());

	AvgStdvMinMaxType &bg = mas5_data->GetBgStats();
    if (separateEntries == false)
    {
	    swprintf(buf, 256, L"Avg:%0.2f,Std:%0.2f,Min:%0.1f,Max:%0.1f", bg.avg, bg.stdv, bg.min, bg.max);
	    AddChipSum(L"Background", buf);
    }
    else
    {
	    AddChipSum(L"BG Avg", bg.avg);
	    AddChipSum(L"BG Std", bg.stdv);
	    AddChipSum(L"BG Min", bg.min);
	    AddChipSum(L"BG Max", bg.max);
    }

	AvgStdvMinMaxType &ns = mas5_data->GetNoiseStats();
    if (separateEntries == false)
    {
	    swprintf(buf, 256, L"Avg:%0.2f,Std:%0.2f,Min:%0.1f,Max:%0.1f", ns.avg, ns.stdv, ns.min, ns.max);
	    AddChipSum(L"Noise", buf);
    }
    else
    {
	    AddChipSum(L"Noise Avg", ns.avg);
	    AddChipSum(L"Noise Std", ns.stdv);
	    AddChipSum(L"Noise Min", ns.min);
	    AddChipSum(L"Noise Max", ns.max);
    }

	ControlInformationList &cntrls = mas5_data->GetControlInfo();
	for (ControlInformationList::iterator it=cntrls.begin(); it!=cntrls.end(); ++it)
	{
		ControlInformationType &cntrl = *it;
        if (separateEntries == false)
        {
		    swprintf(buf, 256, L"Avg:%0.0f,Count:%d", cntrl.avg, cntrl.count);
		    if (cntrl.qcType == affxcdf::CheckerboardPositiveQCProbeSetType)
			    AddChipSum(L"Corner+", buf);
		    else if (cntrl.qcType == affxcdf::CheckerboardNegativeQCProbeSetType)
			    AddChipSum(L"Corner-", buf);
		    else if (cntrl.qcType == affxcdf::CentralCrossPositiveQCProbeSetType)
			    AddChipSum(L"Central+", buf);
		    else if (cntrl.qcType == affxcdf::CentralCrossNegativeQCProbeSetType)
			    AddChipSum(L"Central-", buf);
        }
        else
        {
		    if (cntrl.qcType == affxcdf::CheckerboardPositiveQCProbeSetType)
            {
			    AddChipSum(L"Corner+ Avg", (int) (cntrl.avg+0.5f));
			    AddChipSum(L"Corner+ Count", cntrl.count);
            }
            else if (cntrl.qcType == affxcdf::CheckerboardNegativeQCProbeSetType)
            {
			    AddChipSum(L"Corner- Avg", (int) (cntrl.avg+0.5f));
			    AddChipSum(L"Corner- Count", cntrl.count);
            }
		    else if (cntrl.qcType == affxcdf::CentralCrossPositiveQCProbeSetType)
            {
			    AddChipSum(L"Central+ Avg", (int) (cntrl.avg+0.5f));
			    AddChipSum(L"Central+ Count", cntrl.count);
            }
		    else if (cntrl.qcType == affxcdf::CentralCrossNegativeQCProbeSetType)
            {
			    AddChipSum(L"Central- Avg", (int) (cntrl.avg+0.5f));
			    AddChipSum(L"Central- Count", cntrl.count);
            }
        }
	}

	if (mas5_data->GetCompStatResult(0) != NULL)
	{
		AddAlgParam(L"BaselineRawQ", mas5_data->GetBaselineRawQ());
	}
}

/*
 * Store the detection statistics.
 */
void MAS5CHPUtils::StoreReportDetectionStats()
{
    double dval;
	int total = report_data->Results().ProbeSetResults().NumSets();

	int pcount = report_data->Results().ProbeSetResults().PresentCalls().Count();
	double psignal = report_data->Results().ProbeSetResults().PresentCalls().Signal();
	AddChipSum(L"#P", pcount);
    dval = (total > 0 ? 100.0*pcount/total : 0);
	AddChipSum(L"%P", (float)dval);
    dval = (pcount > 0 ? psignal/pcount : 0);
    AddChipSum(L"Signal(P)", (float)dval);

	int mcount = report_data->Results().ProbeSetResults().MarginalCalls().Count();
	double msignal = report_data->Results().ProbeSetResults().MarginalCalls().Signal();
	AddChipSum(L"#M", mcount);
    dval = (total > 0 ? 100.0*mcount/total : 0);
	AddChipSum(L"%M", (float)dval);
    dval = (mcount > 0 ? msignal/mcount : 0);
	AddChipSum(L"Signal(M)", (float)dval);

	int acount = report_data->Results().ProbeSetResults().AbsentCalls().Count();
	double asignal = report_data->Results().ProbeSetResults().AbsentCalls().Signal();
	AddChipSum(L"#A", acount);
    dval = (total > 0 ? 100.0*acount/total : 0);
	AddChipSum(L"%A", (float)dval);
    dval = (acount > 0 ? asignal/acount : 0);
	AddChipSum(L"Signal(A)", (float)dval);

    dval = (total > 0 ? (psignal+msignal+asignal)/total : 0);
	AddChipSum(L"Signal(All)", (float)dval);
}

/*
 * Store the report parameters.
 */
void MAS5CHPUtils::StoreReportParameters()
{
	int total = report_data->Results().ProbeSetResults().NumSets();
	AddChipSum(L"#Probe Sets Exceeding Probe Pair Threshold", total);
	AddChipSum(L"Probe Pair Threshold", report_data->Results().ProbePairThreshold());
	if (report_data->Results().AntiSenseControls() == true)
		AddChipSum(L"Control Direction", L"Antisense");
	else
		AddChipSum(L"Control Direction", L"Sense");
}

/*
 * Store a control statistic.
 */
void MAS5CHPUtils::StoreReportControlStat(const wchar_t *controlType,
                                          const ExpressionControlResult &result)
{
	if (result.HasControlResult(ExpressionControl::THREE_PRIME_PROBE_SET) == false &&
		result.HasControlResult(ExpressionControl::THREE_PRIME_PROBE_SET) == false &&
		result.HasControlResult(ExpressionControl::THREE_PRIME_PROBE_SET) == false)
		return;

	// Store the 3/5/M values.
	const wchar_t *Detection[4] = {L"P", L"M", L"A", L"No Call"};
	const wchar_t *ControlValueString[3] = {L"3", L"M", L"5"};
	wchar_t name[256];
	wchar_t value[256];
	for (int i=0; i<NUM_SETS_PER_CONTROL; i++)
	{
		if (result.HasControlResult((ExpressionControl::ControlValueType)i) == true)
		{
			swprintf(name, sizeof(name), 
               L"%s_%s_%s_signal", controlType,
               StringUtils::ConvertMBSToWCS(result.GetName()).c_str(), ControlValueString[i]);
			AddChipSum(name, result.GetControlSignalResult((ExpressionControl::ControlValueType)i));

			swprintf(name, sizeof(name),
               L"%s_%s_%s_detection", controlType,
               StringUtils::ConvertMBSToWCS(result.GetName()).c_str(), ControlValueString[i]);
			swprintf(value, sizeof(value),
               L"%s", Detection[result.GetControlDetectionResult((ExpressionControl::ControlValueType)i)]);
			AddChipSum(name, value);
		}
	}

	// Store the average value.
	float avg = 0.0f;
	int iTotal = 0;
	for (int i=0; i<NUM_SETS_PER_CONTROL; i++)
	{
		if (result.HasControlResult((ExpressionControl::ControlValueType)i) == true)
		{
			avg += result.GetControlSignalResult((ExpressionControl::ControlValueType)i);
			++iTotal;
		}
	}
	if (iTotal == 0)
		avg = 0.0f;
	else
		avg /= iTotal;

	swprintf(name, 256, L"%s_%s_avg-signal", controlType, StringUtils::ConvertMBSToWCS(result.GetName()).c_str());
	AddChipSum(name, avg);

	// Store the 3/5 ratio.
	float ratio = result.GetThreeFiveRatio();
	if (ratio >= 0)
	{
		swprintf(name, 256, L"%s_%s_3-5-ratio", controlType, StringUtils::ConvertMBSToWCS(result.GetName()).c_str());
		AddChipSum(name, ratio);
	}
}

/*
 * Store the control statistics.
 */
void MAS5CHPUtils::StoreReportControlStats()
{
	ExpressionControlResultList::const_iterator it;
	for (it = report_data->Results().SpikeStats().begin(); it != report_data->Results().SpikeStats().end(); ++it)
	{
		const ExpressionControlResult &result = *it;
		StoreReportControlStat(L"Spike", result);
	}
	for (it = report_data->Results().HousekeepingStats().begin(); it != report_data->Results().HousekeepingStats().end(); ++it)
	{
		const ExpressionControlResult &result = *it;
		StoreReportControlStat(L"Housekeeping", result);
	}
}

/*
 * Store the change statistics.
 */
void MAS5CHPUtils::StoreReportChangeStats()
{
	if (report_data->DataAccessor()->HasComparisonData() == false)
		return;

	int total = report_data->Results().ProbeSetResults().NumSets();

	AddChipSum(L"#I", report_data->Results().IncreaseStats().ChangeCount());
	if (total > 0)
	{
		AddChipSum(L"%I", 100.0f*report_data->Results().IncreaseStats().ChangeCount()/total);
	}

	AddChipSum(L"#MI", report_data->Results().IncreaseStats().ModerateCount());
	if (total > 0)
	{
		AddChipSum(L"%MI", 100.0f*report_data->Results().IncreaseStats().ModerateCount()/total);
	}

	AddChipSum(L"#D", report_data->Results().DecreaseStats().ChangeCount());
	if (total > 0)
	{
		AddChipSum(L"%D", 100.0f*report_data->Results().DecreaseStats().ChangeCount()/total);
	}

	AddChipSum(L"#MD", report_data->Results().DecreaseStats().ModerateCount());
	if (total > 0)
	{
		AddChipSum(L"%MD", 100.0f*report_data->Results().DecreaseStats().ModerateCount()/total);
	}

	AddChipSum(L"#NC", report_data->Results().NoChangeStats().ChangeCount());
	if (total > 0)
	{
		AddChipSum(L"%NC", 100.0f*report_data->Results().NoChangeStats().ChangeCount()/total);
	}

	AddChipSum(L"#(A/M->P, MI/I)", report_data->Results().IncreaseStats().DetectionChangeCount());
	if (total > 0)
	{
		AddChipSum(L"%(A/M->P, MI/I)", 100.0f*report_data->Results().IncreaseStats().DetectionChangeCount()/total);
	}

	AddChipSum(L"#(P->A/M, MD/D)",  report_data->Results().DecreaseStats().DetectionChangeCount());
	if (total > 0)
	{
		AddChipSum(L"%(P->A/M, MD/D)", 100.0f*report_data->Results().DecreaseStats().DetectionChangeCount()/total);
	}

	AddChipSum(L"#(P->P, MI/I)", report_data->Results().IncreaseStats().DetectionPresentCount());
	if (total > 0)
	{
		AddChipSum(L"%(P->P, MI/I)", 100.0f*report_data->Results().IncreaseStats().DetectionPresentCount()/total);
	}

	AddChipSum(L"#(P->P, MD/D)", report_data->Results().DecreaseStats().DetectionPresentCount());
	if (total > 0)
	{
		AddChipSum(L"%(P->P, MD/D)", 100.0f*report_data->Results().DecreaseStats().DetectionPresentCount()/total);
	}

	AddChipSum(L"#(A/M->A/M, MI/I)", report_data->Results().IncreaseStats().DetectionAbsentCount());
	if (total > 0)
	{
		AddChipSum(L"%(A/M->A/M, MI/I)", 100.0f*report_data->Results().IncreaseStats().DetectionAbsentCount()/total);
	}

	AddChipSum(L"#(A/M->A/M, MD/D)", report_data->Results().DecreaseStats().DetectionAbsentCount());
	if (total > 0)
	{
		AddChipSum(L"%(A/M->A/M, MD/D)", 100.0f*report_data->Results().DecreaseStats().DetectionAbsentCount()/total);
	}

	AddChipSum(L"#(P->P, NC)", report_data->Results().NoChangeStats().DetectionPresentCount());
	if (total > 0)
	{
		AddChipSum(L"%(P->P, NC)", 100.0f*report_data->Results().NoChangeStats().DetectionPresentCount()/total);
	}

	AddChipSum(L"#(A/M->A/M, NC)", report_data->Results().NoChangeStats().DetectionAbsentCount());
	if (total > 0)
	{
		AddChipSum(L"%(A/M->A/M, NC)", 100.0f*report_data->Results().NoChangeStats().DetectionAbsentCount()/total);
	}

	wchar_t buf2[256];
	float lower;
	float upper;
	ChangeStats::GetBinRange(0, lower, upper);
	swprintf(buf2, 256, L"#I(%0.0f<=SLR<%0.0f)", lower, upper);
	AddChipSum(buf2, report_data->Results().IncreaseStats().FoldChangeCount(0));
	for (int i=1; i<NUMBER_FOLD_CHANGE_BINS; i++)
	{
		ChangeStats::GetBinRange(i, lower, upper);
		swprintf(buf2, 256, L"#I(SLR>=%0.0f)", lower);
		AddChipSum(buf2, report_data->Results().IncreaseStats().FoldChangeCount(i));
	}

	ChangeStats::GetBinRange(0, lower, upper);
	swprintf(buf2, 256, L"#D(%0.0f<=SLR<%0.0f)", lower, upper);
	AddChipSum(buf2, report_data->Results().DecreaseStats().FoldChangeCount(0));
	for (int i=1; i<NUMBER_FOLD_CHANGE_BINS; i++)
	{
		ChangeStats::GetBinRange(i, lower, upper);
		swprintf(buf2, 256, L"#D(SLR>=%0.0f)", lower);
		AddChipSum(buf2, report_data->Results().DecreaseStats().FoldChangeCount(i));
	}
}

/*
 * Store the report results to the chp object.
 */
void MAS5CHPUtils::StoreReport()
{
	if (report_data == NULL)
		return;

	StoreReportParameters();
	StoreReportDetectionStats();
	StoreReportControlStats();
	StoreReportChangeStats();
}

/*
 * Save the results to the CHP file.
 */
bool MAS5CHPUtils::SaveCHPFile(const std::string &chpFile, CExpressionAlgorithmImplementation &exp, ExpressionProbeSetReporter &report, bool reportPass, bool saveToLegacyFile)
{
	mas5_data = &exp;
	if (reportPass == true)
		report_data = &report;
	else
		report_data = NULL;

	// Save the results to the CHP file.
	if (saveToLegacyFile == true)
		return SaveGCOSCHPFile(chpFile);
	else
		return SaveCommandConsoleCHPFile(chpFile);
}

/*
 * Save the results to a GCOS CHP file.
 */
bool MAS5CHPUtils::SaveGCOSCHPFile(const std::string &chpFile)
{
	try
	{
		affxchpwriter::CCHPFileWriter writer;
		gcos_data = &writer;
		writer.SetFileName(chpFile.c_str());
		writer.SetAlgorithmName(MAS5_ALG_NAME);
		writer.SetAlgorithmVersion(MAS5_ALG_VERSION);
		writer.SetParentCelFileName(mas5_data->GetCellData().GetFileName().c_str());
		writer.SetProgID(MAS5_PROG_ID);
		int nsets=mas5_data->GetNumResults();
		writer.InitializeForWriting(mas5_data->GetRows(), mas5_data->GetCols(),
			nsets, mas5_data->GetChipType().c_str(),
			affxcdf::ExpressionProbeSetType);

		// Save the algorithm parameters and the chip summary parameters
		StoreAlgParams();
		StoreSummaryParams(false);
		StoreReport();

		// Write the background info.
		AllZonesInfoType zoneInfo = mas5_data->GetBackgroundZoneInfo();
		writer.AddBackgroundInfo(zoneInfo.number_zones, zoneInfo.smooth_factor);
		for (int i=0; i<zoneInfo.number_zones; i++)
		{
			ZoneInfo &z = zoneInfo.pZones[i];
			writer.AddBackgroundZone((int)z.center.x, (int)z.center.y, z.background);
		}

		// Write the results.
		bool hasCompData = (mas5_data->GetCompStatResult(0) != NULL);
		for (int i=0; i<nsets; i++)
		{
			AbsStatExpressionProbeSetResultType *result = mas5_data->GetAbsStatResult(i);
			writer.GetExpressionResults(i)->Signal = result->Signal;
			writer.GetExpressionResults(i)->Detection = result->Detection;
			writer.GetExpressionResults(i)->DetectionPValue = result->DetectionPValue;
			writer.GetExpressionResults(i)->NumPairs = result->NumPairs;
			writer.GetExpressionResults(i)->NumUsedPairs = result->NumUsedPairs;
			writer.GetExpressionResults(i)->m_HasCompResults = hasCompData;

			if (hasCompData == true)
			{
				CompStatExpressionProbeSetResultType *cresult = mas5_data->GetCompStatResult(i);
				writer.GetExpressionResults(i)->Change = cresult->Change;
				writer.GetExpressionResults(i)->ChangePValue = cresult->ChangePValue;
				writer.GetExpressionResults(i)->SignalLogRatio = cresult->SignalLogRatio;
				writer.GetExpressionResults(i)->SignalLogRatioLow = cresult->SignalLogRatioLow;
				writer.GetExpressionResults(i)->SignalLogRatioHigh = cresult->SignalLogRatioHigh;
				writer.GetExpressionResults(i)->NumCommonPairs = cresult->NumCommonPairs;
			}
		}

		// Create the new file.
		if (writer.CreateNewFile() == false)
			return false;

		// Save the results.
		if (writer.Save() == false)
			return false;

		return true;
	}
	catch (...)
	{
		return false;
	}
}

/*
 * Store the headers of the cel and baseline to the CHP file.
 */
void MAS5CHPUtils::StoreParentHeaders()
{
	GenericData *gdata = mas5_data->GetCellData().GetGenericData();
	if (gdata != NULL)
		cc_data->GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr());

	bool hasCompData = (mas5_data->GetCompStatResult(0) != NULL);
	if (hasCompData == true)
	{
		gdata = mas5_data->GetBaselineCellData().GetGenericData();
		if (gdata != NULL)
			cc_data->GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr());
	}
}

/*
 * Save the program information to the AGCC file.
 */
void MAS5CHPUtils::StoreProgramInfo()
{
	ParameterNameValueType p;
	p.SetName(PROGRAM_NAME);
	p.SetValueText(programName);
	cc_data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(p);
	p.SetName(PROGRAM_ID);
	p.SetValueText(programId);
	cc_data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(p);
	p.SetName(PROGRAM_COMPANY);
	p.SetValueText(programCompany);
	cc_data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(p);
}

/*
 * Save the results to a Command Console CHP file.
 */
bool MAS5CHPUtils::SaveCommandConsoleCHPFile(const std::string &chpFile)
{
	// Save the file
	try
	{
		// Determine the max probe set name length
		int maxln = 0;
		int nsets=mas5_data->GetNumResults();
		for (int i=0; i<nsets; i++)
			maxln = max(maxln, (int) mas5_data->GetProbeSetName(i).size());
		int nzones=mas5_data->GetBackgroundZoneInfo().number_zones;
		CHPData data(chpFile, CHP_EXPRESSION_ASSAY_TYPE);
		cc_data = &data;
		bool hasCompData = (mas5_data->GetCompStatResult(0) != NULL);
		data.SetEntryCount(nsets, maxln, hasCompData);
		data.SetBackgroundZoneCnt(nzones);

		data.SetAlgName(StringUtils::ConvertMBSToWCS(MAS5_ALG_NAME));
		data.SetAlgVersion(StringUtils::ConvertMBSToWCS(MAS5_ALG_VERSION));
		data.SetArrayType(StringUtils::ConvertMBSToWCS(mas5_data->GetChipType()));
		data.SetCols(mas5_data->GetCols());
		data.SetRows(mas5_data->GetRows());
		data.SetParentCell(StringUtils::ConvertMBSToWCS(mas5_data->GetCellData().GetFileName()));
		data.SetProgId(StringUtils::ConvertMBSToWCS(MAS5_PROG_ID));

		// Save the program info, algorithm parameters and the chip summary parameters
		StoreProgramInfo();
		StoreAlgParams();
		StoreSummaryParams(true);
		StoreReport();
		StoreParentHeaders();

		// Create a writer.
		CHPFileWriter writer(data);

		// Write the results.
		writer.SeekToDataSet();
		for (int i=0; i<nsets; i++)
		{
			AbsStatExpressionProbeSetResultType *result = mas5_data->GetAbsStatResult(i);
			if (hasCompData == false)
			{
				CHPExpressionEntry e(
					mas5_data->GetProbeSetName(i),
					result->Detection, result->DetectionPValue, result->Signal,
					result->NumPairs, result->NumUsedPairs);
				writer.WriteExpressionEntry(e);
			}
			else
			{
				CompStatExpressionProbeSetResultType *comp = mas5_data->GetCompStatResult(i);
				CHPExpressionEntry e(
					mas5_data->GetProbeSetName(i),
					result->Detection, result->DetectionPValue, result->Signal,
					result->NumPairs, result->NumUsedPairs,
					true,
					comp->Change,
					comp->ChangePValue,
					comp->SignalLogRatio,
					comp->SignalLogRatioLow,
					comp->SignalLogRatioHigh,
					comp->NumCommonPairs);
				writer.WriteExpressionEntry(e);
			}
		}

		// Write the background info.
		writer.SeekToBgSet();
		AllZonesInfoType zoneInfo = mas5_data->GetBackgroundZoneInfo();
		float sf = zoneInfo.smooth_factor;
		for (int i=0; i<nzones; i++)
		{
			ZoneInfo &z = zoneInfo.pZones[i];
			CHPBackgroundZone bgz(z.center.x,z.center.y, z.background,sf);
			writer.WriteBackgroundZone(bgz);
		}
		return true;
	}
	catch(...)
	{
		return false;
	}
}

