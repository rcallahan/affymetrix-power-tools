////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

#include "exp_report/workflow/ExpressionReportWorkflow.h"
//
#include "exp_report/src/ExpressionControlsParameterExtraction.h"
#include "exp_report/src/ReportDataAccessor.h"
//
#include "calvin_files/fusion/src/FusionCELData.h"
#include "calvin_files/fusion/src/FusionCHPData.h"
#include "calvin_files/fusion/src/FusionCHPLegacyData.h"
#include "calvin_files/fusion/src/FusionCHPQuantificationData.h"
#include "calvin_files/fusion/src/FusionCHPQuantificationDetectionData.h"
#include "calvin_files/fusion/src/FusionPSIData.h"
#include "calvin_files/parsers/src/CHPFileReader.h"
#include "calvin_files/parsers/src/CHPQuantificationDetectionFileReader.h"
#include "calvin_files/parsers/src/CHPQuantificationFileReader.h"
#include "calvin_files/utils/src/FileUtils.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "calvin_files/writers/src/CalvinCHPFileWriter.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationDetectionFileWriter.h"
#include "calvin_files/writers/src/CalvinCHPQuantificationFileWriter.h"
#include "file/CDFFileData.h"
#include "util/Fs.h"
#include "util/Util.h"
//
#include <cmath>
#include <cstdio>
#include <sstream>

#ifndef max
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#endif

using namespace ExpressionReport;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_data;
using namespace affymetrix_fusion_io;
using namespace affymetrix_calvin_io;
using namespace affxchp;
using namespace affxchpwriter;
using namespace affxcdf;
using namespace std;

/*
 * Convert to a string.
 */
static string ToString(int value)
{
	ostringstream str;
	str << value;
	return str.str();
}

/*
 * Convert to a string.
 */
static string ToString(float value)
{
	ostringstream str;
	str << value;
	return str.str();
}

/*! Provides data access for CHP files (signal or legacy).
 * This class provides reporting of absolute data only, no comparison.
 * This is the case as it's designed to work solely off of one CHP file, and
 * as such there is no baseline data.
 */
class FusionChipAccessor : public ReportDataAccessor
{
private: // Properties

    /*! Legacy chp file data. */
    FusionCHPLegacyData *legchp;

    /*! Signal chp file data. */
    FusionCHPQuantificationData *sigchp;

    /*! Quant/detection chp file data. */
    FusionCHPQuantificationDetectionData *detchp;

    /*! The CEL file data. */
    FusionCELData *cel;

    /*! The CDF file data. */
    CCDFFileData *cdf;

    /*! A map of probe set id to name. */
    ProbeSetFileEntryMap *probeSetMap;

    /*! Flag indicating detection calls exist. */
    bool detectionCalls;

    /*! The detection threshold. */
    float detectionThreshold;

private: // Functions

    /*! Read the CDF file.
     * @param libPath The path to the library files.
     * @return True if successful.
     */
    bool ReadCDFFile(const char *libPath)
    {
        string path = Fs::join(libPath,GetChipType()+".cdf");
        string pgfpath = Fs::join(libPath,GetChipType()+".pgf");
        if (FileUtils::Exists(path.c_str()) == true && FileUtils::Exists(pgfpath.c_str()) == false)
        {
            cdf = new CCDFFileData();
            cdf->SetFileName(path.c_str());
            return cdf->Read();
        }
        return true;
    }

    /*! Read the CEL file.
     * @param fileName The file name.
     * @return True if successful.
     */
    bool ReadCELFile(const char *fileName)
    {
        cel = new FusionCELData();
        cel->SetFileName(fileName);
        return cel->Read();
    }

    /*! Read the CHP file.
     * @param fileName The file name.
     * @return True if successful.
     */
    bool ReadCHPFile(const char *fileName)
    {
        // Read the CHP file. This function will read any type of CHP file whose parsers (from Fusion)
        // have been compiled and linked into the program.
        FusionCHPData *chp = FusionCHPDataReg::Read(fileName);
        if (chp == NULL)
            return false;

        // The following function will determine if the CHP file read contains "legacy" format data. This
        // can be either a GCOS/XDA file or a Command Console file. The "legacy" format data is that
        // which contains a signal, detection call, detection p-value, probe pairs, probe pairs used and
        // comparison results if a comparison analysis was performed. This may be from the MAS5, RMA
        // or PLIER algorithms. For RMA and PLIER the detection call will be No Call and the p-value will be 0.
        // Note: The file may also contain genotyping results from the GTYPE software. The ExtractData function
        // will perform a check to ensure it is an expression CHP file.
        legchp = FusionCHPLegacyData::FromBase(chp);
        if (legchp != NULL)
        {
            detectionCalls = true;
            return true;
        }

        sigchp = FusionCHPQuantificationData::FromBase(chp);
        if (sigchp != NULL)
        {
            detectionCalls = false;
            return true;
        }

        detchp = FusionCHPQuantificationDetectionData::FromBase(chp);
        if (detchp != NULL)
        {
            detectionCalls = true;
            return true;
        }

        delete chp;
        return false;
    }

    /*! Get the array/chip type. */
    string GetChipType()
    {
        if (legchp != NULL)
            return StringUtils::ConvertWCSToMBS(legchp->GetHeader().GetChipType());
        else if (sigchp != NULL)
            return StringUtils::ConvertWCSToMBS(sigchp->GetArrayType());
        else if (detchp != NULL)
            return StringUtils::ConvertWCSToMBS(detchp->GetArrayType());
        return "";
    }

public: // Functions

    /*! Constructor */
    FusionChipAccessor() { legchp = NULL; sigchp = NULL; detchp = NULL; cel = NULL; cdf = NULL; detectionCalls = false; }

    /*! Close the files. */
    ~FusionChipAccessor() { Close(); }

    /*! Get the type of file being accessed. */
    ExpressionReportWorkflow::ChpFileType FileType()
    {
        if (sigchp != NULL)
            return ExpressionReportWorkflow::CC_QUANTIFICATION;

        else if (detchp != NULL)
            return ExpressionReportWorkflow::CC_QUANTIFICATION_DETECTION;

        else if (legchp != NULL)
        {
            if (legchp->GetGenericData() != NULL)
                return ExpressionReportWorkflow::CC_LEGACY;
            else
                return ExpressionReportWorkflow::XDA;
        }

        return ExpressionReportWorkflow::UNKNOWN_FILE;
    }

    /*! Flag indicating detection calls exist. */
    bool HasDetectionCalls() const { return detectionCalls; }

    /*! A map of probe set id to probe set name. */
    void SetProbeSetMap(ProbeSetFileEntryMap *psmap) { probeSetMap = psmap; }

    /*! The detection threshold. */
    float &DetectionThreshold() { return detectionThreshold; }

    /*! Closes the files. */
    void Close()
    {
        if (cdf != NULL)
        {
            cdf->Close();
            delete cdf;
            cdf = NULL;
        }
        if (cel != NULL)
        {
            cel->Close();
            delete cel;
            cel = NULL;
        }
        if (legchp != NULL)
        {
            delete legchp;
            legchp = NULL;
        }
        if (sigchp != NULL)
        {
            delete sigchp;
            sigchp = NULL;
        }
        if (detchp != NULL)
        {
            delete detchp;
            detchp = NULL;
        }
        detectionCalls = false;
    }

    /*! Read the CHP file.
     * @param fileName The CHP file name.
     * @param libPath The path to the library files.
     * @return True if successful.
     */
    bool ReadFiles(const char *fileName, const char *celFile, const char *libPath)
    {
        Close();

        if (ReadCHPFile(fileName) == false)
            return false;

        if (ReadCELFile(celFile) == false)
            return false;

        if (ReadCDFFile(libPath) == false)
            return false;

        return true;
    }

	/*! Determines if the direction is anti-sense.
	 * @return True if anti-sense.
	 */
	bool IsAntiSense()
	{
		int antiSenseCount = 0;
		int n = GetNumProbeSets();
		for (int i=0; i<n; i++)
		{
			if (IsAntiSense(i) == true)
				++antiSenseCount;
		}
		return (antiSenseCount > n/2);
	}

	/*! Gets the number of expression probe sets. */
	int GetNumProbeSets()
    {
        if (legchp != NULL)
            return legchp->GetHeader().GetNumProbeSets();
        else if (sigchp != NULL)
            return sigchp->GetEntryCount();
        else if (detchp != NULL)
            return detchp->GetEntryCount();
        return 0;
    }

	/*! Gets the number of used probe pairs in a probe set.
	 * @param index The probe set index.
	 * @return The number of probe pairs.
	 */
	int GetNumPairs(int index)
    {
        if (legchp != NULL)
        {
            FusionExpressionProbeSetResults res;
            legchp->GetExpressionResults(index, res);
            return res.GetNumPairs();
        }
        else if (cdf != NULL)
        {
            CCDFProbeSetInformation set;
            cdf->GetProbeSetInformation(index, set);
            return set.GetNumLists();
        }
        else if (detchp != NULL)
            return 999; // this is just so the probe pair threshold test will pass.
        return 0;
    }

	/*! Checks if a probe set is targeting an anti-sense target.
	 * @param index The probe set index.
	 * @return True if the probe set is designed to interrogate an anti-sense target.
	 */
	bool IsAntiSense(int index)
    {
        if (cdf != NULL)
        {
            std::string name = cdf->GetProbeSetName(index);
		    std::string ext = name.substr(name.length()-2, 2);
		    if (ext == "AT" || ext == "at")
			    return true;
            else
                return false;
        }
        else
            return true;
    }

	/*! Gets the probe set signal value
	 * @param index The probe set index.
	 * @return The signal value.
	 */
	float GetSignal(int index)
    {
        if (legchp != NULL)
        {
            FusionExpressionProbeSetResults res;
            legchp->GetExpressionResults(index, res);
            return res.GetSignal();
        }
        else if (sigchp != NULL)
        {
            ProbeSetQuantificationData res;
            sigchp->GetQuantificationEntry(index, res);
            return res.quantification;
        }
        else if (detchp != NULL)
        {
            ProbeSetQuantificationDetectionData res;
            detchp->GetQuantificationDetectionEntry(index, res);
            return res.quantification;
        }
        return 0.0f;
    }

	/*! Gets the probe set detection value
	 * @param index The probe set index.
	 * @return The detection value.
	 */
	DetectionCall GetDetection(int index)
    {
        if (legchp != NULL)
        {
            FusionExpressionProbeSetResults res;
            legchp->GetExpressionResults(index, res);
            DetectionCall det = (DetectionCall) res.GetDetection();
            return det;
        }
        else if (detchp != NULL)
        {
            ProbeSetQuantificationDetectionData res;
            detchp->GetQuantificationDetectionEntry(index, res);
            return (res.pvalue <= detectionThreshold ? ReportDataAccessor::DetectionPresent : ReportDataAccessor::DetectionAbsent);
        }
        return ReportDataAccessor::DetectionNoCall;
    }

	/*! Checks if the data has comparison results.
	 * @return True if comparison data exists.
	 */
	bool HasComparisonData() { return false; }

	/*! Gets the probe set detection value for the baseline file.
	 * @param index The probe set index.
	 * @return The detection value.
	 */
	DetectionCall GetBaselineDetection(int index) { return ReportDataAccessor::DetectionNoCall; }

	/*! Gets the change call.
	 * @param index The probe set index.
	 * @return The change call.
	 */
	ChangeCall GetChange(int index) { return ReportDataAccessor::ChangeNoCall; }

	/*! Gets the signal log ratio.
	 * @param index The probe set index.
	 * @return The SLR value.
	 */
	float GetSignalLogRatio(int index) { return 0.0f; }

    /*! Gets the list of intensities for a QC probe set.
     * @param qctype The type of QC probe set.
     * @return The list of intensities.
     */
    std::vector<float> GetIntensities(QCProbeSetType qctype)
    {
        std::vector<float> intensities;

        if (cdf != NULL)
        {
            CCDFQCProbeSetInformation set;
            GeneChipQCProbeSetType gtype;

            if (qctype == ReportDataAccessor::CentralCrossNegativeProbeSetType)
                gtype = affxcdf::CentralCrossNegativeQCProbeSetType;

            else if (qctype == ReportDataAccessor::CentralCrossPositiveProbeSetType)
                gtype = affxcdf::CentralCrossPositiveQCProbeSetType;

            else if (qctype == ReportDataAccessor::CheckerboardNegativeProbeSetType)
                gtype = affxcdf::CheckerboardNegativeQCProbeSetType;

            else if (qctype == ReportDataAccessor::CheckerboardPositiveProbeSetType)
                gtype = affxcdf::CheckerboardPositiveQCProbeSetType;

            else
                return intensities;

            cdf->GetQCProbeSetInformation(gtype, set);
            int n = set.GetNumCells();
            intensities.resize(n);
            CCDFQCProbeInformation probe;
            for (int i=0; i<n; i++)
            {
                set.GetProbeInformation(i, probe);
                intensities[i] = cel->GetIntensity(probe.GetX(), probe.GetY());
            }
        }
        return intensities;
    }

    /*! Gets the name of a probe set at the given index position. This should only be
     * used for quant/detection CHP files.
     * @param index The probe set index.
     * @return The probe set name.
     */
    std::string GetProbeSetName(int index)
    {
        if (detchp != NULL && probeSetMap != NULL)
        {
            ProbeSetQuantificationDetectionData res;
            detchp->GetQuantificationDetectionEntry(index, res);
            return (*probeSetMap)[(res.id != -1 ? ToString(res.id) : res.name)];
        }
        else if (sigchp != NULL && probeSetMap != NULL)
        {
            ProbeSetQuantificationData res;
            sigchp->GetQuantificationEntry(index, res);
            return (*probeSetMap)[(res.id != -1 ? ToString(res.id) : res.name)];
        }
        return std::string("");
    }

    /*! Gets the name of a probe set at the given index position. This should only be
     * used for quant/detection CHP files.
     * @param index The probe set index.
     * @return The probe set name.
     */
    std::string GetRawProbeSetName(int index)
    {
        if (detchp != NULL)
        {
            ProbeSetQuantificationDetectionData res;
            detchp->GetQuantificationDetectionEntry(index, res);
            return (res.id != -1 ? ToString(res.id) : res.name);
        }
        else if (sigchp != NULL)
        {
            ProbeSetQuantificationData res;
            sigchp->GetQuantificationEntry(index, res);
            return (res.id != -1 ? ToString(res.id) : res.name);
        }
        return std::string("");
    }

	std::wstring GetAlgoName()
	{
		if (detchp != NULL)
		{
			return detchp->GetAlgName();
		}
	}

	bool IsLogScale()
	{
		ParameterNameValueType param;
		ParameterNameValueTypeList params;
		bool isLog = false;
		if (legchp != NULL)
		{
			isLog = false;
		}
		else
		{
			if (detchp != NULL)
			{
				params = detchp->GetAlgParams();
			}
			else if (sigchp != NULL)
			{
				params = sigchp->GetAlgParams();
			}
			for (ParameterNameValueTypeList::const_iterator it=params.begin(); it != params.end(); it++)
			{
				param = *it;
				if (param.GetName() == QUANTIFICATION_SCALE)
				{
					if (param.GetValueText() == QUANTIFICATION_SCALE_LOG)
					{
						isLog = true;
					}
				}
			}
		}
		return isLog;
	}
};

/*
 * Store the extra metrics to the chip summary section.
 */
void ExpressionReportWorkflow::StoreExtraMetrics(const ParameterNameValueTypeList &extraMetrics)
{
    if (ccsig_data != NULL)
    {
        ccsig_data->AddSummaryParams(extraMetrics);
    }
    else if (ccdet_data != NULL)
    {
        ccdet_data->AddSummaryParams(extraMetrics);
    }
	else
    {
        for (ParameterNameValueTypeList::const_iterator it=extraMetrics.begin(); it!=extraMetrics.end(); it++)
        {
            if (ccleg_data != NULL)
            {
                switch (it->GetParameterType())
                {
                case ParameterNameValueType::Int32Type:
                    ccleg_data->AddChipSum(it->GetName(), it->GetValueInt32());
                    break;
                case ParameterNameValueType::Int16Type:
                    ccleg_data->AddChipSum(it->GetName(), (int32_t)it->GetValueInt16());
                    break;
                case ParameterNameValueType::Int8Type:
                    ccleg_data->AddChipSum(it->GetName(), (int32_t)it->GetValueInt8());
                    break;
                case ParameterNameValueType::UInt16Type:
                    ccleg_data->AddChipSum(it->GetName(), (int32_t)it->GetValueUInt16());
                    break;
                case ParameterNameValueType::UInt8Type:
                    ccleg_data->AddChipSum(it->GetName(), (int32_t)it->GetValueUInt8());
                    break;
                case ParameterNameValueType::FloatType:
                    ccleg_data->AddChipSum(it->GetName(), it->GetValueFloat());
                    break;
                case ParameterNameValueType::TextType:
                    ccleg_data->AddChipSum(it->GetName(), it->GetValueText());
                    break;
                default:
                    ccleg_data->AddChipSum(it->GetName(), it->ToString());
                    break;
                }
                ccleg_data->AddChipSum(it->GetName(), it->GetValueFloat());
            }
            else if (gcos_data != NULL)
            {
                gcos_data->AddChipSummaryParameter(StringUtils::ConvertWCSToMBS(it->GetName()).c_str(), StringUtils::ConvertWCSToMBS(it->ToString()).c_str());
            }
        }
    }
}

/*
 * Store the summary parameters to the chp object.
 */
void ExpressionReportWorkflow::AddChipSum(const std::wstring &name, const std::wstring &value)
{
    if (ccsig_data != NULL)
    {
        ParameterNameValueTypeList plist;
        ParameterNameValueType p;
        p.SetName(name);
        p.SetValueText(value);
        plist.push_back(p);
        ccsig_data->AddSummaryParams(plist);
    }
    else if (ccdet_data != NULL)
    {
        ParameterNameValueTypeList plist;
        ParameterNameValueType p;
        p.SetName(name);
        p.SetValueText(value);
        plist.push_back(p);
        ccdet_data->AddSummaryParams(plist);
    }
	else if (ccleg_data != NULL)
		ccleg_data->AddChipSum(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddChipSummaryParameter(StringUtils::ConvertWCSToMBS(name).c_str(), StringUtils::ConvertWCSToMBS(value).c_str());
}

/*
 * Store the summary parameters to the chp object.
 */
void ExpressionReportWorkflow::AddChipSum(const std::wstring &name, int value)
{
	if (ccsig_data != NULL)
    {
        ParameterNameValueTypeList plist;
        ParameterNameValueType p;
        p.SetName(name);
        p.SetValueInt32(value);
        plist.push_back(p);
        ccsig_data->AddSummaryParams(plist);
    }
	else if (ccdet_data != NULL)
    {
        ParameterNameValueTypeList plist;
        ParameterNameValueType p;
        p.SetName(name);
        p.SetValueInt32(value);
        plist.push_back(p);
        ccdet_data->AddSummaryParams(plist);
    }
	else if (ccleg_data != NULL)
		ccleg_data->AddChipSum(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddChipSummaryParameter(StringUtils::ConvertWCSToMBS(name).c_str(), ToString(value).c_str());
}

/*
 * Store the summary parameters to the chp object.
 */
void ExpressionReportWorkflow::AddChipSum(const std::wstring &name, float value)
{
	if (ccsig_data != NULL)
    {
        ParameterNameValueTypeList plist;
        ParameterNameValueType p;
        p.SetName(name);
        p.SetValueFloat(value);
        plist.push_back(p);
        ccsig_data->AddSummaryParams(plist);
    }
    else if (ccdet_data != NULL)
    {
        ParameterNameValueTypeList plist;
        ParameterNameValueType p;
        p.SetName(name);
        p.SetValueFloat(value);
        plist.push_back(p);
        ccdet_data->AddSummaryParams(plist);
    }
	else if (ccleg_data != NULL)
		ccleg_data->AddChipSum(name, value);

	else if (gcos_data != NULL)
		gcos_data->AddChipSummaryParameter(StringUtils::ConvertWCSToMBS(name).c_str(), ToString(value).c_str());
}

/*
 * Store the detection statistics.
 */
void ExpressionReportWorkflow::StoreReportDetectionStats(bool includeMarginal)
{
    float fval;
	int total = report.Results().ProbeSetResults().NumSets();

	int pcount = report.Results().ProbeSetResults().PresentCalls().Count();
	float psignal = report.Results().ProbeSetResults().PresentCalls().Signal();
	AddChipSum(L"#P", pcount);
    fval = (total > 0 ? 100.0f*pcount/total : 0);
	AddChipSum(L"%P", fval);
    fval = (pcount > 0 ? psignal/pcount : 0);
    AddChipSum(L"Signal(P)", fval);

    float msignal = 0.0f;
    if (includeMarginal == true)
    {
	    int mcount = report.Results().ProbeSetResults().MarginalCalls().Count();
	    msignal = report.Results().ProbeSetResults().MarginalCalls().Signal();
	    AddChipSum(L"#M", mcount);
        fval = (total > 0 ? 100.0f*mcount/total : 0);
	    AddChipSum(L"%M", fval);
        fval = (mcount > 0 ? msignal/mcount : 0);
	    AddChipSum(L"Signal(M)", fval);
    }

	int acount = report.Results().ProbeSetResults().AbsentCalls().Count();
	float asignal = report.Results().ProbeSetResults().AbsentCalls().Signal();
	AddChipSum(L"#A", acount);
    fval = (total > 0 ? 100.0f*acount/total : 0);
	AddChipSum(L"%A", fval);
    fval = (acount > 0 ? asignal/acount : 0);
	AddChipSum(L"Signal(A)", fval);

    fval = (total > 0 ? (psignal+msignal+asignal)/total : 0);
	AddChipSum(L"Signal(All)", fval);
}

/*
 * Store the report parameters.
 */
void ExpressionReportWorkflow::StoreReportParameters()
{
	int total = report.Results().ProbeSetResults().NumSets();
	AddChipSum(L"#Probe Sets Exceeding Probe Pair Threshold", total);
	AddChipSum(L"Probe Pair Threshold", report.Results().ProbePairThreshold());
	if (report.Results().AntiSenseControls() == true)
		AddChipSum(L"Control Direction", L"Antisense");
	else
		AddChipSum(L"Control Direction", L"Sense");
}

/*
 * Store a control statistic.
 */
void ExpressionReportWorkflow::StoreReportControlStat(const wchar_t *controlType,
                                                      const ExpressionControlResult &result,
                                                      bool hasDetectionCalls)
{
	if (result.HasControlResult(ExpressionControl::THREE_PRIME_PROBE_SET) == false &&
		result.HasControlResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == false &&
        result.HasControlResult(ExpressionControl::MIDDLE_PROBE_SET) == false)
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
			swprintf(name, 256, L"%s_%s_%s_signal", controlType, StringUtils::ConvertMBSToWCS(result.GetName()).c_str(), ControlValueString[i]);
			AddChipSum(name, result.GetControlSignalResult((ExpressionControl::ControlValueType)i));

            if (hasDetectionCalls == true)
            {
			    swprintf(name, 256, L"%s_%s_%s_detection", controlType, StringUtils::ConvertMBSToWCS(result.GetName()).c_str(), ControlValueString[i]);
			    swprintf(value, 256, L"%s", Detection[result.GetControlDetectionResult((ExpressionControl::ControlValueType)i)]);
			    AddChipSum(name, value);
            }
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
	if (result.HasControlResult(ExpressionControl::THREE_PRIME_PROBE_SET) == true &&
				result.HasControlResult(ExpressionControl::FIVE_PRIME_PROBE_SET) == true)  //if (ratio >= 0)
	{
		swprintf(name, 256, L"%s_%s_3-5-ratio", controlType, StringUtils::ConvertMBSToWCS(result.GetName()).c_str());
		AddChipSum(name, ratio);
	}
}

/*
 * Store the control statistics.
 */
void ExpressionReportWorkflow::StoreReportControlStats(bool hasDetectionCalls)
{
	ExpressionControlResultList::const_iterator it;
	for (it = report.Results().SpikeStats().begin(); it != report.Results().SpikeStats().end(); ++it)
	{
		const ExpressionControlResult &result = *it;
		StoreReportControlStat(L"Spike", result, hasDetectionCalls);
	}
	for (it = report.Results().HousekeepingStats().begin(); it != report.Results().HousekeepingStats().end(); ++it)
	{
		const ExpressionControlResult &result = *it;
		StoreReportControlStat(L"Housekeeping", result, hasDetectionCalls);
	}
}

/*
 * Store the control probe sets into the CHP file.
 */
void ExpressionReportWorkflow::StoreReportControlProbeSets()
{
    wstring name;
    NameAvgCountList::iterator it;
    for (it = report.Results().ControlStats().begin(); it != report.Results().ControlStats().end(); it++)
    {
        name = StringUtils::ConvertMBSToWCS(it->name) + L" Avg";
        AddChipSum(name, it->avg);
        name = StringUtils::ConvertMBSToWCS(it->name) + L" Count";
        AddChipSum(name, it->count);
    }
}

/*
 * Store the report probe set values.
 */
void ExpressionReportWorkflow::StoreReportProbeSetValues()
{
    NameFloatValuePairList::iterator it;
    for (it=report.Results().ProbeSetValues().begin(); it!=report.Results().ProbeSetValues().end(); it++)
    {
        AddChipSum(StringUtils::ConvertMBSToWCS(it->name), it->value);
    }
}

/*
 * Store the report contents into the CHP file.
 */
void ExpressionReportWorkflow::StoreReport(bool hasDetectionCalls, bool includeReportParameters, bool includeMarginal)
{
    if (includeReportParameters == true)
        StoreReportParameters();
    if (hasDetectionCalls == true)
        StoreReportDetectionStats(includeMarginal);
    StoreReportControlProbeSets();
	StoreReportControlStats(hasDetectionCalls);
    StoreReportProbeSetValues();
}

/*
 * Update the CHP file with the results of the report.
 */
bool ExpressionReportWorkflow::UpdateXDAFileWithReport(const char *fileName, bool hasDetectionCalls, const ParameterNameValueTypeList &extraMetrics)
{
    // Read the input file.
    CCHPFileData chp;
    chp.SetFileName(fileName);
    if (chp.Read() == false)
        return false;
    CCHPFileHeader &header =  chp.GetHeader();

    // Create a temp file name.
    string tmpFile = fileName + string(".tmp_report");
    try
    {
        // Create a new file to store the original contents plus the report results.
        gcos_data = new CCHPFileWriter();
        gcos_data->SetAlgorithmName(header.GetAlgName().c_str());
        gcos_data->SetAlgorithmVersion(header.GetAlgVersion().c_str());
        gcos_data->SetParentCelFileName(header.GetParentCellFile().c_str());
        gcos_data->SetProgID(header.GetProgID().c_str());
        for (TagValuePairTypeList::iterator it=header.AlgorithmParameters().begin(); it!=header.AlgorithmParameters().end(); it++)
        {
            gcos_data->AddAlgorithmParameter(it->Tag.c_str(), it->Value.c_str());
        }
        for (TagValuePairTypeList::iterator it=header.SummaryParameters().begin(); it!=header.SummaryParameters().end(); it++)
        {
            gcos_data->AddChipSummaryParameter(it->Tag.c_str(), it->Value.c_str());
        }
        StoreReport(hasDetectionCalls, true, true);
        StoreExtraMetrics(extraMetrics);
        BackgroundZoneInfo &zinfo = header.GetBackgroundZoneInfo();
        gcos_data->AddBackgroundInfo(zinfo.number_zones, zinfo.smooth_factor);
        for (BackgroundZoneTypeList::iterator it=zinfo.zones.begin(); it!=zinfo.zones.end(); it++)
        {
            gcos_data->AddBackgroundZone((int) it->centerx, (int) it->centery, it->background);
        }
        gcos_data->InitializeForWriting(header.GetRows(), header.GetCols(), header.GetNumProbeSets(), header.GetChipType().c_str(), affxcdf::ExpressionProbeSetType, false);
        gcos_data->SetFileName(tmpFile.c_str());
	    if (gcos_data->CreateNewFile() == false)
        {
            throw;
        }
	    if (gcos_data->SaveHeader() == false)
        {
            throw;
        }
        int n = header.GetNumProbeSets();
	    CExpressionProbeSetResults *entry;
	    for (int i = 0; i < n; i++)
	    {
            entry = chp.GetExpressionResults(i);
		    if (gcos_data->SaveExpressionEntry(entry) == false)
                throw;
	    }
	    if (gcos_data->Close() == false)
        {
            throw;
        }
        delete gcos_data;
        gcos_data = NULL;

        // Move the temp file to the original input file.
        if (remove(fileName) != 0)
            return false;
        if (rename(tmpFile.c_str(), fileName) != 0)
            return false;

        return true;
    }
    catch (...)
    {
        if (gcos_data != NULL)
        {
            gcos_data->Close();
            delete gcos_data;
            gcos_data = NULL;
        }
        remove(tmpFile.c_str());
        return false;
    }
}

/*
 * Update the CHP file with the results of the report.
 */
bool ExpressionReportWorkflow::UpdateCCLegacyFileWithReport(const char *fileName, bool hasDetectionCalls, const ParameterNameValueTypeList &extraMetrics)
{
    // Create a temp file name.
    string tmpFile = fileName + string(".tmp_report");
    CHPData *chp = new CHPData;
    try
    {
        // Read the input file.
        CHPFileReader reader;
        reader.SetFilename(fileName);
        reader.Read(*chp);
        int nsets = chp->GetEntryCount();
        int nz = chp->GetBackgroundZoneCnt();

        // Determine the max probe set name length
        int maxln = -1;
        bool hasCompData = false;
        CHPExpressionEntry entry;
        for (int iset=0; iset<nsets; iset++)
		{
            chp->GetEntry(iset, entry);
            hasCompData = entry.GetHasComparisonData();
            maxln = max(maxln, (int) entry.GetProbeSetName().length());
        }

        // Create a new file to store the report results.
        ccleg_data = new CHPData(tmpFile, chp->GetAssayType());
        ccleg_data->SetBackgroundZoneCnt(nz);
        ccleg_data->SetEntryCount(nsets, maxln, hasCompData);
        int nparents = chp->GetGenericData().Header().GetGenericDataHdr()->GetParentCnt();
        for (int iparent=0; iparent<nparents; ++iparent)
        {
            ccleg_data->GetGenericData().Header().GetGenericDataHdr()->AddParent(
                chp->GetGenericData().Header().GetGenericDataHdr()->GetParent(iparent));
        }
        int nparams = chp->GetGenericData().Header().GetGenericDataHdr()->GetNameValParamCnt();
        for (int iparam=0; iparam<nparams; ++iparam)
        {
            ccleg_data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(
                chp->GetGenericData().Header().GetGenericDataHdr()->GetNameValParam(iparam));
        }

        // Add the report results to the chp object.
        StoreReport(hasDetectionCalls, true, true);
        StoreExtraMetrics(extraMetrics);

        // Create a new file.
        CHPFileWriter *writer = new CHPFileWriter(*ccleg_data);
        writer->SeekToDataSet();
		for (int iset=0; iset<nsets; iset++)
		{
            chp->GetEntry(iset, entry);
			writer->WriteExpressionEntry(entry);
        }

		// Write the background info.
		writer->SeekToBgSet();
        CHPBackgroundZone zone;
        for (int iz=0; iz<nz; iz++)
        {
            chp->GetBackgroundZone(iz, zone);
			writer->WriteBackgroundZone(zone);
        }
        ccleg_data->Clear();
        delete ccleg_data;
        ccleg_data = NULL;

        chp->Clear();
        delete chp;
        chp = NULL;

        delete writer;
        writer = NULL;

        // Move the temp file to the original input file.
        if (remove(fileName) != 0)
            return false;
        if (rename(tmpFile.c_str(), fileName) != 0)
            return false;

        return true;
    }
    catch (...)
    {
        if (ccleg_data != NULL)
        {
            ccleg_data->Clear();
            delete ccleg_data;
            ccleg_data = NULL;
        }
        if (chp != NULL)
        {
            chp->Clear();
            delete chp;
            chp = NULL;
        }
        remove(tmpFile.c_str());
        return false;
    }
}

/*
 * Update the CHP file with the results of the report.
 */
bool ExpressionReportWorkflow::UpdateCCQuantificationFileWithReport(const char *fileName, const ParameterNameValueTypeList &extraMetrics)
{
    // Create a temp file name.
    string tmpFile = fileName + string(".tmp_report");
    CHPQuantificationData *chp = new CHPQuantificationData;
    try
    {
        // Read the input file.
        CHPQuantificationFileReader reader;
        reader.SetFilename(fileName);
        reader.Read(*chp);
        int nsets = chp->GetEntryCount();

        // Determine if probe set names or id's were set
        int maxln = -1;
        ProbeSetQuantificationData entry;
        for (int iset=0; iset<nsets; iset++)
		{
            chp->GetQuantificationEntry(iset, entry);
            if (entry.name.length() == 0)
                break;
            else
                maxln = max(maxln, (int) entry.name.length());
        }

        // Create a new file to store the report results.
        ccsig_data = new CHPQuantificationData(tmpFile);
        int nparents = chp->GetGenericData().Header().GetGenericDataHdr()->GetParentCnt();
        for (int iparent=0; iparent<nparents; ++iparent)
        {
            ccsig_data->GetGenericData().Header().GetGenericDataHdr()->AddParent(
                chp->GetGenericData().Header().GetGenericDataHdr()->GetParent(iparent));
        }
        int nparams = chp->GetGenericData().Header().GetGenericDataHdr()->GetNameValParamCnt();
        for (int iparam=0; iparam<nparams; ++iparam)
        {
            ccsig_data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(
                chp->GetGenericData().Header().GetGenericDataHdr()->GetNameValParam(iparam));
        }
        if (maxln == -1)
            ccsig_data->SetEntryCount(nsets);
        else
            ccsig_data->SetEntryCount(nsets, maxln);

        // Add the report results to the chp object.
        StoreReport(false, false, false);
        StoreExtraMetrics(extraMetrics);

        // Create a new file.
        CHPQuantificationFileWriter *writer = new CHPQuantificationFileWriter(*ccsig_data);
        writer->SeekToDataSet();
		for (int iset=0; iset<nsets; iset++)
		{
            chp->GetQuantificationEntry(iset, entry);
            writer->WriteEntry(entry);
        }
        ccsig_data->Clear();
        delete ccsig_data;
        ccsig_data = NULL;

        chp->Clear();
        delete chp;
        chp = NULL;

        delete writer;
        writer = NULL;

        // Move the temp file to the original input file.
        if (remove(fileName) != 0)
            return false;
        if (rename(tmpFile.c_str(), fileName) != 0)
            return false;

        return true;
    }
    catch (...)
    {
        if (ccsig_data != NULL)
        {
            ccsig_data->Clear();
            delete ccsig_data;
            ccsig_data = NULL;
        }
        if (chp != NULL)
        {
            chp->Clear();
            delete chp;
            chp = NULL;
        }
        remove(tmpFile.c_str());
        return false;
    }
}

/*
 * Update the CHP file with the results of the report.
 */
bool ExpressionReportWorkflow::UpdateCCQuantificationDetectionFileWithReport(const char *fileName, const ParameterNameValueTypeList &extraMetrics)
{
    // Create a temp file name.
    string tmpFile = fileName + string(".tmp_report");
    CHPQuantificationDetectionData *chp = new CHPQuantificationDetectionData;
    try
    {
        // Read the input file.
        CHPQuantificationDetectionFileReader reader;
        reader.SetFilename(fileName);
        reader.Read(*chp);
        int nsets = chp->GetEntryCount();

        // Determine if probe set names or id's were set
        int maxln = -1;
        ProbeSetQuantificationDetectionData entry;
        for (int iset=0; iset<nsets; iset++)
		{
            chp->GetQuantificationDetectionEntry(iset, entry);
            if (entry.name.length() == 0)
                break;
            else
                maxln = max(maxln, (int) entry.name.length());
        }

        // Create a new file to store the report results.
        ccdet_data = new CHPQuantificationDetectionData(tmpFile);
        int nparents = chp->GetGenericData().Header().GetGenericDataHdr()->GetParentCnt();
        for (int iparent=0; iparent<nparents; ++iparent)
        {
            ccdet_data->GetGenericData().Header().GetGenericDataHdr()->AddParent(
                chp->GetGenericData().Header().GetGenericDataHdr()->GetParent(iparent));
        }
        int nparams = chp->GetGenericData().Header().GetGenericDataHdr()->GetNameValParamCnt();
        for (int iparam=0; iparam<nparams; ++iparam)
        {
            ccdet_data->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(
                chp->GetGenericData().Header().GetGenericDataHdr()->GetNameValParam(iparam));
        }
        if (maxln == -1)
            ccdet_data->SetEntryCount(nsets);
        else
            ccdet_data->SetEntryCount(nsets, maxln);

        // Add the report results to the chp object.
        StoreReport(true, false, false);
        StoreExtraMetrics(extraMetrics);

        // Create a new file.
        CHPQuantificationDetectionFileWriter *writer = new CHPQuantificationDetectionFileWriter(*ccdet_data);
        writer->SeekToDataSet();
		for (int iset=0; iset<nsets; iset++)
		{
            chp->GetQuantificationDetectionEntry(iset, entry);
            writer->WriteEntry(entry);
        }
        ccdet_data->Clear();
        delete ccdet_data;
        ccdet_data = NULL;

        chp->Clear();
        delete chp;
        chp = NULL;

        delete writer;
        writer = NULL;

        // Move the temp file to the original input file.
        if (remove(fileName) != 0)
            return false;
        if (rename(tmpFile.c_str(), fileName) != 0)
            return false;

        return true;
    }
    catch (...)
    {
        if (ccdet_data != NULL)
        {
            ccdet_data->Clear();
            delete ccdet_data;
            ccdet_data = NULL;
        }
        if (chp != NULL)
        {
            chp->Clear();
            delete chp;
            chp = NULL;
        }
        remove(tmpFile.c_str());
        return false;
    }
}

/*
 * Update the CHP file with the results of the report.
 */
bool ExpressionReportWorkflow::UpdateFileWithReport(const char *fileName, ChpFileType fileType, bool hasDetectionCalls, const ParameterNameValueTypeList &extraMetrics)
{
    if (fileType == ExpressionReportWorkflow::XDA)
        return UpdateXDAFileWithReport(fileName, hasDetectionCalls, extraMetrics);

    else if (fileType == ExpressionReportWorkflow::CC_LEGACY)
        return UpdateCCLegacyFileWithReport(fileName, hasDetectionCalls, extraMetrics);

    else if (fileType == ExpressionReportWorkflow::CC_QUANTIFICATION)
        return UpdateCCQuantificationFileWithReport(fileName, extraMetrics);

    else if (fileType == ExpressionReportWorkflow::CC_QUANTIFICATION_DETECTION)
        return UpdateCCQuantificationDetectionFileWithReport(fileName, extraMetrics);

    return false;
}

/*
 * Compute the probe set statistics.
 */
bool ExpressionReportWorkflow::Run(const char *fileName, const char *celFile, const char *libPath, const char *controlFile, const char *probeSetFile, const ParameterNameValueTypeList &extraMetrics, bool includeAllProbesets)
{
    // Read the CHP file
    FusionChipAccessor dataAccessor;
    if (dataAccessor.ReadFiles(fileName, celFile, libPath) == false)
        return false;
    dataAccessor.DetectionThreshold() = detectionThreshold;

    // Read the controls file.
    if (paramFilesRead == false)
    {
        paramFilesRead = true;
        probePairThr = 1;
        controls.Clear();
        if (controlFile != NULL && FileUtils::Exists(controlFile) == true)
        {
            if (ExpressionControlsParameterExtraction::ExtractParameters(controlFile, probePairThr, controls) == false)
                return false;
        }

        // Read the probe set indicies from the probe set file. This file is a TSV file with ID and probe set names.
        probeSetIndicies.clear();
        if (probeSetFile != NULL && FileUtils::Exists(probeSetFile) == true)
        {
            if (ExpressionProbeSetFileExtraction::ExtractParameters(probeSetFile, probeSetMap) == false)
                return false;

            // Create a list of probe set indicies for the probe id's found in the file.
            string name;
            int nSets = dataAccessor.GetNumProbeSets();
            ProbeSetFileEntryMap::iterator pos;
            for (int iSet=0; iSet<nSets; iSet++)
            {
                name = dataAccessor.GetRawProbeSetName(iSet);
                pos = probeSetMap.find(name);
                if (pos != probeSetMap.end())
                {
                   probeSetIndicies.push_back(iSet);
                }
            }
        }
    }
    dataAccessor.SetProbeSetMap(&probeSetMap);

    // Run the report.
    report.Clear();
    ChpFileType fileType = dataAccessor.FileType();
    bool hasDetectionCalls = dataAccessor.HasDetectionCalls();
	bool isLog = dataAccessor.IsLogScale();
    bool status = report.Run(dataAccessor.IsAntiSense(), probePairThr, &dataAccessor, &controls, probeSetIndicies, isLog, includeAllProbesets);
    dataAccessor.Close();

    // Update the CHP file with the report results.
    if (status == true)
        status = UpdateFileWithReport(fileName, fileType, hasDetectionCalls, extraMetrics);

    return status;
}

/*
 * Initialize the memebers.
 */
ExpressionReportWorkflow::ExpressionReportWorkflow()
{
    paramFilesRead = false;
    gcos_data = NULL;
    ccleg_data = NULL;
    ccsig_data = NULL;
    ccdet_data = NULL;
    detectionThreshold = 0.0f;
}
