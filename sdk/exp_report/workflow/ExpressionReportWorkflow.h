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

#ifndef _ExpressionReportWorkflow_HEADER_
#define _ExpressionReportWorkflow_HEADER_

/*! \file ExpressionReportWorkflow.h Defines classes to execute the expression report on a CHP file (updating it with the results). */

#include "exp_report/src/ExpressionProbeSetFileExtraction.h"
#include "exp_report/src/ExpressionProbeSetReporter.h"
#include "exp_report/src/ExpressionReportControls.h"
#include "exp_report/src/ExpressionReportData.h"
//
#include "calvin_files/data/src/CHPData.h"
#include "calvin_files/data/src/CHPQuantificationData.h"
#include "calvin_files/data/src/CHPQuantificationDetectionData.h"
#include "calvin_files/parameter/src/ParameterNameValueType.h"
#include "file/CHPFileWriter.h"
//
#include <cstring>
#include <string>
//

/*! Defines the quantification scale parameter name. */
#define QUANTIFICATION_SCALE L"quantification-scale"

/*! Defines the quantification scale parameter value log2. */
#define QUANTIFICATION_SCALE_LOG L"log2"

namespace ExpressionReport
{

/*! Interfaces to execute the expression report on a CHP file (updating it with the results). */
class ExpressionReportWorkflow
{
public:

    /*! The type of CHP file. */
    typedef enum _ChpFileType { XDA, CC_LEGACY, CC_QUANTIFICATION, CC_QUANTIFICATION_DETECTION, UNKNOWN_FILE } ChpFileType;

private:
	/*! The command console CHP data object. */
	affymetrix_calvin_io::CHPData *ccleg_data;

	/*! The command console CHP data object. */
	affymetrix_calvin_io::CHPQuantificationData *ccsig_data;

	/*! The command console CHP data object. */
    affymetrix_calvin_io::CHPQuantificationDetectionData *ccdet_data;

	/*! The GCOS CHP data object. */
	affxchpwriter::CCHPFileWriter *gcos_data;

    /*! The report results. */
    ExpressionProbeSetReporter report;

    /*! The detection threshold. */
    float detectionThreshold;

    /*! The probe pair threshold. */
    int probePairThr;

    /*! The controls. */
    ExpressionControls controls;

    /*! A map of probe set id to name. */
    ProbeSetFileEntryMap probeSetMap;

    /*! A list of probe sets to include in the report. */
    std::list<int> probeSetIndicies;

    /*! A flag to indicate that the param files have been read. */
    bool paramFilesRead;

    /*! Update the CHP file with the results of the report.
     * @param fileName The name of the CHP file.
     * @param fileType The type of CHP file.
     * @param hasDetectionCalls Flag indicating that detection calls exist.
     * @param extraMetrics A list of extra report metrics to add to the CHP file.
     * @return True if successful.
     */
    bool UpdateFileWithReport(const char *fileName, ChpFileType fileType, bool hasDetectionCalls, const ParameterNameValueTypeList &extraMetrics);

    /*! Update the CHP file with the results of the report.
     * @param fileName The name of the CHP file.
     * @param hasDetectionCalls Flag indicating detection calls exist.
     * @param extraMetrics A list of extra report metrics to add to the CHP file.
     * @return True if successful.
     */
    bool UpdateCCLegacyFileWithReport(const char *fileName, bool hasDetectionCalls, const ParameterNameValueTypeList &extraMetrics);

    /*! Update the CHP file with the results of the report.
     * @param fileName The name of the CHP file.
     * @param extraMetrics A list of extra report metrics to add to the CHP file.
     * @return True if successful.
     */
    bool UpdateCCQuantificationFileWithReport(const char *fileName, const ParameterNameValueTypeList &extraMetrics);

    /*! Update the CHP file with the results of the report.
     * @param fileName The name of the CHP file.
     * @param extraMetrics A list of extra report metrics to add to the CHP file.
     * @return True if successful.
     */
    bool UpdateCCQuantificationDetectionFileWithReport(const char *fileName, const ParameterNameValueTypeList &extraMetrics);

    /*! Update the CHP file with the results of the report.
     * @param fileName The name of the CHP file.
     * @param hasDetectionCalls Flag indicating detection calls exist.
     * @param extraMetrics A list of extra report metrics to add to the CHP file.
     * @return True if successful.
     */
    bool UpdateXDAFileWithReport(const char *fileName, bool hasDetectionCalls, const ParameterNameValueTypeList &extraMetrics);

    /*! Stores the report contents into the CHP file.
     * @param hasDetectionCalls Flag indicating detection calls exist.
     * @param includeReportParameters True if the report parameters are to be included.
     * @param includeMarginal True if marginal call stats are to be included.
     */
    void StoreReport(bool hasDetectionCalls, bool includeReportParameters, bool includeMarginal);

    /*! Stores the extra metrics to the CHP file header.
     * @param extraMetrics A list of extra report metrics to add to the CHP file.
     */
    void StoreExtraMetrics(const ParameterNameValueTypeList &extraMetrics);

    /*! Stores the detection statistics.
     * @param includeMarginal True if marginal call stats are to be included.
     */
    void StoreReportDetectionStats(bool includeMarginal);

    /*! Stores the control statistics.
     * @param hasDetectionCalls Flag indicating detection calls exist.
     */
    void StoreReportControlStats(bool hasDetectionCalls);

    /*! Stores the probe set values from the report. */
    void StoreReportProbeSetValues();

    /*! Stores the control probe sets statistics. */
    void StoreReportControlProbeSets();

    /*! Stores the report parameters. */
    void StoreReportParameters();

    /*! Store a control statistic.
     * @param controlType The name of the control.
     * @param result The control result.
     * @param hasDetectionCalls Flag indicating detection calls exist.
     */
    void StoreReportControlStat(const wchar_t *controlType,
                                const ExpressionControlResult &result,
                                bool hasDetectionCalls);

    /*! Store the summary parameters to the chp object.
     * @param name The parameter name.
     * @param value The parameter value.
     */
    void AddChipSum(const std::wstring &name, const std::wstring &value);

    /*! Store the summary parameters to the chp object.
     * @param name The parameter name.
     * @param value The parameter value.
     */
    void AddChipSum(const std::wstring &name, int value);

    /*! Store the summary parameters to the chp object.
     * @param name The parameter name.
     * @param value The parameter value.
     */
    void AddChipSum(const std::wstring &name, float value);

public:

    /*! The detection threshold. */
    float &DetectionThreshold() { return detectionThreshold; }

	/*! Run the report on the input file.
	 * @param fileName The file name of the CHP file.
     * @param celFile The name of the CEL file.
     * @param libPath The path to the library files.
     * @param controlFile The file name of the report controls file.
     * @param probeSetFile The file name of the file containing a list of probe sets to include in the report.
     * @param extraMetrics A list of extra report metrics to add to the CHP file.
	 * @param includeAllProbesets Include all probe sets in the report
	 * @return True if successful.
	 */
	bool Run(const char *fileName, const char *celFile, const char *libPath, const char *controlFile, const char *probeSetFile, const ParameterNameValueTypeList &extraMetrics, bool includeAllProbesets);

    /*! Constructor */
    ExpressionReportWorkflow();
};

}

#endif
