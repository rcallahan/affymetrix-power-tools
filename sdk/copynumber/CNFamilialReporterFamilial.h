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

#ifndef _CNFamilialReporterFamilial_H_
#define _CNFamilialReporterFamilial_H_
/**
 * @file CNFamilialReporterFamilial.h
 *
 * @brief This header contains the CNFamilialReporterFamilial class definition.
 */

#include "copynumber/CNFamilialReporter.h"
//
#include "util/PgOptions.h"
//

/**
 * @brief  The Familial file Familial Reporter.
 *
 */
class CNFamilialReporterFamilial : public CNFamilialReporter
{
public:
    static std::string getType() {return "familial-output";}
    static std::string getDescription() {return "Familial Output";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNFamilialReporterFamilial();
    virtual ~CNFamilialReporterFamilial();

    virtual void run();

protected:
    void setupHeader(const AffxString& strFileName, affymetrix_calvin_io::GenericData& data);
    void addParentHeaders(affymetrix_calvin_io::GenericData& data);
    void addStandardHeader(affymetrix_calvin_io::GenericDataHeader* pHeader, BaseEngine* pEngine);
    void addPgOptions(affymetrix_calvin_io::GenericDataHeader* pHeader, PgOptions& opt, const AffxString& strPrefix);
    void loadParam(const AffxString& strName, PgOpt::PgOptType type, const AffxString& strValue, affymetrix_calvin_parameter::ParameterNameValueType& param);

    void setupDataSetHeaderSamples(affymetrix_calvin_io::GenericData& data, int iGroupIndex);
    void setupDataSetHeaderPaternityTesting(affymetrix_calvin_io::GenericData& data, int iGroupIndex);
    void setupDataSetHeaderMIE(affymetrix_calvin_io::GenericData& data, int iGroupIndex);
    void setupDataSetHeaderSegmentOverlaps(affymetrix_calvin_io::GenericData& data, CNFamilialAnalysisMethod* pMethod, int iGroupIndex);
    void setupDataSetHeaderSegments(affymetrix_calvin_io::GenericData& data, CNFamilialAnalysisMethod* pMethod, int iGroupIndex);

    void writeDataGroupAndDataSetHeaders(affymetrix_calvin_io::GenericFileWriter* writer);
    int setupDataSet(affymetrix_calvin_io::GenericFileWriter* writer, int iGroupIndex, int iSetIndex);
    void writeDataSetData(affymetrix_calvin_io::GenericFileWriter* writer);
    void writeDataSetSamples(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex);
    void writeDataSetPaternityTesting(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex);
    void writeDataSetMIE(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex);
    void writeDataSetSegmentOverlaps(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex);
    void writeDataSetSegments(affymetrix_calvin_io::GenericFileWriter* writer, CNFamilialAnalysisMethod* pMethod, int iGroupIndex, int iSetIndex);
    void writeStandardSegment(affymetrix_calvin_io::GenericFileWriter* writer, CNSegment* p, int iGroupIndex, int iSetIndex);
};

#endif


