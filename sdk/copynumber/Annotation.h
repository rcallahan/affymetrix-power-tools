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

/**
 * @file   Annotation.h
 * @brief  Wrapper for the SQLite Annotation database.
 *
 */

#ifndef _ANNOTATION_H_
#define _ANNOTATION_H_

#include "copynumber/CNLog2RatioData.h"
//
#include "chipstream/ChipLayout.h"
#include "util/AffxSplitArray.h"
//
#include "util/SQLite.h"
//

/**
 *  Annotation - Wrapper class for the SQLite Annotation database.
 */
class Annotation : public SQLiteDatabase
{
public:
    Annotation();
    virtual ~Annotation();

public:
    static void test();
	static bool isCytoScanHD(const std::string& strFileName);

	void loadAnnotationForCopynumber(BaseEngine& engine, bool bNewProbeSets, AffxSplitArray<CNProbeSet> &vProbeSets, AffxArray<CNProbeSet>& arProbeSets);

    int loadAnnotationForAdapterTypeNormTran(const std::string& strSnpFileName, const std::string& strCnFileName, ChipLayout& layout, AffxMultiDimensionalArray<char>& mxAdapterProbeTypes, AffxMultiDimensionalArray<char>& mxAdapterNspTypes, AffxMultiDimensionalArray<char>& mxAdapterStyTypes);
    int loadAnnotationForAdapterTypeNormTran(const std::string& strFileName, ChipLayout& layout, AffxMultiDimensionalArray<char>& mxAdapterProbeTypes, AffxMultiDimensionalArray<char>& mxAdapterNspTypes, AffxMultiDimensionalArray<char>& mxAdapterStyTypes);

protected:
    void newCNProbeSetsSQLite(BaseEngine& engine, AffxSplitArray<CNProbeSet> &vProbeSets);
    void newCNProbeSetsCsv(BaseEngine& engine, AffxSplitArray<CNProbeSet> &vProbeSets);
    void newCNProbeSetsHDF5(BaseEngine& engine, AffxSplitArray<CNProbeSet> &vProbeSets);

    void loadAnnotationSQLite(BaseEngine& engine, AffxArray<CNProbeSet>& arProbeSets);
    void loadAnnotationCsv(BaseEngine& engine, AffxArray<CNProbeSet>& arProbeSets);
    void loadAnnotationHDF5(BaseEngine& engine, AffxArray<CNProbeSet>& arProbeSets);

    char parseAdapterCode(AffxString& str, const AffxString& strEnzyme, AffxArray<AffxString>& vAdapters);
    int parseFragmentLength(AffxString& str);
    void loadProbeSetNamesFromRestrictList(BaseEngine& engine, AffxArray<AffxString>& ar);
    AffxString runTime(const AffxString& strAction, time_t runTime);

protected:
    static std::vector<affymetrix_calvin_parameter::ParameterNameValueType> m_vParams;

public:
    static std::vector<affymetrix_calvin_parameter::ParameterNameValueType>* getParams() {return &m_vParams;}

};

#endif


