////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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
 * @file   CNAnalysisMethodCovariateLRAdjuster.h
 *
 * @brief  Class for doing covariate-based log2ratio adjustment
 */
#ifndef _CNAnalysisMethodCovariateLRAdjuster_H_
#define _CNAnalysisMethodCovariateLRAdjuster_H_

#include "copynumber/CNAnalysisMethod.h"
//
#include "util/AffxByteArray.h"
#include "util/AffxConv.h"
#include "util/AffxMultiDimensionalArray.h"
#include "util/Convert.h"
#include "util/Util.h"
//
#include <vector>
//


/**
 * @brief CNAnalysisMethodCovariateLRAdjuster - Class for doing covariate-based log2ratio adjustment
 */
class CNAnalysisMethodCovariateLRAdjuster : public CNAnalysisMethod
{
public:
    static std::string getType() {return "covariate-lr-adjuster";}
    static std::string getDescription() {return "Covariate log2ratio Adjuster";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNAnalysisMethodCovariateLRAdjuster();
    ~CNAnalysisMethodCovariateLRAdjuster();

    virtual bool isSegmentTypeAnalysis() {return false;}

    virtual AffxString getName() {return getType();}

    virtual void run();

    static int determineBins(vector<float>& covValues, vector<int>& binning, int covariateIndex);

protected:
    std::vector<int>& getLRCovariates();
};

#endif  //_CNAnalysisMethodCovariateLRAdjuster_H_
