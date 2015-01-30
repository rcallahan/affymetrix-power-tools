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

#ifndef _CNIntensityAdjustmentMethodHighPassFilter_h_
#define _CNIntensityAdjustmentMethodHighPassFilter_h_
/**
 * @file CNIntensityAdjustmentMethodHighPassFilter.h
 *
 * @brief This header contains the CNIntensityAdjustmentMethodHighPassFilter class definition.
 */

#include "copynumber/CNAnalysisMethod.h"
//
#include <utility>
#include <vector>
//

/**
 * @brief  A class for calulating background.
 *
 */
class CNIntensityAdjustmentMethodHighPassFilter : CNAnalysisMethod
{
private:
    int data_block_rows;  // Height / number of rows in each band
    int data_block_cols; // The band spans the entire image
    int mini_block_rows;
    int mini_block_cols;
    double local_smooth_weight;
    double global_smooth_weight;
    double converged;
    bool m_bUseSingleBlock;

public:
    CNIntensityAdjustmentMethodHighPassFilter();
    virtual ~CNIntensityAdjustmentMethodHighPassFilter();

    static std::string getType() { return "high-pass-filter-intensity-adjustment-method"; }
    static std::string getDescription() { return "CopyNumber HighPassFilter"; }
    static std::string getVersion() { return "1.0"; }

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate * newObject(std::map<std::string, std::string> & param);

    virtual AffxString getName() { return getType(); }
    virtual void run();
    virtual bool isSegmentTypeAnalysis() { return false; }

    void pid_lookup(std::vector<std::pair<int,int> > & lookup, const int row_start,
        const int row_stop, const int image_width);

};

#endif // _CNIntensityAdjustmentMethodHighPassFilter_h_


