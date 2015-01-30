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

#ifndef _CNLog2RatioAdjustmentMethodHighPassFilter_h_
#define _CNLogwRatioAdjustmentMethodHighPassFilter_h_
/**
 * @file CNLog2RatioAdjustmentMethodHighPassFilter.h
 * @brief This header contains the CNLog2RatioAdjustmentMethodHighPassFilter
 * class definition.
 *
 */

#include "copynumber/CNAnalysisMethod.h"
//
#include <utility>
#include <vector>
//

// Local to this module
typedef struct probe_cell_info_t {
    probe_cell_info_t(): row(0), col(0), probe_id(0),
    probe_set_idx(0), autosome(false) { }

    int row, col, probe_id, probe_set_idx;
    bool autosome;
    } ProbeCellInfo;

/**
 * @brief  A class for calulating background.
 *
 */
class CNLog2RatioAdjustmentMethodHighPassFilter : CNAnalysisMethod
{
private:
    int mini_block_rows;
    int mini_block_cols;
    double local_smooth_weight;
    double global_smooth_weight;
    double converged;
	bool use_hpf;
public:
    CNLog2RatioAdjustmentMethodHighPassFilter();
    virtual ~CNLog2RatioAdjustmentMethodHighPassFilter();

    static std::string getType() {
		return "log2ratio-adjustment-method-high-pass-filter";
		}

    static std::string getDescription() {
		return "CopyNumber HighPassFilter Log2Ratio Adjustment";
		}

    static std::string getVersion() { return "1.0"; }

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate * newObject(std::map<std::string, std::string> & param);

    virtual AffxString getName() { return getType(); }
    virtual void run();
    virtual bool isSegmentTypeAnalysis() { return false; }

};

#endif // _CNLog2RatioAdjustmentMethodHighPassFilter_h_


