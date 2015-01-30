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

#ifndef _CNLog2RatioAdjustmentMethodWaveCorrection_H_
#define _CNLog2RatioAdjustmentMethodWaveCorrection_H_
/**
 * @file CNLog2RatioAdjustmentMethodWaveCorrection.h
 *
 * @brief This header contains the CNLog2RatioAdjustmentMethodWaveCorrection class definition.
 */

#include "copynumber/CNAnalysisMethodLog2Ratio.h"

/**
 * @brief  The WaveCorrection log2 ratio adjustment method.
 *
 */
class CNLog2RatioAdjustmentMethodWaveCorrection : public CNAnalysisMethodLog2Ratio
{
private:
    AffxString m_strGroup;
    int m_iBandwidth;
    int m_iBinCount;
    int m_iWaveCount;
    bool m_bWaveSmooth;

public:
    static std::string getType() {return "wave-correction-log2ratio-adjustment-method";}
    static std::string getDescription() {return "Copynumber Wave Correction Log2Ratio Adjustment";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

public:
    CNLog2RatioAdjustmentMethodWaveCorrection();
    virtual ~CNLog2RatioAdjustmentMethodWaveCorrection() {}

    virtual bool isSegmentTypeAnalysis() {return false;}

    virtual AffxString getName() {return getType();}

    virtual void run();

protected:
    void adjustLog2Ratios(int iAutosomeCount, int iLastAutosomeChromosome, int iType);
    void applyRunningMean(int iType);
};

#endif


