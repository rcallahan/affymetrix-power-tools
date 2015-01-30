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

#ifndef _CNAnalysisMethodLog2RatioCyto2_H_
#define _CNAnalysisMethodLog2RatioCyto2_H_
/**
 * @file CNAnalysisMethodLog2RatioCyto2.h
 *
 * @brief This header contains the CNAnalysisMethodLog2RatioCyto2 class definition.
 */

#include "copynumber/CNAnalysisMethodLog2Ratio.h"

/**
 * @brief  The Log2RatioCyto2 analysis method.
 *
 */
class CNAnalysisMethodLog2RatioCyto2 : public CNAnalysisMethodLog2Ratio
{

private:
    int m_iWaveCorrectionBandwidth;
    int m_iWaveCorrectionBinCount;

public:
    static std::string getType() {return "log2-ratio-cyto2";}
    static std::string getDescription() {return "Copynumber Log2RatioCyto2";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    virtual void run();


protected:
    void applyWaveCorrection();
    void adjustLog2Ratios(int iAutosomeCount, int iLastAutosomeChromosome, std::vector<float>& v);
    void applyRunningMean(std::vector<float>& v);
};

#endif


