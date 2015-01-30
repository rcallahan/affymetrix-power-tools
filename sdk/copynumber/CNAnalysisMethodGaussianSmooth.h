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

#ifndef _CNAnalysisMethodGaussianSmooth_H_
#define _CNAnalysisMethodGaussianSmooth_H_
/**
 * @file CNAnalysisMethodGaussianSmooth.h
 *
 * @brief This header contains the CNAnalysisMethodGaussianSmooth class definition.
 */

#include "copynumber/CNAnalysisMethod.h"
//
#include "util/AffxMultiDimensionalArray.h"
//

#define GSERR_CHECK(x, y) if ( !(x) ) throw(Except(y))

/**
 * @brief  The log2 ration gaussian smooth analysis method.
 *
 */
class CNAnalysisMethodGaussianSmooth : public CNAnalysisMethod
{
private:
    //Gaussian Smooth
    int m_iExpSmoothSignal;
    float m_fSmoothSigmaMultiplier;
    int m_iXChromosome;
  //
  int m_smooth_bw;

public:
    CNAnalysisMethodGaussianSmooth();
    virtual ~CNAnalysisMethodGaussianSmooth() {}

    static std::string getType() {return "gaussian-smooth";}
    static std::string getDescription() {return "CopyNumber GaussianSmooth";}
    static std::string getVersion() {return "1.0";}

    static SelfDoc explainSelf();
    static std::vector<SelfDoc::Opt> getDefaultDocOptions();
    static SelfCreate* newObject(std::map<std::string, std::string>& param);

    virtual AffxString getName() {return getType();}
    virtual void run();
    virtual bool isSegmentTypeAnalysis() {return false;}

private:
    void calculateGaussianSmooth(int iChromosome,
            std::vector<int> & vPositions,
            std::vector<float> & vLog2Ratios);

    void smooth( std::vector<float>& lrvals,
            const std::vector<int>& position,
            int smooth_bw, float sigmaMultiplier, int startidx, int sz);

    void normalizeGaussianSmooth();
};

#endif


