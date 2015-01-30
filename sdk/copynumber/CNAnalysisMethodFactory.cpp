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
 * @file CNAnalysisMethodFactory.cpp
 *
 * @brief This file contains the CNAnalysisMethodFactory class members.
 */

#include "copynumber/CNAnalysisMethodFactory.h"
//
#include "copynumber/CNAnalysisMethodAllelePeaks.h"
#include "copynumber/CNAnalysisMethodAllelicDifference.h"
#include "copynumber/CNAnalysisMethodAllelicDifferenceCytoScan.h"
#include "copynumber/CNAnalysisMethodCN.h"
#include "copynumber/CNAnalysisMethodCNCyto2.h"
#include "copynumber/CNAnalysisMethodCNCyto2Gender.h"
#include "copynumber/CNAnalysisMethodCNGender.h"
#include "copynumber/CNAnalysisMethodCNNeutralLOH.h"
#include "copynumber/CNAnalysisMethodCancer.h"
#include "copynumber/CNAnalysisMethodGaussianSmooth.h"
#include "copynumber/CNAnalysisMethodGenotype.h"
#include "copynumber/CNAnalysisMethodKernelSmooth.h"
#include "copynumber/CNAnalysisMethodLOH.h"
#include "copynumber/CNAnalysisMethodLOHCyto2.h"
#include "copynumber/CNAnalysisMethodLOHCytoScan.h"
#include "copynumber/CNAnalysisMethodLog2Ratio.h"
#include "copynumber/CNAnalysisMethodLog2RatioCyto2.h"
#include "copynumber/CNAnalysisMethodMosaicism.h"
#include "copynumber/CNAnalysisMethodNoCall.h"
#include "copynumber/CNAnalysisMethodNormalDiploid.h"
#include "copynumber/CNAnalysisMethodSegmentCN.h"
#include "copynumber/CNAnalysisMethodSegmentLOH.h"
#include "copynumber/CNIntensityAdjustmentMethodHighPassFilter.h"
#include "copynumber/CNIntensityAdjustmentMethodPDNN.h"
#include "copynumber/CNLog2RatioAdjustmentMethodWaveCorrection.h"
#include "copynumber/CNLog2RatioAdjustmentMethodHighPassFilter.h"
#include "copynumber/CNReferenceMethodAdditionalWaves.h"
#include "copynumber/CNReferenceMethodPDNN.h"
#include "copynumber/CNReferenceMethodWaveCorrection.h"
#include "copynumber/CNAnalysisMethodCovariateSignalAdjuster.h"
#include "copynumber/CNAnalysisMethodCovariateLRAdjuster.h"
//


/**
 * @brief Constructor. Registers the objects we know how to create and how
 * they describe themselves.
 */
CNAnalysisMethodFactory::CNAnalysisMethodFactory(bool bAdditionalWaves)
{
    if (bAdditionalWaves) // apt-copynumber-wave
    {
        registerMethod(CNReferenceMethodAdditionalWaves::explainSelf(), &CNReferenceMethodAdditionalWaves::newObject);
        return;
    }

    // Reference methods
    registerMethod(CNReferenceMethodPDNN::explainSelf(), &CNReferenceMethodPDNN::newObject);
    registerMethod(CNReferenceMethodWaveCorrection::explainSelf(), &CNReferenceMethodWaveCorrection::newObject);
    registerMethod(CNReferenceMethodAdditionalWaves::explainSelf(), &CNReferenceMethodAdditionalWaves::newObject);

    //Intensity Adjustment methods
    registerMethod(CNIntensityAdjustmentMethodPDNN::explainSelf(), &CNIntensityAdjustmentMethodPDNN::newObject);
    registerMethod(CNIntensityAdjustmentMethodHighPassFilter::explainSelf(), &CNIntensityAdjustmentMethodHighPassFilter::newObject);

    // Log2Ratio Adjustment methods
    registerMethod(CNLog2RatioAdjustmentMethodWaveCorrection::explainSelf(), &CNLog2RatioAdjustmentMethodWaveCorrection::newObject);
    registerMethod(CNLog2RatioAdjustmentMethodHighPassFilter::explainSelf(), &CNLog2RatioAdjustmentMethodHighPassFilter::newObject);

    // Analysis methods
    registerMethod(CNAnalysisMethodCN::explainSelf(), &CNAnalysisMethodCN::newObject);
    registerMethod(CNAnalysisMethodCNCyto2::explainSelf(), &CNAnalysisMethodCNCyto2::newObject);
    registerMethod(CNAnalysisMethodLog2Ratio::explainSelf(), &CNAnalysisMethodLog2Ratio::newObject);
    registerMethod(CNAnalysisMethodLog2RatioCyto2::explainSelf(), &CNAnalysisMethodLog2RatioCyto2::newObject);
    registerMethod(CNAnalysisMethodAllelicDifference::explainSelf(), &CNAnalysisMethodAllelicDifference::newObject);
    registerMethod(CNAnalysisMethodAllelicDifferenceCytoScan::explainSelf(), &CNAnalysisMethodAllelicDifferenceCytoScan::newObject);
    registerMethod(CNAnalysisMethodGaussianSmooth::explainSelf(), &CNAnalysisMethodGaussianSmooth::newObject);
    registerMethod(CNAnalysisMethodGenotype::explainSelf(), &CNAnalysisMethodGenotype::newObject);
    registerMethod(CNAnalysisMethodKernelSmooth::explainSelf(), &CNAnalysisMethodKernelSmooth::newObject);
    registerMethod(CNAnalysisMethodLOH::explainSelf(), &CNAnalysisMethodLOH::newObject);
    registerMethod(CNAnalysisMethodLOHCyto2::explainSelf(), &CNAnalysisMethodLOHCyto2::newObject);
    registerMethod(CNAnalysisMethodLOHCytoScan::explainSelf(), &CNAnalysisMethodLOHCytoScan::newObject);
    registerMethod(CNAnalysisMethodCNNeutralLOH::explainSelf(), &CNAnalysisMethodCNNeutralLOH::newObject);
    registerMethod(CNAnalysisMethodNormalDiploid::explainSelf(), &CNAnalysisMethodNormalDiploid::newObject);
    registerMethod(CNAnalysisMethodMosaicism::explainSelf(), &CNAnalysisMethodMosaicism::newObject);
//    registerMethod(CNAnalysisMethodNoCall::explainSelf(), &CNAnalysisMethodNoCall::newObject);
//    registerMethod(CNAnalysisMethodCancer::explainSelf(), &CNAnalysisMethodCancer::newObject);
    registerMethod(CNAnalysisMethodCNGender::explainSelf(), &CNAnalysisMethodCNGender::newObject);
    registerMethod(CNAnalysisMethodCNCyto2Gender::explainSelf(), &CNAnalysisMethodCNCyto2Gender::newObject);
    registerMethod(CNAnalysisMethodSegmentCN::explainSelf(), &CNAnalysisMethodSegmentCN::newObject);
    registerMethod(CNAnalysisMethodSegmentLOH::explainSelf(), &CNAnalysisMethodSegmentLOH::newObject);
    registerMethod(CNAnalysisMethodAllelePeaks::explainSelf(), &CNAnalysisMethodAllelePeaks::newObject);

        // Chipstream method
        registerMethod(CNAnalysisMethodChipstream::explainSelf(),  &CNAnalysisMethodChipstream::newObject);

    registerMethod(CNAnalysisMethodCovariateSignalAdjuster::explainSelf(), &CNAnalysisMethodCovariateSignalAdjuster::newObject);
    registerMethod(CNAnalysisMethodCovariateLRAdjuster::explainSelf(), &CNAnalysisMethodCovariateLRAdjuster::newObject);

}

/**
 * @brief Constructor
 */
CNAnalysisMethodFactory::~CNAnalysisMethodFactory()
{
}

/**
 * @brief Create a pointer to a new CNAnalysisMethod object as described
 * in the string specification.
 * @param spec - Specification string
 * @return Pointer to new CNAnalysisMethod objects, must be deleted when
 * finished.
 */
CNAnalysisMethod *CNAnalysisMethodFactory::CNAnalysisMethodForString(const std::string &spec)
{
  CNAnalysisMethod *method = NULL;
  SelfCreate *create = NULL;
  create = SelfCreate::selfCreateFromString(spec, m_Docs, m_Creators, "CNAnalysisMethod");
  /* Check class type. */
  if(InstanceOf(create, CNAnalysisMethod))
  {
    method = static_cast<CNAnalysisMethod *>(create);
  }
  else {
    Err::errAbort("Class doesn't appear to be of type CNAnalysisMethod.");
  }
  return method;
}
