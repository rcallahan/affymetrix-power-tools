////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

#include "copynumber/CNCytoEngine.h"
#include "copynumber/CNAnalysisMethodCovariateParams.h"
#include "copynumber/CPPTest/Setup.h" 
//
//#include "util/AffxArray.h"
#include "util/Util.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <iostream>
//
using namespace std;
/**
 * @class CovariateParamsTest
 * @brief cppunit class for testing CovariateParams functions.
 */

class CovariateParamsTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(CovariateParamsTest);
    CPPUNIT_TEST(determineCovariateMapTest);
    CPPUNIT_TEST(checkParamsTest);
    CPPUNIT_TEST(isReservedCovariateNameTest);
    CPPUNIT_TEST_SUITE_END();

public:  
   
   void determineCovariateMapTest();
   void checkParamsTest();
   void isReservedCovariateNameTest();
     
};
CPPUNIT_TEST_SUITE_REGISTRATION(CovariateParamsTest );


void CovariateParamsTest::determineCovariateMapTest()
{
    // signal-adjustment-covariates only
    CNCytoEngine cnCyto1;
    cnCyto1.setOpt("signal-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16,14");

    CovariateParams::determineCovariateMap(cnCyto1);
    CPPUNIT_ASSERT(CovariateParams::m_allCovariateMap.size() == 3);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-adapter-type") == 1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-length") == 0);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("probe-gc") == 2);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("local-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates.size() == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[2] == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(0) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(1) == true);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(2) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(2) == false);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(2) == true);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(0) == 20);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(1) == 16);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(2) == 14);


    // lr-adjustment-covariates only
    CNCytoEngine cnCyto2;
    cnCyto2.setOpt("lr-adjustment-covariates","somename.order=probe-gc,fragment-length,fragment-adapter-type,something-else.bin-type=equally-populated,discrete,equally-spaced,discrete.bin-count=21,17,23,0.iqr-scaling=off,on,off,off.subtract-from-XY=on,off,on,on");

    CovariateParams::determineCovariateMap(cnCyto2);
    CPPUNIT_ASSERT(CovariateParams::m_allCovariateMap.size() == 4);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-adapter-type") == 2);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-length") == 1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("probe-gc") == 0);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("local-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("something-else") == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates.size() == 4);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates[2] == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates[3] == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedLRBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedLRBins(0) == true);
    CPPUNIT_ASSERT(CovariateParams::isLRCovariateDiscrete(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedLRBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedLRBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::isLRCovariateDiscrete(1) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedLRBins(2) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedLRBins(2) == false);
    CPPUNIT_ASSERT(CovariateParams::isLRCovariateDiscrete(2) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedLRBins(3) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedLRBins(3) == false);
    CPPUNIT_ASSERT(CovariateParams::isLRCovariateDiscrete(3) == true);
    CPPUNIT_ASSERT(CovariateParams::getNumLRBins(0) == 21);
    CPPUNIT_ASSERT(CovariateParams::getNumLRBins(1) == 17);
    CPPUNIT_ASSERT(CovariateParams::getNumLRBins(2) == 23);
    CPPUNIT_ASSERT(CovariateParams::getNumLRBins(3) == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling.size() == 4);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[0] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[1] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[2] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[3] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY.size() == 4);
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[0] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[1] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[2] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[3] == "on");

    // Signal and LR of covariates.
    CNCytoEngine cnCyto3;
    cnCyto3.setOpt("signal-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16,14");
    cnCyto3.setOpt("lr-adjustment-covariates","somename.order=fragment-gc,probe-gc.bin-type=equally-populated,discrete.bin-count=21,17.iqr-scaling=on,on,off.subtract-from-XY=on,off,off");

    CovariateParams::determineCovariateMap(cnCyto3);
    CPPUNIT_ASSERT(CovariateParams::m_allCovariateMap.size() == 4);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-adapter-type") == 1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-length") == 0);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-gc") == 3);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("probe-gc") == 2);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("local-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("something-else") == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates.size() == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[2] == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates[0] == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates[1] == 2);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(0) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(1) == true);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(2) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(2) == false);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(2) == true);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(0) == 20);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(1) == 16);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(2) == 14);

    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedLRBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedLRBins(0) == true);
    CPPUNIT_ASSERT(CovariateParams::isLRCovariateDiscrete(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedLRBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedLRBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::isLRCovariateDiscrete(1) == true);
    CPPUNIT_ASSERT(CovariateParams::getNumLRBins(0) == 21);
    CPPUNIT_ASSERT(CovariateParams::getNumLRBins(1) == 17);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling.size() == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[0] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[1] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[2] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY.size() == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[0] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[1] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[2] == "off");

    // LR m_vIQRScaling and m_vLRSubractFromXY are not checked for values.  No exception expected.
    CNCytoEngine cnCyto4;
    cnCyto4.setOpt("lr-adjustment-covariates","somename.order=fragment-gc,probe-gc.bin-type=equally-populated,discrete.bin-count=21,17.iqr-scaling=on,on,nonsense.subtract-from-XY=on,nonsense,off");
    CovariateParams::determineCovariateMap(cnCyto4);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[2] == "nonsense");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[1] == "nonsense");

    // order values are free text.  No exception expected.
    CNCytoEngine cnCyto5;
    cnCyto5.setOpt("signal-adjustment-covariates","somename.order=nonsense,fragment-adapter-type.bin-type=equally-spaced,equally-populated.bin-count=20,16");
    CovariateParams::determineCovariateMap(cnCyto5);

    // Expect an exception - non-numeric bin-count (nonsense)
    CNCytoEngine cnCyto6;
    cnCyto6.setOpt("signal-adjustment-covariates","somename.order=fragment-adapter-type,probe-gc.bin-type=equally-spaced,nonsense.bin-count=20,nonsense");
    NEGATIVE_TEST(CovariateParams::determineCovariateMap(cnCyto6),Except);

    // allele peaks adjustment covariates only
    CNCytoEngine cnCyto7;
    cnCyto7.setOpt("allele-peaks-adjustment-covariates","someothername.order=fragment-length,probe-gc.bin-type=discrete,equally-spaced.bin-count=17,23.coarse-allele-peak-adjustment=off,off");

    CovariateParams::determineCovariateMap(cnCyto7);
    CPPUNIT_ASSERT(CovariateParams::m_allCovariateMap.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-adapter-type") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-length") == 0);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("probe-gc") == 1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("local-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("something-else") == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedAPBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedAPBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::isAPCovariateDiscrete(0) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedAPBins(1) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedAPBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::isAPCovariateDiscrete(1) == false);
    CPPUNIT_ASSERT(CovariateParams::getNumAPBins(0) == 17);
    CPPUNIT_ASSERT(CovariateParams::getNumAPBins(1) == 23);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust[0] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust[1] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustStep.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustWindow.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustPointCount.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustBandwidth.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCutoff.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCleanthreshold.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustOutliertrim.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_masterPeakPointCount == -1);
    CPPUNIT_ASSERT(CovariateParams::m_masterPeakBandwidth == -1);
    CPPUNIT_ASSERT(CovariateParams::m_covariatePeakPointCount == -1);
    CPPUNIT_ASSERT(CovariateParams::m_covariatePeakBandwidth == -1);
    CPPUNIT_ASSERT(CovariateParams::m_kernelFunctionSelection == "");

    // allele peaks adjustment covariates only
    CNCytoEngine cnCyto8;
    cnCyto8.setOpt("allele-peaks-adjustment-covariates","someothername.order=fragment-length,probe-gc.bin-type=discrete,equally-spaced.bin-count=17,23.coarse-allele-peak-adjustment=off,on"
      ".coarse-allele-peak-adjustment-step=-1,2.coarse-allele-peak-adjustment-window=0,100.coarse-allele-peak-adjustment-point-count=-1,1.coarse-allele-peak-adjustment-bandwidth=0,4"
      ".coarse-allele-peak-adjustment-cutoff=-1,1.coarse-allele-peak-adjustment-clean-threshold=0,1.coarse-allele-peak-adjustment-outlier-trim=-1,1"
      ".master-peaks-point-count=1.master-peaks-bandwidth=4.covariate-bin-peaks-point-count=1.covariate-bin-peaks-bandwidth=4.kernel-function-selection=gaussian");

    CovariateParams::determineCovariateMap(cnCyto8);
    CPPUNIT_ASSERT(CovariateParams::m_allCovariateMap.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-adapter-type") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-length") == 0);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("probe-gc") == 1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("local-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("something-else") == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedAPBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedAPBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::isAPCovariateDiscrete(0) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedAPBins(1) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedAPBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::isAPCovariateDiscrete(1) == false);
    CPPUNIT_ASSERT(CovariateParams::getNumAPBins(0) == 17);
    CPPUNIT_ASSERT(CovariateParams::getNumAPBins(1) == 23);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust[0] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust[1] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustStep.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustStep[0] == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustStep[1] == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustWindow.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustWindow[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustWindow[1] == 100);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustPointCount.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustPointCount[0] == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustPointCount[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustBandwidth.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustBandwidth[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustBandwidth[1] == 4);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCutoff.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCutoff[0] == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCutoff[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCleanthreshold.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCleanthreshold[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCleanthreshold[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustOutliertrim.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustOutliertrim[0] == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustOutliertrim[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_masterPeakPointCount == 1);
    CPPUNIT_ASSERT(CovariateParams::m_masterPeakBandwidth == 4);
    CPPUNIT_ASSERT(CovariateParams::m_covariatePeakPointCount == 1);
    CPPUNIT_ASSERT(CovariateParams::m_covariatePeakBandwidth == 4);
    CPPUNIT_ASSERT(CovariateParams::m_kernelFunctionSelection == "gaussian");

    // Signal, LR and allele peak covariates
    CNCytoEngine cnCyto9;
    cnCyto9.setOpt("signal-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16,14");
    cnCyto9.setOpt("lr-adjustment-covariates","somename.order=fragment-gc,probe-gc.bin-type=equally-populated,discrete.bin-count=21,17.iqr-scaling=on,on,off.subtract-from-XY=on,off,off");
    cnCyto9.setOpt("allele-peaks-adjustment-covariates","someothername.order=fragment-length,probe-gc.bin-type=discrete,equally-spaced.bin-count=17,23.coarse-allele-peak-adjustment=off,off");

    CovariateParams::determineCovariateMap(cnCyto9);
    CPPUNIT_ASSERT(CovariateParams::m_allCovariateMap.size() == 4);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-adapter-type") == 1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-length") == 0);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("fragment-gc") == 3);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("probe-gc") == 2);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("local-gc") == -1);
    CPPUNIT_ASSERT(CovariateParams::mapCovariateNameToIndex("something-else") == -1);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates.size() == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[1] == 1);
    CPPUNIT_ASSERT(CovariateParams::m_vSignalCovariates[2] == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates[0] == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vLRCovariates[1] == 2);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(0) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(1) == true);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedSignalBins(2) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedSignalBins(2) == false);
    CPPUNIT_ASSERT(CovariateParams::isSignalCovariateDiscrete(2) == true);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(0) == 20);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(1) == 16);
    CPPUNIT_ASSERT(CovariateParams::getNumSignalBins(2) == 14);

    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedLRBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedLRBins(0) == true);
    CPPUNIT_ASSERT(CovariateParams::isLRCovariateDiscrete(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedLRBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedLRBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::isLRCovariateDiscrete(1) == true);
    CPPUNIT_ASSERT(CovariateParams::getNumLRBins(0) == 21);
    CPPUNIT_ASSERT(CovariateParams::getNumLRBins(1) == 17);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling.size() == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[0] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[1] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vIQRScaling[2] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY.size() == 3);
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[0] == "on");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[1] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vLRSubractFromXY[2] == "off");

    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates[0] == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vAPCovariates[1] == 2);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedAPBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedAPBins(0) == false);
    CPPUNIT_ASSERT(CovariateParams::isAPCovariateDiscrete(0) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallySpacedAPBins(1) == true);
    CPPUNIT_ASSERT(CovariateParams::useEquallyPopulatedAPBins(1) == false);
    CPPUNIT_ASSERT(CovariateParams::isAPCovariateDiscrete(1) == false);
    CPPUNIT_ASSERT(CovariateParams::getNumAPBins(0) == 17);
    CPPUNIT_ASSERT(CovariateParams::getNumAPBins(1) == 23);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust[0] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust[1] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust.size() == 2);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust[0] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjust[1] == "off");
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustStep.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustWindow.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustPointCount.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustBandwidth.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCutoff.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustCleanthreshold.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_vCoarseAPAdjustOutliertrim.size() == 0);
    CPPUNIT_ASSERT(CovariateParams::m_masterPeakPointCount == -1);
    CPPUNIT_ASSERT(CovariateParams::m_masterPeakBandwidth == -1);
    CPPUNIT_ASSERT(CovariateParams::m_covariatePeakPointCount == -1);
    CPPUNIT_ASSERT(CovariateParams::m_covariatePeakBandwidth == -1);
    CPPUNIT_ASSERT(CovariateParams::m_kernelFunctionSelection == "");
}

void CovariateParamsTest::checkParamsTest()
{
    CNCytoEngine cnCyto1;
    CovariateParams::checkParams(cnCyto1); // Doesn't throw exception

    CNCytoEngine cnCyto2;
    cnCyto2.setOpt("signal-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16,14");
    CovariateParams::checkParams(cnCyto2); // Doesn't throw exception

    CNCytoEngine cnCyto3;
    cnCyto3.setOpt("lr-adjustment-covariates","somename.order=fragment-gc,probe-gc.bin-type=equally-populated,discrete.bin-count=21,17.iqr-scaling=off,on.subtract-from-XY=on,off");
    CovariateParams::checkParams(cnCyto3); // Doesn't throw exception

    CNCytoEngine cnCyto4;
    cnCyto4.setOpt("signal-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16,14");
    cnCyto4.setOpt("lr-adjustment-covariates","somename.order=fragment-gc,probe-gc.bin-type=equally-populated,discrete.bin-count=21,17.iqr-scaling=off,on.subtract-from-XY=on,off");
    cnCyto4.setOpt("allele-peaks-adjustment-covariates", "someothername.order=fragment-length,probe-gc.bin-type=discrete,equally-spaced.bin-count=0,23.coarse-allele-peak-adjustment=off,off");
    CovariateParams::checkParams(cnCyto4); // Doesn't throw exception

    // Expect an exception - signal order and bin-type element count mismatch
    CNCytoEngine cnCyto5;
    cnCyto5.setOpt("signal-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16,14");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto5), Except);

    // Expect an exception - signal order and bin-count element count mismatch
    CNCytoEngine cnCyto6;
    cnCyto6.setOpt("signal-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type.bin-type=equally-spaced,equally-populated.bin-count=20,16,14");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto6), Except);

    // Expect an exception - log2 ratio order and bin-type element count mismatch
    CNCytoEngine cnCyto7;
    cnCyto7.setOpt("lr-adjustment-covariates","somename.order=fragment-adapter-type.bin-type=equally-spaced,uniform,equally-populated.bin-count=20,16,14.iqr-scaling=on,on,off.subtract-from-XY=on,off,off");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto7), Except);

    // Expect an exception - log2 ratio order and bin-count element count mismatch
    CNCytoEngine cnCyto8;
    cnCyto8.setOpt("lr-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16.iqr-scaling=on,on,off.subtract-from-XY=on,off,off");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto8), Except);

    // Expect an exception - log2 ratio order and iqr-scaling element count mismatch
    CNCytoEngine cnCyto9;
    cnCyto9.setOpt("lr-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16.iqr-scaling=on,on.subtract-from-XY=on,off,off");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto9), Except);

    // Expect an exception - log2 ratio order and subtract-from-XY element count mismatch
    CNCytoEngine cnCyto10;
    cnCyto10.setOpt("lr-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,16.iqr-scaling=on,on,off.subtract-from-XY=off,off");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto10), Except);

    // Expect an exception - file-based covariates in lr-adjustment-covariates option but no covariates-file option on the command line.
    CNCytoEngine cnCyto11;
    cnCyto11.setOpt("signal-adjustment-covariates","somename.order=nonsense,fragment-adapter-type.bin-type=equally-spaced,equal-number.bin-count=20,16");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto11),Except);

    // Unknown bin-type does not generate an exception (nonsense)
    CNCytoEngine cnCyto12;
    cnCyto12.setOpt("signal-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type.bin-type=equally-spaced,nonsense.bin-count=20,16");
    CovariateParams::checkParams(cnCyto12); // Doesn't throw exception

    // Expect an exception - non-numeric bin-count (nonsense)
    CNCytoEngine cnCyto13;
    cnCyto13.setOpt("signal-adjustment-covariates","somename.order=fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-spaced.bin-count=20,nonsense");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto13),Except);

    // Expect an exception - negative non-discrete signal bin-count
    CNCytoEngine cnCyto14;
    cnCyto14.setOpt("signal-adjustment-covariates","somename.order=fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-spaced.bin-count=20,-2");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto14),Except);

    // Expect an exception - negative non-discrete LR bin-count
    CNCytoEngine cnCyto15;
    cnCyto15.setOpt("lr-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-populated,discrete.bin-count=20,-1,0.iqr-scaling=on,on,off.subtract-from-XY=on,off,off");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto15),Except);

    // No exception - negative discrete signal bin-count
    CNCytoEngine cnCyto16;
    cnCyto16.setOpt("signal-adjustment-covariates","somename.order=fragment-adapter-type,probe-gc.bin-type=equally-spaced,discrete.bin-count=20,-2");
    CovariateParams::checkParams(cnCyto16);

    // No exception - negative discrete LR bin-count
    CNCytoEngine cnCyto17;
    cnCyto17.setOpt("lr-adjustment-covariates","somename.order=fragment-length,fragment-adapter-type,probe-gc.bin-type=equally-spaced,discrete,discrete.bin-count=20,-1,0.iqr-scaling=on,on,off.subtract-from-XY=on,off,off");
    CovariateParams::checkParams(cnCyto17);

    // Expect an exception - zero non-discrete signal bin-count
   	CNCytoEngine cnCyto18;
    cnCyto18.setOpt("signal-adjustment-covariates","somename.order=fragment-adapter-type,probe-gc.bin-type=equally-spaced,equally-spaced.bin-count=20,0");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto18),Except);

    // Expect an exception - allele peak order and bin-type element count mismatch
    CNCytoEngine cnCyto19;
    cnCyto19.setOpt("allele-peaks-adjustment-covariates","someothername.order=fragment-length,probe-gc,fragment-adapter-type.bin-type=discrete,equally-spaced.bin-count=0,23.coarse-allele-peak-adjustment=off,off");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto19), Except);

    // Expect an exception - allele peak order and bin-count element count mismatch
    CNCytoEngine cnCyto20;
    cnCyto20.setOpt("allele-peaks-adjustment-covariates","someothername.order=fragment-length,probe-gc.bin-type=discrete,equally-spaced.bin-count=0.coarse-allele-peak-adjustment=off,off");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto20), Except);

    // Expect an exception - allele peak order and adjustment window element count mismatch
    CNCytoEngine cnCyto21;
    cnCyto21.setOpt("allele-peaks-adjustment-covariates","someothername.order=fragment-length,probe-gc.bin-type=discrete,equally-spaced.bin-count=17,23.coarse-allele-peak-adjustment=off,on"
      ".coarse-allele-peak-adjustment-step=-1,2.coarse-allele-peak-adjustment-window=0,100,4.coarse-allele-peak-adjustment-point-count=-1,1.coarse-allele-peak-adjustment-bandwidth=0,4"
      ".coarse-allele-peak-adjustment-cutoff=-1,1.coarse-allele-peak-adjustment-clean-threshold=0,1.coarse-allele-peak-adjustment-outlier-trim=-1,1"
      ".master-peaks-point-count=1.master-peaks-bandwidth=4.covariate-bin-peaks-point-count=1.covariate-bin-peaks-bandwidth=4.kernel-function-selection=gaussian");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto21), Except);

    // Expect an exception - negative non-discrete allele peak bin-count
    CNCytoEngine cnCyto22;
    cnCyto22.setOpt("allele-peaks-adjustment-covariates","someothername.order=fragment-length,probe-gc.bin-type=discrete,equally-spaced.bin-count=17,-2.coarse-allele-peak-adjustment=off,off");
    NEGATIVE_TEST(CovariateParams::checkParams(cnCyto22),Except);

    // No exception - negative discrete allele peak bin-count
    CNCytoEngine cnCyto23;
    cnCyto23.setOpt("allele-peaks-adjustment-covariates","someothername.order=fragment-length,probe-gc.bin-type=equally-spaced,discrete.bin-count=17,-2.coarse-allele-peak-adjustment=off,off");
    CovariateParams::checkParams(cnCyto23);
}

void CovariateParamsTest::isReservedCovariateNameTest()
{
    CPPUNIT_ASSERT(CovariateParams::isReservedCovariateName("fragment-adapter-type"));
    CPPUNIT_ASSERT(CovariateParams::isReservedCovariateName("fragment-length"));
    CPPUNIT_ASSERT(CovariateParams::isReservedCovariateName("fragment-gc"));
    CPPUNIT_ASSERT(CovariateParams::isReservedCovariateName("probe-gc"));
    CPPUNIT_ASSERT(CovariateParams::isReservedCovariateName("local-gc"));
    CPPUNIT_ASSERT(CovariateParams::isReservedCovariateName("median-signal"));
    CPPUNIT_ASSERT(CovariateParams::isReservedCovariateName("marker-class"));
    CPPUNIT_ASSERT(CovariateParams::isReservedCovariateName("nonsense") == false);
}
