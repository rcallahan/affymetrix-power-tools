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


#include "mas5-stat/workflow/MAS5ParameterExtraction.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class MAS5ParameterExtractionTest : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( MAS5ParameterExtractionTest );

	CPPUNIT_TEST( testExtractParameters );

	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	void tearDown();
    void testExtractParameters();
};

CPPUNIT_TEST_SUITE_REGISTRATION( MAS5ParameterExtractionTest );

void MAS5ParameterExtractionTest::setUp()
{
}

void MAS5ParameterExtractionTest::tearDown()
{
}

void MAS5ParameterExtractionTest::testExtractParameters()
{
    CExpStatAlgSettings p;

    CPPUNIT_ASSERT(MAS5ParameterExtraction::ExtractParameters("nopath", p) == false);
    CPPUNIT_ASSERT(MAS5ParameterExtraction::ExtractParameters("../data/Hu6800.default.mas5_parameters", p) == true);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Alpha1, 0.04, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Alpha2, 0.06, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.BaseScaleFactor, 1, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.BHCoef, 7, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.BiasCorrect, 0, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.BLCoef, 0.8, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.CMultiplier, .2, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.ContrastTau, .03, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Delta, 0.0000009536743, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Epsilon, 0.5, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.EpsilonAvgLogInten, 0.0001, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.EpsilonAvgLogRatio, 0.0001, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.EpsilonGammas, 0.0001, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.EpsilonSB, 0.0001, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Gamma1H, 0.0025, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Gamma1L, 0.0025, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Gamma2H, 0.003, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Gamma2L, 0.003, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.IntensityHighPercent, 2, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.IntensityLowPercent, 2, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.NoiseFrac, 0.5, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.NormFactor, 1, 0.00001f);
    CPPUNIT_ASSERT(p.NormMaskFile == "normmaskfile");
    CPPUNIT_ASSERT(p.NormMethod == CExpStatAlgSettings::DEFINED_NORMALIZATION_FACTOR);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.NumberBGCells, 2, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.NumberHorZones, 4, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.NumberVertZones, 4, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Perturbation, 1.1, 0.00001f);
    CPPUNIT_ASSERT(p.ProbeMaskFile == "probemaskfile");
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.RelConfInterval, 0.975, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.SaturatedIntensity, 65000, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.HPSaturatedIntensity, 48000, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.ScaleFactor, 1, 0.00001f);
    CPPUNIT_ASSERT(p.ScaleMaskFile == "scalemaskfile");
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.ScaleTau, 10, 0.00001f);
    CPPUNIT_ASSERT(p.SFMethod == CExpStatAlgSettings::SCALE_TO_ALL_PROBE_SETS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.SmoothFactorBG, 100, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.STP, 3, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.Tau, 0.015, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.TGT, 500, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.TuningConstantCAvgLogInten, 5, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.TuningConstantCAvgLogRatio, 5, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.TuningConstantCGammas, 5, 0.00001f);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p.TuningConstantCSB, 5, 0.00001f);
}
