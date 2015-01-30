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

#include "mas5-stat/src/ExpStatAlgSettings.h"
//
#include <cmath>
//

/////////////////////////////////////////////////////////////////

static const int DefaultNumberHorizonalZones = 4;
static const int DefaultNumberVerticalZones = 4;
static const float DefaultNumBGCells = 2.0f;
static const float DefaultNormFactor = 1.0f;
static const float DefaultScaleFactor = 1.0f;
static const float DefaultTGT = 500;
static const float DefaultAlpha1 = 0.05f;
static const float DefaultAlpha2 = 0.065f;
static const float DefaultTau = 0.015f;
static const float DefSaturatedIntensity = 65000.0f;
static const float DefHPSaturatedIntensity = 48000.0f;
static const float DefEpsilon = 0.5f;
static const float DefTuningConstantCSB = 5.0f;
static const float DefTuningConstantCAvgLogInten = 5.0f;
static const float DefTuningConstantCAvgLogRatio = 5.0f;
static const float DefTuningConstantCGammas = 5.0f;
static const float DefEpsilonSB = 0.0001f;
static const float DefEpsilonAvgLogInten = 0.0001f;
static const float DefEpsilonAvgLogRatio = 0.0001f;
static const float DefEpsilonGammas = 0.0001f;
static const float DefContrastTau = 0.03f;
static const float DefSmoothFactorBG = 100.0f;
static const float DefNoiseFrac = 0.5f;
static const float DefScaleTau = 10.0f;
static const float log2DELTA = -20;
static const float DefDelta = (float) pow(2.0, (double) log2DELTA);
static const char  DefAlgorithmName[] = "ExpressionStat";
static const char  DefVersionNumber[] = "5.0";
static const float DefIntensityLowPercent = 2.0f;
static const float DefIntensityHighPercent = 2.0f;
static const float DefaultGamma1H = 0.004500f;
static const float DefaultGamma1L = 0.004500f;
static const float DefaultGamma2H = 0.006000f;
static const float DefaultGamma2L = 0.006000f;
static const float DefaultPerturbation = 1.10f;
static const float DefCMultiplier = 0.2f;
static const float DefBHCoef = 7.0f;
static const float DefBLCoef = 0.8f;
static const float DefBiasCorrect = 0.0f;
static const float DefRelConfInterval = 0.975f;
static const float DefSTP = 3.0f;

/////////////////////////////////////////////////////////////////

CExpStatAlgSettings::CExpStatAlgSettings()
{
	SetDefaults();
}

/////////////////////////////////////////////////////////////////

void CExpStatAlgSettings::SetDefaults()
{
	// Initialize parameters
	ScaleGenes.erase(ScaleGenes.begin(), ScaleGenes.end());
	NormGenes.erase(NormGenes.begin(), NormGenes.end());
	NormMethod = DEFINED_NORMALIZATION_FACTOR;
	SFMethod = DEFINED_SCALING_FACTOR;
	NormFactor = DefaultNormFactor;
	ProbeMaskFile = "";
	ScaleMaskFile = "";
	NormMaskFile = "";
	TGT = DefaultTGT;
	ScaleFactor = DefaultScaleFactor;
	NumberBGCells = DefaultNumBGCells;
	NumberHorZones = DefaultNumberHorizonalZones;
	NumberVertZones = DefaultNumberVerticalZones;
	Alpha1 = DefaultAlpha1;
	Alpha2 = DefaultAlpha2;
	Tau = DefaultTau;
	Gamma1H = DefaultGamma1H;
	Gamma1L = DefaultGamma1L;
	Gamma2H = DefaultGamma2H;
	Gamma2L = DefaultGamma2L;
	Perturbation = DefaultPerturbation;
	SaturatedIntensity = DefSaturatedIntensity;
	HPSaturatedIntensity = DefHPSaturatedIntensity;
	Epsilon = DefEpsilon;
	TuningConstantCSB = DefTuningConstantCSB;
	TuningConstantCAvgLogInten = DefTuningConstantCAvgLogInten;
	TuningConstantCAvgLogRatio = DefTuningConstantCAvgLogRatio;
	TuningConstantCGammas = DefTuningConstantCGammas;
	EpsilonSB = DefEpsilonSB;
	EpsilonAvgLogInten = DefEpsilonAvgLogInten;
	EpsilonAvgLogRatio = DefEpsilonAvgLogRatio;
	EpsilonGammas = DefEpsilonGammas;
	ContrastTau = DefContrastTau;
	SmoothFactorBG = DefSmoothFactorBG;
	NoiseFrac = DefNoiseFrac;
	ScaleTau = DefScaleTau;
	Delta = DefDelta;
	CMultiplier = DefCMultiplier;
	BHCoef = DefBHCoef;
	BLCoef = DefBLCoef;
	BiasCorrect = DefBiasCorrect;
	RelConfInterval = DefRelConfInterval;
	STP = DefSTP;
	IntensityLowPercent = DefIntensityLowPercent;
	IntensityHighPercent = DefIntensityHighPercent;
}

/////////////////////////////////////////////////////////////////
