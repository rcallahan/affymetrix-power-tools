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

#ifndef AlgSettings_HEADER
#define AlgSettings_HEADER

/////////////////////////////////////////////////////////////////

#include <cstring>
#include <list>
#include <string>
//
using namespace std;

/////////////////////////////////////////////////////////////////

class CExpStatAlgSettings {
public:
	enum NormalizationOptionsEnum
	{
		NORM_TO_ALL_PROBE_SETS,
		NORM_TO_SELECTED_PROBE_SETS,
		DEFINED_NORMALIZATION_FACTOR
	};

	enum ScalingOptionsEnum
	{
		SCALE_TO_ALL_PROBE_SETS,
		SCALE_TO_SELECTED_PROBE_SETS,
		DEFINED_SCALING_FACTOR
	};

	NormalizationOptionsEnum NormMethod;
	ScalingOptionsEnum SFMethod;

	string ProbeMaskFile;
	string ScaleMaskFile;
	string NormMaskFile;

	list<string> ScaleGenes;
	list<string> NormGenes;

	int NumberHorZones;
	int NumberVertZones;
	float NumberBGCells;
	float NormFactor;
	float ScaleFactor;
	float TGT;
	float IntensityLowPercent;
	float IntensityHighPercent;
 
	/* For Absolute Call */

	// User Parameters
	float	Alpha1;
	float	Alpha2;
	float	Tau;

	// Other algorithm parameters
	float	HPSaturatedIntensity;
	float	SaturatedIntensity;
	float	Epsilon;
	float	ContrastTau;
	float	TuningConstantCSB;
	float	TuningConstantCAvgLogInten;
	float	TuningConstantCAvgLogRatio;
	float	TuningConstantCGammas;
	float	EpsilonSB;
	float	EpsilonAvgLogInten;
	float	EpsilonAvgLogRatio;
	float	EpsilonGammas;
	float	SmoothFactorBG;
	float	NoiseFrac;
	float	ScaleTau;
	float	Delta;

	/* For Comparative Call */

	// User Parameters
	float	Gamma1H;
	float	Gamma1L;
	float	Gamma2H;
	float	Gamma2L;
	float	Perturbation;

	// Other algorithm parameters
	float	CMultiplier;
	float	BHCoef;
	float	BLCoef;
	float	BiasCorrect;
	float	RelConfInterval;
	float	STP;

	float	BaseScaleFactor;

	void SetDefaults();
	CExpStatAlgSettings();
};

/////////////////////////////////////////////////////////////////

#endif
