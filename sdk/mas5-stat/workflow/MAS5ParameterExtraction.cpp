////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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
#include "calvin_files/parameter/src/ParameterFileData.h"
#include "calvin_files/parsers/src/ParameterFileReader.h"
#include "calvin_files/utils/src/StringUtils.h"
//
#include <cstdlib>
#include <cstring>
#include <string>
//

using namespace affymetrix_calvin_parameter;
using namespace affymetrix_calvin_utilities;
using namespace affymetrix_calvin_io;
using namespace std;

/*
 * Convert a wide string to a float.
 */
static float toFloat(const wstring &str)
{
    string mbs = StringUtils::ConvertWCSToMBS(str);
    double d = strtod(mbs.c_str(),NULL);
    float rv = d;
    return rv;
}

/*
 * Convert a wide string to a float.
 */
static int toInt(const wstring &str)
{
    return (int) toFloat(str);
}

/*
 * Extract the parameters from the parameter file.
 */
bool MAS5ParameterExtraction::ExtractParameters
(
    const char *fileName, 
    CExpStatAlgSettings &parameters
)
{
    // Read the file
    ParameterFileData paramData;
    ParameterFileReader reader;
    if (reader.Read(fileName, paramData) == false)
        return false;

    // Map the parameters from the file to the MAS5 parameter object.
    for (ParameterTypeList::iterator it=paramData.Parameters().begin(); it!=paramData.Parameters().end(); it++)
    {
        if (it->name == L"Alpha1")
            parameters.Alpha1 = toFloat(it->currentValue);

        else if (it->name == L"Alpha2")
            parameters.Alpha2 = toFloat(it->currentValue);

        else if (it->name == L"ScaleFactor")
        {
            parameters.BaseScaleFactor = toFloat(it->currentValue);
            parameters.ScaleFactor = toFloat(it->currentValue);
        }

        else if (it->name == L"BHCoef")
            parameters.BHCoef = toFloat(it->currentValue);

        else if (it->name == L"BiasCorrect")
            parameters.BiasCorrect = toFloat(it->currentValue);

        else if (it->name == L"BLCoef")
            parameters.BLCoef = toFloat(it->currentValue);

        else if (it->name == L"CMultiplier")
            parameters.CMultiplier = toFloat(it->currentValue);

        else if (it->name == L"ContrastTau")
            parameters.ContrastTau = toFloat(it->currentValue);

        else if (it->name == L"Delta")
            parameters.Delta = toFloat(it->currentValue);

        else if (it->name == L"Epsilon")
            parameters.Epsilon = toFloat(it->currentValue);

        else if (it->name == L"EpsilonAvgLogInten")
            parameters.EpsilonAvgLogInten = toFloat(it->currentValue);

        else if (it->name == L"EpsilonAvgLogRatio")
            parameters.EpsilonAvgLogRatio = toFloat(it->currentValue);

        else if (it->name == L"EpsilonGammas")
            parameters.EpsilonGammas = toFloat(it->currentValue);

        else if (it->name == L"EpsilonSB")
            parameters.EpsilonSB = toFloat(it->currentValue);

        else if (it->name == L"Gamma1H")
            parameters.Gamma1H = toFloat(it->currentValue);

        else if (it->name == L"Gamma1L")
            parameters.Gamma1L = toFloat(it->currentValue);

        else if (it->name == L"Gamma2H")
            parameters.Gamma2H = toFloat(it->currentValue);

        else if (it->name == L"Gamma2L")
            parameters.Gamma2L = toFloat(it->currentValue);

        else if (it->name == L"IntensityHighPercent")
            parameters.IntensityHighPercent = toFloat(it->currentValue);

        else if (it->name == L"IntensityLowPercent")
            parameters.IntensityLowPercent = toFloat(it->currentValue);

        else if (it->name == L"NoiseFrac")
            parameters.NoiseFrac = toFloat(it->currentValue);

        else if (it->name == L"NormFactor")
            parameters.NormFactor = toFloat(it->currentValue);

        //else if (it->name == L"")
        //    parameters.NormGenes = ;

        else if (it->name == L"NormMaskFile")
            parameters.NormMaskFile = StringUtils::ConvertWCSToMBS(it->currentValue);

        else if (it->name == L"NFMethod")
        {
            if (it->currentValue == L"DEFINED_NORMALIZATION_FACTOR")
                parameters.NormMethod = CExpStatAlgSettings::DEFINED_NORMALIZATION_FACTOR;

            else if (it->currentValue == L"NORM_TO_ALL_PROBE_SETS")
                parameters.NormMethod = CExpStatAlgSettings::NORM_TO_ALL_PROBE_SETS;

            else if (it->currentValue == L"NORM_TO_SELECTED_PROBE_SETS")
                parameters.NormMethod = CExpStatAlgSettings::NORM_TO_SELECTED_PROBE_SETS;

            else
                parameters.NormMethod = CExpStatAlgSettings::DEFINED_NORMALIZATION_FACTOR;
        }

        else if (it->name == L"NumberBGCells")
            parameters.NumberBGCells = toFloat(it->currentValue);

        else if (it->name == L"NumberHorZones")
            parameters.NumberHorZones = toInt(it->currentValue);

        else if (it->name == L"NumberVertZones")
            parameters.NumberVertZones = toInt(it->currentValue);

        else if (it->name == L"Perturbation")
            parameters.Perturbation = toFloat(it->currentValue);

        else if (it->name == L"ProbeMaskFile")
            parameters.ProbeMaskFile = StringUtils::ConvertWCSToMBS(it->currentValue);

        else if (it->name == L"RelConfInterval")
            parameters.RelConfInterval = toFloat(it->currentValue);

        else if (it->name == L"SaturatedIntensity")
            parameters.SaturatedIntensity = toFloat(it->currentValue);

        else if (it->name == L"HPSaturatedIntensity")
            parameters.HPSaturatedIntensity = toFloat(it->currentValue);

        //else if (it->name == L"")
        //    parameters.ScaleGenes = ;

        else if (it->name == L"ScaleMaskFile")
            parameters.ScaleMaskFile = StringUtils::ConvertWCSToMBS(it->currentValue);

        else if (it->name == L"ScaleTau")
            parameters.ScaleTau = toFloat(it->currentValue);

        else if (it->name == L"SFMethod")
        {
            if (it->currentValue == L"DEFINED_SCALING_FACTOR")
                parameters.SFMethod = CExpStatAlgSettings::DEFINED_SCALING_FACTOR;

            else if (it->currentValue == L"SCALE_TO_ALL_PROBE_SETS")
                parameters.SFMethod = CExpStatAlgSettings::SCALE_TO_ALL_PROBE_SETS;

            else if (it->currentValue == L"SCALE_TO_SELECTED_PROBE_SETS")
                parameters.SFMethod = CExpStatAlgSettings::SCALE_TO_SELECTED_PROBE_SETS;

            else
                parameters.SFMethod = CExpStatAlgSettings::DEFINED_SCALING_FACTOR;
        }

        else if (it->name == L"SmoothFactorBG")
            parameters.SmoothFactorBG = toFloat(it->currentValue);

        else if (it->name == L"STP")
            parameters.STP = toFloat(it->currentValue);

        else if (it->name == L"Tau")
            parameters.Tau = toFloat(it->currentValue);

        else if (it->name == L"TGT")
            parameters.TGT = toFloat(it->currentValue);

        else if (it->name == L"TuningConstantCAvgLogInten")
            parameters.TuningConstantCAvgLogInten = toFloat(it->currentValue);

        else if (it->name == L"TuningConstantCAvgLogRatio")
            parameters.TuningConstantCAvgLogRatio = toFloat(it->currentValue);

        else if (it->name == L"TuningConstantCGammas")
            parameters.TuningConstantCGammas = toFloat(it->currentValue);

        else if (it->name == L"TuningConstantCSB")
            parameters.TuningConstantCSB = toFloat(it->currentValue);
    }

    return true;
}

