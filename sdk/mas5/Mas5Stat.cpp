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

/// @file   Mas5Stat.cpp
/// @brief  Chipstream compatible implementation of the MAS 5.0
///         statistical algorithm.

#include "mas5/Mas5Stat.h"
//
#include "portability/affy-base-types.h"
#include "util/Convert.h"
//
#include <cassert>
#include <cstring>
#include <string.h>
//

using namespace std;

/**
 * @brief Constructor.
 *
 * @param param Map of key/value pairs to set user-defined parameters.
 *
 * Errors: The conversion utilities call errAbort() if any of the user
 * entered values is of the wrong type.
*/
mas5Stat::mas5Stat (map<string,string>& param)
  : log_2 (log (2.0))
{
  // Zero chips with background count.
  m_ChipsWithBackground = 0;
  // Convert, set user specified parameters.
  convertParams (param);
}

/**
 *  @brief Destructor.
 *
*/
mas5Stat::~mas5Stat()
{
  Clear();
}

/**
 * @brief Convert user defined parameters, set the corresponding
 * class parameters.
 *
 * @param param Map of key/value pairs to set user-defined parameters.
 *
 * Default values are set by the SetDefaults() call in the
 * CExpStatAlgSettings constructor.
 *
 * Errors: The conversion utilities will call errAbort() if any of the user
 * entered value strings cannot be converted into the desired type.
*/
void mas5Stat::convertParams (const map<string,string>& param)
{
  map<string,string>::const_iterator iter;
  const map<string,string>::const_iterator paramEnd = param.end();

  iter = param.find ("Alpha1");
  if (iter != paramEnd)
    m_Params.Alpha1 = Convert::toFloat (iter->second.c_str());

  iter = param.find ("Alpha2");
  if (iter != paramEnd)
    m_Params.Alpha2 = Convert::toFloat (iter->second.c_str());

  iter = param.find ("Tau");
  if (iter != paramEnd)
    m_Params.Tau = Convert::toFloat (iter->second.c_str());

  iter = param.find ("TGT");
  if (iter != paramEnd)
    m_Params.TGT = Convert::toFloat (iter->second.c_str());

  iter = param.find ("Gamma1H");
  if (iter != paramEnd)
    m_Params.Gamma1H = Convert::toFloat (iter->second.c_str());

  iter = param.find ("Gamma1L");
  if (iter != paramEnd)
    m_Params.Gamma1L = Convert::toFloat (iter->second.c_str());

  iter = param.find ("Gamma2H");
  if (iter != paramEnd)
    m_Params.Gamma2H = Convert::toFloat (iter->second.c_str());

  iter = param.find ("Gamma2L");
  if (iter != paramEnd)
    m_Params.Gamma2L = Convert::toFloat (iter->second.c_str());

  iter = param.find ("Perturbation");
  if (iter != paramEnd)
    m_Params.Perturbation = Convert::toFloat (iter->second.c_str());

  iter = param.find ("SFMethod");
  if (iter != paramEnd)
    m_Params.SFMethod
      = (CExpStatAlgSettings::ScalingOptionsEnum) Convert::toInt (iter->second.c_str());

  iter = param.find ("NormMethod");
  if (iter != paramEnd)
    m_Params.NormMethod
      = (CExpStatAlgSettings::NormalizationOptionsEnum) Convert::toInt (iter->second.c_str());

  iter = param.find ("NormFactor");
  if (iter != paramEnd)
    m_Params.NormFactor = Convert::toFloat (iter->second.c_str());

  iter = param.find ("ScaleFactor");
  if (iter != paramEnd)
    m_Params.ScaleFactor = Convert::toFloat (iter->second.c_str());

  m_NumberZones = m_Params.NumberVertZones * m_Params.NumberHorZones;
}

/**
 * @brief Set layout (cdf or pgf/clf combination) specific parameters.
 *
 * @param cols Number of feature columns in the array.
 * @param rows Number of feature rows in the array.
*/
void mas5Stat::setLayoutParams (int cols, int rows, const map<probeid_t, bool>& maskedProbes,
                                const int probeSetCount, ChipLayout &layout)

{
  m_Cols = cols;
  m_Rows = rows;
  m_MaskedProbes = maskedProbes;
  m_MaskedProbesEnd = m_MaskedProbes.end();
  m_ProbeSetCount = probeSetCount;
  m_ProbesetProbes = layout.getProbesetProbes();

  // Determine the number of cells per zone in the vertical direction.
  int CellsRemaining = cols % m_Params.NumberVertZones;
  if (CellsRemaining != 0)
    m_ZoneX = (cols + m_Params.NumberVertZones - CellsRemaining) / m_Params.NumberVertZones;
  else
    m_ZoneX = cols / m_Params.NumberVertZones;

  // Determine the number of cells per zone in the horizontal direction.
  CellsRemaining = rows % m_Params.NumberHorZones;
  if (CellsRemaining != 0)
    m_ZoneY = (rows + m_Params.NumberHorZones - CellsRemaining) / m_Params.NumberHorZones;
  else
    m_ZoneY = rows / m_Params.NumberHorZones;

  // Ensure that there are a match and mismatch cell in the same zone.
  m_ZoneY += m_ZoneY % 2; // EXPRESSION_ATOMS_PER_CELL;

  m_ZoneXF = (float) m_ZoneX;
  m_ZoneYF = (float) m_ZoneY;
}

/**
 * @brief Compute background intensity, noise for each zone.
 *
 * @param data Chip data as vector<float>.
 * @param chipsWithBackgroundCount Count of chips for which background
 *        has been calculated.
*/
void mas5Stat::computeBackground (const vector<float>& data, const int chipsWithBackgroundCount)
{
  // Force the caller to keep track of the number of chips for which background
  // has been calculated.
  assert (chipsWithBackgroundCount == ++m_ChipsWithBackground);
  assert (data.size() == (unsigned int)(m_Cols * m_Rows));

  // Global information that is saved.
  CSAllZonesInfoType* pZonesInfo = new CSAllZonesInfoType;
  m_VecZonesInfo.push_back (pZonesInfo);

  // Allocate space for all atoms intensities and ID's.
  int* NumberCellsPerZone = new int[m_NumberZones];

  // Clear arrays.
  memset (NumberCellsPerZone, 0, sizeof(int)*m_NumberZones);

  // Loop over all units to determine the zone ID's and intensities.
  vector<vector<float> > ZoneCells (m_NumberZones);

  // If there is no probeset information, loop over all probes.
  if (m_ProbeSetCount == 0)
    for (int colIx = 0; colIx < m_Cols; ++colIx)
      for (int rowIx = 0; rowIx < m_Rows; ++rowIx)
      {
        const bool bMasked = isMasked (colIx, rowIx);
        if (! bMasked)
        {
	  const int nZone = DetermineZone (colIx, rowIx);
	  ZoneCells[nZone].resize (ZoneCells [nZone].size() + 1);
	  ZoneCells[nZone] [NumberCellsPerZone [nZone] ] =
	    data [colRowToIndex (colIx, rowIx)];

	  if (nZone >= 0 && nZone < m_NumberZones)
	    NumberCellsPerZone [nZone]++;
        }
      }
  // Else loop over the probesets.
  else
  {
    for(uint32_t posIx = 0; posIx < m_ProbesetProbes.size(); posIx++) {
      if(m_ProbesetProbes[posIx] == true) {
        const int colIx = posIx % m_Cols;
        const int rowIx = posIx / m_Cols;
        const bool bMasked = isMasked (colIx, rowIx);
        if (! bMasked)
	  {
	    const int nZone = DetermineZone (colIx, rowIx);
	    ZoneCells[nZone].resize (ZoneCells [nZone].size() + 1);
	    ZoneCells[nZone] [NumberCellsPerZone [nZone] ] =
	      data [colRowToIndex (colIx, rowIx)];
            
	    if (nZone >= 0 && nZone < m_NumberZones)
	      NumberCellsPerZone [nZone]++;
	  }
      }
    }
  }


  // Allocate zones, set smooth factor and set num zones
  pZonesInfo->pZones = new CSZoneInfo[m_NumberZones];
  pZonesInfo->number_zones = m_NumberZones;
  pZonesInfo->smooth_factor = m_Params.SmoothFactorBG;

  // compute background for each zone
  const float lowBG = 0.0;
  const float highBG = m_Params.NumberBGCells / 100.0f;
  for (int iZone = 0; iZone < m_NumberZones; iZone++)
  {
    // Compute the center coordinates of each zone.
    // (x1,y1) is the upper left corner
    // (x2,y2) is the lower right corner
    float x1 = (float)( ((int) (iZone % m_Params.NumberVertZones)) * m_ZoneX );
    float y1 = (float)( ((int) (iZone / m_Params.NumberVertZones)) * m_ZoneY );
    float x2 = x1 + m_ZoneX;
    float y2 = y1 + m_ZoneY;
    pZonesInfo->pZones[iZone].center.x = (x1 + x2) / 2;
    pZonesInfo->pZones[iZone].center.y = (y1 + y2) / 2;

    int iCell = 0;
    int numCell = NumberCellsPerZone [iZone];
    pZonesInfo->pZones [iZone].numCell = numCell;

    vector<float> zoneI (numCell);
    vector<int> rank (numCell);

    for (int i = 0; i < numCell; i++)
    {
      float inten = ZoneCells[iZone][i];
      zoneI[iCell] = modifyIntensitySlightly (inten);
      iCell++;
    }
    FloatPair fp = trimMeanAndStd (zoneI, lowBG, highBG);
    pZonesInfo->pZones[iZone].background = fp.value1;
    pZonesInfo->pZones[iZone].noise = fp.value2;
  }
  // End of computing background intensity and noise for each zone.
  // Carried zones as the required information.

  delete [] NumberCellsPerZone;
}

/**
 * @brief Compute background adjusted intensity.
 *
 * @param probeIx Probe index from the cel file.
 * @param chipIx Microarray index.
 * @param intensity Original intensity.
 *
 * @return float Adjusted intensity.
 *
 * Following ExpressionAlgorithmImplementation.cpp, we calculate the adjusted
 * intensity whether or not the probe is masked.  That code sets UseAtom[iUnit][iAtom]
 * to false if the probe is masked.
 *
 * Errors: assert chipIx is in range.
*/
float mas5Stat::computeBGAdjustedIntensity (const int probeIx, const int chipIx, const float intensity)
{
  assert (chipIx < m_ChipsWithBackground);
  CSAllZonesInfoType* pZonesInfo = m_VecZonesInfo [chipIx];

  // Compute b(x,y), n(x,y)
  const float smoothF = m_Params.SmoothFactorBG;

  // Set the background (with smoothing adjustment) for matchcell.
  const int cellX = indexToX (probeIx);
  const int cellY = indexToY (probeIx);
  float WeightedSumBg = 0.0f;
  float WeightedSumNoise = 0.0f;
  float WeightedSumDenom = 0.0f;
  float background = 0.0f;
  float noise = 0.0f;

  for (int k = 0; k < m_NumberZones; k++)
  {
    const float weightAtXY = computeWeightAtXY ((float)cellX, (float)cellY, pZonesInfo->pZones[k].center.x,
      pZonesInfo->pZones[k].center.y, smoothF);
    WeightedSumBg    += weightAtXY * pZonesInfo->pZones[k].background;
    WeightedSumNoise += weightAtXY * pZonesInfo->pZones[k].noise;
    WeightedSumDenom += weightAtXY;
  }
  if (WeightedSumDenom != 0.0f)
  {
    background = WeightedSumBg / WeightedSumDenom;
    noise = WeightedSumNoise / WeightedSumDenom;
  }

  float modifiedI = modifyIntensitySlightly (intensity);
  return computeAdjustedIntensity (modifiedI, background, noise);
}

/**
 * @brief Clear data.
 *
 */
void mas5Stat::Clear()
{
  // m_Params are set by the constructor, not cleared here.
  const int vecZoneCount = (int) m_VecZonesInfo.size();
  for (int i = 0; i < vecZoneCount; ++i)
  {
    CSAllZonesInfoType* pZonesInfo = m_VecZonesInfo[i];
    delete [] pZonesInfo->pZones;
    delete pZonesInfo;
  }
  m_VecZonesInfo.clear();
  m_MaskedProbes.clear();
  m_MaskedProbesEnd = m_MaskedProbes.end();

  m_Cols = m_Rows = m_ZoneX = m_ZoneY = m_ChipsWithBackground
   = m_NumberZones = 0;
  m_ZoneXF = m_ZoneYF = 0.0;
}

/**
 * @brief Compute detection p-values.
 *
 * @param pmI Vector of perfect match intensities.
 * @param mmI Vector of mismatch intensities.
 * @param pmIndex Vector of perfect match probe ids.
 * @param mmIndex Vector of mismatch probe ids.
 * @param pValue P-value for detection.
 * @param call Expression call (output).
 * @param pairs Probe pairs (output).
 * @param pairsUsed Probe pairs used (output).
 */
void mas5Stat::computeExpressionStat (const vector<float>& pmI, const vector<float>& mmI,
                                      const vector<probeid_t>& pmId,const vector<probeid_t>& mmId,
                                      double& pValue, unsigned char& call,
                                      int& pairs, int& pairsUsed)
{
  const size_t pmSize = pmI.size();
  assert (pmSize != 0);
  assert (pmSize == mmI.size());
  assert (pmSize == pmId.size());
  assert (pmSize == mmId.size());

  vector<double> discMinusTau (pmSize);
  const double tau1 = m_Params.Tau;

  // Loop over the cells in the unit
  int ctr = 0;
  for (unsigned int i = 0; i < pmSize; ++i)
  {
    const int pmX = pmId[i] % m_Cols;
    const int pmY = pmId[i] / m_Cols;
    const int mmX = mmId[i] % m_Cols;
    const int mmY = mmId[i] / m_Cols;

    if (isMasked (pmX, pmY) == false && isMasked (mmX, mmY) == false)
      // Exclude saturated probe pair.
      if (mmI[i] < m_Params.SaturatedIntensity)
      {
	const double sum = pmI[i] + mmI[i];
	if (sum > 0.0f)
	  discMinusTau[ctr] = ((pmI[i] - mmI[i]) / sum) - tau1;
	else
	  discMinusTau[ctr] =  -tau1;
	ctr ++;
      }
  }

  // Compute the absolute call.
  pairs = pmSize;
  if (ctr <= 0)
  {
    call = EXP_NO_ABS_CALL_TYPE;
    pValue = 0.0f;
    pairsUsed = 0;
    return;
  }

  discMinusTau.resize (ctr);
  ExpResults result = newSignRank (discMinusTau, m_Params.Alpha1, m_Params.Alpha2);

  CSAbsStatExpressionProbeSetResultType unit;
  ComputeAbsoluteCall (&unit, &result);
  call = unit.Detection;

  pairs = pmSize;
  pValue = result.p_value;
  pairsUsed = ctr;
}

/**
 * @brief Compute average measurement.
 *
 * @param pmI Vector of perfect match intensities.
 * @param mmI Vector of mismatch intensities.
 * @param pmIndex Vector of perfect match probe ids.
 * @param mmIndex Vector of mismatch probe ids.
 * @param avgMeasurement Average measurement (output).
 */
void mas5Stat::computeMeasurement (const vector<float>& pmI, const vector<float>& mmI,
                                   const vector<probeid_t>& pmId, const vector<probeid_t>& mmId,
                                   double& avgMeasurement)
{
  const size_t pmSize = pmI.size();
  assert (pmSize != 0);
  assert (pmSize == mmI.size());
  assert (pmSize == pmId.size());
  assert (pmSize == mmId.size());

  // Compute Contrast Value.
  vector<float> CT (pmSize);
  ComputeContrastValue (pmI, mmI, CT);

  // Compute Probe Value.
  vector<float> PV (pmSize);
  ComputeProbeValue (pmI, CT, PV);

  // Compute Average Log Intensity for probe set and its confidence interval.
  const float c = m_Params.TuningConstantCAvgLogInten;
  const float epsilon = m_Params.EpsilonAvgLogInten;

  vector<float> PVused (pmSize);
  GetUsedSet (PV, PVused, pmId, mmId);

  avgMeasurement = OneStepBiweightAlgorithm (PVused, c, epsilon);
  // Multiply by user-configurable scale factor (default 1.0).
  avgMeasurement = antiLog (avgMeasurement) * m_Params.ScaleFactor;
}

/**
 * @brief Compute contrast value.
 *
 * @param pmI Perfect match intensities.
 * @param mmI Mismatch intensities.
 * @param CT Contrast value (output).
 */
void mas5Stat::ComputeContrastValue (const vector<float>& pmI, const vector<float>& mmI,
  vector<float>& CT)
{
  float SB;
  ComputeTypicalDifference (pmI, mmI, SB);

  const float ContrastTau = m_Params.ContrastTau;
  const int nProbePair = (int) pmI.size();
  for (int j = 0; j < nProbePair; j++)
  {
    if (mmI[j] < pmI[j])
      CT[j] = mmI[j];
    else if ((mmI[j] >= pmI[j]) &&
         (SB > ContrastTau))
      CT[j] = pmI[j] / (float)antiLog (SB);
    else if ((mmI[j] >= pmI[j]) &&
         (SB <= ContrastTau))
      CT[j] = pmI[j] / (float)antiLog (ContrastTau / (1.0 + (ContrastTau - SB) / m_Params.ScaleTau) );
  }
}

/**
 * @brief Compute typical difference.
 *
 * @param pmI Perfect match intensities.
 * @param mmI Mismatch intensities.
 * @param SB Step biweight (output).
 */
void mas5Stat::ComputeTypicalDifference (const vector<float>& pmI, const vector<float>& mmI, float& SB)
{
  const float c = m_Params.TuningConstantCSB;
  const float epsilon = m_Params.EpsilonSB;

  const int nProbePair = (int) pmI.size(); // Number of Probe Pairs in Probe Set
  vector<float> logPM_minus_logMM (nProbePair);
  for (int j = 0; j < nProbePair; j++)
    logPM_minus_logMM[j] = logtwo (pmI[j]) - logtwo (mmI[j]);

  SB = OneStepBiweightAlgorithm (logPM_minus_logMM, c, epsilon);
}

/**
 * @brief Compute one-step biweight estimate.
 *
 * @param x Vector of log differences.
 * @param c Tuning constant.
 * @param epsilon Epsilon.
 * @return Estimate.
 */
float mas5Stat::OneStepBiweightAlgorithm (const vector<float> & x, const float& c, const float& epsilon)
{
  const int n = (int) x.size();
  if (n == 0)
    return 0.0f;

  const float medianValue = median (x);
  const float MAD = medianAbsoluteDeviation (x) * c + epsilon;
  float value = 0.0f;
  float weightedSumNumer = 0.0f;
  float weightedSumDenom = 0.0f;

  for (int i = 0; i < n; i++)
  {
    const float diff = x[i] - medianValue;
    const float u = diff / MAD;
    const float uSquare = u * u;
    const float oneMinusUSqaure = 1.0f - uSquare;
    if (fabs (u) < 1.0f)
    {
      weightedSumNumer += diff * oneMinusUSqaure * oneMinusUSqaure;
      weightedSumDenom += oneMinusUSqaure * oneMinusUSqaure;
    }
  }

  if (weightedSumDenom != 0.0f)
    value = medianValue + weightedSumNumer / weightedSumDenom;

  return value;
}

/**
 * @brief Compute probe value.
 *
 * @param pmI Perfect match intensities.
 * @param CT Contrast value.
 * @param PV Probe value (output).
 */
void mas5Stat::ComputeProbeValue (const vector<float>& pmI, const vector<float>& CT, vector<float>& PV)
{
  const float delta = m_Params.Delta;
  const float correction = 1.0f + m_Params.BiasCorrect;
  const int nProbePair = (int) pmI.size();
  for (int j = 0; j < nProbePair; j++)
  {
    const float v = pmI[j] - CT[j];
    if (v < delta)
      PV[j] = correction * (float)logtwo (delta);
    else
      PV[j] = correction * (float)logtwo (v);
  }
}

/**
 * @brief Compute used set
 *
 * @param InputSet Input set.
 * @param pmId Perfect match probe ids.
 * @param mmId Mismatch probe ids.
 * @return vector<float> UsedSet Used set (output).
 */
void mas5Stat::GetUsedSet (const vector<float>& InputSet, vector<float>& UsedSet,
                           const vector<probeid_t>& pmId, const vector<probeid_t>& mmId)
{
  const int size = (int) InputSet.size();
  int used = 0;
  for (int j = 0; j < size; j++)
  {
    const int pmX = pmId[j] % m_Cols;
    const int pmY = pmId[j] / m_Cols;
    const int mmX = mmId[j] % m_Cols;
    const int mmY = mmId[j] / m_Cols;

    if (isMasked (pmX, pmY) == false && isMasked (mmX, mmY) == false)
    {
      UsedSet[used] = InputSet[j];
      used++;
    }
  }
  UsedSet.resize (used);
}
