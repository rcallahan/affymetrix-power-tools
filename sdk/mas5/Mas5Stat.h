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
 * @class   mas5Stat
 * @brief  ChipStream and AnalysisStream compatible implementation
 * of the MAS 5.0 statistical algorithms.
 *
 * An implementation of the methods described in Affymetrix
 * Technical Notes "Fine Tuning Your Data Analysis" (2001),
 * "New Statistical Algorithms for Monitoring Gene Expression
 * on GeneChip® Probe Arrays" (2001), "Statistical Algorithms
 * Reference Guide." (2001), and in "Robust estimators for
 * expression analysis", E. Hubbell, W.-m. Liu and R. Mei,
 * Bioinformatics 2002 18:1585-1592, and "Analysis of high
 * density expression microarrays with signed-rank call algorithms",
 * W.-m. Liu, R. Mei, X. Di, T.B. Ryder, E. Hubbell, S. Dee,
 * T.A. Webster, C.A. Harrington, M.-h. Ho, J. Baid and
 * S.P. Smeekens, Bioinformatics 2002 18:1593-1599.
 *
 * This is based on the implementation in class
 * CExpressionAlgorithmImplementation.
 *
 * This class provides ChipStream compatible background
 * correction and AnalysisStream compatible detection
 * of expressed genes and signal estimation.
 *
 * Detection calls use discrimination scores based on
 * Wilcoxon's signed-rank test.  P-values are calculated,
 * giving a confidence level for the pertinent hypothesis.
 * Probe sets for which the p-value for expression is below
 * the Alpha1 threshold are called present, those with
 * p-value above the Alpha2 threshold are called absent,
 * and those with p-value between Alpha1 and Alpha2 are
 * called marginal.
 *
 * In background correction, the microarray is divided
 * into N x N zones (default N = 4); the background in
 * a zone is the average of the lowest 2% of probe set cell
 * intensities within that zone, using a smoothing
 * method to avoid the discontinuity of background
 * on the boundary of zones.
 */

#ifndef _MAS5STAT_H_
#define _MAS5STAT_H_

//
#include "chipstream/ChipLayout.h"
#include "mas5-stat/src/ExpStatAlgSettings.h"
#include "mas5-stat/src/ExpTmpl.h"
#include "mas5-stat/src/mathLib.h"
#include "mas5-stat/src/pTable.h"
//
#include <map>
#include <vector>
//

#define EXP_PRESENT_CALL_TYPE  0
#define EXP_MARGINAL_CALL_TYPE 1
#define EXP_ABSENT_CALL_TYPE   2
#define EXP_NO_ABS_CALL_TYPE   3

////////////////////////////////////////////////////////////////////

typedef struct
{
  float x;
  float y;
} CSCoordinate;

////////////////////////////////////////////////////////////////////

typedef struct
{
  CSCoordinate center;
  int	  numCell;
  float background;
  float noise;
} CSZoneInfo;

////////////////////////////////////////////////////////////////////

typedef struct
{
  int number_zones;
  float smooth_factor;
  CSZoneInfo *pZones;
} CSAllZonesInfoType;

//////////////////////////////////////////////////////////////////////

typedef struct _CSAbsStatExpressionProbeSetResultType
{
  float DetectionPValue;
  float Signal;
  unsigned short NumPairs;
  unsigned short NumUsedPairs;
  unsigned char Detection;
  _CSAbsStatExpressionProbeSetResultType *operator=(_CSAbsStatExpressionProbeSetResultType &src)
  {
    memcpy (this, &src, sizeof(_CSAbsStatExpressionProbeSetResultType));
    return this;
  };
} CSAbsStatExpressionProbeSetResultType;

/////////////////////////////////////////////////////////////////////////////
// Template Functions
/////////////////////////////////////////////////////////////////////////////
template<class T> ExpResults oneSidedSignRank2 (const vector<T> & x, const double alpha)
{
  const int len = (int) x.size();
  // 1. Ignore all zero differences.
  vector <T> newdiff = x;
  int n = 0;
  int i;
  for (i = 0; i < len; ++i)
    if (x[i] != 0.0)
    {
      newdiff[n] = x[i];
      ++n;
    }
  if (n == 0) // No non-zero differences.  Output 0.5 as the one-sided p-value and detection is absent.
    return ExpResults (0.5, 0);
  newdiff.resize (n);

  // 2.  Assign integer ranks to the differences.
  vector<T> ranks (n);
  for (i = 0; i < n; ++i)
    ranks[i] = (float)(i + 1);

  // 3. Convert differences to absolute values and sort in ascending order.
  vector <struct ExpResults> absdiff (n);
  for (i = 0; i < n; ++i)
  {
    absdiff[i].p_value = fabs (newdiff[i]);
    absdiff[i].call = i;
  }
  sort (absdiff.begin(), absdiff.end());

  // 4. If there are ties among absolute differences, all differences in a tie
  //    group are assigned to a rank equal to the average of the integer ranks.
  int nTies = 0;

  // Avoid cross-platform compatibility problems from approximate floating point arithmetic.
  const double tiny = 2e-09;
  for (i = 0; i < n - 1; ++i)
    if (fabs (absdiff[i].p_value - absdiff[i+1].p_value) < tiny)
    {
      nTies ++;
      break;
    }

  int tieGroup = 0;
  double doubleVarMod = 0; // modification of variance due to ties.
  if (nTies)
  {
    i = 0;
    while ( i < n - 1)
    {
      const double initElement = absdiff[i].p_value;
      int tieGroupSize = 1;
      for (int j = i + 1; j < n; ++j)
      {
        if (fabs (absdiff[j].p_value - initElement) < tiny)
	{
          tieGroupSize ++;
          if (j == n - 1)
	  {
            i = j;
            tieGroup ++;
            for (int m = j - tieGroupSize + 1; m <= j; ++m)
              ranks[m] = (2*j - tieGroupSize + 3) / 2.0f;
            doubleVarMod += tieGroupSize *
              ((double)tieGroupSize * tieGroupSize - 1);
            break;
          }
        }
        else
	{
          i = j;
          if (tieGroupSize > 1)
	  {
            tieGroup ++;
            for (int m = j - tieGroupSize; m <= j - 1; ++m)
              ranks[m] = (2 * j - tieGroupSize + 1) / 2.0f;
            doubleVarMod += tieGroupSize *
              ((double)tieGroupSize * tieGroupSize - 1);
          }
          break;
        }
      }
    }
  }
  vector <T> invr (n);
  for (i = 0; i < n; ++i)
    invr[absdiff[i].call] = ranks[i];

  double w = 0;
  for (i = 0; i < n; ++i)
    if (newdiff[i] > 0)
      w += invr[i];

  struct ExpResults ans;
  if (n > 11)
  {
    // Use the asymptotic approximation:
    // S' = [S - n (n+1)/4]/sqrt[n(n+1)(2n+1)/24 - c]
    // where S is the sum of all positive signed ranks
    // and c = sum (b(b*b - 1)/48 is the modification of variance due to ties.
    //
    // p-value = 1 - f (S')
    // where f (S') is the cumulative distribution function of
    // standard normal distribution.
    double dw = w - ((double) n) * (n + 1) / 4.0;
    double denom2 = (((double)n)*(n+1)*(2*n+1) - 0.5*doubleVarMod)/24.0;
    if (denom2 <=0)
      return ExpResults (0, 0);

    double z = dw / sqrt (denom2);
    ans.p_value = 1 - normalCDF(z);
  }
  else
  {
    if (nTies == 0)
    {
      int iCode = 0;
      for (i=0; i < n; ++i)
        if (newdiff[i] > 0)
          iCode += 1 << ((int) invr[i] - 1);
      ans.p_value = fGetPValue (n-1, iCode);
    }
    else
    {
      int twoToN = 1 << n;
      vector<int> mask (n);
      for (i = 0; i < n; ++i)
        mask[i] = 1 << i;
      vector <int> posRanks (twoToN);
      for (i = 0; i < twoToN; ++i)
      {
        double sum = 0;
        for (int j = 0; j < n; ++j)
          if (i & mask[j])
            sum += ranks[j];
        posRanks[i] = (int) sum;
      }
      double tail = 0;
      for (i = 0; i < twoToN; ++i)
      {
        if (posRanks[i] > w)
          tail ++;
        else if (posRanks[i] == w)
          tail += 0.5;
      }
      ans.p_value = tail / (double) twoToN;
    }
  }
  ans.call = (ans.p_value < alpha) ? 1 : 0;
  return ans;
}

//////////////////////////////////////////////////////////////////////

/**
 * @brief Calculate trimmed mean and standard deviation.
 *
 * @param v Vector of data.
 * @param p1 Low threshold (unused).
 * @param p2 High threshold.
 * @return FloatPair Mean, Standard deviation.
 */
template <class T> FloatPair trimMeanAndStd (vector<T> & v, const double p1, const double p2) {
  int total = (int) v.size();
  FloatPair fp (0, 0);
  if (total > 0)
  {
    sort (v.begin(), v.end());
    int n1 = 0;
    int n2 = (int)floor (total * p2);
    double subtotal = n2;
    if (subtotal >= 2)
    {
      double sum = 0;
      for (int i = n1; i < n2; ++i)
        sum += v[i];
      double tMean = sum / subtotal;

      sum = 0;
      for (int i = n1; i < n2; ++i)
        sum += pow (v[i] - tMean, 2);

      fp.value1 = (float) tMean;
      fp.value2 = (float) sqrt (sum / (subtotal - 1));
    }
    else if (subtotal == 1)
      fp.value1 = v[n1];
  }
  return fp;
}

//////////////////////////////////////////////////////////////////////

/**
 * @brief Calculate median.
 *
 * @param x Vector of data.
 * @return Float result.
 */
template<class T> float median(const vector<T> & v) {
  vector <T> u = v;
  int len = (int) u.size();
  if (len < 1)
    return -1;

  sort(u.begin(), u.end());
  int half = len / 2;
  if (len % 2 == 1)
    return (float)u[half];

  else
    return ((float)u[half - 1] + (float)u[half]) / 2.0;

}

//////////////////////////////////////////////////////////////////////

/**
 * @brief Calculate median absolute deviation.
 *
 * @param x Vector of data.
 * @return Float result.
 */
template<class T> float medianAbsoluteDeviation (const vector<T> & x)
{
  const float medianValue = median (x);
  const int size = (int) x.size();
  vector<T> v (size);
  for (int i = 0; i < size;i++)
    v[i] = fabs (x[i] - medianValue);

  return median (v);
}

//////////////////////////////////////////////////////////////////////

/**
 * @brief Calculate sign rank from differences.
 *
 * @param dif Vector of differences.
 * @param alpha1 Threshold for present call.
 * @param alpha2 Threshold for absent call.
 * @return ExpResults p_value, call.
 */
template<class T> ExpResults newSignRank (const vector<T> & dif,
  const double alpha1, const double alpha2)
{
  struct ExpResults newPH;
  struct ExpResults oldPH = oneSidedSignRank2 (dif, alpha1);

  newPH.p_value = oldPH.p_value;
  if (oldPH.call == 1)
    newPH.call = 2;
  else if (oldPH.p_value < alpha2)
    newPH.call = 1;
  else if (oldPH.p_value > 1 - alpha1)
    newPH.call = -2;
  else if (oldPH.p_value > 1 - alpha2)
    newPH.call = -1;
  else
    newPH.call = 0;

  return newPH;
}

/**
 *  mas5Stat
 *  @brief Object for the MAS 5.0 statistical algorithm.
 *
 */
class mas5Stat
{
public:

  /** Constructor.
   * @param param Map of key/value pairs to set user-defined parameters.
   */
  mas5Stat (std::map<std::string,std::string>& param);

  /** Destructor.
   */
  ~mas5Stat();

  /**
   * @brief Set layout (cdf or pgf/clf combination) specific
   * parameters.
   *
   * @param cols Number of feature columns in the array.
   * @param rows Number of feature rows in the array.
   * @param maskedProbes Map describing probes to ignore.
   */
  void setLayoutParams (int cols, int rows, const std::map<probeid_t, bool>& maskedProbes,
                        const int probeSetCount, ChipLayout &layout);

  /**
   * @brief Compute background intensity, noise for each zone.
   *
   * @param data Chip data as vector<float>.
   * @param chipsWithBackgroundCount Count of chips for which background
   *        has been calculated.
   */
  void computeBackground (const std::vector<float>& data, const int chipsWithBackgroundCount);

  /**
   * @brief Compute background adjusted intensity.
   *
   * @param probeIx Probe index from the cel file.
   * @param chipIx Microarray index.
   * @param intensity Original intensity.
   *
   * @return float Adjusted intensity.
   */
  float computeBGAdjustedIntensity (const int probeIx, const int chipIx, const float intensity);

  /**
   * @brief Compute detection p-values.
   *
   * @param pmI Vector of perfect match intensities.
   * @param mmI Vector of mismatch intensities.
   * @param pmIndex Vector of perfect match probe ids.
   * @param mmIndex Vector of mismatch probe ids.
   * @param pValue P-value for detection (output).
   * @param call Expression call (output).
   * @param pairs Probe pairs (output).
   * @param pairsUsed Probe pairs used (output).
   */
  void computeExpressionStat (const std::vector<float>& pmI, const std::vector<float>& mmI,
                              const std::vector<probeid_t>& pmId, const std::vector<probeid_t>& mmId,
                              double& pValue, unsigned char& call, int& pairs, int& pairsUsed);

  /**
   * @brief Compute average measurement.
   *
   * @param pmI Vector of perfect match intensities.
   * @param mmI Vector of mismatch intensities.
   * @param pmIndex Vector of perfect match probe ids.
   * @param mmIndex Vector of mismatch probe ids.
   * @param avgMeasurement Average measurement (output).
   */
  void computeMeasurement (const std::vector<float>& pmI, const std::vector<float>& mmI,
                           const std::vector<probeid_t>& pmId, const std::vector<probeid_t>& mmId,
                           double& avgMeasurement);

  /**
   * @brief Set saturated intensity parameter.
   *
   * This is to be used for testing only.
   *
   * @param saturatedIntensity Saturated intensity.
   */
  void setSaturatedIntensity (const float& saturatedIntensity)
  {
    m_Params.SaturatedIntensity = saturatedIntensity;
  }

  /**
   * @brief Clear data.
   *
   */
  void Clear();

private:

  /**
   * @brief Convert user defined parameters, set the corresponding
   * class parameters.
   *
   * @param param Map of user-defined key/value pairs.
   */
  void convertParams (const std::map<std::string,std::string>& param);

  /**
   * @brief Determine if the probe with the given coordinates is masked.
   *
   * Currently disabled - always returns false.
   *
   * @param col Column.
   * @param row Row.
   * @return True if masked, else false.
   */
  inline bool isMasked (const int col, const int row)
  {
#ifdef SUPPORT_MAS5_PROBE_MASKING
    return m_MaskedProbes.find (colRowToIndex (col,row)) != m_MaskedProbesEnd;
#else
    return false;
#endif
  }

  /**
   * @brief Convert column, row (x, y) to index.
   *
   * Follow CCELFileData::XYToIndex.
   *
   * @param col Column.
   * @param row Row.
   * @return int Index.
   */
  inline int colRowToIndex (const int col, const int row)
  {
    return row * m_Cols + col;
  }

  /**
   * @brief Get X coordinate from index.
   *
   * Follow CCELFileData::IndexToX().
   *
   * @param index Index.
   * @return int X coordinate.
   */
  inline int indexToX (const int index)
  {
    return index % m_Cols;
  }

  /**
   * @brief Get Y coordinate from index.
   *
   * Follow CCELFileData::IndexToY().
   *
   * @param index Index.
   * @return int Y coordinate.
   */
  inline int indexToY (const int index)
  {
    return index / m_Cols;
  }

  /**
   * @brief Find zone index from the cell x, y coordinates.
   *
   * @param cellx Cell x coordinate (column).
   * @param celly Cell y coordinate (row).
   * @return int Zone.
   */
  int DetermineZone (const int cellx, const int celly)
  {
    const float fZx = (float) cellx / m_ZoneXF;
    const float fZy = (float) celly / m_ZoneYF;

    const int Zx = (int) floor (fZx);
    const int Zy = (int) floor (fZy);

    const int zoneID = Zx + Zy * m_Params.NumberVertZones;
    return zoneID;
  }

  /**
   * @brief Set minimum intensity to Epsilon.
   *
   * @param intensity Intensity.
   * @return float Modified intensity.
   */
  inline float modifyIntensitySlightly (const float& intensity)
  {
    return max (intensity, m_Params.Epsilon);
  }

  /**
   * @brief Adjust intensity based on background and noise.
   *
   * @param intensity Intensity.
   * @param background Background.
   * @param noise Noise.
   * @return float Modified intensity.
   */
  inline float computeAdjustedIntensity (const float& intensity, const float& background, const float& noise)
  {
    const float factoredNoise = noise * m_Params.NoiseFrac;
    const float diff = intensity - background;
    return max (max(diff, factoredNoise), 0.5f);
    // AlexC - 1/22/01
    // Code Comments:
    // if too frequent substitution of the noise value, alert might generated.
    // Refer to the page 4 of the Background and Spatial Variation Adjustement spec.,
    // under eq (16), it said that
    // "Production software should probably alert the user if the noise value is being
    //  substituted too frequently (indicating that too much data is below the noise level),
    //  but an appropriate threshold value is not at present known."
  }

  /**
   * @brief Compute weight at x, y coordinates based on the
   * distance from center x and y coordinates.
   *
   * @param x X coordinate.
   * @param y Y coordinate.
   * @param centerX Center X coordinate.
   * @param centerY Center Y coordinate.
   * @param smoothFactor Smoothing factor.
   * @return float Weight.
   */
  inline float computeWeightAtXY (const float& x, const float& y, const float& centerX,
      const float& centerY, const float& smoothFactor)
  {
    return 1.0f / (computeSquaredDistance (x , y, centerX, centerY) + smoothFactor);
  }

  /**
   * @brief Compute squared distance between (x1,y1) and (x2,y2).
   *
   * @param x1 X1 coordinate.
   * @param y1 Y1 coordinate.
   * @param x2 X2 coordinate.
   * @param y2 Y2 coordinate.
   * @return float Distance.
   */
  inline float computeSquaredDistance (const float& x1, const float& y1, const float& x2, const float& y2)
  {
    const float diffx = x1 - x2;
    const float diffy = y1 - y2;
    return diffx * diffx + diffy * diffy;
  }

  /**
   * @brief Compute the absolute call.
   *
   * @param unit Absolute result.
   * @param result Result.
   */
  inline void ComputeAbsoluteCall (CSAbsStatExpressionProbeSetResultType *unit, ExpResults *res)
  {
    // Determine if the gene is present.
    if ( res->call == 2)
      unit->Detection = EXP_PRESENT_CALL_TYPE;

    // Determine if marginal
    else if ( res->call == 1)
      unit->Detection = EXP_MARGINAL_CALL_TYPE;

    // Otherwise the gene is absent.
    else
      unit->Detection = EXP_ABSENT_CALL_TYPE;
  }

  /**
   * @brief Compute contrast value.
   *
   * @param pmI Perfect match intensities.
   * @param mmI Mismatch intensities.
   * @param CT Contrast value (output).
   */
  void ComputeContrastValue (const vector<float>& pmI, const vector<float>& mmI, vector<float>& CT);

  /**
   * @brief Compute probe value.
   *
   * @param pmI Perfect match intensities.
   * @param CT Contrast value.
   * @param PV Probe value (output).
   */
  void ComputeProbeValue (const vector<float>& pmI, const vector<float>& CT, vector<float>& PV);

  /**
   * @brief Compute typical difference.
   *
   * @param pmI Perfect match intensities.
   * @param mmI Mismatch intensities.
   * @param SB Step biweight (output).
   */
  void ComputeTypicalDifference (const vector<float>& pmI, const vector<float>& mmI, float& SB);

  /**
   * @brief Compute logarithm to the base two.
   *
   * @param value Value to compute logarithm of.
   * @return double Logarithm to the base two.
   */
  inline const double logtwo (const double value)
  {
    return (log (value) / log_2);
  }

  /**
   * @brief Compute the inverse of the logarithm
   * (exponent) to the base two.
   *
   * @param value Value to compute exponent of.
   * @return double Exponent to the base two.
   */
  inline const double antiLog (const double value)
  {
    return pow (2.0, value);
  }

  /**
   * @brief Compute one-step biweight estimate.
   *
   * @param x Vector of log differences.
   * @param c Tuning constant.
   * @param epsilon Epsilon.
   * @return float Estimate.
   */
  float OneStepBiweightAlgorithm (const vector<float>& x, const float& c, const float& epsilon);

  /**
   * @brief Compute used set
   *
   * @param InputSet Input set.
   * @param pmId Perfect match probe ids.
   * @param mmId Mismatch probe ids.
   * @return vector<float> UsedSet Used set (output).
   */
  void GetUsedSet (const vector<float>& InputSet, vector<float>& UsedSet,
                   const vector<probeid_t>& pmId, const vector<probeid_t>& mmId);

  /// private data
  /// parameters
  CExpStatAlgSettings m_Params;
  /// number of feature columns
  int m_Cols;
  /// number of feature rows
  int m_Rows;
  /// number of probesets
  int m_ProbeSetCount;
  /// Bitmask where true indicates that probes are present in a probeset
  std::vector<bool> m_ProbesetProbes;
  /// number of columns of features per zone
  int m_ZoneX;
  /// number of rows of features per zone
  int m_ZoneY;
  /// number of columns of features per zone (float)
  float m_ZoneXF;
  /// number of rows of features per zone (float)
  float m_ZoneYF;
  /// count of chips with background calculated
  int m_ChipsWithBackground;
  /// vector of info for all zones
  std::vector<CSAllZonesInfoType*> m_VecZonesInfo;
  /// masked probes
  std::map<probeid_t, bool> m_MaskedProbes;
  /// end of masked probes
  std::map<probeid_t, bool>::const_iterator m_MaskedProbesEnd;
  /// total number of zones
  int m_NumberZones;
  /// logarithm of two
  const double log_2;

};

#endif /* _MAS5STAT_H_ */
