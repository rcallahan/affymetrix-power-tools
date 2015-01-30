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
 * @file   IterPlier.h
 * @author Chuck Sugnet
 * @date   Thu Aug 25 15:09:01 2005
 * 
 * @brief  Class containing methods to iteratively call
 *         plier to improve results.
 */

#ifndef ITERPLIER_H
#define ITERPLIER_H

//
#include "plier/affyplier.h"
//
#include "stats/stats.h"
//
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

/**
 * @class IterPlier
 * @brief Runs plier() interatively using the "best" features in each
 *        iteration, currently the features most correlated with the initial
 *        target signal estimates.
 */
class IterPlier {

public :
    static const int noData = -10000;
  /** 
   * Constructor. Create our own local memory caches.
   * 
   * @param maxFeature - Maximum initial number of features, if more features
   * are needed internal buffers will automatically be enlarged.
   * @param maxChip - Maximum initial number of chips, if more chips are needed
   * internal buffers will automatically be enlarged.
   */
  IterPlier();
  IterPlier(int chip_cnt,int feature_cnt);

  /** 
   * Clean up the local memory caches.
   */
  ~IterPlier();

private:
  void memZap();
  void memFree();
  void memAlloc(int chip_cnt,int feature_cnt);

public:
  void memEnsureSize(int chip_cnt,int feature_cnt);

  /** 
   * Iteratively call plier decreasing the number of features by selecting those
   * that are most correlated with the target response from the previous
   * iteration. 
   * 
   * @param plier - plier object with parameters already set.
   * @param atomCount - Number of atoms. Indicates feature pairs, length of featureRespons,
   *                    and number columns in pm & mm matrices.
   * @param chipCount - Number of experiments. Indicates the hybridizations, length of
   *                    targetResponse and number of rows in pm & mm matrices.
   * @param pm - Matrix of perfect match feature intensities (chipCount rows and atomCount cols);
   * @param mm - Matrix of background feature intensities (chipCount rows and atomCount cols);
   * @param targetResponse - Where signal estimates will be filled in.
   * @param featureResponse - Where feature effects will be filled in.
   * @param residuals - Residuals left after fitting.
   * @param atomUsed - Array of length atomCount indicating features used for final fit.
   * @param iterCount - Number of additional iterations of plier to do.
   * @param iterFeature - Number of features each iteration after initial plier call. e.g iterCount = 2,
   *                      and iterFeature = 22,11
   * @return long plier errorcode.
   */
  long runCorrIterPlier(caffyplier &plier,
                        int atomCount,
                        int chipCount, 
                        double** pm, 
                        double **mm,
                        double *targetResponse, 
                        double *featureResponse,
                        double **residuals,
                        int *atomUsed, 
                        int iterCount,
                        int *iterFeature);

private :
  /* Try not to allocate memory every time called. */
  double **m_pm;  ///< Perfect match matrix to fill in for calling plier iteratively.
  double **m_mm;  ///< Mismatch matrix to fill in for calling plier iteratively.
  double *m_featureResponse; ///< Feature Response array for calling plier iteratively.
  double *m_targetResponse;  ///< Target response array for calling plier iteratively.
  //
  double **residuals;

  //
  int m_cntChip;
  int m_cntFeature;
  //
  int m_maxChip;    ///< Maximum number of chips that m_pm and m_mm can hold.
  int m_maxFeature; ///< Maximum number of features that m_pm and m_mm can hold;

  /** 
   * Run plier choosing only the featureNumToUse best correlated features
   * for the new target signal estimate. 
   * 
   * @param plier - plier object with parameters already set.
   * @param atomCount - Number of atoms. Indicates feature pairs, length of featureRespons,
   *                    and number columns in pm & mm matrices.
   * @param chipCount - Number of experiments. Indicates the hybridizations, length of
   *                    targetResponse and number of rows in pm & mm matrices.
   * @param pm - Matrix of perfect match feature intensities (chipCount rows and atomCount cols);
   * @param mm - Matrix of background feature intensities (chipCount rows and atomCount cols);
   * @param targetResponse - Where signal estimates will be filled in.
   * @param featureResponse - Where feature effects will be filled in.
   * @param residuals - Residuals left after fitting.
   * @param atomUsed - Array of length atomCount indicating features used for final fit.
   * @param featureNumToUse - Number of features to use for new estimate.
   * 
   * @return plier errorcode, zero on success.
   */
  long iterateCorrPlier(caffyplier &plier,
                        int atomCount,
                        int chipCount,
                        double** pm,
                        double **mm, 
                        double *targetResponse,
                        double *featureResponse, 
                        double **residuals,
                        int *atomUsed,
                        int featureToUse);

  /** Comparator struct for sorting. */
  struct less_pair : public binary_function< pair<int,double>, pair<int,double>, bool> {

    bool operator()(const pair<int, double> &x, const pair<int, double> &y) {     
    return x.second < y.second; }
  };

  /** 
   * Copy a column from the matrix formed by pm - mm to the array.
   * 
   * @param pm - perfect match matrix.
   * @param mm - mismatch (or background) matrix.
   * @param colIx - column to copy to array.
   * @param length - number of rows in matrix and length of array.
   * @param array - array to copy column into.
   */  
  void copyColumnToRow(double **pm, double **mm, int colIx, int length, double *array);
};

#endif /* ITERPLIER_H */
