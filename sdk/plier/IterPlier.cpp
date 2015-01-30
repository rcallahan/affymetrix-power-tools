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
 * @file   IterPlier.cpp
 * @author Chuck Sugnet
 * @date   Thu Aug 25 15:15:45 2005
 * 
 * @brief  Class containing methods to iteratively call
 *         plier to improve results.
 */
#include "plier/IterPlier.h"
//
#include "stats/stats.h"
//
#include <cstring>
#include <string>
#include <vector>

/** 
 * Constructor. Create our own local memory caches.
 * 
 * @param maxFeature - Maximum initial number of features, if more features
 * are needed internal buffers will automatically be enlarged.
 * @param maxChip - Maximum initial number of chips, if more chips are needed
 * internal buffers will automatically be enlarged.
 */

IterPlier::IterPlier() {
  memZap();
}
IterPlier::IterPlier(int chip_cnt, int feature_cnt) {
  memZap();
  memAlloc(chip_cnt,feature_cnt);
}
 
//
IterPlier::~IterPlier() {
  memFree();
}

/// @brief     Zero out the vars to a known free state.
void IterPlier::memZap() {
  m_cntChip=0;
  m_cntFeature=0;
  m_maxChip=0;
  m_maxFeature=0;
  m_pm=NULL;
  m_mm=NULL;
  m_featureResponse=NULL;
  m_targetResponse=NULL;
}  
  
/// @brief     Free all the memory allocated.
void
IterPlier::memFree() {
  for(int i = 0; i < m_maxChip; i++) {
    delete [] m_pm[i];
    delete [] m_mm[i];
  }
  delete [] m_pm;
  delete [] m_mm;
  //
  delete [] m_featureResponse;
  delete [] m_targetResponse;
  // memZap so repeated calles to memFree isnt a double delete.
  memZap();
}

/// @brief     Allocate memory of this size
/// @param     chip_cnt     the number of chips
/// @param     feature_cnt  the number of features
void
IterPlier::memAlloc(int chip_cnt, int feature_cnt) {
  // discard all our memory...
  memFree();

  /// ... and reallocate
  m_maxChip=chip_cnt;
  m_maxFeature=feature_cnt;

  m_pm = new double *[m_maxChip];
  m_mm = new double *[m_maxChip];
  for(int i = 0; i < m_maxChip; i++) {
    m_pm[i] = new double [m_maxFeature];
    m_mm[i] = new double [m_maxFeature];
  }
  m_featureResponse = new double[m_maxFeature];
  m_targetResponse = new double[m_maxChip];
}

/// @brief     ensure our allocated memory is at least this big
/// @param     chip_cnt     min number of chips
/// @param     feature_cnt  min number of features
void
IterPlier::memEnsureSize(int chip_cnt, int feature_cnt) {
  if ((m_maxChip<chip_cnt)||(m_maxFeature<feature_cnt)) {
    memAlloc(chip_cnt,feature_cnt);
  }
}

/** 
 * Copy a column from the matrix formed by pm - mm to the array.
 * 
 * @param pm - perfect match matrix.
 * @param mm - mismatch (or background) matrix.
 * @param colIx - column to copy to array.
 * @param length - number of rows in matrix and length of array.
 * @param array - array to copy column into.
 */  
void IterPlier::copyColumnToRow(double **pm, double **mm, int colIx, int length, double *array) {
  int i = 0;
  for(i = 0; i < length; i++) {
    array[i] = pm[i][colIx] - mm[i][colIx];
  }
}


/** 
 * Run plier choosing only the featureNumToUse best correlated features
 * for the new target signal estimate. If correlation cannot be calculated
 * (i.e. due to zero variance) just use first featureNumToUse probes.
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
long IterPlier::iterateCorrPlier(caffyplier &plier, int atomCount, int chipCount, double** pm, double **mm, 
                                 double *targetResponse, double *featureResponse, 
                                 double **residuals, int *atomUsed, int featureNumToUse) {
  int i = 0, j = 0;
  int passCount = 0;
  long errorCode = 0;
  double minVal = 2;
  int corrValidProbes = 0; // Correlation can return nan, if the are all probes are nan, just use orig order.
  vector< pair<int, double> > corrVec(atomCount);

  // Make sure our mem can fit the data.
  memEnsureSize(chipCount,atomCount);

  /* Fill up the featureResponse vector and residual matrix with bad
    numbers to make sure people aren't trying to use invalid data. */
  for(i = 0; i < atomCount; i++) {
    if(featureResponse != NULL)
      featureResponse[i] = IterPlier::noData; 
    if(residuals != NULL) {
      for(j = 0; j < chipCount; j++) 
        residuals[j][i] = IterPlier::noData;
    }
  }
       
  /* For each feature that is going to be used copy into our pm
     matrix and calculate correlation with the current target estimates. */
  for(i = 0; i < atomCount; i++) {
    corrVec[i].first = i;
    if((atomUsed)[i] == 1) {
      copyColumnToRow(pm,mm,i,chipCount, m_targetResponse);
      corrVec[i].second = correlation_coeff(m_targetResponse, 
                                            m_targetResponse+chipCount, 
                                            targetResponse);
      // we use !(corrVec[i].second==corrVec[i].second) rather than isnan() as isnan is
      // not in the c++ standard
      if(corrVec[i].second <= -1 || corrVec[i].second >= 1 || !(corrVec[i].second==corrVec[i].second)) {
        //cerr << "Warning bad correlation:" <<  corrVec[i].second << endl;
        corrVec[i].second = -2;
      }
      else 
        corrValidProbes++;
    }
    else {
      corrVec[i].second = -2; // safe val of -2 is lower than minimum possible correlation of -1
    }
  }

  /* Sort to find best correlated features. If there aren't enough
   valid correlation measurements just use the original order. Don't
   want to take a chance that different sorts will return different
   orders for same values on different platforms. */
  if(corrValidProbes >= featureNumToUse)
    std::sort(corrVec.rbegin(), corrVec.rend(), IterPlier::less_pair());
  minVal = corrVec[featureNumToUse-1].second;

  /* Walk through and identify best correlated features, copying them
     into our memory caches. */
  for(i = 0; i < atomCount ; i++) {
    if(passCount < featureNumToUse) {
      (atomUsed)[corrVec[i].first] = 1;
      for(j = 0; j < chipCount; j++) {
        m_pm[j][passCount] = pm[j][corrVec[i].first];
        m_mm[j][passCount] = mm[j][corrVec[i].first];
      }
      passCount++;
    }
    else {
      (atomUsed)[corrVec[i].first] = 0;
    }
  }

  /* Run plier only on our best correlated features. */
  plier.setNumFeature(passCount);
  plier.setPM(m_pm);
  plier.setMM(m_mm);
  plier.setResiduals(residuals);
  plier.setFeatureResponse(featureResponse);
  plier.setTargetResponse(targetResponse);
  plier.setNumExp(chipCount);
  plier.run(&errorCode);

  /* Move the residual and feature responses back so they
     are synched with the original arrays. */
  passCount = featureNumToUse -1;
  for(i = atomCount - 1; i >= 0; i--) {
    if(atomUsed[i]) {
      /* This works to copy all feature response and residuals back to
         original position in arrays as i is always greater than or
         equal atomCount so we don't overwrite the actual data by
         mistake. */
      assert(passCount >= 0);
      if(featureResponse != NULL && passCount != i) {
        featureResponse[i] = featureResponse[passCount];
        featureResponse[passCount] = IterPlier::noData;
      }
      if(residuals != NULL && passCount != i) {
        for(j = 0; j < chipCount; j++) {
          residuals[j][i] = residuals[j][passCount];
          residuals[j][passCount] = IterPlier::noData;
        }
      }
      passCount--;
    }
    else { /* Atom not used, fill in bad data. */
      if(featureResponse != NULL)
        featureResponse[i] = IterPlier::noData;
      if(residuals != NULL) {
        for(j = 0; j < chipCount; j++) 
          residuals[j][i] = IterPlier::noData;
      }
    }
  }

  return errorCode;
}
                                 
  
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
long IterPlier::runCorrIterPlier(caffyplier &plier, int atomCount, int chipCount, double** pm, double **mm, 
                      double *targetResponse, double *featureResponse, 
                      double **residuals, int *atomUsed, int iterCount, int *iterFeature) {
  long errorCode; /* plier will fill this in with status. */
  int i = 0;
  /* Some sanity checks. */
  assert(pm != NULL);
  assert(targetResponse != NULL);
  assert(featureResponse != NULL);
  assert(iterCount >= 0);
  assert(iterFeature != NULL);
  assert(atomCount > 0);
  assert(chipCount > 0);
  assert(atomUsed != NULL);
  for(i = 1; i < iterCount; i++) 
    assert(iterFeature[i-1] >= iterFeature[i]);
  /* Initialize atomUsed to true (1), all features used in
     first run. */

  // cerr << "Doing atom used init." << endl;
  for(i = 0; i < atomCount; i++) 
    (atomUsed)[i] = 1;

  /* Always going to run plier at least once. */
  // cerr << "Doing some settings." << endl;
  plier.setNumFeature(atomCount);
  plier.setPM(pm);
  plier.setMM(mm);
  plier.setResiduals(residuals);
  plier.setTargetResponse(targetResponse);
  plier.setFeatureResponse(featureResponse);
  plier.setNumExp(chipCount);

  /* Heavy lifting happens here. */
  // cerr << "Running plier." << endl;
  plier.run(&errorCode);
  /* If plier had a problem or we already are below the minimum number
     of features return here. */
  if(errorCode != 0 || iterCount == 0 || variance(targetResponse, targetResponse+chipCount) == 0)
    return errorCode;

  /* Do the iterations here. */
  for(i = 0; i < iterCount; i++) {
    if(atomCount > iterFeature[i]) {
      // cerr << "Running iterCoorPlier." << endl;
    	errorCode = iterateCorrPlier(plier, atomCount, chipCount, pm, mm, targetResponse, featureResponse,
                                 residuals, atomUsed, iterFeature[i]);
    }
    if(errorCode != 0)
      return errorCode;
  }
  return errorCode;
}
