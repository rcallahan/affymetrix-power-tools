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
 * @file   QuantPlierParams.h
 * @author Pete Klosterman
 * @date   Tue Feb 21 14:26:27 2006
 * 
 * @brief Class for plier paramaters.
 */
#ifndef _QUANTPLIERPARAMS_H_
#define _QUANTPLIERPARAMS_H_

/**
 *  Class for plier parameters.
 */
class QuantPlierParams {

public:

  /**
   * Constructor. 
   */
  QuantPlierParams() {
    m_InitAugmentation = 0;
    m_InitDefaultFeatureResponse = 0;
    m_InitDefaultTargetResponse = 0;
    m_SeaAttenuation = 0;
    m_SeaOptConvergence = 0;
    m_SeaOptIteration = 0;
    m_PlierGmCutoff = 0;
    m_PlierDifferentialFeaturePenalty = 0;
    m_PlierDifferentialTargetPenalty = 0;
    m_PlierUseMMLikelihood = 0;
    m_PlierUseInputModel = 0;
    m_PlierFitFeatureResponse = 0;
    m_PlierOptConvergence = 0;
    m_PlierOptIteration = 0;
    m_PlierOptDropMax = 0;
    m_PlierOptLambdaLimit = 0;
    m_PlierOptBalanceMethod = 0;
    m_PlierOptOptimizationMethod = 0;
    m_FixPrecomputed = 0;
    m_SafetyZero = 0;
    m_NumericalTolerance = 0;
    m_FixFeatureEffect = 0;
  }

  /** 
   * @brief Set InitAugmentation parameter.
   */
  inline void setInitAugmentation (double InitAugmentation) {
    m_InitAugmentation = InitAugmentation;
  }
  /** 
   * @brief Get InitAugmentation parameter.
   */
  inline double getInitAugmentation () {
    return m_InitAugmentation;
  }
  /** 
   * @brief Set InitDefaultFeatureResponse parameter.
   */
  inline void setInitDefaultFeatureResponse (double InitDefaultFeatureResponse) {
    m_InitDefaultFeatureResponse = InitDefaultFeatureResponse;
  }
  /** 
   * @brief Get InitDefaultFeatureResponse parameter.
   */
  inline double getInitDefaultFeatureResponse () {
    return m_InitDefaultFeatureResponse;
  }
  /** 
   * @brief Set InitDefaultTargetResponse parameter.
   */
  inline void setInitDefaultTargetResponse (double InitDefaultTargetResponse) {
    m_InitDefaultTargetResponse = InitDefaultTargetResponse;
  }
  /** 
   * @brief Get InitDefaultTargetResponse parameter.
   */
  inline double getInitDefaultTargetResponse () {
    return m_InitDefaultTargetResponse;
  }
  /** 
   * @brief Set SeaAttenuation parameter.
   */
  inline void setSeaAttenuation (float SeaAttenuation) {
    m_SeaAttenuation = SeaAttenuation;
  }
  /** 
   * @brief Get SeaAttenuation parameter.
   */
  inline float getSeaAttenuation () {
    return m_SeaAttenuation;
  }
  /** 
   * @brief Set SeaOptConvergence parameter.
   */
  inline void setSeaOptConvergence (double SeaOptConvergence) {
    m_SeaOptConvergence = SeaOptConvergence;
  }
  /** 
   * @brief Get SeaOptConvergence parameter.
   */
  inline double getSeaOptConvergence () {
    return m_SeaOptConvergence;
  }
  /** 
   * @brief Set SeaOptIteration parameter.
   */
  inline void setSeaOptIteration (long SeaOptIteration) {
    m_SeaOptIteration = SeaOptIteration;
  }
  /** 
   * @brief Get SeaOptIteration parameter.
   */
  inline long getSeaOptIteration () {
    return m_SeaOptIteration;
  }
  /** 
   * @brief Set PlierGmCutoff parameter.
   */
  inline void setPlierGmCutoff (float PlierGmCutoff) {
    m_PlierGmCutoff = PlierGmCutoff;
  }
  /** 
   * @brief Get PlierGmCutoff parameter.
   */
  inline float getPlierGmCutoff () {
    return m_PlierGmCutoff;
  }
  /** 
   * @brief Set PlierDifferentialFeaturePenalty parameter.
   */
  inline void setPlierDifferentialFeaturePenalty (float PlierDifferentialFeaturePenalty) {
    m_PlierDifferentialFeaturePenalty = PlierDifferentialFeaturePenalty;
  }
  /** 
   * @brief Get PlierDifferentialFeaturePenalty parameter.
   */
  inline float getPlierDifferentialFeaturePenalty () {
    return m_PlierDifferentialFeaturePenalty;
  }
  /** 
   * @brief Set PlierDifferentialTargetPenalty parameter.
   */
  inline void setPlierDifferentialTargetPenalty (float PlierDifferentialTargetPenalty) {
    m_PlierDifferentialTargetPenalty = PlierDifferentialTargetPenalty;
  }
  /** 
   * @brief Get PlierDifferentialTargetPenalty parameter.
   */
  inline float getPlierDifferentialTargetPenalty () {
    return m_PlierDifferentialTargetPenalty;
  }
  /** 
   * @brief Set PlierUseMMLikelihood parameter.
   */
  inline void setPlierUseMMLikelihood (bool PlierUseMMLikelihood) {
    m_PlierUseMMLikelihood = PlierUseMMLikelihood;
  }
  /** 
   * @brief Get PlierUseMMLikelihood parameter.
   */
  inline bool getPlierUseMMLikelihood () {
    return m_PlierUseMMLikelihood;
  }
  /** 
   * @brief Set PlierUseInputModel parameter.
   */
  inline void setPlierUseInputModel (bool PlierUseInputModel) {
    m_PlierUseInputModel = PlierUseInputModel;
  }
  /** 
   * @brief Get PlierUseInputModel parameter.
   */
  inline bool getPlierUseInputModel () {
    return m_PlierUseInputModel;
  }
  /** 
   * @brief Set PlierFitFeatureResponse parameter.
   */
  inline void setPlierFitFeatureResponse (bool PlierFitFeatureResponse) {
    m_PlierFitFeatureResponse = PlierFitFeatureResponse;
  }
  /** 
   * @brief Get PlierFitFeatureResponse parameter.
   */
  inline bool getPlierFitFeatureResponse () {
    return m_PlierFitFeatureResponse;
  }
  /** 
   * @brief Set PlierOptConvergence parameter.
   */
  inline void setPlierOptConvergence (double PlierOptConvergence) {
    m_PlierOptConvergence = PlierOptConvergence;
  }
  /** 
   * @brief Get PlierOptConvergence parameter.
   */
  inline double getPlierOptConvergence () {
    return m_PlierOptConvergence;
  }
  /** 
   * @brief Set PlierOptIteration parameter.
   */
  inline void setPlierOptIteration (long PlierOptIteration) {
    m_PlierOptIteration = PlierOptIteration;
  }
  /** 
   * @brief Get PlierOptIteration parameter.
   */
  inline long getPlierOptIteration () {
    return m_PlierOptIteration;
  }
  /** 
   * @brief Set PlierOptDropMax parameter.
   */
  inline void setPlierOptDropMax (double PlierOptDropMax) {
    m_PlierOptDropMax = PlierOptDropMax;
  }
  /** 
   * @brief Get PlierOptDropMax parameter.
   */
  inline double getPlierOptDropMax () {
    return m_PlierOptDropMax;
  }
  /** 
   * @brief Set PlierOptLambdaLimit parameter.
   */
  inline void setPlierOptLambdaLimit (double PlierOptLambdaLimit) {
    m_PlierOptLambdaLimit = PlierOptLambdaLimit;
  }
  /** 
   * @brief Get PlierOptLambdaLimit parameter.
   */
  inline double getPlierOptLambdaLimit () {
    return m_PlierOptLambdaLimit;
  }
  /** 
   * @brief Set PlierOptOptimizationMethod parameter.
   */
  inline void setPlierOptOptimizationMethod (long PlierOptOptimizationMethod) {
    m_PlierOptOptimizationMethod = PlierOptOptimizationMethod;
  }
  /** 
   * @brief Get PlierOptOptimizationMethod parameter.
   */
  inline long getPlierOptOptimizationMethod () {
    return m_PlierOptOptimizationMethod;
  }
  /** 
   * @brief Set PlierOptBalanceMethod parameter.
   */
  inline void setPlierOptBalanceMethod (long PlierOptBalanceMethod) {
    m_PlierOptBalanceMethod = PlierOptBalanceMethod;
  }
  /** 
   * @brief Get PlierOptBalanceMethod parameter.
   */
  inline long getPlierOptBalanceMethod () {
    return m_PlierOptBalanceMethod;
  }
  /** 
   * @brief Set FixPrecomputed parameter.
   */
  inline void setFixPrecomputed (bool FixPrecomputed ) {
    m_FixPrecomputed = FixPrecomputed;
  }
  /** 
   * @brief Get FixPrecomputed parameter.
   */
  inline bool getFixPrecomputed () {
    return m_FixPrecomputed;
  }
  /** 
   * @brief Set NumericalTolerance parameter.
   */
  inline void setNumericalTolerance (double NumericalTolerance) {
    m_NumericalTolerance = NumericalTolerance;
  }
  /** 
   * @brief Get NumericalTolerance parameter.
   */
  inline double getNumericalTolerance () {
    return m_NumericalTolerance;
  }
  /** 
   * @brief Set SafetyZero parameter.
   */
  inline void setSafetyZero (double SafetyZero) {
    m_SafetyZero = SafetyZero;
  }
  /** 
   * @brief Get SafetyZero parameter.
   */
  inline double getSafetyZero () {
    return m_SafetyZero;
  }
  /** 
   * @brief Set FixFeatureEffect parameter.
   */
  inline void setFixFeatureEffect (bool FixFeatureEffect) {
    m_FixFeatureEffect = FixFeatureEffect;
  }
  /** 
   * @brief Get FixFeatureEffect parameter.
   */
  inline bool getFixFeatureEffect() {
    return m_FixFeatureEffect;
  }

private:

  /// Value of the InitAugmentation parameter.
  double m_InitAugmentation;
  /// Value of the InitDefaultFeatureResponse parameter.
  double m_InitDefaultFeatureResponse;
  /// Value of the InitDefaultTargetResponse parameter.
  double m_InitDefaultTargetResponse;
  /// Value of the SeaAttenuation parameter.
  /// What attenuation to use for SEA background adjustment (parameter
  /// for how much to adjust for background signal). 1 is like using
  /// PM only, 0 is like using a hard threshold.
  float m_SeaAttenuation;
  /// Value of the SeaOptConvergence parameter.
  double m_SeaOptConvergence;
  /// Value of the SeaOptIteration parameter.
  long m_SeaOptIteration;
  /// Value of the PlierGmCutoff parameter.
  float m_PlierGmCutoff;
  /// Value of the PlierDifferentialFeaturePenalty parameter.
  float m_PlierDifferentialFeaturePenalty;
  /// Value of the PlierDifferentialTargetPenalty parameter.
  float m_PlierDifferentialTargetPenalty;
  /// Value of the PlierUseMMLikelihood parameter.
  bool m_PlierUseMMLikelihood;
  /// Value of the PlierUseInputModel parameter.
  bool m_PlierUseInputModel;
  /// Value of the PlierFitFeatureResponse parameter.
  bool m_PlierFitFeatureResponse;
  /// Value of the PlierOptConvergence parameter.
  double m_PlierOptConvergence;
  /// Value of the PlierOptIteration parameter.
  long m_PlierOptIteration;
  /// Value of the PlierOptDropMax parameter.
  double m_PlierOptDropMax;
  /// Value of the PlierOptLambdaLimit parameter.
  double m_PlierOptLambdaLimit;
  /// Value of the PlierOptOptimizationMethod parameter.
  /// What optimization level are we using for plier (0 - full plier, 1 - SEA)
  long m_PlierOptOptimizationMethod;
  /// Value of the PlierOptBalanceMethod parameter.
  long m_PlierOptBalanceMethod;
  /// Value of FixPrecomputed
  bool m_FixPrecomputed;
  /// Value of SafetyZero
  double m_SafetyZero;
  /// Value of NumericalTolerance
  double m_NumericalTolerance;
  /// Value of FixFeatureEffect
  bool m_FixFeatureEffect;


};

#endif /* _QUANTPLIERPARAMS_H_ */
