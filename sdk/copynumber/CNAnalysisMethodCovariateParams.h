////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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
 * @file   CNAnalysisMethodCovariateParams.h
 *
 * @brief  Class containing common data for covariate-based signal and log2ratio adjustment
 */
#ifndef _CNAnalysisMethodCovariateParams_H_
#define _CNAnalysisMethodCovariateParams_H_

//
#include "copynumber/CNAnalysisMethod.h"
//
#include "util/AffxByteArray.h"
#include "util/AffxConv.h"
#include "util/AffxMultiDimensionalArray.h"
#include "util/Convert.h"
#include "util/Util.h"
//
#include <vector>
//

/**
  * @brief A class to process the covariate command line
  */
class CovariateParams
{
public:

  static const std::string m_CNReferencePrefix;

  static std::map<string, int> m_allCovariateMap;       // all CovariateEnum-s from the cmd line listed in the order of
                                                        // CNProbeSet::m_vAllCovariates (which is the same as the order
                                                        // of covariates in the reference)

  static std::vector<int> m_vSignalCovariates;          // indices into CNProbeSet::m_vAllCovariates showing the cmd line order
  static std::vector<int> m_vLRCovariates;              // indices into CNProbeSet::m_vAllCovariates showing the cmd line order
  static std::vector<int> m_vAPCovariates;              // indices into CNProbeSet::m_vAllCovariates showing the cmd line order

public:
  static void checkParams(BaseEngine& engine);

  static void determineCovariateMap(BaseEngine& engine);
  static void restoreCovariateParams(std::string str);
  static int mapCovariateNameToIndex(const std::string& str);
  static std::string translateToParamString(const std::map<string, int>& covMap);
  static void translateFromParamString(const std::string& str, std::map<std::string, int>& retMap);
  static bool isReservedCovariateName(const std::string& str);

  // Signal methods
  static bool useEquallyPopulatedSignalBins(int covariateIndex);
  static bool useEquallySpacedSignalBins(int covariateIndex);
  static bool isSignalCovariateDiscrete(int covariateIndex);
  static int getNumSignalBins(int covariateIndex);

  // LR methods
  static bool useEquallyPopulatedLRBins(int covariateIndex);
  static bool useEquallySpacedLRBins(int covariateIndex);
  static bool isLRCovariateDiscrete(int covariateIndex);
  static int getNumLRBins(int covariateIndex);

  // AP methods
  static bool useEquallyPopulatedAPBins(int covariateIndex);
  static bool useEquallySpacedAPBins(int covariateIndex);
  static bool isAPCovariateDiscrete(int covariateIndex);
  static int getNumAPBins(int covariateIndex);

private:
  static void initializeSignalBinTypeAndNum(const std::string& signalStr);
  static void initializeLRBinTypeAndNum(const std::string& lrParameter);
  static void initializeAPBinTypeAndNum(const std::string& apParameter);
  static void initializeBinTypeAndNum(const std::string& parameterString, std::vector<std::string>& vType, std::vector<int>& vNum);
  static void initializeAdditionalLRParams(const string& lrParameter);
  static void initializeAdditionalAPParams(const string& apParameter);

  static bool doesBinningTypeMatch(int index, std::vector<std::string>& vBinTypes, const string& type);
  static int getNumBins(int index, const std::vector<int>& vBinSizes);
  static void clear();
  static bool anyFileBasedCovariates();

  static void addCovariatesToMap(vector<string>& cov, vector<int>& arr);
  static void paramValueToIntVector(map<string, string>& params, const string& paramName, vector<int>& vInt);
  static void paramValueToDoubleVector(map<string, string>& params, const string& paramName, vector<double>& vDouble);
  static int paramValueToInt(map<string, string>& params, const string& paramName, int defaultValue);
  static double paramValueToDouble(map<string, string>& params, const string& paramName, double defaultValue);
  static void elementCountMustMatch(vector<int>& numElements, const string& errmsg);
  static bool anyCoarseAPAdjustments();

public:
  static std::vector<std::string> m_vSignalBinningType;
  static std::vector<int> m_vSignalBinSizes;
  static std::vector<std::string> m_vLRBinningType;
  static std::vector<int> m_vLRBinSizes;
  static std::vector<std::string> m_vIQRScaling;
  static std::vector<std::string> m_vLRSubractFromXY;
  static std::vector<std::string> m_vAPBinningType;
  static std::vector<int> m_vAPBinSizes;
  static std::vector<std::string> m_vCoarseAPAdjust;
  static std::vector<int> m_vCoarseAPAdjustStep;
  static std::vector<int> m_vCoarseAPAdjustWindow;
  static std::vector<int> m_vCoarseAPAdjustPointCount;
  static std::vector<double> m_vCoarseAPAdjustBandwidth;
  static std::vector<double> m_vCoarseAPAdjustCutoff;
  static std::vector<double> m_vCoarseAPAdjustCleanthreshold;
  static std::vector<double> m_vCoarseAPAdjustOutliertrim;
  static int m_masterPeakPointCount;
  static double m_masterPeakBandwidth;
  static int m_covariatePeakPointCount;
  static double m_covariatePeakBandwidth;
  static std::string m_kernelFunctionSelection;
  static std::string m_signalCommand;
  static std::string m_LRCommand;
  static std::string m_APCommand;

};

#endif  //_CNAnalysisMethodCovariateParams_H_
