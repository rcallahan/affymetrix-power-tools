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
 * @file   CNAnalysisMethodCovariateParams.cpp.cpp
 *
 * @brief  Class containing common data for covariate-based signal and log2ratio adjustment
 */

//
#include "copynumber/CNAnalysisMethodCovariateParams.h"
#include "copynumber/Annotation.h"
//
#include "stats/stats.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Verbose.h"
#include "util/Fs.h"
//

using namespace std;

const std::string CovariateParams::m_CNReferencePrefix("#%affymetrix-algorithm-param-covariates-");

std::map<std::string, int> CovariateParams::m_allCovariateMap;
std::vector<int> CovariateParams::m_vSignalCovariates;
std::vector<int> CovariateParams::m_vLRCovariates;
std::vector<int> CovariateParams::m_vAPCovariates;

std::vector<std::string> CovariateParams::m_vSignalBinningType;
std::vector<int> CovariateParams::m_vSignalBinSizes;
std::vector<std::string> CovariateParams::m_vLRBinningType;
std::vector<int> CovariateParams::m_vLRBinSizes;
std::vector<std::string> CovariateParams::m_vIQRScaling;
std::vector<std::string> CovariateParams::m_vLRSubractFromXY;
std::vector<std::string> CovariateParams::m_vAPBinningType;
std::vector<int> CovariateParams::m_vAPBinSizes;
std::vector<std::string> CovariateParams::m_vCoarseAPAdjust;
std::vector<int> CovariateParams::m_vCoarseAPAdjustStep;
std::vector<int> CovariateParams::m_vCoarseAPAdjustWindow;
std::vector<int> CovariateParams::m_vCoarseAPAdjustPointCount;
std::vector<double> CovariateParams::m_vCoarseAPAdjustBandwidth;
std::vector<double> CovariateParams::m_vCoarseAPAdjustCutoff;
std::vector<double> CovariateParams::m_vCoarseAPAdjustCleanthreshold;
std::vector<double> CovariateParams::m_vCoarseAPAdjustOutliertrim;
int CovariateParams::m_masterPeakPointCount;
double CovariateParams::m_masterPeakBandwidth;
int CovariateParams::m_covariatePeakPointCount;
double CovariateParams::m_covariatePeakBandwidth;
string CovariateParams::m_kernelFunctionSelection;
std::string CovariateParams::m_signalCommand;
std::string CovariateParams::m_LRCommand;
std::string CovariateParams::m_APCommand;

/******************************************************************/
/*              CovariateParams class                             */

void CovariateParams::clear()
{
    m_allCovariateMap.clear();
    m_vSignalCovariates.clear();
    m_vLRCovariates.clear();
    m_vAPCovariates.clear();
    m_vSignalBinningType.clear();
    m_vSignalBinSizes.clear();
    m_vLRBinningType.clear();
    m_vLRBinSizes.clear();
    m_vIQRScaling.clear();
    m_vLRSubractFromXY.clear();
    m_vAPBinningType.clear();
    m_vAPBinSizes.clear();

    m_vCoarseAPAdjust.clear();
    m_vCoarseAPAdjustStep.clear();
    m_vCoarseAPAdjustWindow.clear();
    m_vCoarseAPAdjustPointCount.clear();
    m_vCoarseAPAdjustBandwidth.clear();
    m_vCoarseAPAdjustCutoff.clear();
    m_vCoarseAPAdjustCleanthreshold.clear();
    m_vCoarseAPAdjustOutliertrim.clear();
    m_masterPeakPointCount = -1;
    m_masterPeakBandwidth = -1;
    m_covariatePeakPointCount = -1;
    m_covariatePeakBandwidth = -1;
    m_kernelFunctionSelection.clear();
}

void CovariateParams::checkParams(BaseEngine& engine)
{
    // Use the determineCovariateMap method to run the ususal parsing.
    // This may throw an exception.
    determineCovariateMap(engine);

    // Check that the map, bin-type, and bin-count dimensions agree.
    //if (m_vSignalCovariates.size() != m_vSignalBinningType.size())
    //{
    //    Err::errAbort("signal-adjustment-covariates option error - the number elements for order, bin-type and bin-count must match");
    //}

    //if (m_vSignalBinningType.size() != m_vSignalBinSizes.size())
    //{
    //    Err::errAbort("signal-adjustment-covariates option error - the number elements for order, bin-type and bin-count must match");
    //}

    vector<int> numElements;
    numElements.push_back(m_vSignalCovariates.size());
    numElements.push_back(m_vSignalBinningType.size());
    numElements.push_back(m_vSignalBinSizes.size());
    elementCountMustMatch(numElements, "signal-adjustment-covariates option error - the number elements for order, bin-type and bin-count must match");

    //if (m_vLRCovariates.size() != m_vLRBinningType.size())
    //{
    //    Err::errAbort("lr-adjustment-covariates option error - the number elements for order, bin-type and bin-count must match");
    //}

    //if (m_vLRBinningType.size() != m_vLRBinSizes.size())
    //{
    //    Err::errAbort("lr-adjustment-covariates option error - the number elements for order, bin-type and bin-count must match");
    //}

    //if (m_vLRBinSizes.size() != m_vIQRScaling.size())
    //{
    //    Err::errAbort("lr-adjustment-covariates option error - the number elements for order, bin-type, bin-count and iqr-scaling must match");
    //}

    //if (m_vIQRScaling.size() != m_vLRSubractFromXY.size())
    //{
    //    Err::errAbort("lr-adjustment-covariates option error - the number elements for order, bin-type, bin-count, iqr-scaling and subtract-from-XY must match");
    //}

    numElements.clear();
    numElements.push_back(m_vLRCovariates.size());
    numElements.push_back(m_vLRBinningType.size());
    numElements.push_back(m_vLRBinSizes.size());
    numElements.push_back(m_vIQRScaling.size());
    numElements.push_back(m_vLRSubractFromXY.size());
    elementCountMustMatch(numElements, "lr-adjustment-covariates option error - the number elements for order, bin-type, bin-count, iqr-scaling and subtract-from-XY must match");


    //if (m_vAPCovariates.size() != m_vAPBinningType.size())
    //{
    //    Err::errAbort("allele-peaks-adjustment-covariates option error - the number elements for order, bin-type and bin-count must match");
    //}

    //if (m_vAPBinningType.size() != m_vAPBinSizes.size())
    //{
    //    Err::errAbort("allele-peaks-adjustment-covariates option error - the number elements for order, bin-type and bin-count must match");
    //}

    numElements.clear();
    numElements.push_back(m_vAPCovariates.size());
    numElements.push_back(m_vAPBinningType.size());
    numElements.push_back(m_vAPBinSizes.size());
    elementCountMustMatch(numElements, "allele-peaks-adjustment-covariates option error - the number elements for order, bin-type and bin-count must match");

    if (anyCoarseAPAdjustments())
    {
        numElements.clear();
        numElements.push_back(m_vAPCovariates.size());
        numElements.push_back(m_vCoarseAPAdjust.size());
        numElements.push_back(m_vCoarseAPAdjustStep.size());
        numElements.push_back(m_vCoarseAPAdjustWindow.size());
        numElements.push_back(m_vCoarseAPAdjustPointCount.size());
        numElements.push_back(m_vCoarseAPAdjustBandwidth.size());
        numElements.push_back(m_vCoarseAPAdjustCutoff.size());
        numElements.push_back(m_vCoarseAPAdjustCleanthreshold.size());
        numElements.push_back(m_vCoarseAPAdjustOutliertrim.size());

        // Err::errAbort("allele-peaks-adjustment-covariates option error - the number elements for any of the coarse-allele-peak-adjustment sub-options, if used, must match the number of elements for order, bin-type and bin-count");
        elementCountMustMatch(numElements, "allele-peaks-adjustment-covariates option error - the number elements for any of the coarse-allele-peak-adjustment sub-options, if used, must match the number of elements for order, bin-type and bin-count");
    }

    for (int i = 0; i < m_vSignalBinSizes.size(); ++i)
    {
        if (getNumSignalBins(i) <= 0 && isSignalCovariateDiscrete(i) == false) {Err::errAbort("signal-adjustment-covariates option error - bin-count may not be less than or equal to zero for a non-discrete bin type.");}
    }

    for (int i = 0; i < m_vLRBinSizes.size(); ++i)
    {
        if (getNumLRBins(i) <= 0 && isLRCovariateDiscrete(i) == false) Err::errAbort("lr-adjustment-covariates option error - bin-count may not be less than or equal to zero for a non-discrete bin type.");
    }

    for (int i = 0; i < m_vAPBinSizes.size(); ++i)
    {
        if (getNumAPBins(i) <= 0 && isAPCovariateDiscrete(i) == false) Err::errAbort("allele-peaks-adjustment-covariates option error - bin-count may not be less than or equal to zero for a non-discrete bin type.");
    }

    if (anyFileBasedCovariates() == true)
    {
        if (engine.getOpt("covariates-file") == "")
        {
            Err::errAbort("The lr-adjustment-covariates or signal-adjustment-covariates options refer to one or more file-based covariates but the covariates-file option is missing or empty.");
        }
    }

    if ((engine.getOpt("covariates-file") != "") && (!Fs::fileExists(engine.getOpt("covariates-file"))))
    {
        Err::errAbort("Must specify a valid covariates-file.");
    }

    // Don't leave any residual state
    clear();
}

bool CovariateParams::anyCoarseAPAdjustments()
{
    bool result = false;
    for (vector<string>::iterator ii = m_vCoarseAPAdjust.begin(); ii != m_vCoarseAPAdjust.end(); ++ii)
    {
        if (*ii != "off") { result = true; break; }
    }
    return result;
}

void CovariateParams::elementCountMustMatch(vector<int>& numElements, const string& errmsg)
{
    int num = 0;
    vector<int>::iterator ii = numElements.begin();
    if (ii != numElements.end()) { num = *ii; ++ii; }
    for (; ii != numElements.end(); ++ii)
    {
        if (*ii != num) { Err::errAbort(errmsg); }
    }
}

/** Initializes the signal and LR covariate map members.
    The value represents the index in CNProbeSet::getCovariateValue().
    The index is the CovariateEnum.
    Side effect: this method also initializes the number of bins and binning type paramters.
  */
void CovariateParams::determineCovariateMap(BaseEngine& engine)
{
    clear();

    string signalStr;
    vector<int> signalCovariateIndexes;
    if (engine.isOptDefined("signal-adjustment-covariates")) signalStr = engine.getOpt("signal-adjustment-covariates");
    if (signalStr.empty() == false)
    {
        string name;
        map<string, string> param;
        SelfCreate::fillInNameParam(signalStr, name, param);

        map<string, string>::iterator iter;
        if ((iter = param.find("order")) != param.end())
        {
            vector<string> order = StringUtils::Split(iter->second, ",");
            addCovariatesToMap(order, m_vSignalCovariates);
        }
    }

    string lrString;
    vector<int> lrCovariateIndexes;
    if (engine.isOptDefined("lr-adjustment-covariates")) lrString = engine.getOpt("lr-adjustment-covariates");
    if (lrString.empty() == false)
    {
        string name;
        map<string, string> param;
        SelfCreate::fillInNameParam(lrString, name, param);

        map<string, string>::iterator iter;
        if ((iter = param.find("order")) != param.end())
        {
            vector<string> order = StringUtils::Split(iter->second, ",");
            addCovariatesToMap(order, m_vLRCovariates);
        }
    }

    string apString;
    vector<int> apCovariateIndexes;
    if (engine.isOptDefined("allele-peaks-adjustment-covariates")) apString = engine.getOpt("allele-peaks-adjustment-covariates");
    if (apString.empty() == false)
    {
        string name;
        map<string, string> param;
        SelfCreate::fillInNameParam(apString, name, param);

        map<string, string>::iterator iter;
        if ((iter = param.find("order")) != param.end())
        {
            vector<string> order = StringUtils::Split(iter->second, ",");
            addCovariatesToMap(order, m_vAPCovariates);
        }
    }

    initializeSignalBinTypeAndNum(signalStr);
    initializeLRBinTypeAndNum(lrString);
    initializeAPBinTypeAndNum(apString);
    initializeAdditionalLRParams(lrString);
    initializeAdditionalAPParams(apString);

    Verbose::out(4, "CovariateParams::determineCovariateMap - num of covariates to load: " + Convert::toString((int)m_allCovariateMap.size()));
}

/**  Restore covariate parameter from the Parameters section of the CN reference to the
  *  corresponding class member:
  *
  * "all"                           -> m_allCovariateMap
  * "signal"                        -> m_vSignalCovariates
  * "log2ratio                      -> m_vLRCovariates
  * "signal-binning-types"          -> m_vSignalBinningType
  * "signal-bin-sizes"              -> m_vSignalBinSizes
  * "log2ratio-binning-types"       -> m_vLRBinningType
  * "log2ratio-bin-sizes"           -> m_vLRBinSizes
  * "log2ratio-iqr-scaling"         -> m_vIQRScaling
  * "log2ratio-subtract-from-XY"    -> m_vLRSubractFromXY
  */
void CovariateParams::restoreCovariateParams(std::string str)
{
    str = str.substr(m_CNReferencePrefix.size());   // chop the prefix off
    int iFindIndex = str.find('=');
    string covVal  = str.substr(iFindIndex + 1);
    if (covVal.empty()) {
        return;
    }

    string covName = str.substr(0, iFindIndex);
    if (covName == "all") {
        translateFromParamString(covVal, m_allCovariateMap);
    } else if (covName == "signal") {
        Convert::strToIntVec(covVal, ',', m_vSignalCovariates);
    } else if (covName == "log2ratio") {
        Convert::strToIntVec(covVal, ',', m_vLRCovariates);
    } else if (covName == "allele-peaks") {
        Convert::strToIntVec(covVal, ',', m_vAPCovariates);
    } else if (covName == "signal-bin-sizes") {
        Convert::strToIntVec(covVal, ',', m_vSignalBinSizes);
    } else if (covName == "log2ratio-bin-sizes") {
        Convert::strToIntVec(covVal, ',', m_vLRBinSizes);
    } else if (covName == "allele-peaks-bin-sizes") {
        Convert::strToIntVec(covVal, ',', m_vAPBinSizes);
    } else if (covName == "signal-binning-types") {
        Convert::strToStrVec(covVal, ',', m_vSignalBinningType);
    } else if (covName == "log2ratio-binning-types") {
        Convert::strToStrVec(covVal, ',', m_vLRBinningType);
    } else if (covName == "allele-peaks-binning-types") {
        Convert::strToStrVec(covVal, ',', m_vAPBinningType);
    } else if (covName == "log2ratio-iqr-scaling") {
        Convert::strToStrVec(covVal, ',', m_vIQRScaling);
    } else if (covName == "log2ratio-subtract-from-XY") {
        Convert::strToStrVec(covVal, ',', m_vLRSubractFromXY);
    } else if (covName == "allele-peaks-coarse-adjustment") {
        Convert::strToStrVec(covVal, ',', m_vCoarseAPAdjust);
    } else if (covName == "coarse-allele-peak-adjustment-step") {
        Convert::strToIntVec(covVal, ',', m_vCoarseAPAdjustStep);
    } else if (covName == "coarse-allele-peak-adjustment-window") {
        Convert::strToIntVec(covVal, ',', m_vCoarseAPAdjustWindow);
    } else if (covName == "coarse-allele-peak-adjustment-point-count") {
        Convert::strToIntVec(covVal, ',', m_vCoarseAPAdjustPointCount);
    } else if (covName == "coarse-allele-peak-adjustment-bandwidth") {
        Convert::strToDoubleVec(covVal, ',', m_vCoarseAPAdjustBandwidth);
    } else if (covName == "coarse-allele-peak-adjustment-cutoff") {
        Convert::strToDoubleVec(covVal, ',', m_vCoarseAPAdjustCutoff);
    } else if (covName == "coarse-allele-peak-adjustment-clean-threshold") {
        Convert::strToDoubleVec(covVal, ',', m_vCoarseAPAdjustCleanthreshold);
    } else if (covName == "coarse-allele-peak-adjustment-outlier-trim") {
        Convert::strToDoubleVec(covVal, ',', m_vCoarseAPAdjustOutliertrim);
    } else if (covName == "master-peaks-point-count") {
        m_masterPeakPointCount = Convert::toInt(covVal);
    } else if (covName == "master-peaks-bandwidth") {
        m_masterPeakBandwidth = Convert::toDouble(covVal);
    } else if (covName == "covariate-bin-peaks-point-count") {
        m_covariatePeakPointCount = Convert::toInt(covVal);
    } else if (covName == "covariate-bin-peaks-bandwidth") {
        m_covariatePeakBandwidth = Convert::toDouble(covVal);
    } else if (covName == "kernel-function-selection") {
        m_kernelFunctionSelection = covVal;
    } else if (covName == "signal-command") {
        m_signalCommand = covVal;
    } else if (covName == "log2ratio-command") {
        m_LRCommand = covVal;
    } else if (covName == "allele-peaks-command") {
        m_APCommand = covVal;
    }
}

int CovariateParams::mapCovariateNameToIndex(const std::string& str)
{
    map<string, int>::iterator ii = m_allCovariateMap.find(str);
    if (ii == m_allCovariateMap.end()) return -1;
    else return ii->second;
}

void CovariateParams::addCovariatesToMap(vector<string>& cov, vector<int>& arr)
{
    arr.clear();

    for (vector<string>::iterator ii = cov.begin(); ii != cov.end(); ++ii)
    {
        // Check if the covariate is already in the list
        if (m_allCovariateMap.find(*ii) == m_allCovariateMap.end())
        {
            int index = m_allCovariateMap.size();
            m_allCovariateMap[*ii] = index;
        }
        arr.push_back(m_allCovariateMap[*ii]);
    }
}

bool CovariateParams::useEquallyPopulatedSignalBins(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vSignalBinningType, "equally-populated");
}

bool CovariateParams::useEquallySpacedSignalBins(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vSignalBinningType, "equally-spaced");
}

bool CovariateParams::isSignalCovariateDiscrete(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vSignalBinningType, "discrete");
}

bool CovariateParams::useEquallyPopulatedLRBins(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vLRBinningType, "equally-populated");
}

bool CovariateParams::useEquallySpacedLRBins(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vLRBinningType, "equally-spaced");
}

bool CovariateParams::isLRCovariateDiscrete(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vLRBinningType, "discrete");
}

bool CovariateParams::useEquallyPopulatedAPBins(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vAPBinningType, "equally-populated");
}

bool CovariateParams::useEquallySpacedAPBins(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vAPBinningType, "equally-spaced");
}

bool CovariateParams::isAPCovariateDiscrete(int covariateIndex)
{
    return doesBinningTypeMatch(covariateIndex, m_vAPBinningType, "discrete");
}

bool CovariateParams::doesBinningTypeMatch(int index, std::vector<std::string>& vBinTypes, const string& type)
{
    if (vBinTypes.size() <= index)
    {
        Err::errAbort("Bin type index is out of bounds - " + Convert::toString(index));
    }
    return (vBinTypes[index].compare(type) == 0);
}

int CovariateParams::getNumSignalBins(int covariateIndex)
{
    return getNumBins(covariateIndex, m_vSignalBinSizes);
}

int CovariateParams::getNumLRBins(int covariateIndex)
{
    return getNumBins(covariateIndex, m_vLRBinSizes);
}

int CovariateParams::getNumAPBins(int covariateIndex)
{
    return getNumBins(covariateIndex, m_vAPBinSizes);
}

int CovariateParams::getNumBins(int index, const std::vector<int>& vBinSizes)
{
    if (vBinSizes.size() <= index)
    {
        Err::errAbort("Bin size index is out of bounds - " + Convert::toString(index));
    }
    return vBinSizes[index];
}

void CovariateParams::initializeSignalBinTypeAndNum(const string& signalStr)
{
    initializeBinTypeAndNum(signalStr, m_vSignalBinningType, m_vSignalBinSizes);
}
    
void CovariateParams::initializeLRBinTypeAndNum(const string& lrParameter)
{
    initializeBinTypeAndNum(lrParameter, m_vLRBinningType, m_vLRBinSizes);
}

void CovariateParams::initializeAPBinTypeAndNum(const string& apParameter)
{
    initializeBinTypeAndNum(apParameter, m_vAPBinningType, m_vAPBinSizes);
}

void CovariateParams::initializeAdditionalLRParams(const string& lrParameter)
{
    string name;
    map<string, string> param;
    SelfCreate::fillInNameParam(lrParameter, name, param);

    map<string, string>::iterator iter;
    if ((iter = param.find("iqr-scaling")) != param.end())
    {
        m_vIQRScaling = StringUtils::Split(iter->second, ",");
    }

    if ((iter = param.find("subtract-from-XY")) != param.end())
    {
        m_vLRSubractFromXY = StringUtils::Split(iter->second, ",");
    }
}

void CovariateParams::initializeAdditionalAPParams(const string& apParameter)
{
    string name;
    map<string, string> param;
    SelfCreate::fillInNameParam(apParameter, name, param);

    map<string, string>::iterator iter;
    if ((iter = param.find("coarse-allele-peak-adjustment")) != param.end())
    {
        m_vCoarseAPAdjust = StringUtils::Split(iter->second, ",");
    }
    else
    {
      // set the default to no
      for (int i = 0; i < m_vAPCovariates.size(); ++i) { m_vCoarseAPAdjust.push_back("off"); }
    }

    paramValueToIntVector(param, "coarse-allele-peak-adjustment-step", m_vCoarseAPAdjustStep);
    paramValueToIntVector(param, "coarse-allele-peak-adjustment-window", m_vCoarseAPAdjustWindow);
    paramValueToIntVector(param, "coarse-allele-peak-adjustment-point-count", m_vCoarseAPAdjustPointCount);
    paramValueToDoubleVector(param, "coarse-allele-peak-adjustment-bandwidth", m_vCoarseAPAdjustBandwidth);
    paramValueToDoubleVector(param, "coarse-allele-peak-adjustment-cutoff", m_vCoarseAPAdjustCutoff);
    paramValueToDoubleVector(param, "coarse-allele-peak-adjustment-clean-threshold", m_vCoarseAPAdjustCleanthreshold);
    paramValueToDoubleVector(param, "coarse-allele-peak-adjustment-outlier-trim", m_vCoarseAPAdjustOutliertrim);

    m_masterPeakPointCount = paramValueToInt(param, "master-peaks-point-count", -1);
    m_masterPeakBandwidth = paramValueToDouble(param, "master-peaks-bandwidth", -1.0);
    m_covariatePeakPointCount = paramValueToInt(param, "covariate-bin-peaks-point-count", -1);
    m_covariatePeakBandwidth = paramValueToDouble(param, "covariate-bin-peaks-bandwidth", -1.0);

    if ((iter = param.find("kernel-function-selection")) != param.end())
    {
        m_kernelFunctionSelection = iter->second;
    }
}

void CovariateParams::initializeBinTypeAndNum(const string& parameterString, vector<string>& vType, vector<int>& vNum)
{
    string name;
    map<string, string> param;
    SelfCreate::fillInNameParam(parameterString, name, param);

    map<string, string>::iterator iter;
    if ((iter = param.find("bin-type")) != param.end())
    {
        vType = StringUtils::Split(iter->second, ",");
    }
 
    vNum.clear();
    if ((iter = param.find("bin-count")) != param.end())
    {
        vector<string> numBins = StringUtils::Split(iter->second, ",");
        for (vector<string>::iterator ii = numBins.begin(); ii != numBins.end(); ++ii)
        {
            vNum.push_back(Convert::toInt(*ii));
        }
    }
}

void CovariateParams::paramValueToIntVector(map<string, string>& params, const string& paramName, vector<int>& vInt)
{
    map<string, string>::iterator iter;
    if ((iter = params.find(paramName)) != params.end())
    {
        vector<string> values = StringUtils::Split(iter->second, ",");
        for (vector<string>::iterator ii = values.begin(); ii != values.end(); ++ii)
        {
            vInt.push_back(Convert::toInt(*ii));
        }
    }
}

void CovariateParams::paramValueToDoubleVector(map<string, string>& params, const string& paramName, vector<double>& vDouble)
{
    map<string, string>::iterator iter;
    if ((iter = params.find(paramName)) != params.end())
    {
        vector<string> values = StringUtils::Split(iter->second, ",");
        for (vector<string>::iterator ii = values.begin(); ii != values.end(); ++ii)
        {
            vDouble.push_back(Convert::toDouble(*ii));
        }
    }
}

int CovariateParams::paramValueToInt(map<string, string>& params, const string& paramName, int defaultValue)
{
    map<string, string>::iterator iter;
    if ((iter = params.find(paramName)) != params.end())
    {
        vector<string> values = StringUtils::Split(iter->second, ",");
        if (values.size() == 0) { Err::errAbort(paramName + " sub-option error - one value expected, none found"); }
        if (values.size() > 1) { Err::errAbort(paramName + " sub-option error - only one value expected"); }
        return Convert::toInt(values[0]);
    }
    else { return defaultValue; }
}

double CovariateParams::paramValueToDouble(map<string, string>& params, const string& paramName, double defaultValue)
{
    map<string, string>::iterator iter;
    if ((iter = params.find(paramName)) != params.end())
    {
        vector<string> values = StringUtils::Split(iter->second, ",");
        if (values.size() == 0) { Err::errAbort(paramName + " sub-option error - one value expected, none found"); }
        if (values.size() > 1) { Err::errAbort(paramName + " sub-option error - only one value expected"); }
        return Convert::toDouble(values[0]);
    }
    else { return defaultValue; }
}

std::string CovariateParams::translateToParamString(const map<string, int>& covMap)
{
    vector<string> temp;
    temp.resize(covMap.size());
    for (map<string, int>::const_iterator ii = covMap.begin(); ii != covMap.end(); ++ii)
    {
        temp[ii->second] = ii->first;
    }

    std::string retStr;
    for (int iIndex = 0; iIndex < temp.size(); iIndex++)
    {
        retStr += ",";
        retStr += temp[iIndex];
    }
    return retStr.substr(1);
}

void CovariateParams::translateFromParamString(const string& str, map<string, int>& retMap)
{
    std::vector<std::string> vec;
    Convert::strToStrVec(str, ',', vec);
    retMap.clear();
    for (int ii = 0; ii < vec.size(); ii++)
    {
        retMap[vec[ii]] = ii;
    }
}

bool CovariateParams::isReservedCovariateName(const std::string& str)
{
    return (str == "fragment-adapter-type" ||
            str == "fragment-length" ||
            str == "fragment-gc" ||
            str == "probe-gc" ||
            str == "local-gc" ||
            str == "median-signal" ||
            str == "marker-class");
}

bool CovariateParams::anyFileBasedCovariates()
{
    // Are there any non-reserved covariates
    for (map<string, int>::iterator ii = m_allCovariateMap.begin(); ii != m_allCovariateMap.end(); ++ii)
    {
        if (isReservedCovariateName(ii->first) == false) {return true;}
    }
    return false;
}
