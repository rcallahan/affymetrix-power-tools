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
 * @file CNAnalysisMethod.cpp
 *
 * @brief This file contains the CNAnalysisMethod class members.
 */

//
#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNExperiment.h"
#include "copynumber/CNAnalysisMethodCovariateParams.h"
//
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "portability/affy-base-types.h"
#include "util/AffxStatistics.h"
#include "util/Fs.h"
//
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace std;

int CNAnalysisMethod::m_iInstanceCount = 0;
std::vector<affymetrix_calvin_parameter::ParameterNameValueType> CNAnalysisMethod::m_vCelFileParams;
std::vector<affymetrix_calvin_parameter::ParameterNameValueType> CNAnalysisMethod::m_vParams;

AffxString CNAnalysisMethod::getPrefix()
{
  return "affymetrix-algorithm-param-";
}
int CNAnalysisMethod::getSegmentType()
{
  return CNSegment::getSegmentType(getName());
}
AffxString CNAnalysisMethod::getSegmentTypeString()
{
  return CNSegment::getSegmentTypeString(getSegmentType());
}

CNSegmentArray& CNAnalysisMethod::getSegments()
{
  return m_vSegments;
}
CNSegmentArray CNAnalysisMethod::getSegments(int iType)
{
  CNSegmentArray v;
  for (int iIndex = 0; (iIndex < m_vSegments.getCount()); iIndex++) {
    CNSegment* p = m_vSegments.getAt(iIndex);
    if (p->getSegmentType() == iType) {
      v.push_back(p);
    }
  }
  return v;
}
int CNAnalysisMethod::getSegmentCount(int iType)
{
  return m_vSegments.getSegmentCount(iType);
}

void CNAnalysisMethod::setEngine(BaseEngine* p)
{
  m_pEngine = p;
  m_iXChromosome = m_pEngine->getOptInt("xChromosome");
  m_iYChromosome = m_pEngine->getOptInt("yChromosome");
}
BaseEngine* CNAnalysisMethod::getEngine()
{
  return m_pEngine;
}
CNExperiment* CNAnalysisMethod::getExperiment()
{
  return m_pobjExperiment;
}
CNExperimentArray* CNAnalysisMethod::getExperiments()
{
  return m_pvExperiments;
}
void CNAnalysisMethod::setProbes(CNProbeArray& vProbes)
{
  m_pvProbes = &vProbes;
}
CNProbeArray* CNAnalysisMethod::getProbes()
{
  return m_pvProbes;
}
std::vector<affymetrix_calvin_parameter::ParameterNameValueType>* CNAnalysisMethod::getCelFileParams()
{
  return &m_vCelFileParams;
}
std::vector<affymetrix_calvin_parameter::ParameterNameValueType>* CNAnalysisMethod::getParams()
{
  return &m_vParams;
}

void CNAnalysisMethod::setup(CNExperiment& objExperiment, CNProbeSetArray& vProbeSets, CNProbeArray* pvProbes)
{
  m_pobjExperiment = &objExperiment;
  m_pvProbeSets = &vProbeSets;
  m_pvProbes = pvProbes;
}

void CNAnalysisMethod::setup(CNExperimentArray& vExperiments, int experimentIndex, CNProbeSetArray& vProbeSets)
{
  m_pvExperiments = &vExperiments;
  m_pobjExperiment = vExperiments.getAt(experimentIndex);
  m_pvProbeSets = &vProbeSets;
}

void CNAnalysisMethod::setup(CNExperimentArray& vExperiments, CNProbeSetArray& vProbeSets)
{
  m_pvExperiments = &vExperiments;
  m_pvProbeSets = &vProbeSets;
}

void CNAnalysisMethod::isSetup()
{
  if ((m_pobjExperiment == NULL) && (m_pvProbeSets == NULL)) {
    Err::errAbort("CNAnalysisMethod " + getName() + " is not setup properly.");
  }
}

void CNAnalysisMethod::setDataSetOffset(unsigned int ui)
{
  m_uiDataSetOffset = ui;
}
unsigned int CNAnalysisMethod::getDataSetOffset()
{
  return m_uiDataSetOffset;
}

void CNAnalysisMethod::memory(const AffxString& str)
{
  /*
  static uint64_t init = 0;
  uint64_t free = 0;
  uint64_t total = 0;
  uint64_t swapAvailable = 0;
  uint64_t memAvailable = 0;
  Util::memInfo(free, total, swapAvailable, memAvailable, false);
  if (init == 0) {init = total - free;}
  int iMB = (int)((total - free - init) / 1048576.0);
  if (iMB < 0) {iMB = 0;}
  Verbose::out(1, str + ": Memory usage = " + ::getInt(iMB) + " MB");
  */
}

/**
 * Constructor
 */
CNAnalysisMethod::CNAnalysisMethod()
{
    m_iInstanceCount++;
    m_uiDataSetOffset = 0;
    m_pEngine = NULL;
    m_pobjExperiment = NULL;
    m_pvExperiments = NULL;
    m_pvProbeSets = NULL;
    m_iXChromosome = 24;
    m_iYChromosome = 25;
    //m_boundsFilled = false;
    m_pvProbes = NULL;

    m_iStep = 0;
    m_iWindow = 0;
    m_iPointCount = 0;
    m_dBandwidth = 0;
    m_dCutoff = 0;
    m_dCleanThreshold = 0;
    m_bSymmetry = false;
    m_boundsFilled=false;
}

/**
 * Destructor
 */
CNAnalysisMethod::~CNAnalysisMethod()
{
  m_iInstanceCount--;
  if (m_iInstanceCount == 0) {
    m_vCelFileParams.clear();
    m_vParams.clear();
  }
}

/**
 * Return SelfDoc option associated with a specified name.
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc options to search.
 * @param const std::string& - The specified name to search for.
 * @return SelfDoc::Opt* - The SelfDoc option pointer or NULL if not found.
 */
SelfDoc::Opt* CNAnalysisMethod::getSelfDocOpt(std::vector<SelfDoc::Opt>& opts, const std::string& strName)
{
  for (int iIndex = 0; (iIndex < opts.size()); iIndex++) {
    if (opts[iIndex].name == strName) {
      return &opts[iIndex];
    }
  }
  return NULL;
}

/**
 * Setup a bool parameter for this analysis method.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return bool - The value of the parameter as setup by this function.
 */
bool CNAnalysisMethod::setupBoolParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
  SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
  if (popt == NULL) {
    throw(Except("SelfDoc::Opt not found: " + strName));
  }
  bool b = AffxByteArray(popt->defaultVal).parsebool();
  fillInValue(b, std::string(strName), params, doc);
  std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
  affymetrix_calvin_parameter::ParameterNameValueType param;
  param.SetName(wstr);
  param.SetValueInt8(b);
  m_vParams.push_back(param);
  return b;
}

/**
 * Setup an int parameter for this analysis method, including min amd max validation.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return int - The value of the parameter as setup by this function.
 */
int CNAnalysisMethod::setupIntParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
  SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
  if (popt == NULL) {
    throw(Except("SelfDoc::Opt not found: " + strName));
  }
  int i = ::getInt(popt->defaultVal);
  fillInValue(i, std::string(strName), params, doc);
  if ((popt->minVal != "NA") && (popt->minVal != "") && (i < ::getInt(popt->minVal))) {
    throw(Except("SelfDoc::Opt " + strName + " below minimum value of " + popt->minVal));
  }
  if ((popt->maxVal != "NA") && (popt->maxVal != "") && (i > ::getInt(popt->maxVal))) {
    throw(Except("SelfDoc::Opt " + strName + " above maximum value of " + popt->maxVal));
  }
  std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
  affymetrix_calvin_parameter::ParameterNameValueType param;
  param.SetName(wstr);
  param.SetValueInt32(i);
  m_vParams.push_back(param);
  return i;
}

/**
 * Setup a float parameter for this analysis method, including min amd max validation.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return float - The value of the parameter as setup by this function.
 */
float CNAnalysisMethod::setupFloatParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
  SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
  if (popt == NULL) {
    throw(Except("SelfDoc::Opt not found: " + strName));
  }
  float f = (float)::getDouble(popt->defaultVal);
  fillInValue(f, std::string(strName), params, doc);
  if ((popt->minVal != "NA") && (popt->minVal != "") && (f < (float)::getDouble(popt->minVal))) {
    throw(Except("SelfDoc::Opt " + strName + " below minimum value of " + popt->minVal));
  }
  if ((popt->maxVal != "NA") && (popt->maxVal != "") && (f > (float)::getDouble(popt->maxVal))) {
    throw(Except("SelfDoc::Opt " + strName + " above maximum value of " + popt->maxVal));
  }
  std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
  affymetrix_calvin_parameter::ParameterNameValueType param;
  param.SetName(wstr);
  param.SetValueFloat(f);
  m_vParams.push_back(param);
  return f;
}

/**
 * Setup a double parameter for this analysis method, including min amd max validation.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return double - The value of the parameter as setup by this function.
 */
double CNAnalysisMethod::setupDoubleParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
  SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
  if (popt == NULL) {
    throw(Except("SelfDoc::Opt not found: " + strName));
  }
  double d = (float)::getDouble(popt->defaultVal);
  fillInValue(d, std::string(strName), params, doc);
  if ((popt->minVal != "NA") && (popt->minVal != "") && (d < ::getDouble(popt->minVal))) {
    throw(Except("SelfDoc::Opt " + strName + " below minimum value of " + popt->minVal));
  }
  if ((popt->maxVal != "NA") && (popt->maxVal != "") && (d > ::getDouble(popt->maxVal))) {
    throw(Except("SelfDoc::Opt " + strName + " above maximum value of " + popt->maxVal));
  }
  std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
  affymetrix_calvin_parameter::ParameterNameValueType param;
  param.SetName(wstr);
  param.SetValueFloat(d);
  m_vParams.push_back(param);
  return d;
}

/**
 * Setup a string parameter for this analysis method, including min amd max validation.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return std::string - The value of the parameter as setup by this function.
 */
std::string CNAnalysisMethod::setupStringParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
{
  SelfDoc::Opt* popt = getSelfDocOpt(opts, strName);
  if (popt == NULL) {
    throw(Except("SelfDoc::Opt not found: " + strName));
  }
  std::string str = popt->defaultVal;
  fillInValue(str, std::string(strName), params, doc);
  if ((popt->minVal != "") && (str < popt->minVal)) {
    throw(Except("SelfDoc::Opt " + strName + " below minimum value of " + popt->minVal));
  }
  if ((popt->maxVal != "") && (str > popt->maxVal)) {
    throw(Except("SelfDoc::Opt " + strName + " above maximum value of " + popt->maxVal));
  }
  std::wstring wstr = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(strPrefix + strName);
  affymetrix_calvin_parameter::ParameterNameValueType param;
  param.SetName(wstr);
  param.SetValueAscii(str);
  m_vParams.push_back(param);
  return str;
}


/**
 * Cache information about starts and stops on chromosomes.  As in STL,
 * the stops are one past the last position.
 */


void CNAnalysisMethod::fillChrBounds()
{
  fillChrBoundsImpl(m_pvProbeSets, m_chrBounds);
}


pair<int, int> CNAnalysisMethod::getChrBounds(const int chr)
{
  return getChrBounds(chr, getProbeSets());
}


/**
 * Return the bounds of chromosome.
 * @param const int - The chromosome
 * @return pair<int,int> the start and stop of the chromosome.
 */
pair<int, int> CNAnalysisMethod::getChrBounds(const int chr, CNProbeSetArray* pvProbeSets)
{
  // Make sure that the desired map is filled.
  this->fillChrBoundsImpl(pvProbeSets, m_chrBounds);

  // If the chromosome is not there, here is what you get.
  if (m_chrBounds.find(chr) == m_chrBounds.end()) return pair<int, int>(0, 0);

  // If the chromosoem is there, here is what you get.
  return m_chrBounds[chr];
}

/**
 * Returns a vector of chromosomes. Need to make sure that lvalue conforms
 * to the the assigment from this method call.
 * @return vector<int> the list of chromosomes as integers
 */
vector<int> CNAnalysisMethod::getChromosomes(CNProbeSetArray * pvProbeSetArray)
{
  // Make sure that the desired map is filled.
  this->fillChrBoundsImpl(pvProbeSetArray, m_chrBounds);

  // A vector with chromosomes as ints.
  vector<int>chr_vec;

  for (std::map<int, pair<int, int> >::iterator iter = m_chrBounds.begin();
       iter != m_chrBounds.end(); iter++) {
    chr_vec.push_back(iter->first);
  }

  // Return the vector of chromosomes.
  return chr_vec;
}


/**
 * Return the number of probesets for a specified chromosome.
 * @param int - The specified chromosome
 * @return int - The number of probe sets.
 */

int CNAnalysisMethod::getProbeSetCount(int chr, CNProbeSetArray* pvProbeSets)
{
  // Make sure that the desired map is filled.
  this->fillChrBoundsImpl(pvProbeSets, m_chrBounds);

  // If the chromosome is not there, it has size 0.  If this is an
  // error for some context then let the calling context handle the
  // problem.
  if (m_chrBounds.find(chr) == m_chrBounds.end()) return 0;

  // If the chromosoem is there, here is what you get.
  pair<int, int> chr_span = m_chrBounds[chr];
  return chr_span.second - chr_span.first;
}

void CNAnalysisMethod::bin( std::vector<float>& vValues,
                            std::vector<int>& vBinIndexes,
                            int iBinCount)
{
  vBinIndexes.resize(vValues.size());
  AffxMultiDimensionalArray<float> vTemp(vValues.size());
  for (int iIndex = 0; (iIndex < vValues.size()); iIndex++) {
    vTemp.set(iIndex, vValues[iIndex]);
  }
  vTemp.quickSort();
  AffxMultiDimensionalArray<float> v(iBinCount);
  float fPercentile = 0;
  for (int iBinIndex = 0; (iBinIndex < iBinCount); iBinIndex++) {
    fPercentile = vTemp.percentile((1.0 / (double)iBinCount) * (iBinIndex + 1), false);
    v.set(iBinIndex, fPercentile);
  }

  // Do binning.
  unsigned int uiLowIndex = 0;
  unsigned int uiHighIndex = 0;
  for (int iIndex = 0; (iIndex < vValues.size()); iIndex++) {
    vBinIndexes[iIndex] = -1;
    v.binarySearch(vValues[iIndex], uiLowIndex, uiHighIndex);
    vBinIndexes[iIndex] = uiHighIndex;
    if (vBinIndexes[iIndex] == -1) {
      throw(Except("Binning failed."));
    }
  }
}


/**
 * Adapted from AffxMultiDimensionalArray::percentile().
 * NB: dPercentile is forced <= 1.0 as STL operator[] is unchecked.
 */
float CNAnalysisMethod::getPercentile(const vector<float>& vec, double dPercentile)
{
    dPercentile = min(dPercentile, 1.0);
    double dIndex = (vec.size() - 1) * dPercentile;
    double dMultiplier = dIndex - floor(dIndex);
    int low = (int)floor(dIndex);
    int high = (int)ceil(dIndex);

    return vec[low] + dMultiplier*(vec[high] - vec[low]);
}

/**
 * Bin vValues[] into iBinCount bins containing equal number
 * of elements.
 * @param const vector<float> - the vector to be binned
 * @param vector<int> - bin assignments (0-based)
 */
void CNAnalysisMethod::binEqualNumber(
                                    const std::vector<float>& vValues,
                                    std::vector<int>& vBinIndexes,
                                    int iBinCount)
{
    vector<float> vTemp(vValues);
    std::sort(vTemp.begin(), vTemp.end());
    vector<float> v(iBinCount);
    for (int iBinIndex = 0; iBinIndex < iBinCount; iBinIndex++) {
        v[iBinIndex] = getPercentile(vTemp, (iBinIndex + 1)/(double)iBinCount);
    }
    // Do binning.
    for (int iIndex = 0; iIndex < vValues.size(); iIndex++) {
        vector<float>::iterator it = std::lower_bound(v.begin(), v.end(), vValues[iIndex]);
        vBinIndexes[iIndex] = it - v.begin();
    }

    //Verbose::out(4, "CNAnalysisMethod::binEqualNumber - vTemp range = [" + Convert::toString(vTemp[0]) + ", " + 
    //    Convert::toString(vTemp[vTemp.size()-1]) + "]");
    //for (int i =0; i < iBinCount; ++i)
    //{
    //    Verbose::out(4, "   v[" + Convert::toString(i) + "]= " + Convert::toString((double)v[i]));
    //}
}

/**
 * Bin vValues[] into iBinCount equally spaced bins.
 * @param const vector<float> - the vector to be binned
 * @param vector<int> - bin assignments (0-based)
 */
void CNAnalysisMethod::binEqualSpacing(
                                    const std::vector<float>& vValues,
                                    std::vector<int>& vBinIndexes,
                                    int iBinCount)
{
    vector<float> v(iBinCount);
    vector<float>::const_iterator itMin = std::min_element(vValues.begin(), vValues.end());
    vector<float>::const_iterator itMax = std::max_element(vValues.begin(), vValues.end());
    float delta = (*itMax - *itMin)/iBinCount;
    for (int iBinIndex = 0; iBinIndex < iBinCount; iBinIndex++) {
        v[iBinIndex] = *itMin + (iBinIndex + 1)*delta;
    }
    // Do binning.
    for (int iIndex = 0; iIndex < vValues.size(); iIndex++) {
        vector<float>::iterator it = std::lower_bound(v.begin(), v.end(), vValues[iIndex]);
        vBinIndexes[iIndex] = it - v.begin();
    }

    //Verbose::out(4, "CNAnalysisMethod::binEqualSpacing");
    //for (int i =0; i < iBinCount; ++i)
    //{
    //    Verbose::out(4, "   v[" + Convert::toString(i) + "]= " + Convert::toString((double)v[i]));
    //}
}

/** Bin vValues[] (presumed discrete) into bins containing equal values.
  * At the end the number of bins equals the number of distinct values.
  * @param const vector<float> - the vector to be binned
  * @param vector<int> - bin assignments (0-based)
  * @return int - number of created bins
  */
int CNAnalysisMethod::covariateIsBinAssignment(const std::vector<float>& vValues, std::vector<int>& vBinIndexes)
{
    vector<float> vTemp(vValues);
    std::sort(vTemp.begin(), vTemp.end());
    vTemp.erase(std::unique(vTemp.begin(), vTemp.end()), vTemp.end());

    // Do binning
    for (int iIndex = 0; iIndex < vValues.size(); iIndex++) {
        vector<float>::iterator it = std::lower_bound(vTemp.begin(), vTemp.end(), vValues[iIndex]);
        vBinIndexes[iIndex] = it - vTemp.begin();
    }
    return vTemp.size();
}

// This was moved here from CNAnalysisMethodAllelePeaks since it is also used in Chipstream for normal-diploid determination.
double CNAnalysisMethod::bwnrd(const vector<double> &x, double fac)
{
  int n = x.size();

  if (x.size() <= 2) {
    throw(Except("The allele-peaks analysis method has failed. Vector size to small in bwnrd()."));
  }

  // Calculate interquartile range into h
  double *tmpdata = new double [n];
  int prctile25 = (int)(.25 * n - .5);
  if (prctile25 < 0) {
    prctile25 = 0;
  }
  for (int i = 0; i < n; ++i) {
    tmpdata[i] = x[i];
  }
  double val25 = klowest_select(tmpdata, n, prctile25);
  int prctile75 = (int)(.75 * n - .5);
  if (prctile75 < 0) {
    prctile75 = 0;
  }
  for (int i = 0; i < n; ++i) {
    tmpdata[i] = x[i];
  }
  double val75 = klowest_select(tmpdata, n, prctile75);
  double h = (val75 - val25);

  delete [] tmpdata;

  // calculate variance of x into var
  double var = 0;
  double mean = 0;
  for (int i = 0; i < n; ++i) {
    mean += x[i];
  }
  mean = mean / n;
  for (int i = 0; i < n; ++i) {
    var += (x[i] - mean) * (x[i] - mean);
  }
  var = var / (n - 1);

  return (fac * 1.06 * min(sqrt(var), h / 1.34) * pow((double)n, -.2));
}

/**
 * kernel density() with an normal kernel
 * @param dat - data of interest.
 * @param weights (matches dat)
 * @param bandWidth
 *
 * @return - values in density evaluated at values xOut
 */
void CNAnalysisMethod::kdensity(vector<double> &dat,
                                vector<double> &density,
                                vector<double> &xOut,
                                vector<double> &weight,
                                double bandWidth)
{
  double from, to, maximum, minimum;
  // vector<float> xOut(numBins, 0.0);

  if ((dat.size() <= 0) || (density.size() <= 0)) {
    throw(Except("The allele-peaks analysis method has failed. Empty data or density vector in kdensity()."));
  }
  int n = dat.size();

  int m = density.size();

  // bandWidth = findBandWidth(dat);
  maximum = vectorMax(dat);
  minimum = vectorMin(dat);

  /* Set the limits. */
  from = minimum - 3 * bandWidth;
  to = maximum + 3 * bandWidth;
  double xDelta = (to - from) / (m - 1);

  for (int i = 0; i < m ; i++) {
    xOut[i] = from + (i * xDelta);
  }

  // get lazy and allocate a big fat vector, no error checking
  vector<double> z(n*m, 0); // n x m matrix
  int k = 0;
  for (int i = 0; i < n; ++i) { // i-th row
    for (int j = 0; j < m; ++j) {  // j-th row
      // apply the kernel
      z[k] = phi((xOut[j] - dat[i]) / bandWidth);
      k++;
    }
  }
  vector<int> rowIndex(n, 0);
  for (int i = 0; i < n; ++i) { // i-th column
    rowIndex[i] = i * m;
  }
  for (int i = 0; i < m; ++i) { // i-th column
    for (int j = 0; j < n; ++j) {  // j-th row
      density[i] += weight[j] * z[i+rowIndex[j]];
    }
    density[i] = density[i] / bandWidth;
  }
}

double CNAnalysisMethod::trapzoid(vector<double> &x, vector<double> &y)
{
  // x: x values
  // y: y = f(x) values
  // returns the integral of f

  int n = x.size();
  if ((y.size() != x.size()) || (n <= 1)) {
    throw(Except("The allele-peaks analysis method has failed. Unequal vector lengths in trapzoid()."));
  }

  vector<vector<double> > tmpdata(n);

  for (int i = 0; i < n; ++i) {
    vector<double> z(2);
    z[0] = x[i];
    z[1] = y[i];

    tmpdata[i] = z;
  }
  sort(tmpdata.begin(), tmpdata.end(), comparex());

  double integral = 0;
  for (int i = 0; i < (n - 1); ++i) {
    vector<double> z = tmpdata[i];
    vector<double> z1 = tmpdata[i+1];
    integral += (z1[1] + z[1]) * (z1[0] - z[0]) / 2;
  }
  return(integral);
}

/**
 * Find the minimum value in a vector of doubles, should be moved to RMA.cpp
 *
 * @param dat - vector of data.
 *
 * @return minimum value.
 */
double CNAnalysisMethod::vectorMin(vector <double> &dat)
{
  unsigned int i = 0;
  double minimum = -1;
  if (dat.size() <= 0) {
    throw(Except("The allele-peaks analysis method has failed. Data vector is empty in vectorMin()."));
  }
  minimum = dat[0];
  for (i = 0; i < dat.size(); i++) {
    minimum = min(minimum, dat[i]);
  }
  return minimum;
}

/**
 * Find the maximum value in a vector of doubles, should be moved to RMA.cpp.
 *
 * @param dat - vector of data.
 *
 * @return maximum value.
 */
double CNAnalysisMethod::vectorMax(vector <double> &dat)
{
  unsigned int i = 0;
  double maximum = -1;
  if (dat.size() <= 0) {
    throw(Except("The allele-peaks analysis method has failed. Data vector is empty in vectorMax()."));
  }
  maximum = dat[0];
  for (i = 0; i < dat.size(); i++) {
    maximum = max(maximum, dat[i]);
  }
  return maximum;
}


/**
 * Standard normal probability density function.
 * @param x - value of interest.
 * @return - density at value supplied
 */
//double RMA::phi(double x){
double CNAnalysisMethod::phi(double x)
{
  double pi = 3.14159265358979323846;
  return 1 / sqrt(2 * pi)* exp(-0.5 * x * x);
}


/*==================================================================================*/
/* Mean function applied to (running) window. All additions performed using         */
/* addition algorithm which tracks and corrects addition round-off errors (see      */
/*  http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps)*/
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/

/* SumErr - macro calculating error of the summing operation */
#define SumErr(a,b,ab) ((((a)>(b)) == ((a)>-(b))) ?  (b) - ((ab)-(a)) : (a) - ((ab)-(b)) )
/* SUM_1 - macro for calculating Sum+=x; Num+=n; Which is NaN aware and have minimal (single number) overflow error correction */
#define SUM_1(x,n, Sum, Err, Num)   if (x == x){ y=Sum; Err+=x; Sum+=Err; Num+=n; Err=SumErr(y,Err,Sum);  }

void CNAnalysisMethod::runmean(double *In, double *Out, const int *nIn, const int *nWin)
{
  /* medium size version with NaN's and edge calculation, but only one level of round-off correction*/
  int i, k1, k2, Num, n = *nIn, m = *nWin;
  double *in, y, *out, Err, Sum;
  double NaN = std::numeric_limits<double>::quiet_NaN();
  k2  = m >> 1;       /* right half of window size */
  k1  = m - k2 - 1;  /* left half of window size */
  in = In; out = Out;
  Sum = 0;           /* we need to calculate initial 'Sum' */
  Err = 0;
  Num = 0;
  /* step 1 - find mean of elements 0:(k2-1) */
  for (i = 0; i < k2; i++) {
    SUM_1(in[i], 1, Sum, Err, Num)
  }
  /* step 2 - left edge - start expanding the moving window to the right */
  for (i = k2; i < m; i++, out++) {
    SUM_1(in[i], 1, Sum, Err, Num)
    *out = (Num ? (Sum + Err) / Num : NaN);  /* save mean and move window */
  }
  /* step 3: runsum of the rest of the vector. Inside loop is same as:   */
  /* *out = *(out-1) - *in + *(in+m); but with round of error correction */
  for (i = m; i < n; i++, out++, in++) {
    SUM_1(in[m] ,  1, Sum, Err, Num)
    SUM_1(-(*in), -1, Sum, Err, Num)
    *out = (Num ? (Sum + Err) / Num : NaN);  /* save mean and move window */
  }
  /* step 4 - right edge - right side reached the end and left is shrinking  */
  for (i = 0; i < k2; i++, out++, in++) {
    SUM_1(-(*in), -1, Sum, Err, Num)
    *out = (Num ? (Sum + Err) / Num : NaN);  /* save mean and move window */
  }
}

double CNAnalysisMethod::corr(Matrix& mx1, Matrix& mx2)
{
  ColumnVector v1 = mx1.AsColumn();
  ColumnVector v2 = mx2.AsColumn();
  if (v1.Nrows() != v2.Nrows()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  int iLength = v1.Nrows();
  double dSum = 0;
  for (int iElementIndex = 0; (iElementIndex < iLength); iElementIndex++) {
    dSum += v1.element(iElementIndex);
  }
  double m1 = (dSum / (double)iLength);
  dSum = 0;
  for (int iElementIndex = 0; (iElementIndex < iLength); iElementIndex++) {
    dSum += v2.element(iElementIndex);
  }
  double m2 = (dSum / (double)iLength);
  double dNumerator = 0;
  double d1 = 0;;
  double d2 = 0;
  for (int iElementIndex = 0; (iElementIndex < iLength); iElementIndex++) {
    dNumerator += (v1.element(iElementIndex) - m1) * (v2.element(iElementIndex) - m2);
    d1 += (v1.element(iElementIndex) - m1) * (v1.element(iElementIndex) - m1) ;
    d2 += (v2.element(iElementIndex) - m2) * (v2.element(iElementIndex) - m2) ;
  }
  double dDenominator = sqrt(d1 * d2);
  return dNumerator / dDenominator;
}

double CNAnalysisMethod::norm(Matrix& mx)
{
  ColumnVector v = mx.AsColumn();
  int iLength = v.Nrows();
  double dSum = 0;
  for (int iElementIndex = 0; (iElementIndex < iLength); iElementIndex++) {
    dSum += v.element(iElementIndex) * v.element(iElementIndex);
  }
  return sqrt(dSum);
  /*
  DiagonalMatrix D;
  SVD(v, D);
  return D.Maximum();
  */
}

double CNAnalysisMethod::var(Matrix& mx)
{
  ColumnVector v = mx.AsColumn();
  int iLength = v.Nrows();
  double dSum = 0;
  for (int iElementIndex = 0; (iElementIndex < iLength); iElementIndex++) {
    dSum += v.element(iElementIndex);
  }
  double dMean = (dSum / (double)iLength);
  double dVariance = 0.0;
  if (iLength > 1) {
    double dSumOfSquares = 0.0;
    for (int iIndex = 0; (iIndex < iLength); iIndex++) {
      dSumOfSquares += ((v.element(iIndex) - dMean) * (v.element(iIndex) - dMean));
    }
    dVariance = dSumOfSquares / (iLength - 1);
  }
  return dVariance;
}

/**
 * Get the last autosome chromosome. As defined by being the last chromosome before the X chromosome value.
 * @return int - The chromosome number
 */
int CNAnalysisMethod::getLastAutosomeChromosome()
{
  int iLastAutosomeChromosome = 1;
  for (int iIndex = 0; (iIndex < getProbeSets()->getCount()); iIndex++) {
    CNProbeSet* pobjProbeSet = getProbeSets()->getAt(iIndex);
    if (pobjProbeSet->getChromosome() < m_iXChromosome) {
      iLastAutosomeChromosome = Max(iLastAutosomeChromosome, (int)pobjProbeSet->getChromosome());
    }
  }
  return iLastAutosomeChromosome;
}

#define ELEM_TYPE double

#define ELEM_SWAP(a,b) { register ELEM_TYPE t=(a);(a)=(b);(b)=t; }

ELEM_TYPE CNAnalysisMethod::klowest_select(ELEM_TYPE data[], int n, int klowest)
{

  int low, high ;
  int partition, ll, hh;

  low = 0 ; high = n - 1 ;
  for (;;) {
    if (high <= low) /* One element only */
      return data[klowest] ;

    if (high == low + 1) { /* Two elements only */
      if (data[low] > data[high])
        ELEM_SWAP(data[low], data[high]) ;
      return data[klowest] ;
    }

    /* Find the partition of low and high items; swap into position low */
    partition = (low + high) / 2;
    if (data[partition] > data[high]) {
      ELEM_SWAP(data[partition], data[high]) ;
    }
    if (data[low] > data[high]) {
      ELEM_SWAP(data[low], data[high]) ;
    }
    if (data[partition] > data[low]) {
      ELEM_SWAP(data[partition], data[low]) ;
    }

    /* Swap low item (now in position partition) into position (low+1) */
    ELEM_SWAP(data[partition], data[low+1]) ;

    /* Nibble from each end towards pivot, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do ll++; while (data[low] > data[ll]) ;
      do hh--; while (data[hh]  > data[low]) ;

      if (hh < ll)
        break;

      ELEM_SWAP(data[ll], data[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(data[low], data[hh]) ;

    /* Re-set active partition */
    if (hh <= klowest)
      low = ll;
    if (hh >= klowest)
      high = hh - 1;
  }
}

#undef ELEM_SWAP
#undef ELEM_TYPE

double CNAnalysisMethod::percentile(double prctile, double* data, int m)
{
  int tmp_len, nMeds;
  int incr, jump;
  tmp_len = m;
  nMeds = 1;
  incr = 1;
  jump = m;

  int klowest = (int)((prctile / 100) * tmp_len - .5);
  if (klowest < 0) {
    klowest = 0;
  }
  double dResult = 0;
  double* result = &dResult;

  double *tmpdata = new double [tmp_len];
  int indexdata = 0;
  int indexstart = 0;
  for (int j = 0; j < nMeds; ++j) {
    for (int i = 0; i < tmp_len; ++i) {
      tmpdata[i] = data[indexdata];
      indexdata = indexdata + incr;
    }
    indexstart = indexstart + jump;
    indexdata = indexstart;
    // klowest_select returns the k-th zero-based order statistic
    double *jresult = result + j;
    *jresult = klowest_select(tmpdata, tmp_len, klowest);

    // look for any data point with value larger than k-th order statistic
    double next = *jresult;
    int nlower = -1;
    int k = tmp_len;
    for (int i = 0; i < tmp_len; ++i) {
      if (tmpdata[i] > *jresult) {
        next = tmpdata[i];
        k = i;
        break;
      }
      if (tmpdata[i] <= *jresult) {
        nlower++;
      }
    }
    // find the data point larger than k-th order statistic
    for (int i = k; i < tmp_len; ++i) {
      if (tmpdata[i] <= *jresult) {
        nlower++;
      }
      if ((tmpdata[i] > *jresult) && (tmpdata[i] < next)) {
        next = tmpdata[i];
      }
    }
    double incr = next - *jresult;
    if (nlower > klowest) {
      incr = 0;
    }
    if (incr > 0) {
      double factor = (prctile / 100) * tmp_len - (klowest + .5);
      // mexPrintf("p=%f klowest=%d, low=%f, high=%f, incr=%f, factor=%f\n", prctile/100, klowest, *jresult, next, incr, factor);
      *jresult = *jresult + incr * factor;
    }
  }
  delete[] tmpdata;
  return dResult;
}

void CNAnalysisMethod::calculateSummaryLOH()
{
    int iChromosomeStartPosition = m_vSegments.getAt(0)->getStartPosition();
    int iChromosomeEndPosition = m_vSegments.getAt(0)->getEndPosition();
    int iSegmentStartPosition = iChromosomeStartPosition;
    int iSegmentEndPosition = iChromosomeEndPosition;

    float fLengthOfLoh=0;
    float fChromosomeLength=0;
    float fGenomeLength=0;
    float fSumLOH=0;
    float fLOH=0;
    int iNumberOfChromosomes=0;

    char cStartChromosome = m_vSegments.getAt(0)->getChromosome();
    char cPresentChromosome = cStartChromosome;

    if(m_vSegments.getAt(0)->getCall() == 1)
    {
        fLengthOfLoh = iSegmentEndPosition - iSegmentStartPosition;
    }  else
    {
        fLengthOfLoh = 0;
    }
    fChromosomeLength = iSegmentEndPosition - iSegmentStartPosition;

    bool autosomeGenomeLOHSet = false;

    for(int iIndex=1; iIndex < m_vSegments.getCount(); iIndex++){
        CNSegment* pobjSegment = m_vSegments.getAt(iIndex);
        cPresentChromosome = pobjSegment->getChromosome();

        if(cPresentChromosome != cStartChromosome){
            fLOH= fLengthOfLoh / fChromosomeLength;
            pair<char, float>  data(cStartChromosome, fLOH);
            getExperiment()->addChromosomeSummaryData(data);
            fGenomeLength += fChromosomeLength;
            fSumLOH += fLengthOfLoh;

            // Save genomeLOH accumulated so far
            // (i.e. accumulated over all autosomes)
            if (cPresentChromosome == m_iXChromosome)
            {
                getExperiment()->setAutosomeGenomeLOH( fSumLOH / fGenomeLength );
                autosomeGenomeLOHSet = true;
            }

            iNumberOfChromosomes++;

            iChromosomeStartPosition = pobjSegment->getStartPosition();
            iChromosomeEndPosition = pobjSegment->getEndPosition();
            iSegmentStartPosition = iChromosomeStartPosition;
            iSegmentEndPosition = iChromosomeEndPosition;

            if(pobjSegment->getCall() == 1)
            {
                fLengthOfLoh = iSegmentEndPosition - iSegmentStartPosition;
            }  else
            {
                fLengthOfLoh = 0;
            }
            fChromosomeLength = iSegmentEndPosition -iSegmentStartPosition;
            cStartChromosome = cPresentChromosome;

        } else
        {
            iSegmentStartPosition = pobjSegment->getStartPosition();
            iSegmentEndPosition = pobjSegment->getEndPosition();
            iChromosomeEndPosition = pobjSegment->getEndPosition();
            if(pobjSegment->getCall() == 1)
            {
                fLengthOfLoh += (iSegmentEndPosition - iSegmentStartPosition);
            }
            fChromosomeLength = iChromosomeEndPosition - iChromosomeStartPosition;
        }
    }
    fLOH= fLengthOfLoh / fChromosomeLength;
    pair<char, float>  data(cPresentChromosome, fLOH);
    getExperiment()->addChromosomeSummaryData(data);
    fGenomeLength += fChromosomeLength;
    fSumLOH += fLengthOfLoh;
    getExperiment()->setGenomeLOH(  fSumLOH / fGenomeLength );
    if (!autosomeGenomeLOHSet)
    {
        getExperiment()->setAutosomeGenomeLOH(getExperiment()->getGenomeLOH());
    }
    iNumberOfChromosomes++;
    getExperiment()->setNumberOfChromosomesToReport(iNumberOfChromosomes);
}

void CNAnalysisMethod::writeIntensities(std::string intensityFileInfix)
{

  affx::File5_File file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5 = NULL;

  std::string analysisString = "analysis";
  std::string strExperimentName = getExperiment()->getExperimentName();
  std::string fileName = Fs::join(getEngine()->getOpt("out-dir"),
                                  analysisString,
                                  strExperimentName + "." +intensityFileInfix + ".intensities.a5");

  Verbose::out(3, "Writing intensities for " + fileName);

  try{

    file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    group5 = file5.openGroup("Intensities", affx::FILE5_REPLACE);
    tsv5 = group5->openTsv("Intensities", affx::FILE5_REPLACE);
    tsv5->defineColumn(0,0,"ProbeSetID", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0,1,"Log2Intensity", affx::FILE5_DTYPE_FLOAT);

    CNProbe* pProbe;
    for(int iProbeIndex=0;  iProbeIndex < m_pvProbes->getCount(); iProbeIndex++)
    {
      pProbe = m_pvProbes->getAt(iProbeIndex);
      tsv5->set_i(0,0, pProbe->getProbeID());
      tsv5->set_f(0,1, (float) (log((double) pProbe->getIntensity()) / log(2.0)) );
      tsv5->writeLevel(0);
    }

    tsv5->close();
    delete tsv5;
    delete group5;
    file5.close();
  }
  catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump intensities to file."));}
}

void CNAnalysisMethod::writeSignals(std::string signalFileInfix)
{

  affx::File5_File file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5 = NULL;

  std::string analysisString = "analysis";
  std::string strExperimentName = getExperiment()->getExperimentName();
  std::string fileName = Fs::join(getEngine()->getOpt("out-dir"),analysisString,
                                  strExperimentName + "." +signalFileInfix + ".signals.a5");

  Verbose::out(3, "Writing signals for " + fileName);

  try{

    file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    group5 = file5.openGroup("Signals", affx::FILE5_REPLACE);
    tsv5 = group5->openTsv("Signals", affx::FILE5_REPLACE);
    tsv5->defineColumn(0,0,"ProbeSetName", affx::FILE5_DTYPE_STRING, 40);
    tsv5->defineColumn(0,1,"Chromosome", affx::FILE5_DTYPE_CHAR);
    tsv5->defineColumn(0,2,"Position", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0,3,"Log2AMedianIntensity", affx::FILE5_DTYPE_FLOAT);
    tsv5->defineColumn(0,4,"Log2BMedianIntensity", affx::FILE5_DTYPE_FLOAT);
    tsv5->defineColumn(0,5,"Log2AAlleleSignal", affx::FILE5_DTYPE_FLOAT);
    tsv5->defineColumn(0,6,"Log2BAlleleSignal", affx::FILE5_DTYPE_FLOAT);

    CNProbeSet* pProbeSet;
    for(int iProbeSetIndex=0;  iProbeSetIndex < m_pvProbeSets->getCount(); iProbeSetIndex++)
    {
        pProbeSet = m_pvProbeSets->getAt(iProbeSetIndex);
        {
            tsv5->set_string(0,0, pProbeSet->getProbeSetName());
            tsv5->set_c(0,1, pProbeSet->getChromosome());
            tsv5->set_i(0,2, pProbeSet->getPosition());
            tsv5->set_f(0,3, (float) (log((double) pProbeSet->getAMedianIntensity()) / log(2.0)) );
            tsv5->set_f(0,4, (float) (log((double) pProbeSet->getBMedianIntensity()) / log(2.0)) );
            tsv5->set_f(0,5, (float) (log((double) pProbeSet->getAAlleleSignal()) / log(2.0)) );
            tsv5->set_f(0,6, (float) (log((double) pProbeSet->getBAlleleSignal()) / log(2.0)) );
            tsv5->writeLevel(0);
        }
    }

    tsv5->close();
    delete tsv5;
    delete group5;
    file5.close();
  }
  catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump intensities to file."));}
}

void CNAnalysisMethod::writeMediansVectorAndMedian(std::string infix, const std::vector<float>& medians, float grandMedian)
{
  affx::File5_File file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5 = NULL;

  std::string analysisString = "analysis";
  std::string strExperimentName = getExperiment()->getExperimentName();
  std::string fileName = Fs::join(getEngine()->getOpt("out-dir"),
                                  analysisString,
                                  strExperimentName + "." + infix + ".medians.a5");

  Verbose::out(3, "Writing medians for " + fileName);

  try{

    file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    group5 = file5.openGroup("Medians", affx::FILE5_REPLACE);
    tsv5 = group5->openTsv("BinMedians", affx::FILE5_REPLACE);
    tsv5->defineColumn(0, 0, "BinMedians", affx::FILE5_DTYPE_FLOAT);

    for (int ii = 0; ii < medians.size(); ii++)
    {
        tsv5->set_f(0, 0, medians[ii]);
        tsv5->writeLevel(0);
    }
    tsv5->close();
    delete tsv5;

    tsv5 = group5->openTsv("GrandMedian", affx::FILE5_REPLACE);
    tsv5->defineColumn(0, 0, "GrandMedian", affx::FILE5_DTYPE_FLOAT);
    tsv5->set_f(0, 0, grandMedian);
    tsv5->writeLevel(0);

    tsv5->close();
    delete tsv5;
    delete group5;
    file5.close();
  }
  catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump medians to file."));}
}

void CNAnalysisMethod::writeMediansVector(std::string infix, const std::vector<float>& medians)
{
  affx::File5_File file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5 = NULL;

  // Self preservation
  // TODO: Fix this.
  if (getExperiment() == NULL) return;

  std::string analysisString = "analysis";
  std::string strExperimentName = getExperiment()->getExperimentName();
  std::string fileName = Fs::join(getEngine()->getOpt("out-dir"),
                                  analysisString,
                                  strExperimentName + "." + infix + ".medians.a5");

  Verbose::out(3, "Writing medians for " + fileName);

  try{

    file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    group5 = file5.openGroup("Medians", affx::FILE5_REPLACE);
    tsv5 = group5->openTsv("BinMedians", affx::FILE5_REPLACE);
    tsv5->defineColumn(0, 0, "BinMedians", affx::FILE5_DTYPE_FLOAT);

    for (int ii = 0; ii < medians.size(); ii++)
    {
        tsv5->set_f(0, 0, medians[ii]);
        tsv5->writeLevel(0);
    }
    tsv5->close();
    delete tsv5;
    delete group5;
    file5.close();
  }
  catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump medians to file."));}
}

void CNAnalysisMethod::writeSignalsTsv(std::string signalFileInfix, CNProbeSetArray* pvLocalProbeSets)
{
    std::string analysisString = "analysis";
    std::string strExperimentName = getExperiment()->getExperimentName();
    std::string fileName = Fs::join(getEngine()->getOpt("out-dir") , analysisString , strExperimentName)  + "." + signalFileInfix + ".signal.txt";

    affx::TsvFile *tsv = new affx::TsvFile;
    tsv->writeTsv(fileName);

    tsv->defineColumn(0,0,"ProbeSetName", affx::TSV_TYPE_STRING);
    tsv->defineColumn(0,1,"Chromosome", affx::TSV_TYPE_INT);
    tsv->defineColumn(0,2,"Position", affx::TSV_TYPE_INT);
    tsv->defineColumn(0,3,"Log2AMedianIntensity", affx::TSV_TYPE_FLOAT);
    tsv->defineColumn(0,4,"Log2BMedianIntensity", affx::TSV_TYPE_FLOAT);
    tsv->defineColumn(0,5,"Log2AAlleleSignal", affx::TSV_TYPE_FLOAT);
    tsv->defineColumn(0,6,"Log2BAlleleSignal", affx::TSV_TYPE_FLOAT);

    for (int iIndex = 0; (iIndex < pvLocalProbeSets->getCount()); iIndex++)
    {
        CNProbeSet* pProbeSet = getProbeSets()->getAt(iIndex);
        tsv->set(0,0,pProbeSet->getProbeSetName());
        tsv->set(0,1,pProbeSet->getChromosome());
        tsv->set(0,2,pProbeSet->getPosition());
        tsv->set(0,3,(float) (log((double) pProbeSet->getAMedianIntensity()) / log(2.0)));
        tsv->set(0,4,(float) (log((double) pProbeSet->getBMedianIntensity()) / log(2.0)));
        tsv->set(0,5,(float) (log((double) pProbeSet->getAAlleleSignal()) / log(2.0)));
        tsv->set(0,6,(float) (log((double) pProbeSet->getBAlleleSignal()) / log(2.0)));
        tsv->writeLevel(0);
    }
    delete tsv;
}


void CNAnalysisMethod::writeLog2RatiosTsv(std::string l2rFileInfix, CNProbeSetArray* pvLocalProbeSets)
{
    std::string analysisString = "analysis";
    std::string strExperimentName = getExperiment()->getExperimentName();
    std::string fileName = Fs::join(getEngine()->getOpt("out-dir"), analysisString , strExperimentName) + "." + l2rFileInfix + ".l2r.txt";

    affx::TsvFile *tsv = new affx::TsvFile;
    tsv->writeTsv(fileName);

    tsv->defineColumn(0,0,"ProbeSetName", affx::TSV_TYPE_STRING);
    tsv->defineColumn(0,1,"Chromosome", affx::TSV_TYPE_INT);
    tsv->defineColumn(0,2,"Position", affx::TSV_TYPE_INT);
    tsv->defineColumn(0,3,"Log2Ratio", affx::TSV_TYPE_INT);

    for (int iIndex = 0; (iIndex < pvLocalProbeSets->getCount()); iIndex++)
    {
        CNProbeSet* p = pvLocalProbeSets->getAt(iIndex);
        tsv->set(0,0,p->getProbeSetName());
        tsv->set(0,1,p->getChromosome());
        tsv->set(0,2,p->getPosition());
        tsv->set(0,3,p->getLog2Ratio());
        tsv->writeLevel(0);
    }
    delete tsv;
}


void CNAnalysisMethod::writeLog2Ratios(std::string l2rFileInfix)
{
  affx::File5_File file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5 = NULL;

  // Self preservation
  // TODO: Fix this.
  if (getExperiment() == NULL) return;

  std::string analysisString = "analysis";
  std::string strExperimentName = getExperiment()->getExperimentName();
  std::string fileName = Fs::join(getEngine()->getOpt("out-dir"),
                                  analysisString,
                                  strExperimentName + "." + l2rFileInfix + ".l2r.a5");
  
  std::string dirName = Fs::join(getEngine()->getOpt("out-dir"), analysisString);  


  Verbose::out(3, "Writing Log2Ratios for " + fileName);

  try{
    if(!Fs::dirExists(dirName))
    {
        Fs::mkdir(dirName);
    }

    file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    group5 = file5.openGroup("Signals", affx::FILE5_REPLACE);
    tsv5 = group5->openTsv("Signals", affx::FILE5_REPLACE);
    tsv5->defineColumn(0,0,"ProbeSetName", affx::FILE5_DTYPE_STRING,40);
    tsv5->defineColumn(0,1,"Chromosome", affx::FILE5_DTYPE_CHAR);
    tsv5->defineColumn(0,2,"Position", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0,3,"Log2Ratio", affx::FILE5_DTYPE_FLOAT);

    CNProbeSet* pProbeSet;
    for(int iProbeSetIndex=0;  iProbeSetIndex < m_pvProbeSets->getCount(); iProbeSetIndex++)
    {
      pProbeSet = m_pvProbeSets->getAt(iProbeSetIndex);
      tsv5->set_string(0,0, pProbeSet->getProbeSetName());
      tsv5->set_c(0,1, pProbeSet->getChromosome());
      tsv5->set_i(0,2, pProbeSet->getPosition());
      tsv5->set_f(0,3, pProbeSet->getLog2Ratio());
      tsv5->writeLevel(0);
    }

    tsv5->close();
    delete tsv5;
    delete group5;
    file5.close();
  }
  catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump log2 ratios to file."));}
}


void CNAnalysisMethod::writeSmoothedLog2Ratios(std::string l2rFileInfix, CNProbeSetArray* vProbeSets)
{
  affx::File5_File file5;
  affx::File5_Group* group5;
  affx::File5_Tsv* tsv5 = NULL;

  // Self preservation
  // TODO: Fix this.
  if (getExperiment() == NULL) return;

  std::string analysisString = "analysis";
  std::string strExperimentName = getExperiment()->getExperimentName();
  std::string fileName = Fs::join(getEngine()->getOpt("out-dir"),
                                  analysisString,
                                  strExperimentName + "." + l2rFileInfix + ".l2r.a5");
  
  std::string dirName = Fs::join(getEngine()->getOpt("out-dir"), analysisString);  


  Verbose::out(3, "Writing smoothed Log2Ratios for " + fileName);

  try{
    if(!Fs::dirExists(dirName))
    {
        Fs::mkdir(dirName);
    }

    file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
    group5 = file5.openGroup("Signals", affx::FILE5_REPLACE);
    tsv5 = group5->openTsv("Signals", affx::FILE5_REPLACE);
    tsv5->defineColumn(0,0,"ProbeSetName", affx::FILE5_DTYPE_STRING,40);
    tsv5->defineColumn(0,1,"Chromosome", affx::FILE5_DTYPE_CHAR);
    tsv5->defineColumn(0,2,"Position", affx::FILE5_DTYPE_INT);
    tsv5->defineColumn(0,3,"SmoothedLog2Ratio", affx::FILE5_DTYPE_FLOAT);

    CNProbeSet* pProbeSet;
    for(int iProbeSetIndex=0;  iProbeSetIndex < vProbeSets->getCount(); iProbeSetIndex++)
    {
      pProbeSet = vProbeSets->getAt(iProbeSetIndex);
      tsv5->set_string(0,0, pProbeSet->getProbeSetName());
      tsv5->set_c(0,1, pProbeSet->getChromosome());
      tsv5->set_i(0,2, pProbeSet->getPosition());
      tsv5->set_f(0,3, pProbeSet->getSmoothedLog2Ratio());
      tsv5->writeLevel(0);
    }

    tsv5->close();
    delete tsv5;
    delete group5;
    file5.close();
  }
  catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump smoothed log2 ratios to file."));}
}


void CNAnalysisMethod::computeSCAR(std::vector<float>& vSCAR)
{
    float fSCAR;

    for (int iIndex = 0; (iIndex < (int)getProbeSets()->size()); iIndex++)
    {
        CNProbeSet* pobjProbeSet = getProbeSets()->at(iIndex);
        pobjProbeSet->setLoh(numeric_limits<float>::quiet_NaN());
        fSCAR = numeric_limits<float>::quiet_NaN();

        // We attempt to calculate SCAR values for all ProbeSets designated as SNP.
        if(pobjProbeSet->processAsSNP())
        {
            if(pobjProbeSet->getAAlleleSignal() == 0.0 || pobjProbeSet->getBAlleleSignal() == 0.0){
                Verbose::warn(1, "CNAnalysisMethodLOHCyto2:: Invalid A or B allele signals were found for the ProbeSet " + pobjProbeSet->getProbeSetName());
                vSCAR[iIndex] = CN_INVALID_DOUBLE;
                pobjProbeSet->setValidSCARExists(false);
                continue;
            }
            if(pobjProbeSet->getMuAA() == CN_INVALID_DOUBLE || pobjProbeSet->getMuAB() == CN_INVALID_DOUBLE || pobjProbeSet->getMuBB() == CN_INVALID_DOUBLE)
            {
                Verbose::warn(1, "CNAnalysisMethodLOHCyto2:: Invalid snp-reference parameters for probe set " + pobjProbeSet->getProbeSetName() + ". LOH computations will be done without this probe set.");
                vSCAR[iIndex] = CN_INVALID_DOUBLE;
                pobjProbeSet->setValidSCARExists(false);
                continue;
            }

            fSCAR = (
                      2.0*( log2((float)pobjProbeSet->getAAlleleSignal()/(float)pobjProbeSet->getBAlleleSignal()) - pobjProbeSet->getMuAB() )
                    )
                    /
                    (
                      fabs( pobjProbeSet->getMuAA() - pobjProbeSet->getMuAB() )
                      +
                      fabs( pobjProbeSet->getMuBB() - pobjProbeSet->getMuAB() )
                    );
            if(fSCAR < -4.0) fSCAR = -4.0;
            if(fSCAR > 4.0) fSCAR = 4.0;


            vSCAR[iIndex] = fSCAR;
            pobjProbeSet->setValidSCARExists(true);
        } else
        {
            vSCAR[iIndex] = CN_INVALID_DOUBLE;
            pobjProbeSet->setValidSCARExists(false);
        }
        pobjProbeSet->setSCAR(fSCAR);
    }
}

void CNAnalysisMethod::calculateWindows( vector<pair<int, int> >& windowBds,
                                         vector<pair<int, int> >& stepWindowBds,
                                         const pair<int, int>& chrBound,
                                         int iStep,
                                         int iWindow
                                         )
{
    int chrStart = chrBound.first;
    int chrEnd   = chrBound.second - 1;
    int N        = chrEnd - chrStart + 1;
    float fk     = float(N - iWindow)/iStep;
    int numWin   = max(1, int(std::floor(fk)) + 1);    // must have at least one window
    int delta    = N - iWindow - iStep*(numWin - 1);
    int eps      = delta/2;
    int mm       = (iWindow - iStep)/2;
    int nn       = iWindow - mm - iStep;

    int winSize = iWindow - 1;
    if (fk < 1.0 || fk >= 2.0) {
        for (int j = 0; j < numWin; j++) {
            int left = chrStart + j*iStep + eps;
            windowBds.push_back(make_pair(left, left + winSize));
            stepWindowBds.push_back(make_pair(left + mm, left + winSize - nn));
        }
        // Extend first and last windows to cover entire chromosome
        windowBds[0].first = chrStart;
        windowBds[numWin - 1].second = chrEnd;
        stepWindowBds[0].first = chrStart;
        stepWindowBds[numWin - 1].second = chrEnd;
    }
    else {
        // Set up 2 equal windows even if numWin == 1
        windowBds.push_back(make_pair(chrStart, chrStart + winSize + eps));
        windowBds.push_back(make_pair(chrEnd - winSize - eps, chrEnd));
        stepWindowBds.push_back(make_pair(chrStart, chrStart + N/2 - 1));
        stepWindowBds.push_back(make_pair(chrStart + N/2, chrEnd));
    }
}

void CNAnalysisMethod::calculateWeights( CNProbeSetArray* probeSets,
                                         vector<double>& vValsToProcess,
                                         vector<double>& weights,
                                         const pair<int, int>& windowBds,
                                         bool symmetry
                                         )
{
    const bool isCyto2 = m_pEngine->getOptBool("cyto2");
    double dSumReplicateCounts = 0;
    int sizeMax = windowBds.second - windowBds.first + 1;
    if (symmetry) {
        // Twice the size because we are adding mirror images.
        vValsToProcess.resize(2*sizeMax);
        weights.resize(2*sizeMax);

        int k = 0;
        for (int i = windowBds.first; i <= windowBds.second; i++) {
            dSumReplicateCounts += probeSets->getAt(i)->getReplicateCount();

            if (isCyto2) {
                vValsToProcess[k] = probeSets->getAt(i)->getSCAR();
            }
            else {
                vValsToProcess[k] = probeSets->getAt(i)->getAllelicDifference();
            }
            vValsToProcess[k+sizeMax] = -vValsToProcess[k];   // the mirror image
            k++;
        }
        k = 0;
        for (int i = windowBds.first; i <= windowBds.second; i++) {
            if (dSumReplicateCounts == 0) {
                weights[k] = 1.0;
            }
            else {
                weights[k] = probeSets->getAt(i)->getReplicateCount()/dSumReplicateCounts;
            }
            weights[k+sizeMax] = weights[k];
            k++;
        }
    }
    else {
        vValsToProcess.resize(sizeMax);
        weights.resize(sizeMax);

        int k = 0;
        for (int i = windowBds.first; i <= windowBds.second; i++) {
            dSumReplicateCounts += probeSets->getAt(i)->getReplicateCount();

            if (isCyto2) {
                vValsToProcess[k] = probeSets->getAt(i)->getSCAR();
            }
            else {
                vValsToProcess[k] = probeSets->getAt(i)->getAllelicDifference();
            }
            k++;
        }
        k = 0;
        for (int i = windowBds.first; i <= windowBds.second; i++) {
            if (dSumReplicateCounts == 0) {
                weights[k] = 1.0;
            }
            else {
                weights[k] = probeSets->getAt(i)->getReplicateCount()/dSumReplicateCounts;
            }
            k++;
        }
    }
}

int CNAnalysisMethod::findpeaks(    vector<int> &valleys,
                                    vector<int> &peaks,
                                    vector<double> &y,
                                    double delta,
                                    vector<double> &x
                                    )
{

    // findpeaks: Detect peaks and valleys in y = f(x), x is ordered low to high
    //       findpeaks(valleys, peaks, y, delta, x) the indices of x corresponding
    //         to valleys is returned in valleys, the indices of x corresponding to
    //       to peaks is returned in peaks.
    //       A point is considered a maximum peak if it has the maximal
    //       value, and was preceded (to the left) by less than delta

    if ((y.size() <= 3) || (x.size() != y.size()) || (delta <= 0))
    {
        throw(Except("The allele-peaks analysis method has failed. Inconsistent input data found in findpeaks()."));
    }

    double mn = numeric_limits<double>::infinity();
    double mx = -numeric_limits<double>::infinity();

        bool lookformax = true;

    int imax = 0;
    int imin = 0;

    int npeak = 0;

    for (unsigned int i=0; i<y.size(); ++i){
        double v = y[i];
        if (v>mx) {
            mx = v;
            imax = i;
        }
        if (v<mn) {
            mn = v;
            imin = i;
        }
        if (lookformax) {
            if (v < (mx - delta)){
                peaks.push_back(imax);
                npeak++;
                mn = v;
                imin = i;
                lookformax = false;
            }
        } else {
            if (v > (mn + delta)){
                valleys.push_back(imin);
                mx = v;
                imax = i;
                lookformax = true;
            }
        }
    }
    return(npeak);
}


/**
  * For each probe set store maximum number of peaks
  * over all windows containing this probe set.
  */
void CNAnalysisMethod::findMaxPeaks( CNProbeSetArray* probeSets,
                                     const vector<pair<int, int> >& windowBds,
                                     const vector<vector<double> >& peaks
                                     )
{
    int size = windowBds.size();
    for (int i = windowBds[0].first; i <= windowBds[size - 1].second; i++) {
        probeSets->getAt(i)->setMaxPeaks(0);
    }
    for (int i = windowBds[0].first; i <= windowBds[0].second; i++) {
        probeSets->getAt(i)->setMaxPeaks(peaks[0].size());
    }
    for (int win = 1; win < windowBds.size(); win++) {
        int curNumPeaks = peaks[win].size();
        for (int i = windowBds[win].first; i <= windowBds[win].second; i++) {
            if (probeSets->getAt(i)->getMaxPeaks() < curNumPeaks) {
                probeSets->getAt(i)->setMaxPeaks(curNumPeaks);
            }
        }
    }
}


int CNAnalysisMethod::findcutoffInd(const vector<float>& snpqc)
{
    if (snpqc.empty()) {
        Err::errAbort("CNAnalysisMethodAllelePeaks: empty vector snpqc in findcutoffInd()");
    }
    float SNPQC = getExperiment()->getSNPQC();

    vector<float>::const_iterator iter = std::upper_bound(snpqc.begin(), snpqc.end(), SNPQC);

    if (iter == snpqc.begin()) {
        iter++;
    }
    return iter - snpqc.begin() - 1;
}

void CNAnalysisMethod::shrinkToPeaks(CNProbeSetArray* probeSets, CNAnalysisMethod::PeakShrinkOverride *paramOverride)
{
    const bool keepIntermediateData = m_pEngine->getOptBool("keep-intermediate-data");
    const bool isCytoScanHD = m_pEngine->getOptBool("cytoscan-hd");
    const bool isCyto2 = m_pEngine->getOptBool("cyto2");
    const bool useKdensity = m_pEngine->getOptBool("use-old-kdensity-function");
    const int noCovariateIndex = -1;

    // set the parameters to the member values (the default)
    int iStep = m_iStep;
    int iWindow = m_iWindow;
    int iPointCount = m_iPointCount;
    double dBandwidth = m_dBandwidth;
    double dCutoff = m_dCutoff;
    double dCleanThreshold = m_dCleanThreshold;
    bool bSymmetry = m_bSymmetry;

    int iCovariateIndex = noCovariateIndex;   // used for intermediate data output

    bool saveAllelePeaksFlag = false;

    // reset to the override if supplied
    if (paramOverride)
    {
        iStep           = paramOverride->m_iStep_override;
        iWindow         = paramOverride->m_iWindow_override;
        iPointCount     = paramOverride->m_iPointCount_override;
        dBandwidth      = paramOverride->m_dBandwidth_override;
        dCutoff         = paramOverride->m_dCutoff_override;
        dCleanThreshold = paramOverride->m_dCleanThreshold_override;
        bSymmetry       = paramOverride->m_bSymmetry_override;
        iCovariateIndex = paramOverride->m_iCovariateIndex;
        saveAllelePeaksFlag = true;
        m_vCoarseAllelePeaks.resize(probeSets->size());
    }

    // For intermediate testing
    vector<pair<int, int> > peakWindows;
    vector<pair<int, int> > stepWindows;
    vector<pair<int, double> > closestPeak;
    vector<pair<pair<int, int>, vector<double> > > peakSets;
    vector<pair<int, double> > processedValues;
    apParameters apParametersValues;
    cutoffShrinkage cutoffShrinkageValues;
    //

    vector<double> vValsToProcess;
    vector<double> weights;
    vector<double> density(iPointCount);
    vector<double> xi(iPointCount);
    vector<double> vPeakLocations(iPointCount);
    vector<int> chromosomes;

    int cutoffIndFLD3 = findcutoffInd(m_vThreePeakFLD_X);
    int cutoffIndFLD4 = findcutoffInd(m_vFourPeakFLD_X);
    int cutoffIndShrink3 = findcutoffInd(m_vThreePeakShrink_X);
    int cutoffIndShrink4 = findcutoffInd(m_vFourPeakShrink_X);

    // For deep regression
    cutoffShrinkageValues.FLD3threshold = m_vThreePeakFLD_Y[cutoffIndFLD3];
    cutoffShrinkageValues.FLD4threshold = m_vFourPeakFLD_Y[cutoffIndFLD4];
    cutoffShrinkageValues.shrinkage3 = m_vThreePeakShrink_Y[cutoffIndShrink3];
    cutoffShrinkageValues.shrinkage4 = m_vFourPeakShrink_Y[cutoffIndShrink4];
    cutoffShrinkageValues.SNPQC = getExperiment()->getSNPQC();

    apParametersValues.window = iWindow;
    apParametersValues.step = iStep;
    apParametersValues.densityPointCount = iPointCount;  // presumed >= CNAnalysisMethod::allelePeakCount
    apParametersValues.peaksCutoff = dCutoff;
    apParametersValues.cleanThreshold = dCleanThreshold;
    apParametersValues.bandwidthFactor = dBandwidth;
    apParametersValues.symmetry = bSymmetry;
    apParametersValues.covariateIndex = iCovariateIndex;
    //

    vector<bool> filteredProbeSets(probeSets->getCount(), false);

    // Collect the chromosomes and sort them
    for (map<int, pair<int, int> >::iterator it = m_chrBounds.begin();
         it != m_chrBounds.end();
         ++it)
    {
        chromosomes.push_back(it->first);
    }
    sort(chromosomes.begin(), chromosomes.end());

    const bool useOldDensity = isCyto2 || useKdensity;
    for (int chrInd = 0; chrInd < chromosomes.size(); chrInd++) {
        vector<pair<int, int> > windowBds;
        vector<pair<int, int> > stepWindowBds;
        calculateWindows(windowBds, stepWindowBds, m_chrBounds[chromosomes[chrInd]], iStep, iWindow);

        vector<vector<double> > peaks(windowBds.size());

        for (int winInd = 0; winInd < windowBds.size(); winInd++)
        {
            bool symmetrySetting;
            if (useOldDensity) {
                symmetrySetting = bSymmetry;
                for (int z = 0; z < m_iPointCount; z++) {
                    density[z] = 0;
                }
            }
            else {
                // don't symmetrise by brute force, it's done more efficiently in fitDensityCurve()
                symmetrySetting = false;
            }
            calculateWeights(probeSets, vValsToProcess, weights, windowBds[winInd], symmetrySetting);

            double wsum = 0;
            for (int k=0; k<weights.size(); ++k) {
                wsum += weights[k];
            }
            for (int k=0; k<weights.size(); ++k) {
                weights[k] = weights[k]/wsum;
            }

            // call bandwidth calculation
            double bw = bwnrd(vValsToProcess, dBandwidth);

            if (useOldDensity) {
                kdensity(vValsToProcess, density, xi, weights, bw);
            }
            else {
                fitDensityCurve(vValsToProcess, weights, xi, density, iPointCount, bw, bSymmetry);
            }

            // compute overall signal
            vector<double> sqdensity(density.size());
            for (unsigned int k=0; k<density.size(); ++k) {
                sqdensity[k] = density[k] * density[k];
            }
            double overallSignal = trapzoid(xi, sqdensity);

            // find the peaks in (xi, density) and return in peaks
            vector<int> peakIndex;
            vector<int> valleyIndex;
            int npeak1 = findpeaks(valleyIndex, peakIndex, density, 0.046*overallSignal, xi);

            int nData = vValsToProcess.size();
            int npeak = 0;

            // trim off any peaks in extremes of data
            if (npeak1 > 0) {
                // calculate lower percentile
                double *tmpdata = new double[nData];
                int nLow = (int)(dCutoff*nData - 0.5);
                if (nLow < 0){
                    nLow = 0;
                }
                for (int i=0; i<nData; ++i){
                    tmpdata[i] = vValsToProcess[i];
                }
                double percentileLower = klowest_select(tmpdata, nData, nLow);

                // calculate upper percentile
                int nHigh = (int)((1-dCutoff)*nData - 0.5);
                if (nHigh < 0){
                    nHigh = 0;
                }
                for (int i=0; i<nData; ++i){
                    tmpdata[i] = vValsToProcess[i];
                }
                double percentileUpper = klowest_select(tmpdata, nData, nHigh);
                delete [] tmpdata;

                // ignore any peaks outside this range
                for (unsigned int i=0; i<peakIndex.size(); ++i){
                    if ((xi[peakIndex[i]]>=percentileLower) && (xi[peakIndex[i]]<=percentileUpper)){
                        vPeakLocations[npeak++] = xi[peakIndex[i]];
                    }
                }
            }
            int iPeakCount = npeak;
            int maxAllelePeakCount = allelePeakCount;
            if (isCytoScanHD) {
                maxAllelePeakCount = CytoScanHDAllelePeakCount;
            }

            if (iPeakCount > maxAllelePeakCount) {
                iPeakCount = maxAllelePeakCount;
            }
            for (int iAllelePeakIndex = 0; (iAllelePeakIndex < iPeakCount); iAllelePeakIndex++)
            {
                if (vPeakLocations[iAllelePeakIndex] < -iMaxPeak) {
                    vPeakLocations[iAllelePeakIndex] = -iMaxPeak;
                }
                if (vPeakLocations[iAllelePeakIndex] > iMaxPeak) {
                    vPeakLocations[iAllelePeakIndex] = iMaxPeak;
                }
            }

            for (int i = 0; i < iPeakCount; i++) {
                peaks[winInd].push_back(vPeakLocations[i]);
            }
        }
        findMaxPeaks(probeSets, windowBds, peaks);

        int chrStart = m_chrBounds[chromosomes[chrInd]].first;
        int chrEnd   = m_chrBounds[chromosomes[chrInd]].second - 1;
        float cutoff3 = m_vThreePeakFLD_Y[cutoffIndFLD3];
        float cutoff4 = m_vFourPeakFLD_Y[cutoffIndFLD4];

        filterCutOff(probeSets, chrStart, chrEnd, cutoff3, cutoff4, filteredProbeSets);
        filterNoMansLand(probeSets, filteredProbeSets, stepWindowBds, peaks);
        filterShrinkTowardPeaks(
                        probeSets,
                        filteredProbeSets,
                        stepWindowBds,
                        peaks,
                        m_vThreePeakShrink_Y[cutoffIndShrink3],
                        m_vFourPeakShrink_Y[cutoffIndShrink4],
                        closestPeak,
                        processedValues,
                        saveAllelePeaksFlag
                        );

        if (keepIntermediateData) {
            for (int winInd = 0; winInd < windowBds.size(); winInd++) {
                peakWindows.push_back(make_pair(windowBds[winInd].first, windowBds[winInd].second));
                stepWindows.push_back(make_pair(stepWindowBds[winInd].first, stepWindowBds[winInd].second));

                pair<int, int> window = make_pair(windowBds[winInd].first, windowBds[winInd].second);
                peakSets.push_back(make_pair(window, peaks[winInd]));
            }
        }
    }

    if (keepIntermediateData) {
                        writeIntermediateData(
                                            probeSets,
                                            filteredProbeSets,
                                            peakWindows,
                                            stepWindows,
                                            closestPeak,
                                            peakSets,
                                            processedValues,
                                            apParametersValues,
                                            cutoffShrinkageValues
                        );
    }
}

void CNAnalysisMethod::writeIntermediateData(
                                        CNProbeSetArray* probeSets,
                                        const vector<bool>& filteredProbeSets,
                                        const vector<pair<int, int> >& peakWindows,
                                        const vector<pair<int, int> >& stepWindows,
                                        const vector<pair<int, double> >& closestPeak,
                                        const vector<pair<pair<int, int>, vector<double> > >& peakSets,
                                        const vector<pair<int, double> >& processedValues,
                                        const apParameters& apParametersValues,
                                        const cutoffShrinkage& cutoffShrinkageValues
                                        )
{
    const string covarPrefix = ".crude";   // for "crude allele peaks"

    // insert covariate index into the file name if the index is valid
    const bool useCovariateInFileName = apParametersValues.covariateIndex >= 0;

    affx::File5_File file5;
    affx::File5_Group* group5;
    affx::File5_Tsv* tsv5 = NULL;

    std::string analysisString = "analysis";
    std::string strExperimentName = getExperiment()->getExperimentName();
    std::string fileName;
    if (useCovariateInFileName)
    {
        std::string covariateStr = ".covar" + Convert::toString(apParametersValues.covariateIndex);
        fileName = Fs::join(getEngine()->getOpt("out-dir"), analysisString,
                                        strExperimentName + covariateStr + covarPrefix + ".allelepeaks.a5");
    }
    else {
        fileName = Fs::join(getEngine()->getOpt("out-dir"), analysisString,
                                        strExperimentName + ".allelepeaks.a5");
    }

    Verbose::out(3, "Writing intermediate data for " + fileName);

    try {
        file5.open(fileName, affx::FILE5_CREATE | affx::FILE5_REPLACE);
        group5 = file5.openGroup("AllelePeaks", affx::FILE5_REPLACE);

        tsv5 = group5->openTsv("AllelePeaks", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "ProbeSetName", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT);
        if (useCovariateInFileName)
        {
            tsv5->defineColumn(0, 3, "AllelePeak", affx::FILE5_DTYPE_DOUBLE);
        }
        else {
            tsv5->defineColumn(0, 3, "AllelePeaks0", affx::FILE5_DTYPE_INT);
            tsv5->defineColumn(0, 4, "AllelePeaks1", affx::FILE5_DTYPE_INT);
        }

        for (int i = 0; i < probeSets->getCount(); i++) {
            if (filteredProbeSets[i]) {
                tsv5->set_string(0, 0, probeSets->getAt(i)->getProbeSetName());
                tsv5->set_i(0, 1, probeSets->getAt(i)->getChromosome());
                tsv5->set_i(0, 2, probeSets->getAt(i)->getPosition());
                if (useCovariateInFileName)
                {
                    tsv5->set_d(0, 3, m_vCoarseAllelePeaks[i]);
                }
                else {
                    tsv5->set_i(0, 3, probeSets->getAt(i)->getAllelePeaks1());
                    tsv5->set_i(0, 4, probeSets->getAt(i)->getAllelePeaks2());
                }
                tsv5->writeLevel(0);
            }
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("PeakWindows", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "StartProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "StartChromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "StartPosition", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "EndProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 4, "EndChromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 5, "EndPosition", affx::FILE5_DTYPE_INT);

        for (int i = 0; i < peakWindows.size(); i++) {
            int startInd = peakWindows[i].first;
            int endInd   = peakWindows[i].second;
            tsv5->set_string(0, 0, probeSets->getAt(startInd)->getProbeSetName());
            tsv5->set_i(0, 1, probeSets->getAt(startInd)->getChromosome());
            tsv5->set_i(0, 2, probeSets->getAt(startInd)->getPosition());
            tsv5->set_string(0, 3, probeSets->getAt(endInd)->getProbeSetName());
            tsv5->set_i(0, 4, probeSets->getAt(endInd)->getChromosome());
            tsv5->set_i(0, 5, probeSets->getAt(endInd)->getPosition());
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("StepWindows", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "StartProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "StartChromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "StartPosition", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "EndProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 4, "EndChromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 5, "EndPosition", affx::FILE5_DTYPE_INT);

        for (int i = 0; i < stepWindows.size(); i++) {
            int startInd = stepWindows[i].first;
            int endInd   = stepWindows[i].second;
            tsv5->set_string(0, 0, probeSets->getAt(startInd)->getProbeSetName());
            tsv5->set_i(0, 1, probeSets->getAt(startInd)->getChromosome());
            tsv5->set_i(0, 2, probeSets->getAt(startInd)->getPosition());
            tsv5->set_string(0, 3, probeSets->getAt(endInd)->getProbeSetName());
            tsv5->set_i(0, 4, probeSets->getAt(endInd)->getChromosome());
            tsv5->set_i(0, 5, probeSets->getAt(endInd)->getPosition());
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("ClosestPeak", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "ProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "ClosestPeak", affx::FILE5_DTYPE_DOUBLE);

        for (int i = 0; i < closestPeak.size(); i++) {
            int index = closestPeak[i].first;
            tsv5->set_string(0, 0, probeSets->getAt(index)->getProbeSetName());
            tsv5->set_i(0, 1, probeSets->getAt(index)->getChromosome());
            tsv5->set_i(0, 2, probeSets->getAt(index)->getPosition());
            tsv5->set_d(0, 3, closestPeak[i].second);
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("CutoffShrinkages", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "FLD3Threshold", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 1, "FLD4Threshold", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 2, "Shrinkage3", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 3, "Shrinkage4", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 4, "SNPQC", affx::FILE5_DTYPE_DOUBLE);

        tsv5->set_d(0, 0, cutoffShrinkageValues.FLD3threshold);
        tsv5->set_d(0, 1, cutoffShrinkageValues.FLD4threshold);
        tsv5->set_d(0, 2, cutoffShrinkageValues.shrinkage3);
        tsv5->set_d(0, 3, cutoffShrinkageValues.shrinkage4);
        tsv5->set_d(0, 4, cutoffShrinkageValues.SNPQC);
        tsv5->writeLevel(0);

        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("MaxPeak", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "ProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "MaxPeak", affx::FILE5_DTYPE_INT);

        for (int i = 0; i < probeSets->getCount(); i++) {
            tsv5->set_string(0, 0, probeSets->getAt(i)->getProbeSetName());
            tsv5->set_i(0, 1, probeSets->getAt(i)->getChromosome());
            tsv5->set_i(0, 2, probeSets->getAt(i)->getPosition());
            tsv5->set_i(0, 3, probeSets->getAt(i)->getMaxPeaks());
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("Parameters", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "Window", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 1, "Step", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "BandwidthFactor", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 3, "DensityPointCount", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 4, "PeaksCutoff", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 5, "CleanThreshold", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 6, "Symmetry", affx::FILE5_DTYPE_INT);

        tsv5->set_i(0, 0, apParametersValues.window);
        tsv5->set_i(0, 1, apParametersValues.step);
        tsv5->set_d(0, 2, apParametersValues.bandwidthFactor);
        tsv5->set_i(0, 3, apParametersValues.densityPointCount);
        tsv5->set_d(0, 4, apParametersValues.peaksCutoff);
        tsv5->set_d(0, 5, apParametersValues.cleanThreshold);
        int sFlag = apParametersValues.symmetry ? 1 : 0;
        tsv5->set_i(0, 6, sFlag);
        tsv5->writeLevel(0);

        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("PeakSets", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "StartProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "StartChromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "StartPosition", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "EndProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 4, "EndChromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 5, "EndPosition", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 6, "Peak1", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 7, "Peak2", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 8, "Peak3", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 9, "Peak4", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 10, "Peak5", affx::FILE5_DTYPE_DOUBLE);
        tsv5->defineColumn(0, 11, "Peak6", affx::FILE5_DTYPE_DOUBLE);

        for (int j = 0; j < peakSets.size(); j++) {
            int startInd = peakSets[j].first.first;
            int endInd   = peakSets[j].first.second;
            tsv5->set_string(0, 0, probeSets->getAt(startInd)->getProbeSetName());
            tsv5->set_i(0, 1, probeSets->getAt(startInd)->getChromosome());
            tsv5->set_i(0, 2, probeSets->getAt(startInd)->getPosition());
            tsv5->set_string(0, 3, probeSets->getAt(endInd)->getProbeSetName());
            tsv5->set_i(0, 4, probeSets->getAt(endInd)->getChromosome());
            tsv5->set_i(0, 5, probeSets->getAt(endInd)->getPosition());

            int peakInd;
            for (peakInd = 0; peakInd < peakSets[j].second.size(); peakInd++) {
                tsv5->set_d(0, 6 + peakInd, peakSets[j].second[peakInd]);
            }
            for ( ; peakInd < 6; peakInd++) {
                tsv5->set_d(0, 6 + peakInd, numeric_limits<double>::quiet_NaN());
            }
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("ProcessedValues", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "ProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "Value", affx::FILE5_DTYPE_DOUBLE);

        for (int i = 0; i < processedValues.size(); i++) {
            int index = processedValues[i].first;
            tsv5->set_string(0, 0, probeSets->getAt(index)->getProbeSetName());
            tsv5->set_i(0, 1, probeSets->getAt(index)->getChromosome());
            tsv5->set_i(0, 2, probeSets->getAt(index)->getPosition());
            tsv5->set_d(0, 3, processedValues[i].second);
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        tsv5 = group5->openTsv("SCARs", affx::FILE5_REPLACE);
        tsv5->defineColumn(0, 0, "ProbesetID", affx::FILE5_DTYPE_STRING, 20);
        tsv5->defineColumn(0, 1, "Chromosome", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 2, "Position", affx::FILE5_DTYPE_INT);
        tsv5->defineColumn(0, 3, "SCAR", affx::FILE5_DTYPE_DOUBLE);

        for (int i = 0; i < probeSets->getCount(); i++) {
            tsv5->set_string(0, 0, probeSets->getAt(i)->getProbeSetName());
            tsv5->set_i(0, 1, probeSets->getAt(i)->getChromosome());
            tsv5->set_i(0, 2, probeSets->getAt(i)->getPosition());
            tsv5->set_d(0, 3, probeSets->getAt(i)->getSCAR());
            tsv5->writeLevel(0);
        }
        tsv5->close();
        delete tsv5;

        delete group5;
        file5.close();
    }
    catch(...){ throw(Except("Cannot open file " + fileName + " while attempting to dump allele peaks data to file."));}
}




double CNAnalysisMethod::assignConfidenceSigmoidal(     int iHetsInSegment,
                                                        int iNumberInSegment,
                                                        double dHetRate,
                                                        double dPError,
                                                        int iShift,
                                                        int iMinSegmentSize)
{


    if(iNumberInSegment < iMinSegmentSize)
    {
        return 0.0;
    }
    else
    {
        double dScaleMin = 1.0/(exp(shiftValue(0.0, iShift)) + 1);
        double dScaleMax = 1.0/(exp(shiftValue(1.0, iShift)) + 1);

        double dLinearValue = max(min(1.0,( 1.0-( (((double)iHetsInSegment/(double)iNumberInSegment) - dPError) / (dHetRate- dPError) )) ),0.0);

        double dTempValue1 =  1.0/(exp( shiftValue(dLinearValue, iShift) ) + 1.0) ;
        double dTempValue2 = dTempValue1 - dScaleMin;
        double dTempValue3 = dTempValue2 * ( 1.0/ (dScaleMax -dScaleMin) );

        return dTempValue3;

    }
}



double CNAnalysisMethod::shiftValue(    double ecks,
                                        int iShift)
{
    double dValue=0.0;
    dValue = 10 - (ecks*20) + iShift;
    return dValue;
}

float CNAnalysisMethod::getConfidenceThreshold(const std::string& brlmmpStr)
{
    const float defaultConfThreshold = 0.05;   // JPB perhaps this should go somewhere else

    if (brlmmpStr.empty()) {
        return defaultConfThreshold;
    }
    string name;
    map<string, string> param;
    SelfCreate::fillInNameParam(brlmmpStr, name, param);
    map<string, string>::iterator iter;
    if ((iter = param.find("MS")) != param.end()) {
        return getDouble(iter->second);
    }

    return defaultConfThreshold;
}

void CNAnalysisMethod::resetAllelePeakInitialValues()
{
    int iNumberOfProbeSets=CNAnalysisMethod::getProbeSets()->getCount();
    for (int iIndex = 0; iIndex<iNumberOfProbeSets; iIndex++)
    {
        CNAnalysisMethod::getProbeSets()->getAt(iIndex)->setAllelePeaks1(0);
        CNAnalysisMethod::getProbeSets()->getAt(iIndex)->setAllelePeaks2(0);
    }
}


// Fit density curve to specified values[] and weights[] using Gaussian kernel
// with specified bandwidth.
//
// Rather than calculating the density curve directly (which is a time-consuming
// operation) the function uses a known fact that the density curve formula
// happens to be equal to the Fourier transform of an easily computed expression
// involving the given input data. This speeds things up considerably thanks to
// the FFT.
//
// Input:
//      values - sample data to fit density curve to,
//      weights - their multiplicities (must sum up to 1),
//      numPoints - number of points to evaluate the density at,
//      bandwidth - bandwidth of the (Gaussian) kernel,
//      symmetry - whether to symmetrise the data:
//               = true, means symmetrise,
//               = false, means leave the data as is.
// Output:
//      xi - the domain coordinates of the density curve,
//      density - corresponding values at xi points.
//
void CNAnalysisMethod::fitDensityCurve( const std::vector<double>& values,
                                        const std::vector<double>& weights,
                                        std::vector<double>& xi,
                                        std::vector<double>& density,
                                        int numPoints,
                                        double bandwidth,
                                        bool symmetry
                                        )
{
    if (numPoints <= 1) {
        Err::errAbort("Number of density coordinate points must be greater than 1");
    }

    const double pi = 3.14159265358979323846;
    const double enlargeBy = 7.0 * bandwidth;
    const double margin = 4.0 * bandwidth;

    // pretend all values[] are shifted inside the interval [0, RR]
    // with extra headroom (enlargeBy) on both sides
    double minv = *std::min_element(values.begin(), values.end());
    double maxv = *std::max_element(values.begin(), values.end());
    double offset;
    double RR;
    if (symmetry) {
        if (std::fabs(minv) < std::fabs(maxv)) {
            offset = maxv + enlargeBy;
            RR = maxv + enlargeBy + offset;
        } else {
            offset = std::fabs(minv) + enlargeBy;
            RR = std::fabs(minv) + enlargeBy + offset;
        }
    } else {
        offset = -(minv - enlargeBy);
        RR = maxv + enlargeBy + offset;
    }
    
    // find nearest power of 2 larger than numPoints (for the FFT)
    int nearPow2 = 1;
    while (nearPow2 < numPoints) {
        nearPow2 *= 2;
    }

    double delta = (2.0*RR)/(2*nearPow2 - 1);

    // distribute weights[] proportionally along nearPow2 points
    // spread uniformly across [0, RR]
    std::vector<double> densityRe(2*nearPow2, 0.0);
    double delta1 = RR/(nearPow2 - 1);
    for (int iIndex = 0; iIndex < values.size(); iIndex++)
    {
        double curValue = values[iIndex] + offset;
        int intval = (int)std::floor(curValue/delta1);
 
        // NB: intval < nearPow2 - 1 (always)
        double tt = (curValue - intval*delta1)/delta1;
        densityRe[intval] += (1.0 - tt) * weights[iIndex];
        densityRe[intval + 1] += tt * weights[iIndex];
    }
    // pad with zeros to double the sampling rate
    for (int ii = nearPow2; ii < 2*nearPow2; ii++)
    {
        densityRe[ii] = 0.0;
    }
    std::vector<double> densityIm(2*nearPow2, 0.0);

    // calculate sum_i(weights[i] * exp(2 pi i t values[i]),
    // for t = m/(nearPow2 * delta), m = 0 ... (nearPow2 - 1).
    // Trick: this amounts to doing the inverse FFT on densityRe[]
    //
    doFourierTransform(densityRe, densityIm, +1);

    // plug the desired kernel function into the inverse FFT ,
    // Remember that what's needed in the end is not f(t) but
    // f(t*bandwidth). This means the kernel here must be supplied
    // as kernel(x/bandwidth)/bandwidth
    //
    std::vector<double> kernel(2*nearPow2);
    double arg = 0.0;
    for (int ii = 0; ii < nearPow2; ii++) {
        kernel[ii] = std::exp(-(arg/bandwidth)*(arg/bandwidth)/2.0) / std::sqrt(2.0 * pi)/bandwidth;
        arg += delta;
    }
    arg = -nearPow2*delta;
    for (int ii = nearPow2; ii < 2*nearPow2; ii++) {
        kernel[ii] = std::exp(-(arg/bandwidth)*(arg/bandwidth)/2.0) / std::sqrt(2.0 * pi)/bandwidth;
        arg += delta;
    }
    std::vector<double> kernelIm(2*nearPow2, 0.0);    // dummy

    doFourierTransform(kernel, kernelIm, +1);

    for (int ii = 0; ii < 2*nearPow2; ii++) {
        densityRe[ii] *= kernel[ii];
        densityIm[ii] *= kernel[ii];
    }

    // FFT the result - this gives the desired distribution
    // densityIm[] should come out identically zero
    // densityRe[] should come out nonnegative, with second half identically zero
    //
    doFourierTransform(densityRe, densityIm, -1);

    // multiply the result by 1/(2 * nearPow2) to match the true
    // Fourier transform (so the density integrates to 1)
    for (int i = 0; i < nearPow2; i++) {
        densityRe[i] /= 2 * nearPow2;
    }

    // interpolate from nearPow2 points to numPoints over the interval
    // shrunk by margin on both sides
    double eps = (RR - 2.0*margin)/(numPoints - 1);
    density.resize(numPoints);
    double curValue = margin;
    for (int ii = 0; ii < numPoints; ii++)
    {
        int intval = (int)std::floor(curValue/delta1);
        double tt = (curValue - intval*delta1)/delta1;

        // note intval+1 is within range because of the shrinking
        density[ii] = (1.0 - tt) * densityRe[intval] + tt * densityRe[intval + 1];
        curValue += eps;
    }

    // symmetrise the density if needed
    if (symmetry) {
        for (int ii = 0; ii < numPoints/2; ii++) {
            double tmp = (density[ii] + density[numPoints - ii - 1])/2.0;
            density[ii] = density[numPoints - ii - 1] = tmp;
        }
    }

    // set the xi coordinates
    xi.resize(numPoints);
    for (int ii = 0; ii < numPoints; ii++) {
        xi[ii] = margin + ii*eps - offset;
    }
}

int CNAnalysisMethod::reverseBits(int ii, int nbits)
{
    int k = 0;
    for (int i = 0; i < nbits; i++) {
        int j = ii/2;
        k += k;
        k += ii - 2*j;
        ii = j;
    }
    return k;
}

void CNAnalysisMethod::doBitReverseOrder(vector<double>& vec)
{
    int nbits = 1;
    int N = vec.size();
    while ((N /= 2) != 1)
    {
        nbits++;
    }
    for (int i = 1; i < vec.size() - 1; i++) {
        int rIndex = reverseBits(i, nbits);
        if (rIndex > i) {
            std::swap( vec[i], vec[rIndex] );
        }
    }
}

// Fast Fourier Transform, requires input vectors' length = 2^k.
// Performs the transform in-place, overwriting re[] and im[].
// Input:
//      re - real part of the input
//      im - imagianry part of the input
//      expSign - sign of the exponent of the FFT kernel:
//              = -1 means do the forward transform,
//              = +1 means do the inverse transform.
// Output:
//      re - real part of the output
//      im - imaginary part of the output
//
void CNAnalysisMethod::doFourierTransform(vector<double>& re, vector<double>& im, int expSign)
{
    const double pi = 3.14159265358979323846;

    doBitReverseOrder(re);
    doBitReverseOrder(im);

    int N = re.size();       // must be a power of 2, equal to im.size()

    int iStep;
    double exponent = expSign * pi;
    for (int len = 1; len < N; len = iStep) {
        iStep = 2*len;
        double twiddle_re = cos(exponent/len);
        double twiddle_im = sin(exponent/len);
        for (int i = 0; i < N; i += iStep) {
            double arg_re = 1.0;
            double arg_im = 0.0;
            for (int j = i; j < i + len; j++) {
                double oddRe = arg_re * re[j + len] - arg_im * im[j + len];
                double oddIm = arg_re * im[j + len] + arg_im * re[j + len];

                // butterfly
                re[j + len] = re[j] - oddRe;
                im[j + len] = im[j] - oddIm;
                re[j] += oddRe;
                im[j] += oddIm;

                // twiddle update
                double tmp = arg_re;
                arg_re = arg_re * twiddle_re - arg_im * twiddle_im;
                arg_im = tmp    * twiddle_im + arg_im * twiddle_re;
            }
        }
    }
}
