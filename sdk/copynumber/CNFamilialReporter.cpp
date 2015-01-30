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
 * @file CNFamilialReporter.cpp
 *
 * @brief This file contains the CNFamilialReporter class members.
 */

#include "copynumber/CNFamilialReporter.h"
#include "portability/affy-base-types.h"
#include "util/AffxStatistics.h"
#include "util/Fs.h"
//
#include <limits>
//

int CNFamilialReporter::m_iInstanceCount = 0;
std::vector<affymetrix_calvin_parameter::ParameterNameValueType> CNFamilialReporter::m_vParams;

/**
 * @brief Constructor
 */
CNFamilialReporter::CNFamilialReporter()
{
  m_iInstanceCount++;
  m_pEngine = NULL;
  m_pcychpIndex = NULL;
  m_pcychpMother = NULL;
  m_pcychpFather = NULL;

  m_pvMethods = NULL;
}

/**
 * @brief Destructor
 */
CNFamilialReporter::~CNFamilialReporter()
{
  m_iInstanceCount--;
  if (m_iInstanceCount == 0) {
    m_vParams.clear();
  }
}

AffxString CNFamilialReporter::getPrefix()
{
  return "affymetrix-algorithm-param-";
}

std::vector<affymetrix_calvin_parameter::ParameterNameValueType>& CNFamilialReporter::getParams()
{
  return m_vParams;
}

AffxString CNFamilialReporter::getName()
{
  return "";
}

/**
 * @brief Setup the reporter to be run
 * @param BaseEngine& - The engine associated with this reporter
 * @param CNCychp& - The Index's cychp data
 * @param CNCychp& - The Mother's cychp data
 * @param CNCychp& - The Father's cychp data
 * @param AffxArray<CNFamilialAnalysisMethod>& - The analysis methods to extract data from
 */
void CNFamilialReporter::setup(BaseEngine& engine, CNCychp& cychpIndex, CNCychp& cychpMother, CNCychp& cychpFather, AffxArray<CNFamilialAnalysisMethod>& vMethods)
{
  m_pEngine = &engine;
  m_pcychpIndex = &cychpIndex;
  m_pcychpMother = &cychpMother;
  m_pcychpFather = &cychpFather;

  m_pvMethods = &vMethods;
}

CNCychp& CNFamilialReporter::getCychpIndex()
{
  return *m_pcychpIndex;
}
CNCychp& CNFamilialReporter::getCychpMother()
{
  return *m_pcychpMother;
}
CNCychp& CNFamilialReporter::getCychpFather()
{
  return *m_pcychpFather;
}

/**
 * Return SelfDoc option associated with a specified name.
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc options to search.
 * @param const std::string& - The specified name to search for.
 * @return SelfDoc::Opt* - The SelfDoc option pointer or NULL if not found.
 */
SelfDoc::Opt* CNFamilialReporter::getSelfDocOpt(std::vector<SelfDoc::Opt>& opts, const std::string& strName)
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
bool CNFamilialReporter::setupBoolParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
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
int CNFamilialReporter::setupIntParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
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
float CNFamilialReporter::setupFloatParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
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
 * Setup a string parameter for this analysis method, including min amd max validation.
 * @param const std::string& - The name of the parameter
 * @param const std::string& - The prefix to use when setting up the header
 * @param std::map<std::string,std::string>& - A map of parameters to fill in
 * @param SelfDoc& - The SelfDoc object needed by the fillInValue function
 * @param std::vector<SelfDoc::Opt>& - The SelfDoc option to setup.
 * @return std::string - The value of the parameter as setup by this function.
 */
std::string CNFamilialReporter::setupStringParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts)
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
 * @brief Return the analysis associated with the specified name
 * @param const AffxString& - The specified name
 * @return CNFamilialAnalysisMethod* - The associated analysis pointer
 */
CNFamilialAnalysisMethod* CNFamilialReporter::getAnalysis(const AffxString& strName)
{
  for (int iIndex = 0; (iIndex < (int)m_pvMethods->size()); iIndex++) {
    if (m_pvMethods->at(iIndex)->getName() == strName) {
      return m_pvMethods->at(iIndex);
    }
  }
  return NULL;
}

void CNFamilialReporter::loadParam(const AffxString& strName, PgOpt::PgOptType type, const AffxString& strValue, affymetrix_calvin_parameter::ParameterNameValueType& param)
{
  param.SetName(StringUtils::ConvertMBSToWCS(strName));
  switch (type) {
  case PgOpt::BOOL_OPT: param.SetValueInt8(((strValue == "true") ? (char)1 : (char)0)); break;
  case PgOpt::DOUBLE_OPT: param.SetValueFloat((float)::getDouble(strValue)); break;
  case PgOpt::INT_OPT: param.SetValueInt32(::getInt(strValue)); break;
  case PgOpt::STRING_OPT:
    if ((strName.indexOf("command-line") != -1) || (strName.indexOf("analysis") != -1) || (strName.indexOf("program-cvs-id") != -1) || (strName.indexOf("version-to-report") != -1) || (strName.endsWith("-dir"))) {
      param.SetValueText(StringUtils::ConvertMBSToWCS(strValue));
    } else {
      param.SetValueText(StringUtils::ConvertMBSToWCS(Fs::basename(strValue)));
    }
    break;
  default: throw(Except("Cannot find PgOpt type for: " + strName));
  }
}
