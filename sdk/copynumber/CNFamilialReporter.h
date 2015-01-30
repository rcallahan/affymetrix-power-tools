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

#ifndef _CNFamilialReporter_H_
#define _CNFamilialReporter_H_
/**
 * @file CNFamilialReporter.h
 *
 * @brief This header contains the CNFamilialReporter class definition.
 */

#include "copynumber/CNCychp.h"
#include "copynumber/CNFamilialAnalysisMethod.h"
//
#include "calvin_files/utils/src/StringUtils.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
#include "util/BaseEngine.h"
#include "util/CalvinToText.h"
#include "util/Err.h"
#include "util/Guid.h"
#include "util/Util.h"
#include "util/Verbose.h"
//

/**
 * @brief  A base class for copy number quant methods.
 *
 */
class CNFamilialReporter : public SelfDoc, public SelfCreate
{
protected:
  static int m_iInstanceCount;
  static std::vector<affymetrix_calvin_parameter::ParameterNameValueType> m_vParams;

protected:
  static AffxString getPrefix();
  static SelfDoc::Opt* getSelfDocOpt(std::vector<SelfDoc::Opt>& opts, const std::string& strName);
  static bool setupBoolParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
  static int setupIntParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
  static float setupFloatParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);
  static std::string setupStringParameter(const std::string& strName, const std::string& strPrefix, std::map<std::string, std::string>& params, SelfDoc& doc, std::vector<SelfDoc::Opt>& opts);

public:
  static std::vector<affymetrix_calvin_parameter::ParameterNameValueType>& getParams();

protected:
  BaseEngine* m_pEngine;
  CNCychp* m_pcychpIndex;
  CNCychp* m_pcychpMother;
  CNCychp* m_pcychpFather;
  AffxArray<CNFamilialAnalysisMethod>* m_pvMethods;

public:
  CNFamilialReporter();
  virtual ~CNFamilialReporter();

  virtual AffxString getName();

  /**
   * @brief Setup the reporter to be run
   * @param BaseEngine& - The engine associated with this reporter
   * @param CNCychp& - The Index's cychp data
   * @param CNCychp& - The Mother's cychp data
   * @param CNCychp& - The Father's cychp data
   * @param AffxArray<CNFamilialAnalysisMethod>& - The analysis methods to extract data from
   */
  virtual void setup(BaseEngine& engine, CNCychp& cychpIndex, CNCychp& cychpMother, CNCychp& cychpFather, AffxArray<CNFamilialAnalysisMethod>& vMethods);

  CNCychp& getCychpIndex();
  CNCychp& getCychpMother();
  CNCychp& getCychpFather();

  virtual void run() = 0;

protected:
  CNFamilialAnalysisMethod* getAnalysis(const AffxString& strName);

  /**
   * @brief Load a header parameter from individual comppnents
   * @param const AffxString& - The parameter name
   * @param PgOpt::PgOptType - The parameter type
   * @param const AffxString& - The parameter value
   * @param affymetrix_calvin_parameter::ParameterNameValueType& - The parameter to load
   */
  void loadParam(const AffxString& strName, PgOpt::PgOptType type, const AffxString& strValue, affymetrix_calvin_parameter::ParameterNameValueType& param);
};

#endif


