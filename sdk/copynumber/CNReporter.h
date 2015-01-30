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

#ifndef _CNReporter_H_
#define _CNReporter_H_
/**
 * @file CNReporter.h
 *
 * @brief This header contains the CNReporter class definition.
 */

#include "copynumber/CNAnalysisMethod.h"
#include "copynumber/CNExperiment.h"
#include "copynumber/CNSegment.h"
//
#include "util/AffxArray.h"
#include "util/AffxConv.h"
#include "util/BaseEngine.h"
#include "util/CalvinToText.h"
#include "util/Err.h"
#include "util/Guid.h"
#include "util/Util.h"
#include "util/Verbose.h"
//

#ifndef _MSC_VER
#include <sys/types.h>
#include <netinet/in.h>
#include <inttypes.h>
#endif

using namespace std;

class CNReporterMethod : public CNAnalysisMethod
{
protected:
  affymetrix_calvin_io::DataSetWriter* m_pset;

public:
  CNReporterMethod();
  void setDataSetWriter(affymetrix_calvin_io::DataSetWriter& set);
  virtual int getRowCount();
};

/**
 * @brief  A base class for copy number reporters.
 *
 */
class CNReporter
{
protected:
  BaseEngine* m_pEngine;
  CNExperiment* m_pobjExperiment;
  CNProbeSetArray* m_pvProbeSets;
  AffxArray<CNAnalysisMethod>* m_pvMethods;

  AffxString m_strARRFileName;

  void wavinessSegCounts(int& segCountLoss, int& segCountGain, float& sd);
  std::string wavinessAmplitudes();

public:
  CNReporter();
  virtual ~CNReporter();

  CNExperiment* getExperiment();
  CNProbeSetArray* getProbeSets();
  AffxArray<CNAnalysisMethod>* getMethods();

  virtual void defineOptions(BaseEngine& e);
  virtual void checkOptions(BaseEngine& e);

  /**
   * @brief Setup the reporter to be run
   * @param BaseEngine& - The engine associated with this reporter
   * @param CNExperiment& - The experiment
   * @param CNProbeSetArray& - The probe set vector associated with the experiment
   * @param AffxArray<CNAnalysisMethod>& - The analysis methods to extract data from
   */
  void setup(BaseEngine& engine, CNExperiment& objExperiment, CNProbeSetArray& vProbeSets, AffxArray<CNAnalysisMethod>& vMethods);

  /**
   * @brief Is the reporter setup to be rrun
   */
  void isSetup();
  virtual void run() = 0;

  AffxString getArrayName();

  /**
   * @brief Load a header parameter from individual comppnents
   * @param const AffxString& - The parameter name
   * @param PgOpt::PgOptType - The parameter type
   * @param const AffxString& - The parameter value
   * @param affymetrix_calvin_parameter::ParameterNameValueType& - The parameter to load
   */
  void loadParam(const AffxString& strName, PgOpt::PgOptType type, const AffxString& strValue, affymetrix_calvin_parameter::ParameterNameValueType& param);

  static void loadBuffer(char* pBuffer, int& iIndex, AffxString& str, int iLength = -1);

  static void loadBuffer(char* pBuffer, int& iIndex, unsigned char uc);

  static void loadBuffer(char* pBuffer, int& iIndex, char c);

  static void loadBuffer(char* pBuffer, int& iIndex, unsigned int ui);

  static void loadBuffer(char* pBuffer, int& iIndex, float f);

  static AffxString prepareAscii(const AffxString& str, int iLength);

};

struct wavSortCriterion : binary_function<pair<int, float>, pair<int, float>, bool>
{
    bool operator()(const pair<int, float>& lhs, const pair<int, float>& rhs) const {
        return lhs.first < rhs.first;
    }
};

#endif


