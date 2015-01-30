////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////
/**
 * @file   RegressionCheck.h
 * @author Chuck Sugnet
 * @date   Tue Apr 25 16:02:08 2006
 * 
 * @brief  Virtual class for checks to be done after a regression test.
 */

#ifndef REGRESSIONCHECK_H
#define REGRESSIONCHECK_H

#include "util/Util.h"
#include "util/Verbose.h"
//
#include "portability/affy-base-types.h"
//
#include <cstring>
#include <string>
//

/**
 * Abstract base class (i.e. interface) for checks to be run during a regression
 * test.
 */
class RegressionCheck {

public:
 
  /** Constructor */
  RegressionCheck() : m_Name("test"), m_NegTest(false), m_MaxErrorsReport(30), m_CurrentErrorsReported(0)
  { }

  /** Virtual destructor for virtual class. */
  virtual ~RegressionCheck() {}

  /** 
   * Check the condition expected.
   * @param msg - Message to be optionally set while checking.
   * @return - true if passed, false otherwise.
   */
  virtual bool check(std::string &msg) = 0;

protected:
  /** 
   * Utility function to see if a condition is true and if not to add
   * the message to our summary of messages.
   * @param condition - boolean to be tested.
   * @param msg - informative message for this condition (if fails)
   * @param summary - ongoing collection of messages.
   * @return - condition test results.
   */
  bool checkMsg(bool condition, const std::string &msg, std::string &summary) {
    if(!condition) {
      summary += msg;
      reportError(msg);
    }
    return condition;
  }

  /** 
   * Utility function to accumulate error messages
   * @param err - error message
   * @return - void
   */
  void reportError(const std::string &err) {
    // APT-992
#ifndef _WIN32    
	if(m_CurrentErrorsReported++ < m_MaxErrorsReport) {
      Verbose::out(1, "Error encountered: " + err);
      if(m_CurrentErrorsReported >= m_MaxErrorsReport) {
        Verbose::out(1, "Maximum number of errors reported.");
      }
    }
#endif
  }

  /** 
   * Utility function to set the max number of errors to report
   * @param max - maximum number of errors to report (-1 for no limit)
   * @return - void
   */
  void setMaxError(int max) {
    m_MaxErrorsReport = max;
  }

  /** 
   * Are two floats within epsilon (small value) of eachother?
   * 
   * @param gold - correct value.
   * @param gen - value to be tested.
   * @param eps - epsilon, small number they can be different but still equivalent.
   *            i.e. if |generated-gold| >= eps then there is a difference.
   * @param success - boolean set to false if floats not equivalent.
   * @param maxDiff - Maximum difference seen so far.
   * @param log - boolean set to true if differences to be reported.
   * @param frac - Maximum accepted fractional difference in numeric values (not used by default).
   *            i.e. if |generated-gold| >= frac*max(|generated|,|gold|) then there is a difference.
   * @return - true if floats equivalent, false otherwise.
   */
  bool checkFloat(float gold, float gen, double eps, bool &success, double &maxDiff, bool log = false, double frac = 0) {
    bool localSuccess = true;
    // bug fix: proper equivalence test when one or both values to compare are non-finite (see ClearQuest AFFY00024123)
      if (!Util::isFinite(gold) && !Util::isFinite(gen)) {
        // test pair of non-finite values for equivalence (both NAN or both INF with same sign)
        if ( (gold>0)==(gen>0) && (gold<0)==(gen<0) ) {
            return success;
        }
        success = localSuccess = false;
        if (log) {
            reportError("Non-finite non-equivalent floats. gold: '" + ToStr(gold) + "' test: '" + ToStr(gen) + "'");
        }
        return localSuccess;
      }
      if (!Util::isFinite(gold) || !Util::isFinite(gen)) {
        // finite and non-finite values can never be equivalent
        success = localSuccess = false;
        if (log) {
            reportError("One float non-finite difference. gold: '" + ToStr(gold) + "' test: '" + ToStr(gen) + "'");
        }
        return localSuccess;
      }
      
      double diff = fabs(gold - gen);
      maxDiff = Max(diff, maxDiff);
      
      // allowed absolute difference from fractional tolerance (zero by default)
      double eps2 = frac*Max( fabs(gold), fabs(gen) );
      // absolute difference is acceptable if it satisfies either (least restrictive) tolerance
      if (diff > Max(eps,eps2)) {
        success = localSuccess = false;
        if (log) {
            reportError("Different floats. gold: '" + ToStr(gold) + "' test: '" + ToStr(gen) + "'");
        }
        return localSuccess;
      }
      return localSuccess;
  }
 
public:
  std::string m_Name;
  bool m_NegTest;

private:

  int m_MaxErrorsReport;        /// The maximum number of errors to report
  int m_CurrentErrorsReported;  /// The current number of reported errors

};


#endif /* REGRESSIONCHECK_H */
