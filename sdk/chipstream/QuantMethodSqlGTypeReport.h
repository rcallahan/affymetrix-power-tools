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
 * @file   QuantMethodSqlGTypeReport.h
 * @author Chuck Sugnet
 * @date   Mon Nov 30 14:36:41 PST 2009
 *
 * @brief  Class for reporting results of quantification methods into an sql database
 */

#ifndef _QUANTMETHODSQLGTYPEREPORT_H_
#define _QUANTMETHODSQLGTYPEREPORT_H_

//
#include "chipstream/QuantExprMethod.h"
#include "chipstream/QuantMethodReport.h"
#include "chipstream/TsvReport.h"
//
#include "util/Util.h"
//
#include <cstdlib>
#include <iostream>

#include "sqlite3.h"
//

/**
 *   Class for reporting results of quantification methods.
 */
class QuantMethodSqlGTypeReport : public QuantMethodReport {

public:
  
  QuantMethodSqlGTypeReport(const std::string &databaseFile, const std::string &callName, 
                            const std::string &confTable, int maxPsName=30);

  /** Destructor. */
  ~QuantMethodSqlGTypeReport();

  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param qMethod - Quantification method to be used.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  bool prepare(QuantMethod &qMethod, const IntensityMart &iMart);

  /** 
   * After every probeset computation this function is called an is an opportunity
   * to query the quantification method for results, residuals, etc.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  bool report(ProbeSetGroup& psGroup,
              QuantMethod& qMethod, 
              const IntensityMart& iMart, 
              std::vector<ChipStream *>& iTrans, 
              PmAdjuster& pmAdjust);
    
  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * 
   * @param qMethod - Quantification method that was used.
   * 
   * @return true if success, false otherwise.
   */
  bool finish(QuantMethod &qMethod);

  virtual void addStdHeaders(QuantMethodReport *qReport,
                     const std::string& execGuid, 
                     const std::string& reportGuid,
                     const std::string& timeStr,
                     const std::string& commandLine,
                     const std::string& execVersion,
                     const AnalysisInfo& info);

private:

  void createConfTable(sqlite3 *db, const std::string &tableName, const IntensityMart &iMart);

  void createCallTable(sqlite3 *db, const std::string &tableName, const IntensityMart &iMart);

  void createMetaTable(sqlite3 *db, const std::string &metaTableName);

  /*
   * Currently two separate databases as dbs are locked during write since bulk
   * commits are much faster than individual commits. Could easily do one by
   * caching results and then bulk loading them before commiting to keep
   * performance high and lock time minimized. 
   */
  sqlite3 *m_CallDb;
  sqlite3 *m_ConfDb;
  int m_MaxPsNameLength;
  std::string m_CallDbFileName;
  std::string m_ConfDbFileName;
  std::string m_CallTable;
  std::string m_CallMetaTable;
  std::string m_ConfTable;
  std::string m_ConfMetaTable;
};

#endif /* _QUANTMETHODSQLGTYPEREPORT_H_ */
 
