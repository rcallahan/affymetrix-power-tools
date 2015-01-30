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

#ifndef _QUANTMETHODSQLEXPRREPORT_H_
#define _QUANTMETHODSQLEXPRREPORT_H_

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
class QuantMethodSqlExprReport : public QuantMethodReport {

public:

  /**
   * Constructor
   */
  QuantMethodSqlExprReport(const std::string &databaseFile, const std::string &tableName, int psNameLength=30);

  /**
   * Virtual destructor for a virtual class.
   */
  virtual ~QuantMethodSqlExprReport();

  /** 
   * Get set up for a run of reporting probesets. Often used to open file
   * streams and print headers to files etc.
   * 
   * @param qMethod - Quantification method to be used.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool prepare(QuantMethod &qMethod, const IntensityMart &iMart);

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
  virtual bool report(ProbeSetGroup &psGroup, 
                      QuantMethod &qMethod, 
                      const IntensityMart &iMart, 
                      std::vector<ChipStream *> &iTrans, 
                      PmAdjuster &pmAdjust);

  /** 
   * If a probeset fails to compute for whatever reason then this method is
   * called rather than the normal report call above. By default does nothing.
   * 
   * @param psGroup - List of probesets from which probes were used.
   * @param qMethod - Quantification method with compute method called.
   * @param layout - Where the probesets, probes, etc are on the chip.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool reportFailure(ProbeSetGroup &psGroup, 
                             QuantMethod &qMethod,
                             const IntensityMart &iMart, 
                             std::vector<ChipStream *> &iTrans, 
                             PmAdjuster &pmAdjust) { return true; }


  /** 
   * No more probesets will be processed, this is a chance to finish outputting
   * results and clean up.
   * 
   * @param qMethod - Quantification method that was used.
   * 
   * @return true if success, false otherwise.
   */
  virtual bool finish(QuantMethod &qMethod);
    
  virtual void addStdHeaders(QuantMethodReport *qReport,
                             const std::string& execGuid, 
                             const std::string& reportGuid,
                             const std::string& timeStr,
                             const std::string& commandLine,
                             const std::string& execVersion,
                             const AnalysisInfo& info);

  static void execute(sqlite3 *db, const std::string &sql);

private:

  void createSummaryTable(sqlite3 *db, const std::string &tableName, const IntensityMart &iMart);

  void createMetaTable(sqlite3 *db, const std::string &metaTableName);

  /// Size of longest name
  int m_MaxPsNameLength;
  /// Database that we are writing to
  sqlite3 *m_Db;
  /// Name of the file that we're using for database.
  std::string m_DbFileName;
  /// Name of the database table that we are writing into
  std::string m_TableName;
  /// Name of the database table that contains the meta information usually found in headers
  std::string m_TableMetaName;

};

#endif /* _QUANTMETHODSQLEXPRREPORT_H_ */
 
