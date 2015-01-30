////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
//
#include "chipstream/QuantMethodSqlGTypeReport.h"
//
#include "chipstream/QuantGTypeMethod.h"
#include "chipstream/QuantMethodSqlExprReport.h"

using namespace std;

QuantMethodSqlGTypeReport::QuantMethodSqlGTypeReport(const std::string &databaseFile, const std::string &callTable, const std::string &confTable, int maxNameLength) {
  m_CallDbFileName = databaseFile + ".call" + ".db";
  m_ConfDbFileName = databaseFile + ".conf" + ".db";
  m_CallTable = callTable;
  m_CallMetaTable = callTable + "_meta";
  m_ConfTable = confTable;
  m_ConfMetaTable = confTable + "_meta";
  m_MaxPsNameLength = maxNameLength;
  int status = 0;
  if ((status =  sqlite3_open(m_CallDbFileName.c_str(), &m_CallDb)) != SQLITE_OK) {
    Err::errAbort("Couldn't open database: " + m_CallDbFileName + " got error code: " + ToStr(status));
  }
  if ((status =  sqlite3_open(m_ConfDbFileName.c_str(), &m_ConfDb)) != SQLITE_OK) {
    Err::errAbort("Couldn't open database: " + m_ConfDbFileName + " got error code: " + ToStr(status));
  }
}

QuantMethodSqlGTypeReport::~QuantMethodSqlGTypeReport() {
  int status = 0;
  // Close the sqlite database 
  if ((status = sqlite3_close(m_CallDb)) != SQLITE_OK) {
   Verbose::warn(0, "Error trying to close the database with error code: " + ToStr(status));
  }
  if ((status = sqlite3_close(m_ConfDb)) != SQLITE_OK) {
   Verbose::warn(0, "Error trying to close the database with error code: " + ToStr(status));
  }
}

bool QuantMethodSqlGTypeReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
  string sql = "begin transaction;";
  createCallTable(m_CallDb, m_CallTable, iMart);
  QuantMethodSqlExprReport::execute(m_CallDb, sql);
  createConfTable(m_ConfDb, m_ConfTable, iMart);
  QuantMethodSqlExprReport::execute(m_ConfDb, sql);
  return true;
}

bool QuantMethodSqlGTypeReport::report(ProbeSetGroup &psGroup, 
                                       QuantMethod &qMethod, 
                                       const IntensityMart &iMart, 
                                       std::vector<ChipStream *> &iTrans, 
                                       PmAdjuster &pmAdjust) {
  static char buf[256];
  static string delim = ",";
  QuantGTypeMethod *gMethod = dynamic_cast<QuantGTypeMethod *>(&qMethod);
  if(gMethod == NULL) {
    Err::errAbort("Can only use a QuantMethodGTypeReport with QuantGTypeMethods.");
  }
  std::string name = gMethod->getProbeSetName();
  int targetCount = gMethod->getNumCalls();

  string insertCalls = "insert into " + m_CallTable + " values ( '" + psGroup.name + "' ";
  string insertConf = "insert into " + m_ConfTable + " values ( '" + psGroup.name + "' ";
  for(int i = 0; i < targetCount; i++) {
    snprintf(buf, ArraySize(buf), "%.3g", gMethod->getConfidence(i));
    insertConf += delim + buf;
    insertCalls += delim + ToStr((int)gMethod->getCall(i));
  }
  insertCalls += ");";
  insertConf  += ");";
  QuantMethodSqlExprReport::execute(m_CallDb, insertCalls);
  QuantMethodSqlExprReport::execute(m_ConfDb, insertConf);
  return true;
}

void QuantMethodSqlGTypeReport::createConfTable(sqlite3 *db, const std::string &tableName, const IntensityMart &iMart) {
  string createSql = "create table " + tableName + " ( probeset_id varchar(" + ToStr(m_MaxPsNameLength) + "), ";
  for(int i = 0; i < iMart.getCelFileNames().size(); i++) {
    createSql += " '" + iMart.getCelFileNames()[i] + "' real ";
    if(i + 1 != iMart.getCelFileNames().size()) {
      createSql += ",";
    }
  }
  createSql += " );";
  QuantMethodSqlExprReport::execute(db, createSql);
}

void QuantMethodSqlGTypeReport::createCallTable(sqlite3 *db, const std::string &tableName, const IntensityMart &iMart) {
  string createSql = "create table " + tableName + " ( probeset_id varchar(" + ToStr(m_MaxPsNameLength) + ") , ";
  for(int i = 0; i < iMart.getCelFileNames().size(); i++) {
    createSql += " '" + iMart.getCelFileNames()[i] + "' integer(1) ";
    if(i + 1 != iMart.getCelFileNames().size()) {
      createSql += ",";
    }
  }
  createSql += " );";
  QuantMethodSqlExprReport::execute(db, createSql);
}

void QuantMethodSqlGTypeReport::createMetaTable(sqlite3 *db, const std::string &metaTableName) {
  string createSql = "create table " + metaTableName + " ( key varchar(256), value text ); ";
  QuantMethodSqlExprReport::execute(db, createSql);
}

void QuantMethodSqlGTypeReport::addStdHeaders(QuantMethodReport *qReport,
                                             const std::string& execGuid, 
                                             const std::string& reportGuid,
                                             const std::string& timeStr,
                                             const std::string& commandLine,
                                              const std::string& execVersion,
                                             const AnalysisInfo& info) {
  std::vector<std::string>::const_iterator keyIx, paramIx;
  createMetaTable(m_CallDb, m_CallMetaTable);
  createMetaTable(m_ConfDb, m_ConfMetaTable);
  //
  for(keyIx = info.m_ParamNames.begin(), paramIx = info.m_ParamValues.begin();
      keyIx != info.m_ParamNames.end() && paramIx != info.m_ParamValues.end();
      ++keyIx, ++paramIx) {
    string sql = "insert into " + m_CallMetaTable + " values ( '" + *keyIx + "', '" + *paramIx + "');";
    QuantMethodSqlExprReport::execute(m_CallDb, sql);
    sql = "insert into " + m_ConfMetaTable + " values ( '" + *keyIx + "', '" + *paramIx + "');";
    QuantMethodSqlExprReport::execute(m_ConfDb, sql);
  }
}

bool QuantMethodSqlGTypeReport::finish(QuantMethod &qMethod) {
  string sql = "commit;";
  QuantMethodSqlExprReport::execute(m_CallDb, sql);
  QuantMethodSqlExprReport::execute(m_ConfDb, sql);

  string indexCalls = "create index call_index on call (probeset_id);";
  QuantMethodSqlExprReport::execute(m_CallDb,  indexCalls);

  string indexConf = "create index confidence_index on confidence (probeset_id);";
  QuantMethodSqlExprReport::execute(m_ConfDb,  indexConf);
  return true;
}
