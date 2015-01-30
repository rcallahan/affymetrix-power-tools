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
#include "chipstream/QuantMethodSqlExprReport.h"

using namespace std;

QuantMethodSqlExprReport::QuantMethodSqlExprReport(const std::string &databaseFile, const std::string &tableName, int maxLength) {
  m_DbFileName = databaseFile;
  m_TableName = tableName;
  m_TableMetaName = tableName + "_meta";
  m_MaxPsNameLength = maxLength;
  int status = 0;
  if ((status =  sqlite3_open(databaseFile.c_str(), &m_Db)) != SQLITE_OK) {
    Err::errAbort("Couldn't open database: " + databaseFile + " got error code: " + ToStr(status));
  }
}

QuantMethodSqlExprReport::~QuantMethodSqlExprReport() {
  int status = 0;
  // Close the sqlite database 
  if ((status = sqlite3_close(m_Db)) != SQLITE_OK) {
   Verbose::warn(0, "Error trying to close the database with error code: " + ToStr(status));
  }
}

bool QuantMethodSqlExprReport::prepare(QuantMethod &qMethod, const IntensityMart &iMart) {
  createSummaryTable(m_Db, m_TableName, iMart);
  string sql = "begin transaction;";
  execute(m_Db, sql);
  return true;
}

bool QuantMethodSqlExprReport::report(ProbeSetGroup &psGroup, 
                                      QuantMethod &qMethod, 
                                      const IntensityMart &iMart, 
                                      std::vector<ChipStream *> &iTrans, 
                                      PmAdjuster &pmAdjust) {
  static char buf[256];
  static string delim = ",";
  if(psGroup.probeSets[0]->psType != ProbeSet::Expression && psGroup.probeSets[0]->psType != ProbeSet::Copynumber) {
    return false;
  }
  QuantExprMethod* qeMethod = dynamic_cast<QuantExprMethod *>(&qMethod);
  if(qeMethod == NULL) {
    Err::errAbort("Can't call QuantMethodSqlExprReport::report() with something other than a QuantExprMethod.");
  }
  int targetCount = qeMethod->getNumTargets();

  string insert = "insert into " + m_TableName + " values ( '" + psGroup.name + "' ";
  for(int i = 0; i < targetCount; i++) {
    snprintf(buf, ArraySize(buf), "%.4g", qeMethod->getSignalEstimate(i));
    insert += delim + buf;
  }
  insert += ");";
  execute(m_Db, insert);

  return true;
}


void QuantMethodSqlExprReport::createSummaryTable(sqlite3 *db, const std::string &tableName, const IntensityMart &iMart) {
  string createSql = "create table " + tableName + " ( probeset_id varchar(" + ToStr(m_MaxPsNameLength) + ") , ";
  for(int i = 0; i < iMart.getCelFileNames().size(); i++) {
    createSql += " '" + iMart.getCelFileNames()[i] + "' real ";
    if(i + 1 != iMart.getCelFileNames().size()) {
      createSql += ",";
    }
  }
  createSql += " );";
  execute(db, createSql);
}

void QuantMethodSqlExprReport::createMetaTable(sqlite3 *db, const std::string &metaTableName) {
  string createSql = "create table " + metaTableName + " ( key varchar(256), value text ); ";
  execute(db, createSql);
}



void QuantMethodSqlExprReport::execute(sqlite3 *db, const string &sql) {
  
  sqlite3_stmt *stmt = NULL;  
  int status = 0;
    if ((status = sqlite3_prepare_v2(db, sql.c_str(), sql.size(), &stmt, NULL)) != SQLITE_OK) {
      Err::errAbort("Error error code: " + ToStr(status) + " preparing statement: " + sql);
    }
    
    if ((status = sqlite3_step(stmt)) != SQLITE_DONE) {
      Err::errAbort("Error error code: " + ToStr(status) + " executing statement: " + sql);
    }
    
    if ((status = sqlite3_finalize(stmt)) != SQLITE_OK) {
      Err::errAbort("Error error code: " + ToStr(status) + " finalizing statement: " + sql);
    }
}

bool QuantMethodSqlExprReport::finish(QuantMethod &qMethod) {
    string sql = "commit;";
    execute(m_Db, sql);
    string index = "create index summary_index on summary (probeset_id);";
    execute(m_Db, index);
    return true;
}

void QuantMethodSqlExprReport::addStdHeaders(QuantMethodReport *qReport,
                                             const std::string& execGuid, 
                                             const std::string& reportGuid,
                                             const std::string& timeStr,
                                             const std::string& commandLine,
                                              const std::string& execVersion,
                                             const AnalysisInfo& info) {
  std::vector<std::string>::const_iterator keyIx, paramIx;
  createMetaTable(m_Db, m_TableMetaName);
  string sql = "begin transaction;";
  execute(m_Db, sql);
  //
  for(keyIx = info.m_ParamNames.begin(), paramIx = info.m_ParamValues.begin();
      keyIx != info.m_ParamNames.end() && paramIx != info.m_ParamValues.end();
      ++keyIx, ++paramIx) {
    string sql = "insert into " + m_TableMetaName + " values ( '" + *keyIx + "', '" + *paramIx + "');";
    execute(m_Db, sql);
  }
  sql = "commit;";
  execute(m_Db, sql);
}
