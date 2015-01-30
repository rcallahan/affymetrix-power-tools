////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
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

//
#include "file5/File5_Sqlite.h"
//
#include "file5/File5.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include "sqlite3.h"
//
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace affx;

void File5_Sqlite::makeCreateStatement(File5_Tsv *tsv, const string &primaryColumn,
                                       const string &tableName, const string &numeric,
                                       string &createSql) {
  createSql = "create table " + tableName + " ( ";
  int colCount = tsv->getColumnCount(0);
  for ( int i = 0; i < colCount; i++) {
    string name;
    tsv->getColumnName(0,i,&name);
    bool isPrimaryKey = (name == primaryColumn);
    affx::File5_TsvColumn* column = tsv->getColumnPtr(0,i);
    File5_dtype_t type = column->file5_dtype();
    createSql += " '" + name + "' ";
    if (type == FILE5_DTYPE_STRING) {
      int stringSize = column->getOptStringSize();
      createSql += " varchar(" + ToStr(stringSize) + ") ";
    }
    else if (type == FILE5_DTYPE_CHAR) {
      createSql += " int(1) ";
    }
    else if (type == FILE5_DTYPE_SHORT) {
      createSql += " int(2) ";
    }
    else if (type == FILE5_DTYPE_INT) {
      createSql += " int ";
    }
    else if (type == FILE5_DTYPE_FLOAT) {
      createSql += numeric;
    }
    else if (type == FILE5_DTYPE_DOUBLE) {
      createSql += numeric;
    }
    else {
      Err::errAbort("Don't recognize type: " + ToStr(type));
    }

    if (isPrimaryKey) {
      createSql += " primary key ";
    }
    if ((i + 1) != colCount)
      createSql += ", ";
  }
  createSql += ");";
}

void File5_Sqlite::copySchema(sqlite3 *db, File5_Tsv *tsv, File5_Sqlite_Config &config) {
  string createSql;
  makeCreateStatement(tsv, config.primaryColumn, config.tableName, config.numeric, createSql);
  Verbose::out(2, "Creating table: " + config.tableName + " from file: " + tsv->m_name);
  execute(db, createSql);
}

void File5_Sqlite::execute(sqlite3 *db, const string &sql) {
  
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

void File5_Sqlite::copyData(sqlite3 *db, File5_Tsv *tsv, File5_Sqlite_Config &config) {

  string tableName = config.tableName;
  int colCount = tsv->getColumnCount(0);
  string insert;
  int count = 0;
  int interval = 0;
  int lineCount = tsv->getLineCount();
  interval = lineCount / 20;
  int dotMod = max(int(lineCount/20),1);
  Verbose::progressBegin(2, "Inserting " + ToStr(lineCount) + " rows", 20, dotMod, lineCount);
  execute(db, "begin transaction;");
  while(tsv->nextLevel(0) == FILE5_OK) {
    count++;
    Verbose::progressStep(2);
    string insert = "insert into " + tableName + " values ( ";
    for(int i = 0; i < colCount; i++) { 
      affx::File5_TsvColumn* column = tsv->getColumnPtr(0,i);
      File5_dtype_t type = column->file5_dtype();
      if (type == FILE5_DTYPE_STRING) {
        string s;
        tsv->get(0, i, &s);
        insert += " '" + s + "' ";
      }
      else if (type == FILE5_DTYPE_CHAR) {
        char c;
        tsv->get(0, i, &c);
        insert += ToStr((int)c);
      }
      else if (type == FILE5_DTYPE_SHORT) {
        Err::errAbort("No File5_Tsv::get() for shorts.");
      }
      else if (type == FILE5_DTYPE_INT) {
        int integer;
        tsv->get(0, i, &integer);
        insert += ToStr(integer);
      }
      else if (type == FILE5_DTYPE_FLOAT) {
        float f;
        tsv->get(0, i, &f);
        insert += ToStr(f);
      }
      else if (type == FILE5_DTYPE_DOUBLE) {
        double d;
        tsv->get(0, i, &d);
        insert += ToStr(d);
      }
      else {
        Err::errAbort("Don't recognize type: " + ToStr(type));
      }
      if ((i + 1) != colCount)
        insert += ", ";
    }
    insert += ");";
    execute(db, insert);
  }
  
  execute(db, "commit;");
  Verbose::progressEnd(2, "Done.");
}


