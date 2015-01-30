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

#ifndef _FILE5_SQLITE_H_
#define _FILE5_SQLITE_H_

//
#include "file5/File5.h"
//
#include "sqlite3.h"
//
#include <string>

/**
 * Helper class to hold the various options for
 * converting file5 files to sqlite databases.
 */
class File5_Sqlite_Config {

public:

  File5_Sqlite_Config() {
    numeric = "real";
  }

  /// If this column is seen in the file5 table make it the primary column
  std::string primaryColumn;
  /// Name of the database table to be created. Generally file5_tsv name by default.
  std::string tableName;
  /// Precision to use for floating point, default is "real" but could do "numeric(4,5)" as example
  std::string numeric;

};

/**
 * Class with some utility functions for converting File5 files
 * into sqlite databases for analysis. Currently static libraries to
 * be used by other programs.
 */
class File5_Sqlite {
  
public:
  
  /** 
   * Use the column definitions in the File5 definition to create a table in sqlite database
   * 
   * @param db - Open connection to the sqlite database where table will be created.
   * @param tsv5 - File5_Tsv pointer which will be read to create column definitions.
   * @param config - Various configurations to make the table like table name, and primary column name
   */
  static void copySchema(sqlite3 *db, affx::File5_Tsv *tsv5, File5_Sqlite_Config &config);
  

  /** 
   * Fill in an sqlite databse table with the data in the File5_Tsv table
   * 
   * @param db - Open connection to the sqlite database data will be inserted
   * @param tsv5 - File5_Tsv pointer which will be read for data
   * @param config - Various configurations to insert the data like table name etc.
   */
  static void copyData(sqlite3 *db, affx::File5_Tsv *tsv5, File5_Sqlite_Config &config);
  
private:

  /** 
   * Helper function to make an sql statement that will replicated the File5_Tsv table in 
   * sql database.
   * 
   * @param tsv5 - File with data definition.
   * @param primaryColumn - What is the name, if any, of the primary column
   * @param tableName - Name of the table to be created
   * @param numeric - Sql specification for real numbers (eg numeric(4,3))
   * @param createSql - String where sql will be filled in.
   */
  static void makeCreateStatement(affx::File5_Tsv *tsv5, const std::string &primaryColumn,
                                  const std::string &tableName, const std::string &numeric,
                                  std::string &createSql);
  
  /** 
   * Helper function to exectue a particular piece of sql in the sqlite database. Not as efficient
   * as prepared statements as every call prepares a statement, steps through it and finalizes it.
   * Usefule for insert, create table, etc. where no results are expected.
   * 
   * @param db - Open connection to sqlite database.
   * @param sql - Sql to be run by the database.
   */
  static void execute(sqlite3 *db, const std::string &sql);
    
};

#endif /* _FILE5_SQLITE_H_ */
