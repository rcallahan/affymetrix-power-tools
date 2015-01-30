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


#include <iostream>
#include <string>
#include <vector>

#include "sqlite3.h"

#include "util/AptVersionInfo.h"
#include "util/Err.h"
#include "util/PgOptions.h"
#include "util/Verbose.h"
#include "file5/File5.h"
#include "file5/File5_Sqlite.h"

using namespace std;
using namespace affx;

void defineDataSubsetOptions(PgOptions* opts) {
  opts->setUsage("apt-file5-to-sqlite - Create and populate a sqlite database table from a file5 file."
                 "usage:\n"
                 "   apt-file5-to-sqlite\n");

  opts->defineOption("h", "help", PgOpt::BOOL_OPT,
                     "Display program options and extra documentation about usage.",
                     "false");
  opts->defineOption("v", "verbose", PgOpt::INT_OPT,
                     "How verbose to be with status messages 0 - quiet, 1 - usual messages, 2 - more messages.",
                     "1");
  opts->defineOption("", "version", PgOpt::BOOL_OPT,
                     "Display version information.",
                     "false");

  opts->defineOption("", "file5-file", PgOpt::STRING_OPT,
                     "File5 format file to be read from.",
                     "");
  opts->defineOption("", "file5-tsv", PgOpt::STRING_OPT,
                     "File5_Tsv file with parent file5 file to create database table with.",
                     "");
  opts->defineOption("", "database", PgOpt::STRING_OPT,
                     "Sqlite3 database file or name of new file where new table will be created.",
                     "");
  opts->defineOption("", "tablename", PgOpt::STRING_OPT,
                     "Name of the database table to be inserted.",
                     "");
  opts->defineOption("", "primary-column", PgOpt::STRING_OPT,
                     "Name of column in file5 to be primary key.",
                     "");
  opts->defineOption("", "numeric-def", PgOpt::STRING_OPT,
                     "Sql provides a number of options for floating point numbers (eg numeric(4,2)). Default is 'real'",
                     "real");

}

/* Everybody's favorite function. */
int main(int argc, char *argv[]) {
  const string version = AptVersionInfo::versionToReport();

  // Parse options
  PgOptions *opts = NULL;
  opts = new PgOptions();
  defineDataSubsetOptions(opts);
  opts->parseArgv(argv);
  Verbose::setLevel(opts->getInt("verbose"));

  // Check for help or version
  if(opts->getBool("help") || argc == 1) {
    opts->usage();
    cout << "version: " << version << endl;
    exit(0);
  }
  else if(opts->getBool("version")) {
    cout << "version: " << version << endl;
    exit(0);
  }
  
  // Do the file5 import into sql
  File5_Sqlite_Config config;
  string file5File = opts->get("file5-file");
  string file5Tsv = opts->get("file5-tsv");
  string database = opts->get("database");
  config.tableName = opts->get("tablename");
  config.primaryColumn = opts->get("primary-column");
  config.numeric = opts->get("numeric-def");
  sqlite3 *db;
  int status = 0;

  // Open the file5 data table
  File5_File file5;
  File5_Tsv *tsv5 = NULL;
  if ((status= file5.open(file5File, affx::FILE5_OPEN_RO)) != FILE5_OK) {
    Err::errAbort("Couldn't open file: " + file5File);
  }
  tsv5 = file5.openTsv(file5Tsv, affx::FILE5_OPEN_RO);
  if(tsv5 == NULL) {
    Err::errAbort("Couldn't open file5_tsv: " + file5Tsv + " in file5: " + file5File);
  }
  
  // Open the sqlite database connection
  if ((status =  sqlite3_open(database.c_str(), &db)) != SQLITE_OK) {
    Err::errAbort("Couldn't open database: " + database + " got error code: " + ToStr(status));
  }

  // Create the table in the database
  File5_Sqlite::copySchema(db, tsv5, config);

  // Fill in the data
  File5_Sqlite::copyData(db, tsv5, config);

  // Close the sqlite database 
  if ((status = sqlite3_close(db)) != SQLITE_OK) {
    Err::errAbort("Error trying to close the database with error code: " + ToStr(status));
  }
  Verbose::out(1, "Closed database");
  
  // Cleanup the file5 data tables
  tsv5->close();
  delete tsv5;
  file5.close();
  delete opts;
  return 0;
}
