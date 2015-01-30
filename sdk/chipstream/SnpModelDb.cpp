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

#include "util/Err.h"
#include "util/Verbose.h"
#include "chipstream/SnpModelDb.h"

using namespace std;

#define NUM_DOTS 20

SnpModelDb::SnpModelDb() {
    m_Stmt = NULL;
}

SnpModelDb::SnpModelDb(const std::string &fileName) {
    m_Stmt = NULL;
    open(fileName);
}

SnpModelDb::~SnpModelDb() {
    m_Db.close();
    if (m_Stmt != NULL) {
        sqlite3_finalize(m_Stmt);
    }
}

void SnpModelDb::open(const std::string &fileName) {
    if (m_Stmt != NULL) {
        sqlite3_finalize(m_Stmt);
    }
    m_Db.close();
    m_Db.open(fileName);
}

void SnpModelDb::close() {
    m_Db.close();
}

bool SnpModelDb::isSqlModelFile(const std::string &file) {
    bool isSqlModel = false;
    try {
        SQLiteDatabase db;
        db.open(file);
        SQLiteRecordset rset(db);
        string sql = "select count(*) from snp_model;";
        rset.open(sql);
        isSqlModel = true;
        rset.close();
        db.close();
        isSqlModel = true;
    }
    //catch (SQLiteException &e) {    
    catch (SQLiteException &) {
        // Nothing to do here we've already guessed false
    }
    return isSqlModel;
}

void SnpModelDb::getMetaValue(const std::string &key, std::vector<std::string> &values) const {
    values.clear();
    SQLiteRecordset rset(m_Db);
    string sql = "select value from snp_model_meta where key = '" + key + "';";
    rset.open(sql);
    while (rset.fetch()) {
        values.push_back(rset.getString(0));
    }
    rset.close();
}

void SnpModelDb::snpDistFromRSet(SQLiteRecordset &rset, snp_labeled_distribution &snp) const {
    int index = 0;
    snp.probeset_id = rset.getString(index++);
    snp.Dist.aa.m = rset.getDouble(index++);
    snp.Dist.aa.k = rset.getDouble(index++);
    snp.Dist.aa.ss = rset.getDouble(index++); 
    snp.Dist.aa.v = rset.getDouble(index++);
    snp.Dist.aa.ym = rset.getDouble(index++);
    snp.Dist.aa.yss = rset.getDouble(index++);
    snp.Dist.aa.xyss = rset.getDouble(index++);

    snp.Dist.ab.m = rset.getDouble(index++);
    snp.Dist.ab.k = rset.getDouble(index++);
    snp.Dist.ab.ss = rset.getDouble(index++); 
    snp.Dist.ab.v = rset.getDouble(index++);
    snp.Dist.ab.ym = rset.getDouble(index++);
    snp.Dist.ab.yss = rset.getDouble(index++);
    snp.Dist.ab.xyss = rset.getDouble(index++);

    snp.Dist.bb.m = rset.getDouble(index++);
    snp.Dist.bb.k = rset.getDouble(index++);
    snp.Dist.bb.ss = rset.getDouble(index++); 
    snp.Dist.bb.v = rset.getDouble(index++);
    snp.Dist.bb.ym = rset.getDouble(index++);
    snp.Dist.bb.yss = rset.getDouble(index++);
    snp.Dist.bb.xyss = rset.getDouble(index++);
}

bool SnpModelDb::getSnpDistribution(const std::string &probesetName, int copynumber, snp_labeled_distribution &snp) const {
    SQLiteRecordset rset(m_Db);
    snp.probeset_id.clear();
    snp.Dist.Clear();
    string key = copynumber == 1 ? probesetName + ":" + ToStr(copynumber) : probesetName;
    // Better be in right order...
    string sql = "select * from snp_model where probeset_name = '" + key + "'";
    rset.open(sql);
    bool found = false;
    while (rset.fetch()) {
        if (found) {
            Err::errAbort("Multiple models for snp: " + probesetName + " copynumber: " + ToStr(copynumber));
        }
        snpDistFromRSet(rset, snp);
        found = true;
    }
    return found;
}

void SnpModelDb::addMeta(const std::string &key, const std::string &value) {

    std::string escapedValue = value;
	std::string escapedKey = key;
    std::string quote("'");
    std::string twoQuotes("''");
    Util::replaceString(escapedValue, quote, twoQuotes);
	Util::replaceString(escapedKey, quote, twoQuotes);
	string sql = "insert into snp_model_meta values ('" + escapedKey + "','" + escapedValue + "');";
    try {
        m_Db.execute(sql);
    }
    catch (SQLiteException &e) {
        APT_ERR_ABORT("Got Exception at inserting meta - sql: " + sql + ". message: " + e.getMessage());
    }
}

void SnpModelDb::prepareToWriteModels(int numExpected) {
    unsigned int dotMod = Max(int(numExpected/NUM_DOTS), 1);
    try {
        m_Db.execute("begin;");
        Verbose::progressBegin(1, "Preparing models", NUM_DOTS, dotMod, numExpected);
        const char *sql = "insert into snp_model values ( :id, :aamean , :aak, :aavar, :aavark, :aaym, :aayvar, :aacov, "
            " :abmean , :abk, :abvar, :abvark, :abym, :abyvar, :abcov,  " 
            " :bbmean , :bbk, :bbvar, :bbvark, :bbym, :bbyvar, :bbcov  );";
        int sqlMsg = sqlite3_prepare_v2(&m_Db.getConnection(), sql, -1, &m_Stmt, NULL);
        if (m_Stmt == NULL || sqlMsg != SQLITE_OK) {
            APT_ERR_ABORT("Couldn't create prepared statement (code: " + ToStr(sqlMsg) + ")");
        }
    }
    catch (SQLiteException &e) {
        APT_ERR_ABORT("Got Exception at preparing statement. message: " + e.getMessage());
    }
}

void SnpModelDb::writeSnpDistribution(const snp_labeled_distribution &p) {
    int index = 1;
    try {
        Verbose::progressStep(1);
        sqlite3_reset(m_Stmt);
        sqlite3_bind_text(m_Stmt, index++, p.probeset_id.c_str(), p.probeset_id.length(), SQLITE_TRANSIENT);
        
        sqlite3_bind_double(m_Stmt, index++, p.Dist.aa.m);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.aa.k);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.aa.ss); 
        sqlite3_bind_double(m_Stmt, index++, p.Dist.aa.v);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.aa.ym);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.aa.yss);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.aa.xyss);
        
        sqlite3_bind_double(m_Stmt, index++, p.Dist.ab.m);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.ab.k);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.ab.ss); 
        sqlite3_bind_double(m_Stmt, index++, p.Dist.ab.v);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.ab.ym);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.ab.yss);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.ab.xyss);
        
        sqlite3_bind_double(m_Stmt, index++, p.Dist.bb.m);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.bb.k);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.bb.ss); 
        sqlite3_bind_double(m_Stmt, index++, p.Dist.bb.v);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.bb.ym);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.bb.yss);
        sqlite3_bind_double(m_Stmt, index++, p.Dist.bb.xyss);
        int sqlMsg = sqlite3_step(m_Stmt);
        APT_ERR_ASSERT( sqlMsg == SQLITE_DONE, "Error executing prepared statement (code " + ToStr(sqlMsg) +")");
    }
    catch (SQLiteException &e) {
        APT_ERR_ABORT("Got Exception at inserting snp: " + p.probeset_id + " message: " + e.getMessage());
    }
}

void SnpModelDb::doneWriteModels() {
    try {
        int sqlMsg = sqlite3_finalize(m_Stmt);
        APT_ERR_ASSERT( sqlMsg == SQLITE_OK, "Error finalizing statement (code " + ToStr(sqlMsg) +")");
        m_Stmt = NULL;
        m_Db.execute("end;");
        string createIndex = "create index snp_models_idx on snp_model (probeset_name);";
        m_Db.execute(createIndex);
        Verbose::progressEnd(1, "Done.");
    }
    catch (SQLiteException &e) {
        APT_ERR_ABORT("Got Exception at preparing statement. message: " + e.getMessage());
    }
}

void SnpModelDb::setupTables() {
    try {
        string dropExisting = "drop table if exists snp_model_meta;";
        m_Db.execute(dropExisting);
        dropExisting = "drop table if exists snp_model;";
        m_Db.execute(dropExisting);
        string createMetaTable = 
            "create table snp_model_meta ( "
            "key text, "
            "value text "
            ");";
        m_Db.execute(createMetaTable);
        string createIndex = "create index snp_models_meta_idx on snp_model_meta (key);";
        m_Db.execute(createIndex);
        
        string createDataTable = 
            "create table snp_model ( "
            "probeset_name text, " 
            "aa_x_mean numeric, "
            "aa_x_mean_k numeric, "
            "aa_x_var  numeric, "
            "aa_x_var_k numeric, "
            "aa_y_mean numeric, "
            "aa_y_var numeric, "
            "aa_xy_covar numeric, "
            "ab_x_mean numeric, "
            "ab_x_mean_k numeric, "
            "ab_x_var  numeric, "
            "ab_x_var_k numeric, "
            "ab_y_mean numeric, "
            "ab_y_var numeric, "
            "ab_xy_covar numeric, "
            "bb_x_mean numeric, "
            "bb_x_mean_k numeric, "
            "bb_x_var  numeric, "
            "bb_x_var_k numeric, "
            "bb_y_mean numeric, "
            "bb_y_var numeric, "
            "bb_xy_covar numeric "
            ");";
        m_Db.execute(createDataTable);
    }
    catch (SQLiteException &e) {
        APT_ERR_ABORT("Got Exception at creating databases. message: " + e.getMessage());
    }
}
