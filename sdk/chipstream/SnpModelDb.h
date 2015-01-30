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
 * @file   SnpModelDb.h
 * @author Chuck Sugnet
 * @date   Mon Mar 15 10:28:28 2010
 * 
 * @brief Class for reading and writing snp clustering models (eg
 * gaussian elipses of mean and variance). Uses sqlite as a background
 * to provide convenient fast random access to large numbers of snps
 * (e.g. Millions) with minimal memory impact.
 * 
 */

#ifndef SNP_MODEL_DB_H
#define SNP_MODEL_DB_H

#include "util/SQLite.h"
#include "chipstream/QuantLabelZ.h"
#include <string>

/**
 * @brief Class for reading and writing snp clustering models (eg
 * gaussian elipses of mean and variance). Uses sqlite as a background
 * to provide convenient fast random access to large numbers of snps
 * (e.g. Millions) with minimal memory impact.
 */
class SnpModelDb {

public:
    SnpModelDb(const std::string &fileName);
    SnpModelDb();
    ~SnpModelDb();

    static bool isSqlModelFile(const std::string &file);
    void open(const std::string &fileName);
    void close();

    // Reading tables
    void getMetaValue(const std::string &key, std::vector<std::string> &values) const;
    bool getSnpDistribution(const std::string &probesetName, int copynumber, snp_labeled_distribution &snp) const;

    // Writing tables
    void setupTables();
    void addMeta(const std::string &key, const std::string &value);
    void prepareToWriteModels(int numExpected);
    void writeSnpDistribution(const snp_labeled_distribution &p);
    void doneWriteModels();

private:
    void snpDistFromRSet(SQLiteRecordset &rset, snp_labeled_distribution &snp) const;
    sqlite3_stmt *m_Stmt;
    mutable SQLiteDatabase m_Db;
};

#endif
