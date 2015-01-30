////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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

#ifndef ANNOTATIONTABLE_H
#define ANNOTATIONTABLE_H


#include <map>
#include <string>
#include <algorithm>
using namespace std;

#include "sqlite3.h"

#include "util/BaseEngine.h"

#include "AnnotationTsv.h"
using namespace affx;

#include "Configuration.h"
#include "LocalizationTable.h"
#include "ChromosomeTable.h"

class AnnotationTable
{
    public:
    AnnotationTable(sqlite3* AnnotationDb, Configuration* ConfigPtr, 
                                 LocalizationTable* LocalizationTablePtr, ChromosomeTable* ChromosomeTablePtr);
    void LoadColumntTypes();
    bool InsertTsvFiles(std::vector<std::string> ListOfTsvFiles);
    // http://jira.ev.affymetrix.com:8080/browse/CSVTOSQLITE-12
    void AnalyzeFiles(std::vector<std::string> ListOfTsvFiles); 
    void ReportNumberOfRecordsToInsertUpdate(int TotalLinesCount);
    void ReportRequiredColumns();
    
            
    //Check if column exists and adds it into the column list to add
    void CheckIfNeedToAddColumn(string CsvColumnName);
    void CheckIfNeedToAddColumn(string DbColumnName, string ColumnType);
    
    void AddColumns();
    static map<string, string > m_columnNameToColumnTypeMap;
    
    private:
    sqlite3* m_annotationDb;
    bool m_stopDbInsertion;
    
    Configuration* m_configPtr;
    LocalizationTable* m_localizationTablePtr;
    ChromosomeTable* m_chromosomeTablePtr;
    string GetRowCount(string ProbeSet_Id);
    static int RowCountCallback(void *NotUsed, int argc, char **argv, char **azColName);
    static string m_rowCount;

    std::map<string, string> m_listOfColumnsToAdd; // <DbColumn Name, Column Type>
    std::vector<int> m_identicalLinesToSkip; //  Only for the first file
    map<int, string > m_csvColumnIndexToColumnNameMap;
    std::vector<string> m_csvColumns;
    int m_probeSetIdIndex;
    
    map<int64_t, int> m_linesWithDupProbeIds; // <IntsToLong(FileIndex, FileLine), foo int>
    
    bool m_updateRowMode;
    void PrepareInsertStatement(AnnotationTsv* TsvFilePtr);
    void PrepareUpdateStatement(AnnotationTsv* TsvFilePtr);
    vector<string> GetListOfDbColumnsForTsvFile(AnnotationTsv* TsvFilePtr);
    string GetDbColumnNameFromLocalizationTable(string CsvColumnName);        
    
    
    sqlite3_stmt * m_currentStatement;    
    sqlite3_stmt * m_insertStatement;    
    map<string, int> m_insertBindingMap;
    sqlite3_stmt * m_updateStatement;
    map<string, int> m_updateBindingMap;
    void BuildCsvColumnIndexToColumnNameMap(AnnotationTsv* TsvFilePtr);
    void BuildCsvColumnList(AnnotationTsv* TsvFilePtr);
    void GetCsvColumnName(int ColumnIndex, string &CsvColumnName);
    void BindFieldValue(AnnotationTsv* TsvFilePtr, std::string &DbColumnName, std::string &CsvColumnName, std::string &FieldValue);
    
    void ReportTsvErrorStopDbInsertion(std::string Message, int LineNumber, std::string FieldValue, 
                                             std::string ColumnName, std::string ColumnType);
    void FinalizeStatements();
};

#endif
