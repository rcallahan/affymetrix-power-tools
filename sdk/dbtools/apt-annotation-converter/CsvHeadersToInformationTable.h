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
#ifndef CSVHEADERSTOINFORMATIONTABLE_H
#define CSVHEADERSTOINFORMATIONTABLE_H

#include <cstring>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;


#include "HeaderVariableToMatchValidator.h"
#include "Configuration.h"
#include "../../util/Except.h"

typedef std::map<string, map<string, string > > CsvFilesHeadersMap; // map<File Name, map<Header Key, Header Value > >
typedef std::map<string, HeaderVariableToMatchValidator> HeaderValidationMap;  // <Var Name, HeaderVariableToMatchValidator>
#include "sqlite3.h"

class CsvHeadersToInformationTable
{
    public:
    CsvHeadersToInformationTable(sqlite3* AnnotationDb, Configuration* ConfigPtr);
    void LoadInformationTable();
    void UpdateInformationTable();
    void AddVariable(string FileName, string VarName, string VarValue);
    ValidationResult AnalyzeHeaders();
    void Update(CsvFilesHeadersMap HeadersMap);
    
    void (*ValueMismatchReport)(string VarName, map<string, string> FileValueMap);
    void (*MissingValueReport)(string VarName, vector<string> Files);
        
        
    private:
    void UpdateInformationTableMap();
    void BuildValidationMap();
    void DeleteAllFromInformationTable();
    string GetAnnotationTableCount();
    static int AnnotationTableCountCallback(void *NotUsed, int argc, char **argv, char **azColName);
    static string m_annotationTableCount;
    
    sqlite3_stmt * PrepareInsertStatement();
    static int ReadInformationTableCallback(void *NotUsed, int argc, char **argv, char **azColName);
    CsvFilesHeadersMap m_csvFilesHeaders;
    std::vector<string> m_filesOrderKeeper;
    HeaderValidationMap m_headerValidationMap;
    
    static map<string,string> m_informationTableMap; // Loaded from Information table.
    
    Configuration* m_configPtr;
    sqlite3* m_annotationDb;
};

#endif
