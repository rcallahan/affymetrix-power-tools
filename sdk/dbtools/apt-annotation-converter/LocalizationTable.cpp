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

#include "LocalizationTable.h"

#include "StringUtilities.h"


map<string, string > LocalizationTable::m_csvColumnNameToDbColumnNameMap;

LocalizationTable::LocalizationTable(sqlite3* AnnotationDb, Configuration* ConfigPtr)
{
    m_annotationDb = AnnotationDb;
    m_configPtr = ConfigPtr;
}

string LocalizationTable::GetDbColumnName(string CsvColumnName)
{
    string dbColumnName;
    map<string,string>::iterator csvColumnNameToDbColumn = m_csvColumnNameToDbColumnNameMap.find(CsvColumnName);
    if (csvColumnNameToDbColumn != m_csvColumnNameToDbColumnNameMap.end())
    { 
        dbColumnName = csvColumnNameToDbColumn->second;
    }    
    return dbColumnName;
}    
void  LocalizationTable::AppendColumnToAddList(string CsvColumnName)
{    
    std::vector<std::string>::iterator columnToAdd = find(m_listOfRecordsToAdd.begin(), m_listOfRecordsToAdd.end(), CsvColumnName);
    if ( columnToAdd == m_listOfRecordsToAdd.end()) 
    {
        m_listOfRecordsToAdd.push_back(CsvColumnName);
    }     
}

int LocalizationTable::ReadTableCallback(void *NotUsed, int argc, char **argv, char **azColName)
{   

    string csvColumn = argv[1];
    string dbColumn = argv[0];
    if (csvColumn == "Chromosome")
    {// Chromosome as a specisl case
        m_csvColumnNameToDbColumnNameMap.insert(std::make_pair(csvColumn,"Chr_id"));
    }
    else
    {
        m_csvColumnNameToDbColumnNameMap.insert(std::make_pair(csvColumn, dbColumn));    
    }
        
    return 0;
}



void LocalizationTable::Load()
{
    
    /*
        If we ever going to use L10N table for actual translation to local languages, it has to be changes to
        Colname, localName, languageCode.

        Now language code  is the name if the column! 
        "en" is hard coded for now
    */
    char *errorMessage = NULL;
    
    m_csvColumnNameToDbColumnNameMap.clear();
    
      if(sqlite3_exec(m_annotationDb, "Select colname, en from Localization", LocalizationTable::ReadTableCallback, 0, &errorMessage) != SQLITE_OK)
      {
        std::string message = "Cannot read data from the Localization table! ";
        if (errorMessage != NULL)
        {
            message += errorMessage;
            delete(errorMessage);
        }
        throw Except(message.c_str());
      }   
      
      /// Update m_csvColumnNameToDbColumnNameMap with m_alternativeColumnMappings
      map<string, string>::iterator alternMapping;
      map<string,string>::iterator csvColumnNameToDbColumn;
      for (alternMapping = m_configPtr->m_alternativeColumnMappings.begin(); alternMapping != m_configPtr->m_alternativeColumnMappings.end(); alternMapping++)
      {        
        for (csvColumnNameToDbColumn = m_csvColumnNameToDbColumnNameMap.begin(); csvColumnNameToDbColumn != m_csvColumnNameToDbColumnNameMap.end(); csvColumnNameToDbColumn++)
        {
            if (alternMapping->second == csvColumnNameToDbColumn->second) // Same dbColumn
            {
               m_csvColumnNameToDbColumnNameMap.erase(csvColumnNameToDbColumn);                
               break;
            }            
        }
        
        m_csvColumnNameToDbColumnNameMap.insert(std::make_pair(alternMapping->first,alternMapping->second));         
      }
      ///
}

void LocalizationTable::AddRecords()
{    
    sqlite3_stmt *insertStatement;         
    int result = sqlite3_prepare_v2(m_annotationDb, "INSERT INTO Localization (colname, en) VALUES (?,?)", -1, &insertStatement, NULL);
    
    if (result != SQLITE_OK)
    {                     
        throw Except("Internal error. Failed to prepare INSERT INTO Localization statement. ");
    }
    
    std::vector<std::string>::iterator columnsIterator;
    for (columnsIterator = m_listOfRecordsToAdd.begin(); columnsIterator != m_listOfRecordsToAdd.end(); columnsIterator++)
    {
        std::string csvColumnName = *columnsIterator;
        std::string dbColumnName = StringUtilities::MakeColumnName(csvColumnName);
        m_csvColumnNameToDbColumnNameMap.insert(std::make_pair(csvColumnName,dbColumnName));
        
        
        if ((sqlite3_bind_text(insertStatement, 1, dbColumnName.c_str(), dbColumnName.length(), SQLITE_TRANSIENT) != SQLITE_OK) ||
            (sqlite3_bind_text(insertStatement, 2, csvColumnName.c_str(), csvColumnName.length(), SQLITE_TRANSIENT) != SQLITE_OK) 
           )
        {
            std::stringstream message;            
            message << "Could not bind parameters for Localization table update! column name: '" << dbColumnName << "' en: '" << csvColumnName << "'";
            throw Except(message.str().c_str());
        }
                
        
        if (sqlite3_step(insertStatement) != SQLITE_DONE) 
        {
            std::stringstream message;            
            message << "Could not execute SQL statement during Localization table update!  colname: '" << dbColumnName << "' en: '" << csvColumnName << "'";
            throw Except(message.str().c_str());
        }
        sqlite3_reset(insertStatement);
        sqlite3_clear_bindings(insertStatement);
    }
    
    sqlite3_finalize(insertStatement);
}

