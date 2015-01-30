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
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS // To disable localtime function warning
#endif

#include "CsvHeadersToInformationTable.h"
#include "util/Verbose.h"
#include "AptVerboseWrapper.h"
#include "StringUtilities.h"

map<string, string> CsvHeadersToInformationTable::m_informationTableMap;
string CsvHeadersToInformationTable::m_annotationTableCount;

CsvHeadersToInformationTable::CsvHeadersToInformationTable(sqlite3* AnnotationDb, Configuration* ConfigPtr)
{
    m_annotationDb = AnnotationDb;
    ValueMismatchReport = NULL;
    MissingValueReport = NULL; 
    m_configPtr = ConfigPtr;
}

int CsvHeadersToInformationTable::ReadInformationTableCallback(void *NotUsed, int argc, char **argv, char **azColName)
{
    m_informationTableMap.insert(std::make_pair(argv[0],argv[1]));
    return 0;
}

string CsvHeadersToInformationTable::GetAnnotationTableCount()
{
    char *errorMessage = NULL;
    if(sqlite3_exec(m_annotationDb, "Select count(*) from Annotations", CsvHeadersToInformationTable::AnnotationTableCountCallback, 0, &errorMessage) != SQLITE_OK)
      {
        std::string message = "Cannot get row count from Annotations table! ";
        if (errorMessage != NULL)
        {
            message += errorMessage;
            delete(errorMessage);
        }
        AptVerboseWrapper::out(0,message, true);
      }
      return m_annotationTableCount;
}
int CsvHeadersToInformationTable::AnnotationTableCountCallback(void *NotUsed, int argc, char **argv, char **azColName)
{    
    m_annotationTableCount = argv[0];
    return 0;
}

void CsvHeadersToInformationTable::UpdateInformationTableMap()
{
    //The simple case: Update with m_varsToMatch. Values are identical for all files    
    HeaderValidationMap::iterator headerValidationMapIterator;
    for (headerValidationMapIterator = m_headerValidationMap.begin(); headerValidationMapIterator != m_headerValidationMap.end(); headerValidationMapIterator++)    
    {
       std::string key = headerValidationMapIterator->first;
       HeaderVariableToMatchValidator requiredHeaderVariableValidator = headerValidationMapIterator->second;
       if (requiredHeaderVariableValidator.GetStatus() != ValueMismatch)
       {
            std::string valueInCsv = requiredHeaderVariableValidator.GetValue();
            if (false == valueInCsv.empty())
            {
                m_informationTableMap[key] = valueInCsv;
            }
       }
    }
    //////
    
    /// Simple Single Values    
    m_informationTableMap["annotation_row_count"] = GetAnnotationTableCount();
    
    time_t now = time(NULL);
    tm *timeinfo = localtime (&now);
    char buf[256];	
	strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M", timeinfo);
	string db_create_date = buf;
	m_informationTableMap["db_create_date"] = db_create_date;	
	/// 
	
	
	
	/// Add other variables from all txt files
	    std::vector<string> allVariableNames;
	    string imported_files;	    
	    for (unsigned int i = 0; i < m_filesOrderKeeper.size(); i++)
	    {//Collect all variable names	       
	        
	        map<string, string> fileVarList = m_csvFilesHeaders[m_filesOrderKeeper[i]];
	        map<string, string>::iterator varKeyPair;
	        for ( varKeyPair = fileVarList.begin(); varKeyPair != fileVarList.end(); varKeyPair++)
	        {
	            string varName = varKeyPair->first;
	            if (false == varName.empty())
	            {                                        
                    std::vector<string>::iterator varToMatch  = find(m_configPtr->m_varsToMatch.begin(), m_configPtr->m_varsToMatch.end(), varName);
                    if (varToMatch == m_configPtr->m_varsToMatch.end())
                    {
                        if (find(allVariableNames.begin(), allVariableNames.end(), varName) == allVariableNames.end())
                        {
                            allVariableNames.push_back(varKeyPair->first);
                        }
                    }                                     
                }	            
	        }
	        
	        if (imported_files != "")
	        {
	            imported_files += ", ";
	        }
	        imported_files += "'";
	        imported_files += m_filesOrderKeeper[i];
	        imported_files += "'";
	    }
	    
	    m_informationTableMap["imported-files"] = imported_files;
	    
	    for (unsigned int n = 0; n < allVariableNames.size(); n++)
	    {
	        string varName = allVariableNames[n];
	        std::vector<string> varValues;
	        for (unsigned int i = 0; i < m_filesOrderKeeper.size(); i++)
	        {	           
	            map<string, string> fileVarList = m_csvFilesHeaders[m_filesOrderKeeper[i]];
	            map<string, string>::iterator varKeyValue = fileVarList.find(varName);
	            if (varKeyValue != fileVarList.end())
	            {
	                varValues.push_back(varKeyValue->second); // varValue	                
	            }
	            else
	            {
	                varValues.push_back("");
	            }	           	            
	        }
	        bool allValuesIdentical = true;
	        for (unsigned int v = 1; v < varValues.size(); v++)
	        {
	            if (varValues[v-1] != varValues[v])
	            {
	                allValuesIdentical = false;
	                break;
	            }
	        }
	        if (allValuesIdentical)
	        {// Do not use comma separated array is all values are the same
	            m_informationTableMap[varName] = varValues[0];
	        }
	        else
	        {
	            string combinedValues;
	            for (unsigned int v = 0; v < varValues.size(); v++)
	            {
	                if (combinedValues != "")
                    {
                        combinedValues += ", ";
                    }
                    combinedValues += "'";
                    combinedValues += varValues[v];
                    combinedValues += "'";
	            }
	            m_informationTableMap[varName] = combinedValues;
	        }
	    }
	    
}
/* Since table is small, do lazy update: delete all and reinsert
*/
void CsvHeadersToInformationTable::UpdateInformationTable()
{

    UpdateInformationTableMap();
    
    DeleteAllFromInformationTable();
    
    sqlite3_stmt *insertStatement = PrepareInsertStatement();
    
    //char *errorMessage = NULL;
    static map<string,string>::iterator informationTableMapIterator;
    
    for (informationTableMapIterator = m_informationTableMap.begin(); informationTableMapIterator != m_informationTableMap.end(); informationTableMapIterator++)
    {
        std::string key = informationTableMapIterator->first;
        std::string value = informationTableMapIterator->second;
        
        if ((sqlite3_bind_text(insertStatement, 1, key.c_str(), key.length(), SQLITE_TRANSIENT) != SQLITE_OK) ||
            (sqlite3_bind_text(insertStatement, 2, value.c_str(), value.length(), SQLITE_TRANSIENT) != SQLITE_OK) 
           )
        {
            std::stringstream message;            
            message << "Could not bind parameters during Information table update! Key: '" << key << "' Value: '" << value << "'";
            throw Except(message.str().c_str());
        }
                
        
        if (sqlite3_step(insertStatement) != SQLITE_DONE) 
        {
            std::stringstream message;            
            message << "Could not execute SQL statement during Information table update! Key: '" << key << "' Value: '" << value << "'";
            throw Except(message.str().c_str());
        }
        sqlite3_reset(insertStatement);            
        sqlite3_clear_bindings(insertStatement);
    }
    sqlite3_finalize(insertStatement);
}
void CsvHeadersToInformationTable::DeleteAllFromInformationTable()
{
    char *errorMessage = NULL;
    if(sqlite3_exec(m_annotationDb, "Delete from Information", NULL, 0, &errorMessage) != SQLITE_OK)
    {
        
        std::string message = " Cannot delete records from the Information table! ";
        if (errorMessage != NULL)
        {
            message += errorMessage;
            delete(errorMessage);
        }
        throw Except(message.c_str());
    }
}

sqlite3_stmt * CsvHeadersToInformationTable::PrepareInsertStatement()
{
    sqlite3_stmt *compiledStatement;         
    int result = sqlite3_prepare_v2(m_annotationDb, "INSERT INTO Information (key, value) VALUES (?,?)", -1, &compiledStatement, NULL);
    
    if (result != SQLITE_OK)
    {                     
        throw Except("Internal error. Failed to prepare INSERT INTO Information statement. ");
    }

    return compiledStatement;
            
}

void CsvHeadersToInformationTable::LoadInformationTable()
{
    char *errorMessage = NULL;
    
    m_informationTableMap.clear();
    
      if(sqlite3_exec(m_annotationDb, "select key, value from Information", ReadInformationTableCallback, 0, &errorMessage) != SQLITE_OK)
      {
        std::string message = "Cannot read data from the Localization table! ";
        if (errorMessage != NULL)
        {
            message += errorMessage;
            delete(errorMessage);
        }
        throw Except(message.c_str());
      }
}
void CsvHeadersToInformationTable::AddVariable(string FileName, string VarName, string VarValue)
{    
    StringUtilities::FixHeaderVarName(VarName);
    
    CsvFilesHeadersMap::iterator fileIterator = m_csvFilesHeaders.find(FileName);
    if ( fileIterator == m_csvFilesHeaders.end())
    {
        m_filesOrderKeeper.push_back(FileName);
        map<string,string> csvHeaderLine;
        csvHeaderLine.insert(make_pair(VarName, VarValue));       
        m_csvFilesHeaders.insert(make_pair(FileName, csvHeaderLine));		                 
    }
    else
    {
        m_csvFilesHeaders[FileName].insert(make_pair(VarName, VarValue));
    }
    
}

void CsvHeadersToInformationTable::BuildValidationMap()
{
    std::vector<string>::iterator varsToMatchIterator;
    for ( varsToMatchIterator = m_configPtr->m_varsToMatch.begin(); varsToMatchIterator != m_configPtr->m_varsToMatch.end(); varsToMatchIterator++)
    {
        string varName = *varsToMatchIterator;        
        CsvFilesHeadersMap::iterator fileIterator;
        for ( fileIterator = m_csvFilesHeaders.begin(); fileIterator != m_csvFilesHeaders.end(); fileIterator++)
        {
            string fileName = fileIterator->first;
            map<string, string> header = fileIterator->second;
            map<string, string>::iterator headerIterator = header.find(varName);
            
            string varValue;
            if (headerIterator != header.end())
            {
                varValue = headerIterator->second;
            }                                    
            m_headerValidationMap[varName].AddVar(fileName, varValue);   
        }                
    }
}

ValidationResult CsvHeadersToInformationTable::AnalyzeHeaders()
{
    ValidationResult result = Ok;
    
    BuildValidationMap();

   CsvFilesHeadersMap::iterator file;
   for ( file = m_csvFilesHeaders.begin(); file != m_csvFilesHeaders.end(); file++)
   {  
        map<string, string> fileVarList = file->second;  
        map<string, string>::iterator varKeyPair;
        for ( varKeyPair = fileVarList.begin(); varKeyPair != fileVarList.end(); varKeyPair++)
        {
            string varName = varKeyPair->first;
            if (false == varName.empty())
            {
                std::vector<string>::iterator reservedVar = find(m_configPtr->m_reservedVars.begin(), m_configPtr->m_reservedVars.end(), varName);
                if (reservedVar != m_configPtr->m_reservedVars.end())
                { 
                    std::string message = "Variable '";
                    message += varName;
                    message += "' is reserved for internal use and cannot be used in the file header!";                
                    throw Except(message.c_str());
                }
            }	            
        }
    }
	        
    HeaderValidationMap::iterator validationMapIterator;
    for (validationMapIterator = m_headerValidationMap.begin(); validationMapIterator != m_headerValidationMap.end(); validationMapIterator++)
    {
        ValidationResult singleVarResult = validationMapIterator->second.Validate();
        if (singleVarResult  == ValueMismatch)
        {
            result = ValueMismatch; // Error
        }
        else if (singleVarResult  == MissingValue && result != ValueMismatch)
        {
            result = MissingValue; // Warning
        }
        
        validationMapIterator->second.CopyValueToFilesWithEmptyValues();
    }
    
    if (result == ValueMismatch && ValueMismatchReport != NULL)
    {
        for (validationMapIterator = m_headerValidationMap.begin(); validationMapIterator != m_headerValidationMap.end(); validationMapIterator++)
        {
            if (validationMapIterator->second.GetStatus() ==  ValueMismatch)
            {
                ValueMismatchReport(validationMapIterator->first, validationMapIterator->second.GetFileValueMap());
            }              
        }
    }
    
    if (result == MissingValue && MissingValueReport != NULL)
    {
        for (validationMapIterator = m_headerValidationMap.begin(); validationMapIterator != m_headerValidationMap.end(); validationMapIterator++)
        {
            if (validationMapIterator->second.GetStatus() ==  MissingValue)
            {
                MissingValueReport(validationMapIterator->first, validationMapIterator->second.GetFilesWithEmptyValues());
            }              
        }
    }
    
    return result;
}

void CsvHeadersToInformationTable::Update(CsvFilesHeadersMap HeadersMap)
{

}
