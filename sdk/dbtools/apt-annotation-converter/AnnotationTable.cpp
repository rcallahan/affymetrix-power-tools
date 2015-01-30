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

#include <algorithm>
#include <locale>
#include "util/Util.h"
#include "util/Fs.h"
#include "util/Verbose.h"


#include "AnnotationTable.h"
#include "ResourceStrings.h"
#include "AptVerboseWrapper.h"
#include "StringUtilities.h"

static const char* PROBE_SET_ID_IN_CSV = "Probe Set ID";


map<string, string > AnnotationTable::m_columnNameToColumnTypeMap;
string AnnotationTable::m_rowCount;

    
AnnotationTable::AnnotationTable(sqlite3* AnnotationDb, Configuration* ConfigPtr, 
                                 LocalizationTable* LocalizationTablePtr, ChromosomeTable* ChromosomeTablePtr)
{
    m_annotationDb = AnnotationDb;
    m_configPtr = ConfigPtr;
    m_localizationTablePtr = LocalizationTablePtr;
    m_chromosomeTablePtr = ChromosomeTablePtr;
    m_stopDbInsertion = false;
}
void AnnotationTable::LoadColumntTypes()
{
    sqlite3_stmt *compiledStatement;         
    string sqlStatement = "Select * from Annotations";
    int result = sqlite3_prepare_v2(m_annotationDb, sqlStatement.c_str(), -1, &compiledStatement, NULL);
    if (result != SQLITE_OK)
    {    
        std::string message = INTERNAL_ERROR_FAILED_TO_PREPARE_SQL_STATEMENT + sqlStatement;
        throw Except(message.c_str());
    }
    
    m_columnNameToColumnTypeMap.clear();
    for (int columnIndex = 0;  columnIndex < sqlite3_column_count(compiledStatement); columnIndex++)
    {
       std::string columnName = sqlite3_column_name (compiledStatement, columnIndex);
       std::string  columnType = sqlite3_column_decltype (compiledStatement, columnIndex);
       Util::upcaseString_inplace(columnType);
       //std::transform(columnType.begin(), columnType.end(), columnType.begin(), toupper);
       
       //See http://www.sqlite.org/datatype3.html       
       
       if (columnType.find("INT") != std::string::npos)
       {
            columnType = DB_TYPE_INTEGER;
       }
       else if ((columnType.find("REAL") != std::string::npos) || (columnType.find("FLOA") != std::string::npos)
                || (columnType.find("DOUB") != std::string::npos))
       {
            columnType = DB_TYPE_REAL;
       }
       else
       {
            columnType = DB_TYPE_TEXT;
       }
       m_columnNameToColumnTypeMap.insert(std::make_pair(columnName, columnType));
    }
    
    sqlite3_finalize(compiledStatement);
    
}

vector<string> AnnotationTable::GetListOfDbColumnsForTsvFile(AnnotationTsv* TsvFilePtr)
{
    vector<string> result;
    for(int columnIndex = 0; columnIndex < TsvFilePtr->GetColumnCount(); columnIndex++)
    {           
        string csvColumnName;
        GetCsvColumnName(columnIndex, csvColumnName);
        if (m_configPtr->IsCloneableColumn(csvColumnName))
        {
            std::vector<string> columnClones = m_configPtr->GetColumnClones(csvColumnName);            
            for (unsigned int i = 0; i < columnClones.size(); i++)
            {
               result.push_back(columnClones[i]);
            }
        }
        else if (csvColumnName == m_configPtr->m_par1CsvColumnName || csvColumnName == m_configPtr->m_par2CsvColumnName)
        {
            result.push_back(m_configPtr->m_parDbColumnName);
        }
        else
        {// "Normal" column                            
            result.push_back(GetDbColumnNameFromLocalizationTable(csvColumnName));
        }                        
    }
    return result;
}



void AnnotationTable::PrepareUpdateStatement(AnnotationTsv* TsvFilePtr)
{
    std::stringstream  commandText;
    commandText << "UPDATE Annotations SET ";    

    m_updateBindingMap.clear();
    int columnCounter = 1;
    
    vector<string> listOfDbColumns = GetListOfDbColumnsForTsvFile(TsvFilePtr);   
    for(unsigned int columnIndex = 0; columnIndex < listOfDbColumns.size(); columnIndex++)
    {
         string dbColumnName = listOfDbColumns[columnIndex];
                    
        if (dbColumnName != "ProbeSet_ID")
        {
            m_updateBindingMap.insert(std::make_pair(dbColumnName, columnCounter));
            
            if (columnCounter > 1)
            {
                commandText << ", ";            
            }
            
            commandText << "[" << dbColumnName << "] = " << "?";                        
            columnCounter++;
        }
    }

    commandText << " WHERE ProbeSet_ID = ?";
    m_updateBindingMap.insert(std::make_pair("ProbeSet_ID", columnCounter));
        
    if (m_updateStatement != NULL)
    {
        sqlite3_finalize(m_updateStatement);
    }
    int result = sqlite3_prepare_v2(m_annotationDb, commandText.str().c_str(), -1, &m_updateStatement, NULL);
    
    if (result != SQLITE_OK)
    {                     
        std::string message = INTERNAL_ERROR_FAILED_TO_PREPARE_SQL_STATEMENT + commandText.str();
        throw Except(message.c_str());
    }
}

void AnnotationTable::PrepareInsertStatement(AnnotationTsv* TsvFilePtr)
{
    std::stringstream  commandTextColumns;
    commandTextColumns << "INSERT INTO Annotations (";
    
    std::stringstream  commandTextValues;
    commandTextValues << " VALUES (";

    m_insertBindingMap.clear();
    vector<string> listOfDbColumns = GetListOfDbColumnsForTsvFile(TsvFilePtr);   
    for(unsigned int columnIndex = 0; columnIndex < listOfDbColumns.size(); columnIndex++)
    {   
                
        m_insertBindingMap.insert(std::make_pair(listOfDbColumns[columnIndex], columnIndex + 1));
        
        if (columnIndex > 0)
        {
            commandTextColumns << ",";
            commandTextValues << ",";
        }
                
        commandTextColumns << "[" << listOfDbColumns[columnIndex] << "]";        
        commandTextValues << "?";                
    }   
     
    commandTextColumns << ")";
    commandTextValues << ")";         
     
    std::string insertSqlStatement =  commandTextColumns.str() + commandTextValues.str();
    
    if (m_insertStatement != NULL)
    {
        sqlite3_finalize(m_insertStatement);
    }
    
    int result = sqlite3_prepare_v2(m_annotationDb, insertSqlStatement.c_str(), -1, &m_insertStatement, NULL);
    
    if (result != SQLITE_OK)
    {                     
        std::string message = INTERNAL_ERROR_FAILED_TO_PREPARE_SQL_STATEMENT + insertSqlStatement;
        throw Except(message.c_str());
    }
}

void AnnotationTable::CheckIfNeedToAddColumn(string DbColumnName, string ColumnType)
{    
    map<string,string>::iterator columnNameColumnType = m_columnNameToColumnTypeMap.find(DbColumnName);
    if ( columnNameColumnType == m_columnNameToColumnTypeMap.end()) 
    {        
        std::map<string, string>::iterator columnToAdd = m_listOfColumnsToAdd.find(DbColumnName);
        if ( columnToAdd == m_listOfColumnsToAdd.end()) 
        {            
            m_listOfColumnsToAdd.insert(std::make_pair(DbColumnName, ColumnType));
        } 
    } 
}
void AnnotationTable::CheckIfNeedToAddColumn(string CsvColumnName)
{
    string dbColumnName =  m_localizationTablePtr->GetDbColumnName(CsvColumnName);    
    if (dbColumnName.empty())
    { 
       dbColumnName = StringUtilities::MakeColumnName(CsvColumnName);
    }
         
    map<string,string>::iterator columnNameColumnType = m_columnNameToColumnTypeMap.find(dbColumnName);
    if ( columnNameColumnType == m_columnNameToColumnTypeMap.end()) 
    {        
        std::map<string, string>::iterator columnToAdd = m_listOfColumnsToAdd.find(dbColumnName);
        if ( columnToAdd == m_listOfColumnsToAdd.end()) 
        {
            string columnType = DB_TYPE_TEXT;
            
            map<string, string>::iterator columnTypes = m_configPtr->m_columnTypes.find(CsvColumnName);        
            if (columnTypes != m_configPtr->m_columnTypes.end())
            {
                columnType = columnTypes->second; // Column type from config
            }
            m_listOfColumnsToAdd.insert(std::make_pair(dbColumnName, columnType));
        } 
    } 
}

void AnnotationTable::AddColumns()
{    
    
    char *errorMessage = NULL;
    
    std::map<string, string>::iterator columnsIterator;
    for (columnsIterator = m_listOfColumnsToAdd.begin(); columnsIterator != m_listOfColumnsToAdd.end(); columnsIterator++)
    {
        std::string dbColumnName = columnsIterator->first;
        std::string columnType = columnsIterator->second;
        
        m_columnNameToColumnTypeMap.insert(std::make_pair(dbColumnName, columnType)); 
                
        std::stringstream sqlStatement; 
        sqlStatement << "ALTER TABLE Annotations ADD [" << dbColumnName << "] " << columnType;
        
        if(sqlite3_exec(m_annotationDb, sqlStatement.str().c_str(), NULL, 0, &errorMessage))
        {
            std::stringstream message;
            message << "Cannot add column '" << dbColumnName << "' ! ";
            if (errorMessage != NULL)
            {
                message << errorMessage;
                delete(errorMessage);
            }
            throw Except(message.str().c_str());
        }
    }
    
    AptVerboseWrapper::out(2, LIST_OF_COLUMNS_FROM_TSV_FILES, true);
    for (columnsIterator = m_listOfColumnsToAdd.begin(); columnsIterator != m_listOfColumnsToAdd.end(); columnsIterator++)
    {
        AptVerboseWrapper::out(2, columnsIterator->first, true);
    }            
}




int64_t IntsToLong(int High, int Low)
{
    return (((int64_t)High) << 32) | ((unsigned int)Low);
}
void LongToInts(int64_t Long, int &High, int &Low)
{
    High = (int)(Long >> 32);
    Low = (int)Long;
}

void AnnotationTable::AnalyzeFiles(std::vector<std::string> ListOfTsvFileNames)
{
    if (ListOfTsvFileNames.size() == 0) 
    {
        return;
    }
    else if (ListOfTsvFileNames.size() == 1)
    {
        AptVerboseWrapper::out(1, ANALYZING_FILE, true); 
    }
    else
    {
        AptVerboseWrapper::out(1, ANALYZING_FILES, true);
    }    
    
    
    std::map<string, int> probeSetIds; //  <probeSetId, LineNumber>
    m_identicalLinesToSkip.clear();
    m_linesWithDupProbeIds.clear();
    
    int probeSetIdDuplicationCounter = 0;
    const int limitOfErrorsInLogFile = 10;
    static const char* DUPLICATED_PROBSET_ID_ERROR = "File contains duplicated ProbeSet_IDs for records with different values. See the log file for details";
        
    int totalLinesCount = 0;
    bool needToCheckRequiredColumnsSomeValue = true;
    bool needToCheckRequiredColumnsPrefixCheck = true;
    
    for (unsigned int fileIndex = 0; fileIndex < ListOfTsvFileNames.size(); fileIndex++)
    {
        
        std::string fileName = ListOfTsvFileNames[fileIndex];
        AnnotationTsv* tsvFile = new AnnotationTsv();        
        
        if(tsvFile->open(fileName) != TSV_OK)
        {            
            delete(tsvFile);   
            
            std::string message = "Unable to open file '" + fileName + "'";
            throw Except(message.c_str());		
        }
        else
        {
            int64_t fileSize = Fs::fileSize(fileName, false);
            std::string scale = "1% ----------- 50% ----------- 100%";
            int numberOfSteps = scale.size();
            int stepSize = (int)(fileSize / numberOfSteps);
            int step = 0;
            int lastStep = 0;             

            AptVerboseWrapper::progressBegin(1, "Analyzing file " + Fs::basename(fileName) + "\n", numberOfSteps, 1, numberOfSteps);
            AptVerboseWrapper::out(1, scale + "\n", true);
            
            BuildCsvColumnList(tsvFile);            
            
            while (tsvFile->nextLine() == TSV_OK)
            {  
                totalLinesCount++;
                              
                /// Analyze required columns
                    
                    if (needToCheckRequiredColumnsSomeValue && m_configPtr->m_requiredColumnsSomeValue.size() > 0)
                    {   
                        bool allRequiredColumnsSomeValueChecked = true;
                        
                        std::map<string, RequiredColumn>::iterator requiredColunm;                        
                        for (requiredColunm = m_configPtr->m_requiredColumnsSomeValue.begin(); requiredColunm != m_configPtr->m_requiredColumnsSomeValue.end(); requiredColunm++)
                        {
                            if (false == requiredColunm->second.m_verified)                            
                            {
                                allRequiredColumnsSomeValueChecked = false;
                                string columnName = requiredColunm->first;
                                std::vector<std::string>::iterator column = find(m_csvColumns.begin(), m_csvColumns.end(), columnName);
                                if ( column != m_csvColumns.end()) 
                                {// Column exist in this file
                                    string fieldValue;
                                    tsvFile->get(0, columnName, fieldValue);
                                    if (fieldValue != "")
                                    {
                                        requiredColunm->second.m_verified = true;
                                    }
                                }                                                                 
                            }
                        }
                        if (allRequiredColumnsSomeValueChecked)
                        {
                            needToCheckRequiredColumnsSomeValue = false;
                        }
                    }
                    
                    if (needToCheckRequiredColumnsPrefixCheck && m_configPtr->m_requiredColumnsPrefixCheck.size() > 0)
                    {   
                        bool allRequiredColumnsPrefixCheckChecked = true;
                        
                        std::map<string, RequiredColumn>::iterator prefixCheck;                        
                        for (prefixCheck = m_configPtr->m_requiredColumnsPrefixCheck.begin(); 
                                prefixCheck != m_configPtr->m_requiredColumnsPrefixCheck.end(); prefixCheck++)
                        {
                            if (false == prefixCheck->second.m_verified)                            
                            {
                                allRequiredColumnsPrefixCheckChecked = false;
                                string columnName = prefixCheck->first;
                                std::vector<std::string>::iterator column = find(m_csvColumns.begin(), m_csvColumns.end(), columnName);
                                if ( column != m_csvColumns.end()) 
                                {// Column exist in this file
                                    string fieldValue;                                                                        
                                    tsvFile->get(0, columnName, fieldValue);
                                    
                                    for(unsigned int i = 0; i < prefixCheck->second.Prefixes.size(); i++)
                                    {
                                        string prefix = prefixCheck->second.Prefixes[i];
                                        size_t found = fieldValue.find(prefix);
                                        if (found != std::string::npos && 0 == (int)found)
                                        {
                                            prefixCheck->second.Prefixes.erase(find(prefixCheck->second.Prefixes.begin(), 
                                                                                    prefixCheck->second.Prefixes.end(),prefix));
                                        }
                                    }                                    
                                    
                                    if (prefixCheck->second.Prefixes.size() == 0)
                                    {
                                        prefixCheck->second.m_verified = true;
                                    }                                    
                                }                                                                 
                            }
                        }
                        if (allRequiredColumnsPrefixCheckChecked)
                        {
                            needToCheckRequiredColumnsPrefixCheck = false;
                        }
                    }
                    
                    std::map<string, RequiredColumn>::iterator  requiredColumn;
                    for (requiredColumn = m_configPtr->m_requiredColumnsNoEmptyValues.begin(); 
                            requiredColumn != m_configPtr->m_requiredColumnsNoEmptyValues.end(); requiredColumn++)
                    {
                        string columnName = requiredColumn->first;
                        std::vector<std::string>::iterator column = find(m_csvColumns.begin(), m_csvColumns.end(), columnName);
                        if ( column != m_csvColumns.end()) 
                        { // Column exist in this file
                            string fieldValue;                                                
                            tsvFile->get(0, columnName, fieldValue);
                            if (fieldValue.empty())
                            {
                                requiredColumn->second.m_error = true;
                                switch(requiredColumn->second.Severity)
                               {
                                    case Error:
                                        requiredColumn->second.ErrorCounter++;
                                        if (requiredColumn->second.ErrorCounter <= RequiredColumn::MaxNumberOfMessages) 
                                        {
                                            std::stringstream message;
                                            message << "Error! Required column '" << requiredColumn->first << "' is missing its value! File '" << Fs::basename(fileName) << "' Data line " << tsvFile->lineNum() << " \n";
                                            message <<  requiredColumn->second.Message;
                                            AptVerboseWrapper::out(0, message.str());
                                        }
                                    break;
                                    case Warning:
                                        requiredColumn->second.WarningCounter++;
                                        if (requiredColumn->second.WarningCounter <= RequiredColumn::MaxNumberOfMessages)
                                        {                            
                                            std::stringstream message;
                                            message << "WARNING! Required column '" << requiredColumn->first << "' is missing its value! File '" << Fs::basename(fileName) << "' Data line " << tsvFile->lineNum() << " \n";
                                            message <<  requiredColumn->second.Message;
                                            AptVerboseWrapper::out(1, message.str());
                                        }
                                    break;
                                    case Information:                                        
                                        if (requiredColumn->second.WarningCounter <= RequiredColumn::MaxNumberOfMessages)
                                        {                            
                                            std::stringstream message;
                                            message << "Column '" << requiredColumn->first << "' is missing its value! File '" << Fs::basename(fileName) << "' Data line " << tsvFile->lineNum() << " \n";
                                            message <<  requiredColumn->second.Message;
                                            AptVerboseWrapper::out(1, message.str());
                                        }
                                    break;
                               }                                                            
                            }
                        }                       
                    }                    
                //
                
                /// Analyze ProbeSet Uniqueness
                string probeSetId;
                tsvFile->get(0, PROBE_SET_ID_IN_CSV, probeSetId);
                if (probeSetId.empty())
                {
                     std::stringstream message;                                                
                     message << "Probe Set ID column in empty at file '" << fileName <<" data line "<< tsvFile->lineNum();
                     throw Except(message.str().c_str());
                }                
                map<string, int>::iterator find = probeSetIds.find(probeSetId);                    
                if (find != probeSetIds.end())
                { // Found duplicated probe set id
                    if (fileIndex == 0)
                    { //  Analyze Probe Set Id uniqueness only for the first file, for other files records will be updated 
                                   
                        //First, find out if dup Probe Set ID lines are exactly the same.
                        std::vector<string> firstRecord;
                        std::vector<string> currentRecord;
                        for(int columnIndex = 0; columnIndex < tsvFile->GetColumnCount(); columnIndex++)
                        { // Get current dup record fields
                            string fieldValue;
                            tsvFile->get(0, columnIndex, fieldValue);
                            currentRecord.push_back(fieldValue);
                        }
                        
                        ////// Get previous dup record fields
                        
                       int saveCurrentLinePos = tsvFile->lineNumber();
                       tsvFile->gotoLine(find->second - 1);
                       for(int columnIndex = 0; columnIndex < tsvFile->GetColumnCount(); columnIndex++)
                       {
                            string fieldValue;
                            tsvFile->get(0, columnIndex, fieldValue);
                            firstRecord.push_back(fieldValue);
                       }
                       tsvFile->gotoLine(saveCurrentLinePos);
                        
                        bool identicalRecords = true;
                        if (firstRecord.size() == currentRecord.size())
                        {
                            for (unsigned int r = 0; r < firstRecord.size(); r++)
                            {
                                if (firstRecord[r] != currentRecord[r])
                                {
                                    affx::dequote(firstRecord[r]);
                                    affx::dequote(currentRecord[r]);
                                    if (firstRecord[r] != currentRecord[r])
                                    {                                
                                        identicalRecords = false;
                                        break;            
                                    }
                                }
                            }                        
                        }
                        else
                        {
                            identicalRecords = false;                        
                        }
                        
                        if (identicalRecords)
                        {// If the whole line in the txt files is duplicated, we do not treat it as an error, just log a warning.    
                            
                            if (m_identicalLinesToSkip.size() <= 50000)
                            {
                                m_identicalLinesToSkip.push_back(tsvFile->lineNumber());
                                                     
                                std::stringstream message;
                                message << "WARNING! Identical record was detected in '"<< fileName << "'! Probe Set ID '" << probeSetId <<"'";
                                message << " First record: Data line " << find->second; 
                                message << " Second record: Data line " << tsvFile->lineNum();                        
                                message << " The second record will be skipped."; 
                                
                                AptVerboseWrapper::out(0, message.str());
                            }
                            else
                            {// For the performance reason
                                 throw Except("There are too many duplicated identical record in the file(s). See the log file for details.");
                            }
                        }
                        else
                        {                            
                            probeSetIdDuplicationCounter++;
                            if (probeSetIdDuplicationCounter < limitOfErrorsInLogFile)
                            {
                                std::stringstream dupMessage;
                                dupMessage << "Error! The following Probe Set ID '" << probeSetId <<"' is duplicated and fields value(s) are different or their order is different!";
                                dupMessage << " File Name: '" << fileName << "'";
                                dupMessage << " First record: Data line " << find->second; 
                                dupMessage << " Second record: Data line " << tsvFile->lineNum() << " !";
                                AptVerboseWrapper::out(0, dupMessage.str());
                            }
                            else
                            {
                               throw Except(DUPLICATED_PROBSET_ID_ERROR);
                            }                            
                        }
                    }
                    else
                    { // 2, 3 .. n files
                        m_linesWithDupProbeIds.insert(std::make_pair(IntsToLong(fileIndex, tsvFile->lineNum()), 0));
                    }
                }
                else
                {                    
                    probeSetIds.insert(std::make_pair(probeSetId, tsvFile->lineNum()));
                }
                
                // Update progress indicator
                std::fstream::pos_type p = tsvFile->line_fpos();
                // fpos_t currectPosition = p.seekpos();
                
                step = (int) (((int) p) / stepSize);                
                if (step > lastStep)
                {                       
                    int numberOfStepsToDisplay = step - lastStep; 
                    if (numberOfStepsToDisplay > numberOfSteps)
                    {
                        numberOfStepsToDisplay = numberOfSteps;
                    }
                    for (int s = 0; s < numberOfStepsToDisplay; s++)
                    { // this "for" is needed for very small files
                        AptVerboseWrapper::progressStep(1);
                    }
                    lastStep = step;
                }
                //////////////////////////////////////////
                 
            }                        
            
            for (int s = step; s < numberOfSteps; s++)
            {
                AptVerboseWrapper::progressStep(1);
            }
            AptVerboseWrapper::progressEnd(1, "\nDone"); 
        }         
        tsvFile->close();
        delete(tsvFile);
    }
        
    if (probeSetIdDuplicationCounter > 0)
    {
        throw Except(DUPLICATED_PROBSET_ID_ERROR);
    }
    
    ReportRequiredColumns();
            
    if (m_identicalLinesToSkip.size() > 0)
    {
        AptVerboseWrapper::out(0, "WARNING! Some lines with identical record were detected! See the log file for details.\n", true);
    }
    
    ReportNumberOfRecordsToInsertUpdate(totalLinesCount);    
}

void AnnotationTable::ReportRequiredColumns()
{

    std::map<string, RequiredColumn>::iterator requiredColumn;
            
    string requiredColumnErrorMessage;
    string requiredColumnWarningMessage;
    string requiredColumnInfoMessage;
    for (requiredColumn = m_configPtr->m_requiredColumnsNoEmptyValues.begin(); requiredColumn != m_configPtr->m_requiredColumnsNoEmptyValues.end(); requiredColumn++)
        {    
            if (requiredColumn->second.m_error)
            {                 
                switch(requiredColumn->second.Severity)
                {
                    case Error:
                        if (requiredColumnErrorMessage != "")
                        {
                            requiredColumnErrorMessage += ", ";
                        }            
                        requiredColumnErrorMessage += "'" + requiredColumn->first + "'";
                        break;
                    case Warning:
                        if (requiredColumnWarningMessage != "")
                        {
                            requiredColumnWarningMessage += ", ";
                        }            
                        requiredColumnWarningMessage += "'" + requiredColumn->first + "'";
                        break;
                    case Information:
                        if (requiredColumnInfoMessage != "")
                        {
                            requiredColumnInfoMessage += ", ";
                        }            
                        requiredColumnInfoMessage += "'" + requiredColumn->first + "'";
                        break;
                }            
            }
        }
    
    
    if (requiredColumnInfoMessage != "")
    {
        requiredColumnInfoMessage = "The following column(s) have empty values: " + requiredColumnInfoMessage;
        requiredColumnInfoMessage += ". See the log file for details.";
        Verbose::out(1, requiredColumnInfoMessage.c_str());
    }
    
    if (requiredColumnWarningMessage != "")
    {
        requiredColumnWarningMessage = "WARNING! The following recommended column(s) have empty values: " + requiredColumnWarningMessage;
        requiredColumnWarningMessage += ". See the log file for details.";
        Verbose::out(0, requiredColumnWarningMessage.c_str());
    }
    
    if (requiredColumnErrorMessage != "")
    {
        requiredColumnErrorMessage = "The following required column(s) have empty values: " + requiredColumnErrorMessage;
        requiredColumnErrorMessage += ". See the log file for details.";
        throw Except(requiredColumnErrorMessage.c_str());
    }
    
    std::stringstream errorMessage;
    std::stringstream warningMessage;
    std::stringstream infoMessage;
    int errorCounter = 0;
    int warningCounter = 0;
    int infoCounter = 0;
    
    std::map<string, RequiredColumn>::iterator requiredColunm;
                        
    for (requiredColunm = m_configPtr->m_requiredColumnsSomeValue.begin(); requiredColunm != m_configPtr->m_requiredColumnsSomeValue.end(); requiredColunm++)
    {
        if (false == requiredColunm->second.m_verified)
        {
            switch(requiredColunm->second.Severity)
            {
                case Error:
                    if (errorCounter == 0)
                    { 
                        errorMessage << "The following required column(s) are entirely empty: ";
                    }
                    else
                    {
                        errorMessage << ", ";
                    }
                    errorMessage << "'" << requiredColunm->first/*Column Name*/ << "'";
                    errorMessage << " - " << requiredColunm->second.Message << " ";
                    
                    errorCounter++;
                    break;
                case Warning:
                    if (warningCounter == 0)
                    { 
                        warningMessage << "WARNING! The following recommended column(s) are entirely empty: ";
                    }
                    else
                    {
                        warningMessage << ", ";
                    }
                    warningMessage << "'" << requiredColunm->first/*Column Name*/ << "'";
                    warningMessage << " - " << requiredColunm->second.Message << " ";
                                    
                    warningCounter++;
                    break;
                case Information:
                    if (infoCounter == 0)
                    { 
                        infoMessage << "The following column(s) are entirely empty: ";
                    }
                    else
                    {
                        infoMessage << ", ";
                    }
                    infoMessage << "'" << requiredColunm->first/*Column Name*/ << "'";
                    infoMessage << " - " << requiredColunm->second.Message << " ";
                                    
                    infoCounter++;
                    break;
            }            
        }
    }
    if (infoCounter > 0)
    {
        infoMessage << ". Should have at least a single value.";                    
        Verbose::out(0, infoMessage.str());
    }
    
    if (warningCounter > 0)
    {
        warningMessage << ". Should have at least a single value.";                    
        Verbose::out(0, warningMessage.str());
    }               
  
    if (errorCounter > 0)
    {
        errorMessage << ". Should have at least a single value.";                    
        throw Except(errorMessage.str().c_str());
    }
    
    errorMessage.clear();
    warningMessage.clear();
    infoMessage.clear();
    errorCounter = 0;
    warningCounter = 0;
    infoCounter = 0;
    
                        
    for (requiredColunm = m_configPtr->m_requiredColumnsPrefixCheck.begin(); requiredColunm != m_configPtr->m_requiredColumnsPrefixCheck.end(); requiredColunm++)
    {
        if (false == requiredColunm->second.m_verified)
        {
            switch(requiredColunm->second.Severity)
            {
                case Error:
                    if (errorCounter == 0)
                    { 
                        errorMessage << "The following column(s) is missing the required prefix(es): ";
                    }
                    else
                    {
                        errorMessage << "; ";
                    }
                    errorMessage << "Column '" << requiredColunm->first/*Column Name*/ << "'";
                    
                    for(unsigned int i = 0; i < requiredColunm->second.Prefixes.size(); i++)
                    {
                        if (i == 0)
                        {
                            errorMessage << " Prefix(es) - ";
                        }                 
                        else
                        {
                            errorMessage << ", ";
                        }                       
                        errorMessage << "'" << requiredColunm->second.Prefixes[i] << "' ";
                    }
                    errorMessage << " " << requiredColunm->second.Message << " ";
                    
                    errorCounter++;
                break;
                case Warning:
                    if (warningCounter == 0)
                    { 
                        warningMessage << "WARNING! The following column(s) is missing the required prefix(es): ";
                    }
                    else
                    {
                        warningMessage << ", ";
                    }
                    
                    warningMessage << "Column '" << requiredColunm->first/*Column Name*/ << "'";
                    
                    for(unsigned int i = 0; i < requiredColunm->second.Prefixes.size(); i++)
                    {
                        if (i == 0)
                        {
                            warningMessage << " Prefix(es) - ";
                        }                 
                        else
                        {
                            warningMessage << ", ";
                        }                       
                        warningMessage << "'" << requiredColunm->second.Prefixes[i] << "' ";
                    }
                    warningMessage << " " << requiredColunm->second.Message << " ";                                
                                    
                    warningCounter++;
                break;
                case Information:
                    if (infoCounter == 0)
                    { 
                        infoMessage << "The following column(s) is missing the prefix(es): ";
                    }
                    else
                    {
                        infoMessage << ", ";
                    }
                    
                    infoMessage << "Column '" << requiredColunm->first/*Column Name*/ << "'";
                    
                    for(unsigned int i = 0; i < requiredColunm->second.Prefixes.size(); i++)
                    {
                        if (i == 0)
                        {
                            infoMessage << " Prefix(es) - ";
                        }                 
                        else
                        {
                            infoMessage << ", ";
                        }                       
                        infoMessage << "'" << requiredColunm->second.Prefixes[i] << "' ";
                    }
                    infoMessage << " " << requiredColunm->second.Message << " ";                                
                                    
                    infoCounter++;
                break;
            }            
        }
    }
    
    if (infoCounter > 0)
    {        
        Verbose::out(1, infoMessage.str());
    } 
    
    if (warningCounter > 0)
    {        
        Verbose::out(0, warningMessage.str());
    }               
  
    if (errorCounter > 0)
    {     
        throw Except(errorMessage.str().c_str());
    }
    
}

void AnnotationTable::ReportNumberOfRecordsToInsertUpdate(int TotalLinesCount)
{
    std::stringstream message;
    message << "\n" << (TotalLinesCount - m_linesWithDupProbeIds.size()) << " records will be imported into the database. ";
    if (m_linesWithDupProbeIds.size() > 0)
    {
        message << m_linesWithDupProbeIds.size() << " update(s)";
    }
    message << "\n";
    
    AptVerboseWrapper::out(2, message.str(), true);
}

void AnnotationTable::BindFieldValue(AnnotationTsv* TsvFilePtr, std::string &DbColumnName, std::string &CsvColumnName, std::string &FieldValue)
{   
    std::string columnType = m_columnNameToColumnTypeMap[DbColumnName]; 
        
    int parameterIndex;
    
    if (m_updateRowMode)
    {    
        parameterIndex = m_updateBindingMap[DbColumnName];      
    }
    else
    {
        parameterIndex = m_insertBindingMap[DbColumnName];     
    }
    
    /// Special case for the chromosome column
    if (columnType == DB_TYPE_INTEGER && DbColumnName == "Chr_id")
    {
        int64_t chr_id_Key = m_chromosomeTablePtr->GetKey(FieldValue);
        if (sqlite3_bind_int64(m_currentStatement, parameterIndex, chr_id_Key) != SQLITE_OK) 
        {
            ReportTsvErrorStopDbInsertion("Could not bind integer ", TsvFilePtr->lineNum(), FieldValue, CsvColumnName, columnType);
        }
        return;
    }
    ////////////
    
    if (FieldValue == "---" || FieldValue == "DBNULL")
    {
        return;
    }
    
    if (columnType == DB_TYPE_INTEGER)
    {   
        if (false == FieldValue.empty())
        {
            int64_t int64Val;
            std::stringstream iss(FieldValue);            
            if ((iss >> int64Val).fail())
            {
                ReportTsvErrorStopDbInsertion("Could not convert to integer ",TsvFilePtr->lineNum(), FieldValue, CsvColumnName, columnType);
            }
            else
            {
                stringstream checkConvertion;
                checkConvertion << int64Val;
                
                if (checkConvertion.str() != FieldValue)
                {
                    ReportTsvErrorStopDbInsertion("Could not convert to integer ",TsvFilePtr->lineNum(), FieldValue, CsvColumnName, columnType);
                }
                else
                {
                    if (sqlite3_bind_int64(m_currentStatement, parameterIndex, int64Val) != SQLITE_OK) 
                    {
                        ReportTsvErrorStopDbInsertion("Could not bind integer ", TsvFilePtr->lineNum(), FieldValue, CsvColumnName, columnType);
                    }
                }
            }
        }
    }
    else if ( columnType == DB_TYPE_REAL)
    { 
        if (false == FieldValue.empty())
        {
            double d;                        
            std::istringstream iss(FieldValue);            

            if ((iss >> d).fail())
            {
                ReportTsvErrorStopDbInsertion("Could not convert to float", TsvFilePtr->lineNum(), FieldValue, CsvColumnName, columnType);
            }
            else
            {
                if (sqlite3_bind_double(m_currentStatement, parameterIndex, d) != SQLITE_OK) 
                {
                    ReportTsvErrorStopDbInsertion("Could not bind float", TsvFilePtr->lineNum(), FieldValue, CsvColumnName, columnType);
                }
             }
        }
    }
    else if ( columnType == DB_TYPE_TEXT)
    {    
        int result = sqlite3_bind_text(m_currentStatement, parameterIndex, FieldValue.c_str(), FieldValue.length(), SQLITE_TRANSIENT);
        if (result != SQLITE_OK) 
        {
            ReportTsvErrorStopDbInsertion("Could not bind string", TsvFilePtr->lineNum(), FieldValue, CsvColumnName, columnType);
        }
    }
    else
    {
        ReportTsvErrorStopDbInsertion("Internal error! Cannot get a column type for the column:", TsvFilePtr->lineNum(), FieldValue, CsvColumnName, columnType);                    
    }                
}

void AnnotationTable::BuildCsvColumnList(AnnotationTsv* TsvFilePtr)
{   
    m_csvColumns.clear();
    for(int columnIndex = 0; columnIndex < TsvFilePtr->GetColumnCount(); columnIndex++)
    {
        string columnName = TsvFilePtr->getColumnName(0, columnIndex);
        StringUtilities::CheckColumnName(columnName, columnIndex);
        m_csvColumns.push_back(columnName);
    }
}

void AnnotationTable::BuildCsvColumnIndexToColumnNameMap(AnnotationTsv* TsvFilePtr)
{   
    m_csvColumnIndexToColumnNameMap.clear();
    for(int columnIndex = 0; columnIndex < TsvFilePtr->GetColumnCount(); columnIndex++)
    {
        string columnName = TsvFilePtr->getColumnName(0, columnIndex);
        StringUtilities::CheckColumnName(columnName, columnIndex);
        if (columnName == PROBE_SET_ID_IN_CSV)
        {
            m_probeSetIdIndex = columnIndex;
        }
        m_csvColumnIndexToColumnNameMap.insert(std::make_pair(columnIndex, columnName));
    }
}

void AnnotationTable::GetCsvColumnName(int ColumnIndex, string &CsvColumnName)
{
    map<int,string>::iterator csvNameIterator = m_csvColumnIndexToColumnNameMap.find(ColumnIndex);    
    if ( csvNameIterator != m_csvColumnIndexToColumnNameMap.end())
    { 
        CsvColumnName = csvNameIterator->second; // As it is in CSV
    }
    else
    {
        FinalizeStatements();
        throw Except("Internal error! Cannot map column name!");                                                                    
    }
}
string AnnotationTable::GetRowCount(string ProbeSet_Id)
{
    char *errorMessage = NULL;
    string statement = "Select count(*) from Annotations WHERE ProbeSet_ID = '" + ProbeSet_Id + "'";
    if(sqlite3_exec(m_annotationDb, statement.c_str(), AnnotationTable::RowCountCallback, 0, &errorMessage) != SQLITE_OK)
      {
        std::string message = "WARNING! Cannot get row count from Annotations table! ";
        if (errorMessage != NULL)
        {
            message += errorMessage;
            delete(errorMessage);
        }
        AptVerboseWrapper::out(0, message, true);
      }
      return m_rowCount;
}
int AnnotationTable::RowCountCallback(void *NotUsed, int argc, char **argv, char **azColName)
{    
    m_rowCount = argv[0];
    return 0;
}

const int parConversionArraySize = 8;
const char* parConversionArray[parConversionArraySize][3] = 
{
    {"1","0","1"},
    {"0","1","2"},
    {"0","0","0"},
    {"---","0","0"},
    {"0","---","0"},
    {"---","1","2"},
    {"1","---","1"},
    {"---","---","DBNULL"}
};

bool AnnotationTable::InsertTsvFiles(std::vector<std::string> ListOfTsvFileNames)
{   
    if (ListOfTsvFileNames.size() == 0) 
    {
        return true;
    }
    else if (ListOfTsvFileNames.size() == 1)
    {
        AptVerboseWrapper::out(1, IMPORTING_FILE_INTO_DATABASE, true);
    }
    else
    {
        AptVerboseWrapper::out(1, IMPORTING_FILES_INTO_DATABASE, true);
    }
    
    m_stopDbInsertion = false;    
    
    m_updateStatement = NULL;
    m_insertStatement = NULL;
    
    bool anyDupsInFiles = m_linesWithDupProbeIds.size() > 0;
    
    for (unsigned int fileIndex = 0; fileIndex < ListOfTsvFileNames.size(); fileIndex++)
    {
        int updateMessageCounter = 0;
        std::string fileName = ListOfTsvFileNames[fileIndex];
        
        AnnotationTsv* tsvFile = new AnnotationTsv();                 
        if(tsvFile->open(fileName) != TSV_OK)
        {            
            delete(tsvFile); 
            std::string message = "Unable to open file '" + fileName + "'";
            FinalizeStatements();
            throw Except(message.c_str());		
        }
        else        
        { 
            int64_t fileSize = Fs::fileSize(fileName, false);
            std::string scale = "1% ----------- 50% ----------- 100%";
            int numberOfSteps = scale.size();
            int stepSize = (int)(fileSize / numberOfSteps);
            int step = 0;
            int lastStep = 0; 
                      
            AptVerboseWrapper::progressBegin(1, "Importing file " + Fs::basename(fileName) + "\n" , numberOfSteps, 1, numberOfSteps);
            AptVerboseWrapper::out(1, scale + "\n", true);
            
            // Read tsvFile file                                    
            BuildCsvColumnIndexToColumnNameMap(tsvFile);
            PrepareInsertStatement(tsvFile);
            PrepareUpdateStatement(tsvFile); 
   
            m_updateRowMode = false;
            m_currentStatement = m_insertStatement;
    
            int columnCount = tsvFile->GetColumnCount();                                    
            string csvColumnName;
            std::string fieldValue;                        
            
            while (tsvFile->nextLine() == TSV_OK)
            {         
                
                if (fileIndex == 0)    
                {
                    std::vector<int>::iterator isLineToSkip = find(m_identicalLinesToSkip.begin(), m_identicalLinesToSkip.end(), tsvFile->lineNumber());                
                    if ( isLineToSkip != m_identicalLinesToSkip.end())
                    {
                        continue;
                    }
                }                 
                
                string par1CsvColumnValue = "";
                string par2CsvColumnValue = "";
                int parColumnsCounter = 0;
                
                if (anyDupsInFiles)
                {
                    if (fileIndex == 0)
                    {// Probe Set ID in the FIRST file must be unique, so first file is always in the insert mode
                        m_updateRowMode = false;
                    }
                    else
                    {
                        if (m_linesWithDupProbeIds.find(IntsToLong(fileIndex, tsvFile->lineNum())) != m_linesWithDupProbeIds.end())
                        {
                            m_updateRowMode = true;
                            m_currentStatement = m_updateStatement;
                        }
                        else
                        {
                            m_updateRowMode = false;
                            m_currentStatement = m_insertStatement;
                        }                                                
                    }
                }
                for(int columnIndex = 0; columnIndex < columnCount; columnIndex++)
                {
                    
                    tsvFile->get(0, columnIndex, fieldValue);
                    
                    if (updateMessageCounter < 50)
                    {
                        if (columnIndex == m_probeSetIdIndex && m_updateRowMode)
                        {
                            std::stringstream message;
                            message <<  "Record with Probe Set ID '"  << fieldValue << "' will be updated from file '" << fileName << "' Data line " << tsvFile->lineNum();
                            AptVerboseWrapper::out(3, message.str());
                            updateMessageCounter++;
                        }
                    }
                    else if (updateMessageCounter == 50)
                    {                        
                        AptVerboseWrapper::out(3, "Too many update messages in the log file! Other update messages will be omitted.");
                        updateMessageCounter++;
                    }
                    
                    
                    GetCsvColumnName(columnIndex, csvColumnName);
                    
                    if (m_configPtr->IsCloneableColumn(csvColumnName))
                    {
                        std::vector<string> columnClones = m_configPtr->GetColumnClones(csvColumnName);
                        m_configPtr->ConvertTsvToDbValue(columnClones[0], fieldValue);
                        for (unsigned int i = 0; i < columnClones.size(); i++)
                        {
                           BindFieldValue(tsvFile, columnClones[i], csvColumnName, fieldValue);
                        }
                    }
                    else if (csvColumnName == m_configPtr->m_par1CsvColumnName || csvColumnName == m_configPtr->m_par2CsvColumnName)
                    {// CSVTOSQLITE-74 PAR Columns handling. Hardcoded for now. If we will have similar cases with other columns, 
                     // it should be implemented as a customizable feature like the columns cloning
                        
                        parColumnsCounter++;
                        if (csvColumnName == m_configPtr->m_par1CsvColumnName)
                        {
                            par1CsvColumnValue = fieldValue;
                        }
                        else
                        {
                            par2CsvColumnValue = fieldValue;
                        }
                        if (parColumnsCounter == 2)
                        {
                            string parDbColumnValue = "";
                            for (int i = 0; i < parConversionArraySize; i++)
                            {
                                if (par1CsvColumnValue == parConversionArray[i][0] && par2CsvColumnValue == parConversionArray[i][1])
                                {
                                    parDbColumnValue = parConversionArray[i][2];
                                    break;
                                }
                            }                                                             
                            
                            if (parDbColumnValue != "DBNULL")
                            {
                                if (parDbColumnValue != "")                           
                                {
                                    BindFieldValue(tsvFile, m_configPtr->m_parDbColumnName, csvColumnName, parDbColumnValue);
                                }
                                else
                                {
                                   std::stringstream message;
                                   message << "'" << m_configPtr->m_par1CsvColumnName << "' - '" << m_configPtr->m_par2CsvColumnName << "'" << " columns pair has incorrect values. ";
                                   message << "'" << m_configPtr->m_par1CsvColumnName << "' = '" << par1CsvColumnValue << "', ";
                                   message << "'" << m_configPtr->m_par2CsvColumnName << "' = '" << par2CsvColumnValue << "'. ";
                                   ReportTsvErrorStopDbInsertion(message.str() , tsvFile->lineNum(), "", "", ""); 
                                }
                            }                            
                        }
                    }
                    else
                    {// "Normal" column                
                        string dbColumnName = GetDbColumnNameFromLocalizationTable(csvColumnName);                                                
                                                
                        m_configPtr->ConvertTsvToDbValue(csvColumnName, fieldValue);
                        
                        BindFieldValue(tsvFile, dbColumnName, csvColumnName, fieldValue);                        
                    }
                }
                
                if (false == m_stopDbInsertion)
                {
                    int result = sqlite3_step(m_currentStatement);
                    if (result != SQLITE_DONE) 
                    {  
                        std::stringstream message;
                        if (result == SQLITE_CONSTRAINT)
                        {
                            message << "Constraint violation at the data line  " << tsvFile->lineNum();
                        }
                        else
                        {
                            message << "Could not execute SQL statement. Data line " << tsvFile->lineNum();
                        }
                        
                        FinalizeStatements();
                        
                        throw Except(message.str().c_str());
                    }    
                }
                sqlite3_reset(m_currentStatement);
                sqlite3_clear_bindings(m_currentStatement);
                
                std::fstream::pos_type p = tsvFile->line_fpos();
                //fpos_t currectPosition = p.seekpos();
               
                step = (int) (((int) p) / stepSize);
                
                if (step > lastStep)
                {                       
                    int numberOfStepsToDisplay = step - lastStep; 
                    if (numberOfStepsToDisplay > numberOfSteps)
                    {
                        numberOfStepsToDisplay = numberOfSteps;
                    }
                    for (int s = 0; s < numberOfStepsToDisplay; s++)
                    { // this "for" is needed for very small files
                        AptVerboseWrapper::progressStep(1);
                    }
                    lastStep = step;
                }
             }
             
             tsvFile->close();             
             delete(tsvFile);
             for (int s = step; s < numberOfSteps; s++)
             {
                AptVerboseWrapper::progressStep(1);
             }
             AptVerboseWrapper::progressEnd(1, "\nDone");
         }
    }
    
    FinalizeStatements();
    
    return !m_stopDbInsertion;
}

void AnnotationTable::FinalizeStatements()
{
    if (m_insertStatement != NULL)
    {
        sqlite3_finalize(m_insertStatement);
        m_insertStatement = NULL;
    }
    
    if (m_updateStatement != NULL)
    {
        sqlite3_finalize(m_updateStatement);
        m_updateStatement = NULL;
    }
}

string AnnotationTable::GetDbColumnNameFromLocalizationTable(string CsvColumnName)
{
    string dbColumnName = m_localizationTablePtr->GetDbColumnName(CsvColumnName);
                        
    if (dbColumnName.empty())
    {
        FinalizeStatements();
        std::string errorStr;
        errorStr = "Internal error! Cannot get db column name for '" + CsvColumnName + "'!";
        throw Except(errorStr.c_str());
    } 
     
    return dbColumnName;
}

void AnnotationTable::ReportTsvErrorStopDbInsertion(std::string Message, int LineNumber, std::string FieldValue, 
                                                         std::string ColumnName, std::string ColumnType)
{   
    m_stopDbInsertion = true;
         
    std::stringstream message;
    message << Message;
    
    message << " Data line " << LineNumber;        
    if (ColumnName != "")
    {
        message << ", column '" << ColumnName << "'";
    }
    if (ColumnType != "")
    {
        message << ", column type '" << ColumnType << "'";
    }
    if (FieldValue != "")
    {
        message << ", value '" << FieldValue << "'";
    }
    
    throw Except(message.str().c_str()); 
}
    
 
