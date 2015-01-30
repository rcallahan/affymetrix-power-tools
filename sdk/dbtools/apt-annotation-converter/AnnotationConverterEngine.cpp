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

//

#include "calvin_files/exception/src/ExceptionBase.h"
#include "calvin_files/utils/src/StringUtils.h"
#include "util/Fs.h"
#include "util/Verbose.h"

#include "AnnotationTsv.h"
using namespace affx;


#include "AnnotationConverterEngine.h"
#include "StringUtilities.h"
#include "ResourceStrings.h"
#include "AptVerboseWrapper.h"

//
#include <cassert>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>



AnnotationConverterEngine::Reg AnnotationConverterEngine::reg;

AnnotationConverterEngine * AnnotationConverterEngine::FromBase(BaseEngine *engine)
{
    if (engine != NULL && engine->getEngineName() == AnnotationConverterEngine::EngineName())
    {
        return (AnnotationConverterEngine *)engine;
    }
    return NULL;
}

// constructor
AnnotationConverterEngine::AnnotationConverterEngine() 
{
    defineOptions();
    m_annotationDb = NULL;
    m_csvHeadersToInformationTablePtr = NULL;
    m_annotationTable = NULL;
    m_localizationTable = NULL;
    m_chromosomeTable = NULL;        
}

AnnotationConverterEngine::~AnnotationConverterEngine() 
{
    clear();
}

void AnnotationConverterEngine::clear()
{    
    if (m_csvHeadersToInformationTablePtr!= NULL) delete m_csvHeadersToInformationTablePtr;
    if (m_annotationTable!= NULL) delete m_annotationTable;
    if (m_localizationTable!= NULL) delete m_localizationTable;
    if (m_chromosomeTable!= NULL) delete m_chromosomeTable;
    if (m_cdfInformationTable!= NULL) delete m_cdfInformationTable;        
}

static const char* SWITCH_DB_FILE = "db-file";
static const char* SWITCH_DB_TEMPLATE = "db-template";
static const char* SWITCH_TSV_FILE = "tsv-file";
static const char* SWITCH_TSV_FILES = "tsv-files";
static const char* SWITCH_ARRAY_CONFIG = "array-config";
static const char* SWITCH_ARRAY_SET = "array-set";
static const char* MESSAGE_DELAY = "message-delay";

void AnnotationConverterEngine::defineOptions() 
{
    defineOptionSection("Input Options");
    defineOption("", SWITCH_DB_TEMPLATE, PgOpt::STRING_OPT, "Sqllite Annotation database template file","");
    defineOption("", SWITCH_DB_FILE, PgOpt::STRING_OPT, "Target Sqllite Annotation database file","");    
    defineOption("", SWITCH_TSV_FILE, PgOpt::STRING_OPT, "Delimited text file to import in to Annotation table. Use either TCV or CSV-file format","");
    defineOption("", SWITCH_TSV_FILES, PgOpt::STRING_OPT, "Text file specifying TSV files to process, one per line","");
    defineOption("", SWITCH_ARRAY_CONFIG, PgOpt::STRING_OPT, "Array type specific configuration file.","");
    defineOption("", SWITCH_ARRAY_SET, PgOpt::STRING_OPT, "Array set name.","");
    defineOption("", MESSAGE_DELAY, PgOpt::STRING_OPT, "Message Delay (in milliseconds) for the UI tool.","");
}

void AnnotationConverterEngine::defineStates() { }

/**
* Make sure that our options are sane. Call Err::errAbort if not.
*/

// -db-file GenomeWideSNP_6.na30.annot.db  -tsv-file GenomeWideSNP_6.na30.annot.small.csv  -array-config GenomeWideSNP_6.arrayconfig
// or
// -db-file GenomeWideSNP_6.na30.annot.db  -tsv-files ListOfTsvToLoad.txt  -array-config GenomeWideSNP_6.arrayconfig

void AnnotationConverterEngine::checkOptionsImp()
{
           
     //Add template name and copy into db every time     
     CheckSingleRequiredSwitch(SWITCH_DB_TEMPLATE, m_dbTemplate);
     CheckSingleRequiredSwitch(SWITCH_DB_FILE, m_dbFileName);
     CheckSingleRequiredSwitch(SWITCH_ARRAY_CONFIG, m_arraySpecificConfig);
     
     m_tsvFileName = getOpt(SWITCH_TSV_FILE);
     m_tsvListFileName = getOpt(SWITCH_TSV_FILES);
	 m_arraySetName = getOpt(SWITCH_ARRAY_SET);
     
     if (!Fs::fileExists(m_dbTemplate))
     {
          Err::errAbort("Template file '" + m_dbTemplate + "'" + " doesn't exist");
     }
    
     if (m_tsvFileName.empty() && m_tsvListFileName.empty()) 
     {
        std::stringstream message;
        message << "Either '" << SWITCH_TSV_FILE << "' or '" << SWITCH_TSV_FILES << "' must be specified!";
        Err::errAbort(message.str());
     }
     
     if (m_tsvFileName != "" && m_tsvListFileName != "") 
     {
        std::stringstream message;
        message << "Both '" << SWITCH_TSV_FILE << "' and '" << SWITCH_TSV_FILES << "' are specified. Only one of them can be used at a time!";        
        Err::errAbort(message.str());
     } 
     
     m_listOfTsvFiles.clear();     
     if (m_tsvFileName != "")
     { // "tsv-file" specified - loading a single file 
        m_listOfTsvFiles.push_back(m_tsvFileName);
     }
     else
     {
        GetTsvFiles(); 
     } 
     
     string messageDelay =  getOpt(MESSAGE_DELAY);
     if (messageDelay != "")  
     {
        AptVerboseWrapper::Delay = atoi(messageDelay.c_str());
     }
}


void AnnotationConverterEngine::GetTsvFiles() 
{
    if (!Fs::fileExists(m_tsvListFileName))
    {
         Err::errAbort("List file '" + m_tsvListFileName + "'" + " doesn't exist");
    }
    
    m_listOfTsvFiles.clear();
    
    std::ifstream infile(m_tsvListFileName.c_str(), std::ios_base::in);
    
    std::string fileName;
    while (getline(infile, fileName, '\n'))
    {
        Util::trimString(fileName);
        if (false == fileName.empty())
        {                
            std::vector<string>::iterator pos = find(m_listOfTsvFiles.begin(), m_listOfTsvFiles.end(), fileName);
            if ( pos != m_listOfTsvFiles.end())
            {
                Err::errAbort("Duplicated file name '" + fileName + "' in the list file!" );
            }     
            m_listOfTsvFiles.push_back (fileName);
        }
    }
    
    if (m_listOfTsvFiles.size() == 0)
    {
         Err::errAbort("List file '" + m_tsvListFileName + "'" + " has no files to import!");
    }
    
    std::vector<std::string>::iterator fileIterator;
    for (fileIterator = m_listOfTsvFiles.begin(); fileIterator != m_listOfTsvFiles.end(); fileIterator++)
    {
        std::string fileName = *fileIterator;
        if (!Fs::fileExists(fileName))
        {
            Err::errAbort("File '" + fileName + "' listed in the '" + m_tsvListFileName + "' cannot be found!");
        }
    }
    
    std::stringstream message;
    message << "'" << m_tsvListFileName << "'" << " has " << m_listOfTsvFiles.size() << " files to load";
    
    AptVerboseWrapper::out(2, message.str());    
}

void AnnotationConverterEngine::CheckSingleRequiredSwitch(const char* SwitchName, std::string& SwitchValue)
{
     SwitchValue = getOpt(SwitchName);
         
     if (SwitchValue.empty()) 
     {
        std::stringstream message;
        message  << "Must specify '" << SwitchName << "'!";
        Err::errAbort(message.str());        
     }
}

void AnnotationConverterEngine::LoadConfiguration(string FileName)
{        
    try
    {        
        if (!Fs::fileExists(FileName))
        {
            std::stringstream message;
	    message << CANNOT_FIND_CONFIGURATION_FILE << FileName << "'!";
            Err::errAbort(message.str()); 
        }
        m_configuration.LoadConfigurationFile(FileName);
    }
    catch (std::exception& exception)
    {
        std::string message = 	ERROR_LOADING_CONFIGURATION_FILE; 	
	    message += exception.what();
	    Err::errAbort(message);
    }
    catch (...)
    {        
	    Err::errAbort(ERROR_LOADING_CONFIGURATION_FILE);
    }
}



vector<string> AnnotationConverterEngine::m_valueMismatchReportContent;
vector<string> AnnotationConverterEngine::m_missingValueReportContent;


void AnnotationConverterEngine::ExecuteNonQuery(std::string SqlStatement)
{
      char *errorMessage = NULL;
      if(sqlite3_exec(m_annotationDb, SqlStatement.c_str(), NULL, 0, &errorMessage) != SQLITE_OK)
      {
            std::stringstream message;
            message << " Cannot execute query! " << SqlStatement;
            if (errorMessage != NULL)
            {
                message << errorMessage;
                delete(errorMessage);
            }
            throw Except(message.str().c_str());
       }
}

void AnnotationConverterEngine::BeginTransaction()
{
    ExecuteNonQuery("BEGIN TRANSACTION");
}

void AnnotationConverterEngine::CommitTransaction()
{
    ExecuteNonQuery("COMMIT TRANSACTION");
}

void AnnotationConverterEngine::RollbackTransaction()
{
    ExecuteNonQuery("ROLLBACK TRANSACTION");
}

void AnnotationConverterEngine::CloseDbAndRemoveFile()
{        
        sqlite3_close(m_annotationDb);        
        try
        {
            Fs::rm(m_dbFileName, false);
        }
        catch(...)
        {// Deleting the DB file in case of an error is not critical
            AptVerboseWrapper::out(1, "Error deleting the DB file!");
        }
}

sqlite3* AnnotationConverterEngine::OpenDb(std::string dbFileName)
{      
    std::string dbNameUncPath = Fs::convertToUncPath(dbFileName);

    if (!Fs::fileExists(dbNameUncPath))
    {
        Err::errAbort("Database file '" + dbNameUncPath + "'" + " doesn't exist");	
    }
    
    if (sqlite3_open(dbNameUncPath.c_str(), &m_annotationDb) != SQLITE_OK) 
    {
        if (m_annotationDb != NULL)
        {
            sqlite3_close(m_annotationDb);
        }        
        Err::errAbort("Unable to open the database '" + dbNameUncPath + "'");	
     }		
     return m_annotationDb;
}

 void AnnotationConverterEngine::AnalyzeHeaders()
 {        
    if (m_listOfTsvFiles.size() == 0)
    {
        return;
    }  
    else if (m_listOfTsvFiles.size() == 1)
    {
        AptVerboseWrapper::progressBegin(1, ANALYZING_FILE_HEADER, m_listOfTsvFiles.size(), 1, m_listOfTsvFiles.size());
    }
    else
    {
        AptVerboseWrapper::progressBegin(1, ANALYZING_FILES_HEADERS, m_listOfTsvFiles.size(), 1, m_listOfTsvFiles.size());
    }

    std::vector<string> errorColumns; // error if missing
    std::vector<string> warningColumns; // warning if missing
    std::vector<string> informationColumns;// information if missing
    
    RequiredColumn::AddRequiredRecommended(m_configuration.m_requiredColumnsPresence, errorColumns, warningColumns,informationColumns);
    RequiredColumn::AddRequiredRecommended(m_configuration.m_requiredColumnsSomeValue, errorColumns, warningColumns,informationColumns);
    RequiredColumn::AddRequiredRecommended(m_configuration.m_requiredColumnsNoEmptyValues, errorColumns, warningColumns, informationColumns);
    RequiredColumn::AddRequiredRecommended(m_configuration.m_requiredColumnsPrefixCheck, errorColumns, warningColumns, informationColumns);

                
    std::vector<std::string>::iterator fileIterator;
    for (fileIterator = m_listOfTsvFiles.begin(); fileIterator != m_listOfTsvFiles.end(); fileIterator++)
    {
        std::string fileName = *fileIterator;        
        AnnotationTsv* tsvFile = new AnnotationTsv();        
        
        std::vector<string> columnsRequiredInEachFile;
        RequiredColumn::AddRequiredInEachFile(m_configuration.m_requiredColumnsNoEmptyValues, columnsRequiredInEachFile);
        
        if(tsvFile->open(fileName) != TSV_OK)
        {            
            delete(tsvFile);   
            
            std::string message = "Unable to open file '" + fileName + "'";
            throw Except(message.c_str());		
        }
        else
        {                         
            if (tsvFile->headersCount() == 0)
            {
                AptVerboseWrapper::out(1, "WARNING! File '" + fileName + "' has no header." );
            }
            else
            {                                        
                TsvFileHeaderLine* hdr;	                
	            while ((hdr = tsvFile->nextHeaderPtr()) != NULL)
	            {
	                //trim the value with empty ',' from excel export, such as "array_type,,,,,"
	                Util::trimString(hdr->m_value," ,\r\n\t");
	                m_csvHeadersToInformationTablePtr->AddVariable(fileName, hdr->m_key, hdr->m_value);		            
	            }
            }
            m_csvHeadersToInformationTablePtr->AddVariable(fileName, "array_type", m_arraySetName);

            int columnCount = tsvFile->GetColumnCount();
            
            for(int columnIndex = 0; columnIndex < columnCount; columnIndex++)
            {
                
                std::string csvColumnName = tsvFile->getColumnName(0, columnIndex);
                StringUtilities::CheckColumnName(csvColumnName, columnIndex);
                
                std::vector<string>::iterator columnToFind =  find(errorColumns.begin(),errorColumns.end(), csvColumnName);
                if (columnToFind != errorColumns.end())
                {// Ok found it. Can delete 
                    errorColumns.erase(columnToFind);
                }
                                 
                columnToFind =  find(columnsRequiredInEachFile.begin(), columnsRequiredInEachFile.end(), csvColumnName);
                if (columnToFind != columnsRequiredInEachFile.end())
                {// Ok found it. Can delete 
                    columnsRequiredInEachFile.erase(columnToFind);
                }
                                
                columnToFind =  find(warningColumns.begin(),warningColumns.end(), csvColumnName);
                if (columnToFind != warningColumns.end())
                {
                    warningColumns.erase(columnToFind);
                }
                
                columnToFind =  find(informationColumns.begin(),informationColumns.end(), csvColumnName);
                if (columnToFind != informationColumns.end())
                {
                    informationColumns.erase(columnToFind);
                }
                
                if (m_configuration.IsCloneableColumn(csvColumnName))
                {
                    std::vector<string> columnClones = m_configuration.GetColumnClones(csvColumnName);            
                    for (unsigned int i = 0; i < columnClones.size(); i++)
                    {
                        string columnName = columnClones[i];
                        string cloneType = DB_TYPE_TEXT;
                        
                        map<string, string>::iterator columnTypes = m_configuration.m_cloneColumnTypes.find(columnName);        
                        if (columnTypes != m_configuration.m_cloneColumnTypes.end())
                        {
                            cloneType = columnTypes->second;
                        }
                        m_annotationTable->CheckIfNeedToAddColumn(columnName, cloneType);
                    }
                }
                else if (csvColumnName == m_configuration.m_par1CsvColumnName || csvColumnName == m_configuration.m_par2CsvColumnName)
                {
                    m_annotationTable->CheckIfNeedToAddColumn( m_configuration.m_parDbColumnName,m_configuration.m_parDbColumnType);
                }
                else
                {// "Normal" column                                            
                    std::string dbColumnName = m_localizationTable->GetDbColumnName(csvColumnName);
                    if (dbColumnName.empty())
                    { 
                        m_localizationTable->AppendColumnToAddList(csvColumnName);
                        dbColumnName = StringUtilities::MakeColumnName(csvColumnName);
                    }                     
                    m_annotationTable->CheckIfNeedToAddColumn(csvColumnName);
                }                
             }             
         }
         
         tsvFile->close();
         delete(tsvFile); 
         
         if (columnsRequiredInEachFile.size() > 0)
         {
            std::stringstream  message;
            message << "The following required column(s) are missing in the '" << Fs::basename(fileName) << "': ";
            
            for(unsigned int i = 0; i < columnsRequiredInEachFile.size(); i++)
            {
                if (i > 0) message << ", ";    
                
                message << "'" << columnsRequiredInEachFile[i] << "'";
            }
            throw Except(message.str().c_str());
         }
         
         AptVerboseWrapper::progressStep(1);
     }
     
     AptVerboseWrapper::progressEnd(1, "\nDone");
     
     if (informationColumns.size() > 0)
     {
        std::stringstream  message;
        message << "The following recommended column(s) are missing: ";
        
        for(unsigned int i = 0; i < informationColumns.size(); i++)
        {
            if (i > 0) message << ", ";    
            
            message << "'" << informationColumns[i] << "'";
        }
        
        Verbose::out(0, message.str());
     }
     
     if (warningColumns.size() > 0)
     {
        std::stringstream  message;
        message << "Warning! The following recommended column(s) are missing: ";
        
        for(unsigned int i = 0; i < warningColumns.size(); i++)
        {
            if (i > 0) message << ", ";    
            
            message << "'" << warningColumns[i] << "'";
        }
        
        Verbose::out(0, message.str());
     }
     
     if (errorColumns.size() > 0)
     {
        std::stringstream  message;
        message << "The following required column(s) are missing: ";
        
        for(unsigned int i = 0; i < errorColumns.size(); i++)
        {
            if (i > 0) message << ", ";    
            
            message << "'" << errorColumns[i] << "'";
        }
        throw Except(message.str().c_str());
     }          
             
     
     m_csvHeadersToInformationTablePtr->ValueMismatchReport = &AnnotationConverterEngine::ValueMismatchReport;
     m_csvHeadersToInformationTablePtr->MissingValueReport = &AnnotationConverterEngine::MissingValueReport;
     
     ValidationResult validationResult = m_csvHeadersToInformationTablePtr->AnalyzeHeaders();
     if (validationResult == ValueMismatch)
     {
         std::string message = "Header variable(s) have mismatching values:\n";
         for(unsigned int i = 0; i < m_valueMismatchReportContent.size(); i++)
         {
            message += m_valueMismatchReportContent[i] + "\n";
         }                  
         throw Except(message.c_str());
     }
     if (validationResult == MissingValue)
     {
         AptVerboseWrapper::out(0,"WARNING! Header variable(s) have missing or empty values:");
         for(unsigned int i = 0; i < m_missingValueReportContent.size(); i++)
         {
            AptVerboseWrapper::out(0,m_missingValueReportContent[i]);
         }         
     }     
     
     m_annotationTable->AddColumns();
     m_localizationTable->AddRecords();
 
 }
 
              
void AnnotationConverterEngine::AnalyzeTsvFiles()
{   
    AnalyzeHeaders();
    m_annotationTable->AnalyzeFiles(m_listOfTsvFiles);
}

void AnnotationConverterEngine::ValueMismatchReport(string VarName, map<string, string> FileValueMap)
{
    string message = "Variable Name: '" + VarName + "' ";
    
    map<string, string>::iterator pos;    
    for (pos = FileValueMap.begin(); pos != FileValueMap.end(); pos++)
    {
        if (pos != FileValueMap.begin())
        {
            message += ", ";
        }
        message += "File '" + Fs::basename(pos->first) + "' Value: '" + pos->second + "'";
    }
    m_valueMismatchReportContent.push_back(message);
}

void AnnotationConverterEngine::MissingValueReport(string VarName, vector<string> Files)
{
    string message = "Variable Name: '" + VarName + "' File(s): ";
    
    vector<string>::iterator pos;
    for (pos = Files.begin(); pos != Files.end(); pos++)
    {
        if (pos != Files.begin())
        {
            message += ", ";
        }
        message += "'" + Fs::basename(*pos) + "'";
    }
    m_missingValueReportContent.push_back(message);
}




void CopyFile(string in, string out) 
{
    
	std::ifstream ins(in.c_str(), std::ios::binary);
	std::ofstream outs(out.c_str(), std::ios::binary);
	outs << ins.rdbuf();
	ins.close();
	outs.close();
}

    
void AnnotationConverterEngine::checkDiskSpaceImp() 
{

    std::string dbTemplateUncPath = Fs::convertToUncPath(m_dbTemplate);    
    std::string dbFileUncPath = Fs::convertToUncPath(m_dbFileName);
    
    CopyFile(dbTemplateUncPath, dbFileUncPath);
    
    if (!Fs::fileExists(dbFileUncPath))
    {
        Err::errAbort("Cannot copy the template file '" + dbTemplateUncPath + "' to the DB file '" + dbFileUncPath + "'");        
    }
    
    std::string fullPath = Fs::convertToUncPath(m_dbFileName, true);
    string dirname = Fs::dirname(fullPath);    
    
    int64_t out_disk_available = Fs::getFreeDiskSpace(dirname);      
    
    if (out_disk_available != -1)
    {// If getAvailableDiskSpace could not determine disk space it returns -1. Skip verification
        int64_t requiredDiskSpace = 0;
        std::vector<std::string>::iterator fileIterator;    
        for (fileIterator = m_listOfTsvFiles.begin(); fileIterator != m_listOfTsvFiles.end(); fileIterator++)
        {
            std::string fileName = *fileIterator; 
            requiredDiskSpace += Fs::fileSize(fileName, false);
        }    
        
        requiredDiskSpace = requiredDiskSpace * 2.3; //  Worst case scenario: 2 – ASCII file UTF16 db, plus 0.3 "safety buffer"
        
        if (out_disk_available < requiredDiskSpace)
        {
            std::stringstream message;
            message << "Your disc has insufficient space. You need to have at least " << requiredDiskSpace << " bytes of free space.";

            Err::errAbort(message.str());
        }
    }
}

/*  Commenting out - appears to be windows specific, and doesn't appear to be used anywhere -AK
std::string AnnotationConverterEngine::GetErrorMsg()
{
	void* lpMsgBuf;
	if (!FormatMessage( 
			FORMAT_MESSAGE_ALLOCATE_BUFFER | 
			FORMAT_MESSAGE_FROM_SYSTEM | 
			FORMAT_MESSAGE_IGNORE_INSERTS,
			NULL,
			GetLastError(),
			MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
			(LPTSTR) &lpMsgBuf,
			0,
			NULL ))
	{
		// Handle the error.
		return "";
	}

	std::string message = (char*)lpMsgBuf;

	// Free the buffer.
	LocalFree( lpMsgBuf );

	return message;
}

*/

void AnnotationConverterEngine::runImp() 
    {        
    AptVerboseWrapper::progressBegin(1, "Loading configuration files\n", 2, 1, 2);
    string congigName = programPath + ".config";
    if (Fs::fileExists(congigName)) 
    {
      LoadConfiguration(congigName);
    }
    else
    {
        AptVerboseWrapper::out(1, "Could not find configuration file " + congigName);
    }
    AptVerboseWrapper::progressStep(1);// Found or not we need to report a step since the UI tool expects fixed number of steps
    
    LoadConfiguration(m_arraySpecificConfig);
    AptVerboseWrapper::progressStep(1);
    AptVerboseWrapper::progressEnd(1, "\nDone");
    
    
    sqlite3* m_annotationDb = NULL;
    
    
    try
    {        
        AptVerboseWrapper::progressBegin(1, "Loading the database\n", 1, 1, 1);
        
        m_annotationDb = OpenDb(m_dbFileName);                            
        
        m_localizationTable = new LocalizationTable(m_annotationDb,&m_configuration);
        m_localizationTable->Load();    
        
        m_chromosomeTable = new ChromosomeTable(m_annotationDb);
        m_chromosomeTable->Load();
        
        m_cdfInformationTable = new CdfInformationTable(m_annotationDb,&m_configuration);        
        m_cdfInformationTable->Update();
        
        m_annotationTable = new AnnotationTable(m_annotationDb, &m_configuration, m_localizationTable, m_chromosomeTable);
        m_annotationTable->LoadColumntTypes();
        
        m_csvHeadersToInformationTablePtr = new CsvHeadersToInformationTable(m_annotationDb, &m_configuration);
        m_csvHeadersToInformationTablePtr->LoadInformationTable();
        
        AptVerboseWrapper::progressStep(1);
        AptVerboseWrapper::progressEnd(1, "\nDone");
        
        AnalyzeTsvFiles();
    }
    catch (std::exception& exception)
    {   
        CloseDbAndRemoveFile();
	    Err::errAbort(exception.what());	    
    }           
        
    
    BeginTransaction();    
    bool success = false;
    try
    {
        success = m_annotationTable->InsertTsvFiles(m_listOfTsvFiles);        
    }
    catch(std::exception& exception)
    {
         CloseDbAndRemoveFile();
	     Err::errAbort(exception.what());
    }
    if (success)
    {    
        
        AptVerboseWrapper::progressBegin(1, "Completing the database\n", 1, 1, 1);
        CommitTransaction();        
        try
        {            
        
            BeginTransaction();
            m_csvHeadersToInformationTablePtr->UpdateInformationTable();
            CommitTransaction();
        
            sqlite3_close(m_annotationDb);                    
            AptVerboseWrapper::progressStep(1);            
            AptVerboseWrapper::progressEnd(1, "\nDone");            
            AptVerboseWrapper::out(1, "\nProgram finished", true);            
        }
        catch (std::exception& exception)
        {        
            CloseDbAndRemoveFile();
	        Err::errAbort(exception.what());
	    }
    }
    else
    {
        CloseDbAndRemoveFile();
        Err::errAbort("See log file for details.");
    }            
}

