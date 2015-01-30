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

#include "CdfInformationTable.h"

CdfInformationTable::CdfInformationTable(sqlite3* AnnotationDb, Configuration* ConfigPtr)
{
    m_annotationDb = AnnotationDb;
    m_configPtr = ConfigPtr;
}
    
void CdfInformationTable::Update()
{
    if (m_configPtr->m_cdfInformationTable.size() > 0)
    {        
        sqlite3_stmt *insertStatement;

        std::stringstream  commandTextColumns;
        commandTextColumns << "INSERT INTO CdfInformation (";
        
        std::stringstream  commandTextValues;
        commandTextValues << " VALUES (";

        map<string, string>::iterator cdfTableRow;
        int counter = 0;
        for ( cdfTableRow = m_configPtr->m_cdfInformationTable[0].begin(); cdfTableRow != m_configPtr->m_cdfInformationTable[0].end(); cdfTableRow++)
        {
            if (counter > 0)
            {
                commandTextColumns << ",";
                commandTextValues << ",";
            }

            commandTextColumns << "[" << cdfTableRow->first/*column name*/ << "]";        
            commandTextValues << "?";

            counter++;
        }

        commandTextColumns << ")";
        commandTextValues << ")";

        std::string insertSqlStatement =  commandTextColumns.str() + commandTextValues.str();
        int result = sqlite3_prepare_v2(m_annotationDb,  insertSqlStatement.c_str(), -1, &insertStatement, NULL);

        if (result != SQLITE_OK)
        {       
            sqlite3_finalize(insertStatement);              
            throw Except("Internal error. Failed to prepare INSERT INTO CdfInformation statement. Make sure that table exists, and columns in the table and in the configuration file CdfInformationTableRow node match each other.");
        }
        
        for (unsigned int i = 0; i < m_configPtr->m_cdfInformationTable.size(); i++)
        {
            int bindingIndex = 1;
            for ( cdfTableRow = m_configPtr->m_cdfInformationTable[i].begin(); cdfTableRow != m_configPtr->m_cdfInformationTable[i].end(); cdfTableRow++)
            { 
                if(sqlite3_bind_text(insertStatement, bindingIndex, cdfTableRow->second.c_str(), cdfTableRow->second.length(), SQLITE_TRANSIENT) != SQLITE_OK)
                {                
                    std::stringstream message;            
                    message << "Could not bind parameters for CdfInformation table update! column name: '" << cdfTableRow->first << "' column value: '" << cdfTableRow->second << "'";
                    throw Except(message.str().c_str());
                }
                bindingIndex++;
            }
            int result = sqlite3_step(insertStatement);
            if (result != SQLITE_DONE) 
            {
                sqlite3_finalize(insertStatement);
                std::stringstream message;            
                if (result == SQLITE_CONSTRAINT)   
                {//* Contraint violation 
                    message << "Error executing SQL statement during CdfInformation table update due to constraint violation!\n";
                    message << "Make sure CdfInformationTableRow nodes in both configuration files have no duplicated values for unique fields such as the filename or the guid.\n"; 
                    message << "Check the CdfInformation table for other constraints.\n\n";
                }
                else
                {
                    message << "Could not execute SQL statement during CdfInformation table update! Error code: " << result << "\n";
                }
                message << "CdfInformationTableRow causes the problem:\n";
                for ( cdfTableRow = m_configPtr->m_cdfInformationTable[i].begin(); cdfTableRow != m_configPtr->m_cdfInformationTable[i].end(); cdfTableRow++)
                {
                    message << "Attrib:'" << cdfTableRow->first.c_str() << "' Value'" << cdfTableRow->second.c_str() << "'\n";
                } 
                message  << "\n";
                
                throw Except(message.str().c_str());
            }
            sqlite3_reset(insertStatement);
            sqlite3_clear_bindings(insertStatement);
        }
        sqlite3_finalize(insertStatement);
    }
}
