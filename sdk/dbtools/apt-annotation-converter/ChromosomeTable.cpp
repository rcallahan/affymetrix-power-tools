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

#include "ChromosomeTable.h"

map<string, int64_t> ChromosomeTable::m_shortNameToKey; 

ChromosomeTable::ChromosomeTable(sqlite3* AnnotationDb)
{
    m_annotationDb = AnnotationDb;
}

int ChromosomeTable::ReadTableCallback(void *NotUsed, int argc, char **argv, char **azColName)
{   

    string shortName = argv[0];
    int64_t key;
    std::stringstream iss(argv[1]);
    if ((iss >> key).fail())
    {
        std::stringstream message;
        message << "Error loading Chromosome table. Key '" << argv[1] << "' cannot be converted to integer.";
        throw Except(message.str().c_str());
    }
    
    map<string, int64_t>::iterator shortNameToKeyIterator = m_shortNameToKey.find(shortName);
    if (shortNameToKeyIterator == m_shortNameToKey.end())
    {                      
        m_shortNameToKey.insert(std::make_pair(shortName,key));    
    }
    else
    {
        m_shortNameToKey[shortName] = key;
    }

    return 0;
}



void ChromosomeTable::Load()
{
    
    char *errorMessage = NULL;
    
    m_shortNameToKey.clear();
    m_shortNameToKey.insert(std::make_pair("---", 2147483648)); // "---" in CSV files for Chromosome column represents not a number
    m_shortNameToKey.insert(std::make_pair("2147483648", 2147483648));
    m_shortNameToKey.insert(std::make_pair("DBNULL", 2147483648));
    m_shortNameToKey.insert(std::make_pair("NA", 2147483648));
    
    
      if(sqlite3_exec(m_annotationDb, "Select shortname, key from Chromosome", ChromosomeTable::ReadTableCallback, 0, &errorMessage) != SQLITE_OK)
      {
        std::string message = "Cannot read data from the Chromosome table! ";
        if (errorMessage != NULL)
        {
            message += errorMessage;
            delete(errorMessage);
        }
        throw Except(message.c_str());
      }   
    
}
int64_t ChromosomeTable::GetKey(string ShortName)
{
    map<string, int64_t>::iterator shortNameToKeyIterator = m_shortNameToKey.find(ShortName);
    if (shortNameToKeyIterator != m_shortNameToKey.end())
    {
        return shortNameToKeyIterator->second;
    }
    else
    {   
        if (ShortName == "M") 
        {// CSVTOSQLITE-344 Error occurs for "M" value conversion for Chromosome.
            return GetKey("MT");
        }
        std::string message = "Could not find a key for the Short Name '" + ShortName + "' in the Chromosome Table";
        throw Except(message.c_str());
    }    
}
